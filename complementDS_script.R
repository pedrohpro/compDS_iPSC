
setwd("/mnt/sandbox-SSD-1/poliveira/bruna/scDown")
getwd() # checking

rm(list = ls())     # clears Environment 


########################
### loading library  ###
########################

library(devtools)
library(tidyverse) #do it all
library(Seurat) # single cell analysis
library(data.table)
library(org.Hs.eg.db) #Homo sapiens OrgDb
library(biomaRt) #gene annotations
library(RColorBrewer) # for pretty colors
library(ggrepel)
library(umap)
library(VennDiagram) #venn digrams
library(clusterProfiler) #functional analysis
library(pathview)
library(cowplot)
library(GOplot)
library(enrichplot)
library(gage)
library(vsn)
library(RColorBrewer) # colors pallets
library(viridis) # color pallets
library(scales) # colors pallets
library(metap) # for pvalues
library(rsvd) # for ALRA
library(ALRA) # zero imputation
library(SeuratWrappers) # third party seurat functions (eg.: RunALRA())
library(DESeq2)
library(patchwork)
library(scDblFinder) # doublet finder
library(SingleCellExperiment)
library(AUCell) # get set score
library(GSEABase) # build gene set function
library(sessioninfo) # package versions

source("chromopack.R")


############################
### creating directories ###
############################

if (!file.exists("results")){   
  dir.create("results")
  dir.create("results/exploratory_analysis")
}

############################
###         counts       ###
############################

# Runs names lists (maybe better to get this from meta data)
#runs <- sapply(43:71, function(i) paste0("EGAN000033605", i), simplify = T) # for all samples
runs <- c("EGAN00003360559", "EGAN00003360560", "EGAN00003360561", "EGAN00003360562", "EGAN00003360563", "EGAN00003360564", # young samples
          "EGAN00003360565", "EGAN00003360566", "EGAN00003360567", "EGAN00003360568", "EGAN00003360569", "EGAN00003360570")

# Raw counts
counts <- list()
for (run in runs){ 
  counts[[run]] <- Read10X(paste0("data/counts/",run,"/"))
}

# Separate objects
seurats <- list()
for (run in names(counts)) {
  seurats[[run]] <- CreateSeuratObject(counts = counts[[run]], project = run)
}

# merging seurat objects (without batch correction)
scdown <- merge(
  x = seurats$EGAN00003360559, # change the starting sample according to samples used
  y = seurats[names(seurats) != "EGAN00003360559"], # change the starting sample according to samples used
  add.cell.ids = runs,
  project = "allRuns"
)

# joining layers
scdown <- JoinLayers(scdown)

############################
###      meta data       ###
############################

meta <- read.delim("data/meta.csv", sep = ",") %>% 
  select(accession_id, biological_sex, subject_id, phenotype) %>%  # filtering columns
  rename( # renaming columns
    condition = phenotype,
    sex = biological_sex,
    donor = subject_id
  ) %>% 
  mutate(
    condition = case_when( # renaming elements
      condition == "non-diseased" ~ "CT",
      condition == "Down Syndrome" ~ "DS"
    ),
    age = substr(donor, start = nchar(donor) - 1, stop = nchar(donor) - 1) # adding age column
  ) %>% 
  filter( # getting only young samples to avoid the impact of Alzheimer's Disease neuroinflammation
    age == "Y"
  ) %>% 
  mutate_all(factor) # factoring all columns


# Including the sample meta data into the merged cell meta data
scdown@meta.data <- scdown@meta.data %>% 
  rownames_to_column("cell") %>% 
  left_join(meta, by = c("orig.ident" = "accession_id")) %>% 
  column_to_rownames("cell")

# Creating mitochondrial reads percentage
scdown[["percent.mt"]] <- PercentageFeatureSet(scdown, pattern = "^MT-") # high percentage could indicate either cell death or high metabolism cell type


#############################
###   genes of interest   ###
#############################

# genes of interest
int_genes_df <- read.delim("data/genes_down.tsv", sep = "\t")
int_genes <- int_genes_df$genesofinterest
comp_genes <- int_genes_df$complement[!is.na(int_genes_df$complement)]

# cell type markers
markers_df <- read.delim("data/celltypes_markers.csv", sep = ";")
celltypes_markers <- markers_df$marker


############################
###   filtering cells    ###
############################

# Plotting: features, reads, and percent mt reads before QC
violin <- VlnPlot(scdown, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2, pt.size = 0.1, alpha = 0.05)
ggsave("results/exploratory_analysis/before_violin.png", plot = violin, height = 9, width = 9)

# Plotting: nCount vs nFeatures and nCount vs percent.mt before QC
plot1 <- FeatureScatter(scdown, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot3 <- FeatureScatter(scdown, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot4 <- FeatureScatter(scdown, feature1 = "nFeature_RNA", feature2 = "percent.mt")
combined_plot <- (plot1)/(plot3|plot4)
ggsave("results/exploratory_analysis/before_scatter.png", plot = combined_plot, height = 8, width = 11)


################ Filtering
# Number of reads per cell
scdown <- subset(scdown, subset = nCount_RNA > 200 & nCount_RNA < 60000)

# Number of features per cell
scdown <- subset(scdown, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000)

# Percentage of reads mapped to mitochondrial genes 
scdown <- subset(scdown, subset = percent.mt < 10)

# doublets detection
# Convert Seurat object to SingleCellExperiment object
sce <- as.SingleCellExperiment(scdown)
# Run scDblFinder
sce <- scDblFinder(sce)
# Add doublet information to Seurat object metadata
scdown$doublet <- sce$scDblFinder.class
# Filtering out doublets
scdown <- subset(scdown, subset = doublet == "singlet")


# Plotting: features, reads, and percent mt reads after QC
violin <- VlnPlot(scdown, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2, pt.size = 0.1, alpha = 0.05)
ggsave("results/exploratory_analysis/after_violin.png", plot = violin, height = 9, width = 9)

# Plotting: nCount vs nFeatures and nCount vs percent.mt after QC
plot1 <- FeatureScatter(scdown, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot3 <- FeatureScatter(scdown, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot4 <- FeatureScatter(scdown, feature1 = "nFeature_RNA", feature2 = "percent.mt")
combined_plot <- (plot1)/(plot3|plot4)
ggsave("results/exploratory_analysis/after_scatter.png", plot = combined_plot, height = 8, width = 11)


###############################
###   filtering features    ###
###############################

# Number of features before quality control
print(paste0("Number of features before QC: ", nrow(scdown)))

# Removing mitochondrial genes
scdown <- scdown[!grepl("MT-", rownames(scdown)),]

# features information
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

all_features <- getBM(
  attributes = c("entrezgene_id", "external_gene_name", "gene_biotype", "chromosome_name", "start_position", "end_position"),
  filters = "external_gene_name",
  values = rownames(scdown),
  mart = ensembl
)

all_features <- all_features %>% 
  filter(
    complete.cases(.), # dropping features with NAs
    gene_biotype == "protein_coding", # only protein_coding genes
    chromosome_name %in% c(as.character(seq(1, 22))), # only chromosomal genes # c(as.character(seq(1, 22)), "X", "Y")
    !duplicated(external_gene_name) | !duplicated(external_gene_name, fromLast = TRUE) # choosing one of duplicated gene names
  ) %>% 
  mutate(
    chromosome_name = factor(chromosome_name, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")), #"X","Y"
    gene_length = end_position - start_position,
    avg_position = (end_position + start_position)/2
  )

# filtering features from seurat object based on aforementioned filters
scdown <- scdown[rownames(scdown) %in% all_features$external_gene_name,]

# Removing genes with zero expression
aux <- as.data.frame(as.matrix(scdown[["RNA"]]$counts)) # sparse matrix to normal (zero inputation)
aux <- rownames(aux[rowSums(aux) == 0,]) # instead of using rm() changing aux to be smaller
scdown <- scdown[!(rownames(scdown) %in% aux),]
all_features <- all_features[all_features$external_gene_name %in% rownames(scdown),]

# Number of features after quality control
print(paste0("Number of features after QC: ", nrow(scdown)))



######################
###   Processing   ###
######################

scdown <- NormalizeData(scdown)

scdown <- RunALRA(scdown, assay = "RNA") # input is data (normalized) layer from RNA assay. Output is data layer from alra assay
imputed_data <- GetAssayData(scdown, layer = "data", assay = "alra") # Extract the imputed data from the 'alra' assay
scdown@assays[["RNA"]]@layers[["data"]] <- imputed_data # Update the 'data' layer of RNA assay

scdown <- FindVariableFeatures(scdown) # default 2000 features

top10 <- head(VariableFeatures(scdown), 10) # Highlight the 10 most highly variable genes

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(scdown) + ggtitle("All genes")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + ggtitle("Top 10 highly variable genes")
plot3 <- LabelPoints(plot = plot1, points = int_genes, repel = TRUE) + ggtitle("Genes associated with DS")
ggsave("results/exploratory_analysis/variableFeatures.png", plot = plot1 + plot2 + plot3, width = 15, height = 5)

scdown <- ScaleData(scdown, features = rownames(scdown))

scdown <- RunPCA(scdown, features = VariableFeatures(object = scdown)) 

elbow <- ElbowPlot(scdown, ndims = 50) + # help define number of PCs 
  theme(
    plot.background = element_rect(fill = "white"), 
    panel.background = element_rect(fill = "white"), 
    panel.border = element_blank(), 
    axis.line = element_line(color = "black"),  
    legend.background = element_rect(fill = "white"),  
  )
ggsave("results/exploratory_analysis/elbow.png", plot = elbow, height = 4, width = 4)

scdown <- RunUMAP(scdown, dims = 1:10)

scdown <- FindNeighbors(scdown, dims = 1:10, k.param = 30)

scdown <- FindClusters(scdown, resolution = 0.2)



#########################################
###      Batch effect verification    ###
#########################################

# UMAP
for (i in c("seurat_clusters","donor","sex","condition")) {
  clusters <- DimPlot(scdown, reduction = "umap", group.by = i, pt.size = 0.1, label = TRUE, raster=FALSE) +
    ggtitle("") +
    labs(color = i) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    coord_fixed() +
    theme(
      legend.text = element_text(size = 26),  # Adjust legend text size and style
      legend.title = element_text(size = 28, face = "bold"),  # Adjust legend text size and style
      axis.text = element_text(size = 16),    # Adjust axis text size and style
      axis.title = element_text(size = 22, face = "bold"),   # Adjust axis title size and style
      axis.line = element_line(linewidth = 1.5), # Customize axis lines
      axis.ticks = element_line(linewidth = 1.5),                 # Adjust axis ticks size
    )
  
  ggsave(paste0("results/exploratory_analysis/umap_",i,".png"), plot = clusters, width = 12, height = 12)
}

# Bar plots
for (i in c("donor","sex","condition")) {
  aux <- scdown@meta.data %>% 
    select(seurat_clusters, !!sym(i)) %>% 
    group_by(seurat_clusters, !!sym(i)) %>% 
    summarise(num_cells = n())%>%
    group_by(seurat_clusters) %>%
    mutate(pct = num_cells / sum(num_cells) * 100) %>% 
    ungroup()
  
  # absolute value bar plot
  bar_plot <- barfunc(aux, x_axis = "seurat_clusters", y_axis = "num_cells", fill_groups = i)
  ggsave(paste0("results/exploratory_analysis/bar_",i,".png"), plot = bar_plot, width = 5, height = 4)
  
  # pct value bar plot
  bar_plot <- barfunc(aux, x_axis = "seurat_clusters", y_axis = "pct", fill_groups = i, annot_loc = "center", annot_text = "num_cells")
  ggsave(paste0("results/exploratory_analysis/pctbar_",i,".png"), plot = bar_plot, width = 5, height = 4)
}

#################################
###      Bias verification    ###
#################################

# condition vs. sex
aux <- scdown@meta.data %>% 
  select(condition, sex) %>% 
  group_by(condition, sex) %>% 
  summarise(num_cells = n())%>%
  group_by(condition) %>%
  mutate(pct = num_cells / sum(num_cells) * 100) %>% 
  ungroup()

# absolute value bar plot
bar_plot <- barfunc(aux, x_axis = "condition", y_axis = "num_cells", fill_groups = "sex")
ggsave("results/exploratory_analysis/abs_condition_sex.png", plot = bar_plot, width = 3, height = 3)

# pct value bar plot
bar_plot <- barfunc(aux, x_axis = "condition", y_axis = "pct", fill_groups = "sex", annot_loc = "center", annot_text = "num_cells")
ggsave("results/exploratory_analysis/pct_condition_sex.png", plot = bar_plot, width = 3, height = 3)


####################################
###       Finding markers        ###
####################################

# Verifying assay type. Should be RNA
DefaultAssay(scdown)

# Changing identity to cluster column
Idents(scdown) <- "seurat_clusters"

# hline annotation
auxmark <- markers_df %>% 
  group_by(celltype) %>% 
  summarise(n_markers = n()) %>% 
  mutate(celltype = factor(celltype, levels = unique(markers_df$celltype))) %>% 
  arrange(celltype) %>%  
  mutate(
    y = length(unique(scdown@meta.data$seurat_clusters))*1.05,
    yend = y,
    x = 1, # (just for the loop to work)
    xend = n_markers, # provisory (just for the loop to work)
    colour = "#000000",
    celltype = case_when(
      celltype == "Ast" ~ "Astrocyte",
      celltype == "End" ~ "Endothelial",
      celltype == "Exc" ~ "Excitatory neuron",
      celltype == "Inh" ~ "Inhibitory neuron",
      celltype == "Mic" ~ "Microglia",
      celltype == "Oli" ~ "Oligodendrocyte",
      celltype == "OPC" ~ "Oligodendrocyte Progenitor",
      celltype == "Per" ~ "Pericyte",
      celltype == "NPC" ~ "Neural Progenitor"
    )
  ) %>% 
  as.data.frame()
for (i in 2:nrow(auxmark)) {
  auxmark[i,"x"] <- auxmark[i-1,"xend"] + 1
  auxmark[i,"xend"] <- auxmark[i,"x"] + auxmark[i,"n_markers"] -1
}

# canonical markers
dot_plot <- DotPlot(
  object = scdown, 
  features = celltypes_markers, 
  cols = c("blue","red"),
  cluster.idents = T,
  scale = T 
)+
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_blank(),  # Remove panel borders
    axis.line = element_line(color = "black"),  # Set axis line color
    legend.background = element_rect(fill = "white"),  # Set legend background color to white
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
    axis.text.y = element_text(size = 17, face = "bold"),
    #legend.position = "none" # no legend
  )
for (i in 1:nrow(auxmark)) {
  dot_plot <- dot_plot + 
    annotate(
      "segment", # type of annotation to be added
      x = auxmark[i,"x"],
      y = auxmark[i,"y"],
      xend = auxmark[i,"xend"],
      yend = auxmark[i,"yend"],
      #colour = auxmark[i,"colour"],
      
    ) +
    annotate(
      "text", # type of annotation to be added
      x = (auxmark[i,"x"] + auxmark[i,"xend"])/2, 
      y = auxmark[i,"y"] + 0.15, 
      label = auxmark[i,"celltype"], 
      #hjust = 1.1, 
      #vjust = -0.5, 
      color = "black",
      fontface = "bold",
      size = 4
    )
}

ggsave("results/exploratory_analysis/canonical_markers.png", plot = dot_plot, height = 5, width = 22)



###############################
###       Annotating        ###
###############################
Idents(scdown) <- "seurat_clusters"

# # Without ALRA
# scdown <- RenameIdents(scdown, c(
#   `0` = "Oli",
#   `1` = "Ast",
#   `2` = "Exc",
#   `3` = "OPC",
#   `4` = "Mic",
#   `5` = "Exc",
#   `6` = "Inh",
#   `7` = "Inh",
#   `8` = "Exc",
#   `9` = "End",
#   `10` = "End",
#   `11` = "Ast",
#   `12` = "OPC"
# )
# )
# 
# scdown@meta.data <- scdown@meta.data %>%
#   mutate(
#     seurat_clusters = case_when(
#       seurat_clusters == "0" ~ "Oli",
#       seurat_clusters == "1" ~ "Ast",
#       seurat_clusters == "2" ~ "Exc",
#       seurat_clusters == "3" ~ "OPC",
#       seurat_clusters == "4" ~ "Mic",
#       seurat_clusters == "5" ~ "Exc",
#       seurat_clusters == "6" ~ "Inh",
#       seurat_clusters == "7" ~ "Inh",
#       seurat_clusters == "8" ~ "Exc",
#       seurat_clusters == "9" ~ "End",
#       seurat_clusters == "10" ~ "End",
#       seurat_clusters == "11" ~ "Ast",
#       seurat_clusters == "12" ~ "OPC"
#     )
#   )

# With ALRA
scdown <- RenameIdents(scdown, c(
  `0` = "Oli",
  `1` = "Exc",
  `2` = "Oli",
  `3` = "Ast",
  `4` = "Inh",
  `5` = "OPC",
  `6` = "Mic",
  `7` = "End",
  `8` = "Exc",
  `9` = "Mic",
  `10` = "Mic",
  `11` = "Oli",
  `12` = "Oli",
  `13` = "Ast",
  `14` = "Exc",
  `15` = "Inh"
)
)

scdown@meta.data <- scdown@meta.data %>%
  mutate(
    seurat_clusters = case_when(
      seurat_clusters == "0" ~ "Oli",
      seurat_clusters == "1" ~ "Exc",
      seurat_clusters == "2" ~ "Oli",
      seurat_clusters == "3" ~ "Ast",
      seurat_clusters == "4" ~ "Inh",
      seurat_clusters == "5" ~ "OPC",
      seurat_clusters == "6" ~ "Mic",
      seurat_clusters == "7" ~ "End",
      seurat_clusters == "8" ~ "Exc",
      seurat_clusters == "9" ~ "Mic",
      seurat_clusters == "10" ~ "Mic",
      seurat_clusters == "11" ~ "Oli",
      seurat_clusters == "12" ~ "Oli",
      seurat_clusters == "13" ~ "Ast",
      seurat_clusters == "14" ~ "Exc",
      seurat_clusters == "15" ~ "Inh"
    )
  )

# ordering in decreasing number of cells
scdown@meta.data$seurat_clusters <- factor(scdown@meta.data$seurat_clusters, levels = c("Oli","Exc","Ast","OPC","Inh","Mic","End"))


###################################
###       Annotated UMAP        ###
###################################

umap_plot <- DimPlot(
  scdown, 
  reduction = "umap", 
  label = F, 
  pt.size = 0.1, 
  raster = FALSE, 
  label.size = 12, 
  cols = c("Oli" = brewer.pal(7, "Set3")[1], "Exc" = brewer.pal(8, "Set3")[8], "OPC" = brewer.pal(7, "Set3")[3], 
           "Inh" = brewer.pal(7, "Set3")[4], "Ast" = brewer.pal(7, "Set3")[7], "Mic" = brewer.pal(7, "Set3")[6], 
           "End" = brewer.pal(7, "Set3")[5])
) +
  NoAxes() + # removes axis
  theme(
    legend.position = "none"
  )
ggsave("results/exploratory_analysis/ANNOTATED_umap.png", plot = umap_plot, height = 12, width = 12, dpi = 600)


####################################################
###      Repeating Batch effect verification     ###
####################################################

# Bar plots
for (i in c("donor","sex","condition")) {
  aux <- scdown@meta.data %>% 
    select(seurat_clusters, !!sym(i)) %>% 
    group_by(seurat_clusters, !!sym(i)) %>% 
    summarise(num_cells = n())%>%
    group_by(seurat_clusters) %>%
    mutate(pct = num_cells / sum(num_cells) * 100) %>% 
    ungroup()
  
  # absolute value bar plot
  bar_plot <- barfunc(aux, x_axis = "seurat_clusters", y_axis = "num_cells", fill_groups = i)
  ggsave(paste0("results/exploratory_analysis/ANNOTATED_bar_",i,".png"), plot = bar_plot, width = 5, height = 4)
  
  # pct value bar plot
  bar_plot <- barfunc(aux, x_axis = "seurat_clusters", y_axis = "pct", fill_groups = i, annot_loc = "center", annot_text = "num_cells")
  ggsave(paste0("results/exploratory_analysis/ANNOTATED_pctbar_",i,".png"), plot = bar_plot, width = 5, height = 4)
}


##################################################
###          Repeating markers plots           ###
##################################################

Idents(scdown) <- "seurat_clusters"

# hline annotation
auxmark <- markers_df %>% 
  group_by(celltype) %>% 
  summarise(n_markers = n()) %>% 
  mutate(celltype = factor(celltype, levels = unique(markers_df$celltype))) %>% 
  arrange(celltype) %>%  
  mutate(
    y = length(unique(scdown@meta.data$seurat_clusters)) + 0.25,
    yend = y,
    x = 1, # (just for the loop to work)
    xend = n_markers, # provisory (just for the loop to work)
    colour = "#000000",
    celltype = case_when(
      celltype == "Ast" ~ "Astrocyte",
      celltype == "End" ~ "Endothelial",
      celltype == "Exc" ~ "Excitatory neuron",
      celltype == "Inh" ~ "Inhibitory neuron",
      celltype == "Mic" ~ "Microglia",
      celltype == "Oli" ~ "Oligodendrocyte",
      celltype == "OPC" ~ "Oligodendrocyte Progenitor",
      celltype == "Per" ~ "Pericyte",
      celltype == "NPC" ~ "Neural Progenitor"
    )
  ) %>% 
  as.data.frame()
for (i in 2:nrow(auxmark)) {
  auxmark[i,"x"] <- auxmark[i-1,"xend"] + 1
  auxmark[i,"xend"] <- auxmark[i,"x"] + auxmark[i,"n_markers"] -1
}


# canonical markers
dot_plot <- DotPlot(
  object = scdown, 
  features = celltypes_markers, 
  cols = c("blue","red"),
  cluster.idents = T,
  #cluster.features = T, # não existe!!!
  scale = T # zscore 
)+
  theme(
    plot.background = element_rect(fill = "white"),  # Set background color to white
    panel.background = element_rect(fill = "white"),  # Set panel background color to white
    #panel.grid.major = element_blank(),  # Remove major grid lines
    #panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove panel borders
    axis.line = element_line(color = "black"),  # Set axis line color
    legend.background = element_rect(fill = "white"),  # Set legend background color to white
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
    axis.text.y = element_text(size = 17, face = "bold"),
    #legend.position = "none" # no legend
  )
for (i in 1:nrow(auxmark)) {
  dot_plot <- dot_plot + 
    annotate(
      "segment", # type of annotation to be added
      x = auxmark[i,"x"],
      y = auxmark[i,"y"],
      xend = auxmark[i,"xend"],
      yend = auxmark[i,"yend"],
      #colour = auxmark[i,"colour"],
      
    ) +
    annotate(
      "text", # type of annotation to be added
      x = (auxmark[i,"x"] + auxmark[i,"xend"])/2, 
      y = auxmark[i,"y"] + 0.15, 
      label = auxmark[i,"celltype"], 
      #hjust = 1.1, 
      #vjust = -0.5, 
      color = "black",
      fontface = "bold",
      size = 4
    )
}

ggsave("results/exploratory_analysis/ANNOTATED_dotplot.png", plot = dot_plot, height = 8, width = 22)


###########################################
###        Cell type composition        ###
###########################################

# Absolute value
aux <- scdown@meta.data %>%
 select(orig.ident, seurat_clusters, condition) %>%
 group_by(orig.ident, condition, seurat_clusters) %>%
 summarise(num_cells = n()) %>%
 group_by(orig.ident) %>%
 mutate(pct = num_cells / sum(num_cells) * 100) %>%
 ungroup()

# separate CT and DS horizontally
rowOrder <- meta %>%
 arrange(condition) %>%
 pull(accession_id)
aux$orig.ident <- factor(aux$orig.ident, levels = rowOrder)

# absolute value plot
abs_composition <- barfunc(
 aux, x_axis = "orig.ident", y_axis = "num_cells", fill_groups = "seurat_clusters",
 fill_col = c(brewer.pal(n = 10, name = "Set3"), brewer.pal(n = nlevels(aux$seurat_clusters)-10, name = "Paired")),
 rotate_x = 45, size_xaxis_text = 8)
ggsave("results/exploratory_analysis/ANNOTATED_abs_composition.png", plot = abs_composition, height = 5, width = 6)

# Percentage value
pct_composition <- barfunc(
 aux, x_axis = "orig.ident", y_axis = "pct", fill_groups = "seurat_clusters", annot_loc = "center", annot_text = "num_cells",
 fill_col = c(brewer.pal(n = 10, name = "Set3"), brewer.pal(n = nlevels(aux$seurat_clusters)-10, name = "Paired")),
 rotate_x = 45, size_xaxis_text = 8)
ggsave("results/exploratory_analysis/ANNOTATED_pct_composition.png", plot = pct_composition, height = 5, width = 6)

aux <- aux %>%
 group_by(orig.ident, condition) %>%
 summarise(
   Ast_Mic_sum = sum(num_cells[seurat_clusters %in% c("Ast", "Mic")]),
   Exc_Inh_sum = sum(num_cells[seurat_clusters %in% c("Exc", "Inh")]),
   proportion = Ast_Mic_sum / Exc_Inh_sum
 )

##############################################
###    Differential expression analysis    ###
##############################################

# creating composite column
scdown@meta.data$type_cond <- paste0(scdown@meta.data$seurat_clusters, "_", scdown@meta.data$condition)

# Changing identity to composite column
Idents(scdown) <- "type_cond"

DElist <- list()

# Running differential expression (DS vs CT for each cell type)
for (i in levels(scdown$seurat_clusters)){
  DElist[[i]] <- FindMarkers(
    scdown,
    ident.1 = paste0(i,"_DS"), # find markers for this identity
    ident.2 = paste0(i,"_CT"), # against this identity
    logfc.threshold = 0.0, # to analyse every gene
    min.pct = 0.0, # to analyse every gene
    min.cells.feature = 0
  )
}

# adding DEG and chromosome information and cell type information
DElist <- imap(DElist, function(i, ct){ # i is a df and ct is list element name
  i <- i %>% 
    # adding DEG and alteration information and gene name column
    mutate( # creates new columns
      celltype = ct,
      gene = rownames(i),
      DEG = case_when(
        avg_log2FC >= 1 & p_val_adj < 0.05 ~ "UP",
        avg_log2FC <= -1 & p_val_adj < 0.05 ~ "DOWN",
        TRUE ~ "NO" # else none of the other conditions are true
      ),
      alterat = case_when(
        avg_log2FC > 0 & p_val_adj < 0.05 ~ "UP",
        avg_log2FC < 0 & p_val_adj < 0.05 ~ "DOWN",
        TRUE ~ "NO" # else none of the other conditions are true
      )
    ) %>% 
    # adding gene info from biomart
    left_join(all_features, by = c("gene" = "external_gene_name"))  # gene column from i and external_gene_name column from all_features
  return(i)
})

scdown@meta.data$type_cond <- factor(scdown@meta.data$type_cond, levels = c("Oli_DS","Oli_CT","Exc_DS","Exc_CT",
                                                                            "Inh_DS","Inh_CT","Ast_DS","Ast_CT",
                                                                            "OPC_DS","OPC_CT","Mic_DS","Mic_CT",
                                                                            "End_DS","End_CT"))

###############################
###       Violin plot       ###
###############################

norm_df <- as.data.frame(as.matrix(GetAssayData(scdown, layer = "data")))
features_to_plot <- c("C1R","C1S","C4A","C4B","C5","C8G","CD46","CD59","CFI","FCN2","SLC1A2","SLC1A3")

for (i in c("Ast")) { # names(DElist)
  # preping input
  aux <- norm_df[rownames(norm_df) %in% features_to_plot, scdown@meta.data$seurat_clusters == i]
  aux["condition",] <- scdown@meta.data[scdown@meta.data$seurat_clusters == i, "condition"]
  aux <- as.data.frame(t(aux))
  aux <- aux %>% 
    pivot_longer(
      cols = -condition,
      names_to = "gene",
      values_to = "expression"
    ) %>% 
    as.data.frame() %>% 
    mutate(expression = as.numeric(expression))
  
  # preping anot_df
  annot_df <- DElist[[i]][DElist[[i]]$gene %in% features_to_plot,]
  annot_df$avg_log2FC <- paste0("log2FC: ", round(annot_df$avg_log2FC, 2))
  annot_df$p_val_adj <- paste0("padj: ", sprintf("%.1e", annot_df$p_val_adj))
  annot_df$annot <- paste0(annot_df$p_val_adj, "\n", annot_df$avg_log2FC)
  annot_df$max_expression <- NA
  for(j in 1:nrow(annot_df)){
    annot_df$max_expression[j] <- max(aux[aux$gene == annot_df[j,"gene"],"expression"])
  }
  
  # Plotting
  vln_plot <- ggplot(aux, aes(x = gene, y = expression, fill = condition)) +
    geom_violin(
      position = position_dodge(width = 3),
      alpha = 0.5,
      width = 3.5
    )+ 
    facet_wrap( # separated box plots
      ~gene, 
      scales = "free",
      ncol = 6 #round(length(unique(aux$gene))/3)
    )+ 
    scale_fill_manual(values = c("DS" = "#ffcc00", "CT" = "#3771c8")) +  # Specify fill colors for DS and CT
    labs(
      x = NULL,
      y = NULL
    )+ 
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(), #element_line(color = "#cccccc99", linetype = "dotted"),
      strip.background = element_rect(fill = brewer.pal(7, "Set3")[7]), 
      axis.text.x = element_blank(), #removes bottom gene name
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(size = 12, face = "bold")
    )
  
  ggsave(paste0("results/DE_", i, ".pdf"), plot = vln_plot, height = 4, width = 10)
}

##################################
###       Expression UMAP      ###
##################################

##### calculating gene set score

# get expression matrix
expr_matrix <- GetAssayData(scdown, assay = "RNA", layer = "data")

# Create a GeneSet object
AstCompSet <- GeneSet(c("C1R","C1S","C4A","C4B","C5","C8G","CD46","CD59","CFI","FCN2"), setName = "AstCompSet")

# Calculate AUCell scores
cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE, verbose = FALSE)
cells_AUC_ast <- AUCell_calcAUC(AstCompSet, cells_rankings, aucMaxRank = nrow(expr_matrix)*0.05)

# Extract AUC scores and add to Seurat object metadata
auc_scores_ast <- as.data.frame(t(getAUC(cells_AUC_ast)))

# adding scores to metadata
scdown <- AddMetaData(scdown, auc_scores_ast, col.name = "AstCompSet")

# subseting
Idents(scdown) <- "seurat_clusters"
cells_of_interest <- WhichCells(scdown, ident = "Ast")
sc_subset <- subset(scdown, cells = cells_of_interest)

umap_feature <- FeaturePlot(
  sc_subset,
  features = "AstCompSet",
  pt.size = 0.3
) + 
  scale_color_gradient(low = "white", high = "red") + 
  facet_grid(.~sc_subset$condition) +
  labs(
    title = "C1R+C1S+C4A+C4B+C5+\nC8G+CD46+CD59+CFI+FCN2"
  ) +
  theme(
    strip.background = element_rect(fill = brewer.pal(7, "Set3")[7], color = "black"),  
    strip.text = element_text(color = "black", face = "bold"),   # Text color and style
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none" # no legend
  )
ggsave("results/Ast_feature_umap.png", plot = umap_feature, height = 3.5, width = 5, limitsize = FALSE, dpi = 600)


###############################
###      Volcano plot       ###
###############################

for (i in c("Mic","Ast")) {
  # Identifying top DEGs
  aux2 <- DElist[[i]][DElist[[i]]$p_val_adj < 0.05,]
  aux2 <- aux2[order(-abs(aux2$avg_log2FC)),] # orders in decreasing absolute value
  aux2 <- aux2$gene[1:10] # 10 biggest DEGs
  
  # adding highlight
  aux <- DElist[[i]] %>% 
    mutate(
      highlight = case_when(
        gene %in% c("C1QA","C1QB","C1QC","C2","C3AR1","C5AR1","C5AR2","C1R","C1S","C4A","C4B","C8G","SLC1A2","SLC1A3") ~ gene, # label degs and int_genes ############################## change C4A C4B for other genes!
        TRUE ~ NA # else
      ),
      p_val_adj = case_when( # the volcano plot gets too stretched with padj values that tend to zero, so this puts a ceiling
        p_val_adj < 1e-20 ~ 1e-20,
        TRUE ~ p_val_adj
      )
    )
  
  # Plotting
  volc_plot <- ggplot(data = aux, aes(x = avg_log2FC, y = -log10(p_val_adj), col = DEG, label = highlight)) +
    # Log2FC threshold line
    geom_vline(xintercept = c(-1, 1), col = "black", linetype = 'dashed') +
    # padj threshold line
    geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') +
    # dot size
    geom_point(size = 1.5) +
    # color
    scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"),labels = c("DOWN", "NO", "UP")) + 
    # zoom
    coord_cartesian(ylim = c(0, max(-log10(aux$p_val_adj))), xlim = c(min(aux$avg_log2FC), max(aux$avg_log2FC))) + 
    labs(
      color = "DEG",
      x = expression("log"[2]*"FC"), 
      y = expression("-log"[10]*"padj")
      #title = "Differentially expressed genes",
      #subtitle = "Comparing DS against CT (|log2fc| > 1, padj < 0.05)"
    ) +
    scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis (doted background lines)
    geom_text_repel(max.overlaps = Inf, color = "black")+ # To show all labels 
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_line(color = "#99999977", linetype = "dotted"),
      axis.text = element_text(size = 12, face = "bold")
    )
  
  ggsave(paste0("results/volc_", i, ".png"), plot = volc_plot, height = 5, width = 7)
}


###############################
###     DEGChromo Plot      ###
###############################

# Plotting
for (i in c("Mic","Ast","Exc","Inh")){
  
  degchromo <- DEGchromoPlot(
    DElist[[i]], 
    gene_col = "gene", 
    fc_col = "avg_log2FC", 
    x_axis = "chromosome_name",
    only_expr_features = T,
    color_dot_down = "#3771c8aa",
    color_dot_up = "#ffcc00aa",
    color_dot_no = "#dddddd33",
    color_bar_down = "#3771c866",
    color_bar_up = "#ffcc0066",
    color_text_down = "#3771c8",
    color_text_up = "#aa8800",
    color_score_down = "#3771c8",
    color_score_up = "#aa8800",
    size_xaxis_text = 18
  )
  
  ggsave(paste0("results/degchromo_", i, ".png"), plot = degchromo, height = 5, width = 12)
}

###################################
###      DEGs by cell type      ###
###################################

# combine all dfs
DEdf <- do.call(rbind, DElist)

# preping input
aux <- DEdf %>% 
  filter(DEG != "NO") %>% 
  group_by(celltype, DEG) %>% 
  summarise(total = n()) %>% 
  group_by(celltype) %>% #pct calculation
  mutate(pct = (total / sum(total))* 100) %>%  # new column
  ungroup() %>% 
  as.data.frame()

# ordering in decreasing order
ord <- aux %>% 
  filter(DEG == "UP") %>% 
  arrange(-total) %>%
  pull(celltype)

aux$celltype <- factor(aux$celltype, levels = c("Mic","Ast","Oli","End","OPC","Inh","Exc"))

abs_plot <- barfunc(aux, x_axis = "celltype", y_axis = "total", fill_groups = "DEG", annot_loc = "center")
ggsave("results/abs_DEGcelltype.png", plot = abs_plot, height = 5, width = 7)

pct_plot <- barfunc(aux, x_axis = "celltype", y_axis = "pct", fill_groups = "DEG", annot_loc = "center", annot_text = "total")
ggsave("results/pct_DEGcelltype.png", plot = pct_plot, height = 5, width = 7)


# DEGchromoPlot by celltype
DEdf$celltype <- factor(DEdf$celltype, levels = c("Mic","Ast","Oli","End","OPC","Inh","Exc")) # ALRA

degchromo <- DEGchromoPlot(
  DEdf, 
  gene_col = "gene", 
  fc_col = "avg_log2FC", 
  x_axis = "celltype",
  only_expr_features = T,
  color_dot_down = "#3771c8aa",
  color_dot_up = "#ffcc00aa",
  color_dot_no = "#dddddd33",
  color_bar_down = "#3771c866",
  color_bar_up = "#ffcc0066",
  color_text_down = "#3771c8",
  color_text_up = "#aa8800",
  color_score_down = "#3771c8",
  color_score_up = "#aa8800",
  size_xaxis_text = 18,
  size_score = 5
)
ggsave("results/degchromo_degs_by_celltype.png", plot = degchromo, height = 5, width = 6)


##############################################
###           Functional analysis          ###
##############################################

######################
###      ORA       ###
######################

# Running
ORA_ast <- enrichGO(
  keyType = "SYMBOL", #there is no preference for using SYMBOL, ENSEMBL or ENTREZ!
  gene = DElist$Ast %>% filter(DEG != "NO") %>% arrange(-abs(avg_log2FC)) %>% pull(gene),
  universe = rownames(scdown), # sample genes. IMPORTANT!!!
  OrgDb = org.Hs.eg.db,
  ont = "all",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

######################
###      Plots     ###
######################

# ontologies of interest
int_ont <- c("GO:0098883", "GO:0150062", "GO:1905805", "GO:1905806", "GO:1905807", "GO:1905808", "GO:1905810", "GO:1905811", 
             "GO:0016322", "GO:1904799", "GO:1904800", "GO:1904801", "GO:0006958", "GO:0006956", "GO:0004875", "GO:0002430", 
             "GO:0030449", "GO:0008066", "GO:0050808", "GO:0045088", "GO:0050766", "GO:0019955")

# Ast
bar_plot <- funcenrichplot(
  ORA_ast@result %>% filter(ID %in% int_ont) %>% mutate(Celltype = "Astrocyte"),
  p_col = "p.adjust", 
  highlight = "Celltype", 
  color_list = c("Astrocyte" = paste0(brewer.pal(7, "Set3")[7],"aa"))
)
ggsave("results/ORA_ast.pdf", plot = bar_plot, height = 3, width = 5)



session_info()
