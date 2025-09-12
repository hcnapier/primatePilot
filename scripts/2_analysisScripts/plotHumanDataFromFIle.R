setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/seuratObjs")
noGarbage <- readRDS("speciesObjList_cleanCellTypeLabled.rds")
humanPCs_merged <- readRDS("humanPCsMerged.rds")
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data")
mousePCMarkers <- read.csv("mousePCSubtypeMarkers.csv", row.names = 1)
humanPCMarkers <- read.csv("humanPCSubtypeMarkers.csv", row.names = 1)
setwd("/Users/haileynapier/Work/VertGenLab/Projects/zebrinEvolution/Data/geneLists")
orthologsAll <- read.delim("human_rhesusMouseOrthologs.txt", sep = "\t", header = T)
mouseOrthologs <- orthologsAll %>%
  select(Gene.name, Mouse.gene.name, Mouse.homology.type) %>%
  filter(Mouse.homology.type == "ortholog_one2one")
names(mouseOrthologs) <- c("humanGeneName", "geneName", "homologyType")
mousePCMarkers <- left_join(mousePCMarkers, mouseOrthologs, by = "geneName")
humanPCMarkersFilt <- humanPCMarkers %>%
  filter(p_val_adj <= 0.05) 
mousePCMarkersFilt <- mousePCMarkers %>%
  filter(p_val_adj <= 0.05)
sampleNames

# All clusters ----
Idents(noGarbage[["human2"]]) <- "seurat_clusters"
newIDs_human2 <- c("Granule Neuron", 
                   "MLI Type 2", 
                   "MLI Type 1", 
                   "Granule Neuron", 
                   "MLI Type 1", 
                   "Bergmann Glia", 
                   "Golgi Cell", 
                   "MLI Type 2", 
                   "MLI Type 1", 
                   "UBC",
                   "Oligodendrocyte",
                   "PLI", 
                   "Purkinje Cell", 
                   "Golgi Cell")
newIDs_human2 <- c("Granule Neuron", 
                   "MLI Type 2", 
                   "MLI Type 1", 
                   "Granule Neuron", 
                   "MLI Type 1", 
                   "Bergmann Glia", 
                   "Golgi Cell", 
                   "MLI Type 2", 
                   "MLI Type 1", 
                   "Unipolar Brush Cell",
                   "Oligodendrocyte",
                   "Purkinje Layer Interneuron", 
                   "Purkinje Cell", 
                   "Golgi Cell")
names(newIDs_human2) <- levels(noGarbage[['human2']])
noGarbage[['human2']] <- RenameIdents(noGarbage[['human2']], newIDs_human2)
allCellsPlot <- DimPlot(noGarbage[['human2']], reduction = 'umap', label = F)
LabelClusters(allCellsPlot, id = "ident",fontface = "bold")
FeaturePlot(noGarbage[['human2']], features = "CALB1")

# Integrated human PCs ----
DefaultAssay(humanPCs_merged) <- "SCT"
#humanPCs_merged <- JoinLayers(humanPCs_merged)
DimPlot(humanPCs_merged, pt.size = 5, cols = c("#79A3BC", "#E2B1AC"))
FeaturePlot(humanPCs_merged, features = c("GRID2", "CA8", "KCTD8"), ncol = 2, pt.size = 3)
VlnPlot(humanPCs_merged, features = c("GRID2", "CA8", "RNR2", "COX2"), ncol = 2, pt.size = 1.5, cols = c("#79A3BC", "#E2B1AC"))

# Human and mouse marker comparison ----
ggVennDiagram(list(mousePCMarkersFilt$humanGeneName, humanPCMarkersFilt$gene), 
              category.names = c("Mouse", "Human"), label = "none") +
  coord_flip()
intersect(mousePCMarkersFilt$humanGeneName, humanPCMarkersFilt$gene)


# Cell type markers ----
top5NoGarbage <- list()
for(currSample in sampleNames){
  noGarbageMarkers[[currSample]] %>%
    group_by(cluster) %>%
    dplyr::filter("av_log2FC" > 1) %>%
    dplyr::slice_head(n = 5) %>%
    ungroup() -> top5NoGarbage[[currSample]]
}
top3NoGarbage <- list()
for(currSample in sampleNames){
  noGarbageMarkers[[currSample]] %>%
    group_by(cluster) %>%
    dplyr::filter("av_log2FC" > 1) %>%
    dplyr::slice_head(n = 3) %>%
    ungroup() -> top3NoGarbage[[currSample]]
}
DimPlot(noGarbage[['human2']])
DotPlot(noGarbage[['human2']], features = unique(top3NoGarbage[['human2']]$gene), dot.scale = 5) + 
  theme(axis.text.x = element_text(size = 7, angle = 45, vjust = 1, hjust = 1))

top10NoGarbage[['human2']]$gene %>% unique()
