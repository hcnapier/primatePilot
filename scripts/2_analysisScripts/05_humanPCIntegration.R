# Test human PC integration

# 0.0 Setup ----
## 0.1 Load packages ----
library(rliger)
library(Seurat)
library(ggplot2)
library(dplyr)
## 0.2 Load data ----
#setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/seuratObjs")
#noGarbage <- readRDS("speciesObjList_cleanCellTypeLabled.rds")

## 0.3 Data setup ----
human1PCs <- subset(noGarbage[['human1']], idents = "Purkinje")
human2PCs <- subset(noGarbage[['human2']], idents = "Purkinje")
human1PCs@meta.data$sample <- "human1"
human2PCs@meta.data$sample <- "human2"
humanPCs_merged <- merge(human1PCs, human2PCs)
humanPCs_merged
humanPCs_merged[["SCT"]] <- split(humanPCs_merged[["SCT"]], f = humanPCs_merged$sample)
DefaultAssay(humanPCs_merged) <- "RNA"
humanPCs_merged

# 1.0 Liger integration ----
## 1.1 Preprocessing ----
humanPCs_merged <- humanPCs_merged %>%
  normalize() %>%
  selectGenes() %>%
  scaleNotCenter()
humanPCs_merged

## 1.2 Integration ----
humanPCs_merged <- humanPCs_merged %>%
  runINMF(k = 17) %>%
  quantileNorm()
humanPCs_merged

## 1.2 Plot integration result ----
humanPCs_merged <- RunUMAP(humanPCs_merged, reduction = "inmfNorm", dims = 1:17)
plotBySample <- DimPlot(humanPCs_merged, group.by = "sample")
plotCellType <- DimPlot(humanPCs_merged, group.by = "seurat_clusters")
plotBySample + plotCellType

## 1.3 Re-cluster ----
humanPCs_merged <- humanPCs_merged %>%
  FindNeighbors(reduction = "inmfNorm", dims = 1:17) %>%
  FindClusters()
DimPlot(humanPCs_merged, group.by = "seurat_clusters", pt.size = 4)
DimPlot(humanPCs_merged, group.by = "sample")


# 2.0 Save objects ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/seuratObjs")
saveRDS(humanPCs_merged, "humanPCsMerged.rds")
