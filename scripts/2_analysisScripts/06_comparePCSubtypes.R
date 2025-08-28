# 06 Compare PC Subtypes
# Compare genes in PC subtypes across species 
# Hailey Napier
# August 27, 2025

# 0.0 Setup ----
## 0.1 Load packages ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)

## 0.2 Load data ----
#setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/seuratObjs")
#humanPCs_merged <- readRDS("humanPCsMerged.rds")
#noGarbage <- readRDS("speciesObjList_cleanCellTypeLabled.rds")
setwd("/Users/haileynapier/Work/VertGenLab/Projects/zebrinEvolution/Data/geneLists")
orthologsAll <- read.delim("human_rhesusMouseOrthologs.txt", sep = "\t", header = T)
mouseOrthologs <- orthologsAll %>%
  select(Gene.name, Mouse.gene.name, Mouse.homology.type) %>%
  filter(Mouse.homology.type == "ortholog_one2one")
names(mouseOrthologs) <- c("humanGeneName", "geneName", "homologyType")

rhesusOrthologs <- orthologsAll %>%
  select(Gene.name, Macaque.gene.name, Macaque.homology.type) %>%
  filter(Macaque.homology.type == "ortholog_one2one")
names(rhesusOrthologs) <- c("humanGeneName", "geneName", "homologyType")

rm(orthologsAll)

# 1.0 Subset & Re-cluster PCs ----
## 1.1 Mouse ----
mousePCs <- noGarbage[['mouse']] %>% subset(idents = "Purkinje")
mousePCs <- mousePCs %>%
  RunPCA(verbose = F) %>%
  RunUMAP(dims = 1:30, verbose = F) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters()
DimPlot(mousePCs)
  
## 1.2 Rhesus ----
rhesusPCs <- noGarbage[['rhesus']] %>% subset(idents = "Purkinje")
rhesusPCs <- rhesusPCs %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters()
DimPlot(rhesusPCs, pt.size = 2)

## 1.3 Human ---- 
# already re-clustered
DefaultAssay(humanPCs_merged) <- "SCT"
humanPCs_merged <- JoinLayers(humanPCs_merged)
DimPlot(humanPCs_merged, pt.size = 2.5, cols = c("#79A3BC", "#E2B1AC"))


# 2.0 Find subtype markers ----
mousePCMarkers <- mousePCs %>%
  FindMarkers(ident.1 = "3", ident.2 = "2") %>%
  filter(p_val_adj <= 0.1)
mousePCMarkers$geneName <- rownames(mousePCMarkers)

rhesusPCMarkers <- rhesusPCs %>%
  FindMarkers(ident.1 = "0", ident.2 = "1")  %>%
  filter(p_val_adj <= 0.1)
rhesusPCMarkers$geneName <- rownames(rhesusPCMarkers)

humanPCMarkers <- humanPCs_merged %>%
  FindMarkers(ident.1 = "0", ident.2 = "1")  %>%
  filter(p_val_adj <= 0.1)
humanPCMarkers$humanGeneName <- rownames(humanPCMarkers)

# 3.0 Convert subtype markers to human ortholog ----
mousePCMarkers <- left_join(mousePCMarkers, mouseOrthologs, by = "geneName")
rhesusPCMarkers <- left_join(rhesusPCMarkers, rhesusOrthologs, by = "geneName")

# 4.0 Compare marker genes ----
primateSharedMarkers <- intersect(humanPCMarkers$humanGeneName, rhesusPCMarkers$humanGeneName)
allSharedMarkers <- intersect(primateSharedMarkers, mousePCMarkers$humanGeneName)
allSharedMarkers
rhesusMouseSharedMarkers <- intersect(rhesusPCMarkers$humanGeneName, mousePCMarkers$humanGeneName)
humanMouseSharedMarkers <- intersect(humanPCMarkers$humanGeneName, mousePCMarkers$humanGeneName)

# 5.0 Plot ----
## 5.1 Venn diagrams ----
ggVennDiagram(list(mousePCMarkers$humanGeneName, humanPCMarkers$humanGeneName), 
              category.names = c("Mouse", "Human")) +
  coord_flip()
humanMouseSharedMarkers

ggVennDiagram(list(rhesusPCMarkers$humanGeneName, humanPCMarkers$humanGeneName), 
              category.names = c("Rhesus", "Human")) +
  coord_flip()
primateSharedMarkers

## 5.2 Feature plots ----
FeaturePlot(humanPCs_merged, features = "GRID2", pt.size = 3)

## 5.3 Violin plots ----
VlnPlot(humanPCs_merged, features = c("RNR2", "COX3", "GRID2", "CA8"), ncol = 2, cols = c("#79A3BC", "#E2B1AC"))
        