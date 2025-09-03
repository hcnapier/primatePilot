# 07 Compare MLI Subtypes
# Compare genes in MLI subtypes across species 
# Hailey Napier
# August 28, 2025

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

# 1.0 Subset & Re-cluster MLI1s ----
## 1.1 Mouse ----
mouseMLI1s <- noGarbage[['mouse']] %>% subset(idents = c("MLI1"))
mouseMLI1s <- mouseMLI1s %>%
  RunPCA(verbose = F) %>%
  RunUMAP(dims = 1:30, verbose = F) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters()
DimPlot(mouseMLI1s)

## 1.2 Rhesus ----
rhesusMLI1s <- noGarbage[['rhesus']] %>% subset(idents = "MLI1")
rhesusMLI1s <- rhesusMLI1s %>%
  RunPCA(verbose = F) %>%
  RunUMAP(dims = 1:30, verbose = F) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters()
DimPlot(rhesusMLI1s, pt.size = 2)

## 1.3 Human 1 ---- 
human1MLI1s <- noGarbage[['human1']] %>% subset(idents = "MLI1")
human1MLI1s <- human1MLI1s %>%
  RunPCA(verbose = F) %>%
  RunUMAP(dims = 1:30, verbose = F) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters()
DimPlot(human1MLI1s, pt.size = 2)
