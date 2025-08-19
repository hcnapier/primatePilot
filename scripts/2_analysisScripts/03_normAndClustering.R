# 03_normAndClustering
# Normalization, clustering, and dimensional reduction of primate pilot sequencing data
# Hailey Napier
# August 19, 2025

# 0.0 Setup ----
## 0.1 Load packages ----
library(Seurat)
library(ggplot2)
library(glmGamPoi)

## 0.2 Load data ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/seuratObjs")
speciesList <- readRDS("speciesObjList_qcFilt.rds")


# 1.0 Apply sctransform normalization ----
## This takes ~5-10 minutes and eats up ram
speciesList <- lapply(speciesList, PercentageFeatureSet, pattern = "^MT-", col.name = "percent.mt")
options(future.globals.maxSize = 8000 * 1024^2)
speciesList <- lapply(speciesList, SCTransform, vars.to.regress = "percent.mt", verbose = F)

# 2.0 Dimensional reduction via PCA ----
speciesList <- lapply(speciesList, RunPCA, verbose = FALSE)
## Elbow/Scree plot ----
ElbowPlot(speciesList[['human']])
ElbowPlot(speciesList[['rhesus']])
ElbowPlot(speciesList[['mouse']])
## Heatmap ----
DimHeatmap(speciesList[['human']], dims = 1:20, cells = 500, balanced = T)
DimHeatmap(speciesList[['rhesus']], dims = 1:20, cells = 500, balanced = T)
DimHeatmap(speciesList[['mouse']], dims = 1:20, cells = 500, balanced = T)

# 3.0 UMAP & Clustering ----
## Note that the Seurat authors recommend using more PCs for UMAP reduction after SCTransform
speciesList <- lapply(speciesList, RunUMAP, dims = 1:30, verbose = F)
speciesList <- lapply(speciesList, FindNeighbors, dims = 1:30, verbose = F)
speciesList <- lapply(speciesList, FindClusters, verbose = F)

# 4.0 Plot ----
## 4.1 DimPlot ----
DimPlot(speciesList[['human']], label = T)
DimPlot(speciesList[['mouse']], label = T)
DimPlot(speciesList[['rhesus']], label = T)

## 4.2 DimPlot by HTO ----
Idents(speciesList[["human"]]) <- "MULTI_ID"
Idents(speciesList[["mouse"]]) <- "MULTI_ID"
Idents(speciesList[["rhesus"]]) <- "MULTI_ID"
DimPlot(speciesList[['human']], label = F)
DimPlot(speciesList[['mouse']], label = F)
DimPlot(speciesList[['rhesus']], label = F)
