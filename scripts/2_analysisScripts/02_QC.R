# 02_QC
# QC for primate pilot sequencing analysis
# Hailey Napier
# August 19, 2025

# 0.0 Setup ----
## 0.1 Install Packages ---- 
library(Seurat)
library(ggplot2)
library(dplyr)

## 0.2 Load data
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/seuratObjs")
speciesList <- readRDS("speciesObjList.rds")

# 1.0 GEX QC ----
DefaultAssay(speciesList[['human']]) <- "RNA"
DefaultAssay(speciesList[['rhesus']]) <- "RNA"
DefaultAssay(speciesList[['mouse']]) <- "RNA"

## 1.2 Visualize QC metrics by HTO ----
Idents(speciesList[["human"]]) <- "MULTI_ID"
VlnPlot(speciesList[['human']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = 0.1)
Idents(speciesList[["mouse"]]) <- "MULTI_ID"
VlnPlot(speciesList[['mouse']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = 0.1)
Idents(speciesList[["rhesus"]]) <- "MULTI_ID"
VlnPlot(speciesList[['rhesus']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = 0.1)

## 1.3 Visualize QC metrics ----
Idents(speciesList[['human']]) <- speciesList[['human']]$orig.ident
Idents(speciesList[['rhesus']]) <- speciesList[['rhesus']]$orig.ident
Idents(speciesList[['mouse']]) <- speciesList[['mouse']]$orig.ident
VlnPlot(speciesList[['human']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = 0.1)
VlnPlot(speciesList[['mouse']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = 0.1)
VlnPlot(speciesList[['rhesus']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = 0.1)
FeatureScatter(speciesList[['human']], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(speciesList[['mouse']], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(speciesList[['rhesus']], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

## 1.4 Filter out the cells with low numbers of unique features ----
speciesList_filt <- list()
for(currSpecies in speciesNames){
  speciesList_filt[[currSpecies]] <- subset(speciesList[[currSpecies]], subset = nFeature_RNA > 200 & nFeature_RNA < 8000)
}
### Plot filtered objects ----
VlnPlot(speciesList_filt[['human']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = 0.1)
VlnPlot(speciesList_filt[['mouse']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = 0.1)
VlnPlot(speciesList_filt[['rhesus']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = 0.1)
## HTO
Idents(speciesList[["human"]]) <- "MULTI_ID"
VlnPlot(speciesList[['human']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = 0.1)
Idents(speciesList[["mouse"]]) <- "MULTI_ID"
VlnPlot(speciesList[['mouse']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = 0.1)
Idents(speciesList[["rhesus"]]) <- "MULTI_ID"
VlnPlot(speciesList[['rhesus']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = 0.1)

### How many cells remain?
ncol(speciesList_filt[['human']])
ncol(speciesList_filt[['mouse']])
ncol(speciesList_filt[['rhesus']])

### Replace original seurat object with filtered one -----
speciesList <- speciesList_filt
rm(speciesList_filt)

# 2.0 Save object ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/seuratObjs")
saveRDS(speciesList, "speciesObjList_qcFilt.rds")
