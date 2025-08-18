# 01_preprocessingQC
# Preprocessing and QC of Seurat objects for analysis
# Hailey Napier
# August 4, 2025

# 0.0 Setup ----
## 0.1 Load packages -----
library(Seurat)
library(patchwork)
library(dplyr)
library(deMULTIplex)

## 0.2 Source functions ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/functions")
source("addHTOAssay.R")

## 0.3 Read data ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/countMats")
speciesCountMatList_gex <- readRDS("allSpeciesCountMatList_GEX.rds")
speciesCountMatList_hto <- readRDS("allSpeciesCountMatList_HTO.rds")
speciesNames <- c("human", "rhesus", "mouse")

# 1.0 Setup Seurat Object ----
speciesList <- list()
for(currSpecies in speciesNames){
  speciesList[[currSpecies]] <- CreateSeuratObject(counts = speciesCountMatList_gex[[currSpecies]]) # set up GEX Seurat object
  speciesList[[currSpecies]] <- addHTOAssay(obj = speciesList[[currSpecies]], htoMat = speciesCountMatList_hto[[currSpecies]])
}
rm(speciesCountMatList_gex)
rm(speciesCountMatList_hto)


# 2.0 Demultiplex HTOs ----
## 2.1 Try HTODemux ---- 
testObj <- speciesList[["mouse"]]
testHTO <- HTODemux(testObj, positive.quantile = 0.99)
table(testHTO$HTO_classification)
table(testHTO$HTO_classification.global)
#Idents(testObj) <- "HTO_maxID"
#RidgePlot(testObj, assay = "HTO", features = rownames(testObj[["HTO"]])[1:2], ncol = 2)

## 2.2 Try MULTISeqDemux ----
testMultiSeq <- MULTIseqDemux(speciesList[["mouse"]], autoThresh = T, quantile = .14)
table(testMultiSeq$MULTI_ID)
plot(testMultiSeq$MULTI_ID)

## 2.3 Try deMULTIplex
testObj <- speciesList[["mouse"]]
bar.tsne <- barTSNE(testObj[['HTO']]$counts[,1:8])


# 2.0 GEX QC ----
## 2.1 Visualize QC metrics ----
VlnPlot(speciesList[['human']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(speciesList[['mouse']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(speciesList[['rhesus']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(speciesList[['human']], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(speciesList[['mouse']], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(speciesList[['rhesus']], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

## 2.2 Filter out the cells with low numbers of unique features ----
speciesList_filt <- list()
for(currSpecies in speciesNames){
  speciesList_filt[[currSpecies]] <- subset(speciesList[[currSpecies]], subset = nFeature_RNA > 750 & nFeature_RNA < 9000)
}
### Plot filtered objects ----
VlnPlot(speciesList_filt[['human']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(speciesList_filt[['mouse']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(speciesList_filt[['rhesus']], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
### Replace original seurat object with filtered one -----
speciesList <- speciesList_filt
rm(speciesList_filt)


# 3.0 Filter HTO matrix by GEX QC filtering ----
for(currSpecies in speciesNames){
  filtCellList <- Cells(speciesList[[currSpecies]])
  speciesCountMatList_hto[[currSpecies]] <- speciesCountMatList_hto[[currSpecies]][ ,filtCellList]
}


# 4.0 Find and scale variable GEX features ----


# 3.0 Setup & Normalize Seurat Object ----
speciesList <- list()
for(currSpecies in speciesNames){
  print(currSpecies)
  speciesList[[currSpecies]] <- setupSeuratObj_htogex(htoMat = speciesCountMatList_hto[[currSpecies]], gexMat = speciesCountMatList_gex[[currSpecies]])
}
rm(speciesCountMatList_gex)
rm(speciesCountMatList_hto