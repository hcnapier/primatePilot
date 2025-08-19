# 01_htoDemux
# Demultiplex HTOs
# Hailey Napier
# August 4, 2025

# 0.0 Setup ----
## 0.1 Load packages -----
library(Seurat)
library(patchwork)
library(dplyr)
library(deMULTIplex)
library(ggplot2)

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


# 2.0 Test methods ----
## 2.1 Try HTODemux 
testObj <- speciesList[["mouse"]]
testHTO <- HTODemux(testObj, positive.quantile = 0.99)
table(testHTO$HTO_classification)
table(testHTO$HTO_classification.global)
#Idents(testObj) <- "HTO_maxID"
#RidgePlot(testObj, assay = "HTO", features = rownames(testObj[["HTO"]])[1:2], ncol = 2)

## 2.2 Try MULTISeqDemux
testMultiSeq <- MULTIseqDemux(speciesList[["mouse"]], autoThresh = T, quantile = .14)
table(testMultiSeq$MULTI_ID)
plot(testMultiSeq$MULTI_ID)

## 2.3 Try deMULTIplex
## Can't really figure out how to use this...
testObj <- speciesList[["mouse"]]
bar.tsne <- barTSNE(testObj[['HTO']]$counts[,1:8])


# 3.0 Demultiplex all with MULTISeqDemux ----
speciesList <- lapply(speciesList, MULTIseqDemux, autoThresh = T)


# 4.0 Plot HTO demultiplexing results ----
## 4.1 Ridge plots ----
### Mouse ----
Idents(speciesList[["mouse"]]) <- "MULTI_ID"
RidgePlot(speciesList[["mouse"]], assay = "HTO", features = rownames(speciesList[["mouse"]][["HTO"]]), ncol = 4)
### Rhesus ----
Idents(speciesList[["rhesus"]]) <- "MULTI_ID"
RidgePlot(speciesList[["rhesus"]], assay = "HTO", features = rownames(speciesList[["rhesus"]][["HTO"]]), ncol = 4)
### Human ----
Idents(speciesList[["human"]]) <- "MULTI_ID"
RidgePlot(speciesList[["human"]], assay = "HTO", features = rownames(speciesList[["human"]][["HTO"]]), ncol = 4)

## 4.2 NCountRNA Violin Plots ----
### Mouse ----
VlnPlot(speciesList[["mouse"]], features = "nCount_RNA", pt.size = 0.1, log = T) + 
  labs(title = "Mouse\nnCount RNA")
### Rhesus ----
VlnPlot(speciesList[["rhesus"]], features = "nCount_RNA", pt.size = 0.1, log = T) + 
  labs(title = "Rhesus\nnCount RNA")
### Human ----
VlnPlot(speciesList[["human"]], features = "nCount_RNA", pt.size = 0.1, log = T) + 
  labs(title = "Human\nnCount RNA")

## 4.3 tSNE Embedding ----
### Mouse ----
# Remove negative cells
tmp_subset_mouse <- subset(speciesList[["mouse"]], idents = c("Negative", "Doublet"), invert = TRUE)
# Calculate a tSNE embedding of the HTO data
DefaultAssay(tmp_subset_mouse) <- "HTO"
tmp_subset_mouse <- ScaleData(tmp_subset_mouse, features = rownames(tmp_subset_mouse),
                                 verbose = FALSE)
tmp_subset_mouse <- RunPCA(tmp_subset_mouse, features = rownames(tmp_subset_mouse), approx = FALSE)
tmp_subset_mouse <- RunTSNE(tmp_subset_mouse, dims = 1:8, perplexity = 100)
DimPlot(tmp_subset_mouse) + 
  labs(title = "Mouse")
### Rhesus ----
# Remove negative cells
tmp_subset_rhesus <- subset(speciesList[["rhesus"]], idents = c("Negative", "Doublet"), invert = TRUE)
# Calculate a tSNE embedding of the HTO data
DefaultAssay(tmp_subset_rhesus) <- "HTO"
tmp_subset_rhesus <- ScaleData(tmp_subset_rhesus, features = rownames(tmp_subset_rhesus),
                              verbose = FALSE)
tmp_subset_rhesus <- RunPCA(tmp_subset_rhesus, features = rownames(tmp_subset_rhesus), approx = FALSE)
tmp_subset_rhesus <- RunTSNE(tmp_subset_rhesus, dims = 1:8, perplexity = 100)
DimPlot(tmp_subset_rhesus) + 
  labs(title = "Rhesus")
### Human ----
# Remove negative cells
tmp_subset_human <- subset(speciesList[["human"]], idents = c("Negative", "Doublet"), invert = TRUE)
# Calculate a tSNE embedding of the HTO data
DefaultAssay(tmp_subset_human) <- "HTO"
tmp_subset_human <- ScaleData(tmp_subset_human, features = rownames(tmp_subset_human),
                               verbose = FALSE)
tmp_subset_human <- RunPCA(tmp_subset_human, features = rownames(tmp_subset_human), approx = FALSE)
tmp_subset_human <- RunTSNE(tmp_subset_human, dims = 1:8, perplexity = 100)
DimPlot(tmp_subset_human) + 
  labs(title = "Human")

## 4.5 tSNE Embedding, only species-specific tags ----
### Mouse ----
# Remove negative cells
tmp_subset_mouse <- subset(speciesList[["mouse"]], 
                            idents = c("Negative", 
                                       "Doublet",
                                       "HTO-5",
                                       "HTO-6",
                                       "HTO-7", 
                                       "HTO-8"), 
                            invert = TRUE)
# Calculate a tSNE embedding of the HTO data
DefaultAssay(tmp_subset_mouse) <- "HTO"
tmp_subset_mouse <- ScaleData(tmp_subset_mouse, features = rownames(tmp_subset_mouse),
                              verbose = FALSE)
tmp_subset_mouse <- RunPCA(tmp_subset_mouse, features = rownames(tmp_subset_mouse), approx = FALSE)
tmp_subset_mouse <- RunTSNE(tmp_subset_mouse, dims = 1:8, perplexity = 100)
DimPlot(tmp_subset_mouse) + 
  labs(title = "Mouse")
### Rhesus ----
# Remove negative cells
tmp_subset_rhesus <- subset(speciesList[["rhesus"]], 
                            idents = c("Negative", 
                                       "Doublet",
                                       "HTO-1",
                                       "HTO-2",
                                       "HTO-3", 
                                       "HTO-4"), 
                            invert = TRUE)
# Calculate a tSNE embedding of the HTO data
DefaultAssay(tmp_subset_rhesus) <- "HTO"
tmp_subset_rhesus <- ScaleData(tmp_subset_rhesus, features = rownames(tmp_subset_rhesus),
                               verbose = FALSE)
tmp_subset_rhesus <- RunPCA(tmp_subset_rhesus, features = rownames(tmp_subset_rhesus), approx = FALSE)
tmp_subset_rhesus <- RunTSNE(tmp_subset_rhesus, dims = 1:8, perplexity = 100)
DimPlot(tmp_subset_rhesus) + 
  labs(title = "Rhesus")


# 5.0 Save objects ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/seuratObjs")
saveRDS(speciesList, "speciesObjList.R")
