# 04_cellTypeClassification
# Call cell types in primate pilot sequencing data
# Hailey Napier
# August 20, 2025

# 0.0 Setup ----
## 0.1 Load packages ----
library(Seurat)
library(dplyr)
library(ggplot2)

## 0.2 Load data ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/seuratObjs")
speciesList <- readRDS("speciesObjList_normAndClust.rds")
speciesNames <- c("human", "rhesus", "mouse")


# 1.0 Subset to include true singlets ----
## 1.1 Remove doublets and negative cells ----
singletSpeciesList <- list()
for(currSpecies in speciesNames){
  singletSpeciesList[[currSpecies]] <- subset(speciesList[[currSpecies]], idents = c("Negative", "Doublet"), invert = T)
}
## Re-cluster ----
singletSpeciesList <- lapply(singletSpeciesList, RunPCA, verbose = FALSE)
singletSpeciesList <- lapply(singletSpeciesList, RunUMAP, dims = 1:30, verbose = F)
singletSpeciesList <- lapply(singletSpeciesList, FindNeighbors, dims = 1:30, verbose = F)
singletSpeciesList <- lapply(singletSpeciesList, FindClusters, verbose = F)
### UMAP Plot ----
Idents(singletSpeciesList[["human"]]) <- "MULTI_ID"
Idents(singletSpeciesList[["rhesus"]]) <- "MULTI_ID"
Idents(singletSpeciesList[["mouse"]]) <- "MULTI_ID"
DimPlot(singletSpeciesList[['human']])
DimPlot(singletSpeciesList[['rhesus']])
DimPlot(singletSpeciesList[['mouse']])
### Feature plot by read count ----
FeaturePlot(speciesSpecList[['human1']], features = c("nCount_RNA", "nFeature_RNA"))
FeaturePlot(speciesSpecList[['human2']], features = c("nCount_RNA", "nFeature_RNA"))
FeaturePlot(singletSpeciesList[['rhesus']], features = c("nCount_RNA", "nFeature_RNA"))
FeaturePlot(singletSpeciesList[['mouse']], features = c("nCount_RNA", "nFeature_RNA"))

## 1.2 Subset mouse and rhesus samples to get rid of contaminating cells by HTO ----
### Subset ----
speciesSpecList <- list()
speciesSpecList[['mouse']] <- subset(singletSpeciesList[["mouse"]], 
                                     idents = c("HTO-5", "HTO-6", "HTO-7", "HTO-8"), 
                                     invert = T)
speciesSpecList[['rhesus']] <- subset(singletSpeciesList[["rhesus"]], 
                                      idents = c("HTO-1", "HTO-2","HTO-3", "HTO-4"), 
                                      invert = T)
### Also subset human by individual ----
speciesSpecList[['human1']] <- subset(singletSpeciesList[["human"]], 
                                      idents = c("HTO-1", "HTO-2","HTO-3", "HTO-4"), 
                                      invert = F)
speciesSpecList[['human2']] <- subset(singletSpeciesList[["human"]], 
                                     idents = c("HTO-5","HTO-6","HTO-7","HTO-8"), 
                                     invert = F)

### Re-cluster ----
speciesSpecList <- lapply(speciesSpecList, RunPCA, verbose = FALSE)
speciesSpecList <- lapply(speciesSpecList, RunUMAP, dims = 1:30, verbose = F)
speciesSpecList <- lapply(speciesSpecList, FindNeighbors, dims = 1:30, verbose = F)
speciesSpecList <- lapply(speciesSpecList, FindClusters, verbose = F)
### Plot ----
Idents(speciesSpecList[["human1"]]) <- "MULTI_ID"
Idents(speciesSpecList[["human2"]]) <- "MULTI_ID"
Idents(speciesSpecList[["rhesus"]]) <- "MULTI_ID"
Idents(speciesSpecList[["mouse"]]) <- "MULTI_ID"
DimPlot(speciesSpecList[['human1']])
DimPlot(speciesSpecList[['human2']])
ncol(speciesSpecList[['human1']])
ncol(speciesSpecList[['human2']])
DimPlot(speciesSpecList[['rhesus']])
ncol(speciesSpecList[['rhesus']])
DimPlot(speciesSpecList[['mouse']])
ncol(speciesSpecList[['mouse']])


# 2.0 Plot cell-type markers ----
Idents(speciesSpecList[["human1"]]) <- "seurat_clusters"
Idents(speciesSpecList[["human2"]]) <- "seurat_clusters"
Idents(speciesSpecList[["rhesus"]]) <- "seurat_clusters"
Idents(speciesSpecList[["mouse"]]) <- "seurat_clusters"

## Human 1 ----
### PC Markers ----
DimPlot(speciesSpecList[["human1"]])
FeaturePlot(speciesSpecList[['human1']], features = c("PCP2", "CALB1"))
VlnPlot(speciesSpecList[['human1']], features = c("PCP2", "CALB1"))
table(Idents(speciesSpecList[['human1']]))
### MLI Markers ----
FeaturePlot(speciesSpecList[['human1']], features = c("GAD1"))
### Granule Neuron Markers ----

## Human 2 ----
### PC Markers ----
DimPlot(speciesSpecList[["human2"]])
FeaturePlot(speciesSpecList[['human2']], features = c("PCP2", "CALB1"))
VlnPlot(speciesSpecList[['human2']], features = c("PCP2", "CALB1"))
### MLI Markers ----
FeaturePlot(speciesSpecList[['human2']], features = c("GAD1"))
### Granule Neuron Markers ----

## Rhesus ----
### PC Markers ----
DimPlot(speciesSpecList[['rhesus']])
FeaturePlot(speciesSpecList[['rhesus']], features = c("PCP2", "CALB1"))
VlnPlot(speciesSpecList[['rhesus']], features = c("PCP2", "CALB1"))
table(Idents(speciesSpecList[['rhesus']]))

## Mouse ----
### PC Markers ----
DimPlot(speciesSpecList[['mouse']])
FeaturePlot(speciesSpecList[['mouse']], features = c("Pcp2", "Calb1"))
VlnPlot(speciesSpecList[['mouse']], features = c("Pcp2", "Calb1"))
table(Idents(speciesSpecList[['mouse']]))


# 4.0 Plot erroneous HTO cell-type markers ----
Idents(singletSpeciesList[["human"]]) <- "seurat_clusters"
Idents(singletSpeciesList[["rhesus"]]) <- "seurat_clusters"
Idents(singletSpeciesList[["mouse"]]) <- "seurat_clusters"
## Rhesus ----
DimPlot(singletSpeciesList[['rhesus']])
FeaturePlot(singletSpeciesList[['rhesus']], features = c("PCP2", "CALB1"))
VlnPlot(singletSpeciesList[['rhesus']], features = c("PCP2", "CALB1"))
table(Idents(speciesSpecList[['rhesus']]))

## Mouse ----


# 5.0 Plot negative cells with species specific HTOs ----
## 5.1 Subset and re-cluster ----
negSingList <- list()
negSingList[['human1']] <- subset(speciesList[['human']], idents = c("Negative", 
                                                                            "HTO-1",
                                                                            "HTO-2",
                                                                            "HTO-3", 
                                                                            "HTO-4"), invert = F)
negSingList[['human2']] <- subset(speciesList[['human']], idents = c("Negative", 
                                                                            "HTO-5",
                                                                            "HTO-6",
                                                                            "HTO-7", 
                                                                            "HTO-8"), invert = F)
negSingList[['mouse']] <- subset(speciesList[['mouse']], idents = c("Negative", 
                                                                            "HTO-1",
                                                                            "HTO-2",
                                                                            "HTO-3", 
                                                                            "HTO-4"), invert = F)
negSingList[['rhesus']] <- subset(speciesList[['rhesus']], idents = c("Negative", 
                                                                            "HTO-5",
                                                                            "HTO-6",
                                                                            "HTO-7", 
                                                                            "HTO-8"), invert = F)
negSingList <- lapply(negSingList, RunPCA, verbose = F)
negSingList <- lapply(negSingList, RunUMAP, dims = 1:30, verbose = F)
negSingList <- lapply(negSingList, FindNeighbors, dims = 1:30, verbose = F)
negSingList <- lapply(negSingList, FindClusters, verbose = F)

## 5.2 DimPlot ----
Idents(negSingList[["human1"]]) <- "MULTI_ID"
Idents(negSingList[["human2"]]) <- "MULTI_ID"
Idents(negSingList[["rhesus"]]) <- "MULTI_ID"
Idents(negSingList[["mouse"]]) <- "MULTI_ID"
DimPlot(negSingList[['human1']])
DimPlot(negSingList[['human2']])
DimPlot(negSingList[['rhesus']])
DimPlot(negSingList[['mouse']])

## 5.3 Plot cell-type markers ----
Idents(negSingList[["human1"]]) <- "seurat_clusters"
Idents(negSingList[["human2"]]) <- "seurat_clusters"
Idents(negSingList[["rhesus"]]) <- "seurat_clusters"
Idents(negSingList[["mouse"]]) <- "seurat_clusters"

### Human 1 ----
DimPlot(negSingList[['human1']])
table(Idents(negSingList[["human1"]]))
#### PCs ----
FeaturePlot(negSingList[['human1']], features = c("PCP2", "CALB1", "FOXP2"))
VlnPlot(negSingList[['human1']], features = c("PCP2", "CALB1", "FOXP2")) # 18 putative PCs (cluster 15)

### Human 2 ----
DimPlot(negSingList[['human2']])
table(Idents(negSingList[["human2"]]))
#### PCs ----
FeaturePlot(negSingList[['human2']], features = c("PCP2", "CALB1", "FOXP2"))
VlnPlot(negSingList[['human2']], features = c("PCP2", "CALB1", "FOXP2")) # 70 putative PCs (cluster 14)
#### Astrocytes ----

### Rhesus ----
DimPlot(negSingList[['rhesus']])
table(Idents(negSingList[["rhesus"]]))
#### PCs ----
FeaturePlot(negSingList[['rhesus']], features = c("PCP2", "CALB1", "FOXP2"))
VlnPlot(negSingList[['rhesus']], features = c("PCP2", "CALB1", "FOXP2")) # 45 putative PCS (cluster 12)

### Mouse ----
DimPlot(negSingList[['mouse']])
table(Idents(negSingList[['mouse']]))
#### PCs ----
FeaturePlot(negSingList[['mouse']], features = c("Pcp2", "Calb1", "Foxp2"))
VlnPlot(negSingList[['mouse']], features = c("Pcp2", "Calb1", "Foxp2")) # 392 putative PCs (clusters 3&7)

## 5.4 Subset putative Purkinje cells ----
### Subset & Re-cluster ----
putativePCs <- list()
putativePCs[['mouse']] <- subset(negSingList[["mouse"]], 
                                 idents = c("3","7"), 
                                 invert = F)
# putativePCs[['rhesus']] <- subset(negSingList[['rhesus']], 
#                                   idents = c("12"), 
#                                   invert = F)
# putativePCs[['human1']] <- subset(negSingList[['human1']], 
#                                   idents = c("15"), 
#                                   invert = F)
putativePCs[['human2']] <- subset(negSingList[['human2']], 
                                  idents = c("14"), 
                                  invert = F)
#putativePCs[['human1']] <- RunPCA(putativePCs[["human1"]], npcs = 17, maxit = 10000, verbose = F)
putativePCs[['human2']] <- RunPCA(putativePCs[["human2"]], verbose = F)
#putativePCs[['rhesus']] <- RunPCA(putativePCs[["rhesus"]], npcs = 44, maxit = 10000, verbose = F)
putativePCs[['mouse']] <- RunPCA(putativePCs[["mouse"]], verbose = F)
putativePCs <- lapply(putativePCs, RunUMAP, dims = 1:30, verbose = F)
putativePCs <- lapply(putativePCs, FindNeighbors, dims = 1:30, verbose = F)
putativePCs <- lapply(putativePCs, FindClusters, verbose = F)
## Plot ----
### Rhesus ----
#FeaturePlot(putativePCs[['rhesus']], features = c("ALDOC", "MPPED2"), ncol = 1) # ZII+ Markers ----
#FeaturePlot(putativePCs[['rhesus']], features = c("EBF2", "PLCB4"), ncol = 1) # ZII- Markers ----
#FeaturePlot(putativePCs[['rhesus']], features = c("GRID2", "PRKG1")) # Grid2 Markers ----
### Mouse ----
FeaturePlot(putativePCs[['mouse']], features = c("Aldoc", "Mpped2"), ncol = 1) # ZII+ Markers ----
FeaturePlot(putativePCs[['mouse']], features = c("Ebf2", "Plcb4"), ncol = 1) # ZII- Markers ----
FeaturePlot(putativePCs[['mouse']], features = c("Grid2", "Prkg1")) # Grid2 Markers ----
### Human 2 ----
FeaturePlot(putativePCs[['human2']], features = c("ALDOC", "MPPED2"), ncol = 1, pt.size = 3) # ZII+ Markers ----
FeaturePlot(putativePCs[['human2']], features = c("EBF2", "PLCB4"), ncol = 1, pt.size = 3) # ZII- Markers ----
FeaturePlot(putativePCs[['human2']], features = c("GRID2", "PRKG1"), pt.size = 3) # Grid2 Markers ----


