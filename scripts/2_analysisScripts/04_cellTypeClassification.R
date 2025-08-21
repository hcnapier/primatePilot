# 04_cellTypeClassification
# Call cell types in primate pilot sequencing data
# Hailey Napier
# August 20, 2025

# 0.0 Setup ----
## 0.1 Load packages ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)

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
Idents(negSingList[["human1"]]) <- "seurat_clusters"
DimPlot(negSingList[['human1']])
table(Idents(negSingList[["human1"]]))
#### PCs -- 15
FeaturePlot(negSingList[['human1']], features = c("PCP2", "CALB1"))
VlnPlot(negSingList[['human1']], features = c("PCP2", "CALB1")) # 18 putative PCs (cluster 15)
VlnPlot(negSingList[['human1']], features = c("PCP2", "CALB1", "FOXP2", "EBF1", "CNTNAP4", "ITPR1"))
#### Astrocytes -- 14
FeaturePlot(negSingList[['human1']], features = c("GFAP", "SLC4A4", "GABRB1", "SLC1A2"))
VlnPlot(negSingList[['human1']], features = c("GFAP", "SLC4A4", "GABRB1", "SLC1A2", "LUZP2", "AQP4")) # Cluster 14
#### Bergmann Glia -- 8
FeaturePlot(negSingList[['human1']], features = c("SLC4A4", "GABRB1", "SLC1A2","PREX2"))
VlnPlot(negSingList[['human1']], features = c("SLC4A4", "GABRB1", "SLC1A2","PREX2", "SLC1A3", "PTPRZ1")) 
##### Granule Neurons -- 3 & 0
DimPlot(negSingList[['human1']], label = T)
FeaturePlot(negSingList[['human1']], features = c("CHN2", "TENM1", "CADPS2", "CA10"))
VlnPlot(negSingList[['human1']], features = c("CHN2", "TENM1", "CADPS2", "CA10", "GRIK2"))
FindMarkers(negSingList[['human1']], ident.1 = c("3", "0"), ident.2 = NULL)
#### Golgi Cells -- 7
FeaturePlot(negSingList[['human1']], features = c("GABRG3", "CACNG3", "ELAVL2", "KCNN3"))
VlnPlot(negSingList[['human1']], features = c("GABRG3", "CACNG3", "ELAVL2", "KCNN3", "LGI2"))

### Human 2 ----
DimPlot(negSingList[['human2']])
table(Idents(negSingList[["human2"]]))
#### PCs -- 14
FeaturePlot(negSingList[['human2']], features = c("PCP2", "CALB1"))
VlnPlot(negSingList[['human2']], features = c("PCP2", "CALB1")) # 70 putative PCs (cluster 14)
VlnPlot(negSingList[['human2']], features = c("PCP2", "CALB1", "FOXP2", "EBF1", "CNTNAP4", "ITPR1"))
#### Bergmann Glia -- 5 
FeaturePlot(negSingList[['human2']], features = c("SLC4A4", "GABRB1", "SLC1A2","PREX2"))
VlnPlot(negSingList[['human2']], features = c("SLC4A4", "GABRB1", "SLC1A2","PREX2", "SLC1A3", "PTPRZ1")) 
#### Granule Neurons -- 3 & 0
DimPlot(negSingList[['human2']], label = T)
FeaturePlot(negSingList[['human2']], features = c("CHN2", "TENM1", "CADPS2", "CA10"))
VlnPlot(negSingList[['human2']], features = c("CHN2", "TENM1", "CADPS2", "CA10", "GRIK2"))
#### Golgi Cells -- 7 (minor 15)
FeaturePlot(negSingList[['human2']], features = c("GABRG3", "CACNG3", "ELAVL2", "KCNN3"))
VlnPlot(negSingList[['human2']], features = c("GABRG3", "CACNG3", "ELAVL2", "KCNN3", "LGI2"))

### Rhesus ----
DimPlot(negSingList[['rhesus']])
table(Idents(negSingList[["rhesus"]]))
#### PCs -- 12
FeaturePlot(negSingList[['rhesus']], features = c("PCP2", "CALB1"))
VlnPlot(negSingList[['rhesus']], features = c("PCP2", "CALB1")) # 45 putative PCS (cluster 12)
VlnPlot(negSingList[['rhesus']], features = c("PCP2", "CALB1", "FOXP2", "EBF1", "CNTNAP4", "ITPR1"))
#### Bergmann Glia -- 3
FeaturePlot(negSingList[['rhesus']], features = c("SLC4A4", "GABRB1", "SLC1A2","PREX2"))
VlnPlot(negSingList[['rhesus']], features = c("SLC4A4", "GABRB1", "SLC1A2","PREX2", "SLC1A3", "PTPRZ1")) 
#### Granule Neurons -- 9???
DimPlot(negSingList[['rhesus']], label = T)
FeaturePlot(negSingList[['rhesus']], features = c("CHN2", "TENM1", "CADPS2", "CA10"))
VlnPlot(negSingList[['rhesus']], features = c("CHN2", "TENM1", "CADPS2", "CA10", "GRIK2", "RBFOX3"))
FindMarkers(negSingList[['rhesus']], ident.1 = "9", ident.2 = NULL)
#### Golgi Cells -- 5 (minor 14)
FeaturePlot(negSingList[['rhesus']], features = c("GABRG3", "CACNG3", "ELAVL2", "KCNN3"))
VlnPlot(negSingList[['rhesus']], features = c("GABRG3", "CACNG3", "ELAVL2", "KCNN3", "LGI2"))


### Mouse ----
DimPlot(negSingList[['mouse']])
table(Idents(negSingList[['mouse']]))
#### PCs -- 3&7
FeaturePlot(negSingList[['mouse']], features = c("Pcp2", "Calb1"))
VlnPlot(negSingList[['mouse']], features = c("Pcp2", "Calb1")) # 392 putative PCs (clusters 3&7)
VlnPlot(negSingList[['mouse']], features = c("Pcp2", "Calb1", "Foxp2", "Ebf1", "Cntnap4", "Itpr1"))
#### Astrocytes -- 14
FeaturePlot(negSingList[['mouse']], features = c("Gfap", "Slc4a4", "Gabrb1", "Slc1a2"))
VlnPlot(negSingList[['mouse']], features = c("Gfap", "Slc4a4", "Gabrb1", "Slc1a2", "Luzp2", "Aqp4")) # Cluster 14?
# what distinguishes clusters 4 and 14?
fourvsfourteen <- FindMarkers(negSingList[['mouse']], ident.1 = "4", ident.2 = "14")
ms4 <- FindMarkers(negSingList[['mouse']], ident.1 = "4", ident.2 = NULL)
ms14 <- FindMarkers(negSingList[['mouse']], ident.1 = "14", ident.2 = NULL)
#### Bergmann Glia -- 4
FeaturePlot(negSingList[['mouse']], features = c("Slc4a4", "Gabrb1", "Slc1a2", "Prex2"))
VlnPlot(negSingList[['mouse']], features = c("Slc4a4", "Gabrb1", "Slc1a2", "Prex2", "Slc1a3", "Ptprz1"))
#### Granule Neurons -- 12
FeaturePlot(negSingList[['mouse']], features = c("Chn2", "Tenm1", "Cadps2", "Car10"))
VlnPlot(negSingList[['mouse']], features = c("Chn2", "Tenm1", "Cadps2", "Car10", "Grik2"))
#### Golgi Cells -- 19 & 5
DimPlot(negSingList[['mouse']], label = T)
FeaturePlot(negSingList[['mouse']], features = c("Gabrg3", "Cacng3", "Elavl2", "Kcnn3"))
VlnPlot(negSingList[['mouse']], features = c("Gabrg3", "Cacng3", "Elavl2", "Kcnn3", "Lgi2"))



## 5.4 Compare PC markers ----
h1PCMarkers <- FindMarkers(negSingList[['human1']], ident.1 = "15", ident.2 = NULL)
h2PCMarkers <- FindMarkers(negSingList[['human2']], ident.1 = "14", ident.2 = NULL)
rhesusPCMarkers <- FindMarkers(negSingList[['rhesus']], ident.1 = "12", ident.2 = NULL)
mousePCMarkers <- FindMarkers(negSingList[['mouse']], ident.1 = c("3", "7"), ident.2 = NULL)
### Shared markers 
pcMarkerList <- list()
pcMarkerList[['human1']] <- rownames(h1PCMarkers) %>% tolower()
pcMarkerList[['human2']] <- rownames(h2PCMarkers) %>% tolower()
pcMarkerList[['rhesus']] <- rownames(rhesusPCMarkers) %>% tolower()
pcMarkerList[['mouse']] <- rownames(mousePCMarkers) %>% tolower()
# All shared markers 
Reduce(intersect, pcMarkerList)
ggVennDiagram(list(pcMarkerList[['human1']], 
                   pcMarkerList[['human2']], 
                   pcMarkerList[['rhesus']], 
                   pcMarkerList[['mouse']]), 
              category.names = c("Human 1",
                                 "Human 2", 
                                 "Rhesus",
                                 "Mouse"), 
              force_upset = F)
# Human shared markers 
Reduce(intersect, list(pcMarkerList[['human1']], pcMarkerList[['human2']])) %>% length()
ggVennDiagram(list(pcMarkerList[['human1']], 
                   pcMarkerList[['human2']]), 
              category.names = c("Human 1",
                                 "Human 2")) +
  coord_flip()
# Primate shared markers
Reduce(intersect, list(pcMarkerList[['human1']], pcMarkerList[['human2']], pcMarkerList[['rhesus']])) 
ggVennDiagram(list(pcMarkerList[['human1']], 
                   pcMarkerList[['human2']], 
                   pcMarkerList[['rhesus']]), 
              category.names = c("Human 1", 
                                 "Human 2", 
                                 "Rhesus"))
# Human and mouse shared markers 
Reduce(intersect, list(pcMarkerList[['human1']], pcMarkerList[['human2']], pcMarkerList[['mouse']])) # human and mouse shared markers
ggVennDiagram(list(pcMarkerList[['human1']], 
                   pcMarkerList[['human2']], 
                   pcMarkerList[['mouse']]), 
              category.names = c("Human 1", 
                                 "Human 2", 
                                 "Mouse"))
# Rhesus and mouse shared markers
Reduce(intersect, list(pcMarkerList[['rhesus']], pcMarkerList[['mouse']])) # rhesus and mouse shared markers
ggVennDiagram(list(pcMarkerList[['rhesus']], 
                   pcMarkerList[['mouse']]),
              category.names = c("Rhesus", 
                                 "Mouse")) + 
  coord_flip()

## 5.5 Subset putative Purkinje cells ----
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
### Plot ----
#### Rhesus ----
#FeaturePlot(putativePCs[['rhesus']], features = c("ALDOC", "MPPED2"), ncol = 1) # ZII+ Markers ----
#FeaturePlot(putativePCs[['rhesus']], features = c("EBF2", "PLCB4"), ncol = 1) # ZII- Markers ----
#FeaturePlot(putativePCs[['rhesus']], features = c("GRID2", "PRKG1")) # Grid2 Markers ----
#### Mouse ----
FeaturePlot(putativePCs[['mouse']], features = c("Aldoc", "Mpped2"), ncol = 1) # ZII+ Markers ----
FeaturePlot(putativePCs[['mouse']], features = c("Ebf2", "Plcb4"), ncol = 1) # ZII- Markers ----
FeaturePlot(putativePCs[['mouse']], features = c("Grid2", "Prkg1")) # Grid2 Markers ----
#### Human 2 ----
FeaturePlot(putativePCs[['human2']], features = c("ALDOC", "MPPED2"), ncol = 1, pt.size = 3) # ZII+ Markers ----
FeaturePlot(putativePCs[['human2']], features = c("EBF2", "PLCB4"), ncol = 1, pt.size = 3) # ZII- Markers ----
FeaturePlot(putativePCs[['human2']], features = c("GRID2", "PRKG1"), pt.size = 3) # Grid2 Markers ----


