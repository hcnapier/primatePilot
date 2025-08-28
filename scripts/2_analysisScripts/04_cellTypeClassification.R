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
set.seed(1234)

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
DimPlot(negSingList[['rhesus']], label = T)
MDimPlot(negSingList[['mouse']])

## 5.3 Plot cell-type markers ----
Idents(negSingList[["human1"]]) <- "seurat_clusters"
Idents(negSingList[["human2"]]) <- "seurat_clusters"
Idents(negSingList[["rhesus"]]) <- "seurat_clusters"
Idents(negSingList[["mouse"]]) <- "seurat_clusters"
cluster.markers <- list()
cluster.markers <- lapply(negSingList, FindAllMarkers, only.pos = TRUE)
cluster.markers[['mouse']] %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_mouse
DoHeatmap(negSingList[['mouse']], features = top10_mouse$gene) + 
  NoLegend() + 
  theme(axis.text.y = element_text(size = 5))

cluster.markers[['rhesus']] %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_rhesus
DoHeatmap(negSingList[['rhesus']], features = top10_rhesus$gene) + 
  NoLegend() + 
  theme(axis.text.y = element_text(size = 5))

cluster.markers[['human1']] %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_human1
DoHeatmap(negSingList[['human1']], features = top10_human1$gene) + 
  NoLegend() + 
  theme(axis.text.y = element_text(size = 5))

cluster.markers[['human2']] %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_human2
DoHeatmap(negSingList[['human2']], features = top10_human2$gene) + 
  NoLegend() + 
  theme(axis.text.y = element_text(size = 5))

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
#### MLIs 
#### MLI 1 -- 2, 4, 10
#### MLI 2 -- 1, 6
FeaturePlot(negSingList[['human1']], features = c("INPP4B", "SLC24A3", "KIT", "CACNA2D3"))
VlnPlot(negSingList[['human1']], features = c("INPP4B", "SLC24A3", "KIT", "CACNA2D3", "ADARB2"))
VlnPlot(negSingList[['human1']], features = c("SORCS3", "NXPH1"))
#### PLIs -- 13
VlnPlot(negSingList[['human1']], features = c("NXPH2", "TTBK2", "DNAH14", "VCAN"))
FeaturePlot(negSingList[['human1']], features = c("NXPH2", "TTBK2", "DNAH14", "VCAN"))
#### Fibroblasts -- ??? 
FeaturePlot(negSingList[['human1']], features = c("PDGFRA", "PDGFRB", "THY1", "VIM"))
VlnPlot(negSingList[['human1']], features = c("CPED1", "APOD", "BNC2", "BMP6"))
#### Oligodendrocytes -- 5
VlnPlot(negSingList[['human1']], features = c("ST18", "PRR5L", "MBP", "SEC14L5", "DOCK10"))
FeaturePlot(negSingList[['human1']], features = c("ST18", "PRR5L", "MBP", "SEC14L5", "DOCK10"))
#### Endothelial Stalk -- 11?
VlnPlot(negSingList[['human1']], features = c("FLT1", "SLCO1A4", "ADGRL4", "MECOM", "CXCL12"))
#### Endothelial Mural -- ???
VlnPlot(negSingList[['human1']], features = c("ATP13A5", "PDE8B", "CALD1", "PDGFRB", "PLCL1"))



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
#### MLIs 
#### MLI1 -- 9, 4, 2 
#### MLI2 -- 1, 8
FeaturePlot(negSingList[['human2']], features = c("INPP4B", "SLC24A3", "KIT", "CACNA2D3"))
VlnPlot(negSingList[['human2']], features = c("INPP4B", "SLC24A3", "KIT", "CACNA2D3", "ADARB2"))
VlnPlot(negSingList[['human2']], features = c("SORCS3", "NXPH1"))
#### PLIs -- 13
VlnPlot(negSingList[['human2']], features = c("NXPH2", "TTBK2", "DNAH14", "VCAN"))
FeaturePlot(negSingList[['human2']], features = c("NXPH2", "TTBK2", "DNAH14", "VCAN"))
#### Fibroblasts -- 6????
FeaturePlot(negSingList[['human2']], features = c("CPED1", "APOD", "BNC2", "BMP6"))
#### Oligodendrocytes -- 11
VlnPlot(negSingList[['human2']], features = c("ST18", "PRR5L", "MBP", "SEC14L5", "DOCK10"))
FeaturePlot(negSingList[['human2']], features = c("ST18", "PRR5L", "MBP", "SEC14L5", "DOCK10"))
#### Endothelial Stalk -- 6???
VlnPlot(negSingList[['human2']], features = c("FLT1", "SLCO1A4", "ADGRL4", "MECOM", "CXCL12"))
#### Endothelial Mural --
VlnPlot(negSingList[['human2']], features = c("ATP13A5", "PDE8B", "CALD1", "PDGFRB", "PLCL1"))
#### Choroid


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
VlnPlot(negSingList[['rhesus']], features = c("CHN2", "TENM1", "CADPS2", "CA10", "GRIK2", "RBFOX3", "CNTN2", "ETV1"))
FindMarkers(negSingList[['rhesus']], ident.1 = "9", ident.2 = NULL)
#### Golgi Cells -- 5 (minor 14)
FeaturePlot(negSingList[['rhesus']], features = c("GABRG3", "CACNG3", "ELAVL2", "KCNN3"))
VlnPlot(negSingList[['rhesus']], features = c("GABRG3", "CACNG3", "ELAVL2", "KCNN3", "LGI2"))
#### MLIs
#### MLI1 -- 1, 6
#### MLI2 -- 0, 8
FeaturePlot(negSingList[['rhesus']], features = c("INPP4B", "SLC24A3", "KIT", "CACNA2D3"))
VlnPlot(negSingList[['rhesus']], features = c("INPP4B", "SLC24A3", "KIT", "CACNA2D3", "ADARB2"))
VlnPlot(negSingList[['rhesus']], features = c("SORCS3", "NXPH1"))
# dif between 1 and 6
FindMarkers(negSingList[['rhesus']], ident.1 = "1", ident.2 = "6")
# cluster 10
FindMarkers(negSingList[['rhesus']], ident.1 = "10", ident.2 = NULL)
FindMarkers(negSingList[['rhesus']], ident.1 = "10", ident.2 = c("0","8"))
#### PLIs -- 10
FeaturePlot(negSingList[['rhesus']], features = c("NXPH2", "TTBK2", "DNAH14", "VCAN"))
VlnPlot(negSingList[['rhesus']], features = c("NXPH2", "TTBK2", "DNAH14", "VCAN"))
#### Fibroblasts -- ??? 
FeaturePlot(negSingList[['rhesus']], features = c("CPED1", "APOD", "BNC2", "BMP6"))
#### Oligodendrocytes -- 11
VlnPlot(negSingList[['rhesus']], features = c("ST18", "PRR5L", "MBP", "SEC14L5", "DOCK10"))
FeaturePlot(negSingList[['rhesus']], features = c("ST18", "PRR5L", "MBP", "SEC14L5", "DOCK10"))
#### Endothelial Stalk --????
VlnPlot(negSingList[['rhesus']], features = c("FLT1", "SLCO1A4", "ADGRL4", "MECOM", "CXCL12"))
#### Endothelial Mural -- ???
VlnPlot(negSingList[['rhesus']], features = c("ATP13A5", "PDE8B", "CALD1", "PDGFRB", "PLCL1"))
#### Choroid -- ???
VlnPlot(negSingList[['rhesus']], features = c("OTX2OS1", "HTR2C", "COL8A1", "RBM47", "GMNC"))
#### Cluster 2
FindMarkers(negSingList[['rhesus']], ident.1 = "2", ident.2 = NULL)
VlnPlot(negSingList[['rhesus']], features = c("LOC703224", "KHDRBS3", "HS3ST4", "FRMD4B", "GULP1"))
FeaturePlot(negSingList[['rhesus']], features = "nCount_RNA")
VlnPlot(negSingList[['rhesus']], features = c("GRID2", "SORCS3", # Basket cells
                                              "GAD1", "GAD2", # GABAergic cells
                                              "GRM4", # Granule neurons
                                              "CFH", # vascular cells
                                              "TH", "DBH", # dopaminergic cells
                                              "RBFOX3", "APOA2", "KIR3DL1", "KIR3DL2")) 

#### Cluster 4
FindMarkers(negSingList[['rhesus']], ident.1 = '4', ident.2 = NULL)
#### Cluster 7 -- Garbage
FindMarkers(negSingList[['rhesus']], ident.1 = '7', ident.2 = NULL)
#### Cluster 13 
FindMarkers(negSingList[['rhesus']], ident.1 = '13', ident.2 = NULL)
VlnPlot(negSingList[['rhesus']], features = c("BCAS1", "MAG", "MOBP", "MOG", "NG2"))
FindMarkers(negSingList[['rhesus']], ident.1 = "13", ident.2 ="11")


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
#### MLIs
#### MLI1 -- 0,2
#### MLI2 -- 1, 6
FeaturePlot(negSingList[['mouse']], features = c("Inpp4b", "Slc24a3", "Kit", "Cacna2d3"))
VlnPlot(negSingList[['mouse']], features = c("Inpp4b", "Slc24a3", "Kit", "Cacna2d3", "Adarb2"))
# Within MLIs, Sorcs3 is only in MLI1 and Nxph1 is only in MLI2 (Kozareva et al. 2021, )
VlnPlot(negSingList[['mouse']], features = c("Sorcs3", "Lypd6", "Ptprk", "Nxph1", "Cdh22"))
#### PLIs -- 8, 17 (not sure what's happening with 19 & 18)
FeaturePlot(negSingList[['mouse']], features = c("Grm5", "Robo1", "Lrrc7", "Rtl4"))
VlnPlot(negSingList[['mouse']], features = c("Grm5", "Robo1", "Lrrc7", "Rtl4", "Nhs"))
FindMarkers(negSingList[['mouse']], ident.1 = "8", ident.2 = NULL)
#### Fibroblasts -- 15 
FeaturePlot(negSingList[['mouse']], features = c("Cped1", "Apod", "9530026P05Rik", "Bmp6"))
VlnPlot(negSingList[['mouse']], features = c("Cped1", "Apod", "9530026P05Rik", "Bmp6", "Bnc2"))
FindMarkers(negSingList[['mouse']], ident.1 = "11", ident.2 = NULL)
#### Oligodendrocytes -- 10
FeaturePlot(negSingList[['mouse']], features = c("St18", "Prr5l", "Mbp", "Sec14l5"))
VlnPlot(negSingList[['mouse']], features = c("St18", "Prr5l", "Mbp", "Sec14l5", "Dock10"))
#### Endothelial stalk --
FeaturePlot(negSingList[['mouse']], features = c("Flt1", "Slco1a4", "Adgrl4", "Mecom"))
VlnPlot(negSingList[['mouse']], features = c("Flt1", "Slco1a4", "Adgrl4", "Mecom", "Cxcl12"))
#### Endothelial mural -- 11
VlnPlot(negSingList[['mouse']], features = c("Atp13a5", "Pde8b", "Cald1", "Pdgfrb", "Plcl1"))
FeaturePlot(negSingList[['mouse']], features = c("Atp13a5", "Pde8b", "Cald1", "Pdgfrb"))
#### Choroid -- 16
VlnPlot(negSingList[['mouse']], features = c("Otx2os1", "Htr2c", "Col8a1", "Rbm47", "Gmnc"))
FeaturePlot(negSingList[['mouse']], features = c("Otx2os1", "Htr2c", "Col8a1", "Rbm47"))
#### Cluster 9 
FindMarkers(negSingList[['mouse']], ident.1 = "9", ident.2 = NULL)
FindMarkers(negSingList[['mouse']], ident.1 = "9", ident.2 = c("2", "0"))
#### Deep cerebellar nuclei? -- 13
FindMarkers(negSingList[['mouse']], ident.1 = '13', ident.2 = NULL)
FeaturePlot(negSingList[['mouse']], features = c("Lmx1a", "Pax6", "Eomes", "Snca"))
VlnPlot(negSingList[['mouse']], features = c("Lmx1a", "Pax6", "Eomes", "Snca"))
#### Cluster 18 
FindMarkers(negSingList[['mouse']], ident.1 = '18', ident.2 = NULL)
FeaturePlot(negSingList[['mouse']], features = c("nCount_RNA", "nFeature_RNA"))
VlnPlot(negSingList[['mouse']], features = c("Inhba", "Htr2a", "P3h2"))
#### DCN cells 
VlnPlot(negSingList[['mouse']], features = c("Snap25", "Otx2"))


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

## Try to find junk cells ----
noPCs <- list()
noPCs[['mouse']] <- subset(negSingList[['mouse']],
                           idents = c("3","7"),
                           invert = T)
FeaturePlot(noPCs[['mouse']], feature = "nCount_RNA")

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
DimPlot(putativePCs[['mouse']])
msPCSubtypeMarkers <- FindMarkers(putativePCs[['mouse']], ident.1 = "3", ident.2 = "2") %>% filter(p_val_adj < 0.05)
msPCSubtypeMarkerNames <- msPCSubtypeMarkers %>% row.names() %>% tolower()
FindMarkers(putativePCs[['mouse']], ident.1 = "1", ident.2 = "0")
FeaturePlot(putativePCs[['mouse']], features = c("Aldoc", "Mpped2", "Itpr2"), ncol = 1) # ZII+ Markers ----
FeaturePlot(putativePCs[['mouse']], features = c("Ebf2", "Plcb4"), ncol = 1) # ZII- Markers ----
FeaturePlot(putativePCs[['mouse']], features = c("Grid2", "Prkg1")) # Grid2 Markers ----
VlnPlot(putativePCs[['mouse']], features = c("Grid2", "Prkg1", "Gria4", "Rbfox1")) # Grid2 Markers ----
VlnPlot(putativePCs[['mouse']], features = c("Grid2", "Cdh18", "Dlg2", "Kcnab1", "Kctd8", "Large1")) 

#### Human 2 ----
DimPlot(putativePCs[['human2']], pt.size = 3)
h2PCSubtypeMarkers <- FindMarkers(putativePCs[['human2']], ident.1 = '0', ident.2 = '1') %>% filter(p_val_adj < 0.05) 
h2PCSubtypeMarkerNames <- h2PCSubtypeMarkers %>% row.names() %>% tolower() 
FeaturePlot(putativePCs[['human2']], features = c("ALDOC", "MPPED2"), ncol = 1, pt.size = 3) # ZII+ Markers ----
FeaturePlot(putativePCs[['human2']], features = c("EBF2", "PLCB4"), ncol = 1, pt.size = 3) # ZII- Markers ----
FeaturePlot(putativePCs[['human2']], features = c("GRIK1", "SORCS2", "TENM2", "DGKH"), ncol = 1, pt.size = 3) # ZII- Markers ----
VlnPlot(putativePCs[['human2']], features = c("RNR2", "GRIK1", "SORCS2"))
FeaturePlot(putativePCs[['human2']], features = c("GRID2", "PRKG1"), pt.size = 3) # Grid2 Markers ----
VlnPlot(putativePCs[['human2']], features = c("CDH18", "DLG2", "GRID2", "KCNAB1", "KCTD8", "LARGE1"))
# Compare human and mouse subtype markers 
ggVennDiagram(list(h2PCSubtypeMarkerNames, msPCSubtypeMarkerNames), 
              category.names = c("Human", 
                                 "Mouse")) +
  coord_flip() +
  labs(title = "Purkinje Cell Subtype Marker Genes")
intersect(h2PCSubtypeMarkerNames, msPCSubtypeMarkerNames) %>% sort()
h2PCSubtypeMarkerNames

## Assign new cluster IDs ----
Idents(negSingList[["mouse"]]) <- "seurat_clusters"
newIDs_mouse <- c("MLI Type 1", 
                  "MLI Type 2", 
                  "MLI Type 1", 
                  "Purkinje Cell", 
                  "Bergmann Glia", 
                  "Golgi Cell", 
                  "MLI Type 2", 
                  "Purkinje Cell", 
                  "Purkinje Layer Interneuron", 
                  "Garbage", 
                  "Oligodendrocyte", 
                  "Endothelial Mural", 
                  "Granule Neuron", 
                  "Unknown1",
                  "DCN?", 
                  "Fibroblast", 
                  "Choroid", 
                  "Purkinje Layer Interneuron", 
                  "Garbage", 
                  "Golgi Cell", 
                  "Endothelial Stalk")
names(newIDs_mouse) <- levels(negSingList[['mouse']])
negSingList[['mouse']] <- RenameIdents(negSingList[['mouse']],newIDs_mouse)
DimPlot(negSingList[['mouse']], reduction = 'umap', label = T)
Idents(negSingList[["rhesus"]]) <- "seurat_clusters"
newIDs_rhesus <- c("MLI Type 2",
                   "MLI Type 1", 
                   "Unknown1", 
                   "Bergmann Glia", 
                   "Garbage",
                   "Golgi Cells", 
                   "MLI Type 1", 
                   "Garbage", 
                   "MLI Type 2", 
                   "Granule Neuron", 
                   "Purkinje Layer Interneuron", 
                   "Oligodendrocyte", 
                   "Purkinje Cell", 
                   "Immature Oligodendrocyte??",
                   "Golgi Cells??", 
                   "Unknown2")
names(newIDs_rhesus) <- levels(negSingList[['rhesus']])
negSingList[['rhesus']] <- RenameIdents(negSingList[['rhesus']], newIDs_rhesus)
DimPlot(negSingList[['rhesus']], reduction = 'umap', label = T)
Idents(negSingList[["human1"]]) <- "seurat_clusters"
newIDs_human1 <- c("Granule Neuron", 
                   "MLI Type 2", 
                   "MLI Type 1", 
                   "Granule Neuron", 
                   "MLI Type 1", 
                   "Oligodendrocyte", 
                   "MLI Type 2", 
                   "Golgi Cell",
                   "Bergmann Glia", 
                   "Garbage", 
                   "MLI Type 1", 
                   "Garbage", 
                   "Unkown", 
                   "Purkinje Layer Interneuron", 
                   "Astrocyte", 
                   "Purkinje Cell")
names(newIDs_human1) <- levels(negSingList[['human1']])
negSingList[['human1']] <- RenameIdents(negSingList[['human1']], newIDs_human1)
DimPlot(negSingList[['human1']], reduction = 'umap', label = T)
Idents(negSingList[["human2"]]) <- "seurat_clusters"
newIDs_human2 <- c("Granule Neuron", 
                   "MLI Type 2a", 
                   "MLI Type 1a", 
                   "Granule Neuron", 
                   "MLI Type 1b", 
                   "Bergmann Glia", 
                   "Garbage", 
                   "Golgi Cell", 
                   "MLI Type 2b", 
                   "MLI Type 1c", 
                   "Unknown 1",
                   "Oligodendrocyte", 
                   "Unknown 2",
                   "Purkinje Layer Interneuron", 
                   "Purkinje Cell", 
                   "Golgi Cell")
names(newIDs_human2) <- levels(negSingList[['human2']])
negSingList[['human2']] <- RenameIdents(negSingList[['human2']], newIDs_human2)
DimPlot(negSingList[['human2']], reduction = 'umap', label = T)

## Marker gene heatmaps
cellTypeMarkers <- list()
cellTypeMarkers <- lapply(negSingList, FindAllMarkers, only.pos = TRUE)
cellTypeMarkers[['mouse']] %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_mouse
DoHeatmap(negSingList[['mouse']], features = top10_mouse$gene) + 
  NoLegend() + 
  theme(axis.text.y = element_text(size = 5))

cellTypeMarkers[['rhesus']] %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_rhesus
DoHeatmap(negSingList[['rhesus']], features = top10_rhesus$gene) + 
  NoLegend() + 
  theme(axis.text.y = element_text(size = 5))

cellTypeMarkers[['human1']] %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_human1
DoHeatmap(negSingList[['human1']], features = top10_human1$gene) + 
  NoLegend() + 
  theme(axis.text.y = element_text(size = 5))

cellTypeMarkers[['human2']] %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_human2
DoHeatmap(negSingList[['human2']], features = top10_human2$gene) + 
  NoLegend() + 
  theme(axis.text.y = element_text(size = 5))

# Remove Garbage Cells & Recluster ----
noGarbage <- lapply(negSingList, subset, idents = "Garbage", invert = T)
noGarbage <- lapply(noGarbage, RunPCA, verbose = F)
noGarbage <- lapply(noGarbage, RunUMAP, dims = 1:30, verbose = F)
noGarbage <- lapply(noGarbage, FindNeighbors, dims = 1:30, verbose = F)
noGarbage <- lapply(noGarbage, FindClusters, verbose = F)
# pass 2
noGarbage[['human1']] <- subset(noGarbage[['human1']], idents = "6", invert = T)
noGarbage[['human2']] <- subset(noGarbage[['human2']], idents = "11", invert = T)
noGarbage[['mouse']] <- subset(noGarbage[['mouse']], idents = c("8", "16"), invert = T)
noGarbage[['rhesus']] <- subset(noGarbage[['rhesus']], idents = c("14", "15"), invert = T)
noGarbage <- lapply(noGarbage, RunPCA, verbose = F)
noGarbage <- lapply(noGarbage, RunUMAP, dims = 1:30, verbose = F)
noGarbage <- lapply(noGarbage, FindNeighbors, dims = 1:30, verbose = F)
noGarbage <- lapply(noGarbage, FindClusters, verbose = F)
# pass 3
noGarbage[['human1']] <- subset(noGarbage[['human1']], idents = "10", invert = T)
noGarbage[['human1']] <- RunPCA(noGarbage[['human1']], verbose = F)
noGarbage[['human1']] <- RunUMAP(noGarbage[['human1']], dims = 1:30, verbose = F)
noGarbage[['human1']] <- FindNeighbors(noGarbage[['human1']], dims = 1:30, verbose = F)
noGarbage[['human1']] <- FindClusters(noGarbage[['human1']], verbose = F)

## Assign cell types ----
noGarbageMarkers <- list()
noGarbageMarkers <- lapply(noGarbage, FindAllMarkers, only.pos = TRUE)
top10NoGarbage <- list()
sampleNames <- c("mouse", "rhesus", "human1", "human2")
for(currSample in sampleNames){
  noGarbageMarkers[[currSample]] %>%
    group_by(cluster) %>%
    dplyr::filter("av_log2FC" > 1) %>%
    dplyr::slice_head(n = 10) %>%
    ungroup() -> top10NoGarbage[[currSample]]
}
cellTypeMarkers[['human1']] %>% filter(gene == "GRM8")
DimPlot(noGarbage[['mouse']], label = T)
DoHeatmap(noGarbage[['human1']], features = top10NoGarbage[['human1']]$gene) + 
  NoLegend() + 
  theme(axis.text.y = element_text(size = 4))
top10NoGarbage[['human1']] %>% filter(cluster == "3")

top10NoGarbage[['mouse']] %>% filter(gene == "Eomes")
noGarbageMarkers[['human1']] %>% filter(gene == "SLC1A3")

## Assign new cluster IDs
newIDs_human1 <- c("Granule", 
                   "MLI1", 
                   "MLI2",
                   "Granule", 
                   "MLI1", 
                   "ODC", 
                   "Golgi", 
                   "Bergmann", 
                   "UBC",
                   "Granule", 
                   "PLI", 
                   "Purkinje")
names(newIDs_human1) <- levels(noGarbage[['human1']])
noGarbage[['human1']] <- RenameIdents(noGarbage[['human1']], newIDs_human1)
DimPlot(noGarbage[['human1']], reduction = 'umap', label = F)
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
                   "Unipolar Brush Cell",
                   "Oligodendrocyte",
                   "Purkinje Layer Interneuron", 
                   "Purkinje Cell", 
                   "Golgi Cell")
names(newIDs_human2) <- levels(noGarbage[['human2']])
noGarbage[['human2']] <- RenameIdents(noGarbage[['human2']], newIDs_human2)
DimPlot(noGarbage[['human2']], reduction = 'umap', label = F)
newIDs_rhesus <- c("MLI2", 
                   "MLI1", 
                   "UBC", 
                   "Bergmann", 
                   "MLI1", 
                   "Golgi", 
                   "MLI2", 
                   "Golgi", 
                   "Granule", 
                   "PLI", 
                   "ODC", 
                   "Purkinje", 
                   "OPC")
names(newIDs_rhesus) <- levels(noGarbage[['rhesus']])
noGarbage[['rhesus']] <- RenameIdents(noGarbage[['rhesus']], newIDs_rhesus)
DimPlot(noGarbage[['rhesus']], reduction = 'umap', label = T)
Idents(noGarbage[["mouse"]]) <- "seurat_clusters"
newIDs_mouse <- c("MLI1", 
                  "MLI2", 
                  "MLI1", 
                  "Purkinje", 
                  "Bergmann", 
                  "Golgi", 
                  "PLI", 
                  "Purkinje", 
                  "ODC",
                  "Endo. Mural", 
                  "UBC", 
                  "Granule", 
                  "Fibroblast", 
                  "Astrocyte", 
                  "Choroid", 
                  "Golgi", 
                  "Endo. Stalk")
names(newIDs_mouse) <- levels(noGarbage[['mouse']])
noGarbage[['mouse']] <- RenameIdents(noGarbage[['mouse']], newIDs_mouse)
DimPlot(noGarbage[['mouse']], reduction = 'umap', label = T)

# Save Objects ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/seuratObjs")
saveRDS(noGarbage, "speciesObjList_cleanCellTypeLabled.rds")


