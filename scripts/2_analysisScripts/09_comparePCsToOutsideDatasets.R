# 09 compareToOutsideDatasets.R
# Compare Primate Pilot dataset transcriptional similarity to other datasets 
# October 2, 2025
# Hailey Napier

# 0.0 Setup ----
## 0.1 Load packages & functions ----
require(Seurat)
require(dplyr)
require(ggplot2)
require(reshape2)

setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/functions")
source("pseudobulkHumanOrtholog.R")
source("getOrthologCountMat.R")
source("mergeCountMats.R")
source("getDNR.R")

## 0.2 Setup one to one orthologs ----
setwd("/Users/haileynapier/Work/VertGenLab/Projects/zebrinEvolution/Data/geneLists/orthologs")
allPrimateOrthologs <- read.delim("human_MouseMarmosetMacaque_orthologs.txt", sep = "\t", header = T)
names(allPrimateOrthologs) <- tolower(names(allPrimateOrthologs))
allPrimateOrthologs %>%
  filter(macaque.homology.type == "ortholog_one2one") %>%
  filter(mouse.homology.type == "ortholog_one2one") %>%
  filter(white.tufted.ear.marmoset.homology.type == "ortholog_one2one") %>%
  filter(gene.name != "") %>%
  distinct() -> allPrimateOrthologs

## 0.3 Setup primate pilot data ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/seuratObjs")
noGarbage <- readRDS("speciesObjList_cleanCellTypeLabled.rds")
humanPCs_merged <- readRDS("humanPCsMerged.rds")
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/")
pseudobulkMerged <- readRDS("pseudobulkMerged.rds")
pseudobulkMerged %>% 
  as.data.frame() %>%
  select(contains("Purkinje")) -> pseudobulkMerged_pcs
names(pseudobulkMerged_pcs) <- paste(names(pseudobulkMerged_pcs), "Napier", sep = "_")
pseudobulkMerged_pcs <- pseudobulkMerged_pcs[rowSums(pseudobulkMerged_pcs[, -1])>0, ]
rm(pseudobulkMerged)
pseudobulkMerged_pcs$gene.name <- rownames(pseudobulkMerged_pcs)

## 0.4 Setup Hao rhesus data ----
# load data 
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/haoReanalysis/data/seuratObjs")
haoRhesus <- readRDS("macaqueSingleCellPC.rds")
# pseudobulk
pseudobulk_haoRhesus <- pseudobulkHumanOrtholog(haoRhesus, 
                                                groupByIdent = "doublet_info", 
                                                orthologDF = allPrimateOrthologs, 
                                                species = "macaque", 
                                                datasetID = "Purkinje_rhesus_Hao")
# get count matrix in terms of human orthologs
countMat_haoRhesus_pc <- getOrthologCountMat(haoRhesus,
                                             allPrimateOrthologs, 
                                             "macaque")
rm(haoRhesus)

## 0.5 Setup Hao marmoset data (not clean) ----
# load data 
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/haoReanalysis/data/seuratObjs")
haoMarmoset <- readRDS("marmosetSingleCellPC.rds")
# pseudobulk 
pseudobulk_haoMarmoset <- pseudobulkHumanOrtholog(haoMarmoset, 
                                                  groupByIdent = "doublet_info", 
                                                  orthologDF = allPrimateOrthologs, 
                                                  species = "white.tufted.ear.marmoset",
                                                  datasetID = "Purkinje_marmoset_Hao_NOTCLEAN")
# get count matrix in terms of human orthologs
countMat_haoMarmoset_pc <- getOrthologCountMat(haoMarmoset, 
                                               allPrimateOrthologs, 
                                               "white.tufted.ear.marmoset")
rm(haoMarmoset)

## 0.5 Setup Hao marmoset data (CLEAN) ----
# load data 
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/haoReanalysis/data/seuratObjs/cleaned")
haoMarmoset_clean <- readRDS("haoMarmosetPCs_clean.rds")
# pseudobulk 
pseudobulk_haoMarmoset_clean <- pseudobulkHumanOrtholog(haoMarmoset_clean, 
                                                  groupByIdent = "doublet_info", 
                                                  orthologDF = allPrimateOrthologs, 
                                                  species = "white.tufted.ear.marmoset",
                                                  datasetID = "Purkinje_marmoset_Hao")
# get count matrix in terms of human orthologs
countMat_haoMarmoset_pc_clean <- getOrthologCountMat(haoMarmoset_clean, 
                                                     allPrimateOrthologs, 
                                                     "white.tufted.ear.marmoset")
rm(haoMarmoset_clean)

## 0.6 Setup Bartelt mouse data ----
# load data 
barteltMouse <- readRDS("~/Work/VertGenLab/Projects/zebrinEvolution/Data/sequencingData/LukesData/SCA7.8wk.cleaned.rds")
barteltMousePCs <- subset(barteltMouse, idents = "Purkinje cells")
Idents(barteltMousePCs) <- "Type"
barteltMousePCs <- subset(barteltMousePCs, idents = "WT")
rm(barteltMouse)
# pseudobulk
pseudobulk_barteltMousePCs <- pseudobulkHumanOrtholog(barteltMousePCs, 
                                                      groupByIdent = "Type", 
                                                      orthologDF = allPrimateOrthologs, 
                                                      species = "mouse", 
                                                      datasetID = "Purkinje_mouse_Bartelt")
# get count matrix in terms of human orthologs
countMat_barteltMouse_pc <- getOrthologCountMat(barteltMousePCs,
                                                allPrimateOrthologs, 
                                                "mouse")

## 0.6 Merge all pseudobulked DFs into one matrix ----
pseudobulkMerged_pcs <- inner_join(pseudobulkMerged_pcs, pseudobulk_haoRhesus)
pseudobulkMerged_pcs <- inner_join(pseudobulkMerged_pcs, pseudobulk_haoMarmoset)
pseudobulkMerged_pcs <- inner_join(pseudobulkMerged_pcs, pseudobulk_haoMarmoset_clean)
pseudobulkMerged_pcs <- inner_join(pseudobulkMerged_pcs, pseudobulk_barteltMousePCs)
rownames(pseudobulkMerged_pcs) <- pseudobulkMerged_pcs$gene.name
pseudobulkMerged_pcs$gene.name <- NULL
pseudobulkMerged_pcs %>% as.matrix() -> pseudobulkMerged_pcs
rm(pseudobulk_haoRhesus)
rm(pseudobulk_haoMarmoset)
rm(pseudobulk_haoMarmoset_clean)
rm(pseudobulk_barteltMousePCs)

# 1.0 Pairwise Pearson correlations ----
crossDatasetPCCormat <- cor(pseudobulkMerged_pcs)
melted_crossDatasetPCCormat <- melt(crossDatasetPCCormat)
crossDatasetPC_corPlot <- ggplot(data = melted_crossDatasetPCCormat, aes(Var1, Var2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed()
crossDatasetPC_corPlot


# 2.0 Dynamic range all PCs ----
## 2.1 Across all species
### 2.1.1 Get ortholog count matrices for primate pilot data ----
rhesusPCs <- noGarbage[['rhesus']] %>% subset(idents = "Purkinje")
countMat_pilotRhesus_pc <- getOrthologCountMat(rhesusPCs,
                                             allPrimateOrthologs, 
                                             "macaque")
mousePCs <- noGarbage[['mouse']] %>% subset(idents = "Purkinje")
countMat_pilotMouse_pc <- getOrthologCountMat(mousePCs,
                                              allPrimateOrthologs,
                                              "mouse")
human1PCs <- noGarbage[['human1']] %>% subset(idents = "Purkinje")
countMat_pilotHuman1_pcs <- human1PCs@assays$RNA$counts %>% as.matrix()
human2PCs <- noGarbage[['human2']] %>% subset(idents = "Purkinje")
countMat_pilotHuman2_pcs <- human2PCs@assays$RNA$counts %>% as.matrix()

### 2.1.2 Merge count mats for all species for all PCs ----
# (excluding Hao Marmoset because it's weird)
mergedMats <- mergeCountMats(list(countMat_pilotHuman1_pcs, 
                                  countMat_pilotHuman2_pcs, 
                                  countMat_pilotRhesus_pc, 
                                  countMat_pilotMouse_pc, 
                                  countMat_barteltMouse_pc, 
                                  countMat_haoRhesus_pc))
rm(countMat_pilotHuman1_pcs)
rm(countMat_pilotHuman2_pcs)
rm(countMat_barteltMouse_pc)
rm(countMat_haoMarmoset_pc)
rm(countMat_haoMarmoset_pc_clean)
rm(countMat_haoRhesus_pc)
### 2.1.3 Get DNR for all PCs across all species ----
# Single cells
dnr_AllPCs_AcrossSpecies <- getDNR(mergedMats)
range(dnr_AllPCs_AcrossSpecies$dnr, na.rm = T)
hist(dnr_AllPCs_AcrossSpecies$dnr)
# Pseudobulked
pseudobulkAllPCs_dnr <- getDNR(pseudobulkMerged_pcs)
hist(pseudobulkAllPCs_dnr$dnr)


# 3.0 DNR PC subtypes ----
## 3.1 ID subtypes ----
# Here I'm taking the extreme subtype clusters
# Exclude Hao marmoset because of poor cross-species correlation
# Exclude pilot rhesus because not enough cells 
# Bartelt mouse
Idents(barteltMousePCs) <- "seurat_clusters"
barteltMousePCs <- barteltMousePCs %>%
  RunPCA(verbose = F) %>%
  RunUMAP(dims = 1:30, verbose = F) %>%
  FindNeighbors(dims = 1:30, verbose = F) %>%
  FindClusters(verbose = F)
DimPlot(barteltMousePCs, pt.size = 2)
FeaturePlot(barteltMousePCs, features = c("Grid2", "Plcb4", "Aldoc"))
VlnPlot(barteltMousePCs, features = c("Grid2", "Plcb4", "Aldoc"))
barteltMousePCs_zPos <- barteltMousePCs %>% subset(idents = c(3))
barteltMousePCs_zNeg <- barteltMousePCs %>% subset(idents = c(0,2,5))
pseudobulk_bartelt_zPos <- pseudobulkHumanOrtholog(barteltMousePCs_zPos, 
                                                   orthologDF = allPrimateOrthologs, 
                                                   species = "mouse", 
                                                   datasetID = "barteltMouseZPos")
pseudobulk_bartelt_zNeg <- pseudobulkHumanOrtholog(barteltMousePCs_zNeg, 
                                                   orthologDF = allPrimateOrthologs, 
                                                   species = "mouse", 
                                                   datasetID = "barteltMouseZNeg", 
                                                   groupByIdent = "Type")

# Hao rhesus
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/haoReanalysis/data/seuratObjs")
haoRhesus <- readRDS("macaqueSingleCellPC.rds")
haoRhesus_zPos <- haoRhesus %>% 
  subset(idents = c(3))
haoRhesus_zNeg <- haoRhesus %>% 
  subset(idents = c(1,2))
rm(haoRhesus)
pseudobulk_haoRhesus_zNeg <- pseudobulkHumanOrtholog(haoRhesus_zNeg, 
                                                     orthologDF = allPrimateOrthologs, 
                                                     species = "macaque", 
                                                     datasetID = "haoRhesusZNeg", 
                                                     groupByIdent = "doublet_info")
pseudobulk_haoRhesus_zPos <- pseudobulkHumanOrtholog(haoRhesus_zPos, 
                                                     orthologDF = allPrimateOrthologs, 
                                                     species = "macaque", 
                                                     datasetID = "haoRhesusZPos", 
                                                     groupByIdent = "doublet_info")
# Pilot human
humanPCs_zPos <- humanPCs_merged %>%
  subset(idents = 1)
pseudobulk_humanPCs_zPos <- AggregateExpression(humanPCs_zPos, return.seurat = F)
pseudobulk_humanPCs_zPos <- pseudobulk_humanPCs_zPos$RNA %>% as.matrix() %>% as.data.frame()
pseudobulk_humanPCs_zPos$gene.name <- rownames(pseudobulk_humanPCs_zPos)
rownames(pseudobulk_humanPCs_zPos) <- NULL
names(pseudobulk_humanPCs_zPos)[1] <- "pilotHumanZPos"
pseudobulk_humanPCs_zPos$pilotHumanZPos <- pseudobulk_humanPCs_zPos$pilotHumanZPos %>% as.numeric()
humanPCs_zNeg <- humanPCs_merged %>% 
  subset(idents = 0)
pseudobulk_humanPCs_zNeg <- AggregateExpression(humanPCs_zNeg)
pseudobulk_humanPCs_zNeg <- pseudobulk_humanPCs_zNeg$RNA %>% as.matrix() %>% as.data.frame()
pseudobulk_humanPCs_zNeg$gene.name <- rownames(pseudobulk_humanPCs_zNeg) 
rownames(pseudobulk_humanPCs_zNeg) <- NULL
names(pseudobulk_humanPCs_zNeg)[1] <- "pilotHumanZNeg"
pseudobulk_humanPCs_zNeg$pilotHumanZNeg <- pseudobulk_humanPCs_zNeg$pilotHumanZNeg %>% as.numeric()
# Pilot mouse
Idents(mousePCs) <- "seurat_clusters"
mousePCs <- mousePCs %>%
  RunPCA(verbose = F) %>%
  RunUMAP(dims = 1:30, verbose = F) %>%
  FindNeighbors(dims = 1:30, verbose = F) %>%
  FindClusters(verbose = F)
DimPlot(mousePCs, pt.size = 3)
FeaturePlot(mousePCs, features = c("Aldoc", "Grid2", "Plcb4"))
VlnPlot(mousePCs, features = c("Aldoc", "Grid2", "Plcb4"))
mousePCs_zPos <- mousePCs %>%
  subset(idents = 3)
mousePCs_zNeg <- mousePCs %>%
  subset(idents = c(0,1,2))
pseudobulk_pilotMousePCs_zPos <- pseudobulkHumanOrtholog(mousePCs_zPos, 
                                                         orthologDF = allPrimateOrthologs, 
                                                         species = "mouse", 
                                                         datasetID = "pilotMouseZPos")
pseudobulk_pilotMousePCs_zNeg <- pseudobulkHumanOrtholog(mousePCs_zNeg, 
                                                         orthologDF = allPrimateOrthologs, 
                                                         species = "mouse", 
                                                         datasetID = "pilotMouseZNeg")
### 3.1.2 Merge ----
# Z-
pseudobulk_zNeg <- inner_join(pseudobulk_bartelt_zNeg, pseudobulk_haoRhesus_zNeg)
pseudobulk_zNeg <- inner_join(pseudobulk_zNeg, pseudobulk_pilotMousePCs_zNeg)
pseudobulk_zNeg <- inner_join(pseudobulk_zNeg, pseudobulk_humanPCs_zNeg)

# Z+ 
pseudobulk_zPos <- inner_join(pseudobulk_bartelt_zPos, pseudobulk_haoRhesus_zPos)
pseudobulk_zPos <- inner_join(pseudobulk_zPos, pseudobulk_pilotMousePCs_zPos)
pseudobulk_zPos <- inner_join(pseudobulk_zPos, pseudobulk_humanPCs_zPos)

# All PCs, separated by subtype
pseudobulk_PCSubtypes <- inner_join(pseudobulk_zPos, pseudobulk_zNeg)

# Finish matrix construction
rownames(pseudobulk_zPos) <- pseudobulk_zPos$gene.name
pseudobulk_zPos$gene.name <- NULL
pseudobulk_zPos <- pseudobulk_zPos %>% as.matrix()
rownames(pseudobulk_zNeg) <- pseudobulk_zNeg$gene.name
pseudobulk_zNeg$gene.name <- NULL
pseudobulk_zNeg <- pseudobulk_zNeg %>% as.matrix()
rownames(pseudobulk_PCSubtypes) <- pseudobulk_PCSubtypes$gene.name
pseudobulk_PCSubtypes$gene.name <- NULL
#names <- pseudobulk_PCSubtypes %>% colnames() %>% sort()
#pseudobulk_PCSubtypes <- pseudobulk_PCSubtypes[,names]
pseudobulk_PCSubtypes <- pseudobulk_PCSubtypes %>% as.matrix()

## 3.2 Compare PC subtypes ----
# Pearson correlation to compare similarity of PC subtypes
### 3.2.1 Z+ ----
zPosCorMat <- cor(pseudobulk_zPos)
melted_zPosCorMat <- melt(zPosCorMat)
zPos_corPlot <- ggplot(data = melted_zPosCorMat, aes(Var1, Var2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed()
zPos_corPlot

### 3.2.2 All PCs, separated by subtype ----
pcSubtypeCorMat <- cor(pseudobulk_PCSubtypes)
melted_pcSubtypeCorMat <- melt(pcSubtypeCorMat)
pcSubtype_corPlot <- ggplot(data = melted_pcSubtypeCorMat, aes(Var1, Var2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed()
pcSubtype_corPlot 

## 3.3 DNR Z+ ----
### 3.3.1 Get single cell count matrices
human1PCs_zPos_countMat <- humanPCs_zPos@assays$RNA$counts.1 
human2PCs_zPos_countMat <- humanPCs_zPos@assays$RNA$counts.2
humanPCs_zPos_countMat <- mergeCountMats(list(human1PCs_zPos_countMat, human2PCs_zPos_countMat)) 
haoRhesus_zPos_countMat <- getOrthologCountMat(haoRhesus_zPos, 
                                               allPrimateOrthologs,
                                               "macaque")
barteltMousePCs_zPos_countMat <- getOrthologCountMat(barteltMousePCs_zPos, 
                                            allPrimateOrthologs, 
                                            "mouse")
mousePCs_zPos_countMat <- getOrthologCountMat(mousePCs_zPos, 
                                              allPrimateOrthologs, 
                                              "mouse")
zPos_countMat <- mergeCountMats(list(humanPCs_zPos_countMat, 
                                     haoRhesus_zPos_countMat, 
                                     barteltMousePCs_zPos_countMat, 
                                     mousePCs_zPos_countMat))
zPos_dnr <- getDNR(zPos_countMat)
hist(zPos_dnr$dnr)

pseudobulked_zPos_dnr <- getDNR(pseudobulk_zPos)
hist(pseudobulked_zPos_dnr$dnr)

## 3.4 DNR Z- ----
human1PCs_zNeg_countMat <- humanPCs_zNeg@assays$RNA$counts.1 
human2PCs_zNeg_countMat <- humanPCs_zNeg@assays$RNA$counts.2
humanPCs_zNeg_countMat <- mergeCountMats(list(human1PCs_zNeg_countMat, human2PCs_zNeg_countMat)) 
haoRhesus_zNeg_countMat <- getOrthologCountMat(haoRhesus_zNeg, 
                                               allPrimateOrthologs,
                                               "macaque")
barteltMousePCs_zNeg_countMat <- getOrthologCountMat(barteltMousePCs_zNeg, 
                                                     allPrimateOrthologs, 
                                                     "mouse")
mousePCs_zNeg_countMat <- getOrthologCountMat(mousePCs_zNeg, 
                                              allPrimateOrthologs, 
                                              "mouse")
zNeg_countMat <- mergeCountMats(list(humanPCs_zNeg_countMat, 
                                     haoRhesus_zNeg_countMat, 
                                     barteltMousePCs_zNeg_countMat, 
                                     mousePCs_zNeg_countMat))
zNeg_dnr <- getDNR(zNeg_countMat)
hist(zNeg_dnr$dnr)

pseudobulked_zNeg_dnr <- getDNR(pseudobulk_zNeg)
hist(pseudobulked_zNeg_dnr$dnr)

# 4.0 Standard deviation by subtype across species ----
## 4.1 Z- ----
zNegSD <- getNonZeroSD(zNeg_countMat)
zNegSD[zNegSD != 0] %>% hist(breaks = 700)
range(zNegSD, na.rm=T)
sort(zNegSD, decreasing = T) %>% head()

## 4.2 Z+ ----
zPosSD <- getNonZeroSD(zPos_countMat)
zPosSD[zPosSD != 0] %>% hist(breaks = 700)
range(zPosSD, na.rm=T)
sort(zPosSD, decreasing = T) %>% head()


# 5.0 Mean normalized variance by subtype across species ----
## 5.1 Z- ----
zNegVar <- meanNormVar(zNeg_countMat)
zNegVar[zNegVar != 0] %>% hist(breaks = 700)
range(zNegVar, na.rm = T)
sort(zNegVar, decreasing = T) %>% head()
zNegVar_log <- log(zNegVar)
hist(zNegVar_log)
sort(zNegVar_log, decreasing = T) %>% head()

## 5.2 Z+ ----
zPosVar <- meanNormVar(zPos_countMat)
zPosVar[zPosVar != 0] %>% hist(breaks = 700)
range(zPosVar, na.rm = T)
sort(zPosVar, decreasing = T) %>% head()
zPosVar_log <- log(zPosVar)
hist(zPosVar_log)
sort(zPosVar_log, decreasing = T) %>% head()
