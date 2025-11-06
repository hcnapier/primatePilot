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
require(clusterProfiler)
require(org.Hs.eg.db)

setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/functions")
source("pseudobulkHumanOrtholog.R")
source("getOrthologCountMat.R")
source("mergeCountMats.R")
source("getDNR.R")
source("meanNormVar.R")
source("getNonZeroSD.R")
source("getNonZeroVar.R")

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
allPrimateOrthologs$gene...gc.content <- NULL

## 0.3 Setup primate pilot data ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/seuratObjs")
noGarbage <- readRDS("speciesObjList_cleanCellTypeLabled.rds")
humanPCs_merged <- readRDS("humanPCsMerged.rds")

## 0.4 Setup Hao rhesus data ----
# load data 
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/haoReanalysis/data/seuratObjs")
haoRhesus <- readRDS("macaqueSingleCellPC.rds")
# get count matrix in terms of human orthologs
countMat_haoRhesus_pc <- getOrthologCountMat(haoRhesus,
                                             allPrimateOrthologs, 
                                             "macaque")
rm(haoRhesus)
countMat_haoRhesus_pc %>% colSums() %>% hist(main = "Hao Rhesus", xlab = "Read counts per cell", breaks = 14)
countMat_haoRhesus_pc %>% colSums() %>% mean()

## 0.6 Setup Bartelt mouse data ----
# load data 
barteltMouse <- readRDS("~/Work/VertGenLab/Projects/zebrinEvolution/Data/sequencingData/LukesData/SCA7.8wk.cleaned.rds")
barteltMousePCs <- subset(barteltMouse, idents = "Purkinje cells")
Idents(barteltMousePCs) <- "Type"
barteltMousePCs <- subset(barteltMousePCs, idents = "WT")
rm(barteltMouse)
# get count matrix in terms of human orthologs
countMat_barteltMouse_pc <- getOrthologCountMat(barteltMousePCs,
                                                allPrimateOrthologs, 
                                                "mouse")
countMat_barteltMouse_pc %>% colSums() %>% hist(main = "Bartelt Mouse", xlab = "Read counts per cell", breaks = 14)
countMat_barteltMouse_pc %>% colSums() %>% mean()


# 2.0 Dynamic range all PCs ----
## 2.1 Across all species
### 2.1.1 Get ortholog count matrices for primate pilot data ----
rhesusPCs <- noGarbage[['rhesus']] %>% subset(idents = "Purkinje")
countMat_pilotRhesus_pc <- getOrthologCountMat(rhesusPCs,
                                             allPrimateOrthologs, 
                                             "macaque", 
                                             normMethod = "shiftedLog", 
                                             useForL = 7208)
countMat_pilotRhesus_pc %>% colSums() %>% hist(main = "Pilot Rhesus", xlab = "Read counts per cell", breaks = 14)
countMat_pilotRhesus_pc %>% colSums() %>% mean()

mousePCs <- noGarbage[['mouse']] %>% subset(idents = "Purkinje")
countMat_pilotMouse_pc <- getOrthologCountMat(mousePCs,
                                              allPrimateOrthologs,
                                              "mouse", 
                                              normMethod = "shiftedLog", 
                                              useForL = 7208)
countMat_pilotMouse_pc %>% colSums() %>% hist(main = "Pilot Mouse", xlab = "Read counts per cell", breaks = 14)
countMat_pilotMouse_pc %>% colSums() %>% mean()

humanGeneDF <- allPrimateOrthologs %>%
  dplyr::select(c("gene.name", "gene.size"))
human1PCs <- noGarbage[['human1']] %>% subset(idents = "Purkinje")
countMat_pilotHuman1_pcs <- human1PCs@assays$RNA$counts %>% as.matrix() %>% shiftedLogNorm(pseudocount = 1, useForL = 7208)
countMat_pilotHuman1_pcs %>% colSums() %>% hist(main = "Pilot Human1", xlab = "Read counts per cell")
countMat_pilotHuman1_pcs %>% colSums() %>% mean()

human2PCs <- noGarbage[['human2']] %>% subset(idents = "Purkinje")
countMat_pilotHuman2_pcs <- human2PCs@assays$RNA$counts %>% as.matrix() %>% shiftedLogNorm(pseudocount = 1, useForL = 7208)
countMat_pilotHuman2_pcs %>% colSums() %>% hist(main = "Pilot Human2", xlab = "Read counts per cell")
countMat_pilotHuman2_pcs %>% colSums() %>% mean()

### 2.1.2 Merge count mats for all species for all PCs ----
# (excluding Hao Marmoset because it's weird)
mergedMats <- mergeCountMats(list(countMat_pilotHuman1_pcs, 
                                  countMat_pilotHuman2_pcs, 
                                  countMat_pilotRhesus_pc, 
                                  countMat_pilotMouse_pc, 
                                  countMat_barteltMouse_pc, 
                                  countMat_haoRhesus_pc))
mergedMats %>% colSums() %>% hist(main = "All Species Merged", xlab = "Read counts per cell", breaks = 100)
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

# Hao rhesus
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/haoReanalysis/data/seuratObjs")
haoRhesus <- readRDS("macaqueSingleCellPC.rds")
haoRhesus_zPos <- haoRhesus %>% 
  subset(idents = c(3))
haoRhesus_zNeg <- haoRhesus %>% 
  subset(idents = c(1))
rm(haoRhesus)

# Pilot human
humanPCs_zPos <- humanPCs_merged %>%
  subset(idents = 1)

humanPCs_zNeg <- humanPCs_merged %>% 
  subset(idents = 0)

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
  subset(idents = c(2))

## 3.3 DNR Z+ ----
### 3.3.1 Get single cell count matrices
human1PCs_zPos_countMat <- humanPCs_zPos@assays$RNA$counts.1 %>% as.matrix() %>% shiftedLogNorm(pseudocount = 1, removeLowCounts = T, lowCountThreshold = 0.5)
human2PCs_zPos_countMat <- humanPCs_zPos@assays$RNA$counts.2 %>% as.matrix() %>% shiftedLogNorm(pseudocount = 1)
humanPCs_zPos_countMat <- mergeCountMats(list(human1PCs_zPos_countMat, human2PCs_zPos_countMat)) 
humanPCs_zPos_countMat %>% colSums() %>% hist(main = "Human Z+ PCs", xlab = "Read counts per cell")
haoRhesus_zPos_countMat <- getOrthologCountMat(haoRhesus_zPos, 
                                               allPrimateOrthologs,
                                               "macaque", 
                                               normMethod = "shiftedLog",
                                               useForL = 6687)
haoRhesus_zPos_countMat %>% colSums() %>% hist(main = "Hao Rhesus Z+ PCs", xlab = "Read counts per cell")
barteltMousePCs_zPos_countMat <- getOrthologCountMat(barteltMousePCs_zPos, 
                                            allPrimateOrthologs, 
                                            "mouse", 
                                            normMethod = "shiftedLog", 
                                            useForL = 6687)
barteltMousePCs_zPos_countMat %>% colSums() %>% hist(main = "Bartelt Mouse Z+ PCs", xlab = "Read counts per cell")
mousePCs_zPos_countMat <- getOrthologCountMat(mousePCs_zPos, 
                                              allPrimateOrthologs, 
                                              "mouse", 
                                              normMethod = "shiftedLog", 
                                              useForL = 6687)
mousePCs_zPos_countMat %>% colSums() %>% hist(main = "Pilot Mouse Z+ PCs", xlab = "Read counts per cell")
zPos_countMat <- mergeCountMats(list(humanPCs_zPos_countMat, 
                                     haoRhesus_zPos_countMat, 
                                     barteltMousePCs_zPos_countMat, 
                                     mousePCs_zPos_countMat))
zPos_countMat %>% colSums() %>% hist(main = "All Species Merged", xlab = "Read counts per cell", breaks = 100)
zPos_dnr <- getDNR(zPos_countMat)
hist(zPos_dnr$dnr)

pseudobulked_zPos_dnr <- getDNR(pseudobulk_zPos)
hist(pseudobulked_zPos_dnr$dnr)

## 3.4 DNR Z- ----
human1PCs_zNeg_countMat <- humanPCs_zNeg@assays$RNA$counts.1 %>% as.matrix() %>% shiftedLogNorm(pseudocount = 1, useForL = 7861)
human2PCs_zNeg_countMat <- humanPCs_zNeg@assays$RNA$counts.2 %>% as.matrix() %>% shiftedLogNorm(pseudocount = 1, useForL = 7861)
humanPCs_zNeg_countMat <- mergeCountMats(list(human1PCs_zNeg_countMat, human2PCs_zNeg_countMat)) 
humanPCs_zNeg_countMat %>% colSums() %>% hist(main = "Human Z- PCs", xlab = "Read counts per cell")

haoRhesus_zNeg_countMat <- getOrthologCountMat(haoRhesus_zNeg, 
                                               allPrimateOrthologs,
                                               "macaque", 
                                               normMethod = "shiftedLog", 
                                               useForL = 10127)
haoRhesus_zNeg_countMat %>% colSums() %>% hist(main = "Hao Rhesus Z- PCs", xlab = "Read counts per cell")

barteltMousePCs_zNeg_countMat <- getOrthologCountMat(barteltMousePCs_zNeg, 
                                                     allPrimateOrthologs, 
                                                     "mouse", 
                                                     normMethod = "shiftedLog", 
                                                     useForL = 10127)
barteltMousePCs_zNeg_countMat %>% colSums() %>% hist(main = "Bartelt Mouse Z- PCs", xlab = "Read counts per cell")

mousePCs_zNeg_countMat <- getOrthologCountMat(mousePCs_zNeg, 
                                              allPrimateOrthologs, 
                                              "mouse", 
                                              normMethod = "shiftedLog", 
                                              useForL = 10127)
mousePCs_zNeg_countMat %>% colSums() %>% hist(main = "Pilot Mouse Z- PCs", xlab = "Read counts per cell")

zNeg_countMat <- mergeCountMats(list(humanPCs_zNeg_countMat, 
                                     haoRhesus_zNeg_countMat, 
                                     barteltMousePCs_zNeg_countMat, 
                                     mousePCs_zNeg_countMat))
zNeg_countMat %>% colSums() %>% hist(main = "All Species Merged", xlab = "Read counts per cell", breaks = 100)

zNeg_dnr <- getDNR(zNeg_countMat)
hist(zNeg_dnr$dnr)

pseudobulked_zNeg_dnr <- getDNR(pseudobulk_zNeg)
hist(pseudobulked_zNeg_dnr$dnr)


# 4.0 Standard deviation by subtype across species ----
## 4.1 Z- ----
zNegSD <- getNonZeroSD(zNeg_countMat)
zNegSD[zNegSD != 0] %>% hist(breaks = 200, main = "Z- PC SD Across Species", xlab = "SD of gene expression in Z- PCs across species")
zNegSD[zNegSD != 0] %>% log() %>% hist(breaks = 200, main = "Z- PC SD Across Species, log transformed",  xlab = "log(SD) of gene expression in Z- PCs across species")
range(zNegSD, na.rm=T)
sort(zNegSD, decreasing = T) %>% head(n = 10)

## 4.2 Z+ ----
zPosSD <- getNonZeroSD(zPos_countMat)
zPosSD[zPosSD != 0] %>% hist(breaks = 200, main = "Z+ PC SD Across Species", xlab = "SD of gene expression in Z+ PCs across species")
zPosSD[zNegSD != 0] %>% log() %>% hist(breaks = 200, main = "Z+ PC SD Across Species, log transformed",  xlab = "log(SD) of gene expression in Z+ PCs across species")
range(zPosSD, na.rm=T)
sort(zPosSD, decreasing = T) %>% head(n = 10)


# 5.0 Mean normalized variance by subtype across species ----
## 5.1 Z- ----
zNegVar <- getNonZeroVar(zNeg_countMat)
zNegVar[zNegVar != 0] %>% 
  hist(breaks = 200, main = "Variance in Z- PCs Across Species", xlab = "variance in gene expression in Z- PCs across species")
range(zNegVar, na.rm = T)
sort(zNegVar, decreasing = T) %>% head(n = 10)
zNegVar_log <- zNegVar[zNegVar != 0] %>% log()
hist(zNegVar_log, breaks = 200, main = "Variance in Z- PCs Across Species", xlab = "log(variance) in gene expression in Z- PCs across species")
sort(zNegVar_log, decreasing = T) %>% head(n = 10)

## 5.2 Z+ ----
zPosVar <- getNonZeroVar(zPos_countMat)
zPosVar[zPosVar != 0] %>% hist(breaks = 200, main = "Variance in Z+ PCs Across Species", xlab = "variance in gene expression in Z+ PCs across species")
range(zPosVar, na.rm = T)
sort(zPosVar, decreasing = T) %>% head(n = 10)
zPosVar_log <- zPosVar[zPosVar != 0] %>% log()
hist(zPosVar_log, breaks = 200, main = "Variance in Z+ PCs Across Species", xlab = "log(variance) gene expression in Z- PCs across species")
sort(zPosVar_log, decreasing = T) %>% head(n = 10)

# Mean mean normalized variance across species ----
mergedVar_acrSpecies <- mergeVectors(list(zPosVar, zNegVar), list("zPos", "zNeg")) %>%
  na.omit()
meanVar_acrSpecies <- apply(mergedVar_acrSpecies, MARGIN = 1, mean)
hist(log(meanVar_acrSpecies))

# 6.0 Mean normalized variance by species across subtypes ----
## 6.1 Human ----
humanPCs_countMat <- mergeCountMats(list(humanPCs_zNeg_countMat, humanPCs_zPos_countMat))
humanVar <- getNonZeroVar(humanPCs_countMat)
range(humanVar, na.rm = T)
sort(humanVar, decreasing = T) %>% head(n = 10)
humanVar_log <- log(humanVar)
hist(humanVar_log, breaks = 100)
hist(humanVar, breaks = 100)

## 6.2 Rhesus ----
rhesusPCs_countMat <- mergeCountMats(list(haoRhesus_zNeg_countMat,haoRhesus_zPos_countMat))
rhesusVar <- getNonZeroVar(rhesusPCs_countMat)
range(rhesusVar, na.rm = T)
sort(rhesusVar, decreasing = T) %>% head(n = 10)
rhesusVar_log <- log(rhesusVar)
hist(rhesusVar_log, breaks = 100)
hist(rhesusVar, breaks = 100)

## 6.3 Mouse ----
mousePCs_countMat <- mergeCountMats(list(barteltMousePCs_zNeg_countMat, barteltMousePCs_zPos_countMat, mousePCs_zNeg_countMat, mousePCs_zPos_countMat))
mouseVar <- getNonZeroVar(mousePCs_countMat)
range(mouseVar, na.rm = T)
sort(mouseVar, decreasing = T) %>% head(n = 10)
mouseVar_log <- log(mouseVar)
hist(mouseVar_log, breaks = 100)
hist(mouseVar, breaks = 100)

## 6.4 Mean mean normalized variance across subtypes ----
mergedVar_acrSubtypes <- mergeVectors(list(mouseVar, humanVar, rhesusVar), list("mouse", "rhesus", "human")) %>%
  na.omit()
meanVar_acrSubtypes <- apply(mergedVar_acrSubtypes, MARGIN = 1, mean)


# 7.0 Explore average normalized variances across species and subtypes ----
## 7.1 Plot ----
avgVar <- mergeVectors(list(meanVar_acrSpecies, meanVar_acrSubtypes), list("AcrossSpecies", "AcrossSubtypes")) %>%
  filter(AcrossSpecies > 0) %>%
  filter(AcrossSubtypes > 00)
ggplot(data = avgVar, aes(x = AcrossSpecies, y = AcrossSubtypes)) + 
  geom_point() + 
  theme_minimal() +
  ggtitle("Average Variance in Gene Expression")
  
# Log transformed 
avgVar %>% log() %>% arrange(desc(AcrossSpecies)) %>% head(n = 20)
avgVar %>% log() %>% arrange(desc(AcrossSpecies)) %>% tail(n = 20)
ggplot(data = log(avgVar), aes(x = AcrossSpecies, y = AcrossSubtypes)) + 
  geom_point() + 
  theme_minimal() +
  ggtitle("Log Transformed Average Variance in Gene Expression")

## 7.2 GO terms -----
### Get entrez IDs for all genes ----
pcEntrez <- bitr(rownames(avgVar), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

### Low cross-species variance & low cross-subtype variance ----
lowSpec_lowSub <- avgVar %>% 
  log() %>%
  filter(AcrossSubtypes < -2.25) %>%
  filter(AcrossSpecies < -2)
lowSpec_lowSubEntrez <- bitr(rownames(lowSpec_lowSub), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

lowSpec_lowSub_bp <- enrichGO(as.character(lowSpec_lowSubEntrez$ENTREZID), 
                      OrgDb = org.Hs.eg.db,
                      ont = 'BP',
                      keyType = "ENTREZID", 
                      universe = as.character(pcEntrez$ENTREZID),
                      qvalueCutoff = 0.5, 
                      pvalueCutoff = 0.5)
clusterProfiler::dotplot(lowSpec_lowSub_bp)

lowSpec_lowSub_mf <- enrichGO(as.character(lowSpec_lowSubEntrez$ENTREZID), 
                           OrgDb = org.Hs.eg.db,
                           ont = 'MF',
                           keyType = "ENTREZID", 
                           universe = as.character(pcEntrez$ENTREZID),
                           qvalueCutoff = 0.5, 
                           pvalueCutoff = 0.5)
clusterProfiler::dotplot(lowSpec_lowSub_mf)

### Low cross-species variance, high cross-subtype variance ----
lowSpec_hiSub <- avgVar %>% 
  log() %>%
  filter(AcrossSubtypes > -2.25) %>%
  filter(AcrossSpecies < -2)
lowSpec_hiSubEntrez<- bitr(rownames(lowSpec_hiSub), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

lowSpec_hiSub_bp <- enrichGO(as.character(lowSpec_hiSubEntrez$ENTREZID), 
                           OrgDb = org.Hs.eg.db,
                           ont = 'BP',
                           keyType = "ENTREZID", 
                           universe = as.character(pcEntrez$ENTREZID),
                           qvalueCutoff = 0.5, 
                           pvalueCutoff = 0.5)
clusterProfiler::dotplot(lowSpec_hiSub_bp)

lowSpec_hiSub_mf <- enrichGO(as.character(lowSpec_hiSubEntrez$ENTREZID), 
                        OrgDb = org.Hs.eg.db,
                        ont = 'MF',
                        keyType = "ENTREZID", 
                        universe = as.character(pcEntrez$ENTREZID),
                        qvalueCutoff = 0.5, 
                        pvalueCutoff = 0.5)
clusterProfiler::dotplot(lowSpec_hiSub_mf)

### High cross-species variance, low cross-subtype variance ----
hiSpec_lowSub <- avgVar %>% 
  log() %>%
  filter(AcrossSubtypes < -2.25) %>%
  filter(AcrossSpecies > -2)
hiSpec_lowSubEntrez<- bitr(rownames(hiSpec_lowSub), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
hiSpec_lowSub_bp <- enrichGO(as.character(hiSpec_lowSubEntrez$ENTREZID), 
                        OrgDb = org.Hs.eg.db,
                        ont = 'BP',
                        keyType = "ENTREZID", 
                        universe = as.character(pcEntrez$ENTREZID),
                        qvalueCutoff = 0.5, 
                        pvalueCutoff = 0.5)
clusterProfiler::dotplot(hiSpec_lowSub_bp)

hiSpec_lowSub_mf <- enrichGO(as.character(hiSpec_lowSubEntrez$ENTREZID), 
                       OrgDb = org.Hs.eg.db,
                       ont = 'MF',
                       keyType = "ENTREZID", 
                       universe = as.character(pcEntrez$ENTREZID),
                       qvalueCutoff = 0.5, 
                       pvalueCutoff = 0.5)
clusterProfiler::dotplot(hiSpec_lowSub_mf)

### High cross-species variance & high cross-subtype variance ----
hiSpec_hiSub <- avg_normVars %>% 
  log() %>%
  filter(AcrossSubtypes > -2) %>%
  filter(AcrossSpecies > -5)
hiSpec_hiSubEntrez<- bitr(rownames(hiSpec_hiSub), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
hiSpec_hiSub_bp <- enrichGO(as.character(hiSpec_hiSubEntrez$ENTREZID), 
                       OrgDb = org.Hs.eg.db,
                       ont = 'BP',
                       keyType = "ENTREZID", 
                       universe = as.character(pcEntrez$ENTREZID),
                       qvalueCutoff = 0.5, 
                       pvalueCutoff = 0.5)
clusterProfiler::dotplot(hiSpec_hiSub_bp)

hiSpec_hiSub_mf <- enrichGO(as.character(hiSpec_hiSubEntrez$ENTREZID), 
                            OrgDb = org.Hs.eg.db,
                            ont = 'MF',
                            keyType = "ENTREZID", 
                            universe = as.character(pcEntrez$ENTREZID),
                            qvalueCutoff = 0.5, 
                            pvalueCutoff = 0.5)
clusterProfiler::dotplot(hiSpec_hiSub_mf)
