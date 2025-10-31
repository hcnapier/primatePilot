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
allPrimateOrthologs$gene.size <- allPrimateOrthologs$gene.end..bp. - allPrimateOrthologs$gene.start..bp.
allPrimateOrthologs$macaque.gene.size <- allPrimateOrthologs$macaque.chromosome.scaffold.end..bp. - allPrimateOrthologs$macaque.chromosome.scaffold.start..bp.
allPrimateOrthologs$mouse.gene.size <- allPrimateOrthologs$mouse.chromosome.scaffold.end..bp. - allPrimateOrthologs$mouse.chromosome.scaffold.start..bp.
allPrimateOrthologs$white.tufted.ear.marmoset.gene.size <- allPrimateOrthologs$white.tufted.ear.marmoset.chromosome.scaffold.end..bp. - allPrimateOrthologs$white.tufted.ear.marmoset.chromosome.scaffold.start..bp.

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
                                             "macaque", 
                                             normMethod = "shiftedLog", 
                                             useForL = 7208)
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
# pseudobulk
pseudobulk_barteltMousePCs <- pseudobulkHumanOrtholog(barteltMousePCs, 
                                                      groupByIdent = "Type", 
                                                      orthologDF = allPrimateOrthologs, 
                                                      species = "mouse", 
                                                      datasetID = "Purkinje_mouse_Bartelt")
# get count matrix in terms of human orthologs
countMat_barteltMouse_pc <- getOrthologCountMat(barteltMousePCs,
                                                allPrimateOrthologs, 
                                                "mouse", 
                                                normMethod = "shiftedLog", 
                                                useForL = 7208)
countMat_barteltMouse_pc %>% colSums() %>% hist(main = "Bartelt Mouse", xlab = "Read counts per cell", breaks = 14)
countMat_barteltMouse_pc %>% colSums() %>% mean()

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

# Pseudobulked
#pseudobulkAllPCs_dnr <- getDNR(pseudobulkMerged_pcs)
#hist(pseudobulkAllPCs_dnr$dnr)


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
  subset(idents = c(1))
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
  subset(idents = c(2))
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
human1PCs_zPos_countMat <- humanPCs_zPos@assays$RNA$counts.1 %>% as.matrix()
human2PCs_zPos_countMat <- humanPCs_zPos@assays$RNA$counts.2 %>% as.matrix()
humanPCs_zPos_countMat <- mergeCountMats(list(human1PCs_zPos_countMat, human2PCs_zPos_countMat)) 
humanPCs_zPos_countMat %>% colSums() %>% hist(main = "Human Z+ PCs", xlab = "Read counts per cell")
haoRhesus_zPos_countMat <- getOrthologCountMat(haoRhesus_zPos, 
                                               allPrimateOrthologs,
                                               "macaque", 
                                               normMethod = "")
haoRhesus_zPos_countMat %>% colSums() %>% hist(main = "Hao Rhesus Z+ PCs", xlab = "Read counts per cell")
barteltMousePCs_zPos_countMat <- getOrthologCountMat(barteltMousePCs_zPos, 
                                            allPrimateOrthologs, 
                                            "mouse", 
                                            normMethod = "")
barteltMousePCs_zPos_countMat %>% colSums() %>% hist(main = "Bartelt Mouse Z+ PCs", xlab = "Read counts per cell")
mousePCs_zPos_countMat <- getOrthologCountMat(mousePCs_zPos, 
                                              allPrimateOrthologs, 
                                              "mouse", 
                                              normMethod = "")
mousePCs_zPos_countMat %>% colSums() %>% hist(main = "Pilot Mouse Z+ PCs", xlab = "Read counts per cell")
zPos_countMat <- mergeCountMats(list(humanPCs_zPos_countMat, 
                                     haoRhesus_zPos_countMat, 
                                     barteltMousePCs_zPos_countMat, 
                                     mousePCs_zPos_countMat)) %>% 
  shiftedLogNorm()
zPos_countMat %>% colSums() %>% hist(main = "All Species Merged", xlab = "Read counts per cell", breaks = 100)
zPos_dnr <- getDNR(zPos_countMat)
zPos_dnr %>% filter(dnr > 0) -> zPos_dnr
hist(zPos_dnr$dnr, breaks = 100)
median(zPos_dnr$dnr)
highZPosDNR <- zPos_dnr %>%
  filter(dnr > median(zPos_dnr$dnr)) %>%
  dplyr::select(gene)
zPosEntrez <- bitr(rownames(zPos_countMat), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
highZPosDNREntrez <- bitr(highZPosDNR$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
highZPosDNR_bp <- enrichGO(as.character(highZPosDNREntrez$ENTREZID), 
                             OrgDb = org.Hs.eg.db,
                             ont = 'BP',
                             keyType = "ENTREZID", 
                             universe = as.character(pcEntrez$ENTREZID),
                             qvalueCutoff = 0.5, 
                             pvalueCutoff = 0.5)
clusterProfiler::dotplot(lowSpec_hiSub_bp)

#pseudobulked_zPos_dnr <- getDNR(pseudobulk_zPos)
#hist(pseudobulked_zPos_dnr$dnr)

## 3.4 DNR Z- ----
human1PCs_zNeg_countMat <- humanPCs_zNeg@assays$RNA$counts.1 %>% as.matrix() 
human2PCs_zNeg_countMat <- humanPCs_zNeg@assays$RNA$counts.2 %>% as.matrix()
humanPCs_zNeg_countMat <- mergeCountMats(list(human1PCs_zNeg_countMat, human2PCs_zNeg_countMat)) 
humanPCs_zNeg_countMat %>% colSums() %>% hist(main = "Human Z- PCs", xlab = "Read counts per cell")

haoRhesus_zNeg_countMat <- getOrthologCountMat(haoRhesus_zNeg, 
                                               allPrimateOrthologs,
                                               "macaque", 
                                               normMethod = "")
haoRhesus_zNeg_countMat %>% colSums() %>% hist(main = "Hao Rhesus Z- PCs", xlab = "Read counts per cell")

barteltMousePCs_zNeg_countMat <- getOrthologCountMat(barteltMousePCs_zNeg, 
                                                     allPrimateOrthologs, 
                                                     "mouse", 
                                                     normMethod = "")
barteltMousePCs_zNeg_countMat %>% colSums() %>% hist(main = "Bartelt Mouse Z- PCs", xlab = "Read counts per cell")

mousePCs_zNeg_countMat <- getOrthologCountMat(mousePCs_zNeg, 
                                              allPrimateOrthologs, 
                                              "mouse", 
                                              normMethod = "")
mousePCs_zNeg_countMat %>% colSums() %>% hist(main = "Pilot Mouse Z- PCs", xlab = "Read counts per cell")

zNeg_countMat <- mergeCountMats(list(humanPCs_zNeg_countMat, 
                                     haoRhesus_zNeg_countMat, 
                                     barteltMousePCs_zNeg_countMat, 
                                     mousePCs_zNeg_countMat)) %>%
  shiftedLogNorm(useForL = "meanReadDepth")
zNeg_countMat %>% colSums() %>% hist(main = "Z- All Species Merged", xlab = "Read counts per cell", breaks = 100)

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
mergedVar_acrSubtypes <- mergeVectors(list(mouseVar, humanVar, rhesusVar), list("mouse", "human", "rhesus")) %>%
  na.omit()
meanVar_acrSubtypes <- apply(mergedVar_acrSubtypes, MARGIN = 1, mean)


# 7.0 Explore average normalized variances across species and subtypes ----
## 7.1 Plot ----
avgVar <- mergeVectors(list(meanVar_acrSpecies, meanVar_acrSubtypes), list("AcrossSpecies", "AcrossSubtypes")) %>%
  filter(AcrossSpecies > 0) %>%
  filter(AcrossSubtypes > 0)
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
  filter(AcrossSubtypes < -1.75) %>%
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
  filter(AcrossSubtypes > -1.75) %>%
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
  filter(AcrossSubtypes < -1.75) %>%
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
hiSpec_hiSub <- avgVar %>% 
  log() %>%
  filter(AcrossSubtypes > -1.75) %>%
  filter(AcrossSpecies > -2)
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


# 8.0 Explore count matrices further ----
zNeg_countMat %>% hist()
zPos_countMat %>% hist()
zNeg_countMat %>% rowSums() %>% head(n = 20)
zNegVar %>% sort(decreasing = T) %>% head(n = 20)
zNegVar %>% sort(decreasing = T) %>% tail(n = 20)
zPos_countMat %>% rowSums() %>% head(n = 20)

# 9.0 Linear model approach -----
## 9.1 Make sure each count matrix is normalized ----
humanPCs_zPos_countMat_norm <- humanPCs_zPos_countMat %>% shiftedLogNorm()
humanPCs_zNeg_countMat_norm <- humanPCs_zNeg_countMat %>% shiftedLogNorm()
haoRhesus_zNeg_countMat_norm <- haoRhesus_zNeg_countMat %>% shiftedLogNorm()
haoRhesus_zPos_countMat_norm <- haoRhesus_zPos_countMat %>% shiftedLogNorm()
mergedMousePCs_zNeg_countMat <- mergeCountMats(list(mousePCs_zNeg_countMat, 
                                                    barteltMousePCs_zNeg_countMat))
mergedMousePCs_zNeg_countMat_norm <- mergedMousePCs_zNeg_countMat %>% shiftedLogNorm()
mergedMousePCs_zPos_countMat <- mergeCountMats(list(mousePCs_zPos_countMat, 
                                                    barteltMousePCs_zPos_countMat))
mergedMousePCs_zPos_countMat_norm <- mergedMousePCs_zPos_countMat %>% shiftedLogNorm()

## 9.2 Find variance for each species/subtype combination ----
indvCountMatList <- list(humanPCs_zPos_countMat_norm, 
                         humanPCs_zNeg_countMat_norm, 
                         haoRhesus_zPos_countMat_norm, 
                         haoRhesus_zNeg_countMat_norm, 
                         mergedMousePCs_zPos_countMat_norm, 
                         mergedMousePCs_zPos_countMat_norm)
names(indvCountMatList) <- c("human_zPos", 
                             "human_zNeg", 
                             "rhesus_zPos", 
                             "rhesus_zNeg", 
                             "mouse_zPos", 
                             "mouse_zNeg")

indvVarDF <- data.frame("variance" = numeric(), 
                        "species" = character(), 
                        "subtype" = character(), 
                        "gene" = character())

for(currMatName in names(indvCountMatList)){
  if(str_detect(currMatName, "human")){
    species = "human"
  }else if(str_detect(currMatName, "rhesus")){
    species = "rhesus"
  }else if(str_detect(currMatName, "mouse")){
    species = "mouse"
  }
  
  if(str_detect(currMatName, "zPos")){
    subtype = "zPos"
  }else if(str_detect(currMatName, "zNeg")){
    subtype = "zNeg"
  }
  
  print("Assigned species and subtype")
  
  currMat <- indvCountMatList[[currMatName]]
  
  tmpDF <- getNonZeroVar(currMat) %>% data.frame()
  
  print("Calculated variance")
  
  names(tmpDF) <- 'variance'
  tmpDF$gene <- rownames(tmpDF)
  tmpDF$species <- species
  tmpDF$subtype <- subtype
  
  indvVarDF <- rbind(tmpDF, indvVarDF)
}

indvVarDF <- indvVarDF %>% na.omit()

## 9.3 Linear model ----
varBySpecSub_lm <- lm(variance ~ species + subtype, data = indvVarDF)
summary(varBySpecSub_lm)
