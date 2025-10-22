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
rm(barteltMousePCs)

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


# 2.0 Dynamic range across species ----
## 2.1 All PCs ----
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

### 2.1.2 Merge count mats for all species for all PCs (excluding Hao Marmoset because it's weird)
mergedMats <- mergeCountMats(list(countMat_pilotHuman1_pcs, 
                                  countMat_pilotHuman2_pcs, 
                                  countMat_pilotRhesus_pc, 
                                  countMat_pilotMouse_pc, 
                                  countMat_barteltMouse_pc, 
                                  countMat_haoRhesus_pc))
### 2.1.3 Get DNR for all PCs across all species
dnr_AllPCs_AcrossSpecies <- getDNR(mergedMats)
range(dnr_AllPCs_AcrossSpecies$dnr, na.rm = T)
hist(dnr_AllPCs_AcrossSpecies$dnr)

## 2.2 Z+ PCs


## 2.3 Z- PCs


# 3.0 Dynamic range across species 

