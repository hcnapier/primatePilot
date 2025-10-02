# 09 compareToOutsideDatasets.R
# Compare Primate Pilot dataset transcriptional similarity to other datasets 
# October 2, 2025
# Hailey Napier

# 0.0 Setup ----
## 0.1 Load packages ----
require(Seurat)
require(dplyr)
require(ggplot2)

## 0.2 Setup primate pilot data ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/")
pseudobulkMerged <- readRDS("pseudobulkMerged.rds")
pseudobulkMerged %>% 
  as.data.frame() %>%
  select(contains("Purkinje")) -> pseudobulkMerged_pcs
names(pseudobulkMerged_pcs) <- paste(names(pseudobulkMerged_pcs), "Napier", sep = "_")
pseudobulkMerged_pcs <- pseudobulkMerged_pcs[rowSums(pseudobulkMerged_pcs[, -1])>0, ]
geneList <- pseudobulkMerged_pcs %>% rownames()
rm(pseudobulkMerged)

## 0.3 Setup Hao rhesus data ----
# load data 
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/haoReanalysis/data/seuratObjs")
haoRhesus <- readRDS("macaqueSingleCellPC.rds")
# aggregate expression
pseudobulk_haoRhesus <- AggregateExpression(haoRhesus, features = geneList, group.by = "doublet_info", return.seurat = F, assays = "RNA") #group by doublet_info to lump all cells together
# convert to dense matrix 
pseudobulk_haoRhesus <- pseudobulk_haoRhesus$RNA %>% as.matrix()
pseudobulk_haoRhesus <- pseudobulk_haoRhesus[rowSums(pseudobulk_haoRhesus)>0, ] 
rm(haoRhesus)

## 0.4 Setup Hao marmoset data ----
# load data 
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/haoReanalysis/data/seuratObjs")
haoMarmoset <- readRDS("marmosetSingleCellPC.rds")
# aggregate expression
pseudobulk_haoMarmoset <- AggregateExpression(haoMarmoset, features = geneList, group.by = "doublet_info", return.seurat = F, assays = "RNA") #group by doublet_info to lump all cells together
# convert to dense matrix 
pseudobulk_haoMarmoset <- pseudobulk_haoMarmoset$RNA %>% as.matrix()
pseudobulk_haoMarmoset <- pseudobulk_haoMarmoset[rowSums(pseudobulk_haoMarmoset)>0, ] 
rm(haoMarmoset)

## 0.5 Setup Bartelt mouse data ----
# load data 

