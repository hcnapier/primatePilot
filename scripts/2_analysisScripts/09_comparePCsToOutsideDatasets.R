# 09 compareToOutsideDatasets.R
# Compare Primate Pilot dataset transcriptional similarity to other datasets 
# October 2, 2025
# Hailey Napier

# 0.0 Setup ----
## 0.1 Load packages ----
require(Seurat)
require(dplyr)
require(ggplot2)

## 0.2 Setup one to one orthologs ----
setwd("/Users/haileynapier/Work/VertGenLab/Projects/zebrinEvolution/Data/geneLists/orthologs")
allPrimateOrthologs <- read.delim("human_MouseMarmosetMacaque_orthologs.txt", sep = "\t", header = T)
names(allPrimateOrthologs) <- tolower(names(allPrimateOrthologs))
allPrimateOrthologs %>%
  filter(macaque.homology.type == "ortholog_one2one") %>%
  filter(mouse.homology.type == "ortholog_one2one") %>%
  filter(white.tufted.ear.marmoset.homology.type == "ortholog_one2one") %>%
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

## 0.4 Setup Hao rhesus data ----
# load data 
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/haoReanalysis/data/seuratObjs")
haoRhesus <- readRDS("macaqueSingleCellPC.rds")
# aggregate expression
pseudobulk_haoRhesus <- AggregateExpression(haoRhesus, features = pull(allPrimateOrthologs["macaque.gene.name"]), group.by = "doublet_info", return.seurat = F, assays = "RNA") #group by doublet_info to lump all cells together
# convert to dense matrix 
pseudobulk_haoRhesus <- pseudobulk_haoRhesus$RNA %>% as.matrix()
# remove zero counts
pseudobulk_haoRhesus <- pseudobulk_haoRhesus[rowSums(pseudobulk_haoRhesus)>0, ] 
rm(haoRhesus)


## 0.5 Setup Hao marmoset data ----
# load data 
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/haoReanalysis/data/seuratObjs")
haoMarmoset <- readRDS("marmosetSingleCellPC.rds")
# aggregate expression
pseudobulk_haoMarmoset <- AggregateExpression(haoMarmoset, features = pull(allPrimateOrthologs["white.tufted.ear.marmoset.gene.name"]), group.by = "doublet_info", return.seurat = F, assays = "RNA") #group by doublet_info to lump all cells together
# convert to dense matrix 
pseudobulk_haoMarmoset <- pseudobulk_haoMarmoset$RNA %>% as.matrix()
# remove zero counts
pseudobulk_haoMarmoset <- pseudobulk_haoMarmoset[rowSums(pseudobulk_haoMarmoset)>0, ] 
rm(haoMarmoset)

## 0.6 Setup Bartelt mouse data ----
# load data 
barteltMouse <- readRDS("~/Work/VertGenLab/Projects/zebrinEvolution/Data/sequencingData/LukesData/SCA7.8wk.cleaned.rds")
barteltMousePCs <- subset(barteltMouse, idents = "Purkinje cells")
Idents(barteltMousePCs) <- "Type"
barteltMousePCs <- subset(barteltMousePCs, idents = "WT")
rm(barteltMouse)
# aggregate expression
pseudobulk_barteltMousePCs <- AggregateExpression(barteltMousePCs, features = pull(allPrimateOrthologs["mouse.gene.name"]), group.by = "Type", return.seurat = F, assays = "RNA") #group by Type to lump all PCs together
# convert to dense matrix
pseudobulk_barteltMousePCs <- pseudobulk_barteltMousePCs$RNA %>% as.matrix()
# remove zero counts
pseudobulk_barteltMousePCs <- pseudobulk_barteltMousePCs[rowSums(pseudobulk_barteltMousePCs)>0, ] 

## 0.6 Merge all into one matrix ----
as.data.frame(pseudobulk_barteltMousePCs)
as.data.frame(pseudobulk_haoMarmoset)
