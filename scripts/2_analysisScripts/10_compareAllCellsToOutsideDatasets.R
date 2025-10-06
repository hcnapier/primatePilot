# 10_compareToOutsideDatasets.R
# Compare Primate Pilot dataset transcriptional similarity to other datasets 
# October 3, 2025
# Hailey Napier

# RUN ON DCC

# 0.0 Setup ----
## 0.1 Load packages ----
require(Seurat)
require(dplyr)
require(ggplot2)
require(reshape2)

## 0.2 Source functions ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/functions")
source("pseudobulkHumanOrtholog.R")

## 0.3 Setup one to one orthologs ----
setwd("/Users/haileynapier/Work/VertGenLab/Projects/zebrinEvolution/Data/geneLists/orthologs")
allPrimateOrthologs <- read.delim("human_MouseMarmosetMacaque_orthologs.txt", sep = "\t", header = T)
names(allPrimateOrthologs) <- tolower(names(allPrimateOrthologs))
allPrimateOrthologs %>%
  filter(macaque.homology.type == "ortholog_one2one") %>%
  filter(mouse.homology.type == "ortholog_one2one") %>%
  filter(white.tufted.ear.marmoset.homology.type == "ortholog_one2one") %>%
  distinct() -> allPrimateOrthologs

## 0.4 Setup primate pilot data ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/")
pseudobulkMerged <- readRDS("pseudobulkMerged.rds")
pseudobulkMerged %>% 
  as.data.frame() %>%
  select(contains("Purkinje")) -> pseudobulkMerged_pcs
names(pseudobulkMerged_pcs) <- paste(names(pseudobulkMerged_pcs), "Napier", sep = "_")
pseudobulkMerged_pcs <- pseudobulkMerged_pcs[rowSums(pseudobulkMerged_pcs[, -1])>0, ]
rm(pseudobulkMerged)
pseudobulkMerged_pcs$gene.name <- rownames(pseudobulkMerged_pcs)

## 0.5 Setup Hao rhesus data ----
pseudobulkHumanOrtholog()

