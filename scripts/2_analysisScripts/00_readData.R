# 00_readData.R
# Read & format preprocessed sequencing files
# Hailey Napier
# August 4, 2025

# 0.0 Setup ----
## 0.1 Load packages ----
library(Seurat)
## 0.2 Source functions ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/functions")
source("removeZeroRowsCols.R")

# 1.0 Read in count matrices ----
speciesNames <- c("human", "rhesus", "mouse")
speciesCountMatList_hto <- list()
speciesCountMatList_gex <- list()
setwd("/Users/haileynapier/Work/VertGenLab/Projects/zebrinEvolution/Data/sequencingData/20250423_PrimatePilot/pipseekerOuts/filteredMatrix")
for(currSpecies in speciesNames){
  tmp <- Read10X(currSpecies)
  speciesCountMatList_gex[[currSpecies]] <- tmp$`Gene Expression`
  speciesCountMatList_hto[[currSpecies]] <- tmp$`Cell Hashing`
  rm(tmp)
}

# 2.0 Filter hto and gex matrices for cells that are present in both ----
for(currSpecies in speciesNames){
  joint <- intersect(colnames(speciesCountMatList_gex[[currSpecies]]), colnames(speciesCountMatList_hto[[currSpecies]]))
  speciesCountMatList_gex[[currSpecies]] <- speciesCountMatList_gex[[currSpecies]][ ,joint]
  speciesCountMatList_hto[[currSpecies]] <- speciesCountMatList_hto[[currSpecies]][ ,joint]
}

# 3.0 Filter out unused HTOs ----
speciesCountMatList_htoFilt <- lapply(speciesCountMatList_hto, removeZeroRowsCols)

# 4.0 Save filtered count matrices ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/countMats")
saveRDS(speciesCountMatList_gex, "allSpeciesCountMatList_GEX.rds")
saveRDS(speciesCountMatList_htoFilt, "allSpeciesCountMatList_HTO.rds")
