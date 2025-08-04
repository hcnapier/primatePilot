# 00_readData.R
# Read & format preprocessed sequencing files
# Hailey Napier
# August 4, 2025

# 0.0 Setup ----
## 0.1 Load packages ----
library(Seurat)

# 1.0 Read in count matrices ----
speciesNames <- c("human", "rhesus", "mouse")
speciesCountMatList_hto <- list()
speciesCountMatList_gex <- list()
setwd("/Users/haileynapier/Work/VertGenLab/Projects/zebrinEvolution/Data/20250423_PrimatePilot/cellrangerOuts")
for(currSpecies in speciesNames){
  filename <- paste(currSpecies, "filtered_feature_bc_matrix.h5", sep = "_")
  print(filename)
  tmp <- Read10X_h5(filename)
  speciesCountMatList_gex[[currSpecies]] <- tmp$`Gene Expression`
  speciesCountMatList_hto[[currSpecies]] <- tmp$`Antibody Capture`
  rm(tmp)
}

# 2.0 Filter hto and gex matrices for cells that are present in both ----
for(currSpecies in speciesNames){
  joint <- intersect(colnames(speciesCountMatList_gex[[currSpecies]]), colnames(speciesCountMatList_hto[[currSpecies]]))
  speciesCountMatList_gex[[currSpecies]] <- speciesCountMatList_gex[[currSpecies]][ ,joint]
  speciesCountMatList_hto[[currSpecies]] <- speciesCountMatList_hto[[currSpecies]][ ,joint]
}

# 3.0 Save filtered count matrices ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/countMats")
saveRDS(speciesCountMatList_gex, "allSpeciesCountMatList_GEX.rds")
saveRDS(speciesCountMatList_hto, "allSpeciesCountMatList_HTO.rds")
