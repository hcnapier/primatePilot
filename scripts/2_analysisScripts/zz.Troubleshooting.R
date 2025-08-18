speciesNames <- c("human", "rhesus", "mouse")
rawList_hto <- list()
setwd("/Users/haileynapier/Work/VertGenLab/Projects/zebrinEvolution/Data/20250423_PrimatePilot/cellrangerOuts")
for(currSpecies in speciesNames){
  filename <- paste(currSpecies, "raw_feature_bc_matrix.h5", sep = "_")
  print(filename)
  tmp <- Read10X_h5(filename)
  #rawList[[currSpecies]] <- tmp$`Gene Expression`
  rawList_hto[[currSpecies]] <- tmp$`Antibody Capture`
  rm(tmp)
}

rawList_hto[["human"]] %>% sum()
rawList_hto[["mouse"]] %>% sum()
rawList_hto[["rhesus"]] %>% sum()

rawList_hto[["human"]] %>% removeZeroRowsCols()
rawList_hto[["rhesus"]] %>% removeZeroRowsCols()
