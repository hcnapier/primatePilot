# getOrthologCountMat.R
# Hailey Napier
# October 6, 2025

# Function to get a count matrix with gene names converted to human ortholog from a Seurat object

require(Seurat)
require(dplyr)

source("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/functions/TPMNormalize.R")

getOrthologCountMat <- function(obj, orthologDF, species, tpm = FALSE){
  countMat <- obj@assays$RNA$counts # Extract raw count matrix
  speciesGeneName <- paste(species, "gene.name", sep = ".")
  if(tpm){
    countMat <- countMat %>% as.matrix()
    speciesColumns <- colnames(orthologDF)[which(str_detect(colnames(orthologDF), species))]
    geneSizeDF <- orthologDF %>%
      dplyr::select(all_of(speciesColumns))
    countMat <- TPMNormalize(countMat, geneSizeDF)
  }
  countMat %>% as.data.frame() -> countMat
  cellBxs <- colnames(countMat)
  countMat[speciesGeneName] <- rownames(countMat)
  countMat <- inner_join(orthologDF, countMat)
  countMat %>%
    dplyr::select(all_of(cellBxs), gene.name) -> countMat
  rownames(countMat) <- countMat$gene.name
  countMat$gene.name <- NULL
  countMat %>% as.matrix() -> countMat
  return(countMat)
}