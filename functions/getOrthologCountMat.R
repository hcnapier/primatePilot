# getOrthologCountMat.R
# Hailey Napier
# October 6, 2025

# Function to get a count matrix with gene names converted to human ortholog from a Seurat object

require(Seurat)
require(dplyr)

getOrthologCountMat <- function(obj, orthologDF, species){
  countMat <- obj@assays$RNA$counts
  countMat %>% as.data.frame() -> countMat
  cellBxs <- colnames(countMat)
  speciesGeneName <- paste(species, "gene.name", sep = ".")
  countMat[speciesGeneName] <- rownames(countMat)
  countMat <- inner_join(orthologDF, countMat)
  countMat %>%
    select(all_of(cellBxs), gene.name) -> countMat
  rownames(countMat) <- countMat$gene.name
  countMat$gene.name <- NULL
  countMat %>% as.matrix() -> countMat
  return(countMat)
}