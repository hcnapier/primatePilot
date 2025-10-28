# TPMNormalize 
# Function to TPM normalize a count matrix given a species-specific dataframe containing gene sizes
# October 28, 2025

require(tibble)
require(dplyr)

TPMNormalize <- function(countMat, geneSizes){
  if(!is.matrix(countMat)){
    stop("countMat must by of type matrix")
  }
  geneNameCol <- colnames(geneSizes)[which(str_detect(colnames(geneSizes), "gene.name"))]
  geneSizeCol <- colnames(geneSizes)[which(str_detect(colnames(geneSizes), "gene.size"))]
  genes <- rownames(countMat) %>% as.data.frame()
  names(genes) <- geneNameCol
  genes <- inner_join(genes, geneSizes)
  geneSizeMat <- countMat
  for(currGene in genes[geneNameCol]){
    geneSizeMat[currGene,] <- genes[geneSizeCol][which(genes[geneNameCol] == currGene),]
  }
  rpkMat <- countMat/geneSizeMat
  scaleFactor <- sum(rpkMat, na.rm = T)/1000
  tpmMat <- rpkMat/scaleFactor
  return(tpmMat)
}

