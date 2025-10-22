# mergeCountMats
# Function to merge count matrices for variance analysis
# Assumes each count matrix is organized such that each row is a gene name and each column is a cell

# Input: list of count matrices to merge

# Hailey Napier
# October 22, 2025

require(dplyr)
require(tibble)

mergeCountMats <- function(countMats){
  if(is.list(countMats) != T){
    stop("Input must be a list of count matrix objects.")
  }
  outMat <- countMats[[1]] %>% as.data.frame() %>% rownames_to_column()
  for(listIndex in 2:length(countMats)){
    tmp <- countMats[[listIndex]] %>% as.data.frame() %>% rownames_to_column()
    outMat <- inner_join(outMat, tmp, by = "rowname")
  }
  rownames(outMat) <- outMat$rowname
  outMat$rowname <- NULL
  outMat <- outMat %>% as.matrix()
  return(outMat)
}

