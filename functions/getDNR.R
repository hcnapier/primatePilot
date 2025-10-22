# getDNR
# Function to calculate dynamic range (DNR) of each gene in a count matrix
# October 22, 2025

# DNR = log10(maxExp/minExp)
# Input: count matrix where rows are genes and columns are cells. Column row names must be gene names.
# Output: DNR dataframe with a column for gene names and a column for DNR value

source("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/functions/removeZeroRowsCols.R")
require(dplyr)

getDNR <- function(countMat){
  # remove genes with only zero counts 
  countMat <- countMat %>% removeZeroRowsCols()
  for(currGene in countMat %>% rownames()){
    max <- countMat[currGene,] %>% max()
    min <- countMat[currGene,] %>% min()
    if(min == 0){ # ignore zeros
      min = unique(sort(countMat[currGene,]))[[2]]
    }
    dnr <- log10(min/max)
  }
}
