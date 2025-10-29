# shiftedLogNorm
# Function to normalize a count matrix using the shifted log method 
# See https://www.sc-best-practices.org/preprocessing_visualization/normalization.html#shifted-logarithm
# October 29, 2025

require(tibble)
require(dplyr)

source("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/functions/removeZeroRowsCols.R")

shiftedLogNorm <- function(countMat, pseudocount = 1, useForL = "meanReadDepth"){
  if(!is.matrix(countMat)){
    stop("countMat must by of type matrix")
  }
  s <- countMat
  s[] <- 0
  for(currCell in colnames(countMat)){ # Generate a matrix containing the size factor numerator (read depth per cell)
    s[,currCell] <- countMat[,currCell] %>% sum()
  }
  if(useForL == "meanReadDepth"){
    L = countMat %>% colSums() %>% mean()
  }else if(is.numeric(useForL)){
    L = useForL
  }else{
    stop("Please input a valid option to useForL")
  }
  s = s/L
  scaledMat = log((countMat/s) + pseudocount)
  meanRawReadDepth = countMat %>% colSums() %>% mean()
  print(paste("mean raw read depth = ", meanRawReadDepth))
  return(scaledMat)
}