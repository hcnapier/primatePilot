# shiftedLogNorm
# Function to normalize a count matrix using the shifted log method 
# See https://www.sc-best-practices.org/preprocessing_visualization/normalization.html#shifted-logarithm
# October 29, 2025

require(tibble)
require(dplyr)

source("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/functions/removeZeroRowsCols.R")

shiftedLogNorm <- function(countMat, pseudocount = 1, useForL = "meanReadDepth", removeLowCounts = F, lowCountThreshold){
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
    stop("Please input a numeric option to useForL")
  }
  s = s/L
  if(removeLowCounts){
    if(is.numeric(lowCountThreshold)){
      tmpScaledMat <- countMat/s
      tmpScaledMat <- tmpScaledMat[which(rowMeans(tmpScaledMat) > lowCountThreshold),]
      scaledMat = log(tmpScaledMat + pseudocount)
    }else{
      stop("Please input a numeric option to lowCountThreshold")
    }
  }else{
    scaledMat = log((countMat/s) + pseudocount)
  }
  return(scaledMat)
}