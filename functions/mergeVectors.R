# mergeVectors
# Function to merge multiple vectors into one dataframe. Inner join based on vector element names. 
# October 28, 2025

require(tibble)
mergeVectors <- function(vectorList, namesList){
  if(is.list(vectorList) != T){
    stop("Input must be a list of count matrix objects.")
    }
  outDF <- vectorList[1] %>% as.data.frame() %>% rownames_to_column()
  names(outDF)[2] <- namesList[1]
  for(i in 2:length(vectorList)){
    tmp <- vectorList[i] %>% as.data.frame() %>% rownames_to_column()
    names(tmp)[2] <- namesList[i]
    outDF <- inner_join(outDF, tmp, by = "rowname")
  }
  rownames(outDF) <- outDF$rowname
  outDF$rowname <- NULL
  return(outDF)
}
