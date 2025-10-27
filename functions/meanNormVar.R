# meanNormVar
# Function to get the mean normalized variance for every row in a matrix, excluding zeros
# October 27, 2025

source("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/functions/removeZeroRowsCols.R")

meanNormVar <- function(countMat){
  # First remove any rows or columns with only zero counts
  countMat <- countMat %>% removeZeroRowsCols()
  # Next replace remaining zeros with NA
  countMat[countMat == 0] <- NA
  # Compute variance for each row
  countVar <- apply(countMat, MARGIN = 1, var, na.rm = T)
  # Compute mean for each row
  countMean <- apply(countMat, MARGIN = 1, mean, na.rm = T)
  # Compute mean normalized variance
  meanNormVar <- countVar/countMean
  return(meanNormVar)
}