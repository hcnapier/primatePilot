# getNonZeroVar
# Function to get the variance for every row in a matrix, ignoring zeros
# October 20, 2025

# Input is a count matrix where each row is a gene and each column is a cell

source("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/functions/removeZeroRowsCols.R")

getNonZeroVar <- function(countMat){
  # First remove any columns or rows with only zero counts
  countMat <- countMat %>% removeZeroRowsCols()
  # Next replace remaining zeros with NA
  countMat[countMat == 0] <- NA
  # Compute standard deviation for each row
  varOut <- apply(countMat, MARGIN = 1, var, na.rm = T)
  return(varOut)
}

