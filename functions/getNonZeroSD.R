# getNonZeroSD
# Function to get the SD for every row in a matrix, ignoring zeros
# October 27, 2025

# Input is a count matrix where each row is a gene and each column is a cell

source("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/functions/removeZeroRowsCols.R")

getNonZeroSD <- function(countMat){
  # First remove any columns or rows with only zero counts
  countMat <- countMat %>% removeZeroRowsCols()
  # Next replace remaining zeros with NA
  countMat[countMat == 0] <- NA
  # Compute standard deviation for each row
  sdOut <- apply(countMat, MARGIN = 1, sd, na.rm = T)
  return(sdOut)
}

