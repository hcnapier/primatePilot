removeZeroRowsCols <- function(mat){
  result <- mat[rowSums(mat[, -1])>0, ]
  result <- result[, colSums(result != 0) > 0]
  return(result)
}