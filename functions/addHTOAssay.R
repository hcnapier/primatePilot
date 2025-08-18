# addHTOAssay
# Hailey Napier
# August 4, 2025

# Function to add a normalized HTO assay to a Seurat object with GEX data

require(Seurat)

addHTOAssay <- function(obj, htoMat){
  #htoMat <- htoMat + 1
  
  obj[["HTO"]] <- CreateAssayObject(counts = htoMat) # add in HTO library
  obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR") # normalize HTO library using centered log-ratio transformation
  return (obj)
}