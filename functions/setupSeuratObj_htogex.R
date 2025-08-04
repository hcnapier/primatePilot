# setupSeuratObj_htogex
# Hailey Napier
# August 4, 2025

# Function to set up a normalized seurat object for a pooled sequencing library including an
# endogenous gene expression library and a HTO library.

require(Seurat)

setupSeuratObj_htogex <- function(gexMat, htoMat){
  obj <- CreateSeuratObject(counts = gexMat)
  #print("GEX Obj Created")
  obj <- NormalizeData(obj) # normalize gex library
  #print("GEX Obj Normalized")
  obj <- FindVariableFeatures(obj, selection.method = "mean.var.plot") # find variable features in gex library
  #print("GEX Var Features Found")
  obj <- ScaleData(obj, features = VariableFeatures(obj)) # scale variable features in gex library
  #print("GEX Obj Scaled")
  obj[["HTO"]] <- CreateAssayObject(counts = htoMat) # add in HTO library
  #print("HTO Obj Created")
  obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR") # normalize HTO library using centered log-ratio transformation
  #print("DONE")
  return (obj)
}