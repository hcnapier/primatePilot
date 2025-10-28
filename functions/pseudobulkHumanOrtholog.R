# pseudobulkHumanOrtholog.R
# Hailey Napier
# October 3, 2025

# Function to pseudobulk by a given metadata value and convert to the human one2one ortholog

require(Seurat)
require(dplyr)

pseudobulkHumanOrtholog <- function(obj, groupByIdent = "ident", orthologDF, species, datasetID){
  # aggregate expression
  speciesGeneName <- paste(species, "gene.name", sep = ".")
  pseudobulked <- AggregateExpression(obj, features = pull(orthologDF[speciesGeneName]), group.by = groupByIdent, return.seurat = FALSE, assays = "RNA")
  # convert to dense matrix 
  pseudobulked <- pseudobulked$RNA %>% as.matrix()
  # remove zero counts 
  pseudobulked <- pseudobulked[rowSums(pseudobulked)>0, ] %>% as.data.frame()
  names(pseudobulked) <- datasetID
  # convert to human ortholog
  pseudobulked[speciesGeneName] <- rownames(pseudobulked)
  pseudobulked <- inner_join(pseudobulked, orthologDF) 
  pseudobulked %>%
    dplyr::select(all_of(datasetID), gene.name) -> pseudobulked
  return(pseudobulked)
}
