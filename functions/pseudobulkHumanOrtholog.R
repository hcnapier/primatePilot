# pseudobulkHumanOrtholog.R
# Hailey Napier
# October 3, 2025

# Function to pseudobulk by a given metadata value and convert to the human one2one ortholog

require(Seurat)
require(dplyr)

pseudobulkHumanOrtholog <- function(obj, groupByIdent = "ident", orthologDF, species, returnSeuratBool = FALSE){
  # aggregate expression
  pseudobulked <- AggregateExpression(obj, features = pull(orthologDF[species]), group.by = groupByIdent, return.seurat = returnSeuratBool, assays = "RNA")
  # convert to dense matrix 
  pseudobulked <- pseudobulked$RNA %>% as.matrix()
  # remove zero counts 
  pseudobulked <- pseudobulked[rowSums(pseudobulked)>0, ] 
  # convert to human ortholog
  speciesGeneName <- paste(species, "gene.name", sep = ".")
  orthologDF %>%
    select(gene.name) %>%
    select(all_of(speciesGeneName)) -> orthologDF
  pseudobulked <- pseudobulked %>% as.data.frame()
  pseudobulked[speciesGeneName] <- rownames(pseudobulked)
  pseudobulked <- inner_join(pseudobulked, orthologDF) 
  pseudobulked %>%
    select(-all_of(speciesGeneName))
}