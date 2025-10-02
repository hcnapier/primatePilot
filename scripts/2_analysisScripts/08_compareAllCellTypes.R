# 08 Compare All Cell Types
# Compare cluster similarity across species 
# Hailey Napier
# September 25, 2025

# 0.0 Setup ----
## 0.1 Load packages ----
require(Seurat)
require(dplyr)
require(reshape2)
require(ggplot2)
## 0.2 Load data ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/seuratObjs")
noGarbage <- readRDS("speciesObjList_cleanCellTypeLabled.rds")
setwd("/Users/haileynapier/Work/VertGenLab/Projects/zebrinEvolution/Data/geneLists/orthologs")
orthologsAll <- read.delim("human_MouseMarmosetMacaque_orthologs.txt", sep = "\t", header = T)
all121Orthologs <- orthologsAll %>%
  select(Macaque.gene.name, Macaque.homology.type, Mouse.gene.name, Mouse.homology.type, Gene.name) %>%
  filter(Macaque.homology.type == "ortholog_one2one") %>%
  filter(Mouse.homology.type == "ortholog_one2one")
names(all121Orthologs) <- tolower(names(all121Orthologs))
speciesNames <- c("human1", "human2", "rhesus", "mouse")


# 1.0 Pseudobulk ----
noGarbage_pseudobulk <- list()
for(currSpecies in speciesNames){
  if(grepl("human", currSpecies)){
    geneListName <- "gene.name"
  }else if(currSpecies == "rhesus"){
    geneListName <- "macaque.gene.name"
  }else{
    geneListName <- paste(currSpecies, "gene.name", sep = ".")
  }
  noGarbage_pseudobulk[currSpecies] <- AggregateExpression(noGarbage[[currSpecies]], features = pull(all121Orthologs[geneListName]), return.seurat = F)
}

# 2.0 Select common genes ----
## 2.1 Convert gene names to human ortholog ----
### Rhesus 
rhesus_convertToHuman <- data.frame(row.names(noGarbage_pseudobulk$rhesus))
names(rhesus_convertToHuman) <- "macaque.gene.name"
rhesus_convertToHuman <- inner_join(rhesus_convertToHuman, all121Orthologs) %>% distinct()
row.names(noGarbage_pseudobulk$rhesus) <- rhesus_convertToHuman$gene.name
### Mouse
mouse_convertToHuman <- data.frame(row.names(noGarbage_pseudobulk$mouse))
names(mouse_convertToHuman) <- "mouse.gene.name"
mouse_convertToHuman <- inner_join(mouse_convertToHuman, all121Orthologs) %>% distinct()
row.names(noGarbage_pseudobulk$mouse) <- mouse_convertToHuman$gene.name

## 2.2 Select genes that are common between all species ----
sharedGenes <- Reduce(intersect, list(rhesus_convertToHuman$gene.name, mouse_convertToHuman$gene.name, row.names(noGarbage_pseudobulk$human1), row.names(noGarbage_pseudobulk$human2)))

## 2.3 Subset each matrix to include only shared genes ----
for(currSpecies in speciesNames){
  noGarbage_pseudobulk[currSpecies] <- noGarbage_pseudobulk[[currSpecies]][sharedGenes,]
}


# 3.0 Format data ----
## 3.1 Convert sparse matrix to dense matrix 
noGarbage_pseudobulk <- lapply(noGarbage_pseudobulk, as.matrix)
for(currSpecies in speciesNames){
 # colsOrdered <- sort(colnames(noGarbage_pseudobulk[[currSpecies]]))
 # noGarbage_pseudobulk[[currSpecies]] <- noGarbage_pseudobulk[[currSpecies]][, colsOrdered]
  colnames(noGarbage_pseudobulk[[currSpecies]]) <- paste(colnames(noGarbage_pseudobulk[[currSpecies]]), currSpecies, sep = "_")
}

## 3.2 Merge all matrices 
pseudobulkMerged <- noGarbage_pseudobulk[["human1"]] %>% as.data.frame()
pseudobulkMerged$RowNames <- pseudobulkMerged %>% row.names()

for(currSpecies in speciesNames[2:4]){
  temp <- noGarbage_pseudobulk[[currSpecies]] %>% as.data.frame()
  temp$RowNames <- temp %>% row.names()
  pseudobulkMerged <- merge(pseudobulkMerged, temp)
  # order columns by cell type
  colsOrdered <- sort(colnames(pseudobulkMerged))
  pseudobulkMerged <- pseudobulkMerged[,colsOrdered]
}



# human 1 vs human 2
cormat <- cor(noGarbage_pseudobulk[["human1"]], noGarbage_pseudobulk[["human2"]])
melted_cormat_h1h2 <- melt(cormat)
names(melted_cormat_h1h2) <- c("Human1", "Human2", "value")
h1vsh2_corPlot <- ggplot(data = melted_cormat_h1h2, aes(Human1, Human2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

# human 1 vs mouse
