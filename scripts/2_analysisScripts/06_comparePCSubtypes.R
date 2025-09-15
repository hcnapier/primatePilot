# 06 Compare PC Subtypes
# Compare genes in PC subtypes across species 
# Hailey Napier
# August 27, 2025

# 0.0 Setup ----
## 0.1 Load packages ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)

## 0.2 Load data ----
#setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data/seuratObjs")
#humanPCs_merged <- readRDS("humanPCsMerged.rds")
#noGarbage <- readRDS("speciesObjList_cleanCellTypeLabled.rds")
setwd("/Users/haileynapier/Work/VertGenLab/Projects/zebrinEvolution/Data/geneLists")
orthologsAll <- read.delim("human_rhesusMouseOrthologs.txt", sep = "\t", header = T)
mouseOrthologs <- orthologsAll %>%
  select(Gene.name, Mouse.gene.name, Mouse.homology.type) %>%
  filter(Mouse.homology.type == "ortholog_one2one")
names(mouseOrthologs) <- c("humanGeneName", "geneName", "homologyType")

rhesusOrthologs <- orthologsAll %>%
  select(Gene.name, Macaque.gene.name, Macaque.homology.type) %>%
  filter(Macaque.homology.type == "ortholog_one2one")
names(rhesusOrthologs) <- c("humanGeneName", "geneName", "homologyType")

rm(orthologsAll)

setwd("/Users/haileynapier/Work/VertGenLab/Projects/zebrinEvolution/Code/retinaOUProcess/RetinaOUProcess/data/RAnalysisOutput")
G_list <- readRDS("geneEntrezIDList.rds")

setwd("/Users/haileynapier/Work/VertGenLab/Projects/zebrinEvolution/Data/geneLists")
diseaseGenes <- read.csv("diseaseGenes.txt", header = F)
diseaseGenes <- diseaseGenes$V1

# 1.0 Subset & Re-cluster PCs ----
## 1.1 Mouse ----
mousePCs <- noGarbage[['mouse']] %>% subset(idents = "Purkinje")
mousePCs <- mousePCs %>%
  RunPCA(verbose = F) %>%
  RunUMAP(dims = 1:30, verbose = F) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters()
DimPlot(mousePCs)
  
## 1.2 Rhesus ----
rhesusPCs <- noGarbage[['rhesus']] %>% subset(idents = "Purkinje")
rhesusPCs <- rhesusPCs %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters()
DimPlot(rhesusPCs, pt.size = 2)

## 1.3 Human ---- 
# already re-clustered
DefaultAssay(humanPCs_merged) <- "SCT"
humanPCs_merged <- JoinLayers(humanPCs_merged)
DimPlot(humanPCs_merged, pt.size = 5, cols = c("#79A3BC", "#E2B1AC"))


# 2.0 Find subtype markers ----
mousePCMarkers <- mousePCs %>%
  FindMarkers(ident.1 = "3", ident.2 = "2") 
mousePCMarkersFilt <- mousePCMarkers %>%
  filter(p_val_adj <= 0.1)
mousePCMarkers$geneName <- rownames(mousePCMarkers)
mousePCMarkersFilt$geneName <- rownames(mousePCMarkersFilt)

rhesusPCMarkers <- rhesusPCs %>%
  FindMarkers(ident.1 = "0", ident.2 = "1") 
rhesusPCMarkersFilt <- rhesusPCMarkers %>%
  filter(p_val_adj <= 0.1)
rhesusPCMarkers$geneName <- rownames(rhesusPCMarkers)
rhesusPCMarkersFilt$geneName <- rownames(rhesusPCMarkersFilt)

humanPCMarkers <- humanPCs_merged %>%
  FindMarkers(ident.1 = "0", ident.2 = "1") 
humanPCMarkersFilt <- humanPCMarkers %>%
  filter(p_val_adj <= 0.1)
humanPCMarkers$humanGeneName <- rownames(humanPCMarkers)
humanPCMarkersFilt$humanGeneName <- rownames(humanPCMarkersFilt)
humanPCMarkerList <- humanPCMarkersFilt$gene

# 3.0 Convert subtype markers to human ortholog ----
mousePCMarkers <- left_join(mousePCMarkers, mouseOrthologs, by = "geneName")
mousePCMarkersFilt <- left_join(mousePCMarkersFilt, mouseOrthologs, by = "geneName")
rhesusPCMarkers <- left_join(rhesusPCMarkers, rhesusOrthologs, by = "geneName")
rhesusPCMarkersFilt <- left_join(rhesusPCMarkersFilt, rhesusOrthologs, by = "geneName")


# 4.0 Compare marker genes ----
primateSharedMarkers <- intersect(humanPCMarkersFilt$humanGeneName, rhesusPCMarkersFilt$humanGeneName)
allSharedMarkers <- intersect(primateSharedMarkers, mousePCMarkersFilt$humanGeneName)
allSharedMarkers
rhesusMouseSharedMarkers <- intersect(rhesusPCMarkersFilt$humanGeneName, mousePCMarkersFilt$humanGeneName)
humanMouseSharedMarkers <- intersect(humanPCMarkersFilt$humanGeneName, mousePCMarkersFilt$humanGeneName)

# 5.0 Plot ----
## 5.1 Venn diagrams ----
ggVennDiagram(list(mousePCMarkersFilt$humanGeneName, humanPCMarkersFilt$humanGeneName), 
              category.names = c("Mouse", "Human")) +
  coord_flip()
humanMouseSharedMarkers

ggVennDiagram(list(rhesusPCMarkersFilt$humanGeneName, humanPCMarkersFilt$humanGeneName), 
              category.names = c("Rhesus", "Human")) +
  coord_flip()
primateSharedMarkers

## 5.2 Feature plots ----
FeaturePlot(humanPCs_merged, features = "GRID2", pt.size = 3)

## 5.3 Violin plots ----
VlnPlot(humanPCs_merged, features = c("RNR2", "COX3", "GRID2", "CA8"), ncol = 2, cols = c("#79A3BC", "#E2B1AC"))


# 6.0 GO Terms ----
## 6.1 Get entrez gene ID ----
entrezIDList <- G_list %>% dplyr::select(external_gene_name, entrezgene_id)
rm(G_list)
# If multiple entrez IDs, arbitrarily select first
entrezIDList <- entrezIDList %>% 
  group_by(external_gene_name) %>%
  arrange(entrezgene_id) %>% 
  top_n(1)
names(entrezIDList)[1] <- "gene"
names(mousePCMarkers)[7] <- "gene"
names(mousePCMarkersFilt)[7] <- "gene"
names(humanPCMarkers)[6] <- "gene"
names(humanPCMarkersFilt)[6] <- "gene"
mousePCMarkers_entrez <- inner_join(mousePCMarkers, entrezIDList)
mousePCMarkersFilt_entrez <- inner_join(mousePCMarkersFilt, entrezIDList)
humanPCMarkers_entrez <- inner_join(humanPCMarkers, entrezIDList)
humanPCMarkersFilt_entrez <- inner_join(humanPCMarkersFilt, entrezIDList)

## 6.2 Get GO terms ----
humanPCGO <- enrichGO(as.character(humanPCMarkersFilt_entrez$entrezgene_id), 
                              OrgDb = org.Hs.eg.db,
                              ont = 'BP',
                              keyType = "ENTREZID", 
                              universe = as.character(entrezIDList$entrezgene_id),
                              qvalueCutoff = 0.5, 
                              pvalueCutoff = 0.1)
clusterProfiler::dotplot(humanPCGO)

mousePCGO <- enrichGO(as.character(mousePCMarkersFilt_entrez$entrezgene_id), 
                      OrgDb = org.Hs.eg.db,
                      ont = 'BP',
                      keyType = "ENTREZID", 
                      universe = as.character(mousePCMarkers_entrez$entrezgene_id),
                      qvalueCutoff = 0.5, 
                      pvalueCutoff = 0.5)
clusterProfiler::dotplot(mousePCGO)

mousePCGO <- enrichGO(as.character(mousePCMarkersFilt_entrez$entrezgene_id), 
                      OrgDb = org.Hs.eg.db,
                      ont = 'MF',
                      keyType = "ENTREZID", 
                      universe = as.character(mousePCMarkers_entrez$entrezgene_id),
                      qvalueCutoff = 0.5, 
                      pvalueCutoff = 0.5)
clusterProfiler::dotplot(mousePCGO)


# 7.0 Compare subtype marker genes to disease gene list


# 7.0 Save gene lists ----
setwd("~/Work/VertGenLab/Projects/zebrinEvolution/Code/primatePilot/data")
write.csv(mousePCMarkers, "mousePCSubtypeMarkers.csv")
write.csv(humanPCMarkers, "humanPCSubtypeMarkers.csv")

