# General 
install.packages("dplyr")
install.packages("ggplot2")

# Species Assignment Comparisons
install.packages("ggVennDiagram")

# Analysis
install.packages("Seurat")
install.packages("patchwork")

# HTO demultiplexing
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ShortRead")
devtools::install_github('chris-mcginnis-ucsf/MULTI-seq')

# normalization
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("glmGamPoi")

# integration
install.packages("rliger")

# GO terms
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler", force = T)
BiocManager::install("org.Hs.eg.db", force = T)
