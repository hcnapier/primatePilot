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

