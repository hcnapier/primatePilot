# 0.0 Setup ----
## 0.1 Load packages ----
library(ggVennDiagram)
library(dplyr)
library(ggplot2)

## 0.2 Load species assignments lists ----
setwd("/Users/haileynapier/Work/VertGenLab/Projects/zebrinEvolution/Data/20250423_PrimatePilot/demuxSpeciesOuts/cellAssignmentLists")
nameList <- c("cellID", "species_assignment", "numCells", "likelihoodRatio")
orthoAll5000 <- read.delim("121Ortho_allKmers_5000cellswmostreads_species.assignments", sep = "\t", header = F)
hiExpAll5000 <- read.delim("hiExpGenes_allKmers_5000cellswmostreads_species.assignments", sep = "\t", header = F)
hiExpDefault5000 <- read.delim("hiExpGenes_defaultKmers_5000cellswmostreas_species.assignments", sep = "\t", header = F)
orthoDefault5000 <- read.delim("121Ortho_defaultKmers_5000cellswmostreads_species.assignments", sep = "\t", header = F)
names(orthoAll5000) <- nameList
names(hiExpAll5000) <- nameList
names(hiExpDefault5000) <- nameList
names(orthoDefault5000) <- nameList

# 1.0 Compare high expressed genes with all kmers or default kmers ----
## 1.1 Separate by species ----
### Human ----
hiExpAll5000_human <- hiExpAll5000 %>% filter(species_assignment == "human")
hiExpDefault5000_human <- hiExpDefault5000 %>% filter(species_assignment == "human")
hiExpHumanComp <- list(hiExpAll5000_human$cellID, hiExpDefault5000_human$cellID)
ggVennDiagram(hiExpHumanComp, 
              category.names = c("All Kmers", 
                                 "Default Kmers")) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none") +
  labs(title = "Human Cell Assignments",
       subtitle = "Highly expressed genes and 5000 cells with the highest endogenous reads") +
  coord_flip()
### Rhesus ----
hiExpComp_rhesus <- list(default_rhesus = hiExpDefault5000 %>% filter(species_assignment == "rhesus") %>% pull(cellID), 
                         all_rhesus = hiExpAll5000 %>% filter(species_assignment == "rhesus") %>% pull(cellID))
ggVennDiagram(hiExpComp_rhesus, 
              category.names = c("All Kmers", 
                                 "Default Kmers")) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none") +
  labs(title = "Rhesus Cell Assignments",
       subtitle = "Highly expressed genes and 5000 cells with the highest endogenous reads") +
  coord_flip()

## 1.2 Compare all species in one plot -----
hiExpComp <- list(default_human = hiExpDefault5000 %>% filter(species_assignment == "human") %>% pull(cellID),
                  all_human = hiExpAll5000 %>% filter(species_assignment == "human") %>% pull(cellID),
                  default_rhesus = hiExpDefault5000 %>% filter(species_assignment == "rhesus") %>% pull(cellID), 
                  all_rhesus = hiExpAll5000 %>% filter(species_assignment == "rhesus") %>% pull(cellID), 
                  default_mouse = hiExpDefault5000 %>% filter(species_assignment == "mouse") %>% pull(cellID), 
                  all_mouse = hiExpAll5000 %>% filter(species_assignment == "mouse") %>% pull(cellID)
)
ggVennDiagram(hiExpComp, label = "none",
              category.names = c("Default Kmers, Human", 
                                 "All Kmers, Human", 
                                 "Default Kmers, Rhesus",
                                 "All Kmers, Rhesus", 
                                 "Default Kmers Mouse", 
                                 "All Kmers Mouse")) + 
  scale_x_continuous(expand = expansion(mult = .2)) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  labs(title = "All Species Assignments, All Kmers vs. Default (20M) Kmers",
       subtitle = "Highly expressed genes, 5000 cells with the highest endogenous reads") 
# force upset plot
ggVennDiagram(hiExpComp, force_upset = T,
              category.names = c("Default Kmers, Human", 
                                 "All Kmers, Human", 
                                 "Default Kmers, Rhesus",
                                 "All Kmers, Rhesus", 
                                 "Default Kmers Mouse", 
                                 "All Kmers Mouse")) 

# 2.0 Compare high expressed genes and orthologous genes ----
## All species, all groups ----
hiExpOrthoComp <- list(hiExp_default_human = hiExpDefault5000 %>% filter(species_assignment == "human") %>% pull(cellID),
                  hiExp_all_human = hiExpAll5000 %>% filter(species_assignment == "human") %>% pull(cellID),
                  hiExp_default_rhesus = hiExpDefault5000 %>% filter(species_assignment == "rhesus") %>% pull(cellID), 
                  hiExp_all_rhesus = hiExpAll5000 %>% filter(species_assignment == "rhesus") %>% pull(cellID), 
                  hiExp_default_mouse = hiExpDefault5000 %>% filter(species_assignment == "mouse") %>% pull(cellID), 
                  hiExp_all_mouse = hiExpAll5000 %>% filter(species_assignment == "mouse") %>% pull(cellID), 
                  ortho_all_human = orthoAll5000 %>% filter(species_assignment == "human") %>% pull(cellID), 
                  ortho_all_rhesus = orthoAll5000 %>% filter(species_assignment == "rhesus") %>% pull(cellID), 
                  ortho_all_mouse = orthoAll5000 %>% filter(species_assignment == "mouse") %>% pull(cellID)
                  )
ggVennDiagram(hiExpOrthoComp, force_upset = T, label = "percent",
              category.names = c("Default Kmers, Human (highly expressed genes)", 
                                 "All Kmers, Human (highly expressed genes)", 
                                 "Default Kmers, Rhesus (highly expressed genes)",
                                 "All Kmers, Rhesus (highly expressed genes)", 
                                 "Default Kmers, Mouse (highly expressed genes)", 
                                 "All Kmers, Mouse (highly expressed genes)",
                                 "All Kmers, Human (one-to-one orthologs)",
                                 "All Kmers, Rhesus (one-to-one orthologs)",
                                 "All Kmers, Mouse (one-to-one orthologs)"
                                 )) 
## Just compare the groups with all kmers ----
hiExpOrthoComp_allkmers <- list(hiExp_all_human = hiExpAll5000 %>% filter(species_assignment == "human") %>% pull(cellID),
                                ortho_all_human = orthoAll5000 %>% filter(species_assignment == "human") %>% pull(cellID),
                                hiExp_all_rhesus = hiExpAll5000 %>% filter(species_assignment == "rhesus") %>% pull(cellID), 
                                ortho_all_rhesus = orthoAll5000 %>% filter(species_assignment == "rhesus") %>% pull(cellID), 
                                hiExp_all_mouse = hiExpAll5000 %>% filter(species_assignment == "mouse") %>% pull(cellID),
                                ortho_all_mouse = orthoAll5000 %>% filter(species_assignment == "mouse") %>% pull(cellID))
ggVennDiagram(hiExpOrthoComp_allkmers, label = "none",
              category.names = c("Human, Highly Expressed Genes", 
                                 "Human, One-to-one Orthologs", 
                                 "Rhesus, Highly Expressed Genes",
                                 "Rhesus, One-to-one Orthologs",
                                 "Mouse, Highly Expressed Genes",
                                 "Mouse, One-to-one Orthologs")) + 
  scale_x_continuous(expand = expansion(mult = .2)) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  labs(title = "All Species Assignments, Highly Expressed Genes vs. One-to-one Orthologs",
       subtitle = "All Kmers, 5000 cells with the highest endogenous reads") 

ggVennDiagram(hiExpOrthoComp_allkmers, force_upset = T,
              category.names = c("Human, Highly Expressed Genes", 
                                 "Human, One-to-one Orthologs", 
                                 "Rhesus, Highly Expressed Genes",
                                 "Rhesus, One-to-one Orthologs",
                                 "Mouse, Highly Expressed Genes",
                                 "Mouse, One-to-one Orthologs"))

# 3.0 Compare One-to-one orthologous genes ----
## 3.1 Compare all species ----
orthoComp <- list(ortho_default_human = orthoDefault5000 %>% filter(species_assignment == "human") %>% pull(cellID),
                                ortho_all_human = orthoAll5000 %>% filter(species_assignment == "human") %>% pull(cellID),
                                ortho_default_rhesus = orthoDefault5000 %>% filter(species_assignment == "rhesus") %>% pull(cellID), 
                                ortho_all_rhesus = orthoAll5000 %>% filter(species_assignment == "rhesus") %>% pull(cellID), 
                                ortho_default_mouse = orthoDefault5000 %>% filter(species_assignment == "mouse") %>% pull(cellID),
                                ortho_all_mouse = orthoAll5000 %>% filter(species_assignment == "mouse") %>% pull(cellID))

ggVennDiagram(orthoComp, label = "none",
              category.names = c("Human, 20M Kmers", 
                                 "Human, All Kmers", 
                                 "Rhesus, 20M Kmers",
                                 "Rhesus, All Kmers",
                                 "Mouse, 20M Kmers",
                                 "Mouse, All Kmers")) + 
  scale_x_continuous(expand = expansion(mult = .2)) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  labs(title = "All Species Assignments, 20M Kmers vs All Kmers",
       subtitle = "One-to-one orthologous genes, 5000 cells with the highest endogenous reads") 

ggVennDiagram(orthoComp, force_upset = T,
              category.names = c("Human, 20M Kmers", 
                                 "Human, All Kmers", 
                                 "Rhesus, 20M Kmers",
                                 "Rhesus, All Kmers",
                                 "Mouse, 20M Kmers",
                                 "Mouse, All Kmers")) 

## 3.2 Compare all species and doublets ----
orthoComp_doublets <- list(ortho_default_human = orthoDefault5000 %>% filter(species_assignment == "human") %>% pull(cellID),
                  ortho_all_human = orthoAll5000 %>% filter(species_assignment == "human") %>% pull(cellID),
                  ortho_default_rhesus = orthoDefault5000 %>% filter(species_assignment == "rhesus") %>% pull(cellID), 
                  ortho_all_rhesus = orthoAll5000 %>% filter(species_assignment == "rhesus") %>% pull(cellID), 
                  ortho_default_mouse = orthoDefault5000 %>% filter(species_assignment == "mouse") %>% pull(cellID),
                  ortho_all_mouse = orthoAll5000 %>% filter(species_assignment == "mouse") %>% pull(cellID), 
                  ortho_default_humanRhesus = orthoDefault5000 %>% filter(species_assignment == "human+rhesus") %>% pull(cellID),
                  ortho_all_humanRhesus = orthoAll5000 %>% filter(species_assignment == "human+rhesus") %>% pull(cellID),
                  ortho_default_humanMouse = orthoDefault5000 %>% filter(species_assignment == "human+mouse") %>% pull(cellID),
                  ortho_all_humanMouse = orthoAll5000 %>% filter(species_assignment == "human+mouse") %>% pull(cellID),
                  ortho_default_mouseRhesus = orthoDefault5000 %>% filter(species_assignment == "mouse+rhesus") %>% pull(cellID),
                  ortho_all_mouseRhesus = orthoAll5000 %>% filter(species_assignment == "mouse+rhesus") %>% pull(cellID))
ggVennDiagram(orthoComp_doublets, force_upset = T,
              category.names = c("Human, 20M Kmers", 
                                 "Human, All Kmers", 
                                 "Rhesus, 20M Kmers",
                                 "Rhesus, All Kmers",
                                 "Mouse, 20M Kmers",
                                 "Mouse, All Kmers", 
                                 "Human + Rhesus, 20M Kmers", 
                                 "Human + Rhesus, All Kmers", 
                                 "Human + Mouse, 20M Kmers", 
                                 "Human + Mouse, All Kmers", 
                                 "Mouse + Rhesus, 20M Kmers", 
                                 "Mouse + Rhesus, All Kmers")) 

# 4.0 Compare all cells RefSeq to all cells high exp genes ----



