#The purpose of this script is to generate violin plots for Figure 4 (initial cluster ID)

## This script starts from LUN_Final and an initial ordering determined by 
# LUN_Final_Map1 (for seurat_clusters) and LUN_Final_Order1 for (manual.clusters.1)

# Load libraries
library('Seurat')
library('ggplot2')
library('plyr')
library('dplyr')

#Load relevant data
LUN_Final <- readRDS("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/LUN_Final.rds")

#violin and feature plots for 4A,E,I,M,Q,U

Idents(LUN_Final) <- "seurat_clusters"

# Sim1
VlnPlot(LUN_Final, "Sim1", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Sim1"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = FALSE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))

# Plagl1
VlnPlot(LUN_Final, "Plagl1", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Plagl1"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = FALSE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))


# ttc6
VlnPlot(LUN_Final, "Ttc6", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Ttc6"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = FALSE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))


# Cartpt
VlnPlot(LUN_Final, "Cartpt", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Cartpt"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = TRUE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))

# Slc6a2
VlnPlot(LUN_Final, "Slc6a2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Slc6a2"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = FALSE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))


# Prdm6
VlnPlot(LUN_Final, "Prdm6", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Prdm6"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = FALSE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))











