# This script creates Figure 6, the creation of a Dimplot and Dotplot
# with marker genes for the final clusters.

# Input for this script is LUN1_Final

# Load libraries
library('Seurat')
library('ggplot2')
library('plyr')
library('dplyr')

#Load relevant data
LUN_Final <- readRDS("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/LUN_Final.rds")


#6A - a large overview of the clusters


DimPlot(LUN_Final, group.by = 'manual.clusters.1', label = FALSE, 
        label.size = 40, repel = TRUE, 
        pt.size = 4) + 
  theme(text = element_text(size = 100),
        axis.text = element_text(size = 70)) +
  NoLegend()  + NoAxes() #2000x3000



#6B - Dotplot with marker genes

# Read csv with selected marker genes

LUN_Final_DotPlot2 <- read.csv("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/2_Input_Lists/LUN_Final_DotPlot2.csv",
                               sep = ",")


DefaultAssay(LUN_Final) <- "RNA"
Idents(LUN_Final) <- "manual.clusters.1"

DotPlot(LUN_Final, features = LUN_Final_DotPlot2$Marker, dot.scale = 6, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 20)) #1000 x 1500 then x2 in Adobe