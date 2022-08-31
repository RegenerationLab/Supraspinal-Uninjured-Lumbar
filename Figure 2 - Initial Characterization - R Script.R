# This script creates Figure 2
# The input for this script is LUN1, created in "1_Lumbar_Supraspinal Unified Workflow"

# Load libraries
library('Seurat')
library('ggplot2')
library('plyr')
library('dplyr')

#Input data

LUN1 <- readRDS("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/LUN1.rds")


#2A - a large labeled overview of the clusters

DimPlot(LUN1, group.by = 'seurat_clusters', label = FALSE, 
        label.size = 40, repel = TRUE, 
        pt.size = 4) + 
  theme(text = element_text(size = 100),
        axis.text = element_text(size = 70)) +
  NoLegend()  + NoAxes() #3000x3000


#2B - showing consistency across samples

Idents(LUN1) <- "seurat_clusters"

DimPlot(LUN1, group.by = 'seurat_clusters', label = FALSE, 
        label.size = 15, repel = TRUE, split.by = "orig.ident", 
        pt.size = 2) + 
  theme(text = element_text(size = 40),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #3000x1000

#2C - Create a dotplot with markers for off-target cells

Fig2markers <- read.csv("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/2_Input_Lists/Fig2markers.csv",
                        sep=",")

DefaultAssay(LUN1) <- "RNA"

DotPlot(LUN1, features = Fig2markers$Gene, dot.scale = 6, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        axis.text = element_text(size = 20)) #1000 x 1500 then x2 in Adobe

#2D - Create a series of violin plots for neurotransmitter types

VlnPlot(LUN1, "Slc17a7", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) #1000 x 600

VlnPlot(LUN1, "Slc17a6", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) #1000 x 600

VlnPlot(LUN1, "Gad2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) #1000 x 600

VlnPlot(LUN1, "Tph2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) #1000 x 600

VlnPlot(LUN1, "Slc6a2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) #1000 x 600


# 2E - Make feature plots for various markers

FeaturePlot(LUN1, c("Slc17a7"), label = FALSE, pt.size = 1,
            label.size = 5, repel = TRUE, order = FALSE) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #700 x 1000

FeaturePlot(LUN1, c("Slc17a6"), label = FALSE, pt.size = 1,
            label.size = 5, repel = TRUE, order = FALSE) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #700 x 1000

FeaturePlot(LUN1, c("Gad2"), label = FALSE, pt.size = 1,
            label.size = 5, repel = TRUE, order = FALSE) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #700 x 1000

FeaturePlot(LUN1, c("Tph2"), label = FALSE, pt.size = 1,
            label.size = 5, repel = TRUE, order = FALSE) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #700 x 1000

FeaturePlot(LUN1, c("Slc6a2"), label = FALSE, pt.size = 1,
            label.size = 5, repel = TRUE, order = FALSE) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #700 x 1000




