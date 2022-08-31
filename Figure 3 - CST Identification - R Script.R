#The purpose of this script is to establish CST identity in the dataset (Figure 3)
# The input for this script is LUN1, created in "1_Lumbar_Supraspinal Unified Workflow"

# Load libraries
library('Seurat')
library('ggplot2')
library('plyr')
library('dplyr')

#Load relevant data
LUN1 <- readRDS("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/LUN1.rds")


#Positive markers for layer V cortical / CST: Crym, Fezf2, Bcl11b, 
# Bcl6, Pdlim1, Cacna1h, Mylip, Kcng1
#markers I considered but didn't use:(Slco2a1, Crim1, Gng7, Cthrc1)

#Positive markers for Cortex: Satb2, Iptka

#Negative markers for CST: Sulf1, Rims3, Slc30a3, Hpgd (thalamus projecting, didn't work)



# (Not for figure but helpful for orientation - a dimplot to get oriented. CST is the large group of 0 and 2)

Idents(LUN1) <- "seurat_clusters"
DimPlot(LUN1, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# 3A Bcl11b FeaturePlot

FeaturePlot(LUN1, c("Bcl11b"), label = TRUE, pt.size = 1,
            label.size = 10, repel = TRUE, order = FALSE) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #1000 x 1000

# 3B Crym FeaturePlot
FeaturePlot(LUN1, c("Crym"), label = TRUE, pt.size = 1,
            label.size = 10, repel = TRUE, order = FALSE) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #1000 x 1000

# 3C Fezf2 FeaturePlot
FeaturePlot(LUN1, c("Fezf2"), label = TRUE, pt.size = 1,
            label.size = 10, repel = TRUE, order = TRUE) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #1000 x 1000


# 3D-K Example Violin Plots

VlnPlot(LUN1, "Bcl11b", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(text = element_text(size = 20), 
        axis.text = element_text(size = 20)) #1000 x 300

VlnPlot(LUN1, "Crym", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 20)) #1000 x 300

VlnPlot(LUN1, "Fezf2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 20)) #1000 x 300

VlnPlot(LUN1, "Bcl6", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 20)) #1000 x 300

VlnPlot(LUN1, "Pdlim1", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(text = element_text(size = 20)) #1000 x 300

VlnPlot(LUN1, "Cacna1h", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(text = element_text(size = 20)) #1000 x 300

VlnPlot(LUN1, "Mylip", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(text = element_text(size = 20)) #1000 x 300

VlnPlot(LUN1, "Kcng1", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(text = element_text(size = 20)) #1000 x 300






