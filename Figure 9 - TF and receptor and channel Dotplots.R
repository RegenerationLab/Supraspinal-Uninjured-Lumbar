# This file is to load and visualize dot plot of interesting genes, resulting in creation of Figure 9

# Load libraries
library('Seurat')
library('ggplot2')
library('plyr')
library('dplyr')

#Load relevant data
LUN_Final <- readRDS("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/LUN_Final.rds")

#Export a file that lists the average expression of all transcripts for each cluster

DefaultAssay(LUN_Final) <- "RNA"

avgexp <- AverageExpression(LUN_Final, assay = 'RNA', 
                            features = rownames(LUN_Final@assays$RNA@counts),
                            group.by = c('manual.clusters.1'))
avgexp$RNA

write.csv(file = "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LuN/Fig9_Lists/AvgExpression.csv", x = avgexp)



#Input list of TFs

TFs <- read.csv("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/2_Input_Lists/Tfs_v1.csv",
                sep = ",")


DotPlot(LUN_Final, features = TFs$TFs, dot.scale = 5, dot.min = 0.01) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10)) #600 x 1200 for v1



#Input list of receptors

Receptors <- read.csv("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/2_Input_Lists/Receptors_v1.csv",
                      sep = ",")


Idents(LUN_Final)  <- "manual.clusters.1"

DefaultAssay(LUN_Final) <- "RNA"

DotPlot(LUN_Final, features = Receptors$Receptor, dot.scale = 5, dot.min = 0.01) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10)) #600 x 1200


#Input list of Voltage gated ion channels

VGCs <- read.csv("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/2_Input_Lists/VGCs_v3.csv",
                      sep = ",")


Idents(LUN_Final)  <- "manual.clusters.1"

DefaultAssay(LUN_Final) <- "RNA"

DotPlot(LUN_Final, features = VGCs$VGC, dot.scale = 4, dot.min = 0.01) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10)) #600 x 1200


