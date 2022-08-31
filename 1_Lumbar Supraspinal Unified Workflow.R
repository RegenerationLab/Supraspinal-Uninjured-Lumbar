#The purpose of this script is to input and merge three datasets of uninjured supraspinal neurons that were
#retrogradely traced from lumbar spinal cord, remove spurious objects and contaminating cell types, then 
#assign cell identities to the resulting clusters and find marker genes for each.

#At various points the script references additional scripts that were used to create specific figures for publication.

#This script is intended to be executed manually step-by-step, with occasional user judgments informing the 
#precise variables.

#Note that the path names will need to be adjusted according to your own file structure. The barcode, gene, and 
#matrix files for all three samples can be downloaded from GEO. The input files, a series of CSVs that provide
#lists of marker genes and which help to reorder the data, are available on Github within "2_Input_Lists"


# Load libraries
library('Seurat')
library('ggplot2')
library('dplyr')

#########   Input Data    ###########################################
#Step 1: Input data and remove high- and low-count nuclei


Sample1.data <- Read10X(data.dir = "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/1_SourceFiles/Sample1")
Sample1 <- CreateSeuratObject(counts = Sample1.data, project = "Sample1", 
                              min.cells = 3, min.features = 1000)
Sample1 
# Initially 2710
VlnPlot(Sample1, features = c("nFeature_RNA"))
Sample1 <- subset(Sample1, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000)
Sample1 
# Now 2654

Sample2.data <- Read10X(data.dir = "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/1_SourceFiles/Sample2")
Sample2 <- CreateSeuratObject(counts = Sample2.data, project = "Sample2", 
                              min.cells = 3, min.features = 1000)
Sample2
# Initially 2770
VlnPlot(Sample2, features = c("nFeature_RNA"))
Sample2 <- subset(Sample2, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000)
Sample2 
# Now 2608


Sample3.data <- Read10X(data.dir = "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/1_SourceFiles/Sample3")
Sample3 <- CreateSeuratObject(counts = Sample3.data, project = "Sample3", 
                              min.cells = 3, min.features = 1000)
Sample3
#Initially 3579
VlnPlot(LUN1, features = c("nFeature_RNA"))
Sample3 <- subset(Sample3, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000)
Sample3
#Now 3518

# Save RDS files
saveRDS(Sample1, "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/Sample1_1000_6000.rds")
saveRDS(Sample2, "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/Sample2_1000_8000.rds")
saveRDS(Sample3, "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/Sample3_1000_7000.rds")


#To open the individual sample files if needed
Sample1 <- readRDS("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/Sample1_1000_6000.rds")
Sample2 <- readRDS("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/Sample2_1000_8000.rds")
Sample3 <- readRDS("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/Sample3_1000_7000.rds")

############# Merge Data #####################
#Step 2 - merge the three datasets

LUN_merge.list <- list (Sample1, Sample2, Sample3)

LUN_merge.list <- lapply(X=LUN_merge.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = LUN_merge.list)

LUN_merge.anchors <- FindIntegrationAnchors(object.list = LUN_merge.list,
                                            anchor.features = features)

LUN <- IntegrateData(anchorset = LUN_merge.anchors)

DefaultAssay(LUN) <- "integrated"

LUN <- ScaleData(LUN, verbose = FALSE)

LUN <- RunPCA(LUN, npcs = 50, verbose = FALSE)

LUN <- RunUMAP(LUN, reduction = "pca", dims = 1:30)

LUN <- FindNeighbors(LUN, reduction = "pca", dims = 1:30)

LUN <- FindClusters(LUN, resolution = 0.5)

DefaultAssay(LUN) <- "RNA"

LUN

#Indicates 8780 nuclei


# Some visualizations for initial assessment of clustering and sample bias


DimPlot(LUN, group.by = 'seurat_clusters', split.by = 'orig.ident',
        label = TRUE, label.size = 5) + NoLegend() # 1500 x 500

saveRDS(LUN, "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/LUN.rds")
LUN <- readRDS("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/LUN.rds")

###############################################################################
#                             Save and Restart Point
###############################################################################



#Step 3: Remove the clearly stressed, dying, or artificially low cell count clusters, then recluster

FeaturePlot(LUN, c("Creb5"), label = TRUE, 
            label.size = 5, repel = TRUE, order = TRUE) #500 x 500

FeaturePlot(LUN, c("Atf3"), label = TRUE, 
            label.size = 5, repel = TRUE, order = TRUE) #500 x 500

FeaturePlot(LUN, c("nFeature_RNA"), label = TRUE, 
            label.size = 5, repel = TRUE, order = TRUE) #500 x 500

VlnPlot(LUN, "nFeature_RNA", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(text = element_text(size = 20)) #1500 x 500

####Based on the above QC metrics, cluster 6, 7, 20, and 13 are removed

Idents(LUN) <- "seurat_clusters"

LUN1 <- subset(x = LUN, idents = c("6", "7", "20", "13"), invert = TRUE)

#A "sanity check" visualization to confirm proper removal

DimPlot(LUN1, group.by = 'seurat_clusters', split.by = 'orig.ident', label = TRUE, label.size = 5) + NoLegend()


############  Renormalize the new dataset that lacks the dying and low-count clusters cells   ################

DefaultAssay(LUN1) <- "RNA"

LUN1 <- NormalizeData(LUN1, normalization.method = "LogNormalize", scale.factor = 10000)

LUN1 <- FindVariableFeatures(LUN1, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(LUN1)

LUN1 <- ScaleData(LUN1, features = all.genes)

LUN1 <- RunPCA(LUN1, features = VariableFeatures(object = LUN1), npcs = 50)

LUN1 <- FindNeighbors(LUN1, dims = 1:30)
LUN1 <- FindClusters(LUN1, resolution = 0.5)

LUN1 <- RunUMAP(LUN1, dims = 1:30)

LUN1 

#total object is now 7748

# Various visualizations to decide whether clustering is appropriate

FeaturePlot(LUN1, c("Gad2"), 
            label = TRUE, label.size = 3,
            repel = FALSE, order = TRUE, ) #500 x 1000

DimPlot(LUN1, group.by = 'seurat_clusters',  label = TRUE, label.size = 5) + NoLegend()

#### Determine the number of nuclei in each cluster

Idents(LUN1) <- "seurat_clusters"
LUN1_cell.num <- table(Idents(LUN1))
LUN1_cell.num

LUN1

write.csv(x = LUN1_cell.num, 
          file = "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/3_Output_Lists/LUN1_cell_numbers.csv",
          quote = FALSE)



saveRDS(LUN1, "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/LUN1.rds")
LUN1 <- readRDS("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/LUN1.rds")


############################################################

##Launch point for Figure 2, Initial Characterization    ###

##Also Launch for Figure 3, CST identification          ####

############################################################


### Step 4: Identify CST and also establish markers to exclude non-
#CST cortical cells and determine their number

# Make violin plots of Satb2, Camk2a, and Crym and Fezf2 

Idents(LUN1) <- "seurat_clusters"

VlnPlot(LUN1, "Satb2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15)) #Sized for display

VlnPlot(LUN1, "Camk2a", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15)) #Sized for display

VlnPlot(LUN1, "Crym", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15)) #Sized for display

VlnPlot(LUN1, "Fezf2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15)) #Sized for display

VlnPlot(LUN1, "Slc30a3", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15)) #Sized for display


VlnPlot(LUN1, "nFeature_RNA", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(text = element_text(size = 20)) #1500 x 500

#Based on the plots above, cluster 13 appears cortical but not corticospinal


#As further evidence, check markers on a dot plot: they are pan-cortex, layer V, then intracortical

DotPlot(LUN1, features = c("Mef2c", "Grin2a", "Camk2a", "Ptk2b", "Satb2", 
                           "Crym", "Fezf2", "Bcl6", "Slc30a3"), dot.scale = 4, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5),
        axis.text = element_text(size = 10)) #Sized for display
# Note that cluster 13 has high levels of slc30a3, which is an intracortical marker
# Based on these plots, cluster 0 and 2 are likely CST and 13 is non-CST cortical

### Create lists of differentially expressed genes that mark putative
# CST clusters and any non-layer V cortical neurons


LUN1_13vsCST_Markers <- FindMarkers(LUN1,
                                    ident.1 = c("13"),
                                    ident.2 = c("0", "2"))

write.csv(x = LUN1_13vsCST_Markers, 
          file = "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/3_Output_Lists/LUN1_13vsCST.csv",
          quote = FALSE)

#Create a final set of nuclei from which the non-CST cortical neurons have been removed

LUN_Final <- subset(LUN1, idents = c("13"), invert = TRUE)

DimPlot(LUN_Final, group.by = 'seurat_clusters',  label = TRUE, label.size = 5) + NoLegend()


# Step 5 - Identify genes that mark the remaining clusters (note the command below is just an example
#for brevity) 

Idents(LUN_Final) <- "seurat_clusters"

LUN_Final_6 <- FindMarkers(LUN_Final,
                            ident.1 = c("6"))

write.csv(x = LUN_Final_6, 
          file = "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/3_Output_Lists/LUN_Final_6.csv",
          quote = FALSE)

#By cross-comparing the genes enriched in each cluster with data from the Allen Brain Atlas, key marker genes
#are identified and it is possible to make reasonable guesses about cell identity in each cluster.

#Step 6 - reorder the clusters to approximate rostral-caudal sequence of the putative cell identities

#Read a csv with standard markers

LUN_Final_DotPlot1 <- read.csv("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/2_Input_Lists/LUN_Final_DotPlot1.csv",
                               sep = ",")

#Create a dotplot that makes it easy to see which clusters express each marker gene

Idents(LUN_Final)  <- "seurat_clusters"

DefaultAssay(LUN_Final) <- "RNA"

DotPlot(LUN_Final, features = LUN_Final_DotPlot1$Marker, dot.scale = 9, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5),
        axis.text = element_text(size = 10)) #Sized for display


#Assign (or "map") each seurat cluster ID to an anatomical name, called a "manual.cluster.1" 

#Read a csv that has the seurat/manual.cluster correspondence as two columns with matched terms by row

LUN_Final_Map1 <- read.csv("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/2_Input_Lists/LUN_Final_Map1.csv",
                           sep = ",")

#Use the following commands to create a new slot with the manual.clusters (anatomical names)
Idents(LUN_Final) <- "seurat_clusters"

LUN_Final@meta.data$manual.clusters.1 <- plyr::mapvalues(
  x = LUN_Final@meta.data$seurat_clusters,
  from = LUN_Final_Map1$seurat_clusters,
  to = LUN_Final_Map1$manual.cluster.1)

#Check the new data

Idents(LUN_Final) <- "manual.clusters.1"
levels(LUN_Final)

#From the command above the anatomical terms should appear - but not in a rostral/caudal order
#By having the terms ordered properly in the input csv already loaded, the following command will reorder the seurat
#clusters in the correct Rostral/Caudal sequence

Idents(LUN_Final) <- "seurat_clusters"
LUN_Final@meta.data$seurat_clusters <- factor(
  x = LUN_Final@meta.data$seurat_clusters,
  levels = LUN_Final_Map1$seurat_clusters)

#Check that the seurat cluster reordering worked
Idents(LUN_Final) <- "seurat_clusters"
levels(LUN_Final)


#A further check - now the dotplot markers should be arranged in a diagonal line

DefaultAssay(LUN_Final) <- "RNA"

DotPlot(LUN_Final, features = LUN_Final_DotPlot1$Marker, dot.scale = 9, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5),
        axis.text = element_text(size = 10)) #Sized for display


#Order manual.clusters.1 in a rostral-caudal sequence

#Load a csv that has the anatomical names (manual.clusters.1) in rostral/caudal order

LUN_Final_Order1 <- read.csv("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/2_Input_Lists/LUN_Final_Order1.csv",
                             sep = ",")

#Use the following command to reorder the manual clusters

LUN_Final@meta.data$manual.clusters.1 <- factor(
  x = LUN_Final@meta.data$manual.clusters.1,
  levels = LUN_Final_Order1$Order)

#Check the reordering
Idents(LUN_Final) <- "manual.clusters.1"
levels(LUN_Final)


#Perform an additional check with more putative marker genes

#Load an expanded list of putative marker genes

LUN_Final_DotPlot2 <- read.csv("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/2_Input_Lists/LUN_Final_DotPlot2.csv",
                               sep = ",")

#Create a dotplot to check specific expression

DefaultAssay(LUN_Final) <- "RNA"



Idents(LUN_Final) <- "manual.clusters.1"
DotPlot(LUN_Final, features = LUN_Final_DotPlot2$Marker, dot.scale = 4, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10)) #Sized for display

saveRDS(LUN_Final, "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/LUN_Final.rds")
LUN_Final <- readRDS("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/LUN_Final.rds")

###############################################################################
#                             Save and Restart Point
###############################################################################


###############################################################################

#### Launch Point for Figure 4 - Cluster ID
#### Launch Point for Figure 5 - Hindbrain ID
#### Launch Point for Figure 6 - Final Dimplot and Dotplot
###############################################################################




#Export data to summarize differences between clusters

DefaultAssay(LUN_Final) <- "RNA"

avgexp <- AverageExpression(LUN_Final, assay = 'RNA', 
                            features = rownames(LUN_Final@assays$RNA@counts),
                            group.by = c('manual.clusters.1'))
avgexp$RNA

write.csv(file = "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/3_Output_Lists/AvgExpression.csv", x = avgexp)



#Identify marker genes for manual.clusters.1 and export to csv

LUN_Final_All_Markers <- FindAllMarkers(LUN_Final, only.pos = TRUE)

write.csv(x = LUN_Final_All_Markers, 
          file = "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/3_Output_Lists/LUN_Final_All_Markers.csv",
          quote = FALSE)




#### Determine the number of nuclei in each cluster

Idents(LUN_Final) <- "manual.clusters.1"
LUN_Final_cell.num <- table(Idents(LUN_Final))
LUN_Final_cell.num

write.csv(x = LUN_Final_cell.num, 
          file = "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/3_Output_Lists/LUN_Final_Cell_Numbers.csv",
          quote = FALSE)



### Find DEGs in CST versus hindbrain populations

Idents(LUN_Final) <- "manual.clusters.1"
levels(LUN_Final)
LUN_CST_Markers <- FindMarkers(LUN_Final,
                               ident.1 = c("CST"),
                               ident.2 = c("HB-1",
                                           "HB-2",
                                           "HB-3",
                                           "HB-4",
                                           "HB-5"))


write.csv(x = LUN_CST_Markers, 
          file = "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/3_Output_Lists/LUN_CST_Markers.csv",
          quote = FALSE)



######Final save of all seurat objects

saveRDS(Sample1, "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Githubn/Sample1_1000_6000.rds")
saveRDS(Sample2, "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/GithubN/Sample2_1000_8000.rds")
saveRDS(Sample3, "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/GithubN/Sample3_1000_7000.rds")

saveRDS(LUN, "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/GithubN/LUN.rds")
saveRDS(LUN1, "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/GithubN/LUN1.rds")
saveRDS (LUN_Final, "Z:/Blackmore_Lab/7_Single_Cell_Manuscript/GithubN/LUN_Final.rds")






