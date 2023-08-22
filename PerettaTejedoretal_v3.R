## Analysis pipeline for Peretta-Tejedor et al. to analyse CTN1 and related transcripts in murine nephrotoxic nephritis
## Data derived from 2 x Ctrl vs. 2  x NTN (day 1) vs. 2 x NTN (day 5)  scRNAseq datsets (C57Bl/6 background)
## See paper for methdological details and original dataset: https://jasn.asnjournals.org/content/early/2020/07/10/ASN.2020020220
## Daniyal Jafree & Gideon Pomeranz (UCL) | 22nd August 2023 | Version 3

#----------------------------------------------------------------------------------------------------------------#

# Load up packages and set working directory. Change as required
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(tidyverse)
  library(Matrix)
  library(patchwork)
  library(tidyselect)
  library(harmony)
  library(SoupX)
  library(ggrepel)
  library(DoubletFinder)
  setwd("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/glom_analysis/Nuria_CTF1")

#----------------------------------------------------------------------------------------------------------------#

## PART 1: Script for loading individual samples into Seurat objects
# Load data from file, deletes Ensembl IDs and makes gene names unique and create Seurat object for control number 2
  # Add project to differentiate between control and NTN samples
  ctrl2.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/glom_analysis/Nuria_CTF1/Raw_data/GSM4409508_control_2.txt.gz"),sep="\t", header = T, row.names = NULL)
  ctrl2.data <- subset(ctrl2.data, select = -c(IGIS))
  names <- make.unique(as.character(ctrl2.data$SYMBOL))
  rownames(ctrl2.data) <- names
  ctrl2.data <- ctrl2.data[,-1]
  ctrl2 <- CreateSeuratObject(counts = ctrl2.data, min.cells = 2, min.features = 200, project = "Ctrl")
  ctrl2
# Load data from file, deletes Ensembl IDs and makes gene names unique and create Seurat object for control number 3
  # Add project to differentiate between control and NTN samples
  ctrl3.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/glom_analysis/Nuria_CTF1/Raw_data/GSM4409509_control_3.txt.gz"),sep="\t", header = T, row.names = NULL)
  ctrl3.data <- subset(ctrl3.data, select = -c(IGIS))
  name <- make.unique(as.character(ctrl3.data$SYMBOL))
  rownames(ctrl3.data) <- names
  ctrl3.data <- ctrl3.data[,-1]
  ctrl3 <- CreateSeuratObject(counts = ctrl3.data, min.cells = 2, min.features = 200, project = "Ctrl")
  ctrl3
  # Merge control datasets and add cell IDs for each group
  Ctrl <- merge(ctrl2, y = ctrl3, add.cell.ids = c("Ctrl2", "Ctrl3"))
  Ctrl
# Load data from file, deletes Ensembl IDs and makes gene names unique and create Seurat object for NTN number 1 day 1
  # Add project to differentiate between control and NTN samples
  ntn1day1.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/glom_analysis/Nuria_CTF1/Raw_data/GSM4409510_nephritis_d1_1.txt.gz"),sep="\t", header = T, row.names = NULL)
  ntn1day1.data <- subset(ntn1day1.data, select = -c(IGIS))
  name <- make.unique(as.character(ntn1day1.data$SYMBOL))
  rownames(ntn1day1.data) <- names
  ntn1day1.data <- ntn1day1.data[,-1]
  ntn1day1 <- CreateSeuratObject(counts = ntn1day1.data, min.cells = 2, min.features = 200, project = "Ntnday1")
  ntn1day1
# Load data from file, deletes Ensembl IDs and makes gene names unique and create Seurat object for NTN number 2 day 1
  # Add project to differentiate between control and NTN samples
  ntn2day1.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/glom_analysis/Nuria_CTF1/Raw_data/GSM4409511_nephritis_d1_2.txt.gz"),sep="\t", header = T, row.names = NULL)
  ntn2day1.data <- subset(ntn2day1.data, select = -c(IGIS))
  name <- make.unique(as.character(ntn2day1.data$SYMBOL))
  rownames(ntn2day1.data) <- names
  ntn2day1.data <- ntn2day1.data[,-1]
  ntn2day1 <- CreateSeuratObject(counts = ntn2day1.data, min.cells = 2, min.features = 200, project = "Ntnday1")
  ntn2day1
  # Merge NTN day 1 datasets and add cell IDs for each group
  Ntnday1 <- merge(ntn1day1, y = ntn2day1, add.cell.ids = c("ntn1day1", "ntn2day1"))
  Ntnday1
# Load data from file, deletes Ensembl IDs and makes gene names unique and create Seurat object for NTN number 1 day 5
  # Add project to differentiate between control and NTN samples
  ntn1day5.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/glom_analysis/Nuria_CTF1/Raw_data/GSM4409512_nephritis_d5_1.txt.gz"),sep="\t", header = T, row.names = NULL)
  ntn1day5.data <- subset(ntn1day5.data, select = -c(IGIS))
  name <- make.unique(as.character(ntn1day5.data$SYMBOL))
  rownames(ntn1day5.data) <- names
  ntn1day5.data <- ntn1day5.data[,-1]
  ntn1day5 <- CreateSeuratObject(counts = ntn1day5.data, min.cells = 2, min.features = 200, project = "Ntnday5")
  ntn1day5
# Load data from file, deletes Ensembl IDs and makes gene names unique and create Seurat object for NTN number 2 day 5
  # Add project to differentiate between control and NTN samples
  ntn2day5.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/glom_analysis/Nuria_CTF1/Raw_data/GSM4409513_nephritis_d5_2.txt.gz"),sep="\t", header = T, row.names = NULL)
  ntn2day5.data <- subset(ntn2day5.data, select = -c(IGIS))
  name <- make.unique(as.character(ntn2day5.data$SYMBOL))
  rownames(ntn2day5.data) <- names
  ntn2day5.data <- ntn2day5.data[,-1]
  ntn2day5 <- CreateSeuratObject(counts = ntn2day5.data, min.cells = 2, min.features = 200, project = "Ntnday5")
  ntn2day5
  # Merge NTN day 5 datasets and add cell IDs for each group
  Ntnday5 <- merge(ntn1day5, y = ntn2day5, add.cell.ids = c("ntn1day5", "ntn2day5"))
  Ntnday5
# Merge Seurat objects into one to a single, unified Seurat object, and remove files not required
  Combined <- merge(Ctrl, y = c(Ntnday1, Ntnday5))
  Combined
  rm(ctrl2, ctrl2.data, ctrl3, ctrl3.data, ntn1day1, ntn1day1.data, ntn2day1, ntn2day1.data, ntn1day5, ntn1day5.data, ntn2day5, ntn2day5.data, name, names)

#----------------------------------------------------------------------------------------------------------------#

## PART 2: Combined quality control, PCA and provisional clustering of all samples prior to correction for ambient RNA contamination
# Check mitochondrial content and features per barcode for all droplets
  Combined[["percent.mt"]] <- PercentageFeatureSet(Combined, pattern = "^mt-")
  VlnPlot(Combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot1 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot1
  plot2 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot2
# Check plots above and process data to exclude cells with less than 200 or more than 6000 transcripts and less than 15% mitochondrial transcripts
  Combined <- subset(Combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
# Rerun mitochondrial content and features per barcode after QC
  Combined[["percent.mt"]] <- PercentageFeatureSet(Combined, pattern = "^mt-")
  VlnPlot(Combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot1 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot1
  plot2 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot2
# Normalise, find variable features, scale data (based on all genes) and perform PCA. Empirically determine the number of PCs to use for downstream analyses
  Combined <- NormalizeData(Combined, normalization.method = "LogNormalize", scale.factor = 10000)
  Combined <- FindVariableFeatures(Combined, selection.method = "vst", nfeatures = 2000)
  top20 <- head(VariableFeatures(Combined), 20)
  plot1 <- VariableFeaturePlot(Combined)
  plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
  plot2
  all.genes <- rownames(Combined)
  Combined <- ScaleData(Combined, features = all.genes)
  Combined <- RunPCA(Combined, features = VariableFeatures(object = Combined))
  ElbowPlot(Combined)
  pct <- Combined[["pca"]]@stdev / sum(Combined[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co1
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  co2
  pcs <- min(co1, co2)
  pcs
  rm(co1, co2, pcs, plot1, plot2)
  # KNN and Louvain clustering prior to UMAP. Normally resolution should be between 0.4-1.2 for single cell datasets
  Combined <- FindNeighbors(Combined, dims = 1:17)
  Combined <- FindClusters(Combined, resolution = 0.5)
  Combined <- RunUMAP(Combined, dims = 1:17)
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5)

#----------------------------------------------------------------------------------------------------------------#
  
## PART 3: Implement SoupX before re-normalisation, scaling and PCA
# Save metadata and create raw count matrix from Seurat object
  meta_data <- Combined$orig.ident
  rawcountmatrix <-Combined@assays$RNA@counts
# Create metadata from Seurat object
  Combined_metaData <- as.data.frame(Combined@reductions$umap@cell.embeddings)
  colnames(Combined_metaData) <- c('RD1','RD2')
  Combined_metaData$Cluster <- Combined@meta.data$seurat_clusters
# Start by generating SoupChannel object assuming no knowlege of empty droplets
  scNoDrops <- SoupChannel(rawcountmatrix, rawcountmatrix, calcSoupProfile = FALSE)
  soupProf <- data.frame(row.names = rownames(rawcountmatrix), est = rowSums(rawcountmatrix)/sum(rawcountmatrix), counts = rowSums(rawcountmatrix))
  scNoDrops <- setSoupProfile(scNoDrops, soupProf)
# Then add metadata including cluster informatiton and dimension reduction
  scNoDrops <- setClusters(scNoDrops, setNames(Combined_metaData$Cluster, rownames(Combined_metaData)))
# Add dimension reduction for the data
  scNoDrops <- setDR(scNoDrops, Combined_metaData[colnames(rawcountmatrix), c("RD1", "RD2")])
  # Estimate contamination fraction
  scNoDrops <- autoEstCont(scNoDrops)
# Generated corrected count matrix
  correctedcounts <- adjustCounts(scNoDrops)
# Create Seurat object from SoupX-corrected matrix and add groups
  Combined <- CreateSeuratObject(correctedcounts)
  Combined$group <- meta_data
# Repeat normalisation, scaling and PCA
  Combined <- NormalizeData(Combined, normalization.method = "LogNormalize", scale.factor = 10000)
  Combined <- FindVariableFeatures(Combined, selection.method = "vst", nfeatures = 2000)
  top20 <- head(VariableFeatures(Combined), 20)
  plot1 <- VariableFeaturePlot(Combined)
  plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
  plot2
  all.genes <- rownames(Combined)
  Combined <- ScaleData(Combined, features = all.genes)
  Combined <- RunPCA(Combined, features = VariableFeatures(object = Combined))
  ElbowPlot(Combined)
  pct <- Combined[["pca"]]@stdev / sum(Combined[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co1
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  co2
  pcs <- min(co1, co2)
  pcs
  rm(co1, co2, pcs, plot1, plot2, scNoDrops, soupProf, rawcountmatrix, cumu, pct, top20)

#----------------------------------------------------------------------------------------------------------------#
  
## PART 4: Data integration using Harmony, before clustering pre-doublet correction and depreciated FindAllMarkers analysis
# Create plot for PCA, grouped by "orig.ident", and compare expression levels between  datatsets
  p1 <- DimPlot(object = Combined, reduction = "pca", pt.size = .1, group.by = "orig.ident")
  p2 <- VlnPlot(object = Combined, features = "PC_1", group.by = "orig.ident", pt.size = .1)
  p1
  p2
# Run Harmony and perform batch correction by each sample, returns plot with number of iterations required for convergence
  Combined <- Combined %>% 
   RunHarmony("orig.ident", plot_convergence = TRUE)
# Access and show first five Harmony embeddings for each barcode. The title of each embedding, and thus the column names of the matrix, are harmony_x
  Combined_harmony <- Embeddings(Combined, 'harmony')
  Combined_harmony[1:5, 1:5]
  # Replot PCA and expression comparison, this time after Harmony
  p1 <- DimPlot(object = Combined, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
  p2 <- VlnPlot(object = Combined, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
  p1
  p2
  rm(p1, p2)
# Compute UMAP using Harmony embeddings by setting reduction technique to "harmony". Change number of dimensions and resolution of Louvain / Leiden clustering as desired
  Combined <- FindNeighbors(Combined, reduction = "harmony", dims = 1:17) 
  Combined <- FindClusters(Combined, resolution = 0.4)
  Combined <- RunUMAP(Combined, reduction = "harmony", dims = 1:17)
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5)
# First assessment of differentially expressed genes between clusters
  #Combined.markers <- FindAllMarkers(Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  #Combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
  #write.csv(Combined.markers, file = "PerettaTejedoretalmarkers_unannotated.csv")

#----------------------------------------------------------------------------------------------------------------#
  
## PART 5: Removal of doublets and unwanted clusters, keeping and renaming desired clusters
  # Creates 9139 artificial doublets, merges with dataset, finds artificial k nearest neighbors and estimates doublets in real data.
  # pK value of 0.005 based on preliminary analysis using pararmSweep_v3, summarizeSweep and find.pK functions.
  sweep.res.list_kidney <- paramSweep_v3(Combined, PCs = 1:17, sct = FALSE)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  annotations <- Combined@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           
  nExp_poi <- round(0.05*nrow(Combined@meta.data)) # alternatively use nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  Combined <- doubletFinder_v3(Combined, PCs = 1:17, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "DF.classifications_0.25_0.005_1371")  # amend DF.classifications as required
  # Makes a seperate Seurat object of singlets and check UMAP
  Combined <- subset(Combined, subset = DF.classifications_0.25_0.005_1371 == "Singlet") # amend DF.classifications as required
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5)
  # Remove unwanted clusters including afferent/efferent arteriole (4), TECs (5, 15), proliferating cells (7, 14), smooth muscle cells (8), mesangial doublets (11), rbcs (19)
  Combined <- subset(Combined, idents = c("0", "1", "2", "3", "6", "9", "10", "12"))
  new.cluster.ids <- c("GEC",  # 0
                       "Mesangium",  # 1
                       "Podo",  # 2
                       "GEC",  # 3
                       "Immune",  # 6
                       "PEC",  # 9 
                       "GEC",  # 10 
                       "Immune"  # 12
  )
  names(new.cluster.ids) <- levels(Combined)
  Combined <- RenameIdents(Combined, new.cluster.ids)
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5)
  # Calculate top 10 DEGs per cluster and save as CSV file in the current directory
  #Combined.markers <- FindAllMarkers(Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  #Combined.markers %>% group_by(cluster) %>% top_n(n = 10)
  #top10 <- Combined.markers %>% group_by(cluster) %>% top_n(n = 10)
  #write.csv(Combined.markers,"PerettaTejedoretalmarkers_annotated.csv", row.names = FALSE)
  
#----------------------------------------------------------------------------------------------------------------#

## PART 6: Depreciated demonstration of clustering and cell-type assignment
# UMAPs based on cell type and group
#DimPlot(Combined, reduction = "umap", label = F, pt.size = 1)
#DimPlot(Combined, reduction = "umap", label = F, pt.size = 1, split.by = "group")
# Cell type identification / demonstration
#Podocyte_label <- WhichCells(Combined, idents = c("Podocyte"))
#GEC_label <- WhichCells(Combined, idents = c("GEC"))
#Mesangium_label <- WhichCells(Combined, idents = c("Mesangium"))
#PEC_label <- WhichCells(Combined, idents = c("PEC"))
#DimPlot(Combined, reduction = "umap", label = F, pt.size = 1, cells.highlight = Podocyte_label)
#DimPlot(Combined, reduction = "umap", label = F, pt.size = 1, cells.highlight = GEC_label)
#DimPlot(Combined, reduction = "umap", label = F, pt.size = 1, cells.highlight = Mesangium_label)
#DimPlot(Combined, reduction = "umap", label = F, pt.size = 1, cells.highlight = PEC_label)
#FeaturePlot(Combined, features = c("Wt1", "Nphs2"))
#FeaturePlot(Combined, features = c("Emcn", "Kdr"))
#FeaturePlot(Combined, features = c("Ptn", "Pdgfrb"))
#FeaturePlot(Combined, features = c("Pax8", "Jag1"))
  
#----------------------------------------------------------------------------------------------------------------#

## PART 7: Analysis of CTF1 expression by cell type and comparison between conditions
DimPlot(Combined, reduction = "umap", label = F, pt.size = 0.5) #Fig1h
#DimPlot(Combined, reduction = "umap", label = F, pt.size = 0.5, split.by = "group") #Fig1h
FeaturePlot(Combined, features = "Ctf1", pt.size = 0.5, order = T, min.cutoff = 0.5) #Fig1i
VlnPlot(Combined, features = "Ctf1", split.by = "group", idents = "PEC") #Fig1j

# Draw violinplots for markers of interest for each cell type in the glomerulus
#VlnPlot(Combined, features = "Ctf1", split.by = "group", idents = "Podo")
#VlnPlot(Combined, features = "Il6st", split.by = "group", idents = "Podo")
#VlnPlot(Combined, features = "Lifr", split.by = "group", idents = "Podo")
#VlnPlot(Combined, features = "Ctf1", split.by = "group", idents = "GEC")
#VlnPlot(Combined, features = "Il6st", split.by = "group", idents = "GEC")
#VlnPlot(Combined, features = "Lifr", split.by = "group", idents = "GEC")
#VlnPlot(Combined, features = "Ctf1", split.by = "group", idents = "Mesangium")
#VlnPlot(Combined, features = "Il6st", split.by = "group", idents = "Mesangium")
#VlnPlot(Combined, features = "Lifr", split.by = "group", idents = "Mesangium")
#VlnPlot(Combined, features = "Ctf1", split.by = "group", idents = "PEC")

# Use FindAllMarkers function to compare cell type-specific CTF1 between conditions
PECdata <- subset(Combined, idents = "PEC")
Idents(PECdata) <- PECdata@meta.data$group
VlnPlot(PECdata, features = "Ctf1", split.by = "group", pt.size = 0) #Fig1j
compare_ctf1_PEC_1 <- FindMarkers(PECdata, ident.1 = "Ctrl", ident.2 = "Ntnday1", features = "Ctf1", logfc.threshold = 0)
compare_ctf1_PEC_1
compare_ctf1_PEC_2 <- FindMarkers(PECdata, ident.1 = "Ctrl", ident.2 = "Ntnday5", features = "Ctf1", logfc.threshold = 0)
compare_ctf1_PEC_2
compare_ctf1_PEC_3 <- FindMarkers(PECdata, ident.1 = "Ntnday1", ident.2 = "Ntnday5", features = "Ctf1", logfc.threshold = 0)
compare_ctf1_PEC_3

#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
