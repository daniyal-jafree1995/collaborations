## scRNA-seq analysis for Dr Jennie Chandler
## Cross-disease comparison of podocytes in diabetic kidney disease (DKD), nephrotoxic nephritis (NTN) and CD2AP null mice
## Data derived from: https://jasn.asnjournals.org/content/early/2020/07/10/ASN.2020020220
## Authors: Daniyal Jafree (UCL)
## Version 1: 08/04/2022
## Last updated: 03/10/2022
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
  setwd("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/Data/JASN")

#----------------------------------------------------------------------------------------------------------------#
## PART 1: Script for loading individual samples into Seurat objects
  # Load data from file, deletes Ensembl IDs and makes gene names unique and create Seurat object for each control, before merging. Add project to differentiate between samples
  # Ctrl2
  ctrl2.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/Data/JASN/GSE146912_RAW/GSM4409508_control_2.txt.gz"),sep="\t", header = T, row.names = NULL)
  ctrl2.data <- subset(ctrl2.data, select = -c(IGIS))
  names <- make.unique(as.character(ctrl2.data$SYMBOL))
  rownames(ctrl2.data) <- names
  ctrl2.data <- ctrl2.data[,-1]
  ctrl2 <- CreateSeuratObject(counts = ctrl2.data, min.cells = 2, min.features = 200, project = "Ctrl_ntn")
  ctrl2
  # Ctrl3
  ctrl3.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/Data/JASN/GSE146912_RAW/GSM4409509_control_3.txt.gz"),sep="\t", header = T, row.names = NULL)
  ctrl3.data <- subset(ctrl3.data, select = -c(IGIS))
  names <- make.unique(as.character(ctrl3.data$SYMBOL))
  rownames(ctrl3.data) <- names
  ctrl3.data <- ctrl3.data[,-1]
  ctrl3 <- CreateSeuratObject(counts = ctrl3.data, min.cells = 2, min.features = 200, project = "Ctrl_ntn")
  ctrl3
  # Merge control datasets and add cell IDs for each group
  Ctrl_ntn <- merge(ctrl2, y = c(ctrl3), add.cell.ids = c("Ctrl2", "Ctrl3"))
  Ctrl_ntn
  rm(ctrl2, ctrl2.data, ctrl3, ctrl3.data, names)
  # NTN day 5 repeat 1
  ntn1day5.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/Data/JASN/GSE146912_RAW//GSM4409512_nephritis_d5_1.txt.gz"),sep="\t", header = T, row.names = NULL)
  ntn1day5.data <- subset(ntn1day5.data, select = -c(IGIS))
  names <- make.unique(as.character(ntn1day5.data$SYMBOL))
  rownames(ntn1day5.data) <- names
  ntn1day5.data <- ntn1day5.data[,-1]
  ntn1day5 <- CreateSeuratObject(counts = ntn1day5.data, min.cells = 2, min.features = 200, project = "NTN")
  ntn1day5
  # NTN day 5 repeat 2
  ntn2day5.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/Data/JASN/GSE146912_RAW/GSM4409513_nephritis_d5_2.txt.gz"),sep="\t", header = T, row.names = NULL)
  ntn2day5.data <- subset(ntn2day5.data, select = -c(IGIS))
  names <- make.unique(as.character(ntn2day5.data$SYMBOL))
  rownames(ntn2day5.data) <- names
  ntn2day5.data <- ntn2day5.data[,-1]
  ntn2day5 <- CreateSeuratObject(counts = ntn2day5.data, min.cells = 2, min.features = 200, project = "NTN")
  ntn2day5
  # Merge NTN day 5 datasets and add cell IDs for each group
  NTN <- merge(ntn1day5, y = ntn2day5, add.cell.ids = c("ntn1day5", "ntn2day5"))
  NTN
  rm(ntn1day5.data, ntn1day5, ntn2day5.data, ntn2day5, names)
  Combined <- merge(Ctrl_ntn, y = c(NTN))
  Combined
  #----------------------------------------------------------------------------------------------------------------#
  ## PART 2: Quality control
  # Calculate percentage of reads mapping to mitochondrial genome
  Combined[["percent.mt"]] <- PercentageFeatureSet(Combined, pattern="^mt-")
  # subset cells by mitochondrial contact and number of features
  Combined <- subset(Combined, subset = nFeature_RNA > 200 & percent.mt < 10 & nFeature_RNA < 7500) # can optionally assign upper cutoff for nFeature_RNA ?7500
  Combined
  print(VlnPlot(Combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident",ncol = 3))
  plot1 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")
  plot2 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")
  print(plot1 + plot2)
  #----------------------------------------------------------------------------------------------------------------#
  ## PART 3: Normalisation, scaling and PCA with script to empirically decide number of PCs to carry forward (16 according to preliminary analysis)
  # Normalisation of RNA counts
  Combined <- NormalizeData(Combined, normalization.method = "LogNormalize", scale.factor = 10000)
  Combined <- FindVariableFeatures(Combined, selection.method = "vst", nfeatures = 2000)
  # Scaling
  all_genes <- rownames(Combined)
  Combined <- ScaleData(Combined,features = all_genes)
  # PCA
  Combined <- RunPCA(Combined, features = VariableFeatures(object = Combined))
  VizDimLoadings(Combined, dims = 1:2, reduction = "pca")
  DimPlot(Combined, reduction = "pca")
  # Empirically determine number of PCs to use for downstream analysis using elbow plot and script to calculate pcs
  ElbowPlot(Combined, ndims = 50)
  pct <- Combined[["pca"]]@stdev / sum(Combined[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co1
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  co2
  pcs <- min(co1, co2)
  pcs
  #----------------------------------------------------------------------------------------------------------------#
  ## PART 4: Unsupervised clustering using UMAP. Some commented-out options below that Gideon prefers to use.
  # KNN and Louvain clustering prior to UMAP. Normally resolution should be between 0.4-1.2 for single cell datasets
  Combined <- FindNeighbors(Combined, dims = 1:19)
  Combined <- FindClusters(Combined, resolution = 0.4) # alternative Leiden approach from Gideon: FindClusters(wt1, resolution = 0.9, method = "igraph", algorithm = 4)
  # Run and plot UMAP by cell types
  Combined <- RunUMAP(Combined, dims = 1:19)
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5) # can use group.by and call variable from object@meta.data to group cells in UMAP by variable
  FeaturePlot(Combined, features = "Nphs2")
  #----------------------------------------------------------------------------------------------------------------#
  ## PART 5: Implement SoupX package for estimation and correction of raw count matrix for contaminating RNA. Then run Parts 1-4 again with Harmony
  # Create raw count matrix from Seurat object
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
  # Estimate contamination fraction utomatically
  scNoDrops <- autoEstCont(scNoDrops)
  # Generated corrected count matrix
  correctedcounts <- adjustCounts(scNoDrops)
  plotChangeMap(scNoDrops, correctedcounts, "Nphs2")  # Can optionally examine how SoupX has changed expression values for particular genes
  # Make new Seurat object after SoupX correction and carry forward metadata from Part 1
  Combined <- CreateSeuratObject(correctedcounts, meta.data = meta_data) #N.B. metadata is re-used from earlier!
  rm(Combined_metaData, plot1, plot2, rawcountmatrix, scNoDrops, soupProf, all_genes, co1, co2, cumu, pcs, pct, top20)
  # Essential scripts from Parts 1-4. Up to the completion of PCA.
  Combined <- NormalizeData(Combined, normalization.method = "LogNormalize", scale.factor = 10000)
  Combined <- FindVariableFeatures(Combined, selection.method = "vst", nfeatures = 2000)
  top20 <- head(VariableFeatures(Combined), 20)
  plot1 <- VariableFeaturePlot(Combined)
  plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
  plot2
  all.genes <- rownames(Combined)
  Combined <- ScaleData(Combined, features = all.genes)
  Combined <- RunPCA(Combined, features = VariableFeatures(object = Combined))
  ## Harmony script
  # Run Harmony and perform batch correction, returns plot with number of iterations required for convergence.
  Combined <- Combined %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)
  # Access and show first five Harmony embeddings for each barcode.
  Combined_harmony <- Embeddings(Combined, 'harmony')
  Combined_harmony[1:5, 1:5]
  # Plot PCA and expression comparison after Harmony
  p1 <- DimPlot(object = Combined, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
  p2 <- VlnPlot(object = Combined, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
  p1
  p2
  rm(p1, p2)
  #----------------------------------------------------------------------------------------------------------------#
  ## PART 6: UMAP post-SoupX (+/- Harmony), identify and extract podocytes from all datasets, before differential expression testing and outputs for master script
  # Compute UMAP
  Combined <- FindNeighbors(Combined, reduction = "harmony", dims = 1:19) # use Harmony embeddings by setting reduction technique to "harmony"
  Combined <- FindClusters(Combined, resolution = 0.5)
  Combined <- RunUMAP(Combined, reduction = "harmony", dims = 1:19) # use Harmony embeddings by setting reduction technique to "harmony"
  # Generate UMAP and group by desired variables
  DimPlot(Combined, reduction = "umap", label = F, pt.size = 1)
  FeaturePlot(Combined, features = c("Nphs1", "Nphs2"))
  # Downstream applications...
    # Podocyte subclustering and examining podocyte DEGs
      #NTN_comparison <- subset(Combined, idents = c("3", "15"))
      #Idents(NTN_comparison) <- NTN_comparison@meta.data$orig.ident
      #DefaultAssay(NTN_comparison) <- "RNA"
      #NTN_comparison.Podo.DE <- FindAllMarkers(NTN_comparison)
      #NTN_comparison.Podo.DE
      #write.csv(NTN_comparison.Podo.DE, "NTN_comparison.Podo.DE.csv")
  
    # Generation of psuedobulk dataset
      #NTN_PSEUDOBULK <- subset(Combined, subset = orig.ident == "NTN")
      #Idents(NTN_PSEUDOBULK) <- NTN_PSEUDOBULK@meta.data$orig.ident
      #DefaultAssay(NTN_PSEUDOBULK) <- "RNA"
  
  #----------------------------------------------------------------------------------------------------------------#
  #REPEAT FOR DKD. Clear objects from environment prior to running.
  # BTBR Ctrl2 - AKA Ctrl4
  ctrl4.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/Data/JASN/GSE146912_RAW/GSM4409519_BTBR_control_2.txt.gz"),sep="\t", header = T, row.names = NULL)
  ctrl4.data <- subset(ctrl4.data, select = -c(IGIS))
  names <- make.unique(as.character(ctrl4.data$SYMBOL))
  rownames(ctrl4.data) <- names
  ctrl4.data <- ctrl4.data[,-1]
  ctrl4 <- CreateSeuratObject(counts = ctrl4.data, min.cells = 2, min.features = 200, project = "Ctrl_dkd")
  ctrl4
  # BTBR Ctrl1 - AKA Ctrl1
  ctrl1.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/Data/JASN/GSE146912_RAW/GSM4409518_BTBR_control_1.txt.gz"),sep="\t", header = T, row.names = NULL)
  ctrl1.data <- subset(ctrl1.data, select = -c(IGIS))
  names <- make.unique(as.character(ctrl1.data$SYMBOL))
  rownames(ctrl1.data) <- names
  ctrl1.data <- ctrl1.data[,-1]
  ctrl1 <- CreateSeuratObject(counts = ctrl1.data, min.cells = 2, min.features = 200, project = "Ctrl_dkd")
  ctrl1
  # Merge control datasets and add cell IDs for each group
  Ctrl_dkd <- merge(ctrl1, y = c(ctrl4), add.cell.ids = c("Ctrl1", "Ctrl4"))
  Ctrl_dkd
  rm(ctrl1, ctrl1.data, ctrl4, ctrl4.data, names)
  # DKD week 12 repeat 1
  DKD1.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/Data/JASN/GSE146912_RAW/GSM4409520_BTBR_wk12_1.txt.gz"),sep="\t", header = T, row.names = NULL)
  DKD1.data <- subset(DKD1.data, select = -c(IGIS))
  names <- make.unique(as.character(DKD1.data$SYMBOL))
  rownames(DKD1.data) <- names
  DKD1.data <- DKD1.data[,-1]
  DKD1 <- CreateSeuratObject(counts = DKD1.data, min.cells = 2, min.features = 200, project = "DKD")
  DKD1
  # DKD week 12 repeat 2
  DKD2.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/Data/JASN/GSE146912_RAW/GSM4409521_BTBR_wk12_2.txt.gz"),sep="\t", header = T, row.names = NULL)
  DKD2.data <- subset(DKD2.data, select = -c(IGIS))
  names <- make.unique(as.character(DKD2.data$SYMBOL))
  rownames(DKD2.data) <- names
  DKD2.data <- DKD2.data[,-1]
  DKD2 <- CreateSeuratObject(counts = DKD2.data, min.cells = 2, min.features = 200, project = "DKD")
  DKD2
  # Merge DKD week 12 datasets and add cell IDs for each group
  DKD <- merge(DKD1, y = DKD2, add.cell.ids = c("DKD1", "DKD2"))
  DKD
  rm(DKD1.data, DKD1, DKD2.data, DKD2, names)
  Combined <- merge(Ctrl_dkd, y = c(DKD))
  Combined
  #----------------------------------------------------------------------------------------------------------------#
  ## PART 2: Quality control
  # Calculate percentage of reads mapping to mitochondrial genome
  Combined[["percent.mt"]] <- PercentageFeatureSet(Combined, pattern="^mt-")
  # subset cells by mitochondrial contact and number of features
  Combined <- subset(Combined, subset = nFeature_RNA > 200 & percent.mt < 10 & nFeature_RNA < 6000) # can optionally assign upper cutoff for nFeature_RNA ?7500
  Combined
  print(VlnPlot(Combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident",ncol = 3))
  plot1 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")
  plot2 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")
  print(plot1 + plot2)
  #----------------------------------------------------------------------------------------------------------------#
  ## PART 3: Normalisation, scaling and PCA with script to empirically decide number of PCs to carry forward (16 according to preliminary analysis)
  # Normalisation of RNA counts
  Combined <- NormalizeData(Combined, normalization.method = "LogNormalize", scale.factor = 10000)
  Combined <- FindVariableFeatures(Combined, selection.method = "vst", nfeatures = 2000)
  # Scaling
  all_genes <- rownames(Combined)
  Combined <- ScaleData(Combined,features = all_genes)
  # PCA
  Combined <- RunPCA(Combined, features = VariableFeatures(object = Combined))
  VizDimLoadings(Combined, dims = 1:2, reduction = "pca")
  DimPlot(Combined, reduction = "pca")
  # Empirically determine number of PCs to use for downstream analysis using elbow plot and script to calculate pcs
  ElbowPlot(Combined, ndims = 50)
  pct <- Combined[["pca"]]@stdev / sum(Combined[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co1
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  co2
  pcs <- min(co1, co2)
  pcs
  #----------------------------------------------------------------------------------------------------------------#
  ## PART 4: Unsupervised clustering using UMAP. Some commented-out options below that Gideon prefers to use.
  # KNN and Louvain clustering prior to UMAP. Normally resolution should be between 0.4-1.2 for single cell datasets
  Combined <- FindNeighbors(Combined, dims = 1:22)
  Combined <- FindClusters(Combined, resolution = 0.4) # alternative Leiden approach from Gideon: FindClusters(wt1, resolution = 0.9, method = "igraph", algorithm = 4)
  # Run and plot UMAP by cell types
  Combined <- RunUMAP(Combined, dims = 1:22)
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5) # can use group.by and call variable from object@meta.data to group cells in UMAP by variable
  FeaturePlot(Combined, features = "Nphs2")
  #----------------------------------------------------------------------------------------------------------------#
  ## PART 5: Implement SoupX package for estimation and correction of raw count matrix for contaminating RNA. Then run Parts 1-4 again with Harmony
  # Create raw count matrix from Seurat object
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
  # Estimate contamination fraction utomatically
  scNoDrops <- autoEstCont(scNoDrops)
  # Generated corrected count matrix
  correctedcounts <- adjustCounts(scNoDrops)
  plotChangeMap(scNoDrops, correctedcounts, "Nphs2")  # Can optionally examine how SoupX has changed expression values for particular genes
  # Make new Seurat object after SoupX correction and carry forward metadata from Part 1
  Combined <- CreateSeuratObject(correctedcounts, meta.data = meta_data) #N.B. metadata is re-used from earlier!
  rm(Combined_metaData, plot1, plot2, rawcountmatrix, scNoDrops, soupProf, all_genes, co1, co2, cumu, pcs, pct, top20)
  # Essential scripts from Parts 1-4. Up to the completion of PCA.
  Combined <- NormalizeData(Combined, normalization.method = "LogNormalize", scale.factor = 10000)
  Combined <- FindVariableFeatures(Combined, selection.method = "vst", nfeatures = 2000)
  top20 <- head(VariableFeatures(Combined), 20)
  plot1 <- VariableFeaturePlot(Combined)
  plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
  plot2
  all.genes <- rownames(Combined)
  Combined <- ScaleData(Combined, features = all.genes)
  Combined <- RunPCA(Combined, features = VariableFeatures(object = Combined))
  ## Harmony script
  # Run Harmony and perform batch correction, returns plot with number of iterations required for convergence.
  Combined <- Combined %>%
    RunHarmony("orig.ident", plot_convergence = TRUE)
  # Access and show first five Harmony embeddings for each barcode.
  Combined_harmony <- Embeddings(Combined, 'harmony')
  Combined_harmony[1:5, 1:5]
  # Plot PCA and expression comparison after Harmony
  p1 <- DimPlot(object = Combined, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
  p2 <- VlnPlot(object = Combined, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
  p1
  p2
  rm(p1, p2)
  #----------------------------------------------------------------------------------------------------------------#
  ## PART 6: UMAP post-SoupX (+/- Harmony), identify and extract podocytes from all datasets, before differential expression testing and outputs for master script
  # Compute UMAP
  Combined <- FindNeighbors(Combined, reduction = "harmony", dims = 1:22) # use Harmony embeddings by setting reduction technique to "harmony"
  Combined <- FindClusters(Combined, resolution = 0.5)
  Combined <- RunUMAP(Combined, reduction = "harmony", dims = 1:22) # use Harmony embeddings by setting reduction technique to "harmony"
  # Generate UMAP and group by desired variables
  DimPlot(Combined, reduction = "umap", label = F, pt.size = 1)
  FeaturePlot(Combined, features = c("Nphs1", "Nphs2"))
  # Downstream applications...
    # Podocyte subclustering and examining podocyte DEGs
      #DKD_comparison <- subset(Combined, idents = c("3"))
      #Idents(DKD_comparison) <- DKD_comparison@meta.data$orig.ident
      #DefaultAssay(DKD_comparison) <- "RNA"
      #DKD_comparison.Podo.DE <- FindAllMarkers(DKD_comparison)
      #DKD_comparison.Podo.DE
      #write.csv(DKD_comparison.Podo.DE, "DKD_comparison.Podo.DE.csv")
  
    # Generation of psuedobulk dataset
      #DKD_PSEUDOBULK <- subset(Combined, subset = orig.ident == "DKD")
      #Idents(DKD_PSEUDOBULK) <- DKD_PSEUDOBULK@meta.data$orig.ident
      #DefaultAssay(DKD_PSEUDOBULK) <- "RNA"
  
  #----------------------------------------------------------------------------------------------------------------#
  #REPEAT FOR CD2AP null. Clear objects from environment prior to running.
  # Ctrl_cd2ap
  ctrl5.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/Data/JASN/GSE146912_RAW/GSM4409516_CD2AP_WT.txt.gz"),sep="\t", header = T, row.names = NULL)
  ctrl5.data <- subset(ctrl5.data, select = -c(IGIS))
  names <- make.unique(as.character(ctrl5.data$SYMBOL))
  rownames(ctrl5.data) <- names
  ctrl5.data <- ctrl5.data[,-1]
  Ctrl_cd2ap <- CreateSeuratObject(counts = ctrl5.data, min.cells = 2, min.features = 200, project = "Ctrl_cd2ap", add.cell.ids = c("Ctrl5"))
  Ctrl_cd2ap
  rm(ctrl5.data, names)
  # CD2AP mutant
  CD2AP1.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/Data/JASN/GSE146912_RAW/GSM4409517_CD2AP_KO.txt.gz"),sep="\t", header = T, row.names = NULL)
  CD2AP1.data <- subset(CD2AP1.data, select = -c(IGIS))
  names <- make.unique(as.character(CD2AP1.data$SYMBOL))
  rownames(CD2AP1.data) <- names
  CD2AP1.data <- CD2AP1.data[,-1]
  CD2AP <- CreateSeuratObject(counts = CD2AP1.data, min.cells = 2, min.features = 200, project = "CD2AP", add.cell.ids = c("CD2AP"))
  CD2AP
  rm(CD2AP1.data, names)
  # Merge Seurat objects into one to a single, unified Seurat object, and remove files not required
  Combined <- merge(Ctrl_cd2ap, y = c(CD2AP))
  Combined
  #----------------------------------------------------------------------------------------------------------------#
  ## PART 2: Quality control
  # Calculate percentage of reads mapping to mitochondrial genome
  Combined[["percent.mt"]] <- PercentageFeatureSet(Combined, pattern="^mt-")
  # subset cells by mitochondrial contact and number of features
  Combined <- subset(Combined, subset = nFeature_RNA > 200 & percent.mt < 1 & nFeature_RNA < 5000) # can optionally assign upper cutoff for nFeature_RNA ?7500
  Combined
  print(VlnPlot(Combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident",ncol = 3))
  plot1 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")
  plot2 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")
  print(plot1 + plot2)
  #----------------------------------------------------------------------------------------------------------------#
  ## PART 3: Normalisation, scaling and PCA with script to empirically decide number of PCs to carry forward (16 according to preliminary analysis)
  # Normalisation of RNA counts
  Combined <- NormalizeData(Combined, normalization.method = "LogNormalize", scale.factor = 10000)
  Combined <- FindVariableFeatures(Combined, selection.method = "vst", nfeatures = 2000)
  # Scaling
  all_genes <- rownames(Combined)
  Combined <- ScaleData(Combined,features = all_genes)
  # PCA
  Combined <- RunPCA(Combined, features = VariableFeatures(object = Combined))
  VizDimLoadings(Combined, dims = 1:2, reduction = "pca")
  DimPlot(Combined, reduction = "pca")
  # Empirically determine number of PCs to use for downstream analysis using elbow plot and script to calculate pcs
  ElbowPlot(Combined, ndims = 50)
  pct <- Combined[["pca"]]@stdev / sum(Combined[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co1
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  co2
  pcs <- min(co1, co2)
  pcs
  #----------------------------------------------------------------------------------------------------------------#
  ## PART 4: Unsupervised clustering using UMAP. Some commented-out options below that Gideon prefers to use.
  # KNN and Louvain clustering prior to UMAP. Normally resolution should be between 0.4-1.2 for single cell datasets
  Combined <- FindNeighbors(Combined, dims = 1:16)
  Combined <- FindClusters(Combined, resolution = 0.4) # alternative Leiden approach from Gideon: FindClusters(wt1, resolution = 0.9, method = "igraph", algorithm = 4)
  # Run and plot UMAP by cell types
  Combined <- RunUMAP(Combined, dims = 1:16)
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5) # can use group.by and call variable from object@meta.data to group cells in UMAP by variable
  FeaturePlot(Combined, features = "Nphs2")
  #----------------------------------------------------------------------------------------------------------------#
  ## PART 5: Implement SoupX package for estimation and correction of raw count matrix for contaminating RNA. Then run Parts 1-4 again with Harmony
  # Create raw count matrix from Seurat object
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
  # Estimate contamination fraction utomatically
  scNoDrops <- autoEstCont(scNoDrops)
  # Generated corrected count matrix
  correctedcounts <- adjustCounts(scNoDrops)
  plotChangeMap(scNoDrops, correctedcounts, "Nphs2")  # Can optionally examine how SoupX has changed expression values for particular genes
  # Make new Seurat object after SoupX correction and carry forward metadata from Part 1
  Combined <- CreateSeuratObject(correctedcounts, meta.data = meta_data) #N.B. metadata is re-used from earlier!
  rm(Combined_metaData, plot1, plot2, rawcountmatrix, scNoDrops, soupProf, all_genes, co1, co2, cumu, pcs, pct, top20)
  # Essential scripts from Parts 1-4. Up to the completion of PCA.
  Combined <- NormalizeData(Combined, normalization.method = "LogNormalize", scale.factor = 10000)
  Combined <- FindVariableFeatures(Combined, selection.method = "vst", nfeatures = 2000)
  top20 <- head(VariableFeatures(Combined), 20)
  plot1 <- VariableFeaturePlot(Combined)
  plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
  plot2
  all.genes <- rownames(Combined)
  Combined <- ScaleData(Combined, features = all.genes)
  Combined <- RunPCA(Combined, features = VariableFeatures(object = Combined))
  ## Harmony script
  # Run Harmony and perform batch correction, returns plot with number of iterations required for convergence.
  Combined <- Combined %>%
    RunHarmony("orig.ident", plot_convergence = TRUE)
  # Access and show first five Harmony embeddings for each barcode.
  Combined_harmony <- Embeddings(Combined, 'harmony')
  Combined_harmony[1:5, 1:5]
  # Plot PCA and expression comparison after Harmony
  p1 <- DimPlot(object = Combined, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
  p2 <- VlnPlot(object = Combined, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
  p1
  p2
  rm(p1, p2)
  #----------------------------------------------------------------------------------------------------------------#
  ## PART 6: UMAP post-SoupX (+/- Harmony), identify and extract podocytes from all datasets, before differential expression testing and outputs for master script
  # Compute UMAP
  Combined <- FindNeighbors(Combined, reduction = "harmony", dims = 1:16) # use Harmony embeddings by setting reduction technique to "harmony"
  Combined <- FindClusters(Combined, resolution = 0.5)
  Combined <- RunUMAP(Combined, reduction = "harmony", dims = 1:16) # use Harmony embeddings by setting reduction technique to "harmony"
  # Generate UMAP and group by desired variables
  DimPlot(Combined, reduction = "umap", label = F, pt.size = 1)
  FeaturePlot(Combined, features = c("Nphs1", "Nphs2"))
  # Downstream applications...
    # Podocyte subclustering and examining podocyte DEGs
      #CD2AP_comparison <- subset(Combined, idents = c("1"))
      #Idents(CD2AP_comparison) <- CD2AP_comparison@meta.data$orig.ident
      #DefaultAssay(CD2AP_comparison) <- "RNA"
      #CD2AP_comparison.Podo.DE <- FindAllMarkers(CD2AP_comparison)
      #CD2AP_comparison.Podo.DE
      #write.csv(CD2AP_comparison.Podo.DE, "CD2AP_comparison.Podo.DE.csv")
  
    # Generation of psuedobulk dataset
      #CD2AP_PSEUDOBULK <- subset(Combined, subset = orig.ident == "CD2AP")
      #Idents(CD2AP_PSEUDOBULK) <- CD2AP_PSEUDOBULK@meta.data$orig.ident
      #DefaultAssay(CD2AP_PSEUDOBULK) <- "RNA"
  #----------------------------------------------------------------------------------------------------------------#
  ## FINAL STEP: Load in WT1 glomerulopathy data, generate pseudobulk and perform differential expression across datasets
  Combined <- readRDS("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/WT1_scRNAseq.rds")
  WT1_PSEUDOBULK <- subset(Combined, subset = Conditions == "hetR394W")
  Idents(WT1_PSEUDOBULK) <- WT1_PSEUDOBULK@meta.data$orig.ident
  DefaultAssay(WT1_PSEUDOBULK) <- "RNA"
  Aggregate_Pseudobulk <- merge(WT1_PSEUDOBULK, y = c(CD2AP_PSEUDOBULK, DKD_PSEUDOBULK, NTN_PSEUDOBULK))
  #Aggregate_Pseudobulk <- readRDS("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/PSEUDOBULK_ANALYSIS/Aggregate_Pseudobulk.rds")
  compare_Pseudobulk <- FindAllMarkers(Aggregate_Pseudobulk, min.pct = 0.05, logfc.threshold = 0.5, only.pos = T)
  compare_Pseudobulk
  #write.csv(compare_Pseudobulk, "compare_Pseudobulk.csv")
  DotPlot(Aggregate_Pseudobulk, features = c("Cxcl1", "Cxcl2", "Cxcl13", "Il1b", "Tnfsf9", "Fcer1g", "C1qb"))
  #saveRDS(Aggregate_Pseudobulk, file = "/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/PSEUDOBULK_ANALYSIS/Aggregate_Pseudobulk.rds")
  
  #--------------------------------------------------### END ###--------------------------------------------------#
