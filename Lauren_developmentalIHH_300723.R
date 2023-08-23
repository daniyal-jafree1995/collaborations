# Analysis of E17.5 and E18.5 mouse kidney scRNAseq for Ihh expressing cells
# E17.5 data from lab's own scRNA-seq data
# E18.5 data from https://pubmed.ncbi.nlm.nih.gov/31118232/ 
# Daniyal Jafree | 30th July 2023 | Version 1 | Harmony integration of both datasets

# Download the packages required for analysis of the data. We make use of Seurat v3 by the Sajita lab. Worked tutorials for the use of Seurat are available at https://satijalab.org/seurat/vignettes.html
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(Matrix)
library(ggrepel)
library(patchwork)
library(tidyselect)
library(harmony)
library(SoupX)

#----------------------------------------------------------------------------------------------------------------#

# Set working directory, datasets and create Seurat object
  setwd("~/Desktop/scRNAseq_analyses/scRNAseq_trial/Post_Cell_Ranger/fetal_filtered_feature_bc_matrix")
  e17.data <- Read10X(data.dir = "Compressed")
  e17 <- CreateSeuratObject(counts = e17.data, project = "e17", min.cells = 3, min.features = 200)
  e17
# Re-set working directory, load E18.5 data from Melissa little laboratory and create Seurat object
  setwd("/Users/daniyaljafree/Desktop/Collaborations/Macrophage_kidneydev/Laura_scRNAseq/Mouse/Raw_data/Combes_data_GSE108291")
  e18.data <- Read10X(data.dir = "GSE108291")
  e18 <- CreateSeuratObject(counts = e18.data, project = "e18", min.cells = 3, min.features = 200)
  e18
# Create merged Seurat object with two datasets
  Combined <- merge(e17, y = e18, add.cell.ids = c("e17", "e18"))
  Combined
  
#----------------------------------------------------------------------------------------------------------------#
  
# Create feature and mitochondrial plots before filtering. Change between 'MT' and 'mt' if using human or mouse cells respectively
  Combined[["percent.mt"]] <- PercentageFeatureSet(Combined, pattern = "^mt-")
  VlnPlot(Combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot1 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot1
  plot2 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot2
  # Process data to exclude cells with less than 200 or more than 5000 transcripts and less than 15% mitochondrial transcripts. This is based on examination of the plots above.
  Combined <- subset(Combined, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
  # Displays feature and mitochondrial plots after filtering
  plot1 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot1
  plot2 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot2
  
#----------------------------------------------------------------------------------------------------------------#
  
# Normalisation of RNA counts, scaling and PCA. Determine number of PCs to use for downstream analysis
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
  
# Run Harmony and perform batch integration, returns plot with number of iterations required for convergence.
  Combined <- Combined %>%
  RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)
  # Access and show first five Harmony embeddings for each barcode.
  Combined_harmony <- Embeddings(Combined, 'harmony')
  Combined_harmony[1:5, 1:5]
  
#----------------------------------------------------------------------------------------------------------------#
  
# UMAP post-Harmony
  # Compute UMAP
  Combined <- FindNeighbors(Combined, dims = 1:18, reduction = "harmony") # use Harmony embeddings by setting reduction technique to "harmony"
  Combined <- FindClusters(Combined, resolution = 0.6)
  Combined <- RunUMAP(Combined, dims = 1:18, reduction = "harmony") # use Harmony embeddings by setting reduction technique to "harmony"
  # Plot UMAP and examine Ihh
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5)
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "orig.ident")
  FeaturePlot(Combined, features = "Ihh", order = T, pt.size = 0.5)
  
#----------------------------------------------------------------------------------------------------------------#
  
# Pull metadata from Seurat object and create raw count matrix from it
  metadata <- Combined@meta.data
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
  Combined <- CreateSeuratObject(correctedcounts, meta.data = metadata) #N.B. metadata is re-used from earlier!
  
  #----------------------------------------------------------------------------------------------------------------#
  
# POST-SOUPX - Normalisation of RNA counts, scaling and PCA. Determine number of PCs to use for downstream analysis
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
  
# Run Harmony and perform batch integration, returns plot with number of iterations required for convergence.
  Combined <- Combined %>%
    RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)
  # Access and show first five Harmony embeddings for each barcode.
  Combined_harmony <- Embeddings(Combined, 'harmony')
  Combined_harmony[1:5, 1:5]
  
#----------------------------------------------------------------------------------------------------------------#
  
# UMAP post-Harmony, assignment of cell identity
  # Compute UMAP
  Combined <- FindNeighbors(Combined, dims = 1:18, reduction = "harmony") # use Harmony embeddings by setting reduction technique to "harmony"
  Combined <- FindClusters(Combined, resolution = 0.6)
  Combined <- RunUMAP(Combined, dims = 1:18, reduction = "harmony") # use Harmony embeddings by setting reduction technique to "harmony"
  # Plot UMAP and examine Ihh
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5)
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "orig.ident")
  FeaturePlot(Combined, features = "Ihh", order = T, pt.size = 0.5)
  # Differential expression script and cell type assignment. Change working directory to desired one beforehand
  Combined.markers <- FindAllMarkers(Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  Combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
  write.csv(Combined.markers, file = "Russellmarkers_unannotated.csv")
  Combined <- subset(Combined, idents = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "16", "17", "18"))
  
  #----------------------------------------------------------------------------------------------------------------#
  
  # Normalisation of RNA counts, scaling and PCA. Determine number of PCs to use for downstream analysis
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
  
  # Run Harmony and perform batch integration, returns plot with number of iterations required for convergence.
  Combined <- Combined %>%
    RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)
  # Access and show first five Harmony embeddings for each barcode.
  Combined_harmony <- Embeddings(Combined, 'harmony')
  Combined_harmony[1:5, 1:5]
  
  #----------------------------------------------------------------------------------------------------------------#
  
  # UMAP post-Harmony, assignment of cell identity
  # Compute UMAP
  Combined <- FindNeighbors(Combined, dims = 1:18, reduction = "harmony") # use Harmony embeddings by setting reduction technique to "harmony"
  Combined <- FindClusters(Combined, resolution = 0.6)
  Combined <- RunUMAP(Combined, dims = 1:18, reduction = "harmony") # use Harmony embeddings by setting reduction technique to "harmony"
  # Plot UMAP and examine Ihh
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5)
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "orig.ident")
  # Differential expression script and cell type assignment. Change working directory to desired one beforehand
  Combined2.markers <- FindAllMarkers(Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  Combined2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
  write.csv(Combined2.markers, file = "Russellmarkers_unannotated_2.csv")
  # Rename clusters based on differential expression
  Combined <- subset(Combined, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "20"))
  new.cluster.ids <- c("Interstitium",  # 0
                       "Podocytes",  # 1
                       "Myeloid",  # 2
                       "Endothelium",  # 3
                       "Early PT",  # 4
                       "Distal nephron",  # 5 
                       "Late PT",  # 6 
                       "Proliferating stroma",  # 7
                       "Medullary stroma",  #8
                       "Endothelium", #9
                       "Ureteric bud",  #10
                       "Cap mesenchyme",  # 11
                       "RV/SSB",  # 12
                       "Ren1+ stroma",  # 13
                       "Cortical stroma",  # 14
                       "Pelvic epithelium / CD",  # 15
                       "Ureteric stroma",  # 16
                       "Myeloid",  # 17
                       "Endothelium",  # 18
                       "Myeloid"  # 20
  )
  names(new.cluster.ids) <- levels(Combined)
  Combined <- RenameIdents(Combined, new.cluster.ids)
  DimPlot(Combined, reduction = "umap", label = F, pt.size = 0.5)
  DimPlot(Combined, reduction = "umap", label = F, pt.size = 0.5, split.by = "orig.ident")
  FeaturePlot(Combined, features = c("Ihh", "Shh", "Dhh"), order = T, pt.size = 0.5)
  DotPlot(Combined, features = c("Ptch1",
                                 "Ptch2",
                                 "Smo",
                                 "Sufu",
                                 "Gli1",
                                 "Gli2",
                                 "Gli3"), dot.scale = 10)
  
  table(Combined@active.ident, Combined@meta.data$orig.ident)

  table(Combined@meta.data$orig.ident)
  
  
  
  
  