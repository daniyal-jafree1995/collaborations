## Almalki & Arjun et al. - "Hyperglycaemic exacerbation of myocardial infarction through SGLT1"
## scRNA-seq / snRNA-seq analysis of cardiac expression of SGLT1 and SGLT2 transcripts
## Author: Daniyal Jafree (University College London)
## Version 1: 13/12/22

## Load packages and set working directory. Please change working directory as required.
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(tidyverse)
  library(Matrix)
  library(ggrepel)
  library(patchwork)
  library(tidyselect)
  library(harmony)

### DATASET 1: https://www.nature.com/articles/s44161-022-00028-6
  ## Load dataset 1 before subsetting 'donor' cells, creating UMAP of merged dataset and search genes
  load("~/Desktop/HeartAtlas/GSE183852_DCM_Integrated.Robj")
  DimPlot(RefMerge, pt.size = 0.5, raster= F, label = F)
  RefMerge <- SetIdent(RefMerge, value = RefMerge@meta.data$condition)
  sobj1 <- subset(RefMerge, subset = condition == "Donor")
  sobj1 <- SetIdent(sobj1, value = sobj1@meta.data$Names)
  DimPlot(sobj1, pt.size = 0.5, raster = F, label = F) # FIGURE 3A
  DimPlot(sobj1, pt.size = 0.5, raster = F, label = F, group.by = "orig.ident")
  FeaturePlot(sobj1, features = c("SLC5A1", "SLC5A2"), raster = F, pt.size = 0.5, order = T, keep.scale = "all") # FIGURE 3B
  table(sobj1@meta.data$condition)
  new.cluster.ids <- c("CM+", "CM-", "CM-", "CM-", "CM-", "CM-", "CM-", "CM-", "CM-", "CM-", "CM-", "CM-", "CM-", "CM-")
  names(new.cluster.ids) <- levels(sobj1)
  sobj1 <- RenameIdents(sobj1, new.cluster.ids)
  FindMarkers(sobj1, ident.1 = "CM+", ident.2 = NULL, only.pos = TRUE, logfc.threshold = 0, features = "SLC5A1")
  ## Subset cardiomyocytes from dataset 1 for further analysis
  sobj1 <- SetIdent(sobj1, value = sobj1@meta.data$Names)
  CM <- subset(sobj1, idents = c("Cardiomyocytes"))
  CM <- NormalizeData(CM, normalization.method = "LogNormalize", scale.factor = 10000)
  CM <- FindVariableFeatures(CM, selection.method = "vst", nfeatures = 2000)
  CM <- ScaleData(CM)
  CM <- RunPCA(CM, features = VariableFeatures(object = CM))
   # Empirically determine number of PCs to use for downstream analysis using elbow plot and script to calculate pcs
    ElbowPlot(CM, ndims = 50)
    pct <- CM[["pca"]]@stdev / sum(CM[["pca"]]@stdev) * 100
    cumu <- cumsum(pct)
    co1 <- which(cumu > 90 & pct < 5)[1]
    co1
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    co2
    pcs <- min(co1, co2)
    pcs
  ## Integration and unsupervised clustering of cardiomyocytes from dataset 1 using Harmony
  # KNN and Louvain clustering prior to UMAP. Normally resolution should be between 0.4-1.2 for single cell datasets
  CM <- FindNeighbors(CM, dims = 1:16)
  CM <- FindClusters(CM, resolution = 0.4)
  # Run and plot UMAP by cell types
  CM <- RunUMAP(CM, dims = 1:16)
  DimPlot(CM, reduction = "umap", label = TRUE, pt.size = 0.5) # can use group.by and call variable from object@meta.data to group cells in UMAP by variable
  # Run Harmony and perform batch correction, returns plot with number of iterations required for convergence.
  CM <- CM %>%
    RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)
  # Compute UMAP
  CM <- FindNeighbors(CM, dims = 1:16, reduction = "harmony")
  CM <- FindClusters(CM, resolution = 0.4)
  CM <- RunUMAP(CM, dims = 1:16, reduction = "harmony")
  # Generate UMAP and group by desired variables, and assess expression across subgroups from dataset 1
  DimPlot(CM, reduction = "umap", label = TRUE, pt.size = 0.5)
  DimPlot(CM, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "Sex")
  DimPlot(CM, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "Age_Group_Tertile")
  VlnPlot(CM, features = c("SLC5A1"), group.by = "Sex", pt.size = 0) # FIGURE 3D
  CM <- SetIdent(CM, value = CM@meta.data$Age_Group_Tertile)
  levels(CM) <- c('Young', 'Middle', 'Old')
  VlnPlot(CM, features = c("SLC5A1"), pt.size = 0) # FIGURE 3E
  DimPlot(CM, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "orig.ident")
  table(CM@meta.data$condition)
  CM <- SetIdent(CM, value = CM@meta.data$Sex)
  FindAllMarkers(CM, features = "SLC5A1", logfc.threshold = 0)
  CM <- SetIdent(CM, value = CM@meta.data$Age_Group_Tertile)
  FindAllMarkers(CM, features = "SLC5A1", logfc.threshold = 0)
  
### DATASET 2: https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.119.045401
  ## Create Seurat object from dataset 2
  heartdata2 <- Read10X("/Users/daniyaljafree/Desktop/HeartAtlas/BroadInst/Read10XInputs")
  metadata = read.delim("/Users/daniyaljafree/Desktop/HeartAtlas/BroadInst/meta.data.v3.txt", row.names = 1)
  metadata = metadata[-1,]
  sobj2 <- CreateSeuratObject(counts = heartdata2, project = "heart2", meta.data = metadata, min.cells = 3, min.features = 300)
  ## Subset cardiomyocytes from dataset 2 for further analysis
  sobj2 <- SetIdent(sobj2, value = sobj2@meta.data$Cluster)
  new.cluster.ids <- c("CM+", "CM-", "CM-", "CM-", "CM-", "CM+", "CM-", "CM-", "CM-", "CM-", "CM-", "CM+", "CM-", "CM-", "CM+", "CM+", "CM+")
  names(new.cluster.ids) <- levels(sobj2)
  sobj2 <- RenameIdents(sobj2, new.cluster.ids)
  FindMarkers(sobj2, ident.1 = "CM+", ident.2 = NULL, only.pos = TRUE, logfc.threshold = 0, features = "SLC5A1")
  table(sobj2@meta.data$orig.ident)
  
  CM2 <- subset(sobj2, idents = c("03. Atrial Cardiomyocyte", "04. Ventricular Cardiomyocyte I", "15. Ventricular Cardiomyocyte III", "06. Ventricular Cardiomyocyte II"))
  CM2 <- NormalizeData(CM2, normalization.method = "LogNormalize", scale.factor = 10000)
  CM2 <- FindVariableFeatures(CM2, selection.method = "vst", nfeatures = 2000)
  CM2 <- ScaleData(CM2)
  CM2 <- RunPCA(CM2, features = VariableFeatures(object = CM2))
  # Empirically determine number of PCs to use for downstream analysis using elbow plot and script to calculate pcs
    ElbowPlot(CM2, ndims = 50)
    pct <- CM2[["pca"]]@stdev / sum(CM2[["pca"]]@stdev) * 100
    cumu <- cumsum(pct)
    co1 <- which(cumu > 90 & pct < 5)[1]
    co1
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    co2
    pcs <- min(co1, co2)
    pcs
  ## Integration and unsupervised clustering of cardiomyocytes from dataset 1 using Harmony
  # KNN and Louvain clustering prior to UMAP. Normally resolution should be between 0.4-1.2 for single cell datasets
  CM2 <- FindNeighbors(CM2, dims = 1:12)
  CM2 <- FindClusters(CM2, resolution = 0.4)
  # Run and plot UMAP by cell types
  CM2 <- RunUMAP(CM2, dims = 1:12)
  DimPlot(CM2, reduction = "umap", label = TRUE, pt.size = 0.5) # can use group.by and call variable from object@meta.data to group cells in UMAP by variable
  # Run Harmony and perform batch correction, returns plot with number of iterations required for convergence.
  CM2 <- CM2 %>%
  RunHarmony("biological.individual", plot_convergence = TRUE, max.iter.harmony = 10)
  # Compute UMAP
  CM2 <- FindNeighbors(CM2, dims = 1:12, reduction = "harmony")
  CM2 <- FindClusters(CM2, resolution = 0.4)
  CM2 <- RunUMAP(CM2, dims = 1:12, reduction = "harmony")
  # Generate UMAP and group by desired variables, and assess expression across subgroups from dataset 1
  CM2 <- SetIdent(CM2, value = CM2@meta.data$orig.ident)
  DimPlot(CM2, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 0) # FIGURE 3C
  DimPlot(CM2, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "orig.ident", label.size = 0) # FIGURE 3C part 2
  DotPlot(CM2, features = c("SLC5A1"), group.by = "orig.ident", pt.size = 0)
  VlnPlot(CM2, features = c("SLC5A1"), group.by = "orig.ident", pt.size = 0) # FIGURE 3F
  table(CM2@meta.data$biological.individual)
  FindAllMarkers(CM2, features = "SLC5A1", logfc.threshold = 0)

#--------------------------------------------------### END ###--------------------------------------------------#

  # Redundant code
  #CM <- SetIdent(CM, value = CM@meta.data$orig.ident)
  #CM@meta.data$diabetic_status <- ifelse(test = CM@meta.data$orig.ident %in% c("TWCM-10-68",
                                                                               #"TWCM-11-41",
                                                                               #"TWCM-11-42",
                                                                               #"TWCM-11-78",
                                                                               #"TWCM-13-1",
                                                                               #"TWCM-13-101",
                                                                               #"TWCM-13-192",
                                                                               #"TWCM-13-168"), yes = "diabetes", no = "normoglycaemic")
  #CM@meta.data$valvularproblem<- ifelse(test = CM@meta.data$orig.ident %in% c("TWCM-11-41",
                                                                              #"TWCM-11-82",
                                                                              #"TWCM-13-198",
                                                                              #"TWCM-13-36",
                                                                              #"TWCM-11-78",
                                                                              #"TWCM-14-173",
                                                                              #"TWCM-13-96"), yes = "valveproblem", no = "normal")
  #CM <- SetIdent(CM, value = CM@meta.data$valvularproblem)
  #CMclean <- subset(CM, idents = "normal")
  #CMclean <- SetIdent(CMclean, value = CMclean@meta.data$diabetic_status)
  #VlnPlot(CMclean, features = c("SLC5A1"), group.by = "diabetic_status", pt.size = 0) # UNUSED, SUPPLEMENTARY
  #CM <- SetIdent(CM, value = CM@meta.data$diabetic_status)
  ##DimPlot(CMclean, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "diabetic_status")
  #FindAllMarkers(CM, features = "SLC5A1", logfc.threshold = 0)
  #DotPlot(CMclean, features = c("SLC5A1"), group.by = "diabetic_status") # UNUSED, SUPPLEMENTARY


