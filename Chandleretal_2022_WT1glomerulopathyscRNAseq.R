## scRNA-seq analysis for Wt1(+/R394W) mutant analysis for Dr Jennie Chandler
## Data derived from glomeruli isolated using dynabead technique from n = 2 littermate Ctrl and n = 2 Wt1(+/R394W) nephropathy, followed by 10x Genomics Chromium v3
## Authors: Gideon Pomeranz (UCL), Daniyal Jafree (UCL), Dr Andrew Mason (University of York)
## Version 4: 23/08/2022 - final version pre-submission
## Last updated: 03/10/2022
#----------------------------------------------------------------------------------------------------------------#
## Load packages and set working directory. Please change working directory as required.
  # Load packages
  set.seed(42)
  library(dplyr)
  library(Seurat)
  library(SeuratDisk)
  library(patchwork)
  library(tidyverse)
  library(Matrix)
  library(ggrepel)
  library(patchwork)
  library(tidyselect)
  library(SoupX)
  library(harmony)
  library(DoubletFinder)
  library(EnhancedVolcano)
  library(VennDiagram)
  library(NICHES)
  
  # Set working directory to save outputs to. Please change as required.
  setwd("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie")
#----------------------------------------------------------------------------------------------------------------#
## PART 1: Dr Andrew Mason's script for getting metadata from Seurat object. Assumed you have navigated to the directory
  # read in CellRanger data
  counts <- Read10X("Data")
  dim(counts) # rows = genes, columns = cells
  # define labels
  ids <- c("S1","S2","S3","S4")
  con <- c("wildtype","wildtype","hetR394W","hetR394W")
  # set up sample and conditions labels for the backgrounds
  cell_counts <- table(substring(colnames(counts),18))
  samples <- c()
  for (s in 1:length(ids)) { 
    samples <- c(samples, (rep(paste(ids)[s], cell_counts[s])))
  } 
  conditions <- c()
    for (c in 1:length(con)) { 
  conditions <- c(conditions, (rep(paste(con)[c], cell_counts[c])))
  } 
  # create metadata df
  metadata <- data.frame(Sample=samples, Conditions=conditions, row.names=colnames(counts)) # create a metadata df
  # check cell numbers match up
  metadata %>% group_by(Sample) %>% summarise(Cells=n())
  metadata %>% group_by(Conditions) %>% summarise(Cells=n())
  # Load data into a Seurat object
  Combined <- CreateSeuratObject(counts=counts, project="JC_glom", meta.data=metadata)
  Idents(Combined) <- Combined@meta.data$Conditions
  # Part 1 redundant object clean up
  rm(c, s, cell_counts, con, conditions, ids, samples, counts)
#----------------------------------------------------------------------------------------------------------------#
  ## Gideon depreciated data loading, T2G and SeuratDisk script
    # load dataset and transcript to gene
      #setwd("~//Documents/PhD_2018/R_analysis/wt1_mutant")
      #tr2g <- read_tsv("t2g.txt", col_names = c("transcript", "gene", "gene_symbol")) %>%
      #dplyr::select(-transcript) %>%
      #distinct()
    # use SeuratDisk to change h5ad into Seurat objects
      #Convert("control1.h5ad", dest = "h5seurat", overwrite = TRUE)
      #control1 <- LoadH5Seurat("control1.h5seurat")
    # extract the counts and change to gene names
      #c1_counts <- control1@assays$RNA@counts
      #rownames(c1_counts) <- tr2g$gene_symbol[match(rownames(c1_counts), tr2g$gene)]
      #control1 <- CreateSeuratObject(c1_counts, project = "control1")
    # Convert("control2.h5ad", dest = "h5seurat", overwrite = TRUE)
      #control2 <- LoadH5Seurat("control2.h5seurat")
      #c2_counts <- control2@assays$RNA@counts
      #rownames(c2_counts) <- tr2g$gene_symbol[match(rownames(c2_counts), tr2g$gene)]
      #control2 <- CreateSeuratObject(c2_counts, project = "control2")
    # Convert("mutant1.h5ad", dest = "h5seurat", overwrite = TRUE)
      #mutant1 <- LoadH5Seurat("mutant1.h5seurat")
      #m1_counts <- mutant1@assays$RNA@counts
      #rownames(m1_counts) <- tr2g$gene_symbol[match(rownames(m1_counts), tr2g$gene)]
      #mutant1 <- CreateSeuratObject(m1_counts, project = "mutant1")
    # Convert("mutant2.h5ad", dest = "h5seurat", overwrite = TRUE)
      #mutant2 <- LoadH5Seurat("mutant2.h5seurat")
      #m2_counts <- mutant2@assays$RNA@counts
      #rownames(m2_counts) <- tr2g$gene_symbol[match(rownames(m2_counts), tr2g$gene)]
      #mutant2 <- CreateSeuratObject(m2_counts, project = "mutant2")
      #wt1_data <- merge(control1, y=c(control2,mutant1,mutant2), add.cell.ids = c("wt1", "WT2", "Mut1", "Mut2"), project = "wt1")
      #wt1 <- wt1_data
      #rm(c1_counts,c2_counts, m1_counts, m2_counts, control1, control2,mutant1, mutant2)
      #rm(wt1_data)
#----------------------------------------------------------------------------------------------------------------#
## PART 2: Quality control
  # Calculate percentage of reads mapping to mitochondrial genome
  Combined[["percent.mt"]] <- PercentageFeatureSet(Combined, pattern="^mt-")
  # subset cells by mitochondrial contact and number of features
  Combined <- subset(Combined, subset = nFeature_RNA > 150 & percent.mt < 20) # can optionally assign upper cutoff for nFeature_RNA ?7500
  Combined
  print(VlnPlot(Combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident",ncol = 3))
  plot1 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")
  plot2 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")
  print(plot1 + plot2)
#----------------------------------------------------------------------------------------------------------------#
  ## Gideon depreciated script for cell cycle correction
    # Calculate cell cycle scores These are human gene names, might want to change to mouse for this to work.
      #s.genes <- cc.genes.updated.2019$s.genes
      #g2m.genes <- cc.genes.updated.2019$g2m.genes
      #wt1 <- CellCycleScoring(wt1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    # regressing out all cell cycle effects can blur the distinction between stem and progenitor cells as well.
    # regressing out the difference between the G2M and S phase scores. 
    # This means that signals separating non-cycling cells and cycling cells will be maintained, 
    # but differences in cell cycle phase amongst proliferating cells (which are often uninteresting), will be regressed out of the data
      #wt1$CC.Difference <- wt1$S.Score - wt1$G2M.Score
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
  # save PCA as png
  #pdf(paste("Figures", "PCA.pdf", sep = "/"))
  #print(DimPlot(wt1, reduction = "pca", group.by = "orig.ident"))
  #dev.off()
  #pdf(paste("Figures", "Elbow_plot.pdf", sep = "/"))
  #print(ElbowPlot(wt1,ndims = 50))
  #dev.off()
#----------------------------------------------------------------------------------------------------------------#
  ## Depreciated NCA script from Gideon 
      # X <- as.matrix(t(wt1@assays[["RNA"]]@scale.data))
      # y <- wt1@meta.data[["RNA_snn_res.0.9"]]
      # nca = sk$NeighborhoodComponentsAnalysis(n_components=as.integer(10),random_state=as.integer(42), verbose=as.integer(1))
      # nca_fit_0.6 = nca$fit_transform(X, y)
      # write.csv(nca_fit_0.6, file ="nca_fit_0.6_rna.csv")
    # after saving nca results you can read them back in. 
      # nca_fit_0.6 <- read.csv("nca_fit_0.6.csv")
      # nca_fit_0.6$X <- NULL
      # nca_fit_0.6 <- as.matrix(nca_fit_0.6)
      # colnames(nca_fit_0.6) <- paste0("NCA_", 1:10)
    # transpose results to fit seurat
      # wt1[["nca_0.9"]] <- CreateDimReducObject(embeddings = nca_fit_0.9, key = "nca09_", assay = DefaultAssay(wt1))
      # DimPlot(wt1, reduction = "nca_0.9", pt.size = 0.5,cells = 1:19373)
#----------------------------------------------------------------------------------------------------------------#
## PART 4: Unsupervised clustering using UMAP. Some commented-out options below that Gideon prefers to use.
  # KNN and Louvain clustering prior to UMAP. Normally resolution should be between 0.4-1.2 for single cell datasets
  Combined <- FindNeighbors(Combined, dims = 1:16)
  Combined <- FindClusters(Combined, resolution = 0.4) # alternative Leiden approach from Gideon: FindClusters(wt1, resolution = 0.9, method = "igraph", algorithm = 4)
  # Run and plot UMAP by cell types
  Combined <- RunUMAP(Combined, dims = 1:16)
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5) # can use group.by and call variable from object@meta.data to group cells in UMAP by variable
  #save UMAPS as png
  # pdf(paste("Figures", "UMAP_samples.pdf", sep = "/"))
  # print(DimPlot(wt1, reduction = "umap", group.by = "orig.ident"))
  # dev.off()
  # pdf(paste("Figures", "UMAP_noID.pdf", sep = "/"))
  # print(DimPlot(wt1, reduction = "umap", group.by = "seurat_clusters"))
  # dev.off()
  # Gideon alternatives using python UMAP
    #wt1 <- RunUMAP(wt1, dims=1:18, reduction = "harmony", umap.method = "umap-learn", metric = "correlation")
    #wt1 <- RunUMAP(wt1, dims=1:14, umap.method = "umap-learn", metric = "correlation")
    #wt1 <- RunTSNE(wt1, dims=1:14, reduction = "harmony")
#----------------------------------------------------------------------------------------------------------------#
## PART 5: Implement SoupX package for estimation and correction of raw count matrix for contaminating RNA. Then run Parts 1-4 again.
  # Create raw count matrix from Seurat object
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
#----------------------------------------------------------------------------------------------------------------#
  ## Daniyal depreciated Harmony script. In preliminary analysis, only mutant podocytes clustered seperately, but perhaps this is biologically meaningful to show
    # Run Harmony and perform batch correction, returns plot with number of iterations required for convergence.
      #Combined <- Combined %>%
      #RunHarmony("Conditions", plot_convergence = TRUE)
    # Access and show first five Harmony embeddings for each barcode.
      #Combined_harmony <- Embeddings(Combined, 'harmony')
      #Combined_harmony[1:5, 1:5]
#----------------------------------------------------------------------------------------------------------------#
## PART 6: UMAP post-SoupX (+/- Harmony), further cleanup and assignment of cell identity
  # Compute UMAP
  Combined <- FindNeighbors(Combined, dims = 1:16) # use Harmony embeddings by setting reduction technique to "harmony"
  Combined <- FindClusters(Combined, resolution = 0.5)
  Combined <- RunUMAP(Combined, dims = 1:16) # use Harmony embeddings by setting reduction technique to "harmony"
  # Reorder identities so that control comes before disease in all successive graphs.
  my_levels <- c("wildtype", "hetR394W")
  Combined@meta.data$Conditions <- factor(x = Combined@meta.data$Conditions, levels = my_levels)
  # Generate UMAP and group by desired variables
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5)
  DimPlot(Combined, reduction = "umap", group.by = "Conditions", pt.size = 0.5)
  # Differential expression script by Gideon before cell type assignment
  #Combined.markers <- FindAllMarkers(Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  #Combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
  #write.csv(Combined.markers, file = "Chanderetalmarkers_unannotated.csv")
  # Cleanup to remove cluster 8 (likely dying endothelial cells), cluster 15 (endothelial doublets), cluster 17 (non-specific proliferating cluster)
  Combined <- subset(Combined, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "9", "10", "11", "12", "13", "14", "16"))
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5)
  DimPlot(Combined, reduction = "umap", split.by = "Conditions", pt.size = 0.5)
  # Rename clusters based on differential expression in preliminary analysis
  new.cluster.ids <- c("Podo_R394W",  # 0
                       "Podo_wt",  # 1
                       "Podo_wt",  # 2
                       "Podo_wt",  # 3
                       "GEC",  # 4
                       "Podo_R394W",  # 5 
                       "Podo_R394W",  # 6 
                       "Mesangium",  # 7
                       "EA_Endo",  #9
                       "AA_Endo",  #10
                       "SMC",  # 11
                       "PEC",  # 12
                       "Myeloid",  # 13
                       "B_lymph",  # 14
                       "T_lymph"  # 16
  )
  names(new.cluster.ids) <- levels(Combined)
  Combined <- RenameIdents(Combined, new.cluster.ids)
  DimPlot(Combined, reduction = "umap", label = F, pt.size = 1)
#----------------------------------------------------------------------------------------------------------------#
## PART 7: Doublet identification and removal
  # DoubletFinder script from Daniyal
  # In preliminary analysis this showed bias towards removing 'real' mesangial cells, but also identified a doublet endothelial cluster.
  # Creates 2294 artificial doublets, merges with dataset, finds artificial k nearest neighbours and estimates doublets in real data.
  # pK value of 0.13 based on preliminary analysis using pararmSweep_v3, summarizeSweep and find.pK functions.
  sweep.res.list_kidney <- paramSweep_v3(Combined, PCs = 1:16, sct = FALSE)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  annotations <- Combined@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           
  nExp_poi <- round(0.05*nrow(Combined@meta.data)) # alternatively use nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  Combined <- doubletFinder_v3(Combined, PCs = 1:16, pN = 0.25, pK = 0.13, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "DF.classifications_0.25_0.13_344")  # amend DF.classifications as required
  # Makes a seperate Seurat object of doublets, 
  doublets <- subset(Combined, subset = DF.classifications_0.25_0.13_344 == "Doublet") # amend DF.classifications as required
  doublets.cells <- WhichCells(doublets, idents = c("Mesangium", "SMC"), invert = TRUE)
  Combined <- Combined[,!colnames(Combined) %in% doublets.cells]
#----------------------------------------------------------------------------------------------------------------#
## PART 8: UMAP by different parameters, cell type verification using differential expression and counting
  # UMAPs by cell type and by condition
  DimPlot(Combined, reduction = "umap", label = F, pt.size = 1, cols = c("Podo_R394W" = "dodgerblue4",
                                                                         "Podo_wt" = "dodgerblue",
                                                                         "GEC"  = "red",
                                                                         "Mesangium" = "seagreen",
                                                                         "EA_Endo" = "darkorchid",
                                                                         "AA_Endo" = "tomato1",
                                                                         "SMC" = "brown",
                                                                         "PEC" = "yellow2",
                                                                         "Myeloid" = "gray83",
                                                                         "B_lymph" = "gray46",
                                                                         "T_lymph" = "gray20"
  )) ### GENERATES FIGURE 1E
  DimPlot(Combined, reduction = "umap", split.by = "Conditions", pt.size = 1, cols = c("Podo_R394W" = "dodgerblue4",
                                                                                       "Podo_wt" = "dodgerblue",
                                                                                       "GEC"  = "red",
                                                                                       "Mesangium" = "seagreen",
                                                                                       "EA_Endo" = "darkorchid",
                                                                                       "AA_Endo" = "tomato1",
                                                                                       "SMC" = "brown",
                                                                                       "PEC" = "yellow2",
                                                                                       "Myeloid" = "gray83",
                                                                                       "B_lymph" = "gray46",
                                                                                       "T_lymph" = "gray20"
  )) ### GENERATES FIGURE S1A
  # Differential expression script by Gideon after cell type assignment
  #Combined.markers <- FindAllMarkers(Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  #Combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
  #write.csv(Combined.markers, file = "Chanderetalmarkers_annotated.csv")
  # Dotplot of canonical markers for cell types identified in the dataset
  glomerular.cell.markers <- c("Wt1", "Nphs2", # podocyte markers
                               "Emcn", "Ehd3", # GEC markers
                               "Ptn", "Pdgfrb", # mesangial markers
                               "Sox17", # Arterial markers including...
                               "Plvap", # EA markers
                               "Edn1", # AA markers
                               "Acta2", "Myh11", # SMC markers
                               "Pax8", "Cldn1", # PEC markers
                               "Ptprc", "Lyz2", # Myeloid markers
                               "Cd79a", "Igkc", # B lymphocyte markers
                               "Cd3e", "Trbc2") # T lymphocyte markers
  DotPlot(Combined, features = glomerular.cell.markers, dot.scale = 5)
  #saveRDS(Combined, file = "/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/WT1_scRNAseq.rds")
  Combined <- readRDS("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/WT1_scRNAseq.rds")
  # Assess numbers of each cell type by condition
  cell.numbers <- table(Combined@active.ident, Combined@meta.data$Conditions)
  cell.numbers <- table(Combined@active.ident)
  cell.numbers.by.sample <- table(Combined@active.ident, Combined@meta.data$Sample)
  write.csv(cell.numbers.by.sample, file = "Chanderetalcellnumbersbysample.csv")
#----------------------------------------------------------------------------------------------------------------#
## PART 9: Differential expression, comparing cell type-specific transcriptome between conditions using volcano plots and heatmap
  # Compare podocytes from Ctrl and mutant kidneys
  Podo_subset <- subset(Combined, idents = c("Podo_wt", "Podo_R394W"))
  my_levels <- c("Podo_wt", "Podo_R394W")
  Podo_subset@active.ident <- factor(x = Podo_subset@active.ident, levels = my_levels)
  Podocyte.DE <- FindMarkers(Podo_subset, ident.1 = "Podo_R394W", ident.2 = "Podo_wt")
  Podocyte.DE <- Podocyte.DE %>% 
    as.data.frame() %>% 
    arrange(desc(avg_log2FC))
  Podocyte.DE <- slice(Podocyte.DE, -c(1, 2, 3, 4),)
  write.csv(Podocyte.DE, "compare_podocytes.csv")
  # Draw volcano plot for podocytes using EnhancedVolcano
  Podocyte.DE.lowthreshold <- FindMarkers(Podo_subset, ident.1 = "Podo_R394W", ident.2 = "Podo_wt", min.pct = 0.1, logfc.threshold = 0.05)
  Podocyte.DE.lowthreshold <- Podocyte.DE.lowthreshold %>% 
    as.data.frame() %>% 
    arrange(desc(avg_log2FC))
  Podocyte.DE.lowthreshold <- slice(Podocyte.DE.lowthreshold, -c(1, 2, 3, 4),)
  EnhancedVolcano(Podocyte.DE.lowthreshold,
                  lab = rownames(Podocyte.DE.lowthreshold),
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  FCcutoff = 0.25,
                  pCutoff = 0.05,
                  colCustom = NULL,
                  boxedLabels = F,
                  drawConnectors = F,
                  labSize = 2,
                  col = c('black', 'black', 'blue', 'red'),
                  widthConnectors = 0.5)
  # Compare mesangial cells from Ctrl and mutant kidneys
  Mesangium_subset <- subset(Combined, idents = "Mesangium")
  my_levels <- c("wildtype", "hetR394W")
  Mesangium_subset@meta.data$Conditions <- factor(x = Mesangium_subset@meta.data$Conditions, levels = my_levels)
  Idents(Mesangium_subset) <- "Conditions"
  Mesangium.DE <- FindAllMarkers(Mesangium_subset)
  Mesangium.DE <- Mesangium.DE %>% 
    as.data.frame() %>% 
    arrange(desc(avg_log2FC))
  to_remove <- Mesangium.DE[endsWith(rownames(Mesangium.DE), ".1"),]
  Mesangium.DE <- Mesangium.DE[!(rownames(Mesangium.DE) %in% rownames(to_remove)),]
  write.csv(Mesangium.DE, "Mesangium.DE.csv")
  # Draw volcano plot for mesangial cells using EnhancedVolcano
  Mesangium.DE.lowthreshold <- FindAllMarkers(Mesangium_subset, min.pct = 0.1, logfc.threshold = 0.05)
  Mesangium.DE.lowthreshold <- Mesangium.DE.lowthreshold %>% 
    as.data.frame() %>% 
    arrange(desc(avg_log2FC))
  to_remove <- Mesangium.DE.lowthreshold[endsWith(rownames(Mesangium.DE.lowthreshold), ".1"),]
  Mesangium.DE.lowthreshold <- Mesangium.DE.lowthreshold[!(rownames(Mesangium.DE.lowthreshold) %in% rownames(to_remove)),]
  EnhancedVolcano(Mesangium.DE.lowthreshold,
                  lab = rownames(Mesangium.DE.lowthreshold),
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  FCcutoff = 0.25,
                  pCutoff = 0.05,
                  colCustom = NULL,
                  boxedLabels = F,
                  drawConnectors = F,
                  labSize = 0,
                  col = c('black', 'black', 'blue', 'red'),
                  widthConnectors = 0.5)
  # Compare GECs from Ctrl and mutant kidneys
  GEC_subset <- subset(Combined, idents = "GEC")
  GEC_subset@meta.data$Conditions <- factor(x = GEC_subset@meta.data$Conditions, levels = my_levels)
  Idents(GEC_subset) <- "Conditions"
  GEC.DE <- FindAllMarkers(GEC_subset)
  GEC.DE <- GEC.DE %>% 
    as.data.frame() %>% 
    arrange(desc(avg_log2FC))
  to_remove <- GEC.DE[endsWith(rownames(GEC.DE), ".1"),]
  GEC.DE <- GEC.DE[!(rownames(GEC.DE) %in% rownames(to_remove)),]
  write.csv(GEC.DE, "compare_GEC.csv")
  # Draw volcano plot for GECs using EnhancedVolcano
  GEC.DE.lowthreshold <- FindAllMarkers(GEC_subset, min.pct = 0.1, logfc.threshold = 0.05)
  GEC.DE.lowthreshold <- GEC.DE.lowthreshold %>% 
    as.data.frame() %>% 
    arrange(desc(avg_log2FC))
  to_remove <- GEC.DE.lowthreshold[endsWith(rownames(GEC.DE.lowthreshold), ".1"),]
  GEC.DE.lowthreshold <- GEC.DE.lowthreshold[!(rownames(GEC.DE.lowthreshold) %in% rownames(to_remove)),]
  EnhancedVolcano(GEC.DE.lowthreshold,
                  lab = rownames(GEC.DE.lowthreshold),
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  FCcutoff = 0.25,
                  pCutoff = 0.05,
                  colCustom = NULL,
                  boxedLabels = F,
                  drawConnectors = F,
                  labSize = 0,
                  col = c('black', 'black', 'blue', 'red'),
                  widthConnectors = 0.5)
  # Compare PECs from Ctrl and mutant kidneys
  PEC_subset <- subset(Combined, idents = "PEC")
  PEC_subset@meta.data$Conditions <- factor(x = PEC_subset@meta.data$Conditions, levels = my_levels)
  Idents(PEC_subset) <- "Conditions"
  PEC.DE <- FindAllMarkers(PEC_subset)
  PEC.DE <- PEC.DE %>% 
    as.data.frame() %>% 
    arrange(desc(avg_log2FC))
  to_remove <- PEC.DE[endsWith(rownames(GEC.DE), ".1"),]
  PEC.DE <- PEC.DE[!(rownames(PEC.DE) %in% rownames(to_remove)),]
  write.csv(PEC.DE, "compare_PEC.csv")
  # Draw volcano plot for GECs using EnhancedVolcano
  PEC.DE.lowthreshold <- FindAllMarkers(PEC_subset, min.pct = 0.1, logfc.threshold = 0.05)
  PEC.DE.lowthreshold <- PEC.DE.lowthreshold %>% 
    as.data.frame() %>% 
    arrange(desc(avg_log2FC))
  to_remove <- PEC.DE.lowthreshold[endsWith(rownames(PEC.DE.lowthreshold), ".1"),]
  PEC.DE.lowthreshold <- PEC.DE.lowthreshold[!(rownames(PEC.DE.lowthreshold) %in% rownames(to_remove)),]
  EnhancedVolcano(PEC.DE.lowthreshold,
                  lab = rownames(PEC.DE.lowthreshold),
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  FCcutoff = 0.25,
                  pCutoff = 0.05,
                  colCustom = NULL,
                  boxedLabels = F,
                  drawConnectors = F,
                  labSize = 0,
                  col = c('black', 'black', 'blue', 'red'),
                  widthConnectors = 0.5)
  # Create heatmaps of genes uniquely expressed by podocytes in mutant as compared to other glomerular pathologies, vs conserved across glomerular pathologies (shared across a minimum of three pathologies.)
  Unique_genes <- read.csv("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/Data/JASN/WT1_comparison/Compare/Combined_unique_thresholded.csv")
  Unique_gene_list <- Unique_genes$gene
  DoHeatmap(Podo_subset, features = Unique_gene_list)
  Conserved_genes <- read.csv("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/WT1_model_Jennie/Data/JASN/WT1_comparison/Compare/Combined_conserved_thresholded.csv")
  Conserved_genes_list <- Conserved_genes$gene
  DoHeatmap(Podo_subset, features = c("R3hdml", "Adm", "Dusp14", "Gpx3", "Vamp5", "F2r", "Nr4a1", "Tspan2", "Ankrd1", "Lgals1", "Gm1673", "Gadd45b", "Clu", "Fuca1", "Cxcl1", "Kazn", "Crip1", "Tnfrsf12a", "Csrp1", "Cnn2", "Bok", "Net1", "Prss23", "Fam129b", "Pard3b", "Snx31", "Mafb", "Sulf1"))
  # Violin plots showing top 5 upregulated and top 5 downregulated DEGs (encoding proteins) in mutant podocytes 
  VlnPlot(Podo_subset, features = c("R3hdml", "Cdkn1c", "Gpc3", "Aplp1", "Pkm"), pt.size = 0)
  VlnPlot(Podo_subset, features = c("Astn2", "Col1a2", "H2-D1", "B2m", "Mylk"), pt.size = 0)
  VlnPlot(Combined, features = c("Cxcl12", "Dpp4"), split.by = "Conditions")
  # Expression values for genes of interest requested by Jennie, done by cell type
  podogenesofinterest <- FindAllMarkers(Podo_subset, features = c("Npnt", "Anxa1", "Anxa2", "Col4a3", "Col4a4", "Col4a5", "Rtn4", "Fgfr2", "Adam17", "Timp3"), logfc.threshold = 0, min.pct = 0, return.thresh = 1)
  GECgenesofinterest <- FindAllMarkers(GEC_subset, features = c("Dysf", "Flt1", "Kdr"), logfc.threshold = 0, min.pct = 0, return.thresh = 1)
  mesangialgenesofinterest <- FindAllMarkers(Mesangium_subset, features = c("Itga8", "Itga1"), logfc.threshold = 0, min.pct = 0, return.thresh = 1)
  PECgenesofinterest <- FindAllMarkers(PEC_subset, features = c("Lingo1", "Sdc4", "Erbb4"), logfc.threshold = 0, min.pct = 0, return.thresh = 1)
  # Change directory at this point if required
  write.csv(podogenesofinterest, "podogenesofinterest.csv")
  write.csv(GECgenesofinterest, "GECgenesofinterest.csv")
  write.csv(mesangialgenesofinterest, "mesangialgenesofinterest.csv")
  write.csv(PECgenesofinterest, "PECgenesofinterest.csv")
  VlnPlot(Podo_subset, features = c("Npnt", "Anxa1", "Anxa2", "Col4a3", "Col4a4", "Col4a5", "Rtn4", "Fgfr2", "Adam17", "Timp3"), pt.size = 0)
  VlnPlot(GEC_subset, features = c("Dysf", "Flt1", "Kdr"), pt.size = 0)
  VlnPlot(Mesangium_subset, features = c("Itga8", "Itga1"), pt.size = 0)
  VlnPlot(PEC_subset, features = c("Lingo1", "Sdc4", "Erbb4"), pt.size = 0)
  
#----------------------------------------------------------------------------------------------------------------#
                          
## PART 10: Computation of putative cell-cell interactions
  # Utilises recently released NICHES package, preprint available at:
  Combined_tuft <- subset(Combined, idents = c("Podo_wt", "Podo_R394W", "GEC", "Mesangium", "PEC"))
  Combined_tuft <- RenameIdents(Combined_tuft, "Podo_wt"  = "Podocyte", "Podo_R394W"  = "Podocyte")
  table(Combined_tuft@active.ident, Combined_tuft@meta.data$Conditions)
  Combined_tuft@meta.data$celltype_aggregate = paste(Combined_tuft@active.ident, Combined_tuft@meta.data$Conditions, sep = "_")
  Combined_tuft <- SetIdent(Combined_tuft, value = "celltype_aggregate")
  DimPlot(Combined_tuft, reduction = "umap", label = F, pt.size = 1)
  Tuft_network_niches <- RunNICHES((Combined_tuft),
                                  assay = 'RNA',
                                  species = 'mouse',
                                  LR.database = 'omnipath',
                                  CellToCell = T)
  Tuft_network <- Tuft_network_niches$CellToCell
  # Subsets analysis to only include putative communications with podocytes
  # NICHES analysis split by sender / recipient   
    # Analysis of podocyte-derived signals
    Signal_from_podo_NICHES <- subset(Tuft_network, idents = c("Podocyte_wildtype—GEC_wildtype", "Podocyte_wildtype—Mesangium_wildtype", "Podocyte_wildtype—PEC_wildtype", "Podocyte_hetR394W—GEC_hetR394W", "Podocyte_hetR394W—Mesangium_hetR394W", "Podocyte_hetR394W—PEC_hetR394W"))
    my_levels.cells <- c("Podocyte_wildtype—GEC_wildtype", "Podocyte_hetR394W—GEC_hetR394W", "Podocyte_wildtype—Mesangium_wildtype", "Podocyte_hetR394W—Mesangium_hetR394W", "Podocyte_wildtype—PEC_wildtype", "Podocyte_hetR394W—PEC_hetR394W")
    Signal_from_podo_NICHES@active.ident <- factor(x = Signal_from_podo_NICHES@active.ident, levels = my_levels.cells)
    Signal_from_podo_NICHES <- ScaleData(Signal_from_podo_NICHES)
    Signal_from_podo_NICHES <- FindVariableFeatures(Signal_from_podo_NICHES)
    Signal_from_podo_NICHES <- RunPCA(Signal_from_podo_NICHES)
    ElbowPlot(Signal_from_podo_NICHES,ndims=50)
    PCHeatmap(Signal_from_podo_NICHES,dims = 1:5,balanced = T,cells = 100) # Change dims according to elbow plot
    Signal_from_podo_NICHES <- RunUMAP(Signal_from_podo_NICHES,dims = 1:20)
    # UMAP of predicated cell-cell interactions, first by sending cell type then by recieving cell type
    DimPlot(Signal_from_podo_NICHES,reduction = 'umap',label = F, pt.size = 2)
    DimPlot(Signal_from_podo_NICHES,reduction = 'umap',group.by = 'SendingType',label = T)
    DimPlot(Signal_from_podo_NICHES,reduction = 'umap',group.by = 'ReceivingType',label = F, pt.size = 1)
    #Idents(Signal_from_podo_NICHES) <- genotype #group by genotype
    # Heatmap of pairwise cell-cell interactions with wildtype or mutant podocytes
    podo_signals <- FindAllMarkers(Signal_from_podo_NICHES, min.pct = 0.05, logfc.threshold = 0.05, only.pos = T)
    write.csv(podo_signals, "podocyte_sender.csv") # Save and export for curation
    #podo_signals_curated <- podo_signals %>% filter(row_number() %in% c(1,2,3,4,5,7,10,13,19,20,24,25,26,27,32,36,38,44,46,47,181,182,185,186,190,192,194,196,204,207,211,213,216,223,224,238,243,245,257,258,387,388,389,397,399,400,402,403,407,409,411,412,416,417,419,422,431,432,436,439,679,683,684,689,690,692,695,698,699,708,710,713,724,725,726,727,729,732,735,736,998,999,1000,1001,1002,1004,1007,1013,1014,1016,1017,1018,1023,1024,1029,1033,1035,1036,1037,1044,1240,1244,1247,1254,1255,1256,1258,1260,1263,1264,1265,1267,1269,1271,1274,1275,1276,1277,1280,1286))
    #GOI_niche <- podo_signals_curated %>% group_by(cluster) %>% top_n(13, avg_log2FC)
    Podo_GEC_signals <- podo_signals %>% filter(row_number() %in% c(1,2,5,10,19,26,27,32,36,38,44,46,47,182,192,4,196,204,224,243,181,185,186,190,207,211,213,216,223,238,245,257,258))
    Podo_Mes_signals <- podo_signals %>% filter(row_number() %in% c(387,388,389,397,399,400,402,403,407,409,411,412,416,417,419,422,431,432,436,439,679,683,684,689,690,692,695,698,699,708,710,713,724,725,726,727,729,732,735,736))
    Podo_PEC_signals <- podo_signals %>% filter(row_number() %in% c(998,999,1000,1001,1002,1004,1007,1013,1014,1016,1017,1018,1023,1024,1029,1033,1035,1036,1037,1044,1240,1244,1247,1254,1255,1256,1258,1260,1263,1264,1265,1267,1269,1271,1274,1275,1276,1277,1280,1286))
    #write.csv(Podo_GEC_signals, "Podo_GEC_signals.csv") # Save and export for sanity check
    #write.csv(Podo_Mes_signals, "Podo_Mes_signals.csv") # Save and export for sanity check
    #write.csv(Podo_PEC_signals, "Podo_PEC_signals.csv") # Save and export for sanity check
    # Podocyte-GEC heatmap
    Podo_GEC_signals <- Podo_GEC_signals[c(1,2,4,5,6,7,8,9,10,11,12,13,14,16,20,3,21,22,28,30,15,17,18,19,23,24,25,26,27,29,31,32,33),]
    DoHeatmap(Signal_from_podo_NICHES, features = Podo_GEC_signals$gene)+ 
      scale_fill_gradientn(colors = c("purple","black", "yellow"))
    # Podocyte-Mesangium heatmap
    Podo_Mes_signals <- Podo_Mes_signals[c(1,2,8,11,14,17,18,20,23,24,25,30,32,6,7,36,12,13,15,19,22,26,31,34,35,38,39,40),]
    DoHeatmap(Signal_from_podo_NICHES, features = Podo_Mes_signals$gene)+ 
      scale_fill_gradientn(colors = c("purple","black", "yellow"))
    # Podocyte-PEC heatmap
    Podo_PEC_signals <- Podo_PEC_signals[c(6,7,8,12,16,18,19,20,1,2,3,4,5,9,10,11,13,14,15,17,23,27,34,36,37,38,40),]
    DoHeatmap(Signal_from_podo_NICHES, features = Podo_PEC_signals$gene)+ 
      scale_fill_gradientn(colors = c("purple","black", "yellow"))
    # Depreciated podocyte-GEC crosstalk subcluster analysis
    #Signal_from_podo_GEC <- subset(Tuft_network, idents = c("Podocyte_wildtype—GEC_wildtype", "Podocyte_hetR394W—GEC_hetR394W"))
    #my_levels.cells <- c("Podocyte_wildtype—GEC_wildtype", "Podocyte_hetR394W—GEC_hetR394W")
    #Signal_from_podo_GEC@active.ident <- factor(x = Signal_from_podo_GEC@active.ident, levels = my_levels.cells)
    #Signal_from_podo_GEC <- ScaleData(Signal_from_podo_GEC)
    #Signal_from_podo_GEC <- FindVariableFeatures(Signal_from_podo_GEC)
    #Signal_from_podo_GEC <- RunPCA(Signal_from_podo_GEC)
    #ElbowPlot(Signal_from_podo_GEC,ndims=50)
    #PCHeatmap(Signal_from_podo_GEC,dims = 1:5,balanced = T,cells = 100) # Change dims according to elbow plot
    #Signal_from_podo_GEC <- FindNeighbors(Signal_from_podo_GEC, dims = 1:20) # use Harmony embeddings by setting reduction technique to "harmony"
    #Signal_from_podo_GEC <- FindClusters(Signal_from_podo_GEC, resolution = 0.3)
    #Signal_from_podo_GEC <- RunUMAP(Signal_from_podo_GEC,dims = 1:20)
    # UMAP of predicated cell-cell interactions, first by sending cell type then by recieving cell type
    #DimPlot(Signal_from_podo_GEC,reduction = 'umap',label = F, pt.size = 2)
    #DimPlot(Signal_from_podo_GEC,reduction = 'umap',split.by = 'SendingType',label = T)
    #DimPlot(Signal_from_podo_GEC,reduction = 'umap',split.by = 'ReceivingType',label = F, pt.size = 1)
    #Idents(Signal_from_podo_NICHES) <- genotype #group by genotype
    #podotoGECsignals <- FindAllMarkers(Signal_from_podo_GEC, min.pct = 0.05, logfc.threshold = 0.05, only.pos = T)
    #write.csv(podotoGECsignals, "podocyte_GEC_clustering.csv") # Save and export for curation
    
#--------------------------------------------------### END ###--------------------------------------------------#
    