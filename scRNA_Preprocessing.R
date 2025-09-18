#' @title Read, preprocess, and integrate single-cell RNA-seq data
#' @description Reads multiple 10X datasets, performs QC, normalization, variable feature selection,PCA, integration, and Harmony batch correction.
#' @param samples input Character vector of sample names corresponding to folders containing 10X data
#' @param data_dir input String path to the directory containing sample folders
#' @return integrated The object processed through the Seurat pipeline

library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)

preprocessIntegrateSC <- function(samples, data_dir) {
  
  sce_list <- list()
  for (i in seq_along(samples)) {
    sample_path <- file.path(data_dir, samples[i])
    s <- CreateSeuratObject(counts = Read10X(sample_path), min.cells = 3, min.features = 200)
    
    s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
    s$sample <- samples[i]
    s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000)
    s <- FindVariableFeatures(s, selection.method = "vst")
    s <- ScaleData(s)
    s <- RunPCA(s)
    sce_list[[i]] <- s
  }
  names(sce_list) <- samples
  
  # Integration
  anchors <- FindIntegrationAnchors(object.list = sce_list)
  integrated <- IntegrateData(anchorset = anchors)
  DefaultAssay(integrated) <- "integrated"
  integrated <- ScaleData(integrated)
  integrated <- RunPCA(integrated)
  counts <- integrated@assays$RNA@counts
  rb_genes <- rownames(counts)[grepl("^RP[SL]", rownames(counts), ignore.case = TRUE)]
  percent_ribo <- Matrix::colSums(counts[rb_genes, ]) / Matrix::colSums(counts) * 100
  integrated <- AddMetaData(integrated, percent_ribo, col.name = "percent.ribo")
  
  # Quality control: filter cells
  integrated <- subset(integrated,
                       subset = nFeature_RNA > 300 & nFeature_RNA < 7000 &
                         nCount_RNA < 100000 & percent.mt < 10)
  
  # Harmony batch correction
  integrated <- RunHarmony(integrated, reduction = "pca", group.by.vars = "sample", reduction.save = "harmony")
  integrated <- RunPCA(integrated)
  
  # Clustering and UMAP
  integrated <- FindNeighbors(integrated, dims = 1:50)
  integrated <- FindClusters(integrated, resolution = 0.2)
  integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:50, reduction.name = "umap")
  
  return(integrated)
}

samples<-c("C1_UNS","C2_UNS","C3_UNS","P3_UNS","P4_UNS","C1_LPS","C2_LPS","C3_LPS","P3_LPS","P4_LPS")
sceN1<-preprocessIntegrateSC(samples,datafold)# datafold is the path to 10X genomic files

