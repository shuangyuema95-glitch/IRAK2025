#' @title Run CytoTRACE on single-cell data
#' @description Perform CytoTRACE analysis on a Seurat object for selected cells
#' @param seurat_obj input Seurat object containing RNA assay
#' @param celltype_col input Column name in metadata for cell types (default: "celltype")
#' @param selected_cells Vector of cell types to include (default: all)
#' @param reduction_method Dimensionality reduction for plotting, e.g., "umap" or "pca" (default: "umap")
#' @return A list containing CytoTRACE results and ggplot object


library(CytoTRACE)
library(Seurat)
library(ggplot2)

runCytoTRACE <- function(seurat_obj, 
                         celltype_col = "celltype", 
                         selected_cells = NULL, 
                         reduction_method = "umap") {
  
  if (!is.null(selected_cells)) {
    seurat_obj <- subset(seurat_obj, subset = get(celltype_col) %in% selected_cells)
  }
  
  expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  expr_matrix <- as.data.frame(expr_matrix)
  
  phenot <- as.character(seurat_obj[[celltype_col]][,1])
  names(phenot) <- rownames(seurat_obj@meta.data)
  
  results <- CytoTRACE(expr_matrix)
  
  if (reduction_method %in% names(seurat_obj@reductions)) {
    emb <- Embeddings(seurat_obj, reduction = reduction_method)
  } else {
    stop(paste0("Reduction method '", reduction_method, "' not found in Seurat object"))
  }
  
  cyt_plot <- plotCytoTRACE(results, phenotype = phenot, emb = emb)
  return(list(cyto_results = results, plot = cyt_plot))
}

