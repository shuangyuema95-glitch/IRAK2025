#' @title plotPseudobulkHeatmap
#' @description For each cell type, compute pseudobulk (AverageExpression) per sample and plot heatmap.
#' @param seurat_obj input Seurat object
#' @param celltypes input Vector of cell types to iterate
#' @param gene_list input Vector of genes to plot (will include CXCL8 automatically)
#' @param sample_order input  Character vector of samples (column order in heatmap)
#' @pathway input The pathway names corresponding to the genes used in the heatmap
#' @param assay Assay to use (default "RNA")
#' @param slot Slot to use (default "data")
#' @param colors Color palette for heatmap (default Blue-White-Red)
#' @param breaks Numeric vector of breaks (optional)
#' @param show_gene_names Whether to plot gene names (default TRUE)
#' @param fontsize Font size for heatmap (default 3.5)
#' @return A list with two elements: $with_gene_names and $without_gene_names

library(Seurat)
library(pheatmap)
library(ggplotify)

plotPseudobulkHeatmap <- function(seurat_obj,
                                  celltypes,
                                  gene_list,
                                  sample_order,
                                  pathway,
                                  assay = "RNA",
                                  slot = "data",
                                  colors = colorRampPalette(c("Blue","White","Red"))(101),
                                  breaks = seq(-2, 2, 0.5),
                                  show_gene_names = TRUE,
                                  fontsize = 3.5) {
  genes_use <- unique(gene_list)
  res_with <- list()
  res_without <- list()
  for (i in seq_along(celltypes)) {
    cell_i <- celltypes[i]
    message("Processing cell type: ", cell_i)
    
    # subset
    h <- subset(seurat_obj, subset = CT %in% cell_i)
    h$CT <- factor(h$CT)
    
    # AverageExpression (pseudobulk)
    h_avg <- AverageExpression(h, group.by = "sample", assays = assay, slot = slot)
    h_mat <- as.data.frame(h_avg[[assay]])
    h_mat <- h_mat[match(genes_use, rownames(h_mat)), match(sample_order, colnames(h_mat))]
    res_with[[cell_i]] <- as.ggplot(pheatmap(
      h_mat,
      scale = "row",
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      show_rownames = TRUE,
      border_color = NA,
      main = paste(cell_i, "pathway", sep = ""),
      fontsize = fontsize,
      legend = TRUE,
      treeheight_col = 0,
      treeheight_row = 0,
      legend_breaks = seq(-2, 2, 0.5),
      color = colors,
      annotation_legend = TRUE,
      breaks = breaks,
      cellwidth = 8,
      cellheight = 4,
      clustering_distance_rows = "correlation",
      clustering_method = "average"
    ))
    
    res_without[[cell_i]] <- as.ggplot(pheatmap(
      h_mat,
      scale = "row",
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      show_rownames = FALSE,
      border_color = NA,
      main = paste(cell_i, "pathway", sep = ""),
      fontsize = fontsize,
      legend = TRUE,
      treeheight_col = 0,
      treeheight_row = 0,
      legend_breaks = seq(-2, 2, 0.5),
      color = colors,
      annotation_legend = TRUE,
      breaks = breaks,
      cellwidth = 8,
      cellheight = 4,
      clustering_distance_rows = "correlation",
      clustering_method = "average"
    ))
  }
  
  return(list(with_gene_names = res_with, without_gene_names = res_without))
}
