#' @title SCORPIUS Trajectory Analysis for Seurat Object
#' @description Perform SCORPIUS trajectory inference and pseudotime analysis on a Seurat object, 
#'   and plot density along pseudotime for each cell type.
#' @param sce_obj input Seurat object containing single-cell RNA-seq data.
#' @param celltype_col input Column name in \code{sce_obj@meta.data} specifying cell types. Default is "celltype".
#' @param assay Assay to use for expression matrix. Default is "RNA".
#' @param slot Slot in the assay to extract data. Default is "data".
#' @param save_file Optional file path to save results (Rdata). Default is NULL (no saving).
#' @param plot_title Title for the pseudotime density plot. Default is "SCORPIUS Pseudotime".


library(SCORPIUS)
library(preprocessCore)
library(ggplot2)
library(cowplot)

run_SCORPIUS <- function(sce_obj, celltype_col = "celltype", 
                         assay = "RNA", slot = "data",
                         save_file = NULL, plot_title = "SCORPIUS Pseudotime") {
  

  expr_mat <- GetAssayData(sce_obj, assay = assay, slot = slot)
  expr_mat <- as.matrix(expr_mat)
  
  if(nrow(expr_mat) < ncol(expr_mat)) {
    expr_mat <- t(expr_mat)
  }
  
  if(!celltype_col %in% colnames(sce_obj@meta.data)) {
    stop(paste("Column", celltype_col, "not found in sce_obj@meta.data"))
  }
  celltype_vec <- sce_obj@meta.data[match(rownames(expr_mat), rownames(sce_obj@meta.data)), celltype_col]
  

  space <- reduce_dimensionality(expr_mat, dist_method = "spearman")
  draw_trajectory_plot(space, celltype_vec, contour = TRUE)
 
  traj <- infer_trajectory(space)
  draw_trajectory_plot(space, celltype_vec, traj$path, contour = TRUE)
  

  SCOR_time <- data.frame(trajtime = traj$time, ID = names(traj), celltype = celltype_vec)
  
  pseudotime_plot <- ggplot(SCOR_time, aes(x = trajtime, fill = celltype, color = celltype)) +
    geom_density(color = "black") +
    facet_wrap(~ celltype, ncol = 1, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    labs(x = "Pseudotime", y = "Density", title = plot_title) +
    scale_fill_manual(values = rep(c("#EC706E","#E1C767"), length.out = length(unique(celltype_vec))))
  
  print(pseudotime_plot)
  
 
  if(!is.null(save_file)) {
    save(space, SCOR_time, traj, expr_mat, celltype_vec, file = save_file)
    message("Results saved to: ", save_file)
  }
  
  return(list(space = space, traj = traj, pseudotime_plot = pseudotime_plot, SCOR_time = SCOR_time))
}
