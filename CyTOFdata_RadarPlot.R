#' @title Generate Radar Plot from CyTOF t-SNE Object
#' @description Calculate cell type proportions per sample from a CyTOF t-SNE object and plot radar plot.
#' @param tsne_obj input CyTOF t-SNE object (SingleCellExperiment from CATALYST).
#' @param cell_order input Character vector specifying cell type order for radar plot. Default: c("B cell","CD4T","CD8T","DC","monocyte","Treg").
#' @param grid_min Minimum radar grid value. Default: 0.
#' @param grid_mid Middle radar grid value. Default: 0.05.
#' @param grid_max Maximum radar grid value. Default: 0.12.
#' @param palette Vector of colors for each sample. Default: NULL.
#' @param output_dir Directory to save PDF. Default: ".".
#' @param plot_name File name (without extension) for PDF. Default: "radar_plot".
#' @return List contains cell frequency matrix and radar plot

library(ggradar)
library(ggplot2)
library(dplyr)

plot_cytof_radar <- function(tsne_obj,
                             cell_order = c("B cell","CD4T","CD8T","DC","monocyte","Treg"),
                             grid_min = 0, grid_mid = 0.05, grid_max = 0.12,
                             palette = NULL,
                             output_dir = ".",
                             plot_name = "radar_plot") {
  
  library(dplyr)
  library(ggradar)
  
  CountCell <- as.data.frame(table(tsne_obj$sample_id, tsne_obj$celltype))
  
  bili <- data.frame()
  for(s in unique(CountCell$Var1)) {
    DAT1 <- CountCell %>% dplyr::filter(Var1 == s)
    DAT1$prop <- DAT1$Freq / sum(DAT1$Freq)
    bili <- rbind(bili, DAT1)
  }
  
 
  radar_df <- data.frame()
  for(s in unique(bili$Var1)) {
    D <- bili[bili$Var1 == s, ]
    D <- D[match(cell_order, D$Var2), ] # reorder cells
    radar_df <- rbind(radar_df, D$prop)
  }
  colnames(radar_df) <- cell_order
  radar_df$sample <- unique(bili$Var1)
  radar_df <- radar_df[, c(ncol(radar_df), 1:(ncol(radar_df)-1))] # sample first column
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  radar_plot <- ggradar(radar_df,
                        grid.min = grid_min,
                        grid.mid = grid_mid,
                        grid.max = grid_max,
                        group.point.size = 3,
                        group.line.width = 1.5,
                        gridline.min.colour = "grey",
                        gridline.mid.colour = "black",
                        gridline.max.colour = "black",
                        axis.label.size = 5,
                        axis.line.colour = "grey",
                        legend.text.size = 10,
                        legend.position = "left",
                        background.circle.colour = "white",
                        background.circle.transparency = 0.1,
                        group.colours = palette)
  
  ggsave(filename = file.path(output_dir, paste0(plot_name, ".pdf")),
         plot = radar_plot, width = 8, height = 6)
  
  return(list(radar_df = radar_df, radar_plot = radar_plot))
}