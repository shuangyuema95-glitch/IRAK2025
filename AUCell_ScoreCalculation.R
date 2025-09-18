#' @title Plot response score boxplots for selected cell types
#' @description Subset Seurat object by condition and cell type, and plot response scores per sample/case
#' @param seurat_obj input Seurat object containing metadata with score column
#' @param score_col input Character, name of the column in metadata containing response score
#' @param celltypes input Vector of cell types to plot
#' @param conditions input Vector of conditions to plot separately, e.g. c("basal","LPS")
#' @param sample_map input Named vector to recode sample names to case labels
#' @param colors input Named list of color vectors for each condition
#' @param comparisons input List of comparisons for geom_signif
#' @param ncol input Number of columns in the combined plot
#' @return A named list of ggplot objects per condition and combined annotated figure

library(Seurat)
library(ggplot2)
library(cowplot)
library(ggsignif)
library(dplyr)

plotResponseBox <- function(seurat_obj, score_col, celltypes, conditions, 
                            sample_map, colors, comparisons, ncol = 2){
  
  results <- list()
  for(cond in conditions){
    # Subset metadata for condition
    meta <- seurat_obj@meta.data
    meta <- meta[meta$condition == cond, ]
    
    # Recode sample/case
    meta$case <- recode(meta$sample, !!!sample_map)
    
    # Prepare list to store per-celltype plots
    box_list <- list()
    for(cell in celltypes){
      DATA <- meta[meta$celltype == cell, ]
      
      ymax <- max(DATA[[score_col]], na.rm = TRUE)
      
      box_list[[cell]] <- ggplot(DATA, aes(x = case, y = .data[[score_col]], fill = case)) +
        geom_boxplot(size=0.7, outlier.shape = NA, width=0.55) +
        labs(x = "", y = "") +
        labs(title = cell) +
        scale_fill_manual(values = colors[[cond]]) +
        theme_classic() +
        theme(
          legend.position = "none",
          axis.text.x = element_text(colour="black", size=10),
          axis.text.y = element_text(colour="black", size=10),
          axis.title.y = element_text(size=12),
          plot.title = element_text(size=12, hjust = 0.5, lineheight=0.2),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        geom_signif(comparisons = comparisons, step_increase = 0.1,
                    map_signif_level = TRUE, y_position = ymax, textsize = 3,
                    parse = FALSE, color = "black", tip_length = 0)
    }
    
    # Combine plots
    combined <- plot_grid(plotlist = box_list, ncol = ncol, align = "hv", byrow = FALSE)
    annotated <- annotate_figure(combined, left = text_grob(paste(score_col, "(Z-score)"),
                                                            size = 10, face = "bold", rot = 90))
    
    results[[cond]] <- list(per_celltype = box_list, combined = annotated)
  }
  
  return(results)
}

#-------------------------------
# Example usage:

celltypes <- c("CD8 T","memory CD4 T","naive CD4 T","monocyte","macrophage_like","pDC")
sample_map <- c("C1_LPS"="ctrl","C2_LPS"="ctrl","C3_LPS"="ctrl",
                "P3_LPS"="P3","P4_LPS"="P4")  
colors <- list(
  basal = c("ctrl"="#1F4FA1","C1"="#E8652B","C2"="#E46266"),  
  LPS   = c("ctrl"="#40ABE2","P3"="#E8652B","P4"="#E46266")
)

comparisons <- list(c("ctrl","P3"), c("ctrl","P4"))

# Call function
response_plots <- plotResponseBox(seurat_obj = sceN1,
                                  score_col = "aucScore",
                                  celltypes = celltypes,
                                  conditions = conditions,
                                  sample_map = sample_map,
                                  colors = colors,
                                  comparisons = comparisons,
                                  ncol = 2)

# View combined plot for LPS
response_plots$LPS$combined
# View combined plot for basal
response_plots$basal$combined
