#' @title Plot marker gene expressions
#' @param seurat_obj input Seurat object containing cell type metadata in column 'celltype'
#' @param features input marker genes
#' @return bubble_plot

library(Seurat)
library(ggplot2)

plot_LPS_bubble <- function(seurat_obj, features, cols=c("lightgrey","red")){
  
  bubble_plot <- DotPlot(
    seurat_obj,
    features = unique(features),
    assay = "RNA",
    cols = cols
  ) +
    RotatedAxis() +
    theme(
      axis.line = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(angle = 45, vjust = 0.5, hjust = 1),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      panel.border = element_rect(color = "black", size = 1, fill = NA)
    ) +
    ylab(NULL) + xlab(NULL)
  
  return(bubble_plot)
}

#--------------------------------------
# Example usage:
features2 <- c("CD14","S100A8","S100A12","CLEC12A","TCF7L2","FCGR3A",
               "FCN1","CD300E","CTSC","GLUL",
               "HES1","IL1B","CD68","HIF1A","EGR1","ANPEP","VEGFA",
               "NLRP3","MNDA","CXCL8")

bubble_LPSref <- plot_LPS_bubble(sceN1LPS, features2, cols=cols)
bubble_LPSref
