#' @title Plot cell type proportions
#' @description Calculate and visualize the proportion of each cell type in a Seurat object
#' @param seurat_obj input Seurat object containing cell type metadata in column 'celltype'
#' @param exclude_index default Numeric vector of indices of cell types to exclude (optional)
#' @param fill_color Fill color for bars (default "#EC706E")
#' @return A ggplot object showing cell type proportions

library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)

plotCellProportion <- function(seurat_obj, exclude_index = NULL, fill_color = "#EC706E") {
  
  cell_counts <- table(seurat_obj@meta.data$celltype)
  if(!is.null(exclude_index)){
    cell_counts <- cell_counts[-exclude_index]
  }
  
  CellPRO <- data.frame(
    cell = names(cell_counts),
    prop = round(as.numeric(cell_counts)/sum(as.numeric(cell_counts)), 4)
  ) %>%arrange(desc(prop))
  
  CellPRO$cell <- factor(CellPRO$cell)
  CellPRO$total <- rev(letters[1:nrow(CellPRO)])
  

  p <- ggplot(CellPRO, aes(x = prop, y = total, fill = fill_color)) +
    geom_bar(stat = "identity", width = 0.74) +
    theme_minimal() +
    labs(title = "Proportion of different cells (%)", x = "", y = "") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_y_discrete(labels = rev(as.character(CellPRO$cell))) +
    geom_text(aes(label = prop * 100), vjust = 0.5, hjust = 0.5, color = "black") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(colour = "black", size = 12),
      axis.text.y = element_text(colour = "black", size = 12),
      axis.title.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5, lineheight = 0.2),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

