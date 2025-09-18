#' @title Plot DPT Pseudotime Density for Mono/Macro
#' @description Generate density plot of DPT pseudotime for monocyte and macrophage-like cells.
#' @param data_input input Data frame or CSV file path containing pseudotime results. Must have columns: "celltype" and "dpt_pseudotime".
#' @param monocyte_label input Character, label for monocyte cells. Default: "monocyte".
#' @param macro_label input Character, label for macrophage-like cells. Default: "macrophage".
#' @param fill_colors Vector of two colors for monocyte and macrophage, default c("#EC706E","#E1C767").
#' @return ggplot object of pseudotime density plot.

library(ggplot2)

plot_DPT_density <- function(data_input,
                             monocyte_label = "monocyte",
                             macro_label = "macrophage_like",
                             fill_colors = c("#EC706E","#E1C767")) {
  
  
  if(is.character(data_input)) {
    data <- read.csv(data_input, stringsAsFactors = FALSE)
  } else {
    data <- data_input
  }
  data$celltype <- factor(data$celltype)
  data$chara_order <- ifelse(data$celltype == "monocyte", "A", "B")
  p <- ggplot(data, aes(x = dpt_pseudotime, fill = chara_order, color = chara_order)) +
    geom_density(color = "black") +
    scale_fill_manual(values = fill_colors) +
    facet_wrap(~ chara_order, ncol = 1, scales = "free_y",
               labeller = labeller(chara_order = c("A" = monocyte_label, "B" = macro_label))) +
    theme_minimal() +
    theme(legend.position = "none",
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    labs(x = "DPT", y = "Density")
  
  return(p)
}


library(Seurat)
library(SeuratDisk)

#Save the Seurat object and convert it to h5ad format for DPT calculation in Python
SaveH5Seurat(sceN1, filename = "sceN1.h5Seurat")
Convert("sceN1.h5Seurat", dest = "h5ad")