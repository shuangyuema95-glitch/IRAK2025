#' @title Slingshot trajectory analysis for Mono/Macro
#' @description Perform Slingshot pseudotime analysis on monocyte/macrophage_like cells and plot density along pseudotime
#' @param seurat_obj input Seurat object
#' @param cell_types input Character vector of cell types to include, e.g., c("monocyte","macrophage_like")
#' @param umap_slot Slot of UMAP embeddings, default "umap"
#' @param colors Named vector of colors for cell types, default c("#EC706E","#E1C767")
#' @return List containing Slingshot object and ggplot density plot

library(slingshot)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
runSlingshotMonoMacro <- function(seurat_obj, 
                                  cell_types = c("monocyte","macrophage_like"),
                                  umap_slot = "umap",
                                  colors = c("monocyte"="#EC706E","macrophage_like"="#E1C767")) {

  
  sce_sub <- subset(seurat_obj, subset = celltype %in% cell_types)
    sce_sce <- as.SingleCellExperiment(sce_sub, assay = "RNA")
    sce_slingshot <- slingshot(sce_sce,
                             reducedDim = toupper(umap_slot),
                             clusterLabels = sce_sce$celltype,
                             approx_points = 1500)
  lineages <- slingPseudotime(sce_slingshot)
  lineages_long <- as.data.frame(lineages) %>%
    mutate(cell_id = rownames(lineages)) %>%
    pivot_longer(-cell_id, names_to = "lineage", values_to = "pseudotime")
  
  umap_coords <- as.data.frame(reducedDims(sce_sce)[[toupper(umap_slot)]])
  umap_coords$cell_id <- rownames(umap_coords)
  umap_coords$celltype <- sce_sub@meta.data$celltype
  
  merged_data <- left_join(umap_coords, lineages_long, by = "cell_id")
  merged_data <- merged_data %>% filter(celltype %in% cell_types, !is.na(pseudotime))
  
  density_plot <- ggplot(merged_data, aes(x = pseudotime, fill = celltype, color = celltype)) +
    geom_density(color="black") +
    scale_fill_manual(values = colors) +
    facet_wrap(~ celltype, ncol = 1, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    labs(x = "Pseudotime", y = "Density")
  
  return(list(slingshot_obj = sce_slingshot, density_plot = density_plot))
}
