#' @title Single-cell clustering marker analysis
#' @description Perform cluster marker identification and conserved marker extraction.
#' @param seurat_obj input Integrated Seurat object with clustering (seurat_clusters) and sample info
#' @param marker_table input Data frame containing gene-to-celltype annotation (columns: gene, celltype)
#' @param annotations Optional data frame with gene descriptions (columns: gene_name, description)
#' @retun A list contains different marker gene results by distinct Seurat functions

library(Seurat)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(cowplot)

analyzeClusterMarkers <- function(seurat_obj, marker_table, annotations = NULL) {
  all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  
  get_conserved <- function(cluster_id) {
    if(sum(seurat_obj$seurat_clusters == cluster_id) < 3) return(NULL)
    
    markers <- FindConservedMarkers(
      seurat_obj,
      ident.1 = cluster_id,
      grouping.var = "sample",
      only.pos = TRUE
    ) %>%
      rownames_to_column(var = "gene")
    
    if(!is.null(annotations)) {
      markers <- markers %>%
        left_join(annotations[, c("gene_name", "description")],
                  by = c("gene" = "gene_name"))
    }
    
    markers <- cbind(cluster_id = cluster_id, markers)
    return(markers)
  }
  
  # Build conserved marker list for all clusters
  conserved_markers <- list()
  clusters <- unique(seurat_obj$seurat_clusters)
  for(clu in clusters) {
    conserved_markers[[as.character(clu)]] <- get_conserved(clu)$gene
  }
  
  # Function to count marker frequency by cell type for standard markers
  count_marker_fn <- function(cluster) {
    gene <- all_markers[all_markers$cluster == cluster, 'gene']
    freq_table <- marker_table %>%
      filter(gene %in% gene) %>%
      count(celltype) %>%
      arrange(desc(n))
    return(freq_table)
  }
  
  # Function to count marker frequency by cell type for conserved markers
  count_marker2_fn <- function(cluster) {
    gene <- conserved_markers[[as.character(cluster)]]
    freq_table <- marker_table %>%
      filter(gene %in% gene) %>%
      count(celltype) %>%
      arrange(desc(n))
    return(freq_table)
  }
  
  return(list(
    all_markers = all_markers,
    conserved_markers = conserved_markers,
    count_marker_fn = count_marker_fn,
    count_marker2_fn = count_marker2_fn
  ))
}

markerList<-analyzeClusterMarkers(sceN1,marker_table)


#plot

colormaps<- c("#4466B0","#74B47F","#F1D19F","#A982B2","#40ABE2","#E2A9C9",
           "#0F9296","#F1BAAB","#EC706E","#AF478A","#84CAE8","#E1C767",
           "#DBD8E7","#C1E1DA","#25853A")

#-----------------------------
# UMAP / DimPlot for all cells colored by celltype
p3_basal <- DimPlot(sceN1,
                    label = FALSE,
                    group.by = "celltype",
                    pt.size = 0.001,
                    cols = colormaps,
                    raster = FALSE) +
  theme(plot.title = element_blank()) +
  theme(
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.x = element_text(size = 22, color = "black"),
    axis.title.y = element_text(size = 22, color = "black"),
    legend.text = element_text(size = 24)
  )
p3_basal

#-----------------------------
# DimPlot colored by sample group
DimPlot(sceN1, label = FALSE, group.by = "group", pt.size = 0.001, raster = FALSE)
group_colors <- c("#84CAE8", "#40ABE2", "#EC706E", "#DBD8E7", "#F08619")
DimPlot(sceN1, label = FALSE, group.by = "group", pt.size = 0.001,
        raster = FALSE, cols = group_colors)

#-----------------------------
# Subset cells by condition
sceN1LPS   <- subset(sceN1, subset = condition %in% "LPS")
sceN1basal <- subset(sceN1, subset = condition %in% "basal")
sample_subset <- c("C1", "C2", "C3", "P3", "P4")
sample_plots <- list()
for(i in sample_subset){
  da <- subset(sceN1basal, subset = group %in% i)
  sample_plots[[i]] <- DimPlot(da,
                               label = FALSE,
                               group.by = "celltype",
                               pt.size = 0.1,
                               cols = yanse) +
    labs(title = i) +
    theme(plot.title = element_text(size = 12)) +
    NoLegend()
}
# Combine plots in a grid
plot_grid(plotlist = sample_plots,
          ncol = length(sample_subset),
          align = "hv",
          rel_widths = rep(1, length(sample_subset)))


