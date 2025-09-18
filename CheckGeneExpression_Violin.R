#' @title Generate violin plots for selected genes in specific cell types
#' @description Subset Seurat object by cell type, split by sample/case, and plot expression for multiple genes under different conditions.
#' @param seurat_obj input Seurat object containing RNA assay and metadata columns 'celltype', 'group', 'condition'
#' @param genes input Vector of gene names to plot
#' @param celltypes input Vector of cell types to include
#' @param cases input Named vector to recode sample/group to case labels, e.g. c('ZCH'='P1','MR'='P2')
#' @param output_dir input Directory where plots will be saved
#' @param conditions Vector of conditions to plot separately, e.g. c('basal','LPS')
#' @param colors List of color vectors for each condition (should match number of cases)
#' @return Invisible list of ggplot objects


library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)

plotViolinGenes <- function(seurat_obj, genes, celltypes, cases, output_dir, 
                            conditions = c('basal','LPS'), colors = NULL) {
  # Add recoded case and simplified cell type info
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    mutate(case = recode(group, !!!cases),
           CT = recode(celltype,
                       "naive CD4 T" = "T cell",
                       "CD8 T" = "T cell",
                       "memory CD4 T" = "T cell"))
  
  # Subset object to selected cell types
  seurat_subset <- subset(seurat_obj, subset = CT %in% celltypes)
  seurat_subset$CT <- factor(seurat_subset$CT)
  
  violin_list <- list()
  
  # Loop over genes
  for(gene in genes){
    if(gene %in% rownames(seurat_subset[['RNA']])){
      for(cond in conditions){
        sub_obj <- subset(seurat_subset, subset = condition %in% cond)
        
        col_vector <- if(!is.null(colors) && !is.null(colors[[cond]])) colors[[cond]] else rep("grey50", length(unique(sub_obj$case)))
        
        vln <- VlnPlot(sub_obj, features = gene, group.by = "CT", split.by = "case", pt.size = 0,
                       cols = col_vector) +
          xlab(NULL) + ylab(NULL) + ggtitle(gene) +
          theme(
            legend.position = "none",
            axis.text.y = element_text(size = 10),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.title = element_text(size = 12, hjust = 0.5, lineheight = 0.2),
            panel.border = element_rect(color = "black", size = 0.5, fill = NA)
          )
        
        # Create condition-specific output folder
        cond_dir <- file.path(output_dir, cond)
        if(!dir.exists(cond_dir)) dir.create(cond_dir, recursive = TRUE)
        
        ggsave(vln, file = file.path(cond_dir, paste0(gene, ".pdf")),
               width = 3.56, height = 3.53)
        
        violin_list[[paste0(gene, "_", cond)]] <- vln
      }
    }
  }
  
  return(invisible(violin_list))
}

#-------------------------------
# Example usage:
# Define genes, cell types, and case mapping
Viogene <- c("STAT5A","DDX58","DDX60","IRF2","IFIT2","S100A12","REL","LAMP3","IFIH1","TAP1","IRF1","IFI6")
interested_celltype <- c("macrophage")
colors_list <- list(
  basal = c("#40ABE2","#84CAE8","#DBD8E7","#E8652B","#E46266"),
  LPS   = c("#5DB7D1","#84CAE8","#DBD8E7","#E8652B","#E46266")
)
# Call function
violin_results <- plotViolinGenes(sceN1, genes = Viogene, celltypes = interested_celltype,
                                  cases = cases,
                                  output_dir = "E:/IRAK2/ViolinPlots",
                                  conditions = c("basal","LPS"),
                                  colors = colors_list)