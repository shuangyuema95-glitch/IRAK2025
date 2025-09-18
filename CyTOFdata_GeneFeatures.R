#' @title Plot Gene Expression on t-SNE
#' @description Generate facet t-SNE plots for a list of target genes, colored by expression, and save as PDF.
#' @param tsne_obj input t-SNE object (e.g., from CATALYST runDR).
#' @param genes input Character vector of target gene names or regex patterns to match genes.
#' @param facet_col Column in tsne_obj metadata to facet by. Default is "sample_order".
#' @param ncol Number of columns in facet wrap. Default is 6.
#' @param palette Color palette for expression. Default is rev(hcl.colors(20, "Spectral")).
#' @param output_dir Directory to save PDF files. Default is current working directory.
#' @return A named list of ggplot objects for each gene.

library(cytofkit2)
library(CATALYST)
library(dplyr)
library(flowCore)
library(cytofWorkflow)
library(ggplot2)
library(mvtnorm)
library(patchwork)
library(reshape2)
library(pander)
library(tidyverse)
library(dplyr)
library(do)


plot_gene_tsne <- function(tsne_obj, genes, facet_col = "sample_order", ncol = 6,
                           palette = rev(hcl.colors(20, "Spectral")),
                           output_dir = ".") {
  
  gene_feature <- list()
  
  target_genes <- names(tsne_obj@rowRanges)[grep(paste(genes, collapse="|"), names(tsne_obj@rowRanges))]
  
  if(length(target_genes) == 0) {
    stop("No matching genes found in tsne_obj rowRanges.")
  }
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  for(i in seq_along(target_genes)) {
    gene <- target_genes[i]
    
    p <- plotDR(tsne_obj, "TSNE", color_by = gene, facet_by = facet_col, ncol = ncol,
                a_pal = palette) +
      xlab(NULL) + ylab(NULL) +
      theme(
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.4, 'cm'),
        legend.key.width = unit(0.4, 'cm'),
        legend.position = "right",
        plot.title = element_blank(),
        axis.line.y = element_line(size = 0.45),
        panel.background = element_rect(color = "black", size = 0.45),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()
      ) +
      ylab(gene)
    
    ggsave(filename = file.path(output_dir, paste0(gene, "_gene_tSNE.pdf")),
           plot = p, width = 14.91, height = 8.54, units = "in")
    
    gene_feature[[gene]] <- p
  }
  
  return(gene_feature)
}
