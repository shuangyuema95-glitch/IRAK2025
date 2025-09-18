#' @title GO enrichment for single-cell by cell type
#' @description Perform differential expression between patient and control for a given cell type and visualize top GO terms
#' @param seurat_obj input Seurat object with metadata column 'status' indicating 'patient' or 'control'
#' @param celltypes input Vector of cell types to analyze
#' @param OrgDb OrgDb object or string, default "org.Hs.eg.db"
#' @param padj_filter Adjusted p-value cutoff for GO terms
#' @param topN Number of top GO terms to plot per cell type
#' @return A list containing per-cell type GO result table and ggplot

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)

GO_enrichment_scRNA_patient <- function(seurat_obj, celltypes, 
                                        OrgDb = "org.Hs.eg.db",
                                        padj_filter = 0.05,
                                        topN = 10) {
 
  go_results <- list()
  go_plots <- list()
  
  for(ct in celltypes){
    message("Processing cell type: ", ct)
    
    # Subset by cell type
    sub_obj <- subset(seurat_obj, subset = CT == ct)
    
    # Set status as Idents
    Idents(sub_obj) <- sub_obj$status
    
    # Differential expression: patient vs control
    diff_expr <- FindMarkers(sub_obj, ident.1 = "patient", ident.2 = "control")
    diff_expr <- diff_expr %>% 
      rownames_to_column(var = "symbol") %>%
      filter(!is.na(avg_log2FC)) 
    
    # Convert SYMBOL to ENTREZID
    gene_df <- bitr(diff_expr$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
    diff_expr <- inner_join(diff_expr, gene_df, by = c("symbol" = "SYMBOL"))
    
    # Prepare fold change vector named by ENTREZID
    fc_vec <- diff_expr$avg_log2FC
    names(fc_vec) <- diff_expr$ENTREZID
    fc_vec <- sort(fc_vec, decreasing = TRUE)
    
    # GO enrichment (BP, CC, MF)
    ontologies <- c("BP", "CC", "MF")
    go_list <- lapply(ontologies, function(ont) {
      enrichGO(gene = names(fc_vec),
               OrgDb = OrgDb,
               ont = ont,
               pAdjustMethod = "BH",
               pvalueCutoff = 1,
               qvalueCutoff = 1,
               readable = TRUE)@result %>%
        filter(p.adjust < padj_filter) %>%
        mutate(ONTOLOGY = ont,
               logP = -log10(p.adjust)) %>%
        arrange(desc(logP)) %>%
        slice_head(n = topN)
    })
    
    go_res <- bind_rows(go_list)
    go_res$Description <- str_to_title(go_res$Description)
    go_res$shun <- paste0("c", seq_len(nrow(go_res)))
    
    # Plot
    ont_colors <- c("BP" = "#F1BAAB", "CC" = "#A982B2", "MF" = "#4466B0")
    gobulkPlot <- ggplot(go_res, aes(x = shun, y = logP, fill = ONTOLOGY)) +
      geom_bar(stat = 'identity', width = 0.72, position = position_dodge(0.63)) +
      scale_fill_manual(values = ont_colors) +
      scale_x_discrete(labels = go_res$Description) +
      coord_flip() +
      theme_classic() +
      theme(
        legend.position = "right",
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      labs(title = paste0(ct, " GO Enrichment"), x = "", y = "-log10(p.adjust)")
    
    go_results[[ct]] <- go_res
    go_plots[[ct]] <- gobulkPlot
  }
  
  return(list(result = go_results, plot = go_plots))
}
