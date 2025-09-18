#' @title GSEA_for_CellType
#' @description For a specific cell type, perform differential expression and GSEA.
#' @param seurat_obj input Seurat object
#' @param celltype input Character, specific cell type to analyze
#' @param ident.1 input Character, first group in FindMarkers (e.g., "LPS")
#' @param ident.2 input Character, second group (e.g., "basal")
#' @param assay Assay to use (default "RNA")
#' @param slot Slot to use (default "data")
#' @param minGSSize Minimal gene set size for GSEA
#' @param maxGSSize Maximal gene set size for GSEA
#' @param pvalueCutoff p-value cutoff for GSEA
#' @return A list containing: diff_expr (data.frame), gsea_res (enrichResult), gsea_plots (list of gseaplot2 objects)

library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)

GSEA_for_CellType_All <- function(seurat_obj,
                                  celltype,
                                  ident.1,
                                  ident.2,
                                  assay = "RNA",
                                  slot = "data",
                                  minGSSize = 10,
                                  maxGSSize = 500,
                                  pvalueCutoff = 1) {
  
  sub_obj <- subset(seurat_obj, subset = CT %in% celltype)
  
  Idents(sub_obj) <- sub_obj$condition
  diff_expr <- FindMarkers(sub_obj, ident.1 = ident.1, ident.2 = ident.2, assay = assay, slot = slot)
  diff_expr <- diff_expr %>% rownames_to_column(var = "symbol") %>% filter(!is.na(avg_log2FC))
  
  gene_entrez <- bitr(diff_expr$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_entrez <- distinct(gene_entrez, SYMBOL, .keep_all = TRUE)
  
  diff_expr_all <- diff_expr %>% inner_join(gene_entrez, by = c("symbol" = "SYMBOL")) %>% arrange(desc(avg_log2FC))
  
  fc_list <- diff_expr_all$avg_log2FC
  names(fc_list) <- diff_expr_all$ENTREZID
  
  m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene)
  

  gsea_res <- GSEA(fc_list,
                   TERM2GENE = m_t2g,
                   minGSSize = minGSSize,
                   maxGSSize = maxGSSize,
                   pvalueCutoff = pvalueCutoff,
                   pAdjustMethod = "BH")
  
  gsea_plots <- list()
  if(length(gsea_res) > 0){
    for(i in seq_len(nrow(gsea_res@result))){
      gsea_plots[[i]] <- gseaplot2(gsea_res, geneSetID = i, color = "blue", title = gsea_res@result$Description[i])
    }
  }
  
  return(list(diff_expr = diff_expr_all, gsea_res = gsea_res, gsea_plots = gsea_plots))
}