#' @title plotCelltypeDEVolcano
#' @description Calculate differential expression (patient vs control) per cell type and plot volcano/manhattan style
#' @param seurat_obj input input Seurat object
#' @param condition_col input Column name for condition (patient/control)
#' @param celltype_col input Column name for cell type
#' @param control_ids input Vector of sample IDs for controls
#' @param patient_ids input Vector of sample IDs for patients
#' @param highlight_genes Optional vector of genes to highlight in the plot
#' @param topGeneN Number of top genes for labeling
#' @return List containing diff expression dataframe and plot object

library(dplyr)
library(Seurat)
library(openxlsx)
library(scRNAtoolVis)

plotCelltypeDEVolcano <- function(seurat_obj, 
                                  condition_col = "group", 
                                  celltype_col = "celltype",
                                  control_ids, 
                                  patient_ids, 
                                  highlight_genes = NULL,
                                  topGeneN = 10) {
  
  celltypes <- unique(seurat_obj[[celltype_col]][,1])
  man_CtrPat <- data.frame()
  
  for(ct in celltypes) {
    Ct <- subset(seurat_obj, subset = get(celltype_col) == ct)
    
    ctrl_cells <- rownames(Ct@meta.data %>% filter(get(condition_col) %in% control_ids))
    pat_cells  <- rownames(Ct@meta.data %>% filter(get(condition_col) %in% patient_ids))
    expr_mat <- as.data.frame(GetAssayData(Ct, assay = "RNA")[, c(ctrl_cells, pat_cells)])
    expr_mat$Low <- apply(expr_mat, 1, function(x) sum(x==0))
    expr_mat <- expr_mat[-which(expr_mat$Low > ceiling(ncol(expr_mat)*0.7)), ]
    expr_mat$Low <- NULL
    DifWil <- t(apply(expr_mat, 1, function(x) {
      y1 <- as.numeric(x[1:length(ctrl_cells)])
      y2 <- as.numeric(x[(length(ctrl_cells)+1):ncol(expr_mat)])
      p <- wilcox.test(y1, y2, paired = FALSE)$p.value
      fd <- log2(mean(y2)/mean(y1))
      pct1 <- round(sum(y1>0)/length(y1),2)
      pct2 <- round(sum(y2>0)/length(y2),2)
      return(c(p, fd, pct1, pct2))
    }))
    DifWil <- as.data.frame(DifWil)
    colnames(DifWil) <- c("p_val","avg_log2FC","pct.1","pct.2")
    DifWil$p_val_adj <- p.adjust(DifWil$p_val, method = "BH")
    DifWil$cluster <- ct
    DifWil$gene <- rownames(DifWil)
    man_CtrPat <- rbind(man_CtrPat, DifWil)
  }
  
  man_CtrPat$label <- ifelse(man_CtrPat$gene %in% highlight_genes, 0.0001, NA)
  man_plot <- jjVolcano(diffData = man_CtrPat, 
                        topGeneN = topGeneN,
                        back.col = "white",
                        aesCol = c("#4466B0","#E71F19"),
                        tile.col = rep("grey", length(unique(celltypes))),
                        col.type = "adjustP",
                        myMarkers = highlight_genes,
                        adjustP.cutoff = 0.05,
                        celltypeSize = 4
  ) +
    xlab(NULL) + ylab("Log2FC (patient vs control)") +
    ylim(-3,3) +
    geom_hline(yintercept = 1, linetype="dashed", color="black")
  
  return(list(diff_df = man_CtrPat, plot = man_plot))
}

