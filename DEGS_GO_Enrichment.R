#' @title NorRNA
#' @description Preprocess gene expression matrix by normalization
#' @title GO_enrichment_plot
#' @description Gene Ontology function for bulk RNA
#' @param data Input gene expression matrix of bulk-RNA seq data 
#' @return resdata Normalized gene expression matrix

library(DESeq2)
library(clusterProfiler)
library(dplyr)
library(stringr)
library(ggplot2)

exp_matrix<-fread("path/to/geneexp.txt",header=T)
NorRNA<-function(data){
  database=data[,c(1:ncol(data))];dim(data)
  head(database)
  database=database[complete.cases(database),];dim(database)
  cs1<-length(grep("P",colnames(database)));cs1# label of patients
  cs2<-(ncol(database)-cs1)
  condition<-factor(c(rep("ctrl",cs2),rep("patient",cs1)))
  dds <- DESeqDataSetFromMatrix(database, DataFrame(condition), design= ~ condition)
  dds <- DESeq(dds) 
  res <- results(dds)
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, 
                                                            normalized=TRUE)),by="row.names",sort=FALSE)
  return(resdata)
}
diffExp<-NorRNA(exp_matrix)

diffExp%>%filter(pvalue<0.05,abs(log2FoldChange)>1)->diffGenes#Filter out genes with significantly different expression levels
diffGenes <- diffExp %>% filter(pvalue < 0.05, abs(log2FoldChange) > 1)
gene_symbols <- diffGenes$Row.names 

GO_enrichment_plot <- function(gene_symbols,
                               OrgDb = "org.Hs.eg.db",
                               pCutoff = 0.05,
                               qCutoff = 0.05,
                               padj_filter = 0.05,
                               topN = 10) {
  gene_df <- bitr(gene_symbols,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = OrgDb)
  entrez_ids <- gene_df$ENTREZID
  
  ontologies <- c("BP", "CC", "MF")
  go_list <- lapply(ontologies, function(ont) {
    enrichGO(gene = entrez_ids,
             OrgDb = OrgDb,
             ont = ont,
             pAdjustMethod = "BH",
             pvalueCutoff = pCutoff,
             qvalueCutoff = qCutoff,
             readable = TRUE)@result %>%
      dplyr::filter(p.adjust < padj_filter) %>%
      mutate(ONTOLOGY = ont,
             logP = -log10(p.adjust)) %>%  
      arrange(desc(logP)) %>% 
      slice_head(n = topN)
  })
  
  go_res <- bind_rows(go_list)
  go_res$Description <- str_to_title(go_res$Description)
  go_res$shun <- paste0("c", seq_len(nrow(go_res)))
  
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
    labs(title = "GO Enrichment", x = "", y = "-log10(p.adjust)")
  
  return(list(result = go_res, plot = gobulkPlot))
}
