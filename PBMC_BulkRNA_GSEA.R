#' @title NorRNA
#' @description Preprocess gene expression matrix by normalization
#' @title bulkRNA_GSEA
#' @description Geneset Enrichment Analysis function
#' @title plot_GSEA_barplot
#' @description Barplot for GSEA results with significan
#' @gsea_result input GSEA result by function bulkRNA_GSEA
#' @param data Input gene expression matrix of bulk-RNA seq data 
#' @param signaling Input normalized gene expression matrix by function 'NorRNA',e.g. signaling="Type I interferon"
#' @return resdata Normalized gene expression matrix
#' @return gsea_res2 GSEA results
#' @retutn GSEAbarplot  GSEA bar plot

library(DESeq2)
library(data.table)
library(dplyr)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)

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


bulkRNA_GSEA<-function(diffExp){
diffExp[which(!is.na(diffExp$log2FoldChange)),]->geneFC_df
geneFC_list=as.numeric(geneFC_df$log2FoldChange)
names(geneFC_list)<-geneFC_df$Row.names

entrezid_list <- bitr(geneID = geneFC_df$Row.names, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", #Convert gene IDs into ENTREZ IDs
                     OrgDb = "org.Hs.eg.db")

geneList <- dplyr::distinct(entrezid_list,SYMBOL,.keep_all=TRUE)
colnames(geneList)[1]<-"symbol"
colnames(diffExp)[1]<-"symbol"
geneDF <- diffExp %>% 
  inner_join(geneList,by="symbol") #Add log 2 fold change values to the existing genes
geneDFS<-geneDF %>% arrange(desc(log2FoldChange))#Sort genes by log2 fold change values 
geneFC_list = geneDFS$log2FoldChange 
names(geneFC_list) <- geneDFS$ENTREZID

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

print("====== runing Gene Sets Enrichment Analysis ======")
gsea_res2 <- GSEA(geneFC_list,
                  TERM2GENE = m_t2g,
                  minGSSize = 10,
                  maxGSSize = 500,
                  p.adjustCutoff = 1,
                  pAdjustMethod = "BH"
)
emgt2<-data.frame(gsea_res2)
return(emgt2)
}

emgt2res<-bulkRNA_GSEA(diffExp)#Return a dataframe of GSEA enrichment results.

plot_GSEA_barplot <- function(gsea_result, top_up = 2, top_down_remove = c(1,9), 
                              up_color = '#A982B2', down_color = '#4466B0') {
  gsea_bar_pva <- gsea_result %>%
    filter(p.adjust < 0.05) %>%
    select(Description, NES, p.adjust)
  
  gsea_bar_pva$Description <- gsub("HALLMARK_", "", gsea_bar_pva$Description)
  gsea_bar_pva$Description <- gsub("_", " ", gsea_bar_pva$Description)
  
  gsea_bar_pva$type <- ifelse(gsea_bar_pva$NES > 0, "up", "down")
  
  NESup <- gsea_bar_pva[gsea_bar_pva$type == "up", ] %>%
    arrange(NES)
  NESup$orders <- paste0("c", 1:nrow(NESup))  
  
  NESdown <- gsea_bar_pva[gsea_bar_pva$type == "down", ] %>%
    arrange(desc(NES))
  if(length(top_down_remove) > 0){
    NESdown <- NESdown[-top_down_remove, ]
  }
  NESdown$orders <- paste0("a", 1:nrow(NESdown))
  
  gsea_bar_pva <- rbind(NESup, NESdown)
  
  labelBar <- str_to_title(c(rev(NESdown$Description), NESup$Description))

  GSEAbarplot <- ggplot(gsea_bar_pva, aes(x = orders, y = NES, fill = type)) +
    geom_col(width = 0.72) +
    coord_flip() +
    xlab(NULL) +
    scale_fill_manual(values = c('up' = up_color, 'down' = down_color)) +
    scale_x_discrete(labels = labelBar) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_bw() +
    geom_hline(yintercept = 0, color = "black") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(colour = "black", size = 10),
      axis.text.y = element_text(colour = "black", size = 10),
      axis.title.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(GSEAbarplot)
}

#-------------------------------
# Example usage:
GSEAbarplot <- plot_GSEA_barplot(emgt2res)

