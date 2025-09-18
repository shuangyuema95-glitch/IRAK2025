#' @title HeatmapForBulk
#' @description Plot gene expression heatmap for specific singaling pathway
#' @param data Input gene expression matrix of bulk-RNA seq data 
#' @param signaling Input normalized gene expression matrix by function 'NorRNA',e.g. signaling="Type I interferon"
#' @return hplot Gene expression heatmap

library(pheatmap)
library(ggplotify)

bk = unique(c(seq(-1.5,1.5, length=100)))#Set the color gradient for the heatmap
HeatmapForBulk<-function(data1,signaling){
set.seed(123)
hplot<-pheatmap(data1,scale = "row",cluster_cols=F,cluster_rows = T,
                    show_rownames = T,border_color=NA,main = "Type I interferon (UNS level)",
                    show_colnames = T,
                    fontsize = 6.5,legend = T,treeheight_col = 0,treeheight_row = 0,
                    legend_breaks=seq(-2,2,0.5),
                    color=colorRampPalette(c("#6B76B1","#738AB7","#7998C2","#75B9C2","#88C1CF","#A7D4D9",
                                             "#B5DBE3","#D6EAF0","#FEF9FB","#FEF1E4","#FEEEDE","#F7DEC1",
                                             "#F4DBC3","#F7CB99","#F3BC7F","#EDAC68","#E9A35E","#EE935A",
                                             "#E77051"))(101),
                    annotation_legend = T, breaks = bk,
                    cellwidth = 8.8,cellheight = 7.5,
                    clustering_distance_rows="euclidean",
                    clustering_method = "average"
)
hplot<-as.ggplot(hplot)
return(hplot)
}



