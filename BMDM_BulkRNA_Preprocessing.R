#' @title NorRNA
#' @description Preprocess gene expression matrix by normalization
#' @title remove_Low
#' @description Filter gene with low expressions
#' @param data Input gene expression matrix of bulk-RNA seq data 
#' @param Normdata Input normalized gene expression matrix by function 'NorRNA'
#' @return resdata Normalized gene expression matrix
#' @return The matrix after removing lowly expressed genes

library(DESeq2)

NorRNA<-function(data){
  database=data[,c(1:ncol(data))];dim(data)
  head(database)
  database=database[complete.cases(database),];dim(database)
  
  cs1<-length(grep("12|15|17|56|60",colnames(database)));cs1 #The label of Irak2 ko mice 
  cs2<-(ncol(database)-cs1)
  condition<-factor(c(rep("WT",cs2),rep("KO",cs1)))
   dds <- DESeqDataSetFromMatrix(database, DataFrame(condition), design= ~ condition)
  dds <- DESeq(dds) 
  res <- results(dds)
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, 
                                                            normalized=TRUE)),by="row.names",sort=FALSE)
  return(resdata)
}

remove_Low<-function(Normdata){
  Normdata$fr<-apply(Normdata,1,function(x){
    y<-as.numeric(x);y<-y[8:ncol(Normdata)]
    length(which(y==0))/length(y)
  })
  filter(Normdata,fr<0.75)->Normdata#Remove genes that are expressed as 0 in more than 75% of samples
  Normdata$fr<-NULL
  return(Normdata)
}
