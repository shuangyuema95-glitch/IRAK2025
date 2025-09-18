#' @title NormCorrectBatch
#' @description Preprocess gene expression matrix by normalization and batch correction
#' @param fold input Path to store the gene expression matrix
#' @param data input Gene expression matrix of bulk-RNA seq data 
#' @return adjustedMat gene expression matrix


library(DESeq2)
library(dplyr)
library(sva)
library(tidyverse)

NormCorrectBatch<-function(fold,data){
  file<-paste0(fold,data)
  file<-as.data.frame(fread(file,sep="\t",fill=T,header=T,encoding = "UTF-8"))
  rownames(file)<-file$Geneid;file<-file[,-1]
  file=file[complete.cases(file),];dim(file) 
  print(colnames(file))
  print(rownames(file))
  
  cs1<-length(grep("C",colnames(file)));cs1
  cs2<-(ncol(file)-cs1);cs2
  
  patient_info<-factor(c(rep("ctrl",cs1),rep("patient",cs2)))
  dds <- DESeq(DESeqDataSetFromMatrix(file, DataFrame(patient_info), design= ~ patient_info))
  res <- results(dds)
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
  colnames(resdata)
  rownames(resdata)<-resdata$Row.names
  
  dim(resdata);index1<-grep("_",colnames(resdata))
  resdata$low<-apply(resdata,1,function(x){
    length(which(as.numeric(x[index1])==0))
  });filter(resdata,low<=ceiling(length(index1)*0.75))->resdata;resdata$low<-NULL#Remove genes that are expressed as 0 in more than 75% of samples
  dim(resdata)
  resdata<-resdata[,index1]
  print(rownames(resdata))
  
  if(length(grep("JAK",colnames(resdata)))!=0){#if input data is treatment-related gene expression matrix
    meta_data <- data.frame(
      sample_id = colnames(resdata),
      patient_info =factor(c(rep("Healthy Ctrl",cs1),rep("Patient",cs2))),
      batch_info=c("202304","202304","202311","202311","202411","202411",
                   "202304","202304","202311","202311","202411","202411","202411","202411"))
    rownames(meta_data)<-meta_data$sample_id
    print(meta_data)
    print(meta_data$batch_info)
    model <- model.matrix(~patient_info, data = meta_data)
    adjustedMat<- ComBat(dat = resdata, batch = as.character(meta_data$batch_info), mod = model)
  }else{#if input is not treatment-related gene expression matrix
    meta_data <- data.frame(  
      sample_id = colnames(resdata),
      patient_info =factor(c(rep("Healthy Ctrl",cs1),rep("Patient",cs2))),
      batch_info=rep(c("202304","202304","202311","202311","202411","202411"),2))
    rownames(meta_data)<-meta_data$sample_id
    print(meta_data)
    print(meta_data$batch_info)
    model <- model.matrix(~patient_info, data = meta_data)
    adjustedMat<- ComBat(dat = resdata, batch = as.character(meta_data$batch_info), mod = model)
    
  }
  
  dim(adjustedMat)
  return(adjustedMat)
}








