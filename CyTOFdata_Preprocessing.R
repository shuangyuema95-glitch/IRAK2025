#' @title CyTOF Data Preprocessing and Clustering
#' @description Load FCS files from given directories, prepare data with CATALYST, perform clustering and t-SNE dimensionality reduction, and return plot objects.
#' @param LPS_dir input Directory containing LPS condition FCS files.
#' @param basal_dir input  Directory containing basal/control condition FCS files.
#' @param panel input Data frame specifying marker panel, typically loaded from CATALYST \code{panel}.
#' @param md_file input Path to metadata Rdata file (should contain \code{md1}).
#' @param exclude_samples Character vector of sample IDs to exclude, e.g., c("LL").
#' @param tsne_cells Number of cells to subsample for t-SNE, default 10000.
#' @param cluster_params List of clustering parameters: xdim, ydim, maxK, seed.
#' @param colors Vector of colors for plotting clusters.

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


run_CyTOF_analysis <- function(LPS_dir, basal_dir, panel, md_file,
                               exclude_samples = c("LL"),
                               tsne_cells = 1e4,
                               cluster_params = list(xdim=10, ydim=10, maxK=20, seed=1234),
                               colors = c("#4466B0","#74B47F","#F1D19F","#A982B2","#40ABE2",
                                          "#E2A9C9","#0F9296","#F1BAAB","#EC706E","#AF478A",
                                          "#84CAE8","#E1C767","#DBD8E7","#C1E1DA","#25853A")) {
  load(md_file)  
  if(!exists("md1")) stop("md1 object not found in md_file")
  

  LPS_files <- list.files(LPS_dir, pattern="*fcs", full.names = TRUE)
  basal_files <- list.files(basal_dir, pattern="*fcs", full.names = TRUE)
  
  if(length(exclude_samples) > 0) {
    LPS_files <- LPS_files[!grepl(paste(exclude_samples, collapse="|"), LPS_files)]
    basal_files <- basal_files[!grepl(paste(exclude_samples, collapse="|"), basal_files)]
  }
 
  LPS_flowset <- read.flowSet(files = LPS_files)
  basal_flowset <- read.flowSet(files = basal_files)
  all_flowset <- read.flowSet(files = c(LPS_files, basal_files))

  md_LPS <- md1 %>% dplyr::filter(condition == "LPS", !sample_id %in% exclude_samples)
  md_basal <- md1 %>% dplyr::filter(condition == "ctrl", !sample_id %in% exclude_samples)
  
  sce_LPS <- prepData(LPS_flowset, panel, md_LPS, features = panel$fcs_colname)
  sce_basal <- prepData(basal_flowset, panel, md_basal, features = panel$fcs_colname)
  sce_all <- prepData(all_flowset, panel, md1, features = panel$fcs_colname)
  
  
  sce_LPS <- cluster(sce_LPS, xdim=cluster_params$xdim, ydim=cluster_params$ydim,
                     maxK=cluster_params$maxK, seed=cluster_params$seed)
  sce_basal <- cluster(sce_basal, xdim=cluster_params$xdim, ydim=cluster_params$ydim,
                       maxK=cluster_params$maxK, seed=cluster_params$seed)
  sce_all <- cluster(sce_all, xdim=cluster_params$xdim, ydim=cluster_params$ydim,
                     maxK=cluster_params$maxK+3, seed=cluster_params$seed)
  
  set.seed(cluster_params$seed)
  LPS_tsne <- runDR(sce_LPS, "TSNE", cells = tsne_cells)
  set.seed(cluster_params$seed)
  basal_tsne <- runDR(sce_basal, "TSNE", cells = tsne_cells)
  set.seed(cluster_params$seed)
  all_tsne <- runDR(sce_all, "TSNE", cells = tsne_cells)
  
  
  LPS_map <- plotDR(LPS_tsne, "TSNE", color_by = "meta20", a_pal = colors) +
    theme_classic() +
    theme(legend.position="right", 
          panel.background = element_rect(color="black", size=0.45)) +
    xlab("tsne_1") + ylab("tsne_2")
  
  basal_map <- plotDR(basal_tsne, "TSNE", color_by = "meta20", a_pal = colors) +
    theme_classic() +
    theme(legend.position="right", 
          panel.background = element_rect(color="black", size=0.45)) +
    xlab("tsne_1") + ylab("tsne_2")
 
  return(list(
    sce_LPS = sce_LPS,
    sce_basal = sce_basal,
    sce_all = sce_all,
    LPS_tsne = LPS_tsne,
    basal_tsne = basal_tsne,
    all_tsne = all_tsne,
    LPS_map = LPS_map,
    basal_map = basal_map
  ))
}