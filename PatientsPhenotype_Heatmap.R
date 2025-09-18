#' @title plot_IRAK2_heatmap
#' @description Depict the distribution and frequency of different clinical manifestations in distinct patients
#' @param IRAK2 Input Phenotype status matrix of patients, where each column represents a patient and
#'  each row represents a phenotype; 1 indicates the presence of the phenotype, and 0 indicates its absence
#----------------------example as:
# sys	phenotype	P1	P2	P3	P4	P5	P6	P7	P8	P9	P10	P11	P12
# Immune system	Recurrent infections	1	1	1	0	0	1	1	0	1	0	0	0
# Immune system	Positive autoantibody	0	0	1	1	1	1	1	0	0	1	1	0
# Cutaneous-Mucosal system	Hair loss	0	0	1	1	0	1	0	0	0	0	0	0
# Cutaneous-Mucosal system	Facial erythema	0	0	1	0	0	1	0	0	0	0	0	0
# Cutaneous-Mucosal system	Oral ulcers	0	0	0	0	0	0	1	0	1	1	0	0
# Cutaneous-Mucosal system	Genital ulcers	0	0	0	0	0	0	0	0	1	0	0	0

#' @param count Input Matrix counting the number of patients exhibiting each phenotype
#----------------------example as:
# sum	phenotype	sys
# 6	Recurrent infections	Immune system
# 7	Positive autoantibody	Immune system
# 3	Hair loss	Cutaneous-Mucosal system
# 2	Facial erythema	Cutaneous-Mucosal system
# 3	Oral ulcers	Cutaneous-Mucosal system
# 1	Genital ulcers	Cutaneous-Mucosal system
# 5	Skin rash	Cutaneous-Mucosal system
#' return list contains heatmap and barplot


library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)

plot_IRAK2_heatmap <- function(IRAK2, count) {
  
  Pheno <- data.frame()
  for(i in 3:ncol(IRAK2)){
    data <- data.frame(
      sys = IRAK2$sys,
      phenotype = IRAK2$phenotype,
      value = IRAK2[,i]
    )
    data$Category <- rep(colnames(IRAK2)[i], nrow(IRAK2))
    Pheno <- rbind(Pheno, data)
  }
  
  sys_order <- count %>%group_by(sys) %>%summarise(sys_total = sum(sum)) %>%arrange(desc(sys_total))
  count$sys <- factor(count$sys, levels = sys_order$sys)
  count_sorted <- count %>%arrange(sys, desc(sum))
  sys_order <- levels(count_sorted$sys)
  new_sys_order <- c("Immune system","Hematological system", "Musculoskeletal system", "Cutaneous-Mucosal system",
                     "Urinary system","Digestive system", "Nervous system","Endocrine system")
  
  count_sorted$sys <- factor(count_sorted$sys, levels = new_sys_order)
  count_sorted <- count_sorted %>%arrange(sys, desc(sum))
  total_patients <- sum(count_sorted$sum)
  count_sorted$prop <- count_sorted$sum / total_patients
  count_sorted$phenotype <- factor(count_sorted$phenotype, levels = rev(count_sorted$phenotype))
  Pheno$Category<-factor(Pheno$Category,levels = rev(unique(Pheno$Category)))
  Pheno$phenotype <- factor(Pheno$phenotype, levels =(unique(count_sorted$phenotype)))
  
  statuHeatmap<-ggplot(Pheno, aes(x =phenotype, y = Category, fill = as.factor(value))) +
    geom_tile(color = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("0" = "grey90", "1" = "#7BBFC8")) +
    theme_minimal(base_size = 12) +
    xlab(NULL)+ylab(NULL)+
    scale_x_discrete(expand = expansion(mult = c(0,0.1)))+
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1,colour = "black"),#angle = 45, hjust = 1,
      axis.text.y = element_text(size = 10,colour = "black"),
      panel.grid = element_blank(),
      legend.position = "none"
    )


  count_phenotype<-ggplot(count_sorted, aes(x = phenotype, y = prop)) +
    geom_col(fill="grey90",width=0.67) +
    theme_minimal(base_size = 12) +
    theme_classic()+
    xlab("")+
    ylab(NULL)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x= element_blank(),
      plot.title = element_blank())
  count_phenotype
  
  final_plot <- wrap_plots(count_phenotype, statuHeatmap, ncol = 1) + plot_layout(guides = "collect", heights = c(1,3))
  return(list(heatmap = statuHeatmap, barplot = count_phenotype, final=final_plot))
}
