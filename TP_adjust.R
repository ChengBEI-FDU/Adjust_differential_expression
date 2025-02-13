library(openxlsx)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(pheatmap)
library(rtracklayer)
library(clusterProfiler)
library(reshape2)
library(Rmisc)
library(stringr)
library(infotheo)


#### Set PATH ####

setwd("/Users/beicheng/Desktop/TP_and_Benchmark/")
setwd("/Users/beicheng/Desktop/TP_and_Benchmark/Submit_version/github_code/")
setwd("") # CHANGE HERE



#### Import function for adjusting logFC by TP ####

LoessAdjust <- function(
  de_results, # DE results that contains columns 'gene' and 'logFC'
  genes_tp, # TP data that contains columns 'gene' and 'TP'
  loess_span = 0.5,
  loess_degree = 1
  ) {
  
  # merge DE results and TP list
  de_results_merged <- merge(de_results, 
                             data.frame(gene = genes_tp$gene, TP = genes_tp$TP), 
                             by = 'gene', all.x = T)
  
  # remove genes without TP annotation
  print(paste0(length(de_results_merged$gene[is.na(de_results_merged$TP)]), " gene(s) -- ",
               str_c(de_results_merged$gene[is.na(de_results_merged$TP)], collapse = ", "), " --",
               " is (are) not included due to lack of TP value(s)."))
  de_results_merged <- de_results_merged[!is.na(de_results_merged$TP),]
 
  # up- and down-regulated DEGs should be adjusted separately
  tmp_DE_res_up <- de_results_merged[de_results_merged$logFC > 0,]
  tmp_DE_res_down <- de_results_merged[de_results_merged$logFC < 0,]
  
  # loess fit
  tmp_up_model <- loess(log10(logFC) ~ log10(TP), data = tmp_DE_res_up, span = loess_span, degree = loess_degree)
  tmp_DE_res_up$logFC_residual <- 10**tmp_up_model$residuals
  tmp_DE_res_up$logFC_loess_fit <- 10**tmp_up_model$fitted
  tmp_DE_res_up$logFC_adjusted <- tmp_DE_res_up$logFC_residual + 10**mean(log10(tmp_DE_res_up$logFC_loess_fit))
  
  tmp_down_model <- loess(log10(-logFC) ~ log10(TP), data = tmp_DE_res_down, span = loess_span, degree = loess_degree)
  tmp_DE_res_down$logFC_residual <- 10**tmp_down_model$residuals
  tmp_DE_res_down$logFC_loess_fit <- 10**tmp_down_model$fitted
  tmp_DE_res_down$logFC_adjusted <- -(tmp_DE_res_down$logFC_residual + 10**mean(log10(tmp_DE_res_down$logFC_loess_fit)))
  
  DE_res_adj <- rbind(tmp_DE_res_up, tmp_DE_res_down)
  
  return(DE_res_adj)
}



#### Adjust logFC using TP ####

# 1) prepare two data frames: 
  # a) "de_results" contains the columns "gene and logFC"
  # b) "genes_tp" contains the columns "gene" and "TP"

de_results <- read.csv("./example/DE_index_52.csv", row.names = 1) # PRJNA733783, cholesterol (SRR14690793-95) vs control (SRR14690796-98)
genes_tp <- read.csv("TP_data/TP_mtb.csv", row.names = 1)

# 2) adjust
de_adjusted <- LoessAdjust(de_results = de_results, genes_tp = genes_tp, loess_span = 0.5, loess_degree = 1)

# 3) TP and logFC correlation
ggarrange(
  ggplot(de_adjusted, aes(x = TP, y = abs(logFC))) +
    geom_point(size = 0.5, color = 'grey') +
    geom_smooth(method = 'lm') +
    stat_cor(method = 'spearman') +
    scale_x_log10(breaks = c(0.5,1,2,4)) +
    scale_y_log10("Absolute logFC", breaks = c(0.25,0.5,1,2,4,8)) +
    theme_classic() + ggtitle('Before adjustment'),
  ggplot(de_adjusted, aes(x = TP, y = abs(logFC_adjusted))) +
    geom_point(size = 0.5, color = 'grey') +
    geom_smooth(method = 'lm') +
    stat_cor(method = 'spearman') +
    scale_x_log10(breaks = c(0.5,1,2,4)) +
    scale_y_log10("Absolute adj-logFC", breaks = c(0.25,0.5,1,2,4,8)) +
    theme_classic() + ggtitle('After adjustment')
)
dev.off()


