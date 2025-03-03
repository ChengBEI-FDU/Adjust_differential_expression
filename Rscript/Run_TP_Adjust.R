# Set PATH
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Read files
de_results <- read.csv("input_DE.csv", row.names = 1)
genes_tp <- read.csv("input_TP.csv", row.names = 1)


# parameter
loess_span <- 0.5
loess_degree <- 1


# Adjust
de_results_merged <- merge(de_results, 
                           data.frame(gene = genes_tp$gene, TP = genes_tp$TP), 
                           by = 'gene', all.x = T)

cat(paste0(length(de_results_merged$gene[is.na(de_results_merged$TP)]), " gene(s) -- ",
           paste0(de_results_merged$gene[is.na(de_results_merged$TP)], collapse = ", "), " --",
           " is (are) not included due to lack of TP value(s)."),
    file = "output_log.txt")

de_results_merged <- de_results_merged[!is.na(de_results_merged$TP),]

tmp_DE_res_up <- de_results_merged[de_results_merged$logFC > 0,]
tmp_DE_res_down <- de_results_merged[de_results_merged$logFC < 0,]

tmp_up_model <- loess(log10(logFC) ~ log10(TP), data = tmp_DE_res_up, span = loess_span, degree = loess_degree)
tmp_DE_res_up$logFC_residual <- 10**tmp_up_model$residuals
tmp_DE_res_up$logFC_loess_fit <- 10**tmp_up_model$fitted
tmp_DE_res_up$logFC_adjusted <- tmp_DE_res_up$logFC_residual + 10**mean(log10(tmp_DE_res_up$logFC_loess_fit))

tmp_down_model <- loess(log10(-logFC) ~ log10(TP), data = tmp_DE_res_down, span = loess_span, degree = loess_degree)
tmp_DE_res_down$logFC_residual <- 10**tmp_down_model$residuals
tmp_DE_res_down$logFC_loess_fit <- 10**tmp_down_model$fitted
tmp_DE_res_down$logFC_adjusted <- -(tmp_DE_res_down$logFC_residual + 10**mean(log10(tmp_DE_res_down$logFC_loess_fit)))

DE_res_adj <- rbind(tmp_DE_res_up, tmp_DE_res_down)


# Output
write.csv(DE_res_adj, "output_DE_adjusted.csv")




