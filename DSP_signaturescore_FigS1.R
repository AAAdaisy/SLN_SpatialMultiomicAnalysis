##################################################
## Project: TNBC SLNs DSP
## Script purpose: Signature score FigS1.e-g
## Date: 25-07-18
##################################################

rm(list = ls())
library(dplyr)
library(readxl)
library(tidyverse)
library(ggbeeswarm)
library(ggpubr)
library(latex2exp)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(RColorBrewer)
library(tibble)

##----------------------------------------------##
## Gene signature score calculate-------- 
##----------------------------------------------##
#data prepare
Siglist <-  read.csv("Results/sig_score/signature_list/DSP_signature_list.csv")
CD11c_matrix <- read.csv("Docu/Select_Counts/DSP_Tzone_CD11cEnrich_matrix.csv", row.names = 1,check.names = FALSE)
sample_names <- colnames(CD11c_matrix)
Group.Info <- read.csv("Docu/Select_Counts/DSP_Tzone_Group.Info_CD11cEnrich.csv")

# Calculate signature score：ⁿ√[(X₁+1) × (X₂+1) × ... × (Xₙ+1)]
calculate_signature_score <- function(genes, CD11c_matrix) {
  available_genes <- genes[genes %in% rownames(CD11c_matrix)]
gene_expr <- CD11c_matrix[available_genes, , drop = FALSE]
  signature_scores <- apply(gene_expr, 2, function(x) {
    geometric_mean <- exp(mean(log(x + 1)))
    return(geometric_mean)})
  return(signature_scores)
}

signature_types <- unique(Siglist$Signature)
results_list <- list()

for(sig_type in signature_types) {
  genes_in_signature <- Siglist$Gene[Siglist$Signature == sig_type]
  scores <- calculate_signature_score(genes_in_signature, CD11c_matrix)
  result_df <- data.frame(Sample = names(scores),Signature = sig_type,Score = scores,stringsAsFactors = FALSE)
  results_list[[sig_type]] <- result_df
}
all_scores <- do.call(rbind, results_list)

# Group Info.
all_scores$Group <- NA
for(i in 1:nrow(all_scores)) {
  sample_match <- Group.Info$segement_ROI == all_scores$Sample[i]
  if(any(sample_match)) {
    all_scores$Group[i] <- Group.Info$Group[sample_match][1]
  }
}
all_scores <- all_scores[!is.na(all_scores$Group), ]


##----------------------------------------------##
## Signature score boxplots   --------                 
##----------------------------------------------##
# Folder
output_dir <- "Results/sig_score/Signature_Box"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

plot_list <- list() 
for(sig_type in signature_types) {
  current_data <- all_scores[all_scores$Signature == sig_type, ]
  if(nrow(current_data) > 0 && length(unique(current_data$Group)) == 2) {
    # statistical test
    wilcox_result <- wilcox.test(Score ~ Group, data = current_data)
    p_value_wilcox <- wilcox_result$p.value
    # P-value label
    if(p_value_wilcox < 0.0001) {
      p_label <- sprintf("p = %.2e", p_value_wilcox)  # scientific notation
    } else if(p_value_wilcox < 0.001) {
      p_label <- sprintf("p = %.4f", p_value_wilcox)  # Keep 4 decimal places
    } else {
      p_label <- sprintf("p = %.3f", p_value_wilcox)
    }
    # y-axis range
    y_max <- max(current_data$Score, na.rm = TRUE)
    y_min <- min(current_data$Score, na.rm = TRUE)
    y_range <- y_max - y_min
    p_y_position <- y_max + y_range * 0.15  
    # Boxplot
    p <- ggplot(current_data, aes(x = Group, y = Score, fill = Group)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
      geom_jitter(aes(color = Group), width = 0.2, alpha = 0.9, size = 2) +
      labs(title = sig_type,x = "",y = "Signature Score") +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        axis.title.y = element_text(size = 11),
        axis.text = element_text(size = 11),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.line = element_line(linewidth = 0.7),
        axis.ticks = element_line(linewidth = 0.5),
        legend.position = "none") +
      scale_fill_manual(values = c("pCR" = "#E31A1C", "pNR" = "#1F78B4")) +
      scale_color_manual(values = c("pCR" = "#E31A1C", "pNR" = "#1F78B4")) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      annotate("text", x = 1.5, y = Inf, label = p_label,size = 4, hjust = 0.5, vjust = 1.5, fontface = "bold") 
    plot_list[[sig_type]] <- p  
    
  } else {cat("Warning: Insufficient data for", sig_type, "\n")}
}

# Save picture 
for(i in 1:length(plot_list)) {
  sig_name <- names(plot_list)[i]
  clean_name <- gsub("[^A-Za-z0-9_]", "_", sig_name)# Remove special characters
  pdf_filename <- file.path(output_dir, paste0("Signature_", clean_name, "_Wilcox_boxplot.pdf"))
  ggsave(pdf_filename, plot_list[[i]],width = 3.2, height = 4, bg = "white")
}




