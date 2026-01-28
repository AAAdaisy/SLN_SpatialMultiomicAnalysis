##################################################
## Project: TNBC SLNs nCounter Fig1
## Script purpose: nCounter analysis suppl fig.1
##################################################

rm(list = ls())

library(readxl)
library(tidyverse)
library(latex2exp)
library(limma)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(viridis)
library(RColorBrewer)
library(tibble)
library(reshape2)
library(ggbeeswarm)
library(ggpubr)


##----------------------------------------------##
## nCounter data --------                 
##----------------------------------------------##
Group.Info <- read.table("./docu/nCounter_group_infornation.txt", header = TRUE, stringsAsFactors = FALSE) 
Expression.Matrix <- read.table("./docu/nCounter_normalized_matrix.txt",header = TRUE,row.names = 1,sep = "\t",check.names = FALSE)   
                                

##----------------------------------------------##
## LIMMA DEG --------                 
##----------------------------------------------##
group_factor <- factor(Group.Info$Group, levels = c("pCR", "pNR"))
design <- model.matrix(~0 + group_factor)
colnames(design) <- c("pCR", "pNR")
fit <- lmFit(select_counts, design)
contrast.matrix <- makeContrasts(pCR - pNR, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, robust = TRUE)
DEG_results <- topTable(fit2, n = Inf, adjust = "fdr")

significant_genes <- DEG_results %>%
  rownames_to_column("SYMBOL") %>%
  filter(P.Value < 0.05 & abs(logFC) >= 0.2) %>%
  arrange(P.Value)
write.csv(significant_genes, paste0(docu,"nCounter_Significant_DEGs.csv"), row.names = FALSE)


##----------------------------------------------##
## Boxplot Data prepare --------                 
##----------------------------------------------##
# Genemane list
sig_gene_names <- significant_genes$SYMBOL
group_info_boxplot <- Group.Info
group_info_boxplot$SampleID <- group_info_boxplot$nCounterID
group_info_boxplot <- group_info_boxplot[, c("SampleID", "Group")]
# Check sample matching
sample_names <- colnames(select_counts)
matched_samples <- intersect(sample_names, group_info_boxplot$SampleID)

prepare_expression_data <- function(genes_to_plot, expr_matrix, group_info, matched_samples) {
  expr_subset <- expr_matrix[rownames(expr_matrix) %in% genes_to_plot, matched_samples, drop = FALSE]
  expr_subset$Gene <- rownames(expr_subset)
  expr_subset <- expr_subset[, c("Gene", matched_samples)]
  expr_long <- melt(expr_subset, id.vars = "Gene", variable.name = "SampleID", value.name = "Expression")
  expr_long <- merge(expr_long, group_info, by = "SampleID")
  return(expr_long)
}


##----------------------------------------------##
## Boxplots in pages (25 genes per page, 5×5 layout)--------                 
##----------------------------------------------##
dir.create("Paged_Boxplots", showWarnings = FALSE, recursive = TRUE)
# Layout
genes_per_page <- 25  
total_pages <- ceiling(length(sig_gene_names) / genes_per_page)

for (i in 1:total_pages) {
  start_idx <- (i - 1) * genes_per_page + 1
  end_idx <- min(i * genes_per_page, length(sig_gene_names))
  gene_subset <- sig_gene_names[start_idx:end_idx]
  
  cat(paste("Processing", i, "/", total_pages, "pages,Genes number:", length(gene_subset), "\n"))
  expr_data <- prepare_expression_data(gene_subset, select_counts, 
                                       group_info_boxplot, matched_samples)
  # Plot
  p <- ggplot(data = expr_data, aes(x = Group, y = Expression)) +
    geom_boxplot(aes(color = Group), width = 0.6, size = 0.6, 
                 position = position_dodge(0.6), alpha = 0.5) +
    geom_beeswarm(aes(color = Group, fill = Group), size = 1.2, cex = 4, 
                  dodge.width = 0.8, priority = "ascending") +
    facet_wrap(~Gene, scales = "free_y", ncol = 5) +
    labs(x = "", y = "Normalized Expression", 
         title = paste("Differential Genes - Page", i, "of", total_pages)) +
    theme_few() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 6),
      axis.title = element_text(size = 13),
      strip.text = element_text(size = 13),
      plot.title = element_text(size = 15, hjust = 0.5),
      panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
      aspect.ratio = 1  ) + # small image becomes a square
    scale_color_manual("Group", values = brewer.pal(9, "Set1")[c(1:2)]) +
    scale_fill_manual("Group", values = brewer.pal(9, "Set1")[c(1:2)]) +
    stat_compare_means(method = "wilcox.test", label = "p.format", 
                       label.x = 1.5, paired = FALSE, size = 4)
  
  # Adjust image size dynamically based on the number of genes
  n_genes <- length(gene_subset)
  n_cols <- 5 
  n_rows <- ceiling(n_genes / n_cols)
  base_size <- 1.2  # side length of each small figure
  plot_width <- n_cols * base_size + 2
  plot_height <- n_rows * base_size + 3  
  
  filename_pdf <- paste0("Paged_Boxplots/Page_", sprintf("%02d", i), "_DifferentialGenes.pdf")
  ggsave(p, filename = filename_pdf, width = plot_width, height = plot_height, dpi = 300)

}


##----------------------------------------------##
## Combine Boxplot (Share Y-axis) --------                 
##----------------------------------------------##
genes_combined <- c("VCAM1", "PECAM1", "ICAM1", "ICAM2")
available_combined_genes <- genes_combined[genes_combined %in% rownames(Expression.Matrix)]
expr_combined <- Expression.Matrix[available_combined_genes, common_samples, drop = FALSE]
expr_combined$Gene <- rownames(expr_combined)
expr_combined <- expr_combined[, c("Gene", common_samples)]

expr_long_combined <- melt(expr_combined, id.vars = "Gene", 
                           variable.name = "SampleID", value.name = "NormalizedCount")
expr_long_combined <- merge(expr_long_combined, group_info_boxplot, by = "SampleID")
colnames(expr_long_combined)[colnames(expr_long_combined) == "Gene"] <- "SYMBOL"

# Calculate the Y-axis range
y_max <- max(expr_long_combined$NormalizedCount, na.rm = TRUE)
y_min <- min(expr_long_combined$NormalizedCount, na.rm = TRUE)
y_range <- y_max - y_min
y_limit_max <- y_max + y_range * 0.3  # reserve space for the P value

# Boxplot visualization
p_combined <- ggplot(data = expr_long_combined,
                     aes(x = SYMBOL, y = NormalizedCount)) +
  geom_boxplot(aes(color = Group), width = 0.5, size = 0.4, 
               position = position_dodge(0.6), alpha = 0.5) +
  geom_beeswarm(aes(color = Group, fill = Group), size = 1.2, 
                dodge.width = 0.6, priority = "ascending") +
  labs(x = "", y = "Normalized counts") +
  theme_few() +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    panel.border = element_rect(fill = NA, color = "black", size = 1.2, linetype = "solid")
  ) +
  scale_color_manual("", values = brewer.pal(9, "Set1")[c(1:2)]) +
  scale_fill_manual("", values = brewer.pal(9, "Set1")[c(1:2)]) +
  scale_y_continuous(limits = c(y_min, y_limit_max), expand = c(0, 0)) +
  stat_compare_means(aes(group = Group), method = "wilcox.test", 
                     label = "p.format", size = 4)

dir.create("Combined_Genes_Boxplots", showWarnings = FALSE, recursive = TRUE)
ggsave(p_combined, 
       filename = "Combined_Genes_Boxplots/stromal_cell.pdf", 
       width = 4.2, height = 2.5, dpi = 300)

  




