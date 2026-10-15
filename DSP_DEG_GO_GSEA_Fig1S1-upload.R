##################################################
## Project: TNBC SLNs DSP T-cell zone: CD11c enriched region
## Script purpose:Fig.1c-d DEG-vocano,boxplot; FigS1f-g GO enrichment,GSEA analysis
## By ZQ
##################################################

rm(list = ls())
library(limma)
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
library(tidyr)
library(reshape2)

##==========================================================================
# SLN T cell zone: CD11c enriched region 
# GO Enrichment----------------
##==========================================================================

getwd()
setwd("E:/TNBC_SLN_Program/code_nCounter_DSP_mIF/DSP_code_76109/DSP_code/")

#Data prepare
CD11c_matrix <- read.csv("Docu/DSP_Tzone_CD11cEnrich_matrix.csv", row.names = 1,check.names = FALSE)
head(CD11c_matrix)
Group.Info <- read.csv("Docu/DSP_Tzone_Group.Info_CD11cEnrich.csv")
expr_samples <- colnames(CD11c_matrix)
group_ordered <- sapply(expr_samples, function(sample) {
  row_idx <- which(Group.Info$segement_ROI == sample)
  if(length(row_idx) > 0) {
    return(Group.Info$Group[row_idx])
  } else {return(NA) }
})

for(i in 1:length(expr_samples)) {
  cat(i, ". ", expr_samples[i], " - ", group_ordered[i], "\n")
}


##----------------------------------------------##
## limma DEG --------                 
##----------------------------------------------##
list <- factor(group_ordered, levels = c("pCR", "pNR"))
list <- model.matrix(~factor(list) + 0)
colnames(list) <- c("pCR", "pNR")

expr_log <- log2(as.matrix(CD11c_matrix) + 1)   

group <- factor(group_ordered, levels = c("pCR", "pNR"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_log, design)
contr <- makeContrasts(pCR - pNR, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contr), robust = TRUE, trend = TRUE)
diff_data <- topTable(fit2, number = Inf, adjust.method = "fdr")

##----------------------------------------------##
## Fig1.c CD11c enriched region DEG volcano pathway plot--------                 
##----------------------------------------------##
# Increase, Decrease Gene Label
label_data <- diff_data %>% 
  rownames_to_column("SYMBOL") %>%
  mutate(Changed = if_else(P.Value > 0.05, "N.S.", 
                           if_else(logFC >= 0.1, "Increased",
                                   if_else(logFC <= -0.1, "Decreased", "N.S.")))) %>%
  mutate(Changed = factor(Changed, levels = c("Increased", "Decreased", "N.S.")))
label_data <- label_data %>% mutate(Label = if_else(Changed != "N.S.", SYMBOL, ""))
write.csv(label_data,file="Results/DEG_20261012/DSP_CD11c_DEGlimma_Pvalue_IncDec.csv",row.names = F)

#Increase, Decrease gene pathway-------
tempup <- subset(label_data,Changed=="Increased",select=c(SYMBOL))
tempdown <- subset(label_data,Changed=="Decreased",select=c(SYMBOL))


pathways_list <- list(
  "Antigen process and present"=c("HLA-C","HLA-DRB3","HLA-DRB4","LAMP1","FLT3","UBB","CCL21"),
  "T cell priming/activation"=c("CD8A","CD8B","CD3G","IL2RB","NKG7","GNLY","PRF1"),
  "Metabolism"=c("C1QBP","PKM","NDUFB4","NDUFA12","NDUFB10","NDUFA2","NDUFA13","NDUFA11","NDUFA6","NDUFA1","NDUFA7","NDUFB8","NDUFB1"),
  "Immunosuppression"=c("ENTPD1","DUSP1","IL6","OSM"),
  "Stromal"=c("CD34","ITGA4","ITGB4","FZD1","COL1A2"))

data <- label_data
data$row <- data$SYMBOL
colnames(data)[colnames(data) == "logFC"] <- "log2FoldChange" 

# pathway inform
data$GO_term <- "others"
for (pathway_name in names(pathways_list)) {
  pathway_genes <- pathways_list[[pathway_name]]
  data$GO_term[data$row %in% pathway_genes] <- pathway_name
}
# pathway gene marked
data$label_path <- NA
pathway_genes_all <- unique(unlist(pathways_list))
data$label_path[data$row %in% pathway_genes_all] <- data$row[data$row %in% pathway_genes_all]

head(data)
## Pathway volcano plot--------                 
p <- ggplot(data, aes(x = log2FoldChange, y = -log10(P.Value))) +
  geom_point(data = subset(data, GO_term == "others"),
    aes(color = Changed),alpha = 0.6,size = 2.2,shape = 16)+
  scale_color_manual(name = "Expression",values = c("Increased" = "#FF69B4","Decreased" = "#87CEEB","N.S." = "grey70"))+
  # Marking pathway gene 
  geom_point(data = subset(data, GO_term != "others"),
    aes(fill = GO_term),size = 4,shape = 21,color = "black",stroke = 0.8,alpha = 0.9 )+
  # pathway gene color
  scale_fill_manual(
    name = "Pathway",
    values = c(
      "Antigen process and present" = "#FFC839",
      "T cell priming/activation"="#CE0000",
      "Metabolism" = "#DA70D6",                                          
      "Immunosuppression" = "#009966",
      "Stromal-related" = "#6699FF" ),
    breaks = c( 
      "Antigen process and present",
      "T cell priming/activation",
      "Metabolism",                                          
      "Immunosuppression",
      "Stromal-related"))+
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = 'black', lwd = 0.5) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  # Axis labels
  labs( x = TeX("$Log_2\\, Fold\\,Change$"),y = TeX("$-Log_{10}(P-value)$"),
    title = "SLN T-region DC enriched area pCR vs pNR") +
  # Label text
  geom_text_repel(aes(label = label_path),size = 3.2,box.padding = 0.5,point.padding = 0.3,
    segment.color = "grey40",segment.size = 0.3,segment.alpha = 0.7,min.segment.length = 0.1,
    max.overlaps = Inf,force = 2,max.time = 3,max.iter = 200000,direction = "both",    
    nudge_x = ifelse(data$log2FoldChange > 0, 0.1, -0.1),
    color = "#002060",fontface = "bold") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "grey95"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 1.5, color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18, margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 18),
    axis.title.y = element_text(margin = margin(r = 15)),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.text = element_text(face = "bold", color = "black", size = 14),
    axis.text.x = element_text(face = "bold", color = "black", size = 14),
    axis.text.y = element_text(face = "bold", color = "black", size = 14),
    legend.position = "bottom",legend.box = "vertical",legend.direction = "vertical",
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.key = element_rect(fill = "white", color = NA),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(face = "bold", size = 12),
    plot.margin = margin(15, 15, 15, 15)) +
  # Legend style
  guides(color = "none",
    fill = guide_legend(title = "Immune related functions",override.aes = list(size = 4.5, alpha = 1, stroke = 0.5),
      title.position = "top",title.hjust = 0.5,order = 1)) +
      xlim(c(-1.8, 1.8))+ ylim(c(0,7.5))

ggsave(plot = p,"Results/DEG_20261012/DSP_volcano_pCR.vs.pNR_limma_feature.pdf",width = 8,height = 10,dpi = 600,device = cairo_pdf)


##----------------------------------------------##
## Fig1.d Diff Gene Boxplot--------                 
##----------------------------------------------##
genes_to_plot <- c("HLA-C","HLA-DRB3","HLA-DRB4","LAMP1","UBB","FLT3","CCL21",
                   "CD8A","CD8B","CD3G","IL2RB","NKG7","GNLY","PRF1",
                   "C1QBP","PKM","NDUFB4","NDUFA12","NDUFB10","NDUFA2","NDUFA13","NDUFA11","NDUFA6","NDUFA1","NDUFA7","NDUFB1",
                   "ENTPD1","DUSP1","IL6","OSM",
                   "CD34","ITGA4","ITGB4","FZD1","COL1A2")

# 1. 检查基因是否存在
missing_genes <- setdiff(genes_to_plot, rownames(expr_log))
if (length(missing_genes) > 0) message("不在矩阵中的基因: ", paste(missing_genes, collapse = ", "))
genes_use <- intersect(genes_to_plot, rownames(expr_log))

# 2. 样本分组信息
group_df <- data.frame(
  SampleID = colnames(expr_log),
  Group    = factor(group_ordered, levels = c("pCR", "pNR")),
  Patient  = Group.Info$Sample[match(colnames(expr_log), Group.Info$segement_ROI)])

# 3. 转为长格式
expr_long <- reshape2::melt(expr_log[genes_use, , drop = FALSE],
                            varnames = c("Gene", "SampleID"), value.name = "Expression") %>%
  mutate(Gene = factor(as.character(Gene), levels = genes_use),
         SampleID = as.character(SampleID)) %>%
  left_join(group_df, by = "SampleID")

# 4. 用 limma 结果作为 P 值标注，与火山图保持一致
p_lab <- diff_data[genes_use, c("logFC", "P.Value", "adj.P.Val")] %>%
  rownames_to_column("Gene") %>%
  mutate(Gene = factor(Gene, levels = genes_use),
         label = sprintf("P = %.2g", P.Value)) %>%
  left_join(expr_long %>% group_by(Gene) %>%
              summarise(y = max(Expression) + 0.12 * diff(range(Expression))), by = "Gene")

# 5. 作图
p <- ggplot(expr_long, aes(x = Group, y = Expression)) +
  geom_boxplot(aes(color = Group), width = 0.6, linewidth = 0.6,
               alpha = 0.5, outlier.shape = NA) +
  geom_beeswarm(aes(color = Group), size = 1.2, cex = 4, priority = "ascending") +
  geom_text(data = p_lab, aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE, size = 3) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_color_manual("Group", values = brewer.pal(9, "Set1")[1:2]) +
  labs(x = "", y = TeX("$Log_2$ (Normalized Expression + 1)"),
       title = "Differentially expressed genes: pCR vs pNR") +
  theme_few() +
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 13),
        strip.text = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 15, hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        aspect.ratio = 1)

# 6. 保存
n_cols <- 5
base_size <- 1.2
n_rows <- ceiling(length(genes_use) / n_cols)
plot_width <- n_cols * base_size + 2
plot_height <- n_rows * base_size + 2

# Save
dir.create("./Boxplot", showWarnings = FALSE, recursive = TRUE)
filename <- "./Boxplot/DSP_diffgene_Boxplot.pdf"
ggsave(p, filename = filename, width = plot_width, height = plot_height)


##----------------------------------------------##
## Source data export --------
##----------------------------------------------##
out_dir <- "./Boxplot"

source_long <- expr_long %>%
  dplyr::select(Gene, SampleID, Patient, Group, Expression) %>%
  dplyr::arrange(Gene, Group, Patient, SampleID)

source_wide <- source_long %>%
  tidyr::pivot_wider(names_from = Gene, values_from = Expression) %>%
  dplyr::arrange(Group, Patient, SampleID)

stats_out <- expr_long %>%
  dplyr::group_by(Gene, Group) %>%
  dplyr::summarise(n = dplyr::n(), median = median(Expression),
                   Q1 = quantile(Expression, 0.25), Q3 = quantile(Expression, 0.75),
                   .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Group, values_from = c(n, median, Q1, Q3)) %>%
  dplyr::left_join(p_lab %>% dplyr::select(Gene, log2FC = logFC, P.Value, adj.P.Val), by = "Gene") %>%
  dplyr::mutate(Test = "limma moderated t-test (eBayes, robust = TRUE, trend = TRUE), two-sided",
                Expression_unit = "log2(Q3-normalized count + 1)") %>%
  dplyr::arrange(factor(Gene, levels = genes_use))


# 一个 Excel 文件、每个子图一个 sheet
# install.packages("openxlsx")
openxlsx::write.xlsx(
  list("Fig1d_points" = source_long,
       "Fig1d_wide"   = source_wide,
       "Fig1d_stats"  = stats_out),
  file = file.path(out_dir, "SourceData_Fig1d.xlsx"), overwrite = TRUE)

##==========================================================================
# GO Enrichment----------------
##==========================================================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(RColorBrewer)

head(data)
gene_up= subset(data,Changed =="Increased",select="SYMBOL")
gene_up=gene_up$SYMBOL
gene_up.df <- bitr(gene_up, fromType = "SYMBOL", 
                   toType = c("ENSEMBL","ENTREZID"), 
                   OrgDb = org.Hs.eg.db)
goUP <- enrichGO(gene_up.df$ENTREZID, OrgDb = "org.Hs.eg.db",ont = "ALL", 
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.2, 
                 readable = TRUE)
goUP_result <- goUP@result

gene_down= subset(data,Changed =="Decreased",select="SYMBOL")
gene_down=gene_down$SYMBOL
gene_down.df <- bitr(gene_down, fromType = "SYMBOL", 
                     toType = c("ENSEMBL","ENTREZID"), 
                     OrgDb = org.Hs.eg.db)
goDown <- enrichGO(gene_down.df$ENTREZID, OrgDb = "org.Hs.eg.db",ont = "ALL", 
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.2, 
                   readable = TRUE)
goDown_result <- goDown@result

##----------------------------------------------##
## FigS1.F GO enrichment analysis Lolipop graph --------                 
##----------------------------------------------##

# Data prepare
prepare_lollipop_data <- function(go_up, go_down, top_n = 15) {
  up_sig <- go_up %>%
    filter(p.adjust < 0.05) %>%
    arrange(p.adjust) %>%
    head(top_n) %>%
    mutate(Change = "pCR",type = 1)
 
  down_sig <- go_down %>%
    filter(p.adjust < 0.05) %>%
    arrange(p.adjust) %>%
    head(top_n) %>%
    mutate(Change = "pNR",type = -1)
  
  plot_data <- rbind(up_sig, down_sig)
  plot_data <- plot_data %>%
    mutate(
      neg_log_p = -log10(p.adjust),
      x = neg_log_p * type,
      margin = ifelse(type == 1, -0.2, 0.2), 
      hjust = ifelse(type == 1, 1, 0)) %>%  # Text alignment direction
    arrange(x) %>%
    mutate(Description = factor(Description, levels = Description))
  return(plot_data)
}
plot_data <- prepare_lollipop_data(goUP_result, goDown_result, top_n = 15)

# Color
two_colors <- c("pCR" = "#E31A1C", "pNR" = "#1F78B4")
my_theme <- theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    plot.margin = margin(t = 20, r = 40, b = 20, l = 40),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5))
  
# x-axis range,reserve space for the text
x_range <- range(plot_data$x)
x_limit <- c(floor(x_range[1] - 3), ceiling(x_range[2] + 3))

# Plot
p <- plot_data %>%
  ggplot(aes(x, Description, colour = Change)) +
  geom_segment(aes(xend = 0, yend = Description),linewidth = 1.2, show.legend = FALSE) + 
  geom_point(aes(size = Count), alpha = 0.8) +
  geom_text(aes(x = margin, label = Description, hjust = hjust),
            show.legend = FALSE, colour = "black", size = 4.5, fontface = "bold") +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  labs(
    x = '-log10(p.adjust)',y = NULL, 
    title = "GO Enrichment: pCR vs pNR",
    size = "Gene Count",
    colour = "Regulation") +
  scale_colour_manual(values = two_colors) +
  scale_size_continuous(range = c(3, 8)) +
  scale_x_continuous(labels = abs, limits = x_limit) +
  my_theme

##==========================================================================
# GO GSEA ----------
##==========================================================================
library(clusterProfiler)
library(msigdbr)
library(DOSE)

gene_symbols <- rownames(diff_data)
gene_entrez <- bitr(gene_symbols, 
                    fromType = "SYMBOL",
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)
diff_data$SYMBOL <- rownames(diff_data)
diff_merged <- merge(diff_data, gene_entrez, by = "SYMBOL", all.x = FALSE)

diff_merged <- diff_merged %>%group_by(ENTREZID) %>%slice_max(abs(t), n = 1) %>% ungroup()
gene_list <- diff_merged$t
names(gene_list) <- diff_merged$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)
head(gene_list, 20)

GOgsea<-gseGO(gene_list,
              OrgDb='org.Hs.eg.db',
              ont='ALL',
              pvalueCutoff=0.05,
              minGSSize=10,
              maxGSSize=500,
              keyType="ENTREZID")

##----------------------------------------------##
## FigS1.g GOgesa Comparison Plot --------                 
##----------------------------------------------##
# Pathways prepare
GOgsea_result <- GOgsea@result
significant_gsea <- GOgsea_result %>% filter(p.adjust < 0.05) %>% arrange(desc(abs(NES)))
upregulated_pathways <- significant_gsea %>% filter(NES > 0) %>% arrange(desc(NES))
downregulated_pathways <- significant_gsea %>% filter(NES < 0) %>% arrange(NES)
top_pathways <- rbind(head(upregulated_pathways, 10),head(downregulated_pathways, 10))
  
# Plot  
barplot_comparison <- ggplot(top_pathways, aes(x = reorder(Description, NES), y = NES)) +
  geom_col(aes(fill = ifelse(NES > 0, "pCR enriched", "pNR enriched")),alpha = 0.8) + 
  scale_fill_manual(values = c("pCR enriched" = "#E31A1C","pNR enriched" = "#1F78B4"),name = "Enrichment") + 
  coord_flip() +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  labs(x = "GO Terms",
       y = "Normalized Enrichment Score (NES)",
       title = "GSEA GO term top10 Enriched Pathways",
       subtitle = "pCR vs pNR") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),legend.position = "top")
ggsave("Results/DEG_20261012/FigS1g_GSEA_barplot_comparison.pdf", barplot_comparison, width = 8, height = 6.5)

