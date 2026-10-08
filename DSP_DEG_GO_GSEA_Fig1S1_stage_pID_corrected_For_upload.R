##################################################
## Project: TNBC SLNs DSP T-cell zone: CD11c enriched region
## Script purpose:Fig.1c-d DEG-vocano,boxplot; FigS1c-d GO enrichment,GSEA analysis
## Add cT_cm and Treatment_Class as covariate，using duplicateCorrelation to remove pseudo replication。
## Date: 260522
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
Group.Info <- read.csv("Docu/DSP_Tzone_Group.Info_CD11cEnrich.csv")

df_clinical_2 <- read_xlsx("E:/TNBC_SLN_Program/code_nCounter_DSP_mIF/nCounter_code/nCounter_code/docu/TNBC TDLN cohort 20221003.xlsx")
df_clinical <- read_excel("F:/0_zq/patient_information.xlsx")


library(stringr)

expr_samples <- colnames(CD11c_matrix)

# 基础数据：包含 segement_ROI，并关联 Group.Info 获取 Sample 编号
meta_data <- data.frame(segement_ROI = expr_samples) %>%
  left_join(Group.Info, by = "segement_ROI")

# 通过 Group.Info 的 Sample 匹配 df_clinical_2 的 Sentinel LN Path ID
meta_data <- meta_data %>%
  left_join(df_clinical_2 %>% 
              dplyr::select(`Sentinel LN Path ID`, PatientID_clean = `Patient ID`), 
            by = c("Sample" = "Sentinel LN Path ID"))

# 关联 df_clinical 获取核心临床变量 (patientID, 肿瘤分期，治疗方案)
meta_data <- meta_data %>%
  left_join(df_clinical %>% 
              dplyr::select(PatientID, `cT(stage)`), 
            by = c("PatientID_clean" = "PatientID"))

# 清洗、降维化疗方案并准备模型变量
meta_data_clean <- meta_data %>%
  # 剔除缺失临床信息的无效 ROI
  filter(!is.na(Group) & !is.na(`cT(stage)`) ) %>%
  mutate(
    cT_stage = factor(`cT(stage)`),
    
    Group = factor(Group, levels = c("pCR", "pNR")),
    
    Patient_Block = as.factor(PatientID_clean) 
  )

# 查看一下数据准备情况 (可用于核对最终纳入模型的 ROI 数量和患者数)
print("--- DSP 数据用于 Limma 混合模型的变量分布 ---")
print(table(meta_data_clean$Group))

#pCR pNR 
#8   8 

print(table(meta_data_clean$cT_stage))
#T1 T2 T3 
#10  4  2 
cat("纳入模型的独立患者数量: ", length(unique(meta_data_clean$Patient_Block)), "\n")
# 纳入模型的独立患者数量:  6 

# 确保矩阵的列完全对应清洗好的 meta_data_clean
expr_matrix_clean <- CD11c_matrix[, meta_data_clean$segement_ROI]

# 定义 Limma 的变量
group_factor <- meta_data_clean$Group
cT_stage        <- meta_data_clean$cT_stage
#treat_class  <- meta_data_clean$Treatment_Class
patient_id   <- meta_data_clean$Patient_Block

# 构建设计矩阵 (包含固定效应：分组，肿瘤大小/stage)
design <- model.matrix(~ 0 + group_factor + cT_stage)
colnames(design)[1:2] <- c("pCR", "pNR")

# 计算患者/切片内部的 ROI 间相关性 (解决多个 ROI 的 Pseudo-replication)
corfit <- duplicateCorrelation(expr_matrix_clean, design, block = patient_id)

# 拟合带有混合效应的线性模型 (传入 correlation = corfit$consensus)
df.fit <- lmFit(expr_matrix_clean, design, block = patient_id, correlation = corfit$consensus)

# 构造对比：pCR vs pNR (此时比较结果已完全剔除肿瘤大小/stage和化疗的混杂)
df.matrix <- makeContrasts(pCR - pNR, levels = design)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit, robust = TRUE)

# 提取校正后的差异基因
diff_data_corrected <- topTable(fit, n = Inf, adjust = "fdr")

# 获取显著基因标记 (用于火山图绘制)
label_data_corrected <- diff_data_corrected %>%
  rownames_to_column("SYMBOL") %>%
  mutate(Changed = if_else(P.Value > 0.05, "N.S.", 
                           if_else(logFC >= 1, "Increased",
                                   if_else(logFC <= -1, "Decreased", "N.S.")))) %>%
  mutate(Changed = factor(Changed, levels = c("Increased", "Decreased", "N.S."))) %>%
  mutate(Label = if_else(Changed != "N.S.", SYMBOL, ""))
label_data_ctsage_pID_corrected <- label_data_corrected
# 写入文件，直接对接下游火山图代码
write.csv(label_data_corrected, file="./Results/DSP_CD11c_DEGlimma_Pvalue_IncDec_Corrected_cTstage_patientID.csv", row.names = FALSE)

library(ggplot2)
library(ggrepel)
library(latex2exp)
library(scales)

label_data <- label_data_ctsage_pID_corrected

pathways_list <- list(
  "Antigen process and present"=c("CLU","HLA-DRB3","HLA-DRB4","LAMP1","CTSW","PSMB7","CCL21","CD5","NLRC5"),
  "T cell priming/activation"=c("CD8A","CD8B","CD3E","CD3D","CD2","IL2RB","NKG7","GNLY","PRF1","IL32"),
  "Metabolism"=c("C1QBP","PKM","NDUFB4","NDUFA12","NDUFB10","NDUFA2","NDUFA3","NDUFA2","NDUFA6","NDUFA1","NDUFA7","NDUFB8","NDUFB1","LDHB","SDHA"),
  "Immunosuppression"=c("ENTPD1","DUSP1","IL6","OSM","LILRB1","FCGR1A"),
  "Stromal"=c("CD34","ITGA4","ITGB4","VCAN","FZD1","COL2A1","LOXL2","LAMB3")
)

data <- label_data
data$row <- data$SYMBOL
colnames(data)[colnames(data) == "logFC"] <- "log2FoldChange" 

data$GO_term <- "others"
for (pathway_name in names(pathways_list)) {
  pathway_genes <- pathways_list[[pathway_name]]
  data$GO_term[data$row %in% pathway_genes] <- pathway_name
}

data$label_path <- NA
pathway_genes_all <- unique(unlist(pathways_list))
data$label_path[data$row %in% pathway_genes_all] <- data$row[data$row %in% pathway_genes_all]

data_all<- data

# 1. 过滤：确保只有在通路列表中，且差异显著（P < 0.05）的基因才被标记和染色
# 如果通路基因不显著，将其 GO_term 强制归为 "others"
data <- data %>%
  mutate(GO_term = ifelse(P.Value < 0.05 & GO_term != "others", GO_term, "others"))

# 2. 重新更新 label_path：仅针对显著且在通路内的基因生成标签
data$label_path <- NA
# 找出所有显著且在通路内的基因
significant_pathway_genes <- data$row[data$GO_term != "others"]
data$label_path[data$row %in% significant_pathway_genes] <- data$row[data$row %in% significant_pathway_genes]

# 3. 后续绘图逻辑保持不变
# 此时 ggplot 绘图时，那些非显著的通路基因会自动被归入 "others" 类别，
# 使用灰色点展示，并且不会被添加 label。



p <- ggplot(data, aes(x = log2FoldChange, y = -log10(P.Value))) +
  geom_point(data = subset(data, GO_term == "others"),
             aes(color = Changed), alpha = 0.6, size = 2.2, shape = 16) +
  scale_color_manual(name = "Expression", 
                     values = c("Increased" = "#FF69B4", "Decreased" = "#87CEEB", "N.S." = "grey70")) +
  
  geom_point(data = subset(data, GO_term != "others"),
             aes(fill = GO_term), size = 4, shape = 21, color = "black", stroke = 0.8, alpha = 0.9 ) +
  
  scale_fill_manual(
    name = "Pathway",
    values = c(
      "Antigen process and present" = "#FFC839",
      "T cell priming/activation" = "#CE0000",
      "Metabolism" = "#DA70D6",                                          
      "Immunosuppression" = "#009966",
      "Stromal" = "#6699FF" ),
    breaks = c( 
      "Antigen process and present",
      "T cell priming/activation",
      "Metabolism",                                          
      "Immunosuppression",
      "Stromal")) +
  
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = 'black', lwd = 0.5) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  
  labs( x = TeX("$Log_2\\, Fold\\,Change$"),y = TeX("$-Log_{10}(P-value)$"),
        title = "DC enriched area pCR vs pNR") +
  
  geom_text_repel(aes(label = label_path), size = 3.2, box.padding = 0.5, point.padding = 0.3,
                  segment.color = "grey40", segment.size = 0.3, segment.alpha = 0.7, min.segment.length = 0.1,
                  max.overlaps = Inf, force = 2, max.time = 3, max.iter = 200000, direction = "both",    
                  nudge_x = ifelse(data$log2FoldChange > 0, 0.1, -0.1),
                  color = "#002060", fontface = "bold") +
  
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
    legend.position = "bottom", legend.box = "vertical", legend.direction = "vertical",
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.key = element_rect(fill = "white", color = NA),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(face = "bold", size = 12),
    plot.margin = margin(15, 15, 15, 15)) +
  
  guides(color = "none",
         fill = guide_legend(title = "Immune related functions", override.aes = list(size = 4.5, alpha = 1, stroke = 0.5),
                             title.position = "top", title.hjust = 0.5, order = 1)) +
  
  scale_x_continuous(
    limits = c(-45, 45), 
    oob = scales::squish, 
    breaks = c( -80, -60, -40, -20, 0, 20, 40, 60, 80),
    labels = c("-80","-60","-40","-20", "0", "20", "40", "60", "80") # 坐标轴加上 > 和 < 符号
  ) +
  scale_y_continuous(
    limits = c(0, 7), # 稍微调高 Y 轴上限，容纳极显著基因
    oob = scales::squish,
    breaks = c(0, 2, 4, 6),
    labels = c("0", "2", "4", ">6")
  )

# 3. 保存图片
ggsave(plot = p, "Results/DEG/DSP_volcano_pCR.vs.pNR_limma_ctsage_pID_corrected.pdf", width = 6, height = 8, dpi = 600, device = cairo_pdf)

## Fig1.d Diff Gene Boxplot--------                 
##----------------------------------------------##
#genes_to_plot <- c("CLU","HLA-DRB3","HLA-DRB4","LAMP1","CTSW","PSMB7","CCL21",
#                   "CD8A","CD8B","CD3G","IL2RB","NKG7","GNLY","PRF1",
#                   "C1QBP","PKM","NDUFB4","NDUFA12","NDUFB10","NDUFA2","NDUFA13","NDUFA11","NDUFA6","NDUFA1","NDUFA7","NDUFB8","NDUFB1",
#                   "ENTPD1","DUSP1","IL6","OSM",
#                   "CD34","ITGA4","ITGB4","VCAN","FZD1")
genes_to_plot <- c("CLU","HLA-DRB3","HLA-DRB4","LAMP1","CTSW","PSMB7","CCL21","CD5","NLRC5",
                   "CD8A","CD8B","CD3E","CD3D","CD2","IL2RB","NKG7","GNLY","PRF1","IL32",
                   "NDUFB4","NDUFA12","NDUFB10","NDUFA2","NDUFA3","NDUFA2","LDHB","SDHA",
                   "ENTPD1","DUSP1","OSM","LILRB1","FCGR1A",
                   "VCAN","FZD1","COL2A1","LOXL2","LAMB3")


# Data prepare
expr_matrix <- expr_matrix_clean
#expr_subset <- expr_matrix[genes_to_plot, matched_samples, drop = FALSE]
#expr_subset$Gene <- rownames(expr_subset)
#expr_subset <- expr_subset[, c("Gene", matched_samples)]
#expr_long <- melt(expr_subset, id.vars = "Gene", variable.name = "SampleID", value.name = "Expression")
#expr_long <- merge(expr_long, group_info, by = "SampleID")


library(reshape2)
available_genes <- genes_to_plot[genes_to_plot %in% rownames(expr_matrix)]
expr_subset <- expr_matrix[available_genes, meta_data_clean$segement_ROI, drop = FALSE]
expr_subset$Gene <- rownames(expr_subset)
expr_long <- melt(expr_subset, id.vars = "Gene", variable.name = "segement_ROI", value.name = "Expression")
expr_long <- merge(expr_long, meta_data_clean[, c("segement_ROI", "Group")], by = "segement_ROI")

# Boxplot
p <- ggplot(data = expr_long, aes(x = Group, y = Expression)) +
  geom_boxplot(aes(color = Group), width = 0.6, size = 0.6, 
               position = position_dodge(0.6), alpha = 0.5) +
  geom_beeswarm(aes(color = Group, fill = Group), size = 1.2, cex = 4, 
                dodge.width = 0.8, priority = "ascending") +
  facet_wrap(~Gene, scales = "free_y", ncol = 5) +
  labs(x = "", y = "Normalized Expression", title = "Cell adhesion Genes") +
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
    aspect.ratio = 1) + 
  scale_color_manual("Group", values = brewer.pal(9, "Set1")[c(1:2)]) +
  scale_fill_manual("Group", values = brewer.pal(9, "Set1")[c(1:2)]) +
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     label.x = 1.5, paired = FALSE, size = 4)
# Adjust the image size 
n_genes <- length(genes_to_plot)
n_cols <- 5  
n_rows <- ceiling(n_genes / n_cols) 
base_size <- 1.2  # Side length of each small figure
plot_width <- n_cols * base_size + 2
plot_height <- n_rows * base_size + 2

# Save
dir.create("Results/DEG/Boxplot", showWarnings = FALSE, recursive = TRUE)
filename <- "Results/DEG/Boxplot/DSP_diffgene_Boxplot_STAGE_PID_corrected.pdf"
ggsave(p, filename = filename, width = plot_width, height = plot_height)