# ==============================================
# Odds Ratio Analysis for Spatial Interaction Data
# Adapted from Python version
# Author: ZQ
# Date: 2025-02-13
# ==============================================

# ---------------------------
# 1. 环境设置与包加载
# ---------------------------
setwd("E:/TNBC_SLN_Program/DataAnalysis/SMI/CellposeTest/")

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(glmnet)
library(broom)

# ---------------------------
# 2. 数据加载与预处理
# ---------------------------
data_t <- read_csv('E:/TNBC_SLN_Program/DataAnalysis/SMI/CellposeTest/CCIResultOutput_20241105/mid_t_region_cci_df_updated.csv')

# 转换Group为二分类变量（pCR=1, 非pCR=0）
data_t <- data_t %>%
  mutate(Group = ifelse(Group == "pCR", 1, 0)) %>%
  drop_na()  # 删除缺失值

# ---------------------------
# 3. 核心函数：计算指定细胞类型的Odds Ratio
# ---------------------------
calculate_odds_ratio <- function(data, cell_type) {
  # 过滤指定细胞类型的数据
  cell_data <- data %>%
    filter(Final_cell_type == cell_type)
  
  # 准备模型数据
  x <- as.matrix(dplyr::select(cell_data, -CellID, -Group, -Final_cell_type, -sample_fov))
  y <- cell_data$Group
  
  # 使用glmnet进行正则化逻辑回归
  cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 0.5)
  optimal_lambda <- cv_fit$lambda.min
  final_model <- glmnet(x, y, family = "binomial", alpha = 0.5, lambda = optimal_lambda)
  
  # 提取系数并计算Odds Ratio
  model_coefs <- coef(final_model)
  odds_ratio <- exp(model_coefs[, 1])
  
  # 使用标准glm获取标准误和p值
  cell_data$Group <- as.factor(cell_data$Group)
  glm_data <- cell_data %>%
    dplyr::select(-CellID, -sample_fov, -Final_cell_type)
  
  glm_model <- glm(Group ~ ., data = glm_data, family = "binomial")
  glm_summary <- summary(glm_model)
  
  # 计算置信区间
  se <- sqrt(diag(vcov(glm_model)))
  lower_ci <- exp(model_coefs[, 1] - 1.96 * se)
  upper_ci <- exp(model_coefs[, 1] + 1.96 * se)
  p_values <- glm_summary$coefficients[, 4]
  
  # 整理结果
  results <- data.frame(
    Variable = rownames(model_coefs),
    Odds_Ratio = odds_ratio,
    Lower_CI = lower_ci,
    Upper_CI = upper_ci,
    P_Value = p_values
  )
  
  # 校正p值（BH方法）
  results$P_Adjusted <- p.adjust(results$P_Value, method = "BH")
  
  # 移除截距项
  results <- results %>%
    filter(Variable != "(Intercept)")
  
  return(results)
}

# ---------------------------
# 4. 可视化函数：绘制棒棒糖图
# ---------------------------
plot_odds_ratio <- function(odds_ratio_results, cell_type, color_code, save_path) {
  # 定义颜色映射（基于细胞类型）
  color_map <- c(
    "CD4_Tnaive" = "#99aa00",
    "CD4_Treg" = "#99cc00",
    "CD4_Tmem" = "#99ee00",
    "CD8_Tnaive" = "#0022FF",
    "CD8_Tpex" = "#0044FF",
    "CD8_Tmem" = "#0066FF",
    "Immature_B" = "#FFD799",
    "Naïve_B" = "#FFD777",
    "None_class_switched_memoryB" = "#ffcc00",
    "Class_switched_memoryB" = "#ffaa00",
    "cDC" = "#aa11a6",
    "Macrophage_monocyte" = "#FF00FF",
    "Plasma" = "#5F86C4",
    "NK" = "#0099FF",
    "pDC" = "#BEDADA",
    "Granulocytes" = "#6699BB",
    "Adipose" = "#99776B",
    "Epithelial" = "#f7a899",
    "Fibroblast" = "#e5486e",
    "Endothelial" = "#CD8A66FF"
  )
  
  # 数据预处理
  plot_data <- odds_ratio_results %>%
    mutate(
      Significance = ifelse(P_Adjusted < 0.05, "*", ""),
      Variable = factor(Variable, levels = unique(Variable))
    )
  
  # 创建图形
  p <- ggplot(plot_data, aes(x = Odds_Ratio, y = Variable)) +
    # 置信区间线段
    geom_segment(aes(x = Lower_CI, xend = Upper_CI, yend = Variable), size = 0.5) +
    
    # 中心点（根据细胞类型着色）
    geom_point(size = 4, shape = 21, stroke = 0, fill = color_code) +
    
    # 显著性边框
    geom_point(size = 4, shape = 21, 
               aes(stroke = ifelse(P_Adjusted < 0.05, 1.2, 0)),
               color = "black", fill = NA) +
    
    # 显著性标注
    geom_text(aes(label = Significance), 
              hjust = 0.5, vjust = 0.7, size = 5, 
              fontface = "bold", color = "white") +
    
    # 右侧p值标注
    geom_text(aes(label = ifelse(P_Adjusted < 0.05, 
                                 sprintf("%.2e", P_Adjusted), 
                                 sprintf("%.2f", P_Adjusted))),
              hjust = 0, vjust = 0.5, size = 3.5,
              color = ifelse(plot_data$P_Adjusted < 0.05, "black", "grey"),
              x = max(plot_data$Upper_CI) * 1.15) +
    
    # 左侧颜色标签
    geom_tile(aes(x = -0.5, y = Variable, fill = Variable), 
              width = 0.1, height = 0.8, show.legend = FALSE) +
    
    # 参考线
    geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
    
    # 方向箭头
    annotate("segment", x = 1, xend = 0.5, 
             y = nrow(plot_data) + 0.5, yend = nrow(plot_data) + 0.5,
             arrow = arrow(length = unit(0.02, "npc")), color = "blue") +
    annotate("segment", x = 1, xend = 1.5, 
             y = nrow(plot_data) + 0.5, yend = nrow(plot_data) + 0.5,
             arrow = arrow(length = unit(0.02, "npc")), color = "red") +
    annotate("text", x = 0.5, y = nrow(plot_data) + 0.5, 
             label = "pNR", color = "blue", hjust = 0, vjust = -0.5) +
    annotate("text", x = 1.5, y = nrow(plot_data) + 0.5, 
             label = "pCR", color = "red", hjust = 1, vjust = -0.5) +
    
    # 颜色映射
    scale_fill_manual(values = color_map) +
    
    # 主题设置
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12)
    ) +
    
    # 坐标轴标签
    labs(
      title = paste("Odds ratio of", cell_type, "contact cell types in T region"),
      x = "Odds Ratio",
      y = paste(cell_type, "direct contact cell proportion")
    )
  
  # 保存图形
  ggsave(paste0(save_path, "/T_region_", cell_type, "_lollipop_plot.png"), 
         p, width = 7, height = 6, dpi = 300)
  ggsave(paste0(save_path, "/T_region_", cell_type, "_lollipop_plot.pdf"), 
         p, width = 7, height = 6)
  
  return(p)
}

# ---------------------------
# 5. 主分析流程
# ---------------------------
# 设置输出目录
output_dir <- "./CCIoddsRatio/OddsRatioOutput_sample_fov"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 定义要分析的细胞类型及其颜色代码
cell_analysis <- list(
  cDC = "#aa11a6",
  pDC = "#BEDADA"
  # 可在此添加更多细胞类型
)

# 循环分析每个细胞类型
for (cell in names(cell_analysis)) {
  cat("\n分析细胞类型:", cell, "\n")
  
  # 计算Odds Ratio
  results <- calculate_odds_ratio(data_t, cell)
  
  # 保存结果
  write_csv(results, paste0(output_dir, "/T_region_", cell, "_odds_ratio_results.csv"))
  
  # 绘制图形
  plot_odds_ratio(results, cell, cell_analysis[[cell]], output_dir)
  
  # 打印简要统计
  cat("显著变量数量:", sum(results$P_Adjusted < 0.05), "\n")
  print(results %>% arrange(P_Adjusted) %>% head(5))
}

# ---------------------------
# 6. 结果汇总
# ---------------------------
cat("\n分析完成！所有结果已保存至:", output_dir, "\n")