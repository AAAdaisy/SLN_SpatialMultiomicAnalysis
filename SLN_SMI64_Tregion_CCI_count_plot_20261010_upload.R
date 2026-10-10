## T region cDC/pDC contact Cells  (0, 1, 2, >=3)

library(data.table)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(cowplot)

#setwd("E:/TNBC_SLN_Program/DataAnalysis/SMI/CellposeTest/CCIResultOutput_20241105")

# 1. 设定输出目录
new_output_dir <- "./RegionCCI_output_with_sourcedata/"
if (!dir.exists(new_output_dir)) {
  dir.create(new_output_dir, recursive = TRUE)
}

# 2. 颜色方案定义
colorlist3 <- c("#BC3C29FF", "#0072B5FF")
colorlist4 <- c("#BC3C29FF", "#0072B5FF")

# 3. 封装绘图及 Source Data 导出函数
run_analysis_and_export <- function(data_path, target_celltype, color_pal, output_dir) {
  
  if (!file.exists(data_path)) {
    warning(paste("File not found:", data_path))
    return(NULL)
  }
  
  dat <- fread(data_path)
  count_levels <- c("0", "1", "2", ">=3")
  
  # ----------------------------------------------------
  # Part A: 导出单图及对应 Source Data
  # ----------------------------------------------------
  p1_data <- dat %>% 
    filter(Yaxis == target_celltype) %>%
    mutate(count_group = factor(count_group, levels = count_levels))
  
  if (nrow(p1_data) > 0) {
    # 绘图
    p1 <- ggplot(p1_data, aes(x = Group, y = proportion, fill = Group)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.01) +
      geom_point(aes(color = Group), 
                 position = position_jitterdodge(dodge.width = 1, jitter.width = 0.3, jitter.height = 0), 
                 size = 1.6) +
      facet_wrap(~count_group, scales = "free_y", nrow = 1) +
      scale_color_manual("", values = color_pal[1:2]) +
      scale_fill_manual("", values = color_pal[1:2]) +
      labs(x = "", y = "Proportion") +
      ggthemes::theme_few() +
      theme(
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 9),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1.2, linetype = "solid"),
        legend.title = element_text(size = 9),
        plot.title = element_text(size = 9, hjust = 0.5, colour = "black"),
        strip.text = element_text(size = 9, colour = "black"),
        legend.text = element_text(size = 9)
      ) +
      ggtitle(paste0(target_celltype, "_T_region"))
    
    # 保存图像
    ggsave(file.path(output_dir, paste0(target_celltype, "_T_region.pdf")), p1, height = 2, width = 6)
    ggsave(file.path(output_dir, paste0(target_celltype, "_T_region.jpg")), p1, height = 2, width = 6, dpi = 300)
    
    # 导出 Source Data (点的散点明细值)
    fwrite(p1_data, file.path(output_dir, paste0("SourceData_", target_celltype, "_T_region_facet_points.csv")))
  }
  
  # ----------------------------------------------------
  # Part B: 遍历循环 celltype，导出拼图与逐个小图的 Source Data
  # ----------------------------------------------------
  for (celltype in unique(dat$Yaxis)) {
    plot_list <- list()
    summary_stat_list <- list()
    point_data_list <- list()
    
    for (cnt in count_levels) {
      mid <- dat %>% 
        filter(Yaxis == celltype) %>% 
        filter(count_group == cnt)
      
      # 数据为空时跳过
      if (nrow(mid) == 0) next
      
      # 收集点级别明细数据
      point_data_list[[cnt]] <- mid
      
      # 收集箱线图统计量
      summary_stat_list[[cnt]] <- mid %>%
        group_by(Group, count_group) %>%
        summarise(
          n = n(),
          min = min(proportion, na.rm = TRUE),
          q25 = quantile(proportion, 0.25, na.rm = TRUE),
          median = median(proportion, na.rm = TRUE),
          mean = mean(proportion, na.rm = TRUE),
          q75 = quantile(proportion, 0.75, na.rm = TRUE),
          max = max(proportion, na.rm = TRUE),
          p_value_label = ifelse("p_value_label" %in% colnames(mid), unique(p_value_label)[1], NA),
          .groups = "drop"
        )
      
      # 绘制单个小图
      p_sub <- ggplot(mid, aes(x = Group, y = proportion, fill = Group)) +
        geom_boxplot(alpha = 0.7, outlier.size = 0.01) +
        geom_point(aes(color = Group), 
                   position = position_jitterdodge(dodge.width = 1, jitter.width = 0.3, jitter.height = 0), 
                   size = 1.6) +
        scale_color_manual("", values = color_pal[1:2]) +
        scale_fill_manual("", values = color_pal[1:2]) +
        labs(x = "", y = "Proportion") +
        ggthemes::theme_few() +
        theme(
          axis.text = element_text(color = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 9),
          axis.title = element_text(size = 9),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 1.2, linetype = "solid"),
          legend.title = element_text(size = 9),
          plot.title = element_text(size = 9, hjust = 0.5, colour = "black"),
          strip.text = element_text(size = 9, colour = "black"),
          legend.text = element_text(size = 9)
        ) +
        ggtitle(cnt) +
        annotate("segment", x = 1, xend = 2, 
                 y = max(mid$proportion, na.rm = TRUE), 
                 yend = max(mid$proportion, na.rm = TRUE),
                 arrow = arrow(ends = "both", angle = 90, length = unit(.04, "cm"))) +
        annotate("text", x = 1.5, 
                 y = max(mid$proportion, na.rm = TRUE), 
                 label = unique(mid$p_value_label), size = 4) +
        ylim(NA, max(mid$proportion, na.rm = TRUE) + max(mid$proportion, na.rm = TRUE) / 40)
      
      plot_list[[cnt]] <- p_sub
    }
    
    # 拼图并导出图像
    if (length(plot_list) > 0) {
      p_combined <- cowplot::plot_grid(plotlist = plot_list, nrow = 1, align = 'hv')
      ggsave(file.path(output_dir, paste0(celltype, "_T_region.pdf")), p_combined, height = 2, width = 10)
      ggsave(file.path(output_dir, paste0(celltype, "_T_region.jpg")), p_combined, height = 2, width = 10, dpi = 300)
    }
    
    # 导出各个小图所有散点的 Source Data
    if (length(point_data_list) > 0) {
      all_points_df <- bind_rows(point_data_list)
      fwrite(all_points_df, file.path(output_dir, paste0("SourceData_", celltype, "_dots_all_counts.csv")))
    }
    
    # 导出各个小图箱线图统计量 (Mean/Median/IQR/P-value)
    if (length(summary_stat_list) > 0) {
      all_stats_df <- bind_rows(summary_stat_list)
      fwrite(all_stats_df, file.path(output_dir, paste0("SourceData_", celltype, "_boxplot_statistics.csv")))
    }
  }
}

# -------------------------------------------------------------------------------------

# Task 1: cDC-CD4_Tnaive
run_analysis_and_export(
  data_path = "./RegionCCI_output_count_batch_barplot/t_region_cDC_InteractedCount-data-CD4Tnaive.csv",
  target_celltype = "cDC-CD4_Tnaive",
  color_pal = colorlist3,
  output_dir = new_output_dir
)

# Task 2: pDC-CD4_Tnaive
run_analysis_and_export(
  data_path = "./RegionCCI_output_count_batch_barplot/t_region_pDC_InteractedCount-data-CD4Tnaive.csv",
  target_celltype = "pDC-CD4_Tnaive",
  color_pal = colorlist4,
  output_dir = new_output_dir
)