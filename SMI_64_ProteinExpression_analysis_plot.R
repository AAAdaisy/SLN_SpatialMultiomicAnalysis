################### Next part  ##########

################### Comparing DEGs of cDC in DC_enriched niche vs cDC in other niches --> (compare cDC protein expression) #######
## 0. READ DATA -- Region_df / T_region_CN100 / expression_df
library(tidyverse)
library(data.table)
library(rstatix)
library(ggthemes)
library(ggh4x)
library(ggpubr)
library(aplot)
library(dplyr)

set.seed(123)   # 保证每次运行结果一致

cci_df <- read_csv('./CCI_df_with_sample_fov.csv')
T_region_CN100 <- read.csv("./T_Region_window.100CN.csv")

expression_df <- read.csv("./expression_matrix.csv")
expression_df <- expression_df %>%
  mutate(X = as.character(X)) %>%  
  tibble::column_to_rownames(var = "X") 

## 1. 取 T 区 cDC 
t_region_sample_fov <- Region_df %>% 
  filter(Region_type == "T_region") %>% 
  pull(sample_fov)

cDC_T_region_DCclusterniche <- T_region_CN100 %>%  # DCclusterniche neighborhood100 == 10
  filter(neighborhood100 == "10",
         sample_fov %in% t_region_sample_fov)

cDC_T_region_notDCclusterniche <- T_region_CN100 %>%  # DCclusterniche neighborhood100 != 10
  filter(neighborhood100 != "10",
         sample_fov %in% t_region_sample_fov)

# 从cDC_T_region_notDCclusterniche随机抽取与cDC_T_region_DCclusterniche相同行数的细胞

# 目标行数
n_target <- nrow(cDC_T_region_DCclusterniche)

cDC_T_region_notDCclusterniche_sub <- cDC_T_region_notDCclusterniche %>%
  sample_n(size = n_target)   # 随机抽取相同行数
nrow(cDC_T_region_notDCclusterniche_sub)

library(dplyr)

combined_df <- bind_rows(
  cDC_T_region_DCclusterniche %>% mutate(cluster = "DC_cluster"),
  cDC_T_region_notDCclusterniche_sub %>% mutate(cluster = "notDC_cluster")
)

## 2. 标记“是否处于DC_enrichedniche” 

## DCclusterniche的 cell_id  列表
cDC_T_region_DCclusterniche_id <- combined_df %>% 
  filter(neighborhood100 == "10") %>% 
  distinct(cell_id ) %>% 
  pull(cell_id )

cDC_T_region_notDCclusterniche_id <- combined_df %>% 
  filter(neighborhood100 != "10") %>% 
  distinct(cell_id ) %>% 
  pull(cell_id )   

## 给 T 区全部 cDC 打标签
cDC_T_region_labelled <- combined_df %>% 
  distinct(cell_id, .keep_all = TRUE) %>%   
  mutate(DCclusterniche = if_else(cell_id %in% cDC_T_region_DCclusterniche_id, "DCclusterniche+", "DCclusterniche-"))


## 3. 提取表达矩阵 
cDC_expression <- expression_df %>% 
  ## 只保留 T 区 cDC
  rownames_to_column("cell_id") %>% 
  filter(cell_id %in% cDC_T_region_labelled$cell_id) %>% 
  ## 合并标签
  left_join(cDC_T_region_labelled %>% select(cell_id, DCclusterniche),
            by = c("cell_id" = "cell_id"))

## 4. 背景归一化
## 4.1 去掉非功能蛋白 + 内参
cDC_expr_clean <- cDC_expression %>% 
  select(-c(EpCAM, CD31, Fibronectin, FABP4, CD138, CD45, CD15,
            CD11c, CD123, CD56, CD68, CD19, CD20, CD127, IgD,
            CD3, CD4, CD8, CD45RA, FOXP3, PD.1, TCF7))

## 4.2 Background 计算
## 背景归一化（只处理数值型蛋白列）
cDC_norm <- cDC_expr_clean %>% 
  ## 计算背景
  mutate(Background = (Rb.IgG + Ms.IgG1) / 2) %>% 
  ## 仅对数值型列做归一化
  mutate(across(where(is.numeric), ~ .x / Background)) %>% 
  ## 去掉背景列
  select(-Background) %>% 
  ## 把分组信息粘回来
  relocate(DCclusterniche, .after = last_col())  


## 保存中间文件（可选）
write.csv(cDC_norm,
          file = "./Output/Tregion/cDC_DCclusterniche_vs_notDCclusterniche_40protein-NormalizedBackground.csv",
          row.names = FALSE)

## 5. 差异分析 + 作图  
df <- fread("./Output/Tregion/cDC_DCclusterniche_vs_notDCclusterniche_40protein-NormalizedBackground.csv")

## 5.1 boxplot 数据准备
df_long <- df %>% 
  pivot_longer(-c(cell_id, DCclusterniche),
               names_to = "Protein",
               values_to = "Expression") %>% 
  mutate(DCclusterniche = factor(DCclusterniche, levels = c("DCclusterniche+", "DCclusterniche-")))

## 5.1.1. 统计检验：每个蛋白做 Wilcoxon
stat_df <- df_long %>%
  group_by(Protein) %>%
  wilcox_test(Expression ~ DCclusterniche, detailed = TRUE) %>%   # 返回 estimate 等
  mutate(
    Enrich = case_when(
      p <= 0.05 & estimate > 0 ~ "DCclusterniche",      # DCclusterniche+ > DCclusterniche-
      p <= 0.05 & estimate < 0 ~ "notDCclusterniche",    # DCclusterniche+ < DCclusterniche-
      .default = "n.s."
    )
  ) %>%
  select(Protein, p, Enrich, estimate)

## 5.1.2. 计算均值（lollipop 用）
lipo_mean <- df_long %>%
  reframe(
    Exp = mean(Expression),
    .by = c(Protein, DCclusterniche)
  ) %>%
  pivot_wider(names_from = DCclusterniche, values_from = Exp) %>%
  left_join(stat_df, by = "Protein") %>%
  # 保证列存在
  mutate(across(c(`DCclusterniche+`, `DCclusterniche-`), ~ replace_na(.x, 0))) %>%
  arrange(Enrich, desc(`DCclusterniche+`), p) %>%
  mutate(Protein = factor(Protein, levels = rev(unique(.$Protein))))

## 5.1.3. boxplot 数据（按蛋白排序，和 lollipop 保持一致）
df_plot <- df_long %>%
  left_join(lipo_mean %>% select(Protein, Enrich), by = "Protein") %>%
  mutate(
    Protein = factor(Protein, levels = levels(lipo_mean$Protein)),
    DCclusterniche = factor(DCclusterniche, levels = c("DCclusterniche+", "DCclusterniche-"))
  )

## 6. 画图 
## 6.1 boxplot
p_box <- ggplot(df_plot, aes(DCclusterniche, Expression, fill = DCclusterniche)) +
  geom_boxplot(outlier.size = 0.4, outlier.color = "black", width = 0.6) +
  geom_point(size = 1, alpha = 0.5) +
  facet_wrap(~Protein, scales = "free_y", nrow = 8) +
  ggthemes::theme_few() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.border = element_rect(fill = NA, color = "black", size = 1)
  ) +
  ggpubr::stat_compare_means(comparisons = list(c("DCclusterniche+", "DCclusterniche-")),
                             method = "wilcox", label = "p.format",
                             size = 3, label.y = 0) +
  scale_fill_manual(values = c("#BC3C29FF", "#0072B5FF")) +
  labs(x = "", y = "Expression (background normalized)")

ggsave(p_box,
       filename = "./Output/Tregion/cDC_DCclusterniche_vs_notDCclusterniche_Boxplot.pdf",
       height = 21, width = 18)
ggsave(p_box,
       filename = "./Output/Tregion/cDC_DCclusterniche_vs_notDCclusterniche_Boxplot.png",
       height = 21, width = 18, dpi = 300)

## 6.2 lollipop plot
p_lipo_main <- ggplot(lipo_mean) +
  geom_segment(aes(x = `DCclusterniche+`, xend = `DCclusterniche-`, y = Protein, yend = Protein),
               color = "black", linewidth = 0.8) +
  # 用映射方式指定颜色
  geom_point(aes(x = `DCclusterniche+`, y = Protein, color = "DCclusterniche+"), size = 2.5) +
  geom_point(aes(x = `DCclusterniche-`, y = Protein, color = "DCclusterniche-"), size = 2.5) +
  # originate scale
  #scale_color_manual(values = c("DCclusterniche+" = "#BC3C29FF", "DCclusterniche-" = "#0072B5FF")) +
  # labled scale
  scale_color_manual(
    values = c("DCclusterniche+" = "#da70d6", "DCclusterniche-" = "#4169e1"),
    labels = c("DCclusterniche+" = "cDCs in DCclusterniche",
               "DCclusterniche-" = "cDCs not in DCclusterniche")
  ) +
  ## 经典无网格主题，方便全黑控制
  theme_classic(base_size = 14) +
  theme(
    axis.line    = element_line(color = "black"),   # 轴线黑色 , #size = 1
    axis.text.x  = element_text(color = "black", size = 12),
    axis.text.y  = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12),
    axis.title.y = element_blank(),                          # 已写在 labs()
    legend.position = "right",                               # 独立右侧 legend
    legend.title    = element_blank(),
    legend.text     = element_text(color = "black", size = 10),
    #panel.border    = element_blank()                        # 去掉多余边框
    #加黑色边框
    panel.border    = element_rect(fill = NA, color = "black", size = 1)
  ) +
  labs(x = "Mean expression (background normalized)")

## p 值文字（同你原来，无需改）
p_lipo_text <- ggplot(lipo_mean, aes(x = 0, y = Protein)) +
  geom_text(aes(label = ifelse(p > 0.001,
                               sprintf("%.2g", p),
                               "1e-3"),
                color = Enrich),
            size = 2.5, hjust = 0, nudge_x = 0.02) +
  scale_color_manual(values = c("#da70d6", "grey50", "#4169e1")) +
  theme_void() +
  theme(legend.position = "none")

## 合并
p_lipo <- aplot::insert_right(p_lipo_main, p_lipo_text, width = 0.15)

## 保存
ggsave(p_lipo,
       filename = "./Output/Tregion/cDC_DCclusterniche_vs_notDCclusterniche_Lollipop.pdf",
       height = 8, width = 8)

ggsave(p_lipo,
       filename = "./Output/Tregion/cDC_DCclusterniche_vs_notDCclusterniche_Lollipop.png",
       dpi = 300, height = 8, width = 8)

## ---- 7. 统计结果整理
sum_res <- df_long %>% 
  group_by(Protein, DCclusterniche) %>% 
  wilcox_test(Expression ~ DCclusterniche, detailed = TRUE) %>% 
  summarise(
    mean = mean(Expression),
    median = median(Expression),
    q25  = quantile(Expression, 0.25),
    q75  = quantile(Expression, 0.75),
    .groups = "drop"
  ) %>% 
  left_join(stat_df %>% select(Protein, p), by = "Protein") %>% 
  mutate(q_signif = case_when(
    is.na(p)        ~ "",
    p < 0.001       ~ "***",
    p < 0.01        ~ "**",
    p < 0.05        ~ "*",
    TRUE            ~ ""
  ))

## 宽表：一行一个蛋白
sum_wide <- sum_res %>% 
  pivot_wider(names_from = DCclusterniche,
              values_from = c(mean, median, q25, q75),
              names_sep = "_")

fwrite(sum_wide,
       file = "./Output/Tregion/cDC_DCclusterniche_vs_notDCclusterniche_Statistics_wilcoxtest.csv")

##### select specific proteins to plot #####
## 4. 背景归一化
## 4.1 去掉非功能蛋白 + 内参
cDC_expr_clean <- cDC_expression %>% 
  select(-c(EpCAM, CD31, Fibronectin, FABP4, CD138, CD45, CD15,
            CD11c, CD123, CD56, CD68, CD19, CD20, CD127, IgD,
            CD3, CD4, CD8, CD45RA, FOXP3, PD.1, TCF7,
            # Ms.IgG1, Rb.IgG
            
  ))

## 4.2 Background 计算
## 背景归一化（只处理数值型蛋白列）
cDC_norm <- cDC_expr_clean %>% 
  ## 计算背景
  mutate(Background = (Rb.IgG + Ms.IgG1) / 2) %>% 
  ## 仅对数值型列做归一化（跳过 cell_id 和 Contact_CD4）
  mutate(across(where(is.numeric), ~ .x / Background)) %>% 
  ## 去掉背景列
  select(-Background) %>% 
  ## 把分组信息粘回来
  relocate(DCclusterniche, .after = last_col())  


## 保存中间文件（可选）
write.csv(cDC_norm,
          file = "./Output/Tregion/cDC_DCclusterniche_vs_notDCclusterniche_40protein_select-NormalizedBackground.csv",
          row.names = FALSE)

## 5. 差异分析 + 作图  
df <- fread("./Output/Tregion/cDC_DCclusterniche_vs_notDCclusterniche_40protein_select-NormalizedBackground.csv")

## 5.1 boxplot 数据准备
df_long <- df %>% 
  pivot_longer(-c(cell_id, DCclusterniche),
               names_to = "Protein",
               values_to = "Expression") %>% 
  mutate(DCclusterniche = factor(DCclusterniche, levels = c("DCclusterniche+", "DCclusterniche-")))

## 5.1.1. 统计检验：每个蛋白做 Wilcoxon
stat_df <- df_long %>%
  group_by(Protein) %>%
  wilcox_test(Expression ~ DCclusterniche, detailed = TRUE) %>%   # 返回 estimate 等
  mutate(
    Enrich = case_when(
      p <= 0.05 & estimate > 0 ~ "DCclusterniche",      # DCclusterniche+ > DCclusterniche-
      p <= 0.05 & estimate < 0 ~ "notDCclusterniche",    # DCclusterniche+ < DCclusterniche-
      .default = "n.s."
    )
  ) %>%
  select(Protein, p, Enrich, estimate)

## 5.1.2. 计算均值（lollipop 用）
lipo_mean <- df_long %>%
  reframe(
    Exp = mean(Expression),
    .by = c(Protein, DCclusterniche)
  ) %>%
  # 筛选需要绘制的蛋白
  filter(Protein %in% c("PD.L1","CCR7","HLA.DR","Tim.3")) %>%
  
  pivot_wider(names_from = DCclusterniche, values_from = Exp) %>%
  left_join(stat_df, by = "Protein") %>%
  # 保证列存在
  mutate(across(c(`DCclusterniche+`, `DCclusterniche-`), ~ replace_na(.x, 0))) %>%
  arrange(Enrich, desc(`DCclusterniche+`), p) %>%
  mutate(Protein = factor(Protein, levels = rev(unique(.$Protein))))

## 5.1.3. boxplot 数据（按蛋白排序，和 lollipop 保持一致）
df_plot <- df_long %>%
  left_join(lipo_mean %>% select(Protein, Enrich), by = "Protein") %>%
  # 筛选需要绘制的蛋白
  filter(Protein %in% c("PD.L1","CCR7","HLA.DR","Tim.3")) %>%
  mutate(
    Protein = factor(Protein, levels = levels(lipo_mean$Protein)),
    DCclusterniche = factor(DCclusterniche, levels = c("DCclusterniche+", "DCclusterniche-"))
  )

## 6. 画图 
## 6.1 boxplot
p_box <- ggplot(df_plot, aes(DCclusterniche, Expression, fill = DCclusterniche)) +
  geom_boxplot(outlier.size = 0.4, outlier.color = "black", width = 0.6) +
  geom_point(size = 1, alpha = 0.5) +
  facet_wrap(~Protein, scales = "free_y", nrow = 1) +
  ggthemes::theme_few() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.border = element_rect(fill = NA, color = "black", size = 1)
  ) +
  ggpubr::stat_compare_means(comparisons = list(c("DCclusterniche+", "DCclusterniche-")),
                             method = "wilcox", label = "p.format",
                             size = 3, label.y = 0) +
  scale_fill_manual(values = c("#da70d6", "#4169e1")) +
  labs(x = "", y = "Expression (background normalized)")

ggsave(p_box,
       filename = "./Output/Tregion/cDC_DCclusterniche_vs_notDCclusterniche_select_Boxplot.pdf",
       height = 4, width = 12)
ggsave(p_box,
       filename = "./Output/Tregion/cDC_DCclusterniche_vs_notDCclusterniche_select_Boxplot.png",
       height = 4, width = 12, dpi = 300)

## 6.2 lollipop plot
p_lipo_main <- ggplot(lipo_mean) +
  geom_segment(aes(x = `DCclusterniche+`, xend = `DCclusterniche-`, y = Protein, yend = Protein),
               color = "black", linewidth = 0.8) +
  # 用映射方式指定颜色，才能自动出 legend
  geom_point(aes(x = `DCclusterniche+`, y = Protein, color = "DCclusterniche+"), size = 2.5) +
  geom_point(aes(x = `DCclusterniche-`, y = Protein, color = "DCclusterniche-"), size = 2.5) +
  # originate scale
  #scale_color_manual(values = c("DCclusterniche+" = "#BC3C29FF", "DCclusterniche-" = "#0072B5FF")) +
  # labled scale
  scale_color_manual(
    values = c("DCclusterniche+" = "#da70d6", "DCclusterniche-" = "#4169e1"),
    labels = c("DCclusterniche+" = "cDCs in DCclusterniche",
               "DCclusterniche-" = "cDCs not in DCclusterniche")
  ) +
  ## 经典无网格主题，方便全黑控制
  theme_classic(base_size = 14) +
  theme(
    axis.line    = element_line(color = "black"),   # 轴线黑色 , #size = 1
    axis.text.x  = element_text(color = "black", size = 12),
    axis.text.y  = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12),
    axis.title.y = element_blank(),                          # 已写在 labs()
    legend.position = "right",                               # 独立右侧 legend
    legend.title    = element_blank(),
    legend.text     = element_text(color = "black", size = 10),
    #panel.border    = element_blank()                        # 去掉多余边框
    #加黑色边框
    panel.border    = element_rect(fill = NA, color = "black", size = 1)
  ) +
  labs(x = "Mean expression (background normalized)")

## p 值文字
p_lipo_text <- ggplot(lipo_mean, aes(x = 0, y = Protein)) +
  geom_text(aes(label = ifelse(p > 0.001,
                               sprintf("%.2g", p),
                               "1e-3"),
                color = Enrich),
            size = 2.5, hjust = 0, nudge_x = 0.02) +
  scale_color_manual(values = c("#da70d6", #"grey50", 
                                "#4169e1")) +
  theme_void() +
  theme(legend.position = "none")

## 合并
p_lipo <- aplot::insert_right(p_lipo_main, p_lipo_text, width = 0.15)

## 保存
ggsave(p_lipo,
       filename = "./Output/Tregion/cDC_DCclusterniche_vs_notDCclusterniche_select_Lollipop.pdf",
       height = 3, width = 8)

ggsave(p_lipo,
       filename = "./Output/Tregion/cDC_DCclusterniche_vs_notDCclusterniche_select_Lollipop.png",
       dpi = 300, height = 3, width = 8)
