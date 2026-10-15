rm(list = ls())
install.packages("latex2exp")
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

getwd()
setwd("E:/TNBC_SLN_Program/code_nCounter_DSP_mIF/nCounter_code/nCounter_code/docu")

Group.Info <- read.table("./docu/nCounter_group_infornation.txt", header = TRUE, stringsAsFactors = FALSE) 
Expression.Matrix <- read.table("./docu/nCounter_normalized_matrix.txt",header = TRUE,row.names = 1,sep = "\t",check.names = FALSE)   
exp_df <- read_xlsx("./docu/nCounter_SLN_52sample.xlsx")  
df_clinical_2 <- read_xlsx("./docu/TNBC TDLN cohort.xlsx")
df_clinical <- read_excel("./patient_information.xlsx")
##----------------------------------------------##
## 1. 准备包含协变量的样本信息          --------
##----------------------------------------------##

# 从 exp_df 提取 nCounter ID 与 Sentinel LN Path ID
sample_map <- exp_df %>%
  dplyr::select(nCounterID = `nCounter ID`, 
                SentinelLN = `Sentinel LN Path ID`) %>%
  distinct()

sample_map <- sample_map %>%
  left_join(df_clinical_2 %>%
              dplyr::select(SentinelLN = `Sentinel LN Path ID`, 
                            PatientID = `Patient ID`),
            by = "SentinelLN")


sample_map <- sample_map %>%
  left_join(df_clinical %>%
              dplyr::select(PatientID, `cT(cm)`, `cT(stage)`,Treatment),
            by = "PatientID")


sample_map <- sample_map %>%
  left_join(Group.Info, by = c("nCounterID" = "nCounterID"))


sample_map <- sample_map %>%
  filter(nCounterID %in% colnames(Expression.Matrix),
         !is.na(`cT(cm)`),
         !is.na(`cT(stage)`))


rownames(sample_map) <- sample_map$nCounterID
sample_map <- sample_map[colnames(Expression.Matrix), ]   

sample_map <- sample_map %>%
  mutate(
    Treatment_Class = case_when(
      str_detect(Treatment, "EC|CEF") ~ "Anthracycline_Based",
      str_detect(Treatment, "TP|T$|Paclitaxel") ~ "Taxane_Platinum_Based",
      TRUE ~ "Other"
    ),
    Treatment_Class = as.factor(Treatment_Class)
  )

print(table(sample_map$Treatment_Class))


select_counts <- Expression.Matrix

group_factor <- factor(sample_map$Group, levels = c("pCR", "pNR"))
cT_cm       <- sample_map$`cT(cm)`
cT_stage    <- factor(sample_map$`cT(stage)`)    # 自动处理为虚拟变量
Treatment_Class  <- sample_map$Treatment_Class       


##------------------------------------------------------##
## 2. LIMMA 差异分析（校正 Treatment + cT_cm）  --------##
##------------------------------------------------------##

# 设计矩阵：分组 + 肿瘤大小 
design <- model.matrix(~ 0 + group_factor + Treatment_Class + cT_cm)

colnames(design)[1:2] <- c("pCR", "pNR")   

fit <- lmFit(select_counts, design)

# 对比：pCR vs pNR 
contrast.matrix <- makeContrasts(pCR - pNR, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, robust = TRUE)

DEG_results <- topTable(fit2, n = Inf, adjust = "fdr")

significant_genes <- DEG_results %>%
  rownames_to_column("SYMBOL") %>%
  filter(P.Value < 0.05 & abs(logFC) >= 0.2) %>%
  arrange(P.Value)

write.csv(significant_genes, 
          paste0("./docu/nCounter_Significant_DEGs_cTcm_Treatment_corrected.csv"), 
          row.names = FALSE)

write.csv(DEG_results, 
          paste0("./docu/nCounter_allDEGs_corrected.csv"), 
          row.names = TRUE)