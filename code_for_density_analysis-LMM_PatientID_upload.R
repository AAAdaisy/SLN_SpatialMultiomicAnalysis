#处理数据，提取 PatientID 并合并临床信息
library(dplyr)
library(stringr)
library(readxl)
library(readr)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(tidyr)
library(purrr)


### 20261012 only correct patient_ID proportion #####
df_clinical <- read_excel("F:/0_zq/patient_information.xlsx")

df <- read.csv("E:/TNBC_SLN_Program/DataAnalysis/SMI/CellTypeProportionPlot/RegionAnalysis/Tregion_CellType-Proportion-Sample_fov-Data.csv")

unique(df$sample_fov)
head(df)

df_proportion_clean <- df %>%
  filter(Group %in% c("pCR", "pNR")) %>%
  
  mutate(PatientID_clean = str_replace(sample_fov, "_.*", "")) %>%
  mutate(PatientID_clean = str_replace(PatientID_clean, "^S", "")) %>%
  
  left_join(df_clinical, by = c("PatientID_clean" = "PatientID")) %>%
  
  rename(cT_cm = `cT(cm)`,
         Stage_clean = Stage) %>%
  
  filter(!is.na(cT_cm) & !is.na(Stage_clean)) %>%
  
  mutate(
    Group = factor(Group.x, levels = c("pNR", "pCR")),
    PatientID_clean = as.factor(PatientID_clean),
    Stage_clean = as.factor(Stage_clean)
  )


Proportion_optimal_lmm <- df_proportion_clean %>%
  group_by(Final_cell_type) %>%
  nest() %>%
  mutate(
    model = map(data, ~ {
      tryCatch(
        # 将降维后的 Treatment_Class 作为随机效应
        #lmer(Proportion ~ Group + (1 | Treatment_class) + (1 | PatientID_clean) + cT_stage, data = .x),
        lmer(Proportion ~ Group + (1 | PatientID_clean), data = .x),
        error = function(e) return(NULL)
      )
    }),
    tidied = map(model, ~ if(!is.null(.x)) tidy(.x, effects = "fixed") else NULL)
  ) %>%
  filter(!map_lgl(tidied, is.null)) %>%
  unnest(tidied) %>%
  filter(term == "GrouppCR") %>%
  select(Final_cell_type, estimate, std.error, p.value) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  arrange(p.value)

print("--- 模型结果 ---")
print(Proportion_optimal_lmm)
write_csv(Proportion_optimal_lmm, "E:/TNBC_SLN_Program/DataAnalysis/SMI/CellTypeProportionPlot/RegionAnalysis/Tregion_FinalCellType_Proportion_LMM_PaientID_Result.csv")

### 20261012 only correct patient_ID density #####
df_clinical <- read_excel("F:/0_zq/patient_information.xlsx")

df <- read.csv("E:/TNBC_SLN_Program/DataAnalysis/SMI/CellTypeProportionPlot/RegionAnalysis/Tregion_CellType-Proportion-Sample_fov-Data.csv")

fov_area_mm2 <- 0.2601

unique(df$sample_fov)

df_density_clean <- df %>%
  mutate(Density = Count / fov_area_mm2) %>%
  filter(Group %in% c("pCR", "pNR")) %>%
  
  mutate(PatientID_clean = str_replace(sample_fov, "_.*", "")) %>%
  mutate(PatientID_clean = str_replace(PatientID_clean, "^S", "")) %>%
  
  left_join(df_clinical, by = c("PatientID_clean" = "PatientID")) %>%
  
  rename(cT_cm = `cT(cm)`,
         Stage_clean = Stage) %>%
  
  filter(!is.na(cT_cm) & !is.na(Stage_clean)) %>%
  
  mutate(
    Group = factor(Group.x, levels = c("pNR", "pCR")),
    PatientID_clean = as.factor(PatientID_clean),
    Stage_clean = as.factor(Stage_clean)
  )

density_optimal_lmm <- df_density_clean %>%
  group_by(Final_cell_type) %>%
  nest() %>%
  mutate(
    model = map(data, ~ {
      tryCatch(
        # 将降维后的 Treatment_Class 作为随机效应
        #lmer(Density ~ Group + (1 | Treatment_class) + (1 | PatientID_clean) + cT_stage, data = .x),
        lmer(Density ~ Group + (1 | PatientID_clean), data = .x),
        error = function(e) return(NULL)
      )
    }),
    tidied = map(model, ~ if(!is.null(.x)) tidy(.x, effects = "fixed") else NULL)
  ) %>%
  filter(!map_lgl(tidied, is.null)) %>%
  unnest(tidied) %>%
  filter(term == "GrouppCR") %>%
  select(Final_cell_type, estimate, std.error, p.value) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  arrange(p.value)

print("--- 模型结果 ---")
print(density_optimal_lmm)
write_csv(density_optimal_lmm, "E:/TNBC_SLN_Program/DataAnalysis/SMI/CellTypeProportionPlot/RegionAnalysis/Tregion_FinalCellType_Density_LMM_PaientID_Result.csv")

