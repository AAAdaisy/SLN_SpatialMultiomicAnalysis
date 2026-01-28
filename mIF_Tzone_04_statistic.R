##################################################
## Project: mIF SLN Tregion
## Script purpose: DC 20um neighbors statistic
## Date: 25-11-12
##################################################
rm(list = ls())
gc()
library(dplyr)
library(readr)
library(stringr)


##----------------------------------------------##
## 1.DC as center cell 20um radius Neighboring cells (proportion and counts)
## Proportion of each DC cell type in each FOV, with the denominator being all the cells around it
## 2.The most adjacent cell type of DC, RANK nearestNGB (statistic by FOV)
##----------------------------------------------##
input_dir <- "./Results/neighbor_FOV20um/"
output_dir <- "./Results/neighbor_FOV20um/U20_NGB_byFOV/"  
files <- list.files(input_dir, pattern = "_byFOV_20μm_neighborNGB\\.csv$", full.names = TRUE)
summary_list <- list()
rank_summary_list <- list()
counts_summary_list <- list()

for (file in files) {
  message("Processing: ", basename(file))
  sample_name <- str_extract(basename(file), "^[^_]+")
  data <- read_csv(file, show_col_types = FALSE)
  # Only retain center cell == DC
  dc_data <- data %>% filter(Name == "DC")
  if (nrow(dc_data) == 0) {
    message("Warning! ", sample_name, " NO DC data!")
    next }
  # Calculate the total number of neighboring cells for each DC.
  dc_data <- dc_data %>%
    mutate(total_Cell = `CD4 Tnaive_mem` + `CD4 Th` + `CD8 Teff` + `CD8 Tnaive_mem` + `DC` + `Other`) %>%
    filter(total_Cell > 0)
  
  # === Each FOV:cell type number ===
  fov_counts_summary <- dc_data %>%
    group_by(FOV) %>%
    dplyr::summarise(
      Sample = sample_name,
      Mean_CD4_Tnaive_mem = mean(`CD4 Tnaive_mem`, na.rm = TRUE),
      Mean_CD4_Th = mean(`CD4 Th`, na.rm = TRUE),
      Mean_CD8_Teff = mean(`CD8 Teff`, na.rm = TRUE),
      Mean_CD8_Tnaive_mem = mean(`CD8 Tnaive_mem`, na.rm = TRUE),
      Mean_DC = mean(`DC`, na.rm = TRUE),
      Mean_Other = mean(`Other`, na.rm = TRUE),
      Mean_total_Cell = mean(total_Cell, na.rm = TRUE),
      DC_count = dplyr::n(),
      .groups = 'drop')
  counts_summary_list[[sample_name]] <- fov_counts_summary
  
  # Proportion
  dc_data <- dc_data %>%
    mutate(
      Pro_CD4_Tnaive_mem = `CD4 Tnaive_mem` / total_Cell,
      Pro_CD4_Th = `CD4 Th` / total_Cell,
      Pro_CD8_Teff = `CD8 Teff` / total_Cell,
      Pro_CD8_Tnaive_mem = `CD8 Tnaive_mem` / total_Cell,
      Pro_DC = `DC` / total_Cell,
      Pro_Other = `Other` / total_Cell,
      T_naive = `CD4 Tnaive_mem` + `CD8 Tnaive_mem`,
      Pro_T_naive = T_naive / total_Cell)
  # Group statistics by FOV
  fov_ratio_summary <- dc_data %>%
    group_by(FOV) %>%
    dplyr::summarise(
      Sample = sample_name,
      Mean_Pro_CD4_Tnaive_mem = mean(Pro_CD4_Tnaive_mem, na.rm = TRUE),
      Mean_Pro_CD4_Th = mean(Pro_CD4_Th, na.rm = TRUE),
      Mean_Pro_CD8_Teff = mean(Pro_CD8_Teff, na.rm = TRUE),
      Mean_Pro_CD8_Tnaive_mem = mean(Pro_CD8_Tnaive_mem, na.rm = TRUE),
      Mean_Pro_DC = mean(Pro_DC, na.rm = TRUE),
      Mean_Pro_Other = mean(Pro_Other, na.rm = TRUE),
      Mean_Pro_T_naive = mean(Pro_T_naive, na.rm = TRUE),
      DC_count = dplyr::n(),
      .groups = 'drop')
    
  write_csv(dc_data, file.path(output_dir, paste0(sample_name, "_U20byFOV_neighbor_ratio_perDC.csv")))
  write_csv(fov_ratio_summary, file.path(output_dir, paste0(sample_name, "_U20byFOV_neighbor_ratio.csv")))
  write_csv(fov_counts_summary, file.path(output_dir, paste0(sample_name, "_U20byFOV_neighbor_counts.csv")))
  
  # === DC Nearest neighbor ===
  target_types <- c("CD4 Tnaive_mem", "CD4 Th", "CD8 Teff", "CD8 Tnaive_mem", "DC", "Other")
  dc_nearest <- dc_data %>% filter(nearestNGB %in% target_types)
  if (nrow(dc_nearest) > 0) {

    fov_list <- unique(dc_nearest$FOV)
    fov_rank_results <- list()
    
    for (fov_name in fov_list) {
      fov_subset <- dc_nearest %>% filter(FOV == fov_name)
      
      if (nrow(fov_subset) > 0) {
        tab <- as.data.frame(table(fov_subset$nearestNGB), stringsAsFactors = FALSE)
        colnames(tab) <- c("nearestNGB", "Count")
        tab$Proportion <- tab$Count / sum(tab$Count)
        tab$Sample <- sample_name
        tab$FOV <- fov_name
        tab <- tab[order(-tab$Count), ]
        fov_rank_results[[as.character(fov_name)]] <- tab
      }
    }
    
    if (length(fov_rank_results) > 0) {
      fov_rank_all <- bind_rows(fov_rank_results)
      write_csv(fov_rank_all, file.path(output_dir, paste0(sample_name, "_U20byFOV_nearest_rank.csv")))
      rank_summary_list[[sample_name]] <- fov_rank_all
    }
  }
  
  summary_list[[sample_name]] <- fov_ratio_summary
}

fov_ratio_summary_all <- bind_rows(summary_list)
fov_rank_summary_all <- bind_rows(rank_summary_list)
fov_counts_summary_all <- bind_rows(counts_summary_list)

# Round up cells number
fov_counts_summary_rounded <- fov_counts_summary_all %>%
  mutate(
    Mean_CD4_Tnaive_mem = ceiling(Mean_CD4_Tnaive_mem),
    Mean_CD4_Th = ceiling(Mean_CD4_Th),
    Mean_CD8_Teff = ceiling(Mean_CD8_Teff),
    Mean_CD8_Tnaive_mem = ceiling(Mean_CD8_Tnaive_mem),
    Mean_DC = ceiling(Mean_DC),
    Mean_Other = ceiling(Mean_Other),
    Mean_total_Cell = ceiling(Mean_total_Cell) )
  
write_csv(fov_ratio_summary_all, file.path(output_dir, "U20byFOV_summary_neighbor_ratio_all.csv"))
write_csv(fov_rank_summary_all, file.path(output_dir, "U20byFOV_summary_nearest_rank_all.csv"))
write_csv(fov_counts_summary_all, file.path(output_dir, "U20byFOV_summary_neighbor_counts_all.csv"))
write_csv(fov_counts_summary_rounded, file.path(output_dir, "U20byFOV_summary_neighbor_counts_all_rounded.csv"))

message("All files successfully processed! Output in: ", output_dir)


