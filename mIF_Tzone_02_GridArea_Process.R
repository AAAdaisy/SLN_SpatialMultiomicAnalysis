##################################################
## Project: TNBC SLN
## Script purpose:Rasterized SLN T-regions and visualized
## Date: 25-10-29
##################################################

rm(list = ls())
gc()
library(tidyverse)
library(viridis)
library(ggplot2)
library(scales)


##==========================================================================
# T-cell region rasterization
##==========================================================================

##----------------------------------------------##
# Griding function: Automatically identifies separate Tregions (some LNs have more than one T-region) 
#                 Allocates FOV and visualizes it.
# min_distance: The minimum distance threshold (um) for determining whether two T-regions are separated from each other             
##----------------------------------------------##
assign_fov <- function(df, sample_id, grid_size = 200, min_distance = 500) {
  # 1. Use hierarchical clustering to identify Tregions that are spatially separated
  coords <- df %>% select(Centroid.X.µm, Centroid.Y.µm)
  # Calculate the center point of each grid, reduce calculation.
  df_grid <- df %>%
    mutate(grid_x = floor(Centroid.X.µm / grid_size),grid_y = floor(Centroid.Y.µm / grid_size)) %>%
    group_by(grid_x, grid_y) %>%
    summarise(center_x = mean(Centroid.X.µm),center_y = mean(Centroid.Y.µm),.groups = 'drop')
  # The low number of grid is treated as one area.
  if (nrow(df_grid) <= 2) {
    df <- df %>%
      mutate(grid_x = floor(Centroid.X.µm / grid_size),grid_y = floor(Centroid.Y.µm / grid_size),tregion_id = 1)
  } else {
    # Compute the inter-grid distance matrix.
    dist_matrix <- dist(df_grid %>% select(center_x, center_y))
    # **Crucial:Perform hierarchical clustering to group grids belonging to the same area
    hc <- hclust(dist_matrix, method = "single")
    df_grid$cluster <- cutree(hc, h = min_distance)
    df <- df %>%
      mutate(grid_x = floor(Centroid.X.µm / grid_size),grid_y = floor(Centroid.Y.µm / grid_size)) %>%
      left_join(df_grid %>% select(grid_x, grid_y, cluster),by = c("grid_x", "grid_y")) %>%
      rename(tregion_id = cluster)  
  }
  
  # 2. Each T-region was separately rasterized.
  df_result <- df %>%
    group_by(tregion_id) %>%
    mutate(
      # floor() :coordinates assigned based on the half-open interval [a, b) 
      grid_x_raw = floor(Centroid.X.µm / grid_size),
      grid_y_raw = floor(Centroid.Y.µm / grid_size),
      grid_x_norm = grid_x_raw - min(grid_x_raw),
      grid_y_norm = grid_y_raw - min(grid_y_raw),
      # FOV numbers were calculated using row-major ordering.
      max_x = max(grid_x_norm) + 1,
      fov_num = grid_y_norm * max_x + grid_x_norm + 1,
      # FOV labels
      FOV = paste0(sample_id, "_TR", tregion_id, "_FOV", fov_num)) %>% 
      ungroup() %>%
      select(-grid_x, -grid_y, -grid_x_raw, -grid_y_raw, -grid_x_norm, -grid_y_norm, -max_x, -fov_num, -tregion_id)
  return(df_result)
}


##----------------------------------------------##
###T cell zone visualization ----------
##----------------------------------------------##
plot_fov_distribution <- function(df_fov, sample_id, output_path) {
  n_fov <- n_distinct(df_fov$FOV)
  base_turbo <- viridis(30, option = "turbo") # color template
  #if FOV > 30,add random colors
  if (n_fov <= 30) {
    final_colors <- base_turbo[1:n_fov]
  } else {
    repeats <- ceiling(n_fov / 30)
    final_colors <- rep(base_turbo, repeats)[1:n_fov]
    set.seed(123)
    final_colors <- sample(final_colors)
  }
  
  p <- ggplot(df_fov, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color = FOV)) +
    geom_point(size = 0.4, alpha = 0.6) +
    theme_classic() +
    theme(legend.position = "none",panel.grid = element_blank()) +
    labs(title = paste("FOV distribution -", sample_id),x = "Centroid X (µm)",y = "Centroid Y (µm)") +
    coord_fixed() + scale_color_manual(values = final_colors)
  ggsave(filename = file.path(output_path, paste0(sample_id, "_FOV_plot.pdf")),plot = p,width = 12,height = 10)
  
  return(p)
}


##----------------------------------------------##
## Batch processing --------                 
##----------------------------------------------##
input_dir <- "./Results/annotated_results"
output_table_dir <- "./Results/allFOV_table"
output_picture_dir <- "./Results/allFOV_picture"
dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_picture_dir, recursive = TRUE, showWarnings = FALSE)

csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
# Processing-----
for (i in seq_along(csv_files)) {
  file_path <- csv_files[i]
  file_name <- basename(file_path)
  cat("--------------------------------------------------")
  cat("[", i, "/", length(csv_files), "] Processing:", file_name, "\n")
  
  tryCatch({
    df <- read.csv(file_path, header = TRUE)
    sample_id <- str_extract(df$cell_id[1], "^C\\d+")
    if (is.na(sample_id)) {
      cat("Warning!Cannot choose sample ID,Jump!\n\n")
      next }
    cat("Sample ID:", sample_id, "\n")
    df_fov <- assign_fov(df, sample_id, grid_size = 200, min_distance = 500)
    # Save
    output_csv <- file.path(output_table_dir,str_replace(file_name, "\\.csv$", "_with_FOV.csv")) 
    write.csv(df_fov,output_csv,row.names = FALSE) 
    plot_fov_distribution(df_fov, sample_id, output_picture_dir)
    cat("Picture saved!", paste0(sample_id, "_FOV_plot.pdf"), "\n")

    cat("*****Processing Finish*****")
  }, error = function(e) {
    cat("FALSE!", conditionMessage(e))
  })
}




