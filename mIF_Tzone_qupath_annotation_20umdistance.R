##################################################
## Project: TNBC SLN
## Script purpose: mIF 7-plex Tcell zone 
## Date: 25-10-20
##################################################
rm(list = ls())

library(stringr)
library(readxl)
library(tidyverse)
library(dplyr)
library(plyr)
library(readr)
library(data.table)
library(ggpubr)
library(ggplot2)
library(ggbeeswarm)
library(RColorBrewer)
library(patchwork)

##==========================================================================
# Qupath file batch input & T-cell zone select
##==========================================================================

##----------------------------------------------##
## Qupath file batch input--------  
##----------------------------------------------##
path <- "E:/SLN_mIFqupath/Export"
files <- list.files(path, pattern = "_detec\\.txt$", full.names = TRUE)

for (f in files) {
  name <- gsub("_detec\\.txt$", "", basename(f))
  varname <- paste0("C", name)
  # Add cell_id
  df <- read_tsv(f,
                 locale = locale(encoding = "UTF-8"),
                 col_names = TRUE,
                 trim_ws = TRUE) %>%
    mutate(cell_id = paste0(varname, "_", seq_len(n())))
  assign(varname, df, envir = .GlobalEnv)
  message(paste("loaded", varname, ":", nrow(df), "rows"))
}


##----------------------------------------------##
## Select SLN Tcell-zone data：Parent=="Annotation (SLN_Tregion)"
##----------------------------------------------##
selected_cols <- c("Image", "Name", "Classification", "Parent", "Centroid X µm", "Centroid Y µm", 
                   "Nucleus: Area", "Nucleus: Opal 780 mean", "Cell: Opal 570 mean", 
                   "Cell: Opal 690 mean", "Cell: Opal 480 mean", "Cell: Opal 620 mean", 
                   "Cell: Opal 520 mean", "Nucleus/Cell area ratio", "cell_id")
all_vars <- ls(pattern = "^C\\d+$", envir = .GlobalEnv)

for (varname in all_vars) {
  df <- get(varname, envir = .GlobalEnv)
  df_filtered <- df %>%
    filter(Parent == "Annotation (SLN_Tregion)") %>%
    select(all_of(selected_cols))
  new_varname <- paste0(varname, "_SLN_Tregion")
  assign(new_varname, df_filtered, envir = .GlobalEnv)
  message(paste("Created", new_varname, ":", nrow(df_filtered), "rows"))
}


##==========================================================================
# Cell Identity
##==========================================================================
##----------------------------------------------##
# CD11c+ DC
# CD3+CD8⁺PD1+TCF1- CD8 Teff
# CD3+CD8⁺PD1-TCF1+ CD8 Tnaive/memory
# CD3+CD8-PD1+TCF1- CD4 Th
# CD3+CD8-PD1-TCF1+ CD4 Tnaive/memory
##----------------------------------------------##

##----------------------------------------------##
## based on the threshold add 6 markers Pos/neg --------                 
##----------------------------------------------##
## threshold table
cutoff_df <- read_tsv("F:/Work/SLNs_Work/mIF/mIFR/anno/SLN_mIF_6channel_cutoff.txt")  
marker_mapping <- c(
  "TCF1" = "Nucleus: Opal 780 mean",
  "CD11c" = "Cell: Opal 570 mean",
  "PD1" = "Cell: Opal 690 mean",
  "CD8" = "Cell: Opal 480 mean",
  "CD68" = "Cell: Opal 620 mean",
  "CD3" = "Cell: Opal 520 mean")
## For all "_SLN_Tregion" variables add marker +/- -------------
sln_vars <- ls(pattern = "_SLN_Tregion$", envir = .GlobalEnv)

for (varname in sln_vars) {
  df <- get(varname, envir = .GlobalEnv)
  sample_id <- gsub("_SLN_Tregion$", "", varname)
  cutoff_row <- cutoff_df %>% 
    filter(SampleID == sample_id)
  if (nrow(cutoff_row) == 0) {
    warning(paste("No cutoff found for", sample_id))
    next
  }
  
  for (marker in names(marker_mapping)) {
    opal_col <- marker_mapping[marker]
    cutoff_value <- cutoff_row[[marker]]
    df <- df %>%
      mutate(!!marker := ifelse(.data[[opal_col]] >= cutoff_value, "+", "-"))
  }
  assign(varname, df, envir = .GlobalEnv)
  message(paste("Processed", varname, ":", nrow(df), "cells"))
}


##----------------------------------------------##
## Cell identity based on +/- markers--------                 
##----------------------------------------------##
# Annotate Rules
rules <- read.table("F:/Work/SLNs_Work/mIF/mIFR/anno/SLN_mIF_anno_rule.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
markers <- c("CD11c", "CD3", "CD8", "PD1", "TCF1", "CD68")

match_celltype <- function(cell_row, rules) {
  for (i in 1:nrow(rules)) {
    rule <- rules[i, ]
    match_all <- TRUE
    for (m in markers) {
      val_rule <- rule[[m]]
      if (!is.na(val_rule)) { 
        if (cell_row[[m]] != val_rule) {
          match_all <- FALSE
          break }
      }
    }
    if (match_all) return(rule$CellType)
  }
  return("Other")
}

sln_vars <- ls(pattern = "_SLN_Tregion$")

# Save
out_dir <- "annotated_results"
dir.create(out_dir, showWarnings = FALSE)
for (var_name in sln_vars) {
  message("Processing: ", var_name)
  data <- get(var_name)
  data$Name <- apply(data, 1, match_celltype, rules = rules)
  out_file <- file.path(out_dir, paste0(var_name, "_annotated.csv"))
  write.csv(data, out_file, row.names = FALSE)
  assign(var_name, data, envir = .GlobalEnv)
}


##----------------------------------------------##
## Result summary--------                 
##----------------------------------------------##
library(tibble)

target_vars <- ls(pattern = "_SLN_Tregion$")
cell_types <- c("DC", "CD8 Teff", "CD8 Tnaive_mem", "CD4 Th", "CD4 Tnaive_mem", "Other")
cell_types_safe <- make.names(cell_types)
summary_list <- list()

for (v in target_vars) {
  df <- get(v)
  if (!"Name" %in% colnames(df)) {
    warning(paste("Warning", v, "NO Name cols,skip!"))
    next
  }
  # Count celltype number
  count_table <- table(factor(df$Name, levels = cell_types))
  sample_row <- as.data.frame(as.list(count_table))
  colnames(sample_row) <- make.names(names(sample_row))
  # Add sample name 
  sample_row <- sample_row %>%
    mutate(Sample = v,
           Total = sum(count_table)) %>%
    select(Sample, all_of(cell_types_safe), Total)
  summary_list[[v]] <- sample_row
}

summary_table <- bind_rows(summary_list)
colnames(summary_table) <- c("Sample", cell_types, "Total")
write.csv(summary_table, "annotated_results/SLN_mIF_Tregion_summary_Annota.csv", row.names = FALSE)



##################################################
## Project: TNBC SLN
## Script purpose: mIF LN Tcell-zone center cell Distence r=20um
## Date: 25-10-21
##################################################
rm(list = ls())
gc()
library(dplyr)
library(tidyr)
library(stringr)
library(plyr)
library(RANN)

##==========================================================================
# Cell types within the radius range of the central cell
# Build a KD-tree and find neighbors within a radius of r�?
# Faster than sqrt
##==========================================================================

##----------------------------------------------##
##  Per FOV,center cell 20um radius--------                 
##----------------------------------------------##
r <- 20  #Searching radius
output_dir <- "./Results/neighbor_FOV20um/"  
input_dir <- "./Results/with_FOV/" 


if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

file_list <- list.files(path = input_dir, pattern = ".*_with_FOV\\.csv$", full.names = TRUE)

for (i in seq_along(file_list)) {
  file_path <- file_list[i]
  file_name <- basename(file_path)
  sample_id <- sub("_SLN_Tregion_annotated_with_FOV\\.csv$", "", file_name)
  
  cat("========================================")
  cat("Processing NO.", i, "file:", file_name)
  cat("Sample ID:", sample_id)
  cat("========================================")
  
  tryCatch({
    data <- read.csv(file_path, check.names = FALSE)
    data2 <- data %>%
      group_by(Name) %>%
      mutate(Suf = seq_along(Name)) %>%  # or seq_len(n())
      mutate(CellIndex = paste(Name, Suf, sep = "_")) %>%
      select(-Suf) %>% 
      ungroup()
    
    all_fov_results <- list()
    fov_list <- unique(data2$FOV)
    
    for (fov_idx in seq_along(fov_list)) {
      current_fov <- fov_list[fov_idx]
      cat("  Processing FOV [", fov_idx, "/", length(fov_list), "]: ", current_fov, "\n", sep = "")
      fov_data <- data2 %>% filter(FOV == current_fov)
      if (nrow(fov_data) == 0) {
        cat(" Warning: FOV", current_fov, "has no cells, skipping...\n")
        next } 
      cat(" Cells in this FOV:", nrow(fov_data), "\n")
      
      # note: R will convert spaces into dots
      cell_df <- fov_data %>%
        select(X = Centroid.X.µm, Y = Centroid.Y.µm, CellIndex) %>%
        remove_rownames() %>%
        column_to_rownames("CellIndex") %>%
        as.data.frame(check.names = FALSE)
      
      ##----------------------------------------------##
      # Build KD-tree and find the neighbors within a radius of r
      ##----------------------------------------------##
      coords <- as.matrix(cell_df[, c("X", "Y")])
      max_k <- min(nrow(coords), 1000)  # Limit maximum neighbor
      nn_res <- RANN::nn2(data = coords, query = coords, k = max_k, radius = r)
      neighbor_list <- lapply(seq_len(nrow(coords)), function(j) {
        tryCatch({
          idx_all <- nn_res$nn.idx[j, ] 
          dist_all <- nn_res$nn.dists[j, ]
          
          valid_pos <- which(!is.na(dist_all) & dist_all > 0 & dist_all <= r & idx_all > 0)
          
          if (length(valid_pos) == 0) {
            return(data.frame(CellIndex = rownames(cell_df)[j],
                              nearestNGB = NA,
                              stringsAsFactors = FALSE)) }
          # Neighbor cells info.
          neighbor_idx <- idx_all[valid_pos]
          neighbor_cells <- rownames(cell_df)[neighbor_idx]
          neighbor_dists <- dist_all[valid_pos]
          
          tmp_ngb <- fov_data[match(neighbor_cells, fov_data$CellIndex), , drop = FALSE]
          tmp_ngb$dist <- neighbor_dists
          tmp_ngb <- tmp_ngb %>% filter(!is.na(CellIndex))
          
          if (nrow(tmp_ngb) == 0) {
            return(data.frame(CellIndex = rownames(cell_df)[j],
                              nearestNGB = NA,
                              stringsAsFactors = FALSE)) }
          # Closest cell type
          min_dist <- min(tmp_ngb$dist, na.rm = TRUE)
          nearest_types <- unique(tmp_ngb$Name[which(tmp_ngb$dist == min_dist)])
          nearest_ct <- paste0(nearest_types, collapse = ",")
          # Count the number of cells of each type
          ct_table <- table(tmp_ngb$Name)
          ct_counts <- data.frame(CellIndex = rownames(cell_df)[j],nearestNGB = nearest_ct,
                                  stringsAsFactors = FALSE)
          
          for (cell_type in names(ct_table)) {
            ct_counts[[cell_type]] <- as.integer(ct_table[cell_type]) }
          return(ct_counts)
          
        }, error = function(e) {
          return(data.frame(CellIndex = rownames(cell_df)[j],nearestNGB = NA,stringsAsFactors = FALSE))
        })
      })
      
      # Clean data
      valid_list <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0, neighbor_list)
      
      if (length(valid_list) > 0) {
        mid.data <- dplyr::bind_rows(valid_list)
        mid.data$nearestNGB <- as.character(mid.data$nearestNGB)
        cell_type_cols <- setdiff(names(mid.data), c("CellIndex", "nearestNGB"))
        for (col in cell_type_cols) {
          if (col %in% names(mid.data)) {
            mid.data[[col]] <- as.integer(mid.data[[col]])
            mid.data[[col]][is.na(mid.data[[col]])] <- 0L
          }
        }
        fov_result <- merge(fov_data, mid.data, by = "CellIndex", all.x = TRUE)
        # Find all celltype colum
        all_cell_types <- c("DC", "CD8 Teff", "CD8 Tnaive_mem", "CD4 Th", "CD4 Tnaive_mem","Other")
        for (ct in all_cell_types) {
          if (ct %in% names(fov_result)) {
            fov_result[[ct]][is.na(fov_result[[ct]])] <- 0L }
        }
        
        if ("nearestNGB" %in% names(fov_result)) {
          fov_result$nearestNGB[is.na(fov_result$nearestNGB)] <- "None" }
        
        all_fov_results[[current_fov]] <- fov_result
        cat("Success! Neighbors calculated for", nrow(fov_result), "cells\n")
      } else {
        cat("Warning: No valid neighbors found\n")
      }
    }
    
    ##----------------------------------------------##
    # 4. Combing all Fovs results
    ##----------------------------------------------##
    if (length(all_fov_results) > 0) {
      F.D <- dplyr::bind_rows(all_fov_results)
      all_cell_types <- c("DC", "CD8 Teff", "CD8 Tnaive_mem", "CD4 Th", "CD4 Tnaive_mem","Other")
      for (ct in all_cell_types) {
        if (ct %in% names(F.D)) {
          F.D[[ct]][is.na(F.D[[ct]])] <- 0L
        }
      }
      # Process "NA" to "None"
      if ("nearestNGB" %in% names(F.D)) {
        F.D$nearestNGB <- as.character(F.D$nearestNGB)
        F.D$nearestNGB[is.na(F.D$nearestNGB)] <- "None"
      }
      # Seve
      output_file <- paste0(output_dir, sample_id, "_byFOV_20um_neighborNGB.csv")
      write.csv(F.D, file = output_file, row.names = FALSE)
      
      cat("SUCCESS! Saved to:", output_file, "\n")
      cat("Total FOVs processed:", length(all_fov_results), "\n")
    } else {
      cat("Warning: No results generated for this sample\n")
    }
    
  }, error = function(e) {
    cat("ERROR processing", file_name, ":", conditionMessage(e), "\n")
  })
}

cat("========================================")
cat("All files processed,20um DC neighbor!")
cat("========================================")

