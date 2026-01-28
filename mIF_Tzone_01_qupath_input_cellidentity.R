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




