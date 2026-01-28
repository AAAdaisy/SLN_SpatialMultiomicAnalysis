##################################################
## Project: SMI 64-plex protein 
## Script purpose: After AtoMx QC, use CELESTA tools annotate cells
## Date: 24-09-13
##################################################
rm(list=ls())
gc()

library(base)
library(CELESTA)
library(Rmixmod)
library(spdep)
library(ggplot2)
library(reshape2)
library(zeallot)
library(dplyr)
library(pheatmap)


##==========================================================================
# Cell annotation:using CELESTA (Cell TypE identification with SpaTiAl information) tool.
# A hierarchical annotation strategy was applied: Layer 1 (non-immune cells and major immune cell types), T-cell layer, and B-cell layer.
# **We sincerely thank Dr.Felicia New for her professional advice on SMI proteome annotation!**
##==========================================================================
listname = c("S9895889","S9904003","S9795736","S9904050","T001848336","S9789310")
datapath <- "E:/SMI/SMI_CELESTA/ReQC_09CELESTA/data/COORD/"

for (sample in listname){
  var_test=read.csv(paste0("./data/mergeQCAto_",sample,"_matrix_spatial.csv"),header=TRUE)
  rownames(var_test) = var_test$X
  var_test = var_test[,-1]
  names(var_test)[names(var_test)=="CenterX_global_px"] <- "X"
  names(var_test)[names(var_test)=="CenterY_global_px"] <- "Y"
  min(var_test$X)
  var_test <- var_test %>% select(c("X","Y"), everything())
  write.csv(var_test,file = paste0(datapath,sample,"_matrix_spatial_COORD.csv"),row.names = FALSE)
}


##==========================================================================
# Layer1 AssignCells
##==========================================================================
listname = c("S9895889","S9904003","S9795736","S9904050","T001848336","S9789310")
Layer1_marker_info <- read.csv("./Docu/anno_marker/CELESTA_Layer1_signature.csv")

for (sample in listname){
  var_test <- read.csv(paste0("./data/COORD/",sample,"_matrix_spatial_COORD.csv"))
  CelestaObjL1 <- CreateCelestaObject(project_title = paste0("SMI64_Layer1_",sample),Layer1_marker_info,var_test)
  CelestaObjL2 <- AssignCells(CelestaObjL1,max_iteration=10,cell_change_threshold=0.01,
                            high_expression_threshold_anchor=rep(0.7,15),
                            low_expression_threshold_anchor=rep(0.9,15),
                            high_expression_threshold_index=rep(0.5,15),
                            low_expression_threshold_index=rep(1,15))
  result <- read.csv(paste0("./SMI64_Layer1_",sample,"_final_cell_type_assignment.csv"))
  data_summary <- as.data.frame(table(result$Final.cell.type))
  colnames(data_summary) <- c("celltype", "count")
  write.csv(data_summary,file = paste0(docuL1,sample,"CelestaObjL2_result_data_summary.csv"))
  ###Visual heatmap
  # Order
  row_order <- c( "Epithelial","Endothelial", "Fibroblast","Adipose","Plasma","Granulocytes",
                  "cDC","pDC","NK","Macrophage_monocyte","Bcell","CD4_Tcell","CD8_Tcell")
  col_order <- c("EpCAM","CD31","Fibronectin","FABP4","CD138","CD45","CD15","CD11c","CD123","CD56","CD68","CD20","CD3","CD4","CD8")
  
  prob <- as.data.frame(CelestaObjL2@marker_exp_prob)
  final_celltype <- as.data.frame(result$Final.cell.type)
  prob_result <- cbind(final_celltype,prob)
  colnames(prob_result)[1] <- "Final.cell.type"
  # Group and calculate means
  Prob_split_result <- split(prob_result, prob_result$Final.cell.type)
  Prob_summary_result <- sapply(Prob_split_result, function(group) {
    colMeans(group[, -1])
  })
  Prob_summary_result <- as.data.frame(t(Prob_summary_result))
  Prob_summary_result <- Prob_summary_result[row_order,col_order]
  # Heatmap
  p2 <- pheatmap(Prob_summary_result,
                 show_rownames = TRUE,
                 show_colnames = TRUE,
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 scale = 'none')
  ggsave(p2,filename = paste0(diplotL1,sample,"heatmap_CELESTA_AnnoProp.png"),dpi =600,height =8,width = 9.5,units = "in",limitsize = FALSE)
  dev.off()
  
}


##==========================================================================
# T-cell layer & B-cell layer AssignCells
##==========================================================================
T_marker_info <- read.csv("./Mark_table/CELESTA_T_signature.csv")
B_marker_info <- read.csv("./Mark_table/CELESTA_B_signature.csv")

listname = c("S9895889","S9904003","S9795736","S9904050","T001848336","S9789310")
for (sample in listname){
  result_1Layer <- read.csv(paste0("./Docu/RESULT/SMI64_Layer1_",sample,"_final_cell_type_assignment.csv"))
  var_test <- read.csv(paste0("./data/COORD/",sample,"_matrix_spatial_COORD.csv"))
  data_summary <- as.data.frame(table(result_1Layer$Final.cell.type))
  colnames(data_summary) <- c("celltype", "count")
  ##----------------------------------------------##
  ## T-cell layer annotate--------                 
  ##----------------------------------------------##
  Tclust_temp <- subset(result_1Layer,Final.cell.type %in% c("CD4_Tcell","CD8_Tcell","Tcell"),select = c("Final.cell.type","X","Y"))
  Tcluster <-  merge(Tclust_temp, var_test, by = c("X", "Y"), all.x = TRUE)
  Tcluster_a <- Tcluster[ , !colnames(Tcluster) %in% "Final.cell.type"]
  CelestaObjT1 <- CreateCelestaObject(project_title = paste0("SMI64_Tlayer_", sample),T_marker_info,Tcluster_a)
  
  CelestaObjT2 <- AssignCells(CelestaObjT1,max_iteration=10,cell_change_threshold=0.01,
                              high_expression_threshold_anchor=c(0.7,0.7,0.6,0.7,0.7,0.65,0.7,0.75,0.8,0.7,0.7,0.7,0.7,0.7),
                              low_expression_threshold_anchor=rep(0.9,14),
                              high_expression_threshold_index=c(0.5,0.5,0.4,0.6,0.6,0.4,0.4,0.6,0.65,0.5,0.45,0.5,0.5,0.5),
                              low_expression_threshold_index=c(0.9,0.9,0.8,rep(0.9,14)))
  
  Tclust_result <- read.csv(paste0("./SMI64_Tlayer_",sample, "_final_cell_type_assignment.csv"))
  
  T_summary <- as.data.frame(table(Tclust_result$Final.cell.type))
  colnames(T_summary) <- c("celltype", "count")
  write.csv(T_summary,file = paste0(docuT,sample,"CelestaObjT2_result_T_summary.csv"))
  ##----------------------------------------------##
  ## B-cell layer annotate--------                 
  ##----------------------------------------------##
  Bclust_temp <- subset(result_1Layer,Final.cell.type %in% c("Bcell"),select = c("Final.cell.type","X","Y"))
  Bcluster <-  merge(Bclust_temp, var_test, by = c("X", "Y"), all.x = TRUE)
  Bcluster_a <- Bcluster[ , !colnames(Bcluster) %in% "Final.cell.type"]
  #B_marker_info <- read.csv("./Mark_table/CELESTA_Bsignature.csv")
  CelestaObjB1 <- CreateCelestaObject(project_title = paste0("SMI64_Blayer_",sample),B_marker_info,Bcluster_a)
  
  CelestaObjB2 <- AssignCells(CelestaObjB1,max_iteration=10,cell_change_threshold=0.01,
                              high_expression_threshold_anchor=rep(0.7,6),
                              low_expression_threshold_anchor=rep(0.9,6),
                              high_expression_threshold_index=rep(0.5,6),
                              low_expression_threshold_index=rep(1,6))
  
  Bclust_result <- read.csv(paste0("./SMI64_Blayer_", sample, "_final_cell_type_assignment.csv"))
  B_summary <- as.data.frame(table(Bclust_result$Final.cell.type))
  colnames(B_summary) <- c("celltype", "count")
  write.csv(B_summary,file = paste0(docuB,sample,"CelestaObjB2_result_B_summary.csv"))

}


##==========================================================================
# Cell annotation Tidy and summary; Visual
##==========================================================================
getwd()
rm(list=ls())
gc()

library(CELESTA)
library(Rmixmod)
library(spdep)
library(ggplot2)
library(reshape2)
library(zeallot)
library(dplyr)
library(pheatmap)
library(tidyr)
library(ggpubr)


listname <- c("S9895889", "S9904003","S9904050", "S9795736","T001848336", "S9789310")
for (id in listname) {
  
  ##----------------------------------------------##
  ### Layer1 annotation process -----
  ##----------------------------------------------##
  COORDid <- read.csv(file = paste0("./data/COORD/", id, "_matrix_spatial_COORD.csv"))
  layer1 <- read.csv(paste0("./Docu/RESULT/SMI64_Layer1_", id, "_final_cell_type_assignment.csv"))
  Final_cell <- subset(layer1, select = c("X", "Y", "Final.cell.type"))
  merged_annocell <- merge(COORDid, Final_cell, by = c("X", "Y"), all.x = TRUE)
  merged_annocell <- merged_annocell %>% select("Final.cell.type", everything())
  
  names(merged_annocell)[names(merged_annocell) == "X"] <- "CenterX_global_px"
  names(merged_annocell)[names(merged_annocell) == "Y"] <- "CenterY_global_px"
  space <- read.csv(paste0("./data/mergeQCAto_", id, "_matrix_spatial.csv"))
  id_df <- subset(space, select = c("X", "CenterX_global_px", "CenterY_global_px"))
  names(id_df)[names(id_df) == "X"] <- "cell_id"
  
  ID_annoRE <- merge(id_df, merged_annocell, by = c("CenterX_global_px", "CenterY_global_px"), all.x = TRUE)
  ID_annoRE$Final.cell.type[which(ID_annoRE$Final.cell.type == "Immune")] <- "Unknown"
  ID_annoRE$Final.cell.type[which(ID_annoRE$Final.cell.type == "Tcell")] <- "Unknown"
  idcell <- ID_annoRE[, c("Final.cell.type", "cell_id")]
  write.csv(idcell, file = paste0(docuS, id, "_SMI64_Layer1_CELLID_Final.csv"))
  ##Annotate celltype summary
  summary_idcell <- as.data.frame(table(idcell$Final.cell.type))
  colnames(summary_idcell) <- c("Final.cell.type", "cell_number")
  write.csv(summary_idcell, file = paste0(docuS, id, "_SMI64_Layer1_SUMMARY_Final.csv"))
  
  ##----------------------------------------------##
  ### T-cell layer annotation process -----
  ##----------------------------------------------##
  Tclu <- read.csv(paste0("./SMI64_Tlayer_", id, "_final_cell_type_assignment.csv"))
  Tclu_Final_cell <- subset(Tclu, select = c("X", "Y", "Final.cell.type"))
  Tclu_merged_annocell <- merge(Tclu_Final_cell, COORDid, by = c("X", "Y"), all.x = TRUE)
  Tclu_merged_annocell <- Tclu_merged_annocell %>% select("Final.cell.type", everything())
  
  names(Tclu_merged_annocell)[names(Tclu_merged_annocell) == "X"] <- "CenterX_global_px"
  names(Tclu_merged_annocell)[names(Tclu_merged_annocell) == "Y"] <- "CenterY_global_px"
  
  Tclu_ID_anno <- merge(Tclu_merged_annocell, id_df, by = c("CenterX_global_px", "CenterY_global_px"), all.x = TRUE)
  Tclu_ID_anno <- Tclu_ID_anno %>% select(c("cell_id", "Final.cell.type", "CenterX_global_px", "CenterY_global_px"), everything())
  Tclu_ID_anno$Final.cell.type[which(Tclu_ID_anno$Final.cell.type == "CD8+_TCF1-")] <- "CD8_Teff"
  Tclu_ID_anno$Final.cell.type[which(Tclu_ID_anno$Final.cell.type == "CD8+_TCF1+_CD45RA+")] <- "CD8_Tnaive_mem"
  Tclu_ID_anno$Final.cell.type[which(Tclu_ID_anno$Final.cell.type == "CD8+_TCF1+_CD45RA-_PD1+")] <- "CD8_Tpex"
  Tclu_ID_anno$Final.cell.type[which(Tclu_ID_anno$Final.cell.type == "CD8+_TCF1+_CD45RA-_PD1-")] <- "CD8_Teff"
  Tclu_ID_anno$Final.cell.type[which(Tclu_ID_anno$Final.cell.type == "CD4+_FOXP3+")] <- "CD4_Treg"
  Tclu_ID_anno$Final.cell.type[which(Tclu_ID_anno$Final.cell.type == "CD4+_FOXP3-_TCF1-")] <- "CD4_Teff"
  Tclu_ID_anno$Final.cell.type[which(Tclu_ID_anno$Final.cell.type == "CD4+_TCF1+_CD45RA+")] <- "CD4_Tnaive_mem"
  Tclu_ID_anno$Final.cell.type[which(Tclu_ID_anno$Final.cell.type == "CD4+_TCF1+_CD45RA-")] <- "CD4_Teff"
  Tclu_ID_anno$Final.cell.type[which(Tclu_ID_anno$Final.cell.type %in% c("CD4+_FOXP3-", "CD4+_FOXP3-_TCF1+", "CD4+_T", "CD8+_T", "CD8+_TCF1+", "CD8+_TCF1+_CD45RA-"))] <- "Unknown"
  
  Tclu_idcell <- Tclu_ID_anno[, c("Final.cell.type", "cell_id")]
  write.csv(Tclu_idcell, file = paste0(docuS, id, "_SMI64_Tlayer_CELLID_Final.csv"))
  
  summary_Tclu_idcell <- as.data.frame(table(Tclu_idcell$Final.cell.type))
  colnames(summary_Tclu_idcell) <- c("Final.cell.type", "cell_number")
  write.csv(summary_Tclu_idcell, file = paste0(docuS, id, "_SMI64_Tlayer_SUMMARY_Final.csv"))
  
  ##----------------------------------------------##
  ### B-cell layer annotation process -----
  ##----------------------------------------------##
  Bclu <- read.csv(paste0("./SMI64_Blayer_", id, "_final_cell_type_assignment.csv"))
  
  Bclu_Final_cell <- subset(Bclu, select = c("X", "Y", "Final.cell.type"))
  Bclu_merged_annocell <- merge(Bclu_Final_cell, COORDid, by = c("X", "Y"), all.x = TRUE)
  Bclu_merged_annocell <- Bclu_merged_annocell %>% select("Final.cell.type", everything())
  
  names(Bclu_merged_annocell)[names(Bclu_merged_annocell) == "X"] <- "CenterX_global_px"
  names(Bclu_merged_annocell)[names(Bclu_merged_annocell) == "Y"] <- "CenterY_global_px"
  
  Bclu_ID_anno <- merge(Bclu_merged_annocell, id_df, by = c("CenterX_global_px", "CenterY_global_px"), all.x = TRUE)
  Bclu_ID_anno <- Bclu_ID_anno %>% select(c("cell_id", "Final.cell.type", "CenterX_global_px", "CenterY_global_px"), everything())
  
  Bclu_idcell <- Bclu_ID_anno[, c("Final.cell.type", "cell_id")]
  write.csv(Bclu_idcell, file = paste0(docuS, id, "_SMI64_Blayer_CELLID_Final.csv"))
  
  summary_Bclu_idcell <- as.data.frame(table(Bclu_idcell$Final.cell.type))
  colnames(summary_Bclu_idcell) <- c("Final.cell.type", "cell_number")
  write.csv(summary_Bclu_idcell, file = paste0(docuS, id, "_SMI64_Blayer_SUMMARY_Final.csv"))
  
}







