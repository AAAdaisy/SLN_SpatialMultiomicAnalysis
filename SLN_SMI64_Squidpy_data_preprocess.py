# -*- coding: utf-8 -*-
"""
Project: SMI Protein
Script purpose: Try Squidpy,input data:Slide 23R5087/2000/5353/4439SLZA (4 SLNs)
Created on Thu Aug 29 15:30:36 2024
"""

import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
from scipy import stats
import scanpy.external as sce

from pathlib import Path
import squidpy as sq
import seaborn as sns

pd.set_option('display.max_columns',70)
pd.set_option('display.width', 50)
pd.set_option('display.max_colwidth',50)

plopath = "./SMI-64plex/Step0_Python_dataPreprocess/Plot/"
docupath = "./SMI-64plex/Step0_Python_dataPreprocess/Docu/"
prefix = "SMI_Squidpy_"



"""
Squidpy data import
"""
##Download the Nanostring data, unpack and load to AnnData
adata1 = sq.read.nanostring(
    path="E:/SMI/SMI_Squidpy/SMI_Squidpy/Squidpy_input1", 
    counts_file="23R63134988SLZA_exprMat_file.csv",
    meta_file="23R63134988SLZA_metadata_file.csv",
    fov_file="23R63134988SLZA_fov_positions_file.csv" )

adata2 = sq.read.nanostring(
    path="E:/SMI/SMI_Squidpy/SMI_Squidpy/Squidpy_input2",
    counts_file="23R5087200053534439SLZA_exprMat_file.csv",
    meta_file="23R5087200053534439SLZA_metadata_file.csv",
    fov_file="23R5087200053534439SLZA_fov_positions_file.csv" )


"""
Based on the FOV matching with the sample
"""
adata1.obs.set_index("cell_id", inplace=True)
adata2.obs.set_index("cell_id", inplace=True)
adata1.obs_names = [f"s1_{cell}" for cell in adata1.obs_names]
adata2.obs_names = [f"s2_{cell}" for cell in adata2.obs_names]

#Based on the FOV matching with the sample
#####S1:Slide1#######
adata1.obs['sample'] = ''
adata1.obs['fov'] = adata1.obs['fov'].astype(int)
adata1.obs.loc[adata1.obs['fov'].between(1, 171), 'sample'] = 'T001848336'
adata1.obs.loc[adata1.obs['fov'].between(172, 274), 'sample'] = 'S9789310'
adata1.write_h5ad(os.path.join(docupath, prefix +'Raw_Noqc_2LN.h5ad'))

#####S2:Slide2#########
adata2.obs['sample'] = ''
adata2.obs['fov'] = adata2.obs['fov'].astype(int)
adata2.obs.loc[adata2.obs['fov'].between(1, 39), 'sample'] = 'S9904003'
adata2.obs.loc[adata2.obs['fov'].between(40, 55), 'sample'] = 'S9904050'
adata2.obs.loc[adata2.obs['fov'].between(56, 77), 'sample'] = 'S9895889'
adata2.obs.loc[adata2.obs['fov'].between(78, 390), 'sample'] = 'S9795736'
adata2.write_h5ad(os.path.join(docupath, prefix +'Raw_Noqc_4LN.h5ad'))


"""
Standardized spatial coordinates
"""
adata1.obs["CenterX_global_px"] = adata1.obs["CenterX_global_px"] - adata1.obs["CenterX_global_px"].min()
adata1.obs["CenterY_global_px"] = adata1.obs["CenterY_global_px"] - adata1.obs["CenterY_global_px"].min()

adata2.obs["CenterX_global_px"] = adata2.obs["CenterX_global_px"] - adata2.obs["CenterX_global_px"].min()
adata2.obs["CenterY_global_px"] = adata2.obs["CenterY_global_px"] - adata2.obs["CenterY_global_px"].min()

adata1.obs['CenterX_global_px'].max()
adata2.obs["CenterX_global_px"] = adata2.obs["CenterX_global_px"]+ 108333.333 #13mm

adata1.obsm['spatial'][:, 0] = adata1.obsm['spatial'][:, 0] - adata1.obsm['spatial'][:, 0].min()
adata1.obsm['spatial'][:, 1] = adata1.obsm['spatial'][:, 1] - adata1.obsm['spatial'][:, 1].min()

adata2.obsm['spatial'][:, 0] = adata2.obsm['spatial'][:, 0] - adata2.obsm['spatial'][:, 0].min()
adata2.obsm['spatial'][:, 1] = adata2.obsm['spatial'][:, 1] - adata2.obsm['spatial'][:, 1].min()

adata1.obsm['spatial_fov'][:, 0] = adata1.obsm['spatial_fov'][:, 0] - adata1.obsm['spatial_fov'][:, 0].min()
adata1.obsm['spatial_fov'][:, 1] = adata1.obsm['spatial_fov'][:, 1] - adata1.obsm['spatial_fov'][:, 1].min()

adata2.obsm['spatial_fov'][:, 0] = adata2.obsm['spatial_fov'][:, 0] - adata2.obsm['spatial_fov'][:, 0].min()
adata2.obsm['spatial_fov'][:, 1] = adata2.obsm['spatial_fov'][:, 1] - adata2.obsm['spatial_fov'][:, 1].min()

adata2.obsm['spatial_fov'][:, 0] = adata2.obsm['spatial_fov'][:, 0] + 108333.333

adata1.write_h5ad(os.path.join(docupath, prefix +'OriginAxis_adata1_Noqc.h5ad'))
adata2.write_h5ad(os.path.join(docupath, prefix +'OriginAxis_adata2_Noqc.h5ad'))

#merge data
mergeNoqc = adata1.concatenate(adata2, batch_key='sample_id') #3691566
mergeNoqc.write_h5ad(os.path.join(docupath, prefix +'Raw_Merge_6LNdata_Noqc.h5ad'))


'''
QC:nanostring AtoMx platform
After QC,Unreliable cells will be marked: .obs['remove_flagged_cells']==True
'''
mergeNoqc = sc.read_h5ad("E:/SMI/SMI_Squidpy/SMI_Squidpy/Docu/SMI_Squidpy_Raw_Merge_6LNdata_Noqc.h5ad")
mergeQCAto = mergeNoqc[mergeNoqc.obs['remove_flagged_cells']==False] 
mergeQCAto.write_h5ad(os.path.join(docupath,prefix+'mergeQC_AtoMx.h5ad'))




