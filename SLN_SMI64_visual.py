# -*- coding: utf-8 -*-
"""
Project: SMI protein
Script purpose: S1,S2 data match cells annotation
Created on Tue Oct  8 15:45:18 2024
"""
import os
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

plopath = "./SMI-64plex/Step2_Python_visual/Plot/"
docupath = "./SMI-64plex/Step2_Python_visual/Docu/"
prefix = "SMI_Squidpy_"


"""
Slide1(2LNs)
"""
s1_cell_id = pd.read_csv('./SMI-64plex/Step1_R_cellAnnotate/Result/SMI64_S1_CELESTA_anno_cellid.csv', index=False)
s1_cell_id.set_index('cell_id', inplace=True)
#ln2(Slides NO.1)
ln2_raw = sc.read_h5ad("E:/SMI/SMI_Squidpy/SMI_Squidpy/Docu/SMI_Squidpy_Raw_Noqc_2LN.h5ad")
ln2_raw.obs.set_index("cell", inplace=True)
ln2_raw.raw = ln2_raw.copy() 
ln2_raw.layers["counts"] = ln2_raw.X.copy()
ln2_raw.obs_names = [f"s1_{cell}" for cell in ln2_raw.obs_names]

#mapping CELESTA annotation

ln2_raw.obs['Final.cell.type'] = ln2_raw.obs.index.map(s1_cell_id['Final.cell.type'])
ln2_raw.write_h5ad(os.path.join(docupath, prefix +'S1_ln2_annotation_forSgement.h5ad'))


sc.pp.normalize_total(ln2_raw, inplace=True)
sc.pp.log1p(ln2_raw)
sc.pp.pca(ln2_raw)
sc.pp.neighbors(ln2_raw)
sc.tl.umap(ln2_raw)
sc.tl.leiden(ln2_raw)

color_dict = {
  'Epithelial': '#993366',
  'Endothelial': '#cabbe9',  
  'Fibroblast': '#CD3278', 
  'Adipose': '#CC9999', 
  'Plasma': '#FFD700', 
  'Granulocytes': '#996600',
  'cDC': '#6699CC', 
  'pDC':'#B4CDCD',
  'NK':'#FF8C00',
  'Macrophage_monocyte':'#EE7AE9', 
  'Immature_B':'#FFF8DC',
  'Naïve_B':'#FFEC8B',
  'None_class_switched_memoryB':'#CDBE70',
  'Class_switched_memoryB':'#CDAD00',
  'CD4_Tnaive_mem':'#99CC99',
  'CD4_Teff':'#669933',
  'CD4_Treg':'#99CC00',
  'CD8_Tnaive_mem':'#66CCCC',
  'CD8_Teff':'#009999',
  'CD8_Tpex':'#336666',
  'Unknown':'#B5B5B5'
 }


#spatial_segment
ln2_raw.obs["fov"] = ln2_raw.obs["fov"].astype('str')
ln2_raw.obs["library_id"] = ln2_raw.obs["fov"]
ln2_raw.obs["library_id"] = ln2_raw.obs["library_id"].astype('str')
ln2_raw.obs['Final.cell.type_PDL1'].astype('category')

#spatial_segment per FOVs
for library_id in range(1, 3):
    sq.pl.spatial_segment(
        ln2_raw,
        color="Final.cell.type_PDL1",
        library_key="fov",
        library_id=str(library_id),
        seg_cell_id="cell_ID",
        seg_outline = True, 
        img=False,
        palette = color_dict,
        size=120,
        title = f"slide1 No.{library_id} FOV")
    plt.savefig(os.path.join(plopath, prefix + f'Final.cell.type_PDL1_fov{library_id}_segment.jpg'), bbox_inches='tight', dpi=600)
    plt.close()


"""
Slide2(4LNs)
"""
prefix = "SMI_Squidpy_Slide2_"

s2_cell_id = pd.read_csv('./SMI-64plex/Step1_R_cellAnnotate/Result/SMI64_S2_CELESTA_anno_cellid.csv', index=False)
s2_cell_id.set_index('cell_id', inplace=True)

#ln4(Slides NO.2)
ln4_raw = sc.read_h5ad("E:/SMI/SMI_Squidpy/SMI_Squidpy/Docu/SMI_Squidpy_Raw_Noqc_4LN.h5ad")
ln4_raw.obs.set_index("cell", inplace=True)
ln4_raw.raw = ln4_raw.copy() 
ln4_raw.layers["counts"] = ln4_raw.X.copy()
ln4_raw.obs_names = [f"s2_{cell}" for cell in ln4_raw.obs_names]

#mapping CELESTA annotation
ln4_raw.obs['Final.cell.type'] = ln4_raw.obs.index.map(s2_cell_id['Final.cell.type'])
ln4_raw.obs['Final.cell.type_PDL1'] = ln4_raw.obs.index.map(s2_cell_id['Final.cell.type_PDL1'])
ln4_raw.write_h5ad(os.path.join(docupath, prefix +'ln4_annotation_forSgement.h5ad'))

sc.pp.normalize_total(ln4_raw, inplace=True)
sc.pp.log1p(ln4_raw)
sc.pp.pca(ln4_raw)
sc.pp.neighbors(ln4_raw)
sc.tl.umap(ln4_raw)
sc.tl.leiden(ln4_raw)

sc.pl.umap(ln4_raw, 
           color=["Final.cell.type_PDL1"], 
           palette = color_dict
           )# legend_loc='on data' 
plt.savefig(os.path.join(plopath, prefix + 'umap_Final.cell.type.png'), bbox_inches='tight', dpi=600)
plt.close()

##spatial_segment
ln4_raw.obs["fov"] = ln4_raw.obs["fov"].astype('str')
ln4_raw.obs["library_id"] = ln4_raw.obs["fov"]
ln4_raw.obs["library_id"] = ln4_raw.obs["library_id"].astype('str')
ln4_raw.obs['Final.cell.type_PDL1'].astype('category')

#spatial_segment per FOVs
plopath2 = "./SMI-64plex/Step2_Python_visual/"
for library_id in range(1, 391): 
    sq.pl.spatial_segment(
        ln4_raw,
        color="Final.cell.type_PDL1",
        library_key="fov",
        library_id=str(library_id), 
        seg_cell_id="cell_ID",
        seg_outline = True, 
        img=False,
        #palette = color_dict,
        size=120,
        title = f"slide2 No.{library_id} FOV",
        legend_loc = None ) 
                 
    plt.savefig(os.path.join(plopath2, prefix + f'Final.cell.type_fov{library_id}_segment.jpg'), bbox_inches='tight', dpi=600)
    plt.close()










