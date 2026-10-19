#conda activate scotia_env

import numpy as np
import pandas as pd
from scipy.spatial import KDTree
# pip install POT[all]
import ot
from scipy.stats import mannwhitneyu
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import warnings
from ot.bregman import sinkhorn
import sys


#Recptor -> CD4_c2_Naive
#Sender  -> c6_cDC2_CLEC10A


# 抑制sinkhorn的数值警告
warnings.filterwarnings("ignore", category=RuntimeWarning, module="ot.bregman")
warnings.filterwarnings("ignore", category=UserWarning, module="ot.bregman")

# 保留其他重要警告
warnings.filterwarnings("default")  # 其他警告仍会显示

# 加载数据
adata = sc.read_h5ad("./adata_T.h5ad") # sys.argv[0]
adata.obs["cellid"] = adata.obs.index
radius = 100  # 邻域半径(μm)
adata.obs['x_position'] = adata.obs['x_position'].round(4)
adata.obs['y_position'] = adata.obs['y_position'].round(4)

# 加载LR对
lr_df = pd.read_csv('filtered_LRpairs.csv')
lr_pairs = [{'ligand': row['source'], 'receptor': row['target']}
            for _, row in lr_df.iterrows()]

def compute_reg(C, min_samples=3, default_reg=0.1, quantile=25, scale_factor=5):
    valid_distances = C[C > 0].flatten()
    if len(valid_distances) < min_samples:  # 样本量不足时使用默认值
        return default_reg
    return max(default_reg, np.percentile(valid_distances, quantile) / scale_factor)

def compute_reg_min_distance(C, default_reg=0.1, scale_factor=3):
    valid_distances = C[C > 0].flatten()
    if len(valid_distances) == 0:
        return default_reg
    min_distance = np.min(valid_distances)
    return max(default_reg, min_distance / scale_factor)



all_results = pd.DataFrame()
for fov in adata.obs['fov'].unique().tolist():
    print(f"\nProcessing FOV: {fov}")
    #fov = 14
    adata_fov = adata[adata.obs['fov'] == fov].copy()

    # 提取细胞
    cd4_mask = adata_fov.obs['annotation'] == 'CD4_c2_Naive' ### sys.argv[1]
    cDC_mask = adata_fov.obs['annotation'] == 'c6_cDC2_CLEC10A' ### sys.argv[2]

    cd4_coords = adata_fov[cd4_mask].obs[['x_position', 'y_position']].values
    cd4_ids = adata_fov[cd4_mask].obs_names.tolist()
    cDC_coords = adata_fov[cDC_mask].obs[['x_position', 'y_position']].values
    cDC_ids = adata_fov[cDC_mask].obs_names.tolist()

    # 构建CD4的KDTree
    cd4_tree = KDTree(cd4_coords)

    # 创建cDC→CD4邻域映射
    cDC_to_cd4 = {}
    for i, (cDC_id, coord) in enumerate(zip(cDC_ids, cDC_coords)):
        neighbor_indices = cd4_tree.query_ball_point(coord, r=radius)
        if len(neighbor_indices) > 0:
            cDC_to_cd4[cDC_id] = [cd4_ids[idx] for idx in neighbor_indices]
        if i % 5 == 0:
            print(f"Processed {i}/{len(cDC_ids)} cDCs")

    print(f"cDCs with ≥1 CD4 neighbor: {len(cDC_to_cd4)}/{len(cDC_ids)}")

    # 计算每对LR的OT
    for lr_pair in lr_pairs:
        #lr_pair = {"ligand":"IL4","receptor":"IL4R"}
        print(f"  LR pair: {lr_pair['ligand']}-{lr_pair['receptor']}")
        #lr_pair = {"ligand":"IL4","receptor":"IL4R"}
        # 计算所有cDC的OT得分
        ot_scores = {}
        for cDC_id in cDC_to_cd4:
            #cDC_id = "c_1_14_2965"
            #print(f"\nProcessing cell: {cDC_id}")  # 打印当前细胞ID
            #print(f"Neighbor CD4 count: {len(cd4_ids)}")

            cDC_idx = adata_fov.obs_names.get_loc(cDC_id)
            cd4_ids = cDC_to_cd4[cDC_id]
            cd4_indices = [adata_fov.obs_names.get_loc(id) for id in cd4_ids]
            # 提取表达量
            ligand_expr = adata_fov[cDC_idx, lr_pair['ligand']].X[0, 0]  ### cDC配体
            receptor_expr = adata_fov[cd4_indices, lr_pair['receptor']].X.toarray().flatten()  ### CD4受体
            if ligand_expr == 0 or receptor_expr.sum() == 0:
                #for cd4 in cd4_ids:
                #    ot_scores[cDC_id+"-"+cd4] = np.nan
                continue

            # 归一化概率分布
            mu = np.array([ligand_expr])          # 源分布 (1个cDC)
            nu = receptor_expr#.reshape(-1, 1)       # 目标分布 (多个CD4)
            #nu = nu / nu.sum()                      # 归一化
            nu_normalized = nu * (mu.sum() / nu.sum())
            cDCcoord = np.array([adata_fov.obs.loc[cDC_id, 'x_position'],adata_fov.obs.loc[cDC_id, 'y_position']], dtype=float)
            cd4coords = np.array(adata_fov.obs.loc[cd4_ids, ['x_position', 'y_position']], dtype=float)
            C = np.linalg.norm(cDCcoord - cd4coords, axis=1).reshape(1, -1)
            reg = compute_reg_min_distance(C, scale_factor=3)
            #reg = compute_reg(C, min_samples=5, default_reg=0.1)
            pi = ot.sinkhorn(mu, nu_normalized, C, reg=reg)  # μ.T形状(1,1), ν.T形状(1,n_CD4)
            for ind in range(len(cd4_ids)):
                ot_scores[cDC_id+"-"+cd4_ids[ind]] = pi[0,ind]

        temp_df = pd.DataFrame.from_dict(ot_scores, orient='index', columns=['OT_score'])
        temp_df['lr_pair'] = f"{lr_pair['ligand']}_{lr_pair['receptor']}"
        temp_df['fov'] = fov
        #temp_df = temp_df.join(metadata, how='left')
        all_results = pd.concat([all_results, temp_df])

# 保存最终结果
all_results.to_csv('./output/OT_results_cDC2_to_CD4Tn.csv')
print("Analysis completed. Results saved to ./output/OT_results_cDC2_to_CD4Tn.csv")


