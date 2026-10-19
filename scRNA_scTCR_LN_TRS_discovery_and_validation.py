#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRNA-seq + scTCR-seq Comprehensive Discovery & Validation Pipeline (Complete Version)
- Subcluster UMAP & Clonal Expansion Tier UMAP with Enhanced Point Size (size=22), Exclude CD4_Tstr & Proliferating_T
- CD4 & CD8 Benchmark Tumor Reactivity Boxplots with Mathematically Centered Jitter Alignment
- Identification of LN-specific Tumor-Reactive Gene (LN-TRG) Signature
- Volcano Plot with High-Impact Journal Style Short Leader Lines & Data-Driven Upregulated Gene Annotation
- Internal Dose-Dependent Gradient Validation in pCR CD8_Tem with Pairwise Statistical Testing vs Large
- Cross-Subcluster Validation of LN-TRG Signature (pCR vs pNR across CD4, CD8 and All Subclusters)
- Exporting Gene Lists, Pairwise Statistical Tables & Source Data for downstream CosMx spatial validation
"""

from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import scanpy as sc
import mudata as md

# ==============================================================================
# 0. 画布与路径全局配置
# ==============================================================================
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'

BASE_DIR = Path("/mnt/data1/0_zq/TDLNprogram/data/LN_scRNAseq_scTCRseq/1_DataAnalysis/1_Output")
INPUT_H5MU = BASE_DIR / "processed_data" / "13_TDLN_cd4cd8_res1_1_multisignature_annotated.h5mu"

DIR_FIGURES = BASE_DIR / "figures" / "16_LN_TRG_discovery_and_validation"
DIR_TABLES = BASE_DIR / "tables"
DIR_SOURCE_DATA = BASE_DIR / "source_data" / "16_LN_TRG_discovery"
DIR_DATA = BASE_DIR / "processed_data"

for d in [DIR_FIGURES, DIR_TABLES, DIR_SOURCE_DATA, DIR_DATA]:
    d.mkdir(parents=True, exist_ok=True)

COLOR_PCR = "#BC3C29FF"
COLOR_PNR = "#0072B5FF"
RESPONSE_PALETTE = {"pCR": COLOR_PCR, "pNR": COLOR_PNR}

TIER_ORDER = ["Single (1)", "Small (2-5)", "Medium (6-20)", "Large (21-100)", "Hyperexpanded (>100)"]
TIER_PALETTE = {
    "Single (1)": "#91BFDB", 
    "Small (2-5)": "#E0F3F8", 
    "Medium (6-20)": "#FEE090",
    "Large (21-100)": "#FC8D59", 
    "Hyperexpanded (>100)": "#D73027",
    "No_TCR": "#E0E0E0"
}

CUSTOM_FINE_TYPE_ORDER = [
    "CD8_Naive_T", "CD8_Tisg", "CD8_Tex", "CD4_Naive_T", "CD4_Tcm_1", "CD4_Tcm_2",
    "CD4_Tstr", "CD4_Tem_1", "CD4_Tem_2", "CD4_Treg_1", "CD4_Treg_2", "CD4_Th17",
    "CD8_Tem", "Proliferating_T"
]

def save_dual_figure(fig_stem_name, fig=None):
    png_path = DIR_FIGURES / f"{fig_stem_name}.png"
    pdf_path = DIR_FIGURES / f"{fig_stem_name}.pdf"
    if fig is not None:
        fig.savefig(png_path, dpi=300, bbox_inches="tight")
        fig.savefig(pdf_path, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.savefig(png_path, dpi=300, bbox_inches="tight")
        plt.savefig(pdf_path, bbox_inches="tight")
        plt.close()
    print(f"[SAVED FIG] {fig_stem_name} (.png & .pdf)")

# ==============================================================================
# 1. 加载多模态数据与元数据对齐
# ==============================================================================
print("=" * 70)
print(f"[INITIALIZATION] Loading Annotated Multi-Omics Data from:\n  {INPUT_H5MU}")
print("=" * 70)

if not INPUT_H5MU.exists():
    raise FileNotFoundError(f"Input file not found: {INPUT_H5MU}")

mdata = md.read_h5mu(str(INPUT_H5MU))
adata = mdata["gex"] if "gex" in mdata.mod else (mdata["rna"] if "rna" in mdata.mod else mdata)

if hasattr(mdata, "obs"):
    for col in mdata.obs.columns:
        if col not in adata.obs.columns:
            adata.obs[col] = mdata.obs[col]

fine_col = "fine_cell_type" if "fine_cell_type" in adata.obs.columns else "gex:fine_cell_type"
group_col = "response" if "response" in adata.obs.columns else ("gex:response" if "gex:response" in adata.obs.columns else "Group")
size_col = "sample_clone_size" if "sample_clone_size" in adata.obs.columns else "airr:clone_id_size"

adata.obs["fine_cell_type"] = adata.obs[fine_col].astype(str)
adata.obs["response"] = adata.obs[group_col].astype(str)
adata.obs["sample_clone_size"] = pd.to_numeric(adata.obs[size_col], errors="coerce").fillna(0)

def assign_clone_tier(s):
    try:
        val = float(s)
    except (ValueError, TypeError):
        return "No_TCR"
    if val == 1:
        return "Single (1)"
    elif 2 <= val <= 5:
        return "Small (2-5)"
    elif 6 <= val <= 20:
        return "Medium (6-20)"
    elif 21 <= val <= 100:
        return "Large (21-100)"
    elif val > 100:
        return "Hyperexpanded (>100)"
    else:
        return "No_TCR"

adata.obs["clone_tier"] = adata.obs["sample_clone_size"].apply(assign_clone_tier)
tier_categories = TIER_ORDER + ["No_TCR"]
adata.obs["clone_tier"] = pd.Categorical(adata.obs["clone_tier"], categories=tier_categories)

print(f"Loaded cells: {adata.shape[0]}, genes: {adata.shape[1]}")

# ==============================================================================
# 2. 绘制 UMAP 全景图 (过滤 CD4_Tstr 与 Proliferating_T，增大点尺寸 size=25)
# ==============================================================================
print("\n>>> Step 2: 绘制亚群 UMAP 与 TCR 克隆梯度 UMAP (过滤 CD4_Tstr 与 Proliferating_T)...")

exclude_subsets = ["CD4_Tstr", "Proliferating_T"]
adata_umap = adata[~adata.obs["fine_cell_type"].isin(exclude_subsets)].copy()

fig, axes = plt.subplots(2, 2, figsize=(15, 13))

for i, grp in enumerate(["pCR", "pNR"]):
    ax = axes[0, i]
    sub_adata = adata_umap[adata_umap.obs["response"] == grp]
    sc.pl.umap(
        sub_adata, color="fine_cell_type", ax=ax, show=False, size=25,
        title=f"T Cell Subclusters ({grp})",
        legend_loc="right margin" if i == 1 else "none",
        legend_fontsize=8.5, frameon=False
    )

for i, grp in enumerate(["pCR", "pNR"]):
    ax = axes[1, i]
    sub_adata = adata_umap[adata_umap.obs["response"] == grp]
    sc.pl.umap(
        sub_adata, color="clone_tier", palette=TIER_PALETTE, ax=ax, show=False, size=25,
        title=f"TCR Clonal Expansion Tiers ({grp})",
        legend_loc="right margin" if i == 1 else "none",
        legend_fontsize=8.5, frameon=False
    )

plt.tight_layout()
save_dual_figure("01_UMAP_Subclusters_and_Clonal_Tiers_split_Response", fig)

# ==============================================================================
# 3. CD4 与 CD8 各亚群公共特征集打分对比 (Box + Jitter)
# ==============================================================================
print("\n>>> Step 3: 绘制 CD4 (CRM, NeoTCR4) 与 CD8 (CRM, NeoTCR8) 肿瘤反应打分对比 ...")

cd4_types = [ct for ct in CUSTOM_FINE_TYPE_ORDER if "CD4" in ct and ct in adata.obs["fine_cell_type"].unique() and ct not in exclude_subsets]
cd8_types = [ct for ct in CUSTOM_FINE_TYPE_ORDER if "CD8" in ct and ct in adata.obs["fine_cell_type"].unique()]

crm_col = "tumor_reactivity_crm_score" if "tumor_reactivity_crm_score" in adata.obs.columns else "gex:tumor_reactivity_score"
neotcr4_col = "neotcr4_score" if "neotcr4_score" in adata.obs.columns else "gex:neotcr4_score"
neotcr8_col = "neotcr8_score" if "neotcr8_score" in adata.obs.columns else "gex:neotcr8_score"

def plot_benchmark_box_jitter_perfect_alignment(sub_types, score_col, score_name, fig_name, export_stats_prefix=None):
    df_sub = adata.obs[adata.obs["fine_cell_type"].isin(sub_types)].copy()
    
    stats_records = []
    p_dict = {}
    for st in sub_types:
        v_cr = df_sub.loc[(df_sub["fine_cell_type"] == st) & (df_sub["response"] == "pCR"), score_col].dropna().values
        v_nr = df_sub.loc[(df_sub["fine_cell_type"] == st) & (df_sub["response"] == "pNR"), score_col].dropna().values
        if len(v_cr) > 0 and len(v_nr) > 0:
            u_stat, p_val = stats.mannwhitneyu(v_cr, v_nr, alternative="two-sided")
            p_dict[st] = p_val
            stats_records.append({
                "Score_Name": score_name,
                "Fine_Cell_Type": st,
                "N_pCR": len(v_cr),
                "N_pNR": len(v_nr),
                "Mean_Score_pCR": round(float(np.mean(v_cr)), 4),
                "Median_Score_pCR": round(float(np.median(v_cr)), 4),
                "Mean_Score_pNR": round(float(np.mean(v_nr)), 4),
                "Median_Score_pNR": round(float(np.median(v_nr)), 4),
                "Delta_Mean(pCR_minus_pNR)": round(float(np.mean(v_cr) - np.mean(v_nr)), 4),
                "Mann_Whitney_U_Stat": u_stat,
                "P_Value": p_val
            })
            
    df_stats = pd.DataFrame(stats_records)
    padj_map = {}
    if len(df_stats) > 0:
        df_stats["P_Adj_FDR"] = multipletests(df_stats["P_Value"], method="fdr_bh")[1]
        padj_map = df_stats.set_index("Fine_Cell_Type")["P_Adj_FDR"].to_dict()
        if export_stats_prefix:
            df_stats.to_csv(DIR_TABLES / f"{export_stats_prefix}_stats.csv", index=False)
            df_stats.to_csv(DIR_SOURCE_DATA / f"SourceData_{export_stats_prefix}.csv", index=False)
    
    fig, ax = plt.subplots(figsize=(max(7.0, len(sub_types) * 1.15 + 1.5), 5.6))
    
    np.random.seed(42)
    offset_map = {"pCR": -0.20, "pNR": 0.20}
    color_map = {"pCR": COLOR_PCR, "pNR": COLOR_PNR}
    box_width = 0.32
    
    for idx, st in enumerate(sub_types):
        for grp in ["pCR", "pNR"]:
            offset = offset_map[grp]
            col = color_map[grp]
            center_x = idx + offset
            vals = df_sub.loc[(df_sub["fine_cell_type"] == st) & (df_sub["response"] == grp), score_col].dropna().values
            
            if len(vals) > 0:
                jitter_x = np.random.normal(center_x, 0.045, size=len(vals))
                jitter_x = np.clip(jitter_x, center_x - box_width/2 * 0.75, center_x + box_width/2 * 0.75)
                ax.scatter(jitter_x, vals, color=col, alpha=0.35, s=14, edgecolors='none', zorder=1)
                
                ax.boxplot(
                    vals, positions=[center_x], widths=box_width,
                    patch_artist=True, showfliers=False, zorder=2,
                    boxprops=dict(facecolor=col, alpha=0.85, edgecolor='black', linewidth=1.1),
                    medianprops=dict(color='black', linewidth=1.8),
                    whiskerprops=dict(color='black', linewidth=1.1),
                    capprops=dict(color='black', linewidth=1.1)
                )
                
    y_max = df_sub[score_col].max()
    y_min = df_sub[score_col].min()
    y_span = y_max - y_min
    
    for i, st in enumerate(sub_types):
        q = padj_map.get(st, np.nan)
        if not pd.isna(q):
            txt = f"***\nq={q:.1e}" if q < 0.001 else (f"**\nq={q:.3f}" if q < 0.01 else (f"*\nq={q:.3f}" if q < 0.05 else f"ns\nq={q:.2f}"))
            ax.text(i, y_max + y_span * 0.05, txt, ha="center", va="bottom", fontsize=8, fontweight="bold")
            
    ax.set_xticks(range(len(sub_types)))
    ax.set_xticklabels(sub_types, rotation=35, ha="right", fontsize=9.5)
    ax.set_xlim(-0.6, len(sub_types) - 0.4)
    ax.set_ylim(y_min - y_span * 0.05, y_max + y_span * 0.25)
    ax.set_title(f"{score_name} across Subclusters (pCR vs pNR)", fontweight="bold", fontsize=12)
    ax.set_xlabel("")
    ax.set_ylabel(score_name, fontweight="bold")
    
    legend_elements = [
        mpl.patches.Patch(facecolor=COLOR_PCR, edgecolor='black', label='pCR', alpha=0.85),
        mpl.patches.Patch(facecolor=COLOR_PNR, edgecolor='black', label='pNR', alpha=0.85)
    ]
    ax.legend(handles=legend_elements, title="Response", frameon=False, loc="upper right")
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    save_dual_figure(fig_name, fig)

if crm_col in adata.obs.columns:
    plot_benchmark_box_jitter_perfect_alignment(cd4_types, crm_col, "CRM 50-gene Score (CD4)", "02a_CD4_CRM_Score_Boxplot", "16_CD4_CRM_Score_pCR_vs_pNR")
    plot_benchmark_box_jitter_perfect_alignment(cd8_types, crm_col, "CRM 50-gene Score (CD8)", "02c_CD8_CRM_Score_Boxplot", "16_CD8_CRM_Score_pCR_vs_pNR")

if neotcr4_col in adata.obs.columns:
    plot_benchmark_box_jitter_perfect_alignment(cd4_types, neotcr4_col, "Science NeoTCR4 Score (CD4)", "02b_CD4_NeoTCR4_Score_Boxplot", "16_CD4_NeoTCR4_Score_pCR_vs_pNR")

if neotcr8_col in adata.obs.columns:
    plot_benchmark_box_jitter_perfect_alignment(cd8_types, neotcr8_col, "Science NeoTCR8 Score (CD8)", "02d_CD8_NeoTCR8_Score_Boxplot", "16_CD8_NeoTCR8_Score_pCR_vs_pNR")

# ==============================================================================
# 4. 挖掘 LN 特异性肿瘤反应特征集 (火山图)
# ==============================================================================
print("\n>>> Step 4: 在 pCR CD8_Tem 中提取 Large vs Single 克隆差异基因 (LN-TRG)...")

pcr_cd8_tem_mask = (adata.obs["fine_cell_type"] == "CD8_Tem") & (adata.obs["response"] == "pCR")
adata_tem_pcr = adata[pcr_cd8_tem_mask].copy()

adata_tem_pcr.obs["deg_group"] = np.where(
    adata_tem_pcr.obs["clone_tier"].isin(["Large (21-100)", "Hyperexpanded (>100)"]), "Large",
    np.where(adata_tem_pcr.obs["clone_tier"] == "Single (1)", "Single", "Other")
)

adata_deg = adata_tem_pcr[adata_tem_pcr.obs["deg_group"].isin(["Large", "Single"])].copy()
adata_deg.obs["deg_group"] = adata_deg.obs["deg_group"].astype("category")

n_large = (adata_deg.obs["deg_group"] == "Large").sum()
n_single = (adata_deg.obs["deg_group"] == "Single").sum()
print(f"  Strict Comparison in pCR CD8_Tem: Large/Hyper clones (n={n_large}) vs Singletons (n={n_single})")

sc.tl.rank_genes_groups(adata_deg, groupby="deg_group", groups=["Large"], reference="Single", method="wilcoxon")
deg_df = sc.get.rank_genes_groups_df(adata_deg, group="Large")

tcr_bcr_pattern = r"^(TRA|TRB|TRG|TRD|IGH|IGK|IGL)[VJC]"


CELL_CYCLE_GENES = {
    "MKI67", "TOP2A", "CENPF", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "CDK1", "PCNA",
    "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "STMN1", "TUBB", "TUBA1B",
    "HMGB1", "HMGB2", "HMGN2", "PCLAF", "TYMS", "TK1", "CCNB1", "CCNB2", "CCNA2"
}

deg_df["is_cell_cycle"] = deg_df["names"].isin(CELL_CYCLE_GENES)
deg_df["is_ribosomal"] = deg_df["names"].str.startswith(("RPS", "RPL", "MT-"))
deg_df["is_tcr_bcr"] = deg_df["names"].str.startswith(("TRA","TRB","TRG","TRD","IGH","IGK","IGL","VJC"))


clean_up_genes = deg_df[
    (deg_df["logfoldchanges"] > 0.35) & 
    (deg_df["pvals_adj"] < 0.05) & 
    (~deg_df["is_cell_cycle"]) & 
    (~deg_df["is_tcr_bcr"]) &
    (~deg_df["is_ribosomal"])
]["names"].tolist()

print(f"  [LN-TRG Extracted] Filtered clean upregulated genes in Large clones: {len(clean_up_genes)} genes")

with open(DIR_TABLES / "LN_specific_Tumor_Reactive_Gene_List.txt", "w") as f:
    for g in clean_up_genes:
        f.write(f"{g}\n")

deg_df.to_csv(DIR_TABLES / "16_pCR_CD8_Tem_Large_vs_Single_DEG_all.csv", index=False)
deg_df[deg_df["names"].isin(clean_up_genes)].to_csv(DIR_TABLES / "16_LN_TRG_Signature_Genes_Filtered.csv", index=False)

fig, ax = plt.subplots(figsize=(7.5, 6.2))

deg_df["logFC_plot"] = deg_df["logfoldchanges"].clip(-5.0, 5.0)
deg_df["neg_log10_padj"] = -np.log10(deg_df["pvals_adj"].clip(lower=1e-50))
deg_df["color_cat"] = "Non-significant"
deg_df.loc[deg_df["names"].isin(clean_up_genes), "color_cat"] = "LN-TRG Genes"

sns.scatterplot(
    data=deg_df, x="logFC_plot", y="neg_log10_padj", hue="color_cat",
    palette={"LN-TRG Genes": COLOR_PCR, "Non-significant": "#BDC3C7"},
    alpha=0.75, s=26, edgecolor="none", ax=ax, zorder=2
)

deg_up_clean = deg_df[deg_df["names"].isin(clean_up_genes)].copy()
deg_up_clean["rank_score"] = deg_up_clean["logfoldchanges"] * deg_up_clean["neg_log10_padj"]
top_up_genes_df = deg_up_clean.sort_values("rank_score", ascending=False).head(12).reset_index(drop=True)

print("\n  [Top Upregulated Genes in Large Clones to Annotate]:")
print(top_up_genes_df[["names", "logfoldchanges", "pvals_adj", "neg_log10_padj"]].to_string(index=False))

for idx, row in top_up_genes_df.iterrows():
    gene_name = row["names"]
    gx = row["logFC_plot"]
    gy = row["neg_log10_padj"]
    
    x_offset = 0.28 + (idx % 2) * 0.18
    y_offset = 0.70 + (idx % 3) * 0.55
    
    ax.annotate(
        gene_name,
        xy=(gx, gy),
        xytext=(gx + x_offset, gy + y_offset),
        fontsize=8.5,
        fontweight="bold",
        color="#111111",
        arrowprops=dict(
            arrowstyle="-",
            color="#333333",
            lw=0.75,
            shrinkA=1,
            shrinkB=2
        ),
        zorder=4
    )

ax.axvline(0.35, color="gray", linestyle="--", lw=0.8)
ax.axhline(-np.log10(0.05), color="gray", linestyle=":", lw=0.8)
ax.set_xlim(-5.2, 5.2)
ax.set_title("Identification of LN-specific Tumor Reactive Signature (LN-TRG)\n(pCR CD8_Tem: Large vs Single Clones)", fontweight="bold", fontsize=11)
ax.set_xlabel("Log2 Fold Change (Large vs Single, Clipped at ±5.0)", fontweight="bold")
ax.set_ylabel("-Log10 (FDR q-value)", fontweight="bold")
ax.legend(frameon=False, loc="upper left")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
save_dual_figure("03_pCR_CD8_Tem_Large_vs_Single_Volcano", fig)

# ==============================================================================
# 5. 内生有效性检验：在 pCR CD8_Tem 各克隆梯度中验证剂量递减与两两差异检验 (vs Large)
# ==============================================================================
print("\n>>> Step 5: 计算 LN-TRG Score 并进行 Single/Small/Medium vs Large 两两差异检验与标注...")

valid_clean_genes = [g for g in clean_up_genes if g in (adata.raw.var_names if adata.raw is not None else adata.var_names)]
sc.tl.score_genes(adata, gene_list=valid_clean_genes, score_name="LN_TRG_score", use_raw=True)

pcr_tem_df = adata.obs[(adata.obs["fine_cell_type"] == "CD8_Tem") & (adata.obs["response"] == "pCR")].copy()
valid_tiers = [t for t in TIER_ORDER if t in pcr_tem_df["clone_tier"].unique()]

# 5.1 各梯度基础描述统计
gradient_stats = []
for t in valid_tiers:
    scores = pcr_tem_df.loc[pcr_tem_df["clone_tier"] == t, "LN_TRG_score"].values
    gradient_stats.append({
        "Clone_Tier": t,
        "Cell_Count": len(scores),
        "Mean_LN_TRG_Score": round(float(np.mean(scores)), 4) if len(scores) else np.nan,
        "Median_LN_TRG_Score": round(float(np.median(scores)), 4) if len(scores) else np.nan,
        "Std": round(float(np.std(scores)), 4) if len(scores) else np.nan
    })
df_grad = pd.DataFrame(gradient_stats)
df_grad.to_csv(DIR_TABLES / "16_LN_TRG_Score_Dose_Dependent_Gradient.csv", index=False)

# 5.2 两两差异检验：Single vs Large, Small vs Large, Medium vs Large
target_tier = "Large (21-100)" if "Large (21-100)" in valid_tiers else ("Hyperexpanded (>100)" if "Hyperexpanded (>100)" in valid_tiers else valid_tiers[-1])
pairwise_comparisons = [
    ("Medium (6-20)", target_tier),
    ("Small (2-5)", target_tier),
    ("Single (1)", target_tier)
]

pairwise_test_records = []
scores_target = pcr_tem_df.loc[pcr_tem_df["clone_tier"] == target_tier, "LN_TRG_score"].dropna().values

for group_comp, group_ref in pairwise_comparisons:
    if group_comp in valid_tiers and group_ref in valid_tiers:
        scores_comp = pcr_tem_df.loc[pcr_tem_df["clone_tier"] == group_comp, "LN_TRG_score"].dropna().values
        if len(scores_comp) > 0 and len(scores_target) > 0:
            u_stat, p_val = stats.mannwhitneyu(scores_target, scores_comp, alternative="two-sided")
            pairwise_test_records.append({
                "Comparison": f"{group_comp} vs {group_ref}",
                "Group_1": group_comp,
                "Group_2_Ref": group_ref,
                "N_Group_1": len(scores_comp),
                "N_Group_2_Ref": len(scores_target),
                "Mean_Group_1": round(float(np.mean(scores_comp)), 4),
                "Mean_Group_2_Ref": round(float(np.mean(scores_target)), 4),
                "Delta_Mean(Ref_minus_Comp)": round(float(np.mean(scores_target) - np.mean(scores_comp)), 4),
                "Mann_Whitney_U_Stat": u_stat,
                "P_Value": p_val
            })

df_pairwise_stats = pd.DataFrame(pairwise_test_records)
if len(df_pairwise_stats) > 0:
    df_pairwise_stats["P_Adj_FDR"] = multipletests(df_pairwise_stats["P_Value"], method="fdr_bh")[1]

print("\n  [Pairwise Significance Tests against Large Clones]:")
print(df_pairwise_stats[["Comparison", "N_Group_1", "N_Group_2_Ref", "Delta_Mean(Ref_minus_Comp)", "P_Value", "P_Adj_FDR"]].to_string(index=False))

df_pairwise_stats.to_csv(DIR_TABLES / "16_LN_TRG_Pairwise_Large_vs_Others_stats.csv", index=False)
df_pairwise_stats.to_csv(DIR_SOURCE_DATA / "SourceData_Fig04_pairwise_gradient_tests.csv", index=False)

# 5.3 绘制剂量效应小提琴 + 箱线图并添加两两显著性标注
fig, ax = plt.subplots(figsize=(6.5, 6.0))

sns.violinplot(
    data=pcr_tem_df, x="clone_tier", y="LN_TRG_score", order=valid_tiers,
    palette=TIER_PALETTE, inner=None, linewidth=1.1, alpha=0.6, ax=ax
)
sns.boxplot(
    data=pcr_tem_df, x="clone_tier", y="LN_TRG_score", order=valid_tiers,
    width=0.2, color="white", showfliers=False, linewidth=1.2, ax=ax,
    boxprops=dict(edgecolor="black"), medianprops=dict(color="#D32F2F", linewidth=1.8)
)
sns.stripplot(
    data=pcr_tem_df, x="clone_tier", y="LN_TRG_score", order=valid_tiers,
    color="#2C3E50", size=2.5, alpha=0.35, jitter=0.2, ax=ax, zorder=1
)

tier_rank_map = {t: idx for idx, t in enumerate(valid_tiers)}
pcr_tem_df["tier_rank"] = pcr_tem_df["clone_tier"].map(tier_rank_map)
rho, p_trend = stats.spearmanr(pcr_tem_df["tier_rank"], pcr_tem_df["LN_TRG_score"])

y_max = pcr_tem_df["LN_TRG_score"].max()
y_min = pcr_tem_df["LN_TRG_score"].min()
y_span = y_max - y_min

bracket_levels = [
    ("Medium (6-20)", target_tier, 0.08),
    ("Small (2-5)", target_tier, 0.20),
    ("Single (1)", target_tier, 0.32)
]

for g_comp, g_ref, h_factor in bracket_levels:
    if g_comp in valid_tiers and g_ref in valid_tiers:
        idx_comp = valid_tiers.index(g_comp)
        idx_ref = valid_tiers.index(g_ref)
        
        row_match = df_pairwise_stats[df_pairwise_stats["Group_1"] == g_comp]
        if len(row_match) > 0:
            pval = row_match["P_Value"].values[0]
            padj = row_match["P_Adj_FDR"].values[0]
            
            if padj < 0.001:
                sig_text = f"*** (q={padj:.1e})"
            elif padj < 0.01:
                sig_text = f"** (q={padj:.3f})"
            elif padj < 0.05:
                sig_text = f"* (q={padj:.3f})"
            else:
                sig_text = f"ns (q={padj:.2f})"
                
            y_bar = y_max + y_span * h_factor
            bar_height = y_span * 0.02
            
            ax.plot(
                [idx_comp, idx_comp, idx_ref, idx_ref],
                [y_bar - bar_height, y_bar, y_bar, y_bar - bar_height],
                color="#222222", linewidth=0.9
            )
            ax.text(
                (idx_comp + idx_ref) / 2, y_bar + y_span * 0.015,
                sig_text, ha="center", va="bottom", fontsize=8.0, fontweight="bold", color="#111111"
            )

ax.set_ylim(y_min - y_span * 0.05, y_max + y_span * 0.45)
ax.set_title(f"Internal Dose-Dependent Validation of LN-TRG Signature\n(Overall Trend Spearman rho = {rho:.3f}, p = {p_trend:.2e})", fontweight="bold", fontsize=11)
ax.set_xlabel("TCR Clonal Expansion Tier in pCR CD8_Tem", fontweight="bold")
ax.set_ylabel("LN-specific Tumor Reactive (LN-TRG) Score", fontweight="bold")
ax.set_xticklabels(ax.get_xticklabels(), rotation=25, ha="right", fontsize=9.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
save_dual_figure("04_LN_TRG_Score_Internal_Gradient_Validation", fig)

# ==============================================================================
# 5.4 新增：在 CD4、CD8 及全亚群中系统评估 LN-TRG Score (pCR vs pNR)
# ==============================================================================
print("\n>>> Step 5.4: 计算并绘制各亚群 LN-TRG Score (pCR vs pNR Box+Jitter 对比)...")

# 1. 绘制 CD4 各亚群 LN-TRG Score
plot_benchmark_box_jitter_perfect_alignment(
    cd4_types, "LN_TRG_score", "LN-specific Tumor Reactive (LN-TRG) Score (CD4)",
    "05a_CD4_LN_TRG_Score_Boxplot", "16_CD4_LN_TRG_Score_pCR_vs_pNR"
)

# 2. 绘制 CD8 各亚群 LN-TRG Score
plot_benchmark_box_jitter_perfect_alignment(
    cd8_types, "LN_TRG_score", "LN-specific Tumor Reactive (LN-TRG) Score (CD8)",
    "05b_CD8_LN_TRG_Score_Boxplot", "16_CD8_LN_TRG_Score_pCR_vs_pNR"
)

# 3. 绘制全亚群全景图
all_available_subtypes = [ct for ct in CUSTOM_FINE_TYPE_ORDER if ct in adata.obs["fine_cell_type"].unique() and ct not in exclude_subsets]
plot_benchmark_box_jitter_perfect_alignment(
    all_available_subtypes, "LN_TRG_score", "LN-specific Tumor Reactive (LN-TRG) Score (All Subclusters)",
    "05c_All_Subclusters_LN_TRG_Score_Boxplot", "16_All_Subclusters_LN_TRG_Score_pCR_vs_pNR"
)

# ==============================================================================
# 6. 导出 Source Data 与数据保存
# ==============================================================================
print("\n>>> Step 6: 导出全量 Source Data 并进行数据保存...")

src_export_cols = [
    "response", "fine_cell_type", "sample_clone_size", "clone_tier",
    crm_col, neotcr4_col, neotcr8_col, "LN_TRG_score"
]
src_valid_cols = [c for c in src_export_cols if c in adata.obs.columns]
adata.obs[src_valid_cols].to_csv(DIR_SOURCE_DATA / "SourceData_Figure_scRNA_scTCR_LN_TRG.csv", index=True)

mdata.mod["gex"] = adata
mdata.update()
mdata.write_h5mu(DIR_DATA / "16_TDLN_cd4cd8_LN_TRG_validated.h5mu")
adata.write_h5ad(DIR_DATA / "16_TDLN_cd4cd8_LN_TRG_validated_gex.h5ad")

print("\n" + "=" * 70)
print("[SUCCESS] Refined scRNA-seq + scTCR-seq Pipeline Completed Successfully!")
print(f"1. LN-TRG Gene List        : {DIR_TABLES / 'LN_specific_Tumor_Reactive_Gene_List.txt'}")
print(f"2. Pairwise Stats Table    : {DIR_TABLES / '16_LN_TRG_Pairwise_Large_vs_Others_stats.csv'}")
print(f"3. Subcluster Stats Tables : {DIR_TABLES / '16_CD8_LN_TRG_Score_pCR_vs_pNR_stats.csv'}")
print(f"4. Figures saved to        : {DIR_FIGURES}")
print(f"5. Source Data saved to    : {DIR_SOURCE_DATA}")
print("=" * 70)
