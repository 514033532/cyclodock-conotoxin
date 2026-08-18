#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
整合大图：左图 – 食性小提琴图（带统计标注）→ 标注 (A)
        右图 – 物种排名条形图（高对比度颜色）→ 标注 (B)
输出：一张大图 PNG / PDF / SVG，无遮挡
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
from scikit_posthocs import posthoc_conover
from matplotlib.patches import Patch

# ===================== 路径配置 =====================
BASE_DIR = '/home/zhh2/CONO-PEPTIDE2026/2026result/jianzhan'
OUTPUT_DIR = os.path.join(BASE_DIR, 'combined_figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)

VIOLIN_FILE = os.path.join(BASE_DIR, 'merged_with_target_info2Organism_diet.csv')
RANK_FILE = os.path.join(BASE_DIR, 'fenlei/fig4/species_score_diet.csv')

COL_DIET = 'Organism_diet'
COL_SCORE = 'HADDOCK_score_[a.u.]'

# ================== 1. 读取并准备数据 ==================
df_violin = pd.read_csv(VIOLIN_FILE, sep='\t', encoding='utf-8-sig')
df_violin.columns = df_violin.columns.str.strip()
df_violin = df_violin[[COL_DIET, COL_SCORE]].dropna()
df_violin[COL_DIET] = df_violin[COL_DIET].str.strip()
df_violin = df_violin[df_violin[COL_DIET] != 'unknown']

df_rank = pd.read_csv(RANK_FILE)
df_rank = df_rank.dropna(subset=['mean_score', COL_DIET])
df_rank[COL_DIET] = df_rank[COL_DIET].str.strip()
df_rank = df_rank[df_rank[COL_DIET] != 'unknown']
df_rank_sorted = df_rank.sort_values('mean_score', ascending=True).reset_index(drop=True)
df_rank_sorted['rank'] = range(1, len(df_rank_sorted)+1)

# ================== 2. 统计检验（小提琴部分）==================
groups = df_violin.groupby(COL_DIET)
group_names = df_violin[COL_DIET].unique()
kw_groups = [g[COL_SCORE].values for _, g in groups]
h_stat, p_kw = kruskal(*kw_groups)
posthoc = posthoc_conover(df_violin, val_col=COL_SCORE, group_col=COL_DIET, p_adjust='fdr_bh')

group_order_sorted = sorted(group_names, key=lambda x: df_violin[df_violin[COL_DIET]==x][COL_SCORE].mean())
pos_map = {g:i for i,g in enumerate(group_order_sorted)}

pair_list = []
alpha = 0.05
for i,g1 in enumerate(group_order_sorted):
    for g2 in group_order_sorted[i+1:]:
        p = posthoc.loc[g1, g2]
        if p < alpha:
            if p < 0.001: mark = "***"
            elif p < 0.01: mark = "**"
            else: mark = "*"
            pair_list.append({"x1":pos_map[g1], "x2":pos_map[g2], "p":p, "star":mark})

# ================== 3. 颜色映射（排名图）==================
palette = {
    'piscivorous': '#E41A1C',
    'molluscivorous': '#377EB8',
    'vermivorous': '#4DAF4A',
}
diets = df_rank_sorted[COL_DIET].unique()
extra_colors = ['#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF']
next_color = 0
for d in diets:
    if d not in palette:
        palette[d] = extra_colors[next_color % len(extra_colors)]
        next_color += 1

# ================== 4. 绘制大图（1行2列）==================
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 0.8

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 7))  # 宽度增加，留白更多

# ---------- 左图：小提琴图 ----------
order = group_order_sorted
violin_parts = ax1.violinplot(
    [df_violin[df_violin[COL_DIET]==g][COL_SCORE].values for g in order],
    positions=list(range(len(order))),
    widths=0.7,
    showmeans=False, showmedians=False
)
palette_violin = sns.color_palette("Set2", n_colors=len(order))
for i, pc in enumerate(violin_parts["bodies"]):
    pc.set_facecolor(palette_violin[i])
    pc.set_alpha(0.7)
    pc.set_edgecolor("white")
    pc.set_linewidth(0)

ax1.boxplot(
    [df_violin[df_violin[COL_DIET]==g][COL_SCORE].values for g in order],
    positions=list(range(len(order))),
    widths=0.15,
    patch_artist=True,
    showmeans=False,
    boxprops={"facecolor":"white","edgecolor":"black","linewidth":0.8},
    whiskerprops={"color":"black","linewidth":0.8},
    capprops={"color":"black","linewidth":0.8},
    medianprops={"color":"black","linewidth":1.2}
)

for idx,g in enumerate(order):
    arr = df_violin[df_violin[COL_DIET]==g][COL_SCORE].values
    x_jitter = np.random.normal(idx, scale=0.04, size=len(arr))
    ax1.scatter(x_jitter, arr, s=6, c="black", alpha=0.15, linewidth=0)

mean_vals = df_violin.groupby(COL_DIET)[COL_SCORE].mean().reindex(order)
ax1.scatter(range(len(order)), mean_vals, marker="D", s=30, color="#D62728", zorder=10)

y_min = df_violin[COL_SCORE].min()
y_max = df_violin[COL_SCORE].max()
y_range = y_max - y_min
base_y = y_max + 0.08 * y_range
step_y = 0.06 * y_range
for idx, pair in enumerate(pair_list):
    x1, x2 = pair["x1"], pair["x2"]
    current_y = base_y + idx * step_y
    ax1.plot([x1, x2], [current_y, current_y], lw=1.2, c="black")
    ax1.plot([x1, x1], [current_y-0.02*y_range, current_y], lw=1.2, c="black")
    ax1.plot([x2, x2], [current_y-0.02*y_range, current_y], lw=1.2, c="black")
    ax1.text((x1+x2)/2, current_y + 0.01*y_range, pair["star"], ha="center", va="bottom",
             fontsize=10, fontweight="bold")

ax1.set_xticks(range(len(order)))
ax1.set_xticklabels([g.capitalize() for g in order], fontsize=8)
ax1.set_ylabel("HADDOCK score (a.u.)", fontsize=10, fontweight='bold')
ax1.set_xlabel("Dietary type", fontsize=10, fontweight='bold')
ax1.set_ylim(y_min - 0.5, base_y + len(pair_list)*step_y + 0.12*y_range)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.legend(["Mean"], loc="upper right", frameon=False, fontsize=8)

# ---------- 添加 (A) 标签（左上角，在数据坐标外） ----------
ax1.text(-0.18, 1.02, '(A)', transform=ax1.transAxes,
         fontsize=14, fontweight='bold', va='bottom', ha='left')

# ---------- 右图：物种排名条形图 ----------
bars = ax2.barh(df_rank_sorted['rank'], df_rank_sorted['mean_score'],
                color=[palette[d] for d in df_rank_sorted[COL_DIET]],
                edgecolor='black', linewidth=0.5, height=0.7)

for i, (idx, row) in enumerate(df_rank_sorted.iterrows()):
    score = row['mean_score']
    if score < -25:
        ha = 'left'; x_pos = score + 0.3
    else:
        ha = 'right'; x_pos = score - 0.3
    ax2.text(x_pos, i+1, f'{score:.2f}', va='center', ha=ha, fontsize=6, color='black')

ax2.set_yticks(df_rank_sorted['rank'])
ax2.set_yticklabels(df_rank_sorted['Species'], fontsize=7)
ax2.invert_yaxis()
ax2.set_xlabel('Average HADDOCK Score (a.u.)', fontsize=10, fontweight='bold')
ax2.set_ylabel('Species', fontsize=10, fontweight='bold')
ax2.grid(axis='x', linestyle='--', alpha=0.3)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# ---------- 添加 (B) 标签（左上角） ----------
ax2.text(-0.18, 1.02, '(B)', transform=ax2.transAxes,
         fontsize=14, fontweight='bold', va='bottom', ha='left')

# ---------- 图例位置调整（避免遮挡 (B) 和条形） ----------
# 原图例在 (1.02, 1)，会挡住 (B) 标签，现在放到右图下方或右侧靠外
legend_elements = [Patch(facecolor=palette[d], edgecolor='black', label=d.capitalize()) for d in diets]
ax2.legend(handles=legend_elements,
           loc='lower left',          # 改为左下角，避免与 (B) 重叠
           bbox_to_anchor=(1.02, 0),  # 在图外右侧偏下
           frameon=False, fontsize=8)

# ================== 5. 整体调整与保存 ==================
plt.tight_layout()
plt.subplots_adjust(wspace=0.4)  # 增加左右子图间距

for ext in ['png', 'pdf', 'svg']:
    plt.savefig(os.path.join(OUTPUT_DIR, f'Combined_Violin_Ranking_AB.{ext}'),
                dpi=600 if ext=='png' else None, bbox_inches='tight')
plt.close()

print(f"✅ 大图（含 (A)(B) 标注）已保存至 {OUTPUT_DIR}")
print(f"   Kruskal‑Wallis H={h_stat:.3f}, p={p_kw:.6f}")
print(f"   显著对比数: {len(pair_list)}")
