#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
绘制三张独立图并合并为一张大图（ABC标注，竖排 3×1）
图1：靶蛋白平均 HADDOCK 得分排名 + 数据点
图2：按食性分组的靶蛋白平均得分（前20结合最强）+ 数据点，优化 X 轴排列
图3：按药理分类的平均得分 + 数据点
输出：独立图片 + 大图（PNG, PDF, SVG）
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

# ================== 配置 ==================
INPUT_FILE = '/home/zhh2/CONO-PEPTIDE2026/2026result/jianzhan/merged_with_target_info2Organism_diet.csv'
OUTPUT_DIR = '/home/zhh2/CONO-PEPTIDE2026/2026result/jianzhan/fenlei/fig3'
os.makedirs(OUTPUT_DIR, exist_ok=True)

COL_DIET = 'Organism_diet'
COL_TARGET = 'Target_Protein'
COL_SCORE = 'HADDOCK_score_[a.u.]'
COL_PHARMA = 'Pharmacological_Classification'

# ================== 读取数据（过滤 unknown） ==================
df = pd.read_csv(INPUT_FILE, sep='\t', encoding='utf-8-sig')
df.columns = df.columns.str.strip()
df = df[[COL_DIET, COL_TARGET, COL_SCORE, COL_PHARMA]].dropna()
df[COL_DIET] = df[COL_DIET].str.strip()
df = df[df[COL_DIET] != 'unknown']
print(f"有效数据行数：{len(df)}")
print("食性分布：")
print(df[COL_DIET].value_counts())

# 去除靶蛋白名称中可能的尾随空格，并统一大小写（可选）
df[COL_TARGET] = df[COL_TARGET].str.strip()
# 如有必要，可合并相似名称（这里不操作，仅做清洗）

DIET_ORDER = ['molluscivorous', 'piscivorous', 'vermivorous', 'generalist']
existing_diets = sorted(df[COL_DIET].unique(), key=lambda x: DIET_ORDER.index(x) if x in DIET_ORDER else len(DIET_ORDER))

# ================== 计算各类平均值 ==================
target_mean = df.groupby(COL_TARGET)[COL_SCORE].mean().sort_values()
pharma_mean = df.groupby(COL_PHARMA)[COL_SCORE].mean().sort_values()

# ================== 设置全局绘图风格 ==================
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.major.width'] = 0.8
plt.rcParams['ytick.major.width'] = 0.8

# ================== 图1：靶蛋白排名（水平条形图 + 散点） ==================
fig1, ax1 = plt.subplots(figsize=(6, max(4, len(target_mean)*0.2)))
if len(target_mean) > 30:
    target_mean_plot = target_mean.head(30)
else:
    target_mean_plot = target_mean
y_pos = np.arange(len(target_mean_plot))
ax1.barh(y_pos, target_mean_plot.values, color='#1f4e79', edgecolor='black', linewidth=0.5, alpha=0.6)
for i, (target, mean_val) in enumerate(target_mean_plot.items()):
    scores = df[df[COL_TARGET]==target][COL_SCORE].values
    y_jitter = np.random.normal(i, 0.04, size=len(scores))
    ax1.scatter(scores, y_jitter, s=5, color='black', alpha=0.3, linewidth=0)
ax1.set_yticks(y_pos)
ax1.set_yticklabels(target_mean_plot.index, fontsize=7)
ax1.invert_yaxis()
ax1.set_xlabel('Average HADDOCK Score (a.u.)', fontsize=9, fontweight='bold')
ax1.set_title('Target Proteins Ranked by Average HADDOCK Score', fontsize=10, fontweight='bold')
ax1.grid(axis='x', linestyle='--', alpha=0.3)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig1_Target_Ranking.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig1_Target_Ranking.pdf'), bbox_inches='tight')
plt.close()

# ================== 图2：按食性分组，前20个结合最强的靶蛋白（优化X轴） ==================
top_n = 20
target_overall_mean = df.groupby(COL_TARGET)[COL_SCORE].mean()
target_order = target_overall_mean.sort_values().head(top_n).index

plot_df = df[df[COL_TARGET].isin(target_order)].groupby([COL_DIET, COL_TARGET])[COL_SCORE].mean().reset_index()
plot_df[COL_TARGET] = pd.Categorical(plot_df[COL_TARGET], categories=target_order, ordered=True)

# 动态宽度：每个靶蛋白占用 0.5 英寸宽度，至少 8 英寸
fig_width = max(8, top_n * 0.5)
fig2, ax2 = plt.subplots(figsize=(fig_width, 5))

sns.barplot(data=plot_df, x=COL_TARGET, y=COL_SCORE, hue=COL_DIET,
            hue_order=DIET_ORDER, palette='Set2', edgecolor='black', linewidth=0.5, ax=ax2, alpha=0.7)

sns.stripplot(data=df[df[COL_TARGET].isin(target_order)], x=COL_TARGET, y=COL_SCORE, hue=COL_DIET,
              hue_order=DIET_ORDER, dodge=True, size=3, color='black', alpha=0.3, jitter=0.2, ax=ax2)

ax2.legend_.remove()
palette = sns.color_palette("Set2", n_colors=len(DIET_ORDER))
existing_palette = [palette[DIET_ORDER.index(d)] for d in existing_diets]
legend_elements = [Patch(facecolor=existing_palette[i], edgecolor='black', label=d) for i, d in enumerate(existing_diets)]
ax2.legend(handles=legend_elements, title='Diet', bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)

ax2.set_xlabel('Target Protein', fontsize=9, fontweight='bold')
ax2.set_ylabel('Average HADDOCK Score (a.u.)', fontsize=9, fontweight='bold')
ax2.set_title(f'Top {top_n} Target Proteins with Strongest Binding (by overall mean)', fontsize=10, fontweight='bold')
# 关键优化：x轴标签旋转90°，紧凑排列，字体调小
ax2.tick_params(axis='x', labelrotation=90, labelsize=6.5)
ax2.tick_params(axis='y', labelsize=8)
ax2.grid(axis='y', linestyle='--', alpha=0.3)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# 自动调整y轴边距
y_min, y_max = ax2.get_ylim()
ax2.set_ylim(y_min - 0.5, y_max + 0.5)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig2_Target_by_Diet.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig2_Target_by_Diet.pdf'), bbox_inches='tight')
plt.close()

# ================== 图3：药理分类条形图 + 散点 ==================
fig3, ax3 = plt.subplots(figsize=(6, max(4, len(pharma_mean)*0.2)))
if len(pharma_mean) > 20:
    pharma_mean_plot = pharma_mean.head(20)
else:
    pharma_mean_plot = pharma_mean
y_pos = np.arange(len(pharma_mean_plot))
ax3.barh(y_pos, pharma_mean_plot.values, color='#377eb8', edgecolor='black', linewidth=0.5, alpha=0.6)
for i, (pharma, mean_val) in enumerate(pharma_mean_plot.items()):
    scores = df[df[COL_PHARMA]==pharma][COL_SCORE].values
    y_jitter = np.random.normal(i, 0.04, size=len(scores))
    ax3.scatter(scores, y_jitter, s=5, color='black', alpha=0.3, linewidth=0)
ax3.set_yticks(y_pos)
ax3.set_yticklabels(pharma_mean_plot.index, fontsize=7)
ax3.invert_yaxis()
ax3.set_xlabel('Average HADDOCK Score (a.u.)', fontsize=9, fontweight='bold')
ax3.set_title('Average HADDOCK Score by Pharmacological Family', fontsize=10, fontweight='bold')
ax3.grid(axis='x', linestyle='--', alpha=0.3)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig3_Pharma_Ranking.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig3_Pharma_Ranking.pdf'), bbox_inches='tight')
plt.close()

# ================== 组合大图（3行1列，竖排） ==================
fig, axes = plt.subplots(3, 1, figsize=(10, 18))
ax1, ax2, ax3 = axes

# ----- (A) 图1 -----
target_mean_plot = target_mean.head(30) if len(target_mean)>30 else target_mean
y_pos = np.arange(len(target_mean_plot))
ax1.barh(y_pos, target_mean_plot.values, color='#1f4e79', edgecolor='black', linewidth=0.5, alpha=0.6)
for i, (target, mean_val) in enumerate(target_mean_plot.items()):
    scores = df[df[COL_TARGET]==target][COL_SCORE].values
    y_jitter = np.random.normal(i, 0.04, size=len(scores))
    ax1.scatter(scores, y_jitter, s=5, color='black', alpha=0.3, linewidth=0)
ax1.set_yticks(y_pos)
ax1.set_yticklabels(target_mean_plot.index, fontsize=7)
ax1.invert_yaxis()
ax1.set_xlabel('Average HADDOCK Score (a.u.)', fontsize=9, fontweight='bold')
ax1.set_title('Target Proteins Ranked by Average HADDOCK Score', fontsize=10, fontweight='bold')
ax1.grid(axis='x', linestyle='--', alpha=0.3)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.text(-0.1, 1.02, '(A)', transform=ax1.transAxes, fontsize=14, fontweight='bold', va='bottom', ha='right')

# ----- (B) 图2（优化X轴）-----
sns.barplot(data=plot_df, x=COL_TARGET, y=COL_SCORE, hue=COL_DIET,
            hue_order=DIET_ORDER, palette='Set2', edgecolor='black', linewidth=0.5, ax=ax2, alpha=0.7)
sns.stripplot(data=df[df[COL_TARGET].isin(target_order)], x=COL_TARGET, y=COL_SCORE, hue=COL_DIET,
              hue_order=DIET_ORDER, dodge=True, size=3, color='black', alpha=0.3, jitter=0.2, ax=ax2)
ax2.legend_.remove()
existing_palette = [palette[DIET_ORDER.index(d)] for d in existing_diets]
legend_elements = [Patch(facecolor=existing_palette[i], edgecolor='black', label=d) for i, d in enumerate(existing_diets)]
ax2.legend(handles=legend_elements, title='Diet', bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
ax2.set_xlabel('Target Protein', fontsize=9, fontweight='bold')
ax2.set_ylabel('Average HADDOCK Score (a.u.)', fontsize=9, fontweight='bold')
ax2.set_title(f'Top {top_n} Target Proteins with Strongest Binding (by overall mean)', fontsize=10, fontweight='bold')
ax2.tick_params(axis='x', labelrotation=90, labelsize=6.5)
ax2.tick_params(axis='y', labelsize=8)
ax2.grid(axis='y', linestyle='--', alpha=0.3)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
y_min, y_max = ax2.get_ylim()
ax2.set_ylim(y_min - 0.5, y_max + 0.5)
ax2.text(-0.1, 1.02, '(B)', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='bottom', ha='right')

# ----- (C) 图3 -----
pharma_mean_plot = pharma_mean.head(20) if len(pharma_mean)>20 else pharma_mean
y_pos = np.arange(len(pharma_mean_plot))
ax3.barh(y_pos, pharma_mean_plot.values, color='#377eb8', edgecolor='black', linewidth=0.5, alpha=0.6)
for i, (pharma, mean_val) in enumerate(pharma_mean_plot.items()):
    scores = df[df[COL_PHARMA]==pharma][COL_SCORE].values
    y_jitter = np.random.normal(i, 0.04, size=len(scores))
    ax3.scatter(scores, y_jitter, s=5, color='black', alpha=0.3, linewidth=0)
ax3.set_yticks(y_pos)
ax3.set_yticklabels(pharma_mean_plot.index, fontsize=7)
ax3.invert_yaxis()
ax3.set_xlabel('Average HADDOCK Score (a.u.)', fontsize=9, fontweight='bold')
ax3.set_title('Average HADDOCK Score by Pharmacological Family', fontsize=10, fontweight='bold')
ax3.grid(axis='x', linestyle='--', alpha=0.3)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.text(-0.1, 1.02, '(C)', transform=ax3.transAxes, fontsize=14, fontweight='bold', va='bottom', ha='right')

plt.tight_layout()
plt.subplots_adjust(hspace=0.4)

for ext in ['png', 'pdf', 'svg']:
    plt.savefig(os.path.join(OUTPUT_DIR, f'Combined_ABC_Vertical_withPoints.{ext}'),
                dpi=600 if ext=='png' else None, bbox_inches='tight')
plt.close()

print(f"✅ 所有图表（图2 X轴已优化为垂直标签）已保存至：{OUTPUT_DIR}")
print("实际存在的食性：", existing_diets)
