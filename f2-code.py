#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
整合脚本：绘制唯一值统计图（两个多子图 + 一个旋转90°的水平条形图）
输出目录：/home/zhh2/CONO-PEPTIDE2026/2026result/jianzhan/Stat_Plots2/Combined_Figures
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# ==================== 配置区 ====================
input_csv = "/home/zhh2/CONO-PEPTIDE2026/2026result/jianzhan/merged_with_target_info2.csv"
output_dir = Path("/home/zhh2/CONO-PEPTIDE2026/2026result/jianzhan/fenlei/fig2")
output_dir.mkdir(parents=True, exist_ok=True)

# 全部10列（顺序可自定义）
all_cols = [
    "Pharmacological_Classification",
    "Disulfide_Connectivity",
    "Source_Type",
    "Target_Protein",
    "Protein_Type",
    "Species",
    "Sequence_Type",
    "Molecular_Class",
    "Superfamily",
    "Protein_Name"
]

# 大气风格配置参数
DPI = 300
FONT_FAMILY = 'Arial'
plt.rcParams['font.family'] = FONT_FAMILY
plt.rcParams['font.size'] = 13

# ==================== 读取并计算数据 ====================
def read_csv_robust(filepath):
    strategies = [
        {'sep': '\t', 'engine': None, 'on_bad_lines': 'skip', 'low_memory': False},
        {'sep': ',', 'engine': None, 'on_bad_lines': 'skip', 'low_memory': False},
        {'sep': None, 'engine': 'python', 'on_bad_lines': 'skip'},
    ]
    for kwargs in strategies:
        try:
            df = pd.read_csv(filepath, encoding='utf-8', **kwargs)
            if not df.empty:
                return df
        except:
            continue
    raise ValueError("无法读取数据，请检查文件路径。")

df = read_csv_robust(input_csv)
df = df.replace(r'^\s*$', np.nan, regex=True)

existing_cols = [c for c in all_cols if c in df.columns]

# ==================== 生成唯一值统计并排序 ====================
def generate_unique_counts(df, cols):
    result = []
    for col in cols:
        cnt = df[col].dropna().nunique()
        result.append((col, cnt))
    result.sort(key=lambda x: x[1], reverse=True)
    return result

unique_data = generate_unique_counts(df, existing_cols)

# 固定拆分为前6个 与 后4个
unique_6 = unique_data[:6]
unique_4 = unique_data[6:]

# ==================== 统一风格绘图函数（用于多子图） ====================
def plot_consistent_unique_chart(items, start_letter_idx, layout_shape, figure_number):
    n_rows, n_cols = layout_shape
    fig_width = n_cols * 8.5
    fig_height = n_rows * 7.5
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))
    axes = axes.flatten()
    
    plt.rcParams['font.family'] = FONT_FAMILY
    plt.rcParams['font.size'] = 13

    max_val = max([val for _, val in items]) if items else 0

    for i, (col, cnt) in enumerate(items):
        ax = axes[i]
        cat_names = [col]
        cat_vals = [cnt]
        bars = ax.barh(cat_names, cat_vals, color='#1f4e79', edgecolor='black', linewidth=0.8, height=0.5)
        ax.text(cnt + (max_val * 0.02), 0, f'{int(cnt)}', 
                ha='left', va='center', fontsize=15, color='black', fontweight='bold')
        alpha_idx = start_letter_idx + i
        ax.set_title(f"({chr(65 + alpha_idx)}) {col}", fontsize=16, fontweight='bold', loc='left', pad=15)
        ax.set_xlabel("Unique Count", fontsize=14, fontweight='bold')
        ax.tick_params(axis='y', labelsize=13)
        ax.tick_params(axis='x', labelsize=12)
        ax.set_xlim(0, max_val * 1.15)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.2)
        ax.spines['bottom'].set_linewidth(1.2)
        ax.grid(axis='x', linestyle='--', alpha=0.4, linewidth=0.7)

    for j in range(len(items), len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout(pad=3.0, w_pad=2.5, h_pad=3.5)
    pdf_path = output_dir / f"Combined_Unique_{figure_number}.pdf"
    png_path = output_dir / f"Combined_Unique_{figure_number}.png"
    plt.savefig(pdf_path, dpi=DPI, bbox_inches='tight', format='pdf')
    plt.savefig(png_path, dpi=DPI, bbox_inches='tight', format='png')
    plt.close()
    print(f"✅ 已保存: Combined_Unique_{figure_number} (PDF & PNG)")

# ==================== 生成第一张多子图（前6列） ====================
if len(unique_6) > 0:
    print("正在绘制第一大张图 (前6列的唯一点计)...")
    plot_consistent_unique_chart(unique_6, start_letter_idx=0, layout_shape=(2, 3), figure_number="1")

# ==================== 生成第二张多子图（后4列） ====================
if len(unique_4) > 0:
    print("正在绘制第二大张图 (后4列的唯一点计)...")
    plot_consistent_unique_chart(unique_4, start_letter_idx=6, layout_shape=(2, 2), figure_number="2")

# ==================== 生成第三张图（顺时针旋转90°的水平条形图） ====================
print("正在绘制第三张图 (旋转90°的水平条形图)...")
# 使用与上面相同的唯一值数据，但以单图形式展示全部类别
# 这里我们直接使用之前计算好的 unique_data（全部9个类别）
# 但确保数据与之前一致，可以重新从 df 生成，但 unique_data 已经包含了所有列（可能不足9个）
# 如果 existing_cols 少于9个，则 unique_data 长度会少，我们用全部数据
# 为了与第二段代码数据一致，我们使用硬编码的数据，也可从 unique_data 取全部
# 但第二段代码是固定的9个类别，我们优先使用 unique_data 的全部项（保证一致）
all_items = unique_data  # 已经排序，包含所有列
# 但第二段代码中 Variable 是硬编码的，为了统一，我们采用 unique_data 的列名和计数
# 如果 unique_data 顺序可能与硬编码不同，但没关系，我们直接使用
# 但为了标题 (K)，我们按原有格式，用全部列

# 转换为 DataFrame 用于绘图
plot_df = pd.DataFrame(all_items, columns=["Variable", "Unique_Count"])
# 按计数降序（已经排过）
# 如果希望与硬编码顺序一致，可以手动指定顺序，但这里保持排序

fig, ax = plt.subplots(figsize=(8, 12))
bars = ax.barh(plot_df["Variable"], plot_df["Unique_Count"],
               color='#1f4e79', edgecolor='black', linewidth=0.8, height=0.6)
for bar in bars:
    width = bar.get_width()
    ax.text(width + 0.5, bar.get_y() + bar.get_height()/2,
            f'{int(width)}', va='center', ha='left',
            fontsize=14, fontweight='bold', color='black')
ax.set_xlabel("Unique Count", fontsize=14, fontweight='bold')
ax.set_ylabel("Category", fontsize=14, fontweight='bold')
ax.set_title("(K) Unique Counts per Category", fontsize=16, fontweight='bold', pad=20)
ax.tick_params(axis='y', labelsize=12)
ax.tick_params(axis='x', labelsize=12)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.2)
ax.spines['bottom'].set_linewidth(1.2)
ax.grid(axis='x', linestyle='--', alpha=0.4, linewidth=0.7)

plt.tight_layout()
pdf_path = output_dir / "Combined_Unique_Rotated90.pdf"
png_path = output_dir / "Combined_Unique_Rotated90.png"
plt.savefig(pdf_path, dpi=DPI, bbox_inches='tight')
plt.savefig(png_path, dpi=DPI, bbox_inches='tight')
plt.close()
print(f"✅ 已保存: Combined_Unique_Rotated90 (PDF & PNG)")

print(f"\n🎉 全部唯一点计图生成完毕！文件保存在: {output_dir}")
