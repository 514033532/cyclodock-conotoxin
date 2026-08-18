#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
独立桑基图脚本（仅绘制桑基图）
输入：merged_with_target_info2Organism_diet.csv
输出：sankey_diet_species.html（交互式）保存到 fenlei/fig1/
"""

import os
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import plot
from collections import Counter

# ================== 配置 ==================
INPUT_FILE = '/home/zhh2/CONO-PEPTIDE2026/2026result/jianzhan/merged_with_target_info2Organism_diet.csv'
OUTPUT_DIR = '/home/zhh2/CONO-PEPTIDE2026/2026result/jianzhan/fenlei/fig1/'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 列名（请确认与实际文件一致）
COL_DIET = 'Organism_diet'
COL_SPECIES = 'Species'
COL_SUPERFAMILY = 'Superfamily'
COL_PHARMA = 'Pharmacological_Classification'
COL_TARGET = 'Target_Protein'

# 每层最多保留的类别数（取出现频率最高的前 N 个）
TOP_N = 30

# ================== 读取数据 ==================
df = pd.read_csv(INPUT_FILE, sep='\t', encoding='utf-8-sig')
df.columns = df.columns.str.strip()

# 只保留需要的列，并删除含有缺失值的行
required_cols = [COL_DIET, COL_SPECIES, COL_SUPERFAMILY, COL_PHARMA, COL_TARGET]
df = df[required_cols].dropna()

# 去除食性列的前后空格，过滤掉 'unknown'
df[COL_DIET] = df[COL_DIET].str.strip()
df = df[df[COL_DIET] != 'unknown']
print(f"总有效记录数：{len(df)}")

# ================== 定义层级（顺序固定） ==================
layers = [
    ('diet', COL_DIET),
    ('species', COL_SPECIES),
    ('superfamily', COL_SUPERFAMILY),
    ('pharma', COL_PHARMA),
    ('target', COL_TARGET)
]

# ================== 逐层过滤（保留每层 TOP_N） ==================
filtered_df = df.copy()
for layer_key, col in layers:
    top_vals = filtered_df[col].value_counts().head(TOP_N).index
    filtered_df = filtered_df[filtered_df[col].isin(top_vals)]

print(f"过滤后记录数：{len(filtered_df)}")

# ================== 构建节点索引 ==================
all_nodes = []
node_index = {}
for layer_key, col in layers:
    for val in filtered_df[col].unique():
        if val not in node_index:
            node_index[val] = len(all_nodes)
            all_nodes.append(val)

# ================== 构建边（相邻层之间） ==================
sources, targets = [], []
for _, row in filtered_df.iterrows():
    for i in range(len(layers)-1):
        src = row[layers[i][1]]
        tgt = row[layers[i+1][1]]
        sources.append(node_index[src])
        targets.append(node_index[tgt])

# 聚合相同边的计数
link_counter = Counter(zip(sources, targets))
sources_agg = [s for (s, t), cnt in link_counter.items()]
targets_agg = [t for (s, t), cnt in link_counter.items()]
values_agg = [cnt for (s, t), cnt in link_counter.items()]

print(f"总边数（聚合后）：{len(sources_agg)}")

# ================== 绘制交互式桑基图 ==================
fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=all_nodes,
        color="lightblue"
    ),
    link=dict(
        source=sources_agg,
        target=targets_agg,
        value=values_agg
    )
)])

fig.update_layout(
    title=f"Diet → Species → Superfamily → Pharmacology → Target (Top {TOP_N} per level)",
    font=dict(size=10),
    width=1000,
    height=700
)

# 保存为 HTML 文件
output_path = os.path.join(OUTPUT_DIR, 'sankey_diet_species.html')
plot(fig, filename=output_path, auto_open=False)
print(f"✅ 桑基图已保存至：{output_path}")
