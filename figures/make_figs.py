import json, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, '/home/zhh2/.claude/skills/scipilot-figure-skill')
try:
    from setup_style import setup_style
    setup_style(journal='nature', lang='en')
except Exception as e:
    print('setup_style 不可用, 用内置预设:', e)
    plt.rcParams.update({'font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,
                         'axes.spines.top':False,'axes.spines.right':False,'font.family':'sans-serif'})
D = json.load(open('fig_data.json'))
okabe = ['#0072B2','#D55E00','#009E73','#CC79A7','#56B4E9','#E69F00']
ONE=3.46; TWO=7.09

# ---- Fig 1: 活性矩阵热图 (双栏) ----
fig, ax = plt.subplots(figsize=(TWO, 4.6))
M=D['matrix']; T=D['toxins']; Tg=D['targets']
code={'T':0,'W':1,'N':2,'-':3}
mat=np.full((len(T),len(Tg)),3)
for i,t in enumerate(T):
    for j,g in enumerate(Tg):
        if f'{t}|{g}' in M: mat[i,j]=code[M[f'{t}|{g}']]
cmap=matplotlib.colors.ListedColormap(['#1b7837','#a6dba0','#cccccc','#ffffff'])
im=ax.imshow(mat,cmap=cmap,vmin=-0.5,vmax=3.5,aspect='auto')
ax.set_xticks(range(len(Tg))); ax.set_xticklabels(Tg, rotation=45, ha='right')
ax.set_yticks(range(len(T))); ax.set_yticklabels(T)
for i in range(len(T)):
    for j in range(len(Tg)):
        if mat[i,j]<3: ax.text(j,i,['T','W','N'][mat[i,j]],ha='center',va='center',fontsize=7)
ax.grid(False)
cbar=fig.colorbar(im, ax=ax, ticks=[0,1,2,3], fraction=0.02)
cbar.ax.set_yticklabels(['True positive','Weak','Negative','Not tested'], fontsize=6.5)
ax.set_title('Benchmark set: literature activity matrix', fontsize=9)
ax.set_xlabel('Human target'); ax.set_ylabel('Toxin')
plt.tight_layout(); plt.savefig('Fig1_activity_matrix.png',dpi=300); plt.savefig('Fig1_activity_matrix.pdf'); plt.close()

# ---- Fig 2: 交叉设计分数 (双栏, 3 面板, stripplot 直列每个点) ----
fig, axes = plt.subplots(1,3, figsize=(TWO, 2.9), sharey=False)
for ax,(name,groups) in zip(axes, D['fig2'].items()):
    xs=[]; ys=[]; labels=[]
    for gi,(gn,vals) in enumerate(groups.items()):
        for i,(tox,sc) in enumerate(vals.items()):
            jitter=np.random.default_rng(42+gi).uniform(-0.12,0.12)
            ax.scatter(gi+jitter, sc, s=22, color=okabe[gi], edgecolor='k', linewidth=0.4, zorder=3)
    ax.set_xticks([0,1]); ax.set_xticklabels(list(groups.keys()), fontsize=6.5, rotation=12)
    ax.set_title(name, fontsize=8.5)
    ax.set_ylabel('HADDOCK score (lower = better)')
    # 组内中位线
    for gi,(gn,vals) in enumerate(groups.items()):
        med=np.median(list(vals.values()))
        ax.plot([gi-0.3,gi+0.3],[med,med], color=okabe[gi], lw=1.2, zorder=2)
fig.legend(['Positive','Negative'], loc='lower center', ncol=2, fontsize=7, frameon=False)
plt.tight_layout(rect=[0,0.04,1,1]); plt.savefig('Fig2_crossover.png',dpi=300); plt.savefig('Fig2_crossover.pdf'); plt.close()

# ---- Fig 3: 口袋组成堆叠柱 (单栏) ----
fig, ax = plt.subplots(figsize=(ONE, 2.9))
cats=list(D['fig3']['NaV1.7 pocket-1 (6N4I)'].keys())
x=np.arange(2); bottom=np.zeros(2)
for ci,c in enumerate(cats):
    vals=[D['fig3'][p][c] for p in D['fig3']]
    ax.bar(x, vals, bottom=bottom, color=okabe[ci], label=c, width=0.5)
    bottom+=np.array(vals)
ax.set_xticks(x); ax.set_xticklabels(['NaV1.7\npocket-1 (6N4I)','NaV1.4\npocket-1 (6AGF)'], fontsize=7)
ax.set_ylabel('Atom count in pocket')
ax.legend(fontsize=6, frameon=False, loc='upper left', bbox_to_anchor=(1.0,1.0))
ax.set_title('Pocket composition audit', fontsize=9)
plt.tight_layout(); plt.savefig('Fig3_pocket_audit.png',dpi=300); plt.savefig('Fig3_pocket_audit.pdf'); plt.close()
print('3 图已生成')
