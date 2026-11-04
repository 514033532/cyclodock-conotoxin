import json, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,
                     'axes.spines.top':False,'axes.spines.right':False,'font.family':'sans-serif'})
D=json.load(open('fig_data_extra.json'))
okabe=['#0072B2','#D55E00','#009E73','#CC79A7']
ONE=3.46; TWO=7.09

# Fig 4: 盲接靶内对比 (4 面板, 2x2)
fig,axes=plt.subplots(2,2,figsize=(TWO,5.2))
for ax,(name,groups) in zip(axes.ravel(), D['fig4_blind'].items()):
    for gi,(gn,vals) in enumerate(groups.items()):
        for i,(t,s) in enumerate(vals.items()):
            ax.scatter(gi+np.random.default_rng(7+gi).uniform(-0.1,0.1), s, s=22, color=okabe[gi], edgecolor='k', linewidth=0.4, zorder=3)
        med=np.median(list(vals.values()))
        ax.plot([gi-0.3,gi+0.3],[med,med],color=okabe[gi],lw=1.2)
    ax.set_xticks([0,1]); ax.set_xticklabels(['Positive','Negative'],fontsize=6.5)
    ax.set_title(name,fontsize=8.5)
    ax.set_ylabel('HADDOCK score')
fig.legend(['Positive','Negative'],loc='lower center',ncol=2,fontsize=7,frameon=False)
plt.tight_layout(rect=[0,0.03,1,1]); plt.savefig('Fig4_blind.png',dpi=300); plt.savefig('Fig4_blind.pdf'); plt.close()

# Fig 5: 分数尺度 (分组柱: 毒素 x 表面)
fig,ax=plt.subplots(figsize=(ONE,2.9))
toxins=list(D['fig5_scale'].keys()); surfaces=['CaV2.2','NaV1.7','NaV1.4']
w=0.25; x=np.arange(len(toxins))
for si,surf in enumerate(surfaces):
    vals=[D['fig5_scale'][t][surf] for t in toxins]
    ax.bar(x+(si-1)*w, vals, w, color=okabe[si], label=surf)
ax.set_xticks(x); ax.set_xticklabels(toxins)
ax.set_ylabel('HADDOCK score (lower=better)')
ax.set_title('Same toxin, different surfaces: score scales',fontsize=8.5)
ax.legend(fontsize=6.5, frameon=False, title='Surface', title_fontsize=7,
           loc='center left', bbox_to_anchor=(1.01, 0.5))
plt.tight_layout(); plt.savefig('Fig5_scale.png',dpi=300,bbox_inches='tight'); plt.savefig('Fig5_scale.pdf',bbox_inches='tight'); plt.close()

# Fig 6: 循环 vs 非循环 (2 面板)
fig,axes=plt.subplots(1,2,figsize=(TWO,3.0),sharey=False)
for ax,(name,vals) in zip(axes, D['fig6_circular'].items()):
    for i,(t,s) in enumerate(vals.items()):
        ax.scatter(i, s, s=26, color=okabe[0] if t in ('KIIIA','GIIIA','CnIIIC') else okabe[1],
                   edgecolor='k', linewidth=0.4, zorder=3)
        ax.annotate(t, (i,s), textcoords='offset points', xytext=(0,5), fontsize=6, ha='center')
    ax.set_title(name,fontsize=8.5); ax.set_ylabel('HADDOCK score')
    ax.set_xticks([])
axes[0].annotate('apparent separation',xy=(0.5,0.9),xycoords='axes fraction',ha='center',fontsize=7,color='#D55E00')
axes[1].annotate('no separation',xy=(0.5,0.9),xycoords='axes fraction',ha='center',fontsize=7,color='#D55E00')
plt.tight_layout(); plt.savefig('Fig6_circular.png',dpi=300); plt.savefig('Fig6_circular.pdf'); plt.close()
print('Fig4-6 已生成')
