"""Figure 3: independent MCP-counter scores and conditional associations."""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

ROOT=Path(__file__).resolve().parents[1]
src=ROOT/'results/virtual_validation'
out=ROOT/'manuscript/figures/Figure_3_virtual_immune_validation.png'
d=pd.read_csv(src/'MCPcounter_correlations.csv')
p=pd.read_csv(src/'MCPcounter_partial_rank.csv')
cohorts=['TCGA','CGGA693','CGGA325','GSE16011','GSE149009']
proxies=['Monocytic lineage','T cells','Neutrophils','Endothelial cells','Fibroblasts']
w=d[d.proxy.isin(proxies)].pivot(index='proxy',columns='cohort',values='rho').loc[proxies,cohorts]
fig,ax=plt.subplots(1,2,figsize=(10.6,4.2),gridspec_kw={'width_ratios':[1.5,1]})
sns.heatmap(w,ax=ax[0],cmap='RdBu_r',center=0,vmin=-.8,vmax=.8,annot=True,fmt='.2f',cbar_kws={'label':'Spearman rho'})
ax[0].set_title('A  Published cell-population scores');ax[0].set_xlabel('');ax[0].set_ylabel('')
p=p.set_index('cohort').loc[cohorts]
bars=ax[1].barh(range(len(cohorts)),p.partial_rank_rho,color=['#7b8fa7']*5)
ax[1].set_yticks(range(len(cohorts)),cohorts);ax[1].invert_yaxis();ax[1].axvline(0,color='black',lw=.8)
ax[1].set_xlim(0,.6);ax[1].set_xlabel('Partial rank correlation')
ax[1].set_title('B  Monocytic score after adjustment')
for b,r in zip(bars,p.itertuples()):
    ax[1].text(b.get_width()+.012,b.get_y()+b.get_height()/2,f'P={r.permutation_p:.3g}',va='center',fontsize=9)
fig.tight_layout();fig.savefig(out,dpi=300,bbox_inches='tight');plt.close(fig)
print(out)
