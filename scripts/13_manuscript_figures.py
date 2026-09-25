"""Figures for observational manuscript; all displayed estimates are saved CSVs."""
from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results';F=ROOT/'manuscript/figures';F.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
c=pd.read_csv(R/'xenobiotic_marker_correlations_all_cohorts.csv')
m=pd.read_csv(R/'pathway_cox_all.csv')
fig,ax=plt.subplots(1,2,figsize=(10,3.5),gridspec_kw={'width_ratios':[1.2,1]})
coh=['TCGA','CGGA693','CGGA325']; colors=['#265D73','#BF634A','#6D7C47']
for i,co in enumerate(coh):
    d=c[(c.cohort==co)&c.proxy.isin(['myeloid_identity','macrophage_state','tumor_glial_proxy'])].set_index('proxy')
    x=np.arange(3)+(i-1)*.22
    ax[0].bar(x,d.loc[['myeloid_identity','macrophage_state','tumor_glial_proxy'],'rho'],width=.20,color=colors[i],label=co)
ax[0].axhline(0,color='gray',lw=.8);ax[0].set_xticks(range(3),['Myeloid identity','Macrophage state','Glial proxy'],rotation=15,ha='right')
ax[0].set_ylim(-.8,1);ax[0].set_ylabel('Spearman correlation with xenobiotic score');ax[0].legend(frameon=False,fontsize=8)
xeno=m[(m.scope=='GBM')&(m.model=='age_MGMT')&(m.pathway=='XENOBIOTIC_METABOLISM')].set_index('cohort')
for i,co in enumerate(coh):
    q=xeno.loc[co];ax[1].scatter(q.HR_per_SD,i,color=colors[i],s=50)
    ax[1].plot([q.CI_low,q.CI_high],[i,i],color=colors[i],lw=2)
ax[1].axvline(1,color='gray',lw=.8);ax[1].set_yticks(range(3),['TCGA','CGGA693','CGGA325']);ax[1].invert_yaxis()
ax[1].set_xlabel('Hazard ratio per SD of pathway score');ax[1].set_xlim(.5,2.0)
fig.tight_layout();fig.savefig(F/'Figure_1_bulk_composition_and_survival.png',dpi=220,bbox_inches='tight');plt.close(fig)
s=pd.read_csv(R/'spatial_xenobiotic/cohort_summary.csv')
fig,ax=plt.subplots(figsize=(7,3.6))
for i,(co,color) in enumerate([('GSE237183','#265D73'),('GSE318562','#BF634A')]):
    d=s[(s.dataset==co)&(s.outcome=='myeloid_identity')].set_index('model')
    ax.plot([0,1],[d.loc['rho_unadjusted','median_rho'],d.loc['rho_adjusted','median_rho']],marker='o',lw=2,color=color,label=f'{co} ({int(d.iloc[0].n_patients)} patients)')
ax.set_xticks([0,1],['Unadjusted','Adjusted for hypoxia,\nglial proxy and depth']);ax.set_ylabel('Median patient spatial correlation')
ax.set_ylim(0,.32);ax.legend(frameon=False);fig.tight_layout();fig.savefig(F/'Figure_2_spatial_coenrichment.png',dpi=220,bbox_inches='tight');plt.close(fig)
