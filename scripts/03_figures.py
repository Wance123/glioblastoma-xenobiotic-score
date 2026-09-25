from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parents[1]
out=root/'results'
m=pd.read_csv(out/'survival_models.csv')
r=pd.read_csv(out/'recurrence_unpaired_gene_models.csv')
fig,ax=plt.subplots(figsize=(7,3.5))
models=['unadjusted','age_adjusted','clinical_complete_case']
for j,(cohort,color) in enumerate([('CGGA693','#246a73'),('CGGA325','#a33a35')]):
    d=m[m.cohort.eq(cohort)].set_index('model').loc[models]
    y=[i+(j-.5)*.18 for i in range(3)]
    ax.errorbar(d.HR_per_score_unit,y,xerr=[d.HR_per_score_unit-d.CI_low,d.CI_high-d.HR_per_score_unit],fmt='o',color=color,capsize=3,label=cohort)
ax.axvline(1,color='gray',lw=1)
ax.set_yticks(range(3),['Unadjusted','Age adjusted','Clinical complete case'])
ax.set_xlabel('Hazard ratio per unit of fixed four-gene score (95% CI)')
ax.legend(frameon=False)
fig.tight_layout();fig.savefig(out/'figure1_survival_forest.png',dpi=200);plt.close(fig)
fig,ax=plt.subplots(figsize=(7,3.5))
for j,(cohort,color) in enumerate([('CGGA693','#246a73'),('CGGA325','#a33a35')]):
    d=r[r.cohort.eq(cohort)].set_index('gene').loc[['MPI','PMM2','GMPPA','GMPPB']]
    y=[i+(j-.5)*.18 for i in range(4)]
    ax.errorbar(d.recurrent_minus_primary_log2,y,xerr=[d.recurrent_minus_primary_log2-d.CI_low,d.CI_high-d.recurrent_minus_primary_log2],fmt='o',color=color,capsize=3,label=cohort)
ax.axvline(0,color='gray',lw=1)
ax.set_yticks(range(4),['MPI','PMM2','GMPPA','GMPPB'])
ax.set_xlabel('Recurrent minus primary log2 expression, adjusted for age (95% CI)')
ax.legend(frameon=False)
fig.tight_layout();fig.savefig(out/'figure2_recurrence_forest.png',dpi=200);plt.close(fig)
