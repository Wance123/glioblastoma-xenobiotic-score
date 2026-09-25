import os
"""TCGA composition audit and purity-adjusted sensitivity for six pathways."""
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[1]
SRC=Path(os.environ['JOB46_ROOT'])/'reanalysis_v2'
z=pd.read_csv(SRC/'results/all_cohort_scores.csv',low_memory=False)
z=z[z.cohort.eq('TCGA')&z.eligible_GBM.eq(True)].copy()
expr=pd.read_csv(SRC/'data/TCGA_GBM_log_expression.tsv.gz',sep='\t').set_index('sample')
z=z.set_index('sample')
for gene in ['NDST3','ABCA5']:
    z[gene]=expr.loc[gene].reindex(z.index)
variables=['NDST3','ABCA5']+['HALLMARK_'+x for x in ['GLYCOLYSIS','FATTY_ACID_METABOLISM','OXIDATIVE_PHOSPHORYLATION','XENOBIOTIC_METABOLISM','CHOLESTEROL_HOMEOSTASIS','BILE_ACID_METABOLISM']]
rows=[];models=[]
for v in variables:
    for y in ['purity','immune_estimate','myeloid_identity','tumor_glial_proxy','hypoxia']:
        a=z[[v,y]].apply(pd.to_numeric,errors='coerce').dropna()
        if len(a)>=30:
            rho,p=spearmanr(a[v],a[y]);rows.append({'variable':v,'correlate':y,'n':len(a),'spearman_rho':rho,'p':p})
    if v.startswith('HALLMARK_'):
        d=z[['time','event','age','MGMT','purity',v]].apply(pd.to_numeric,errors='coerce').dropna()
        d['score_z']=(d[v]-d[v].mean())/d[v].std(ddof=0)
        fit=CoxPHFitter(penalizer=.01).fit(d.drop(columns=v),'time','event')
        q=fit.summary.loc['score_z']
        models.append({'pathway':v,'n':len(d),'events':int(d.event.sum()),'HR_per_SD':float(np.exp(q['coef'])),'p':float(q['p'])})
r=pd.DataFrame(rows);r['BH_FDR_per_correlate']=r.groupby('correlate').p.transform(lambda p:multipletests(p,method='fdr_bh')[1])
r.to_csv(ROOT/'results/TCGA_composition_correlations.csv',index=False)
m=pd.DataFrame(models);m['BH_FDR_six']=multipletests(m.p,method='fdr_bh')[1]
m.to_csv(ROOT/'results/TCGA_purity_adjusted_pathway_cox.csv',index=False)
print(r[r.correlate.isin(['purity','immune_estimate'])].to_string(index=False))
print(m.to_string(index=False))
