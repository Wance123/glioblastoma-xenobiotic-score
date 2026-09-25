import os
"""Locked NDST3/ABCA5 test in independent TCGA IDH-WT primary GBM.

The two genes were selected using CGGA before reading TCGA outcomes. This
script uses the already audited TCGA expression/clinical joins from job46.
"""
from pathlib import Path
import numpy as np,pandas as pd
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[1]
SRC=Path(os.environ['JOB46_ROOT'])/'reanalysis_v2'
ledger=pd.read_csv(SRC/'results/TCGA_all_sample_ledger.csv').set_index('sample')
eligible=ledger[ledger.eligible_GBM.eq(True)].copy()
expr=pd.read_csv(SRC/'data/TCGA_GBM_log_expression.tsv.gz',sep='\t').set_index('sample')
out=[]
for gene in ['NDST3','ABCA5']:
    for label,covars in [('age',['age']),('age_MGMT',['age','MGMT']),('age_MGMT_KPS',['age','MGMT','KPS'])]:
        d=eligible[['time','event']+covars].join(expr.loc[gene].rename('expression'),how='inner').dropna()
        d=d[d.time>0].copy()
        d['expression_z']=(d.expression-d.expression.mean())/d.expression.std(ddof=0)
        fit=CoxPHFitter(penalizer=.01).fit(d.drop(columns='expression'),'time','event')
        z=fit.summary.loc['expression_z']
        out.append({'gene':gene,'model':label,'n':len(d),'events':int(d.event.sum()),
                    'HR_per_SD':float(np.exp(z['coef'])),'CI_low':float(np.exp(z['coef lower 95%'])),
                    'CI_high':float(np.exp(z['coef upper 95%'])),'p':float(z['p']),
                    'coef':float(z['coef'])})
r=pd.DataFrame(out)
for label in r.model.unique():
    mask=r.model.eq(label)
    r.loc[mask,'Holm_two_candidates']=multipletests(r.loc[mask,'p'],method='holm')[1]
r.to_csv(ROOT/'results/TCGA_NDST3_ABCA5_validation.csv',index=False)
print('eligible',len(eligible),'events',int(eligible.event.sum()))
print(r.to_string(index=False))
