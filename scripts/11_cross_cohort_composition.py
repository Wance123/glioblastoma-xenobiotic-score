import os
"""Cross-cohort association of xenobiotic score with independent marker panels.

Marker panels were prepared in job46; remove overlapping Hallmark genes from
the original panel definitions before interpreting this as independent evidence.
"""
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from lifelines import CoxPHFitter
ROOT=Path(__file__).resolve().parents[1]
SRC=Path(os.environ['JOB46_ROOT'])/'reanalysis_v2/results/all_cohort_scores.csv'
d=pd.read_csv(SRC,low_memory=False)
rows=[]; models=[]
for cohort in ['TCGA','CGGA693','CGGA325']:
    z=d[(d.cohort==cohort)&d.eligible_GBM.eq(True)].copy()
    score='HALLMARK_XENOBIOTIC_METABOLISM'
    for proxy in ['myeloid_identity','macrophage_state','microglia_state','tumor_glial_proxy','mesenchymal','hypoxia']:
        a=z[[score,proxy]].apply(pd.to_numeric,errors='coerce').dropna()
        rho,p=spearmanr(a[score],a[proxy])
        rows.append({'cohort':cohort,'proxy':proxy,'n':len(a),'rho':rho,'p':p})
    for proxy in ['myeloid_identity','tumor_glial_proxy','hypoxia']:
        a=z[['time','event','age','MGMT',score,proxy]].apply(pd.to_numeric,errors='coerce').dropna()
        for v in [score,proxy]:a[v]=(a[v]-a[v].mean())/a[v].std(ddof=0)
        a=a.rename(columns={score:'xeno_score',proxy:'proxy_score'})
        fit=CoxPHFitter(penalizer=.01).fit(a,'time','event')
        q=fit.summary.loc['xeno_score']
        models.append({'cohort':cohort,'adjusted_for':proxy,'n':len(a),'events':int(a.event.sum()),
                       'HR_xeno_per_SD':float(np.exp(q['coef'])),'p':float(q['p'])})
pd.DataFrame(rows).to_csv(ROOT/'results/xenobiotic_marker_correlations_all_cohorts.csv',index=False)
pd.DataFrame(models).to_csv(ROOT/'results/xenobiotic_marker_adjusted_cox.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(models).to_string(index=False))
