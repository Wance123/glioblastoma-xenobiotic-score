import os
"""New exploratory pathway analysis across TCGA and two CGGA cohorts.

Six fixed Hallmark pathway rank scores from the audited job46 bulk preparation
are tested in IDH-WT GBM and adult IDH-WT diffuse glioma (grade adjusted).
No outcome-driven gene selection or threshold optimization occurs.
"""
from pathlib import Path
import numpy as np,pandas as pd
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[1]
SRC=Path(os.environ['JOB46_ROOT'])/'reanalysis_v2/results/all_cohort_scores.csv'
NAMES=['GLYCOLYSIS','FATTY_ACID_METABOLISM','OXIDATIVE_PHOSPHORYLATION',
       'XENOBIOTIC_METABOLISM','CHOLESTEROL_HOMEOSTASIS','BILE_ACID_METABOLISM']
d=pd.read_csv(SRC,low_memory=False)
rows=[]; audit=[]
for scope in ['GBM','IDHwt_all_grade']:
    for cohort in ['TCGA','CGGA693','CGGA325']:
        z=d[d.cohort.eq(cohort)&d.primary.eq(True)&d.IDH.eq(0)&d.age.ge(18)&d.time.gt(0)&d.event.isin([0,1])].copy()
        if scope=='GBM':z=z[z.GBM.eq(True)]
        else:z=z[z.grade.isin([2,3,4])]
        audit.append({'scope':scope,'cohort':cohort,'n':len(z),'events':int(z.event.sum()),
            'grades':';'.join(f'{int(k)}:{int(v)}' for k,v in z.grade.value_counts().sort_index().items())})
        for name in NAMES:
            col='HALLMARK_'+name
            score=pd.to_numeric(z[col],errors='coerce')
            sd=score.std(ddof=0)
            if not np.isfinite(sd) or sd==0:continue
            z['score_z']=(score-score.mean())/sd
            covars=['age','MGMT']+(['grade','IDH'] if scope=='IDHwt_all_grade' else [])
            # IDH is fixed at WT here, so omit it to avoid singularity.
            covars=[v for v in covars if v!='IDH']
            for model,covs in [('age',['age']+(['grade'] if scope=='IDHwt_all_grade' else [])),('age_MGMT',covars)]:
                f=z[['time','event','score_z']+covs].dropna()
                if len(f)<30:continue
                fit=CoxPHFitter(penalizer=.01).fit(f,'time','event')
                q=fit.summary.loc['score_z']
                rows.append({'scope':scope,'cohort':cohort,'pathway':name,'model':model,'n':len(f),'events':int(f.event.sum()),
                    'coef':float(q['coef']),'HR_per_SD':float(np.exp(q['coef'])),
                    'CI_low':float(np.exp(q['coef lower 95%'])),'CI_high':float(np.exp(q['coef upper 95%'])),
                    'p':float(q['p'])})
r=pd.DataFrame(rows)
r['BH_FDR_six_per_cohort_scope_model']=r.groupby(['scope','cohort','model']).p.transform(lambda p:multipletests(p,method='fdr_bh')[1])
r.to_csv(ROOT/'results/pathway_cox_all.csv',index=False)
pd.DataFrame(audit).to_csv(ROOT/'results/pathway_cohort_audit.csv',index=False)
wide=r[(r.scope=='GBM')&(r.model=='age_MGMT')].pivot(index='pathway',columns='cohort',values=['coef','p','BH_FDR_six_per_cohort_scope_model'])
wide.columns=['_'.join(c) for c in wide.columns];wide=wide.reset_index()
wide['same_direction_all']=((np.sign(wide.coef_TCGA)==np.sign(wide.coef_CGGA693))&(np.sign(wide.coef_TCGA)==np.sign(wide.coef_CGGA325)))
wide['nominal_all']=wide.same_direction_all&(wide.p_TCGA<.05)&(wide.p_CGGA693<.05)&(wide.p_CGGA325<.05)
wide['FDR_all']=wide.same_direction_all&(wide.BH_FDR_six_per_cohort_scope_model_TCGA<.05)&(wide.BH_FDR_six_per_cohort_scope_model_CGGA693<.05)&(wide.BH_FDR_six_per_cohort_scope_model_CGGA325<.05)
wide.to_csv(ROOT/'results/pathway_GBM_three_cohort_replication.csv',index=False)
print(pd.DataFrame(audit).to_string(index=False))
print(wide.to_string(index=False))
