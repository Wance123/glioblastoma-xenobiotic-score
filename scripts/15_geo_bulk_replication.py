import os
"""External GEO composition replication with frozen xenobiotic and marker panels."""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from lifelines import CoxPHFitter

ROOT=Path(__file__).resolve().parents[1]
DATA=Path(os.environ['GEO_LOCAL_ROOT'])
OUT=ROOT/'results'/'geo_external_replication'
OUT.mkdir(parents=True,exist_ok=True)
genes=set(json.loads((ROOT/'input/HALLMARK_XENOBIOTIC_METABOLISM.json').read_text())['HALLMARK_XENOBIOTIC_METABOLISM']['geneSymbols'])
gs=pd.read_csv(ROOT/'config/gene_sets.csv')
panels={m:set(gs.loc[gs.module.eq(m),'gene']) for m in ['myeloid_identity','macrophage_state','tumor_glial_proxy']}
assert all(not genes.intersection(v) for v in panels.values())

def scored(path):
    x=pd.read_csv(path,sep='\t',compression='gzip',index_col=0)
    x=x.loc[:,~x.columns.duplicated()]
    ranks=x.rank(axis=1,pct=True)
    scores={}
    coverage={}
    for name,panel in [('xenobiotic',genes),*panels.items()]:
        present=sorted(panel.intersection(x.columns))
        coverage[name]=f'{len(present)}/{len(panel)}'
        scores[name]=ranks[present].mean(axis=1)
    return pd.DataFrame(scores),coverage

rows=[]; samples=[]; survival=[]; hashes=[]
for accession,filename in [('GSE16011','GSE16011_genesymbol.txt.gz'),('GSE149009','GSE149009_genesymbol.txt.gz')]:
    p=DATA/filename
    scores,cov=scored(p)
    if accession=='GSE16011':
        mapping=json.loads((DATA/'GSE16011_gsm_to_dbid.json').read_text())
        clin=pd.read_csv(DATA/'GSE16011_clinical.csv')
        clin=clin[clin.histology.eq('GBM (grade IV)')&clin.age.ge(18)].copy()
        scores['database_id']=[mapping.get(s,np.nan) for s in scores.index]
        scores=scores.merge(clin,on='database_id',how='inner')
        scores['sample']=scores['expr_label']
        scores['cohort']=accession
        scores['eligibility']='adult historical GBM; IDH unknown'
    else:
        scores['sample']=scores.index
        scores['cohort']=accession
        scores['eligibility']='GBM tumor tissue; IDH not uniformly annotated'
    for proxy in panels:
        z=scores[['xenobiotic',proxy]].dropna()
        rho,pv=spearmanr(z.xenobiotic,z[proxy])
        rows.append(dict(cohort=accession,proxy=proxy,n=len(z),rho=rho,p=pv,
                         xeno_coverage=cov['xenobiotic'],proxy_coverage=cov[proxy]))
    if accession=='GSE16011':
        for adjust in [['age'],['age','gender']]:
            z=scores[['OS_days','OS_status','xenobiotic',*adjust]].copy()
            if 'gender' in z:z['gender']=z.gender.map({'Male':1,'Female':0})
            z=z.apply(pd.to_numeric,errors='coerce').dropna()
            z=z[z.OS_days.gt(0)&z.OS_status.isin([0,1])]
            z['xenobiotic']=(z.xenobiotic-z.xenobiotic.mean())/z.xenobiotic.std(ddof=0)
            fit=CoxPHFitter(penalizer=.01).fit(z,duration_col='OS_days',event_col='OS_status')
            q=fit.summary.loc['xenobiotic']
            survival.append(dict(cohort=accession,adjustment='+'.join(adjust),n=len(z),events=int(z.OS_status.sum()),
                                 HR_per_SD=float(np.exp(q['coef'])),CI_low=float(np.exp(q['coef lower 95%'])),
                                 CI_high=float(np.exp(q['coef upper 95%'])),p=float(q['p'])))
    samples.append(scores)
    for f in [p]+([DATA/'GSE16011_clinical.csv',DATA/'GSE16011_gsm_to_dbid.json'] if accession=='GSE16011' else []):
        hashes.append(dict(file=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))

pd.DataFrame(rows).to_csv(OUT/'correlations.csv',index=False)
pd.DataFrame(survival).to_csv(OUT/'survival.csv',index=False)
pd.concat(samples,ignore_index=True).to_csv(OUT/'sample_scores.csv',index=False)
pd.DataFrame(hashes).to_csv(OUT/'input_sha256.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
print(pd.DataFrame(survival).to_string(index=False))
