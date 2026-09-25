import os
"""Cross-cohort MCP-counter and conditional rank-association checks."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr,rankdata
from sklearn.decomposition import PCA

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/virtual_validation'
BASE=pd.read_csv(Path(os.environ['JOB46_ROOT'])/'reanalysis_v2/results/all_cohort_scores.csv',low_memory=False)
GEO=pd.read_csv(ROOT/'results/geo_external_replication/sample_scores.csv',low_memory=False)
MAP=json.loads(Path(os.environ['GEO_LOCAL_ROOT'])/'GSE16011_gsm_to_dbid.json'.read_text())
CONTROLS=['T cells','Endothelial cells','Fibroblasts','Neutrophils']

def residual_rank(y,X):
    y=rankdata(y)
    X=np.column_stack([np.ones(len(y)),*[rankdata(X[c]) for c in X]])
    return y-X@np.linalg.lstsq(X,y,rcond=None)[0]

rows=[]; controls=[]; axes=[]
for cohort in ['TCGA','CGGA693','CGGA325','GSE16011','GSE149009']:
    m=pd.read_csv(OUT/f'{cohort}_MCPcounter.csv')
    if cohort.startswith('GSE'):
        z=GEO[GEO.cohort.eq(cohort)].copy()
        if cohort=='GSE16011':
            m['database_id']=m['sample'].map(MAP)
            z=m.merge(z[['database_id','xenobiotic','eligibility']],on='database_id',validate='one_to_one')
        else:z=m.merge(z[['sample','xenobiotic','eligibility']],on='sample',validate='one_to_one')
        x='xenobiotic'
    else:
        z=m.merge(BASE[(BASE.cohort.eq(cohort))&BASE.eligible_GBM.eq(True)][['sample','HALLMARK_XENOBIOTIC_METABOLISM','myeloid_identity','hypoxia']],on='sample',validate='one_to_one')
        x='HALLMARK_XENOBIOTIC_METABOLISM'
    for proxy in [c for c in m.columns if c!='sample']:
        if proxy=='database_id':continue
        q=z[[x,proxy]].apply(pd.to_numeric,errors='coerce').dropna()
        if q[proxy].nunique()<3:continue
        rho,p=spearmanr(q[x],q[proxy])
        rows.append(dict(cohort=cohort,proxy=proxy,n=len(q),rho=rho,p=p))
    q=z[[x,'Monocytic lineage',*CONTROLS]].apply(pd.to_numeric,errors='coerce').dropna()
    rx=residual_rank(q[x],q[CONTROLS]);ry=residual_rank(q['Monocytic lineage'],q[CONTROLS])
    rho,_=spearmanr(rx,ry)
    rng=np.random.default_rng(49017)
    null=np.array([spearmanr(rx,rng.permutation(ry)).statistic for _ in range(5000)])
    p_perm=(1+np.sum(np.abs(null)>=abs(rho)))/(1+len(null))
    controls.append(dict(cohort=cohort,n=len(q),adjusted_for='T cells+endothelial+fibroblasts+neutrophils',partial_rank_rho=rho,permutation_p=p_perm))
    cellcols=[c for c in m.columns if c!='sample' and c!='database_id']
    pc=z[[x,*cellcols]].apply(pd.to_numeric,errors='coerce').dropna()
    values=pc[cellcols]
    values=(values-values.mean())/values.std(ddof=0)
    fit=PCA(n_components=1).fit(values)
    pc1=fit.transform(values)[:,0]
    if np.corrcoef(pc1,values['Monocytic lineage'])[0,1]<0:pc1=-pc1
    r,p=spearmanr(pc[x],pc1)
    axes.append(dict(cohort=cohort,n=len(pc),PC1_variance_fraction=fit.explained_variance_ratio_[0],xeno_PC1_rho=r,p=p,
                     monocytic_PC1_loading=np.corrcoef(pc1,values['Monocytic lineage'])[0,1],
                     fibroblast_PC1_loading=np.corrcoef(pc1,values['Fibroblasts'])[0,1]))

pd.DataFrame(rows).to_csv(OUT/'MCPcounter_correlations.csv',index=False)
pd.DataFrame(controls).to_csv(OUT/'MCPcounter_partial_rank.csv',index=False)
pd.DataFrame(axes).to_csv(OUT/'MCPcounter_PCA.csv',index=False)
print(pd.DataFrame(rows).query('proxy in ["Monocytic lineage","T cells","Endothelial cells","Fibroblasts"]')[['cohort','proxy','n','rho','p']].to_string(index=False))
print(pd.DataFrame(controls).to_string(index=False))
print(pd.DataFrame(axes).to_string(index=False))
