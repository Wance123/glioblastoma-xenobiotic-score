import os
"""Exploratory metabolism-wide screen with discovery/replication separation.

Gene universe: union of six official MSigDB Hallmark sets, frozen to JSON in
input before outcomes are read. All tests and non-hits are written to CSV.
"""
from pathlib import Path
import json, zipfile, requests
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

ROOT=Path(__file__).resolve().parents[1]; IN=ROOT/'input'; OUT=ROOT/'results'
LOCAL=Path(os.environ['GEO_LOCAL_ROOT'])
SETS=['HALLMARK_GLYCOLYSIS','HALLMARK_FATTY_ACID_METABOLISM',
'HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_XENOBIOTIC_METABOLISM',
'HALLMARK_CHOLESTEROL_HOMEOSTASIS','HALLMARK_BILE_ACID_METABOLISM']
COHORTS={'CGGA693':(IN/'CGGA.mRNAseq_693_clinical.20200506.txt.zip',IN/'CGGA.mRNAseq_693.RSEM-genes.20200506.txt.zip'),
'CGGA325':(LOCAL/'CGGA.mRNAseq_325_clinical.20200506.txt.zip',LOCAL/'CGGA.mRNAseq_325.RSEM-genes.20200506.txt.zip')}
def load(p):
    with zipfile.ZipFile(p) as z:
        return pd.read_csv(z.open(next(n for n in z.namelist() if n.endswith('.txt') and not n.startswith('__MACOSX'))),sep='\t',low_memory=False)
membership={}
for name in SETS:
    file=IN/(name+'.json')
    if not file.exists():
        resp=requests.get('https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/'+name+'.json',timeout=30)
        resp.raise_for_status(); json.loads(resp.content)[name]; file.write_bytes(resp.content)
    for gene in json.loads(file.read_text())[name]['geneSymbols']:
        membership.setdefault(gene,[]).append(name)
pd.DataFrame([{'gene':g,'sets':';'.join(v)} for g,v in sorted(membership.items())]).to_csv(OUT/'metabolic_gene_universe.csv',index=False)
rows=[]; audit=[]
for cohort,(cp,ep) in COHORTS.items():
    c=load(cp); x=load(ep)
    c['Age']=pd.to_numeric(c.Age,errors='coerce')
    c['time']=pd.to_numeric(c.OS,errors='coerce')
    c['event']=pd.to_numeric(c['Censor (alive=0; dead=1)'],errors='coerce')
    c['MGMT']=c.MGMTp_methylation_status.map({'methylated':1,'un-methylated':0})
    c['TMZ']=pd.to_numeric(c['Chemo_status (TMZ treated=1;un-treated=0)'],errors='coerce')
    c['RT']=pd.to_numeric(c['Radio_status (treated=1;un-treated=0)'],errors='coerce')
    c=c[(c.PRS_type=='Primary')&(c.Histology=='GBM')&(c.IDH_mutation_status=='Wildtype')&(c.Age>=18)&(c.time>0)&c.event.isin([0,1])].set_index('CGGA_ID')
    x=x[x.Gene_Name.isin(membership)].drop_duplicates('Gene_Name').set_index('Gene_Name').T
    x=x.apply(pd.to_numeric,errors='coerce')
    common=c.index.intersection(x.index); c=c.loc[common]; x=x.loc[common]
    audit.append({'cohort':cohort,'n':len(c),'events':int(c.event.sum()),'genes_universe':len(membership),'genes_detected':x.shape[1]})
    base=c[['time','event','Age','MGMT','TMZ','RT']].dropna()
    for gene in sorted(membership):
        if gene not in x: continue
        v=np.log2(x.loc[base.index,gene]+1)
        if v.notna().sum()!=len(v) or v.std(ddof=0)<0.1 or (x.loc[base.index,gene]>0).mean()<0.5: continue
        v=(v-v.mean())/v.std(ddof=0)
        d=base.copy(); d['gene_score']=v
        try:
            fit=CoxPHFitter(penalizer=0.01).fit(d,'time','event')
            z=fit.summary.loc['gene_score']
            rows.append({'cohort':cohort,'gene':gene,'n':len(d),'events':int(d.event.sum()),'coef':float(z['coef']),'HR_per_SD':float(np.exp(z['coef'])),'p':float(z['p'])})
        except Exception:
            pass
d=pd.DataFrame(rows)
d['BH_FDR_within_cohort']=d.groupby('cohort').p.transform(lambda v:multipletests(v,method='fdr_bh')[1])
d.to_csv(OUT/'metabolic_survival_all_genes.csv',index=False)
pd.DataFrame(audit).to_csv(OUT/'metabolic_screen_audit.csv',index=False)
p=d.pivot(index='gene',columns='cohort',values=['coef','p','BH_FDR_within_cohort'])
p.columns=['_'.join(col) for col in p.columns]; p=p.reset_index()
p['same_direction']=np.sign(p.coef_CGGA693)==np.sign(p.coef_CGGA325)
p['replicated_nominal']=(p.p_CGGA693<.05)&(p.p_CGGA325<.05)&p.same_direction
p['replicated_discovery_FDR']=(p.BH_FDR_within_cohort_CGGA693<.05)&(p.p_CGGA325<.05)&p.same_direction
p.sort_values('p_CGGA693').to_csv(OUT/'metabolic_survival_replication.csv',index=False)
print('universe',len(membership));print(pd.DataFrame(audit).to_string(index=False))
print('nominal replicated',int(p.replicated_nominal.sum()),'discovery FDR + external nominal',int(p.replicated_discovery_FDR.sum()))
print(p[p.replicated_nominal].sort_values('p_CGGA693').head(30).to_string(index=False))
