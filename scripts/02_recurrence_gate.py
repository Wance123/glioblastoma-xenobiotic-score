import os
"""Exploratory cross-sectional primary versus recurrent IDH-WT GBM comparison.

These are unpaired specimens; no patient-matched recurrence claim is possible.
"""
from pathlib import Path
import zipfile
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[1]
LOCAL=Path(os.environ['GEO_LOCAL_ROOT'])
GENES=["MPI","PMM2","GMPPA","GMPPB"]
COHORTS={
"CGGA693":(ROOT/"input/CGGA.mRNAseq_693_clinical.20200506.txt.zip",ROOT/"input/CGGA.mRNAseq_693.RSEM-genes.20200506.txt.zip"),
"CGGA325":(LOCAL/"CGGA.mRNAseq_325_clinical.20200506.txt.zip",LOCAL/"CGGA.mRNAseq_325.RSEM-genes.20200506.txt.zip")}
def load(path):
    with zipfile.ZipFile(path) as z:
        return pd.read_csv(z.open(next(n for n in z.namelist() if n.endswith('.txt') and not n.startswith('__MACOSX'))),sep='\t',low_memory=False)
out=[]
for cohort,(cp,ep) in COHORTS.items():
    c=load(cp); e=load(ep)
    c['Age']=pd.to_numeric(c.Age,errors='coerce')
    c=c[(c.PRS_type.isin(['Primary','Recurrent'])) & (c.Histology.isin(['GBM','rGBM'])) &
        (c.IDH_mutation_status=='Wildtype') & (c.Age>=18)].copy()
    c['recurrent']=(c.PRS_type=='Recurrent').astype(int)
    x=e[e.Gene_Name.isin(GENES)].set_index('Gene_Name').T
    c=c.set_index('CGGA_ID').join(x,how='inner')
    print(cohort,'primary',sum(c.recurrent==0),'recurrent',sum(c.recurrent==1))
    for gene in GENES:
        c['logexpr']=np.log2(pd.to_numeric(c[gene],errors='coerce')+1)
        m=smf.ols('logexpr ~ recurrent + Age',data=c).fit(cov_type='HC3')
        out.append({'cohort':cohort,'gene':gene,'n':int(m.nobs),'n_recurrent':int(c.recurrent.sum()),
            'recurrent_minus_primary_log2':float(m.params['recurrent']),
            'CI_low':float(m.conf_int().loc['recurrent',0]),'CI_high':float(m.conf_int().loc['recurrent',1]),
            'p':float(m.pvalues['recurrent'])})
d=pd.DataFrame(out)
d['BH_FDR_within_cohort']=d.groupby('cohort').p.transform(lambda x:multipletests(x,method='fdr_bh')[1])
d.to_csv(ROOT/'results/recurrence_unpaired_gene_models.csv',index=False)
print(d.to_string(index=False))
