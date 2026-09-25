import os
"""Patient-paired lineage comparison of the xenobiotic Hallmark program."""
from pathlib import Path
import gzip,json
import numpy as np,pandas as pd
from scipy.stats import wilcoxon

ROOT=Path(__file__).resolve().parents[1]
DATA=Path(os.environ['GSE131928_ROOT'])
ANN=Path(os.environ['JOB46_ROOT'])/'results/single_cell'
name='HALLMARK_XENOBIOTIC_METABOLISM'
genes=set(json.loads((ROOT/'input'/(name+'.json')).read_text())[name]['geneSymbols'])
rows=[];aud=[];tests=[]
for platform,pattern in [('Smartseq2','*Smartseq2*processed_TPM.tsv.gz'),('10X','*10X*processed_TPM.tsv.gz')]:
    path=next(DATA.glob(pattern))
    with gzip.open(path,'rt') as f:
        barcodes=f.readline().rstrip('\n').split('\t')[1:]
        sums=np.zeros(len(barcodes));found=[]
        for line in f:
            g=line.split('\t',1)[0]
            if g in genes:
                v=np.fromstring(line.split('\t',1)[1],sep='\t')
                assert len(v)==len(barcodes)
                sums+=np.log1p(v);found.append(g)
    scores=pd.DataFrame({'barcode':barcodes,'score':sums/len(found)})
    meta=pd.read_csv(ANN/f'GSE131928_{platform}_scores.csv.gz',usecols=['barcode','patient','lineage'])
    d=meta.merge(scores,on='barcode',validate='one_to_one')
    aud.append({'platform':platform,'genes_set':len(genes),'genes_measured':len(found),'cells_joined':len(d),'patients':d.patient.nunique()})
    for (patient,lineage),g in d.groupby(['patient','lineage']):
        if len(g)>=5:rows.append({'platform':platform,'patient':patient,'lineage':lineage,'n_cells':len(g),'mean_score':g.score.mean()})
r=pd.DataFrame(rows);r.to_csv(ROOT/'results/GSE131928_xenobiotic_patient_lineage.csv',index=False)
for platform in ['Smartseq2','10X']:
    p=r[r.platform.eq(platform)].pivot(index='patient',columns='lineage',values='mean_score')
    if {'Myeloid','Glial_tumor_like'}.issubset(p.columns):
        q=p[['Myeloid','Glial_tumor_like']].dropna()
        stat=wilcoxon(q.Myeloid-q.Glial_tumor_like) if len(q)>=5 else None
        tests.append({'platform':platform,'paired_patients':len(q),'myeloid_minus_glial_median':float((q.Myeloid-q.Glial_tumor_like).median()),
            'positive_patients':int((q.Myeloid>q.Glial_tumor_like).sum()),'wilcoxon_p':None if stat is None else float(stat.pvalue)})
pd.DataFrame(aud).to_csv(ROOT/'results/GSE131928_xenobiotic_coverage.csv',index=False)
pd.DataFrame(tests).to_csv(ROOT/'results/GSE131928_xenobiotic_paired_tests.csv',index=False)
print(pd.DataFrame(aud).to_string(index=False));print(pd.DataFrame(tests).to_string(index=False))
