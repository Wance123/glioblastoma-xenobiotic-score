import os
"""Patient-level descriptive cellular-source audit for NDST3 and ABCA5.

Uses original GSE131928 processed TPM, joined to pre-existing job46 heuristic
lineage labels by barcode. No cells are treated as independent patients.
"""
from pathlib import Path
import gzip
import numpy as np,pandas as pd

ROOT=Path(__file__).resolve().parents[1]
DATA=Path(os.environ['GSE131928_ROOT'])
ANN=Path(os.environ['JOB46_ROOT'])/'results/single_cell'
all_rows=[];coverage=[]
for platform,pattern in [('Smartseq2','*Smartseq2*processed_TPM.tsv.gz'),('10X','*10X*processed_TPM.tsv.gz')]:
    path=next(DATA.glob(pattern))
    meta=pd.read_csv(ANN/f'GSE131928_{platform}_scores.csv.gz',usecols=['barcode','patient','lineage'])
    assert meta.barcode.is_unique
    with gzip.open(path,'rt') as f:
        head=f.readline().rstrip('\n').split('\t')
        targets={}
        for line in f:
            g=line.split('\t',1)[0]
            if g in ['NDST3','ABCA5']:
                v=np.fromstring(line.split('\t',1)[1],sep='\t')
                assert len(v)==len(head)-1,(platform,g,len(v),len(head)-1)
                targets[g]=v
                if len(targets)==2:break
    assert set(targets)=={'NDST3','ABCA5'}
    x=pd.DataFrame({'barcode':head[1:],'NDST3':targets['NDST3'],'ABCA5':targets['ABCA5']})
    d=meta.merge(x,on='barcode',how='inner',validate='one_to_one')
    coverage.append({'platform':platform,'matrix_cells':len(x),'annotated_adult_cells':len(meta),'joined_cells':len(d),'patients':d.patient.nunique()})
    for (patient,lineage),g in d.groupby(['patient','lineage']):
        if len(g)<5:continue
        for gene in ['NDST3','ABCA5']:
            all_rows.append({'platform':platform,'patient':patient,'lineage':lineage,'gene':gene,
                'n_cells':len(g),'detection_fraction':float((g[gene]>0).mean()),
                'mean_log1p_TPM':float(np.log1p(g[gene]).mean())})
res=pd.DataFrame(all_rows)
res.to_csv(ROOT/'results/GSE131928_patient_lineage_gene_expression.csv',index=False)
pd.DataFrame(coverage).to_csv(ROOT/'results/GSE131928_origin_coverage.csv',index=False)
summary=res.groupby(['platform','lineage','gene']).agg(patients=('patient','nunique'),
    median_detection=('detection_fraction','median'),median_mean_log1p_TPM=('mean_log1p_TPM','median')).reset_index()
summary.to_csv(ROOT/'results/GSE131928_cell_origin_summary.csv',index=False)
print(pd.DataFrame(coverage).to_string(index=False));print(summary.to_string(index=False))
