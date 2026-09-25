import os
"""Patient-aware spatial xenobiotic and myeloid-marker association.

Spatial RNA co-localization only; existing job46 spot QC and marker scores are
reused, while xenobiotic score is freshly calculated from original matrices.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,re,itertools
import numpy as np,pandas as pd,h5py
from scipy import sparse,io,stats

ROOT=Path(__file__).resolve().parents[1]
J46=Path(os.environ['JOB46_ROOT'])
LED=pd.read_csv(J46/'results/audit/sample_ledger.csv')
H='HALLMARK_XENOBIOTIC_METABOLISM'
GENES=set(json.loads((ROOT/'input'/(H+'.json')).read_text())[H]['geneSymbols'])
OUT=ROOT/'results/spatial_xenobiotic';OUT.mkdir(exist_ok=True)
def decode(v):return np.array([a.decode() if isinstance(a,bytes) else str(a) for a in v])
def read_h5(p):
    with h5py.File(p) as f:
        g=f['matrix'];x=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g['shape'][:]))
        return x,decode(g['features']['name'][:]),decode(g['barcodes'][:])
def file_for(r):
    base=Path(os.environ['SPATIAL_RAW_ROOT'])/r['dataset']/'raw'
    if r['dataset']=='GSE237183':return next(base.glob(r['gsm']+'*_filtered_feature_bc_matrix.h5'))
    title=r['title'];matches=[]
    norm=lambda s:re.sub(r'(?<=-)0(?=\d)','',s)
    for p in base.glob('*_filtered_feature_bc_matrix'):
        a=p.name.replace('_filtered_feature_bc_matrix','').split('-',1)[1]
        if norm(a)==norm(title):matches.append(p)
    assert len(matches)==1,(title,matches)
    return matches[0]
def partial_rank(x,y,z):
    x=stats.rankdata(x);y=stats.rankdata(y)
    c=np.column_stack([np.ones(len(x))]+[stats.rankdata(z[:,j]) for j in range(z.shape[1])])
    a=x-c@np.linalg.lstsq(c,x,rcond=None)[0];b=y-c@np.linalg.lstsq(c,y,rcond=None)[0]
    den=np.sqrt((a@a)*(b@b));return float((a@b)/den) if den>0 else np.nan
rows=[];coverage=[];nullrows=[];RNG=np.random.default_rng(490912)
for r in LED[(LED.include==True)&LED.dataset.isin(['GSE237183','GSE318562'])].to_dict('records'):
    p=file_for(r);stem=p.name.replace('_filtered_feature_bc_matrix.h5','').replace('_filtered_feature_bc_matrix','')
    spot=pd.read_csv(J46/'results/spatial'/f'{r["dataset"]}_{stem}_scores.csv.gz',usecols=['barcode','macrophage_state','myeloid_identity','microglia_state','tumor_glial_proxy','hypoxia','log_counts'])
    if p.is_file():x,genes,bars=read_h5(p)
    else:
        x=io.mmread(str(p/'matrix.mtx.gz')).tocsc();genes=pd.read_csv(p/'features.tsv.gz',sep='\t',header=None)[1].astype(str).to_numpy();bars=pd.read_csv(p/'barcodes.tsv.gz',sep='\t',header=None)[0].astype(str).to_numpy()
    assert x.shape==(len(genes),len(bars))
    use=[i for i,g in enumerate(genes) if g in GENES]
    assert len(use)>=100,(stem,len(use))
    counts=np.asarray(x.sum(axis=0)).ravel();keep=counts>0
    v=x[use][:,keep].astype(float)@sparse.diags(1e4/counts[keep]);v.data=np.log1p(v.data)
    a=v.toarray().T;sd=a.std(axis=0);a=(a-a.mean(axis=0))/np.where(sd>1e-10,sd,1)
    xx=pd.DataFrame({'barcode':bars[keep],'xeno':a.mean(axis=1)})
    frame=spot.merge(xx,on='barcode',validate='one_to_one',how='inner')
    assert len(frame)==len(spot),(stem,len(frame),len(spot))
    coverage.append({'dataset':r['dataset'],'patient':r['patient'],'section':stem,'genes_defined':len(GENES),'genes_measured':len(use),'n_spots':len(frame)})
    for outcome in ['macrophage_state','myeloid_identity','microglia_state']:
        rho=stats.spearmanr(frame.xeno,frame[outcome]).statistic
        adj=partial_rank(frame.xeno.to_numpy(),frame[outcome].to_numpy(),frame[['hypoxia','tumor_glial_proxy','log_counts']].to_numpy())
        rows.append({'dataset':r['dataset'],'patient':r['patient'],'section':stem,'outcome':outcome,'n_spots':len(frame),'rho_unadjusted':rho,'rho_adjusted':adj})
    # Matched-gene negative controls: choose genes close in mean normalized
    # expression and detection, excluding every target Hallmark gene.
    used_norm=(x[:,keep].astype(float)@sparse.diags(1e4/counts[keep])).tocsr()
    mean=np.asarray(used_norm.mean(axis=1)).ravel()
    det=np.asarray((x[:,keep]>0).mean(axis=1)).ravel()
    features=np.column_stack([np.log1p(mean),det]);features=(features-features.mean(axis=0))/np.maximum(features.std(axis=0),1e-8)
    eligible=np.array([i for i,g in enumerate(genes) if g not in GENES and det[i]>.01 and mean[i]>0])
    pools=[]
    for i in use:
        order=np.argsort(np.sum((features[eligible]-features[i])**2,axis=1))
        pools.append(eligible[order[:100]])
    barcode_pos=pd.Index(bars[keep]).get_indexer(frame.barcode)
    assert (barcode_pos>=0).all()
    for null_id in range(50):
        chosen=[]
        for pool in pools:
            candidates=np.setdiff1d(pool,chosen,assume_unique=False)
            chosen.append(int(RNG.choice(candidates if len(candidates) else pool)))
        vn=used_norm[chosen].copy();vn.data=np.log1p(vn.data)
        aa=vn.toarray().T[barcode_pos]
        ss=aa.std(axis=0);aa=(aa-aa.mean(axis=0))/np.where(ss>1e-10,ss,1)
        score=aa.mean(axis=1)
        nullrows.append({'dataset':r['dataset'],'patient':r['patient'],'section':stem,'null_id':null_id,
            'rho_adjusted':partial_rank(score,frame.myeloid_identity.to_numpy(),frame[['hypoxia','tumor_glial_proxy','log_counts']].to_numpy())})
    print(r['dataset'],r['patient'],stem,len(frame),flush=True)
pd.DataFrame(rows).to_csv(OUT/'section_effects.csv',index=False)
pd.DataFrame(coverage).to_csv(OUT/'coverage.csv',index=False)
pd.DataFrame(nullrows).to_csv(OUT/'matched_gene_null_section.csv',index=False)
sections=pd.DataFrame(rows)
pat=sections.groupby(['dataset','patient','outcome'],as_index=False).agg(n_sections=('section','size'),rho_unadjusted=('rho_unadjusted','median'),rho_adjusted=('rho_adjusted','median'))
pat.to_csv(OUT/'patient_effects.csv',index=False)
summ=[]
for (dataset,outcome),g in pat.groupby(['dataset','outcome']):
    for model in ['rho_unadjusted','rho_adjusted']:
        vals=g[model].to_numpy();p=stats.binomtest(int((vals>0).sum()),len(vals),.5,alternative='two-sided').pvalue
        summ.append({'dataset':dataset,'outcome':outcome,'model':model,'n_patients':len(vals),'n_positive':int((vals>0).sum()),'median_rho':float(np.median(vals)),'patient_sign_p':p})
pd.DataFrame(summ).to_csv(OUT/'cohort_summary.csv',index=False)
print(pd.DataFrame(summ).to_string(index=False))
obs=pat[pat.outcome=='myeloid_identity'].groupby('dataset').rho_adjusted.median()
nr=pd.DataFrame(nullrows).groupby(['dataset','null_id','patient'],as_index=False).rho_adjusted.median().groupby(['dataset','null_id'],as_index=False).rho_adjusted.median()
ns=[]
for dataset,q in nr.groupby('dataset'):
    actual=float(obs[dataset]);p=(1+int((q.rho_adjusted.abs()>=abs(actual)).sum()))/(1+len(q))
    ns.append({'dataset':dataset,'observed_patient_median_rho':actual,'null_panels':len(q),'null_median_rho':float(q.rho_adjusted.median()),'empirical_two_sided_p':p})
pd.DataFrame(ns).to_csv(OUT/'matched_gene_null_summary.csv',index=False)
print(pd.DataFrame(ns).to_string(index=False))
