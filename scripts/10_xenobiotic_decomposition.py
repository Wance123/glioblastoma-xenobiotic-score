import os
"""Functional decomposition of the 200-gene Hallmark xenobiotic program.

MSigDB GO/Reactome intersections are fixed before outcome analysis. Modules
overlap biologically; all tested groups and negative findings are retained.
"""
from pathlib import Path
import json
import numpy as np,pandas as pd
from lifelines import CoxPHFitter
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[1]
SRC=Path(os.environ['JOB46_ROOT'])/'reanalysis_v2'
H='HALLMARK_XENOBIOTIC_METABOLISM'
NAMES={
 'phase_I':'REACTOME_PHASE_I_FUNCTIONALIZATION_OF_COMPOUNDS',
 'phase_II':'REACTOME_PHASE_II_CONJUGATION_OF_COMPOUNDS',
 'transport':'REACTOME_TRANSPORT_OF_SMALL_MOLECULES',
 'GO_xenobiotic':'GOBP_XENOBIOTIC_METABOLIC_PROCESS',
 'detoxification':'GOBP_DETOXIFICATION',
 'oxidative_stress':'GOBP_RESPONSE_TO_OXIDATIVE_STRESS',
 'immune_annotated':'GOBP_IMMUNE_RESPONSE',
}
def genes(name):return set(json.loads((ROOT/'input'/(name+'.json')).read_text())[name]['geneSymbols'])
hall=genes(H)
modules={'whole_hallmark':hall}
modules.update({k:hall&genes(v) for k,v in NAMES.items()})
modules['immune_excluded']=hall-modules['immune_annotated']
modules['enzyme_core']=(modules['phase_I']|modules['phase_II']|modules['GO_xenobiotic']|modules['detoxification'])-modules['immune_annotated']
pd.DataFrame([{'module':k,'gene':g,'source':NAMES.get(k,'derived_from_msigdb_intersections')} for k,v in modules.items() for g in sorted(v)]).to_csv(ROOT/'results/xenobiotic_module_genes.csv',index=False)
ledger=pd.read_csv(SRC/'results/all_cohort_scores.csv',low_memory=False)
results=[];corr=[];coverage=[]; patient_scores=[]
for cohort in ['TCGA','CGGA693','CGGA325']:
    expr=pd.read_csv(SRC/'data'/f'{cohort}_GBM_log_expression.tsv.gz',sep='\t',index_col=0)
    expr=expr.apply(pd.to_numeric,errors='coerce')
    ranks=expr.rank(axis=0,pct=True)
    clin=ledger[(ledger.cohort==cohort)&ledger.eligible_GBM.eq(True)].set_index('sample')
    clin=clin[~clin.index.duplicated(keep='first')]
    for module,gs in modules.items():
        measured=sorted(set(ranks.index)&gs)
        coverage.append({'cohort':cohort,'module':module,'defined':len(gs),'measured':len(measured)})
        if len(measured)<5:continue
        score=ranks.loc[measured].mean(axis=0).rename('score')
        d=clin.join(score,how='inner')
        d['score_z']=(d.score-d.score.mean())/d.score.std(ddof=0)
        patient_scores.extend([{'cohort':cohort,'sample':idx,'module':module,'score':v} for idx,v in d.score.items()])
        for covset,covs in [('age_MGMT',['age','MGMT'])]+([('age_MGMT_purity',['age','MGMT','purity'])] if cohort=='TCGA' else []):
            f=d[['time','event','score_z']+covs].apply(pd.to_numeric,errors='coerce').dropna()
            if len(f)<30:continue
            fit=CoxPHFitter(penalizer=.01).fit(f,'time','event')
            q=fit.summary.loc['score_z']
            results.append({'cohort':cohort,'module':module,'model':covset,'n':len(f),'events':int(f.event.sum()),
                'coef':float(q['coef']),'HR_per_SD':float(np.exp(q['coef'])),'p':float(q['p'])})
        if cohort=='TCGA':
            for v in ['purity','immune_estimate','myeloid_identity','tumor_glial_proxy']:
                a=d[['score',v]].apply(pd.to_numeric,errors='coerce').dropna()
                rho,p=spearmanr(a.score,a[v]);corr.append({'module':module,'correlate':v,'n':len(a),'rho':rho,'p':p})
r=pd.DataFrame(results)
r['BH_FDR_per_cohort_model']=r.groupby(['cohort','model']).p.transform(lambda x:multipletests(x,method='fdr_bh')[1])
r.to_csv(ROOT/'results/xenobiotic_submodule_cox.csv',index=False)
pd.DataFrame(corr).to_csv(ROOT/'results/xenobiotic_submodule_composition.csv',index=False)
pd.DataFrame(coverage).to_csv(ROOT/'results/xenobiotic_module_coverage.csv',index=False)
pd.DataFrame(patient_scores).to_csv(ROOT/'results/xenobiotic_sample_scores.csv',index=False)
print(pd.DataFrame(coverage).to_string(index=False));print(r.to_string(index=False))
