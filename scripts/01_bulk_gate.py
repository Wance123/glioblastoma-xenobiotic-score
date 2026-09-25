import os
"""Prespecified exploratory mannose-utilization gate in adult primary IDH-WT GBM.

Inputs are the official CGGA RSEM and clinical archives. No outcome-driven
feature selection is performed. Run: python scripts/01_bulk_gate.py
"""
from pathlib import Path
import hashlib, json, zipfile
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
LOCAL = Path(os.environ['GEO_LOCAL_ROOT'])
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
GENES = ["MPI", "PMM2", "GMPPA", "GMPPB"]
COHORTS = {
    "CGGA693": (ROOT / "input/CGGA.mRNAseq_693_clinical.20200506.txt.zip", ROOT / "input/CGGA.mRNAseq_693.RSEM-genes.20200506.txt.zip"),
    "CGGA325": (LOCAL / "CGGA.mRNAseq_325_clinical.20200506.txt.zip", LOCAL / "CGGA.mRNAseq_325.RSEM-genes.20200506.txt.zip"),
}

def read_zip(path):
    with zipfile.ZipFile(path) as z:
        name = next(n for n in z.namelist() if n.endswith(".txt") and not n.startswith("__MACOSX"))
        return pd.read_csv(z.open(name), sep="\t", low_memory=False)

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()

summary = []
models = []
genes_out = []
provenance = []
for cohort, (clinical_path, expr_path) in COHORTS.items():
    for p in [clinical_path, expr_path]:
        provenance.append({"cohort": cohort, "file": str(p), "bytes": p.stat().st_size, "sha256": sha256(p)})
    clinical = read_zip(clinical_path)
    expr = read_zip(expr_path)
    clinical["Age"] = pd.to_numeric(clinical["Age"], errors="coerce")
    clinical["OS"] = pd.to_numeric(clinical["OS"], errors="coerce")
    event_col = "Censor (alive=0; dead=1)"
    clinical[event_col] = pd.to_numeric(clinical[event_col], errors="coerce")
    eligible = clinical[(clinical.PRS_type == "Primary") & (clinical.Histology == "GBM") &
        (clinical.IDH_mutation_status == "Wildtype") & (clinical.Age >= 18) &
        (clinical.OS > 0) & clinical[event_col].isin([0,1])].copy()
    assert eligible.CGGA_ID.is_unique
    g = expr.loc[expr.Gene_Name.isin(GENES)].copy()
    assert set(GENES) == set(g.Gene_Name), (cohort, set(GENES)-set(g.Gene_Name))
    assert g.Gene_Name.is_unique
    x = g.set_index("Gene_Name").T.apply(pd.to_numeric, errors="coerce")
    assert (x.to_numpy() >= 0).all()
    eligible = eligible.set_index("CGGA_ID").join(x, how="inner")
    for gene in GENES:
        eligible[gene] = np.log2(eligible[gene] + 1)
    # Within-cohort standardization avoids RNA-seq batch-scale differences.
    for gene in GENES:
        sd = eligible[gene].std(ddof=0)
        eligible[gene+"_z"] = (eligible[gene]-eligible[gene].mean())/sd
    eligible["mannose_core_score"] = eligible[[g+"_z" for g in GENES]].mean(axis=1)
    eligible["event"] = eligible[event_col].astype(int)
    eligible["time"] = eligible.OS
    eligible["MGMT_methylated"] = eligible.MGMTp_methylation_status.map({"methylated":1,"un-methylated":0})
    eligible["TMZ"] = pd.to_numeric(eligible["Chemo_status (TMZ treated=1;un-treated=0)"],errors="coerce")
    eligible["RT"] = pd.to_numeric(eligible["Radio_status (treated=1;un-treated=0)"],errors="coerce")
    eligible.to_csv(OUT/f"{cohort}_eligible_scores.csv")
    summary.append({"cohort":cohort,"all_samples":len(clinical),"eligible":len(eligible),
        "events":int(eligible.event.sum()),"median_OS_days":float(eligible.time.median()),
        "MGMT_known":int(eligible.MGMT_methylated.notna().sum()),
        "TMZ_known":int(eligible.TMZ.notna().sum()),"RT_known":int(eligible.RT.notna().sum())})
    for gene in GENES:
        c = CoxPHFitter().fit(eligible[["time","event",gene+"_z"]],"time","event")
        row = c.summary.loc[gene+"_z"]
        genes_out.append({"cohort":cohort,"gene":gene,"HR_per_SD":float(np.exp(row["coef"])),"p":float(row["p"])})
    for specification, covars in [
        ("unadjusted",["mannose_core_score"]),
        ("age_adjusted",["mannose_core_score","Age"]),
        ("clinical_complete_case",["mannose_core_score","Age","MGMT_methylated","TMZ","RT"]),
    ]:
        d = eligible[["time","event"]+covars].dropna()
        if len(d) < 30 or d.event.sum() < 15:
            models.append({"cohort":cohort,"model":specification,"n":len(d),"events":int(d.event.sum()),"status":"insufficient"})
            continue
        c = CoxPHFitter(penalizer=0.01).fit(d,"time","event")
        row = c.summary.loc["mannose_core_score"]
        models.append({"cohort":cohort,"model":specification,"n":len(d),"events":int(d.event.sum()),
            "HR_per_score_unit":float(np.exp(row["coef"])),"CI_low":float(np.exp(row["coef lower 95%"])),
            "CI_high":float(np.exp(row["coef upper 95%"])),"p":float(row["p"]),
            "concordance":float(c.concordance_index_),"status":"fit"})

gene_df = pd.DataFrame(genes_out)
gene_df["BH_FDR_within_cohort"] = gene_df.groupby("cohort")["p"].transform(lambda v: multipletests(v,method="fdr_bh")[1])
pd.DataFrame(summary).to_csv(OUT/"cohort_audit.csv",index=False)
pd.DataFrame(models).to_csv(OUT/"survival_models.csv",index=False)
gene_df.to_csv(OUT/"gene_exploratory_models.csv",index=False)
pd.DataFrame(provenance).to_csv(OUT/"input_sha256.csv",index=False)
print(pd.DataFrame(summary).to_string(index=False))
print(pd.DataFrame(models).to_string(index=False))
