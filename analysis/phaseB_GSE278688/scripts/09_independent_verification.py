#!/usr/bin/env python
# =============================================================================
# 09_independent_verification.py -- Phase B step 1, GSE278688 (FROZEN RUN)
#
# Independent re-implementation of the primary contrast IN A DIFFERENT LANGUAGE,
# reusing the manuscript's OWN functions by import, to check the R script
# (08_final_analysis.R). Three things are load-bearing and are each tested:
#   (1) ddof=0 vs ddof=1 in the within-sample z-standardisation
#   (2) CP10K+log1p vs log1p(raw counts)
#   (3) that no sample's programme score silently used a single gene
# A positive control is included so that "agrees" is not a vacuous statement.
# =============================================================================
import os
import sys, os, json
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.join(os.environ.get("SIGLEC15_SUBMISSION_ROOT", r"E:\Tex\HMlong\npjPrecOncol\FINAL_SUBMISSION_v6"), "biological_validation", "analysis_F_cross_scale"))
import numpy as np, pandas as pd
import analysis_F_cross_scale as F          # the manuscript's own module

ROOT = os.environ.get("SIGLEC15_GSE278688_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927\phaseB_GSE278688")
OUT  = os.path.join(ROOT, "output")
def log(m): print(m, flush=True)

log(f"imported manuscript module; zscore_columns = {F.zscore_columns}")
log(f"manuscript SEED={F.SEED}  N_RESAMPLES={F.N_RESAMPLES}")

d = pd.read_csv(os.path.join(OUT, "_myeloid_gene_cache.csv"))
R = pd.read_csv(os.path.join(OUT, "primary_contrast_by_sample.csv"))
elig = list(R.loc[R.status == "OK", "sample_id"])
log(f"eligible samples from R output: {elig}")

def build(df, mode):
    out = df.copy()
    if mode == "cp10k":
        out[["SIGLEC15","ACP5","MMP9"]] = np.log1p(
            out[["SIGLEC15","ACP5","MMP9"]].to_numpy(float)
            * 1e4 / out["libsize"].to_numpy(float)[:,None])
    elif mode == "log1p_raw":
        out[["SIGLEC15","ACP5","MMP9"]] = np.log1p(out[["SIGLEC15","ACP5","MMP9"]].to_numpy(float))
    else:
        raise ValueError(mode)
    return out

def contrast(df, mode="cp10k", ddof_mode="manuscript"):
    """ddof_mode='manuscript' -> F.zscore_columns (ddof=0); 'sample' -> ddof=1."""
    dd = build(df, mode)
    rows = []
    for s in elig:
        g = dd[dd.sample_id == s]
        vals = g[["ACP5","MMP9"]].to_numpy(float)
        if ddof_mode == "manuscript":
            z = F.zscore_columns(vals)                 # ddof=0, sd==0 -> nan
        else:
            z = (vals - vals.mean(0)) / vals.std(0, ddof=1)
        score = np.nanmean(z, axis=1)
        pos = g["SIGLEC15"].to_numpy(float) > 0
        rows.append({"sample_id": s, "n_det": int(pos.sum()), "n_non": int((~pos).sum()),
                     "contrast": float(np.nanmean(score[pos]) - np.nanmean(score[~pos]))})
    return pd.DataFrame(rows)

log("\n=== (1)+(2) primary contrast under each convention ===")
a = contrast(d, "cp10k",     "manuscript")   # what R does
b = contrast(d, "cp10k",     "sample")       # ddof=1
c = contrast(d, "log1p_raw", "manuscript")   # no library-size scaling

cmp = a[["sample_id","contrast"]].rename(columns={"contrast":"python_cp10k_ddof0"})
cmp = cmp.merge(b[["sample_id","contrast"]].rename(columns={"contrast":"python_cp10k_ddof1"}), on="sample_id")
cmp = cmp.merge(c[["sample_id","contrast"]].rename(columns={"contrast":"python_log1praw_ddof0"}), on="sample_id")
cmp = cmp.merge(R[["sample_id","contrast"]].rename(columns={"contrast":"R_cp10k_ddof0"}), on="sample_id")
cmp["diff_R_vs_python"] = cmp["R_cp10k_ddof0"] - cmp["python_cp10k_ddof0"]
pd.set_option("display.width", 200)
log(cmp.to_string(index=False))

log("")
log(f"R vs Python (both CP10K, ddof=0) max abs diff : {cmp['diff_R_vs_python'].abs().max():.3e}")
log(f"  -> AGREEMENT: {np.allclose(cmp['R_cp10k_ddof0'], cmp['python_cp10k_ddof0'], atol=1e-9)}")
log(f"mean contrast  R={cmp['R_cp10k_ddof0'].mean():.6f}  Python={cmp['python_cp10k_ddof0'].mean():.6f}")
log(f"mean contrast  ddof=1 (CP10K)      = {cmp['python_cp10k_ddof1'].mean():.6f}")
log(f"mean contrast  log1p(raw), ddof=0  = {cmp['python_log1praw_ddof0'].mean():.6f}")
log("")
log("POSITIVE CONTROL for the test itself: the two alternative conventions must")
log("produce values that DIFFER from the adopted one, otherwise the comparison")
log(f"  ddof=1   vs adopted: max abs diff = {(cmp['python_cp10k_ddof1']-cmp['python_cp10k_ddof0']).abs().max():.4f}")
log(f"  log1praw vs adopted: max abs diff = {(cmp['python_log1praw_ddof0']-cmp['python_cp10k_ddof0']).abs().max():.4f}")

log("\n=== (3) invariance of the programme genes within each sample ===")
dd = build(d, "cp10k")
inv_rows = []
for s in sorted(d.sample_id.unique()):
    g = dd[dd.sample_id == s]
    row = {"sample_id": s, "n_myeloid": len(g), "tissue": g.tissue.iloc[0]}
    for gene in ("ACP5","MMP9"):
        row[f"{gene}_sd_ddof0"] = float(g[gene].to_numpy(float).std(0))
    row["invariant"] = ",".join([g_ for g_ in ("ACP5","MMP9") if row[f"{g_}_sd_ddof0"] == 0]) or "-"
    inv_rows.append(row)
inv = pd.DataFrame(inv_rows)
log(inv.to_string(index=False))
log("")
sd0 = inv[inv.invariant != "-"]
log(f"samples where a programme gene is invariant within sample: {len(sd0)} of {len(inv)}")
log(f"  of which TUMOUR samples: {int((sd0.tissue=='Tumor').sum())} of {int((inv.tissue=='Tumor').sum())}")
log(f"  within the ELIGIBLE set: {len([s for s in sd0.sample_id if s in elig])} of {len(elig)}")

out = {"R_vs_python_agreement": bool(np.allclose(cmp['R_cp10k_ddof0'], cmp['python_cp10k_ddof0'], atol=1e-9)),
       "max_abs_diff_R_python": float(cmp['diff_R_vs_python'].abs().max()),
       "mean_contrast_R": float(cmp['R_cp10k_ddof0'].mean()),
       "mean_contrast_ddof1": float(cmp['python_cp10k_ddof1'].mean()),
       "mean_contrast_log1p_raw": float(cmp['python_log1praw_ddof0'].mean()),
       "samples_with_invariant_gene_all": int(len(sd0)),
       "samples_with_invariant_gene_eligible": int(len([s for s in sd0.sample_id if s in elig])),
       "invariant_detail": sd0[["sample_id","tissue","invariant"]].to_dict("records")}
with open(os.path.join(OUT,"_verification.json"),"w") as f: json.dump(out, f, indent=2)
log("\nwritten _verification.json")
