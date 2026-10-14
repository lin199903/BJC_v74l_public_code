#!/usr/bin/env python
# =============================================================================
# 05_normalisation_fork_test.py -- Phase B step 1, GSE278688 (FROZEN RUN)
#
# The frozen brief says "SIGLEC15 detectability = stored log1p-normalised
# expression > 0". The DEPOSITED matrix is raw integer UMI counts (verified:
# MatrixMarket flag 'integer', all values integral, min 1, max 16277).
# So "stored log1p-normalised" does not exist as such.
#
# This script establishes, BEFORE any analysis choice is made, exactly which
# frozen-rule quantities are invariant to the normalisation choice and which
# are not. It makes NO analysis choice; it measures the size of the fork.
#
#   Candidate A: log1p(count)                     (no library-size scaling)
#   Candidate B: log1p(1e4 * count / libsize)     (standard log-normalisation)
# =============================================================================
import os
import sys, os, gzip, time, json
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np, pandas as pd, scipy.io, scipy.sparse as sp

ROOT = os.environ.get("SIGLEC15_GSE278688_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927\phaseB_GSE278688")
OUT  = os.path.join(ROOT, "output")
t0=time.time()
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

d = pd.read_csv(os.path.join(OUT, "_myeloid_gene_cache.csv"))
log(f"myeloid cache: {d.shape}")

GENES = ["SIGLEC15","ACP5","MMP9"]
lib = d["libsize"].values.astype(float)
raw = {g: d[g].values.astype(float) for g in GENES}

# ---------- Candidate A ----------
A = {g: np.log1p(raw[g]) for g in GENES}
# ---------- Candidate B ----------
scale = 1e4 / lib
B = {g: np.log1p(raw[g] * scale) for g in GENES}

out = {}
log("\n=== 1. Is the DETECTABILITY set invariant? ===")
for g in GENES:
    det_raw = raw[g] > 0
    for nm, M in (("A log1p(count)", A), ("B log1p(CP10K)", B)):
        det = M[g] > 0
        log(f"  {g:9s} {nm:16s} n_detect={det.sum():5d}  identical_to_raw>0: {bool(np.array_equal(det,det_raw))}")
    out[f"detect_invariant_{g}"] = True

log("\n=== 2. How different are A and B for the PROGRAMME genes? ===")
for g in ["ACP5","MMP9"]:
    a, b = A[g], B[g]
    nz = a > 0
    # correlation among cells that express the gene at all
    if nz.sum() > 2:
        r_p = np.corrcoef(a[nz], b[nz])[0,1]
        r_s = pd.Series(a[nz]).corr(pd.Series(b[nz]), method="spearman")
        log(f"  {g:9s} n_expr={nz.sum():5d}  Pearson(A,B)={r_p:.4f}  Spearman(A,B)={r_s:.4f}")
        log(f"            mean A={a[nz].mean():.4f}  mean B={b[nz].mean():.4f}")
        log(f"            max |A-B| = {np.abs(a[nz]-b[nz]).max():.4f}")

log("\n=== 3. Within-sample z-scores of the programme genes, A vs B ===")
sid = d["sample_id"].values
for g in ["ACP5","MMP9"]:
    zg_A = np.full(len(d), np.nan); zg_B = np.full(len(d), np.nan)
    for s in pd.unique(sid):
        i = np.where(sid == s)[0]
        for arr, zg in ((A[g], zg_A), (B[g], zg_B)):
            x = arr[i]; sd = x.std(ddof=1)
            zg[i] = (x - x.mean())/sd if sd > 0 else np.nan
    m = np.isfinite(zg_A) & np.isfinite(zg_B)
    r_p = np.corrcoef(zg_A[m], zg_B[m])[0,1]
    r_s = pd.Series(zg_A[m]).corr(pd.Series(zg_B[m]), method="spearman")
    log(f"  {g:9s} z(A) vs z(B): Pearson={r_p:.6f}  Spearman={r_s:.6f}  max|diff|={np.abs(zg_A[m]-zg_B[m]).max():.4f}")

log("\n=== 4. The sample contrast, A vs B (the quantity that actually matters) ===")
rows = []
for s in sorted(pd.unique(sid)):
    i = np.where(sid == s)[0]
    det = raw["SIGLEC15"][i] > 0
    rec = {"sample_id": s, "n_det": int(det.sum()), "n_nondet": int((~det).sum())}
    for nm, M in (("A", A), ("B", B)):
        zs = []
        inv = []
        for g in ["ACP5","MMP9"]:
            x = M[g][i]; sd = x.std(ddof=1)
            if sd == 0: inv.append(g); zs.append(None)
            else: zs.append((x - x.mean())/sd)
        if any(z is None for z in zs):
            rec[f"prog_{nm}"] = np.nan; rec[f"inv_{nm}"] = ",".join(inv)
        else:
            prog = (zs[0] + zs[1])/2
            rec[f"prog_{nm}"] = prog
            rec[f"inv_{nm}"] = ""
        rec[f"contrast_{nm}"] = (np.nan if rec[f"prog_{nm}"] is np.nan
                                 else float(rec[f"prog_{nm}"][det].mean() - rec[f"prog_{nm}"][~det].mean()))
    rows.append(rec)

# invariances of the programme gene within sample
invtab = []
for r in rows:
    invtab.append((r["sample_id"], r["inv_A"], r["inv_B"]))
log("  sample   invariant_genes_A   invariant_genes_B")
for s,a,b in invtab: log(f"    {s:16s} '{a}'   '{b}'")

cA = np.array([r["contrast_A"] for r in rows], float)
cB = np.array([r["contrast_B"] for r in rows], float)

# restrict to frozen-eligible samples (Tumor + >=50 myeloid + both groups present)
m = pd.read_csv(os.path.join(OUT,"sample_manifest.csv")) if os.path.exists(os.path.join(OUT,"sample_manifest.csv")) else None
elig = [r["sample_id"] for r in rows if r["n_det"]>=1 and r["n_nondet"]>=1]
log(f"\n  samples with both groups (any tissue): {len(elig)}")

log("\n  per-sample contrast A vs B:")
log("    sample            n_det n_nondet   contrast_A   contrast_B   A-B")
for r in rows:
    log(f"    {r['sample_id']:16s} {r['n_det']:5d} {r['n_nondet']:8d}   {r['contrast_A']:10.4f}   {r['contrast_B']:10.4f}  {r['contrast_A']-r['contrast_B']:7.4f}")

ok = np.isfinite(cA) & np.isfinite(cB)
log(f"\n  mean contrast A = {np.nanmean(cA[ok]):.4f}")
log(f"  mean contrast B = {np.nanmean(cB[ok]):.4f}")
log(f"  Pearson(A,B)    = {np.corrcoef(cA[ok],cB[ok])[0,1]:.6f}")
log(f"  Spearman(A,B)   = {pd.Series(cA[ok]).corr(pd.Series(cB[ok]),method='spearman'):.6f}")
log(f"  max |A-B|       = {np.abs(cA[ok]-cB[ok]).max():.4f}")
log(f"  sign agreement  = {int((np.sign(cA[ok])==np.sign(cB[ok])).sum())}/{ok.sum()}")

out.update(dict(n_myeloid=int(len(d)), mean_contrast_A=float(np.nanmean(cA[ok])),
                mean_contrast_B=float(np.nanmean(cB[ok])),
                pearson_contrast=float(np.corrcoef(cA[ok],cB[ok])[0,1]),
                max_abs_diff_contrast=float(np.abs(cA[ok]-cB[ok]).max()),
                sign_agreement=f"{int((np.sign(cA[ok])==np.sign(cB[ok])).sum())}/{int(ok.sum())}"))
with open(os.path.join(OUT,"_normalisation_fork_test.json"),"w") as f:
    json.dump(out, f, indent=2)
log(f"\nwritten _normalisation_fork_test.json  TOTAL {time.time()-t0:.1f}s")
