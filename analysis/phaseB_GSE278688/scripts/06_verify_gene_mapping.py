#!/usr/bin/env python
# =============================================================================
# 06_verify_gene_mapping.py -- Phase B step 1, GSE278688 (FROZEN RUN)
# Independent control that the MMREAD row order equals the features-file order.
# A wrong gene->row map would invalidate everything, so this is checked, not
# assumed: known housekeeping/lineage genes must show the expected pattern.
# =============================================================================
import os
import sys, os, gzip, time
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np, pandas as pd, scipy.io

ROOT = os.environ.get("SIGLEC15_GSE278688_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927\phaseB_GSE278688")
DATA = os.path.join(ROOT, "data")
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

feat = pd.read_csv(os.path.join(DATA,"GSE278688_sc_features.tsv.gz"), sep="\t", header=None)
sym  = feat[1].values
meta = pd.read_csv(os.path.join(DATA,"GSE278688_sc_metadata.csv.gz"), index_col=0)

M = scipy.io.mmread(os.path.join(DATA,"GSE278688_sc_matrix.mtx.gz")).tocsr()
log(f"matrix {M.shape} nnz={M.nnz}")

ct = meta["all_celltype"].values

# Positive controls: each gene should peak in its expected compartment.
controls = {
  "EPCAM":   ("Ductal",       "epithelial"),
  "KRT19":   ("Ductal",       "epithelial"),
  "PTPRC":   ("CD8T",         "pan-immune (CD45)"),
  "CD3D":    ("CD8T",         "T cell"),
  "CD79A":   ("B",            "B cell"),
  "COL1A1":  ("CAF",          "fibroblast"),
  "LYZ":     ("Myeloid",      "myeloid"),
  "CD68":    ("Myeloid",      "macrophage"),
  "CPA3":    ("Myeloid",      "mast (myeloid-labelled)"),
  "SIGLEC15":("Myeloid",      "myeloid checkpoint"),
  "ACP5":    ("Myeloid",      "osteoclast/ macrophage"),
  "MMP9":    ("Myeloid",      "neutrophil/ macrophage"),
}
log(f"\n{'gene':10s} {'row':>6s} {'pct_expr':>9s} {'mean_all':>10s} | mean by compartment (highest 3)")
for g,(exp,desc) in controls.items():
    w = np.where(sym == g)[0]
    if len(w)==0: log(f"{g:10s}  ABSENT"); continue
    r = int(w[0])
    v = np.asarray(M[r,:].todense()).ravel()
    pct = 100*np.mean(v>0); mean_all = v.mean()
    tab = pd.Series(v).groupby(ct).mean().sort_values(ascending=False)
    top = ", ".join(f"{k}={x:.2f}" for k,x in tab.head(3).items())
    flag = "" if tab.index[0]==exp else f"   <-- EXPECTED {exp}"
    log(f"{g:10s} {r:6d} {pct:8.2f}% {mean_all:10.4f} | {top}{flag}   [{desc}]")

log("\nCONTROL INTERPRETATION: if gene->row mapping were scrambled, the positive")
log("controls above would NOT peak in their expected compartment.")
