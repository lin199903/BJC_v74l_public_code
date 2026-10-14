#!/usr/bin/env python
# =============================================================================
# 03_load_and_cache.py -- Phase B step 1, GSE278688 (FROZEN RUN)
# Loads the deposited sparse matrix ONCE, keeps only the 9,102 Myeloid cells
# and the genes needed, caches to RDS-friendly CSV/feather for the R analysis.
# READ-ONLY w.r.t. the downloaded deposit.
#
# CRITICAL: does NOT modify cell-type labels. Uses only `all_celltype`.
# =============================================================================
import os
import sys, os, gzip, time, json
sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")

import numpy as np, pandas as pd, scipy.io, scipy.sparse as sp

ROOT = os.environ.get("SIGLEC15_GSE278688_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927\phaseB_GSE278688")
DATA = os.path.join(ROOT, "data")
OUT  = os.path.join(ROOT, "output")
os.makedirs(OUT, exist_ok=True)
t0 = time.time()
def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

mtx_p = os.path.join(DATA, "GSE278688_sc_matrix.mtx.gz")
log(f"matrix : {mtx_p} ({os.path.getsize(mtx_p):,} bytes)")

# ---- header ----
with gzip.open(mtx_p, "rt") as f:
    for line in f:
        if line.startswith("%"):
            log(f"comment: {line.strip()}"); continue
        dims = [int(x) for x in line.split()]
        break
nrow, ncol, nnz = dims
log(f"dims   : nrow={nrow} ncol={ncol} nnz={nnz}")

feat = pd.read_csv(os.path.join(DATA, "GSE278688_sc_features.tsv.gz"),
                   sep="\t", header=None)
bar  = pd.read_csv(os.path.join(DATA, "GSE278688_sc_barcodes.tsv.gz"),
                   sep="\t", header=None)[0].values
meta = pd.read_csv(os.path.join(DATA, "GSE278688_sc_metadata.csv.gz"), index_col=0)

log(f"features {feat.shape}  barcodes {len(bar)}  metadata {meta.shape}")
assert nrow == len(feat), f"nrow {nrow} != features {len(feat)}"
assert ncol == len(bar),  f"ncol {ncol} != barcodes {len(bar)}"
assert list(meta.index) == list(bar), "metadata index != barcodes order"

# ---- read matrix ----
log("reading MTX (this is the expensive step) ...")
M = scipy.io.mmread(mtx_p)
log(f"mmread done: {M.shape} nnz={M.nnz}  ({time.time()-t0:.1f}s)")
M = M.tocsr()   # rows x cols
log(f"csr done    ({time.time()-t0:.1f}s)")

if M.shape == (nrow, ncol):
    orient = "genes x cells (rows=features)"
elif M.shape == (ncol, nrow):
    orient = "cells x genes (rows=barcodes)"
    M = M.T.tocsr()
else:
    raise SystemExit(f"unexpected shape {M.shape} vs header ({nrow},{ncol})")
log(f"orientation: {orient}")

# ---- value character ----
vals = M.data
log(f"value min/max      : {vals.min():.6g} / {vals.max():.6g}")
log(f"integer-valued     : {bool(np.all(vals == np.round(vals)))}")
log(f"frac == 1          : {float(np.mean(vals == 1)):.4f}")
log(f"frac non-integer   : {float(np.mean(vals != np.round(vals))):.6f}")

# ---- per-cell library size (column sums) ----
lib = np.asarray(M.sum(axis=0)).ravel()
log(f"libsize min/median/max: {lib.min():.4g} / {np.median(lib):.4g} / {lib.max():.4g}")

# ---- myeloid subset ----
my_mask = (meta["all_celltype"].values == "Myeloid")
my_idx = np.where(my_mask)[0]
log(f"Myeloid cells: {len(my_idx)}")

# ---- gene rows of interest (+ housekeeping for orientation sanity) ----
sym = feat[1].values
probe_genes = ["SIGLEC15", "ACP5", "MMP9", "EPCAM", "COL1A1", "PTPRC", "CD3D", "LYZ"]
gidx = {}
for g in probe_genes:
    w = np.where(sym == g)[0]
    if len(w) == 0:
        log(f"  !! {g} ABSENT"); continue
    gidx[g] = int(w[0])
log(f"gene rows: {gidx}")

# ---- cache: myeloid cells x genes of interest, RAW stored values ----
sub = M[np.array(list(gidx.values())), :][:, my_idx]           # genes x myeloid
sub = sub.T.tocsr()                                            # myeloid x genes
df = pd.DataFrame(sub.toarray(), columns=list(gidx.keys()), index=None)
df.insert(0, "cell", bar[my_idx])
df.insert(1, "sample_id", meta["sample_id"].values[my_idx] if "sample_id" in meta
          else [b.rsplit("-GEX_out_", 1)[0].rsplit("_GEX_out_", 1)[0] for b in bar[my_idx]])
df.insert(2, "tissue", meta["tissue"].values[my_idx])
df.insert(3, "patient", meta["patients"].values[my_idx])
df.insert(4, "libsize", lib[my_idx])
df.insert(5, "n_genes", np.asarray((M[:, my_idx] > 0).sum(axis=0)).ravel())
df.to_csv(os.path.join(OUT, "_myeloid_gene_cache.csv"), index=False)
log(f"cache written: {df.shape} -> _myeloid_gene_cache.csv")

# ---- cache 2: SIGLEC15 for ALL cells (needed for specimen-level all-cell agg)
si = gidx["SIGLEC15"]
all_sig = np.asarray(M[si, :].todense()).ravel()
pid = np.array([b.replace("-GEX_out_", "|").replace("_GEX_out_", "|").split("|")[0]
                for b in bar])
adf = pd.DataFrame({
    "cell": bar, "sample_id": pid,
    "tissue": meta["tissue"].values,
    "patient": meta["patients"].values,
    "all_celltype": meta["all_celltype"].values,
    "SIGLEC15": all_sig, "libsize": lib})
adf.to_csv(os.path.join(OUT, "_allcells_siglec15_cache.csv"), index=False)
log(f"all-cell cache written: {adf.shape} -> _allcells_siglec15_cache.csv")

summary = dict(
    nrow=nrow, ncol=ncol, nnz=int(M.nnz), orientation=orient,
    values_integer=bool(np.all(vals == np.round(vals))),
    value_min=float(vals.min()), value_max=float(vals.max()),
    n_myeloid=int(len(my_idx)), genes=list(gidx.keys()),
    libsize_median=float(np.median(lib)), libsize_min=float(lib.min()),
    libsize_max=float(lib.max()),
)
with open(os.path.join(OUT, "_matrix_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
log("summary written -> _matrix_summary.json")
log(f"TOTAL {time.time()-t0:.1f}s")
