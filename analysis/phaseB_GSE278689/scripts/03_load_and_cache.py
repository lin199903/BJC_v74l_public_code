#!/usr/bin/env python
# =============================================================================
# 03_load_and_cache.py -- Phase B step 2, GSE278689
#
# Streams the deposited Matrix Market file ONCE and writes two caches:
#
#   output/_myeloid_gene_cache.csv   every deposited `Myeloid` cell, with raw
#                                    counts for SIGLEC15, ACP5, MMP9 and its
#                                    full library size (sum over ALL genes)
#   output/_allcells_siglec15_cache.csv  every cell, with raw SIGLEC15 count
#                                    and full library size
#
# No filtering, no threshold, no annotation is applied here. Cells are carried
# through with their deposited labels so that eligibility can be decided later,
# explicitly, by the frozen rules.
#
# The matrix is not subset to "Tumor" at load time: the manifest must be able to
# show what was excluded and why.
# =============================================================================
from __future__ import annotations

import os
import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("SIGLEC15_GSE278689_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927\phaseB_GSE278689"))
DATA = ROOT / "data"
OUT = ROOT / "output"
TARGETS = ["SIGLEC15", "ACP5", "MMP9"]
# Marker genes are captured for EVERY cell so that step 06 can verify the
# gene -> row mapping with positive controls. They are used for verification
# only; this run never re-annotates a cell.
MARKERS = ["KRT19", "COL1A1", "CD68", "CPA3", "EPCAM", "PTPRC"]
CAPTURE = TARGETS + MARKERS
CHUNK = 5_000_000


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def read_mtx_header(path: Path):
    """Return (n_rows_genes, n_cols_cells, n_entries, n_leading_lines).

    n_leading_lines counts every physical line before the first data line, so the
    caller can skip exactly that many and never re-parse the dimension line as an
    observation (a silent corruption that would still look dimensionally valid).
    """
    leading = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("%"):
                leading += 1
                continue
            parts = line.split()
            if len(parts) != 3:
                leading += 1
                continue
            return int(parts[0]), int(parts[1]), int(parts[2]), leading + 1
    raise RuntimeError("no dimensions found in matrix header")


def main() -> None:
    log("ROOT", ROOT)
    log("numpy", np.__version__, "pandas", pd.__version__)

    t0 = time.time()
    feats = pd.read_csv(DATA / "GSE278689_sn_features.tsv.gz", sep="\t", header=None)
    genes = feats[1].astype(str).tolist()
    gene_row = {g: i + 1 for i, g in enumerate(genes)}          # 1-based mtx row
    log("features", len(genes), "genes;  columns:", list(feats.columns))

    missing = [g for g in TARGETS if g not in gene_row]
    if missing:
        raise SystemExit(f"FATAL: target gene(s) absent from feature table: {missing}")
    absent_markers = [g for g in MARKERS if g not in gene_row]
    if absent_markers:
        log("marker genes absent from feature table:", absent_markers)
    tgt_row = {g: gene_row[g] for g in CAPTURE if g in gene_row}
    log("target gene rows (1-based):", tgt_row)
    dup = sorted({g for g in genes if genes.count(g) > 1}) if len(set(genes)) != len(genes) else []
    log("duplicated gene symbols:", dup if dup else "none")

    bcs = pd.read_csv(DATA / "GSE278689_sn_barcodes.tsv.gz", sep="\t", header=None)[0].astype(str).tolist()
    md = pd.read_csv(DATA / "GSE278689_sn_metadata.csv.gz", index_col=0)
    md.index = md.index.astype(str)
    md.index.name = "barcode"
    log("barcodes", len(bcs), " unique", len(set(bcs)), " metadata rows", len(md))
    if len(bcs) != len(md) or set(bcs) != set(md.index):
        raise SystemExit("FATAL: barcode file and metadata do not describe the same cells")

    n_genes, n_cells, n_entries, n_skip = read_mtx_header(DATA / "GSE278689_sn_matrix.mtx.gz")
    log(f"matrix header: {n_genes} genes x {n_cells} cells, {n_entries:,} entries, "
        f"{n_skip} leading line(s) to skip")
    if n_genes != len(genes) or n_cells != len(bcs):
        raise SystemExit("FATAL: matrix header does not match features/barcodes")

    ## ---------------------------------------------------------------- pass ---
    libsize = np.zeros(n_cells, dtype=np.int64)
    ngenes = np.zeros(n_cells, dtype=np.int32)
    val = {g: np.zeros(n_cells, dtype=np.float64) for g in CAPTURE}
    seen = 0
    non_integer = 0
    neg = 0
    dim_bad = 0
    max_value = 0.0

    reader = pd.read_csv(DATA / "GSE278689_sn_matrix.mtx.gz", sep=r"\s+", header=None,
                         comment="%", skiprows=n_skip, dtype={0: np.int32, 1: np.int32, 2: np.float64},
                         chunksize=CHUNK, engine="c")
    for k, ch in enumerate(reader):
        ch.columns = ["g", "c", "v"]
        seen += len(ch)
        if seen > n_entries:
            raise SystemExit(f"FATAL: read more than the declared {n_entries} entries -- "
                             "the dimension line was parsed as data")
        g = ch["g"].to_numpy()
        c = ch["c"].to_numpy() - 1
        v = ch["v"].to_numpy()
        if (c < 0).any() or (c >= n_cells).any() or (g < 1).any() or (g > n_genes).any():
            bad = (c < 0) | (c >= n_cells) | (g < 1) | (g > n_genes)
            dim_bad += int(bad.sum())
            g, c, v = g[~bad], c[~bad], v[~bad]
        non_integer += int((v != np.floor(v)).sum())
        neg += int((v < 0).sum())
        if v.size:
            max_value = max(max_value, float(v.max()))
        libsize += np.bincount(c, weights=v, minlength=n_cells).astype(np.int64)
        ngenes += np.bincount(c, minlength=n_cells).astype(np.int32)
        for name in CAPTURE:
            r = tgt_row.get(name)
            if r is None:
                continue
            m = g == r
            if m.any():
                val[name][c[m]] = v[m]
        log(f"  chunk {k}: {len(ch):,} entries  cumulative {seen:,}/{n_entries:,} "
            f"({100*seen/n_entries:.1f}%)")
        if seen >= n_entries:
            break

    log(f"read {seen:,} entries in {time.time()-t0:.1f}s   max stored value {max_value:g}")
    if seen != n_entries:
        raise SystemExit(f"FATAL: read {seen} entries, header declared {n_entries}")
    if non_integer:
        log(f"WARNING: {non_integer} non-integer values (deposit is not pure integer counts)")
    if neg:
        log(f"WARNING: {neg} negative values")
    if dim_bad:
        log(f"WARNING: {dim_bad} out-of-range indices skipped")

    ## --------------------------------------------------------------- caches --
    bar = np.array(bcs, dtype=object)
    sample_id = np.array([b.rsplit("_", 1)[0] for b in bcs], dtype=object)

    allc = pd.DataFrame({
        "barcode": bar,
        "sample_id": sample_id,
        "patient": md["patients"].to_numpy(),
        "tissue": md["tissue"].to_numpy(),
        "all_celltype": md["all_celltype"].to_numpy(),
        "SIGLEC15": val["SIGLEC15"],
        "libsize": libsize,
        "n_genes": ngenes,
    })
    for g in MARKERS:
        if g in val:
            allc[g] = val[g]
    allc.to_csv(OUT / "_allcells_siglec15_cache.csv", index=False)
    log("wrote _allcells_siglec15_cache.csv", allc.shape)

    md = md.loc[bcs]
    my = md["all_celltype"].to_numpy() == "Myeloid"
    mye = pd.DataFrame({
        "barcode": bar[my],
        "sample_id": sample_id[my],
        "patient": md["patients"].to_numpy()[my],
        "tissue": md["tissue"].to_numpy()[my],
        "SIGLEC15": val["SIGLEC15"][my],
        "ACP5": val["ACP5"][my],
        "MMP9": val["MMP9"][my],
        "libsize": libsize[my],
        "n_genes": ngenes[my],
    })
    mye.to_csv(OUT / "_myeloid_gene_cache.csv", index=False)
    log("wrote _myeloid_gene_cache.csv", mye.shape)

    summary = {
        "matrix": "GSE278689_sn_matrix.mtx.gz",
        "n_genes": int(n_genes), "n_cells": int(n_cells), "n_entries": int(n_entries),
        "entries_read": int(seen),
        "non_integer_values": int(non_integer),
        "negative_values": int(neg),
        "out_of_range_indices": int(dim_bad),
        "all_values_integral": bool(non_integer == 0),
        "max_stored_value": max_value,
        "n_leading_header_lines": int(n_skip),
        "libsize_min": int(libsize.min()), "libsize_max": int(libsize.max()),
        "libsize_median": float(np.median(libsize)),
        "libsize_zero_cells": int((libsize == 0).sum()),
        "myeloid_cells": int(my.sum()),
        "siglec15_total_positive_cells": int((val["SIGLEC15"] > 0).sum()),
        "runtime_seconds": round(time.time() - t0, 1),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    (OUT / "_matrix_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    log(json.dumps(summary, indent=1))
    log("DONE")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
