"""Gate 1 (part B) — R3 and R4 from the single-cell objects, with R1 reproduced
as an internal correctness check.

R3 = mean stored SIGLEC15 expression among SIGLEC15-detectable tumour-myeloid
     cells of the sample (positive-cell intensity).
R4 = all-cell pseudo-whole-tumour SIGLEC15 aggregate over every cell type in the
     sample, computed as log1p(mean(exp(x)-1)) to match the archived myeloid
     aggregate convention.

Correctness gate: the myeloid aggregate recomputed here must reproduce the
archived `myeloid_sample_aggregates_scored.csv` SIGLEC15 column.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ANALYSIS = Path(r".\biological_validation"
                r"\analysis_F_cross_scale\outputs")
OUT = Path(r".\data\analysis\V46_validation_extension\Gate1")
OUT.mkdir(parents=True, exist_ok=True)

SOURCES = [
    (r".\data\analysis\scRNA_GSE212966.h5ad", "GSE212966", "tissue", "PDAC"),
    (r".\data\analysis\scRNA_GSE156405_GSE155698_integrated.h5ad", None,
     "tissue_unified", "PDAC"),
]
MIN_MYELOID = 50


def aggregate(expr: np.ndarray) -> float:
    """log1p(mean(exp(x)-1)) over stored log1p values."""
    return float(np.log1p(np.mean(np.expm1(expr.astype(np.float64)))))


def main() -> None:
    rows = []
    for path, forced_ds, tissue_col, tissue_val in SOURCES:
        a = ad.read_h5ad(path, backed="r")
        o = a.obs
        a.file.close()
        keep_tissue = o[tissue_col].astype(str).eq(tissue_val)
        sub = o.loc[keep_tissue]
        if forced_ds:
            sub = sub.assign(_ds=forced_ds)
        else:
            sub = sub.assign(_ds=sub["dataset"].astype(str))
        sub = sub.assign(_sig=pd.to_numeric(sub["SIGLEC15"], errors="coerce").fillna(0.0))
        sub = sub.assign(_mye=sub["celltype"].astype(str).eq("Myeloid"))

        for (ds, sid), g in sub.groupby(["_ds", "sample_id"], observed=True):
            mye = g.loc[g["_mye"], "_sig"].to_numpy(float)
            allc = g["_sig"].to_numpy(float)
            if mye.size < MIN_MYELOID:
                continue
            pos = mye[mye > 0]
            rows.append({
                "dataset": ds,
                "sample_id": str(sid),
                "n_myeloid": int(mye.size),
                "n_all_cells": int(allc.size),
                "n_positive": int(pos.size),
                # R1 reproduced (myeloid aggregate) -- correctness check
                "R1_myeloid_aggregate_recomputed": aggregate(mye),
                # R3 positive-cell intensity
                "R3_mean_positive_intensity": float(pos.mean()) if pos.size else np.nan,
                # R4 all-cell pseudo-whole-tumour aggregate
                "R4_allcell_aggregate": aggregate(allc),
                "n_celltypes": int(g["celltype"].nunique()),
            })
        print(f"  {path.split(chr(92))[-1]}: {len(rows)} cumulative sample rows")

    mine = pd.DataFrame(rows)
    arch = pd.read_csv(ANALYSIS / "myeloid_sample_aggregates_scored.csv")[
        ["dataset", "sample_id", "SIGLEC15", "n_myeloid", "n_siglec15_positive"]]

    chk = mine.merge(arch, on=["dataset", "sample_id"], how="outer",
                     suffixes=("_new", "_arch"), indicator=True)
    diffs = chk.loc[chk["_merge"].eq("both")].copy()
    diffs["abs_diff"] = (diffs["R1_myeloid_aggregate_recomputed"]
                         - diffs["SIGLEC15"]).abs()

    report = {
        "n_samples_recomputed": int(len(mine)),
        "n_samples_archived": int(len(arch)),
        "match_both": int((chk["_merge"] == "both").sum()),
        "only_new": chk.loc[chk["_merge"] == "left_only", ["dataset", "sample_id"]].to_dict("records"),
        "only_archived": chk.loc[chk["_merge"] == "right_only", ["dataset", "sample_id"]].to_dict("records"),
        "R1_reproduction": {
            "max_abs_diff": float(diffs["abs_diff"].max()) if len(diffs) else None,
            "n_within_1e-9": int((diffs["abs_diff"] < 1e-9).sum()),
            "n_compared": int(len(diffs)),
            "verdict": ("EXACT" if len(diffs) and diffs["abs_diff"].max() < 1e-9
                        else "MISMATCH"),
        },
        "n_myeloid_match": int((diffs["n_myeloid_new"] == diffs["n_myeloid_arch"]).sum())
        if "n_myeloid_new" in diffs else None,
    }

    mine.to_csv(OUT / "GATE1_R3_R4_by_sample.csv", index=False)
    (OUT / "_gate1_partB_reproduction.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
