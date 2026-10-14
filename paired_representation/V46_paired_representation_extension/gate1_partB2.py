"""Gate 1 (part B, corrected) — R3/R4 from the exact inputs used by the frozen
analysis, with R1 reproduced as the correctness gate.

Inputs are the three files listed in the archived input_manifest.csv (SHA-256
verified). Correctness gate: recomputed myeloid aggregate SIGLEC15 must equal the
archived `myeloid_sample_aggregates_scored.csv` values.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

D = Path(r".\biological_validation\data\scrna_pdac")
ANALYSIS = Path(r".\biological_validation"
                r"\analysis_F_cross_scale\outputs")
OUT = Path(r".\data\analysis\V46_validation_extension\Gate1")
OUT.mkdir(parents=True, exist_ok=True)

TUMOUR_TISSUES = ["PDAC", "PDAC_tumor", "PDAC_tissue"]
FILES = ["scRNA_GSE212966.h5ad", "scRNA_GSE156405.h5ad", "scRNA_GSE155698_sub.h5ad"]
MIN_MYELOID = 50


def aggregate(expr: np.ndarray) -> float:
    return float(np.log1p(np.mean(np.expm1(expr.astype(np.float64)))))


def main() -> None:
    rows = []
    per_file = {}
    for name in FILES:
        a = ad.read_h5ad(D / name, backed="r")
        o = a.obs.copy()
        a.file.close()
        keep = (o["tissue"].astype(str).isin(TUMOUR_TISSUES)
                & o["celltype"].astype(str).eq("Myeloid"))
        sub = o.loc[keep].copy()
        sub["_sig"] = pd.to_numeric(sub["SIGLEC15"], errors="coerce").fillna(0.0)
        per_file[name] = {
            "n_cells_total": int(len(o)),
            "n_tumour_myeloid": int(len(sub)),
            "n_samples": int(sub["sample_id"].nunique()),
            "datasets": sorted(sub["dataset"].astype(str).unique()),
        }
        # all-cell tumour compartment for R4 (same tissue filter, every cell type)
        allk = o.loc[o["tissue"].astype(str).isin(TUMOUR_TISSUES)].copy()
        allk["_sig"] = pd.to_numeric(allk["SIGLEC15"], errors="coerce").fillna(0.0)
        for sid, g in sub.groupby("sample_id", observed=True):
            mye = g["_sig"].to_numpy(float)
            if mye.size < MIN_MYELOID:
                continue
            ds = str(g["dataset"].iloc[0])
            pos = mye[mye > 0]
            ac = allk.loc[allk["sample_id"].astype(str).eq(str(sid)), "_sig"].to_numpy(float)
            rows.append({
                "dataset": ds,
                "sample_id": str(sid),
                "n_myeloid": int(mye.size),
                "n_all_cells_tumour": int(ac.size),
                "n_positive": int(pos.size),
                "R1_myeloid_aggregate_recomputed": aggregate(mye),
                "R3_mean_positive_intensity": float(pos.mean()) if pos.size else np.nan,
                "R4_allcell_aggregate": aggregate(ac) if ac.size else np.nan,
            })

    mine = pd.DataFrame(rows)
    arch = pd.read_csv(ANALYSIS / "myeloid_sample_aggregates_scored.csv")[
        ["dataset", "sample_id", "SIGLEC15", "n_myeloid", "n_siglec15_positive"]]

    m = mine.merge(arch, on=["dataset", "sample_id"], how="outer",
                   suffixes=("_new", "_arch"), indicator=True)
    both = m.loc[m["_merge"].eq("both")].copy()
    both["abs_diff_R1"] = (both["R1_myeloid_aggregate_recomputed"] - both["SIGLEC15"]).abs()
    both["d_n"] = both["n_myeloid_new"] - both["n_myeloid_arch"]
    both["d_pos"] = both["n_positive"] - both["n_siglec15_positive"]

    report = {
        "per_file": per_file,
        "n_samples_recomputed": int(len(mine)),
        "match_both": int((m["_merge"] == "both").sum()),
        "only_new": m.loc[m["_merge"] == "left_only", ["dataset", "sample_id"]].to_dict("records"),
        "only_archived": m.loc[m["_merge"] == "right_only", ["dataset", "sample_id"]].to_dict("records"),
        "correctness_gate": {
            "n_compared": int(len(both)),
            "max_abs_diff_R1": float(both["abs_diff_R1"].max()) if len(both) else None,
            "n_R1_exact_1e-9": int((both["abs_diff_R1"] < 1e-9).sum()),
            "n_myeloid_count_match": int((both["d_n"] == 0).sum()),
            "n_positive_count_match": int((both["d_pos"] == 0).sum()),
            "verdict": "EXACT" if len(both) and both["abs_diff_R1"].max() < 1e-9 else "MISMATCH",
        },
        "worst_offenders": both.nlargest(5, "abs_diff_R1")[
            ["dataset", "sample_id", "n_myeloid_new", "n_myeloid_arch",
             "R1_myeloid_aggregate_recomputed", "SIGLEC15", "abs_diff_R1"]
        ].to_dict("records") if len(both) else [],
    }

    mine.to_csv(OUT / "GATE1_R3_R4_by_sample.csv", index=False)
    (OUT / "_gate1_partB_corrected.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False, default=str))


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
