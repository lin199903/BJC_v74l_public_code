"""Cross-check: Task5 annotation table (from .raw), SIGLEC15 positive rate vs the positive rate re-read by V04 from .X.

Used only for the completeness statement of Gate 0 provenance, not involved in any biological inference.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(r".")
OUT = ROOT / "data" / "analysis" / "V04_discriminating_analysis"
ANNOTATION = ROOT / "data" / "analysis" / "Task5_myeloid_subtypes" / "myeloid_combined_obs.csv"
TUMOUR_TISSUES = {"PDAC", "PDAC_tumor", "PDAC_tissue"}

lab = pd.read_csv(ANNOTATION, index_col=0, low_memory=False)
lab = lab.loc[lab["tissue"].astype(str).isin(TUMOUR_TISSUES)]

rows = []
for ds, g in lab.groupby("dataset"):
    rows.append({
        "dataset": ds,
        "task5_n_cells": int(len(g)),
        "task5_n_siglec15_positive": int((pd.to_numeric(g["SIGLEC15"], errors="coerce") > 0).sum()),
        "task5_positive_fraction": float((pd.to_numeric(g["SIGLEC15"], errors="coerce") > 0).mean()),
    })
task5 = pd.DataFrame(rows).set_index("dataset")

sm = pd.read_csv(OUT / "_v04_sample_metrics.csv")
v04 = sm.groupby("dataset")[["n_myeloid", "n_positive"]].sum()
v04["v04_positive_fraction"] = v04["n_positive"] / v04["n_myeloid"]

merged = task5.join(v04, how="outer")
merged["positive_fraction_delta"] = merged["v04_positive_fraction"] - merged["task5_positive_fraction"]

print(merged.to_string())
print()
print("total V04 cells:", int(v04["n_myeloid"].sum()),
      " positive:", int(v04["n_positive"].sum()),
      f" ({100 * v04['n_positive'].sum() / v04['n_myeloid'].sum():.3f}%)")

(OUT / "_v04_layer_crosscheck.json").write_text(
    json.dumps({
        "purpose": "Gate 0 completeness check only; not used in any biological inference.",
        "task5_siglec15_layer": ".raw",
        "v04_siglec15_layer": ".X",
        "per_dataset": merged.reset_index().to_dict(orient="records"),
        "v04_total_cells": int(v04["n_myeloid"].sum()),
        "v04_total_positive": int(v04["n_positive"].sum()),
    }, indent=2, ensure_ascii=False), encoding="utf-8")
