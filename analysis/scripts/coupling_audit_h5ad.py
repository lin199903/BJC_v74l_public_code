"""Coupling Audit part 2 (h5ad): per-sample SIGLEC15 intensity and TREM2_TAM fraction."""
import json
import sys

import numpy as np
import pandas as pd

ANALYSIS = r".\biological_validation\analysis_F_cross_scale"
sys.path.insert(0, ANALYSIS)
import run_S20_state_controls as S

out = {}
for path in S.SCRNA_PATHS:
    obs, expr, dataset = S.extract_myeloid(path, ["SIGLEC15"])
    # print obs columns for subtype annotation
    cols = list(obs.columns)
    out.setdefault("obs_columns", {})[dataset] = cols
    sub_col = None
    for cand in ["cell_type2", "fine_subtype", "subtype", "annotation", "cell_type_fine", "myeloid_subtype"]:
        if cand in cols:
            sub_col = cand
            break
    rows = []
    for sample_id, idx in obs.groupby("sample_id", observed=True).groups.items():
        idx = np.asarray(list(idx), dtype=int)
        if len(idx) < 50:
            continue
        o = obs.iloc[idx]
        x = expr.iloc[idx]["SIGLEC15"].to_numpy(float)
        pos = x > 0
        n_pos = int(pos.sum())
        intensity = float(x[pos].mean()) if n_pos else None
        trem2_frac = None
        if sub_col is not None:
            sub = o[sub_col].astype(str)
            trem2_mask = sub.str.contains("TREM2", case=False, na=False)
            if trem2_mask.sum() > 0:
                trem2_frac = float(trem2_mask.sum() / len(o))
        rows.append({"dataset": dataset, "sample_id": str(sample_id),
                     "n_myeloid": int(len(o)), "n_positive": n_pos,
                     "mean_positive_intensity": intensity,
                     "trem2_tam_fraction": trem2_frac})
    out[dataset] = rows

with open(r"./temp\COUPLING_AUDIT_h5ad.json", "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=2)
print("COUPLING_AUDIT_H5AD_BEGIN")
print(json.dumps({k: v for k, v in out.items() if k != "obs_columns"}, indent=2))
print("OBS_COLUMNS", json.dumps(out.get("obs_columns", {})))
print("COUPLING_AUDIT_H5AD_END")
