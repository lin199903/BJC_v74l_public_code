"""Coupling Audit, part 1 (CSV-based): patient-specific K_i distribution and correlates."""
import json

import numpy as np
import pandas as pd
from scipy import stats as st

J = r"./prior_analysis_v2"
RES = J + r"\deposit_v24\results"

core = pd.read_csv(RES + r"\cell_level_sample_effects.csv")
core = core[core.program == "Derosiers_ACP5_MMP9_core"].copy()
matched = pd.read_csv(RES + r"\depth_matched_cell_effects.csv")
matched = matched[matched.program == "Derosiers_ACP5_MMP9_core"].copy()

core["prevalence"] = core.n_positive / core.n_myeloid
core["pos_per_1k"] = core.n_positive / core.n_myeloid * 1000

report = {}
# 1. distribution of K_i (unadjusted)
K = core.difference_positive_minus_negative
report["n_samples"] = int(len(core))
report["n_positive_ki"] = int((K > 0).sum())
report["n_negative_ki"] = int((K < 0).sum())
report["mean_ki"] = float(K.mean())
report["median_ki"] = float(K.median())
report["sd_ki"] = float(K.std(ddof=1))
report["iqr_ki"] = [float(K.quantile(0.25)), float(K.quantile(0.75))]
report["range_ki"] = [float(K.min()), float(K.max())]
report["n_ki_above_0_5"] = int((K > 0.5).sum())
report["n_ki_below_0"] = int((K < 0).sum())

# 2. dataset-specific medians
ds = core.groupby("dataset").difference_positive_minus_negative.agg(["count", "median", "mean", "min", "max"])
report["by_dataset"] = {k: {"n": int(r["count"]), "median": float(r["median"]),
                            "mean": float(r["mean"]), "min": float(r["min"]), "max": float(r["max"])}
                        for k, r in ds.iterrows()}

# 3. rank stability: unadjusted vs depth-adjusted
m = core.dropna(subset=["difference_positive_minus_negative", "depth_adjusted_difference"])
r_s, p_s = st.spearmanr(m.difference_positive_minus_negative, m.depth_adjusted_difference)
report["rank_stability_unadj_vs_depthadj"] = {"n": int(len(m)), "spearman": float(r_s), "p": float(p_s)}

# 4. depth-matched subset rank stability
mm = core.merge(matched[["dataset", "sample_id", "difference_positive_minus_negative"]],
                on=["dataset", "sample_id"], suffixes=("", "_matched"))
mm = mm.dropna(subset=["difference_positive_minus_negative", "difference_positive_minus_negative_matched"])
r_m, p_m = st.spearmanr(mm.difference_positive_minus_negative, mm.difference_positive_minus_negative_matched)
report["rank_stability_unadj_vs_matched"] = {"n": int(len(mm)), "spearman": float(r_m), "p": float(p_m)}

# 5. K_i vs prevalence
r_p, p_p = st.spearmanr(core.prevalence, K)
report["ki_vs_prevalence"] = {"n": int(len(core)), "spearman": float(r_p), "p": float(p_p)}

# 6. K_i vs positive-cell count
r_c, p_c = st.spearmanr(core.n_positive, K)
report["ki_vs_positive_count"] = {"n": int(len(core)), "spearman": float(r_c), "p": float(p_c)}

# 7. high-prevalence patients: coupling direction
hi = core[core.prevalence > core.prevalence.median()]
lo = core[core.prevalence <= core.prevalence.median()]
report["coupling_by_prevalence_median_split"] = {
    "high_prevalence_median_ki": float(hi.difference_positive_minus_negative.median()),
    "low_prevalence_median_ki": float(lo.difference_positive_minus_negative.median()),
    "n_high": int(len(hi)), "n_low": int(len(lo)),
}

# 8. per-sample table
per = core[["dataset", "sample_id", "n_myeloid", "n_positive", "prevalence",
            "difference_positive_minus_negative", "depth_adjusted_difference"]].sort_values(
    ["dataset", "sample_id"])
per = per.merge(matched[["dataset", "sample_id", "difference_positive_minus_negative"]]
                .rename(columns={"difference_positive_minus_negative": "matched_difference"}),
                on=["dataset", "sample_id"], how="left")
per.to_csv(RES + r"\COUPLING_AUDIT_ki_by_patient.csv", index=False)

with open(RES + r"\COUPLING_AUDIT_summary.json", "w", encoding="utf-8") as fh:
    json.dump(report, fh, indent=2, default=str)

print(json.dumps(report, indent=2))
