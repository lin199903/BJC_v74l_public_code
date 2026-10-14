"""S21 ambient-robustness analysis: deep-cell vs shallow-cell contrasts.

The processed objects store only log1p-normalised expression, so SoupX/CellBender
cannot be run (no raw counts). Ambient contamination inflates the transcript
content of low-library-size cells most. We therefore re-estimate the Level 1
core contrast separately within deep cells (total_counts >= sample median) and
shallow cells (below the median), restricting BOTH positive and negative groups
to the same depth stratum. If ambient co-contamination of low-depth cells drove
the contrast, the shallow-stratum contrast should dominate; if the state is
biological, the deep-stratum contrast should persist. Outputs to a writable TEMP
dir; summary printed to stdout.
"""
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ANALYSIS = Path(r".\biological_validation\analysis_F_cross_scale")
sys.path.insert(0, str(ANALYSIS))
import run_S20_state_controls as S
import analysis_F_cross_scale as F

PROG = {"anchor": "SIGLEC15", "genes": ["ACP5", "MMP9"]}


def main():
    out_dir = Path(tempfile.gettempdir()) / "S21_outputs"
    out_dir.mkdir(exist_ok=True)
    genes = sorted(set(PROG["genes"] + [PROG["anchor"]]))
    rows = []
    for path in S.SCRNA_PATHS:
        obs, expr, dataset = S.extract_myeloid(path, genes)
        for sample_id, idx in obs.groupby("sample_id", observed=True).groups.items():
            idx = np.asarray(list(idx), dtype=int)
            if len(idx) < 50:
                continue
            sample_obs = obs.iloc[idx].reset_index(drop=True)
            sample_expr = expr.iloc[idx].reset_index(drop=True)
            if PROG["anchor"] not in sample_expr or "total_counts" not in sample_obs:
                continue
            x = sample_expr[PROG["anchor"]].to_numpy(float)
            pos = x > 0
            if pos.sum() == 0 or (~pos).sum() == 0:
                continue
            tc = pd.to_numeric(sample_obs["total_counts"], errors="coerce").to_numpy(float)
            median_tc = np.nanmedian(tc)
            deep = tc >= median_tc
            shallow = tc < median_tc
            present = [g for g in PROG["genes"] if g in sample_expr.columns]
            z = F.zscore_columns(sample_expr[present].to_numpy(float))
            score = np.nanmean(z, axis=1)
            for name, mask in (("deep", deep), ("shallow", shallow)):
                p, n = pos & mask, (~pos) & mask
                if p.sum() == 0 or n.sum() == 0:
                    continue
                eff = np.nanmean(score[p]) - np.nanmean(score[n])
                rows.append({
                    "dataset": dataset, "sample_id": str(sample_id),
                    "stratum": name, "n_positive": int(p.sum()),
                    "n_negative": int(n.sum()),
                    "difference_positive_minus_negative": float(eff),
                })
    df = pd.DataFrame(rows)
    report = {}
    for stratum, g in df.groupby("stratum"):
        v = g["difference_positive_minus_negative"].to_numpy(float)
        ci = F.stratified_bootstrap_mean(g, "difference_positive_minus_negative")
        wilc = F.stats.wilcoxon(v, alternative="two-sided", zero_method="wilcox")
        report[stratum] = {
            "n_samples": int(len(g)), "n_positive_cells": int(g["n_positive"].sum()),
            "mean_difference": float(v.mean()), "ci_low": float(ci[0]),
            "ci_high": float(ci[1]), "wilcoxon_p": float(wilc.pvalue),
        }
    paired = (
        df.pivot(index=["dataset", "sample_id"], columns="stratum",
                 values="difference_positive_minus_negative")
        .dropna(subset=["deep", "shallow"])
        .reset_index()
    )
    paired["deep_minus_shallow"] = paired["deep"] - paired["shallow"]
    paired_ci = F.stratified_bootstrap_mean(paired, "deep_minus_shallow")
    paired_wilc = F.stats.wilcoxon(
        paired["deep_minus_shallow"].to_numpy(float),
        alternative="two-sided", zero_method="wilcox"
    )
    report["deep_minus_shallow"] = {
        "n_paired_samples": int(len(paired)),
        "mean_difference": float(paired["deep_minus_shallow"].mean()),
        "ci_low": float(paired_ci[0]),
        "ci_high": float(paired_ci[1]),
        "wilcoxon_p": float(paired_wilc.pvalue),
    }
    df.to_csv(out_dir / "S21_depth_strata.csv", index=False)
    paired.to_csv(out_dir / "S21_depth_strata_paired.csv", index=False)
    with open(out_dir / "S21_summary.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)
    print("S21_SUMMARY_JSON_BEGIN")
    print(json.dumps(report, indent=2))
    print("S21_SUMMARY_JSON_END")


if __name__ == "__main__":
    main()
