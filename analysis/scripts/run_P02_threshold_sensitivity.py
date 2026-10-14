#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""P0-2: positivity-threshold sensitivity analysis for the locked S17 pipeline.

Reuses the locked S17 statistical infrastructure (seed, resampling, inferential
frameworks) and varies only the SIGLEC15 positivity definition.

Thresholds are expressed in the log1p-normalised expression space used by S17:
  t = 0.00  locked baseline (any detected expression)
  t = 0.25  remove low-tail positives
  t = 0.50  ~median of positive-cell expression (0.478)
  t = 1.00  ~P75 of positive-cell expression (0.726)
  t = 1.21  ~P90 of positive-cell expression (1.214)

Cell-level contrasts depend directly on the threshold. The patient-level
myeloid aggregate uses the continuous cell-equal mean (log[1 + mean{exp(x)-1}])
and is therefore mathematically invariant to the threshold; it is re-computed
and reported per row as a verification. An alternative aggregate definition
(pos_rate = fraction of cells above threshold) is added as a cross-check.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import run_S18_S19 as s18  # noqa: E402

from analysis_F_cross_scale import (  # noqa: E402
    N_RESAMPLES, PRIMARY, RNG, SEED, TARGET,
    permutation_p, residualize_depth, score_program,
    spearman_bootstrap, stratified_bootstrap_mean, stratified_rank_inference,
)

OUT = HERE / "outputs_P02_threshold_sensitivity"
OUT.mkdir(parents=True, exist_ok=True)

THRESHOLDS = [0.00, 0.25, 0.50, 1.00, 1.21]


def collect_samples():
    """Return per-sample (obs, expr, dataset) for tumour myeloid cells, >=50 cells."""
    samples = []
    for path in s18.SCRNA_PATHS:
        obs, expr, dataset = s18.extract_h5ad_custom(path)
        for sample_id, idx in obs.groupby("sample_id", observed=True).groups.items():
            idx_arr = np.asarray(list(idx), dtype=int)
            if len(idx_arr) < 50:
                continue
            samples.append({
                "dataset": dataset,
                "sample_id": str(sample_id),
                "obs": obs.iloc[idx_arr].reset_index(drop=True),
                "expr": expr.iloc[idx_arr].reset_index(drop=True),
            })
    return samples


def run_cell_level(samples, threshold):
    """Cell-level contrast (patient as unit) for PRIMARY program at a threshold."""
    rows = []
    for s in samples:
        expr = s["expr"]
        positive = expr[TARGET].to_numpy(float) > threshold
        if positive.sum() == 0 or (~positive).sum() == 0:
            continue
        score, present = score_program(expr, PRIMARY)
        if not np.isfinite(score).any():
            continue
        residual = residualize_depth(score, s["obs"])
        rows.append({
            "dataset": s["dataset"], "sample_id": s["sample_id"],
            "n_myeloid": len(expr), "n_positive": int(positive.sum()),
            "n_negative": int((~positive).sum()),
            "n_genes_present": len(present),
            "difference_positive_minus_negative": float(
                np.nanmean(score[positive]) - np.nanmean(score[~positive])
            ),
            "depth_adjusted_difference": float(
                np.nanmean(residual[positive]) - np.nanmean(residual[~positive])
            ),
        })
    effects = pd.DataFrame(rows)
    if effects.empty:
        return None
    values = effects["difference_positive_minus_negative"].to_numpy(float)
    ci_lo, ci_hi = stratified_bootstrap_mean(effects, "difference_positive_minus_negative")
    w = stats.wilcoxon(values, alternative="two-sided", zero_method="wilcox")
    adj = effects["depth_adjusted_difference"].to_numpy(float)
    adj_lo, adj_hi = stratified_bootstrap_mean(effects, "depth_adjusted_difference")
    w_adj = stats.wilcoxon(adj, alternative="two-sided", zero_method="wilcox")
    return {
        "n_samples": len(effects),
        "n_datasets": effects["dataset"].nunique(),
        "n_positive_cells_total": int(effects["n_positive"].sum()),
        "mean_sample_difference": float(values.mean()),
        "ci_low": float(ci_lo), "ci_high": float(ci_hi),
        "wilcoxon_p": float(w.pvalue),
        "depth_adjusted_mean": float(adj.mean()),
        "depth_adjusted_ci_low": float(adj_lo),
        "depth_adjusted_ci_high": float(adj_hi),
        "depth_adjusted_wilcoxon_p": float(w_adj.pvalue),
    }


def run_aggregate(samples, threshold):
    """Patient-level aggregate correlation (primary definition and pos_rate variant)."""
    rows = []
    for s in samples:
        expr = s["expr"]
        aggregate = np.log1p(np.mean(np.expm1(expr.to_numpy(float)), axis=0))
        agg_row = dict(zip(expr.columns, aggregate))
        rows.append({
            "dataset": s["dataset"], "sample_id": s["sample_id"],
            "n_myeloid": len(expr),
            "n_positive": int((expr[TARGET].to_numpy(float) > threshold).sum()),
            "pos_rate": float((expr[TARGET].to_numpy(float) > threshold).mean()),
            "agg_SIGLEC15": float(agg_row[TARGET]),
            "ACP5": float(agg_row["ACP5"]),
            "MMP9": float(agg_row["MMP9"]),
        })
    frame = pd.DataFrame(rows)

    # Score the aggregate values within dataset (identical to S17 aggregate_scored)
    scored = []
    for dataset, group in frame.groupby("dataset"):
        group = group.copy()
        group["agg_ACP5_MMP9"], _ = score_program(group, PRIMARY)
        scored.append(group)
    scored = pd.concat(scored, ignore_index=True)

    rho, lo, hi, p = stratified_rank_inference(scored, "agg_SIGLEC15", "agg_ACP5_MMP9")
    rho_pr, p_pr, lo_pr, hi_pr = spearman_bootstrap(
        scored["pos_rate"].to_numpy(float),
        scored["agg_ACP5_MMP9"].to_numpy(float),
    )
    pp_pr = permutation_p(
        scored["pos_rate"].to_numpy(float),
        scored["agg_ACP5_MMP9"].to_numpy(float),
        rho_pr,
    )
    return {
        "n_samples_aggregate": len(scored),
        "agg_rho": float(rho), "agg_ci_low": float(lo), "agg_ci_high": float(hi),
        "agg_permutation_p": float(p),
        "pos_rate_rho": float(rho_pr), "pos_rate_ci_low": float(lo_pr),
        "pos_rate_ci_high": float(hi_pr), "pos_rate_permutation_p": float(pp_pr),
    }


def main():
    print(f"SEED = {SEED}  N_RESAMPLES = {N_RESAMPLES}")
    print(f"Thresholds (log1p space): {THRESHOLDS}")
    samples = collect_samples()
    print(f"Eligible samples (>=50 myeloid): {len(samples)}")

    rows = []
    for t in THRESHOLDS:
        cell = run_cell_level(samples, t)
        agg = run_aggregate(samples, t)
        if cell is None:
            print(f"\nt = {t}: no estimable sample contrast")
            rows.append({"threshold": t, "cell_n_samples": 0})
            continue
        row = {
            "threshold": t,
            "cell_n_samples": cell["n_samples"],
            "cell_n_datasets": cell["n_datasets"],
            "cell_n_positive_total": cell["n_positive_cells_total"],
            "cell_mean_difference": cell["mean_sample_difference"],
            "cell_ci_low": cell["ci_low"], "cell_ci_high": cell["ci_high"],
            "cell_wilcoxon_p": cell["wilcoxon_p"],
            "cell_depth_adjusted_mean": cell["depth_adjusted_mean"],
            "cell_depth_ci_low": cell["depth_adjusted_ci_low"],
            "cell_depth_ci_high": cell["depth_adjusted_ci_high"],
            "cell_depth_wilcoxon_p": cell["depth_adjusted_wilcoxon_p"],
            "agg_n_samples": agg["n_samples_aggregate"],
            "agg_rho": agg["agg_rho"], "agg_ci_low": agg["agg_ci_low"],
            "agg_ci_high": agg["agg_ci_high"], "agg_permutation_p": agg["agg_permutation_p"],
            "pos_rate_rho": agg["pos_rate_rho"],
            "pos_rate_ci_low": agg["pos_rate_ci_low"],
            "pos_rate_ci_high": agg["pos_rate_ci_high"],
            "pos_rate_permutation_p": agg["pos_rate_permutation_p"],
        }
        rows.append(row)
        print(f"\n--- t = {t} ---")
        for k, v in row.items():
            print(f"  {k}: {v}")

    table = pd.DataFrame(rows)
    table.to_csv(OUT / "threshold_sensitivity_table.csv", index=False)
    (OUT / "threshold_sensitivity_summary.json").write_text(
        json.dumps({
            "analysis": "P0-2 positivity-threshold sensitivity (S17 locked pipeline)",
            "seed": SEED, "n_resamples": N_RESAMPLES,
            "thresholds_log1p": THRESHOLDS,
            "primary_program": PRIMARY,
            "locked_aggregate_definition": "log[1 + mean{exp(x)-1}] over all tumour myeloid cells",
            "results": rows,
        }, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"\n-> {OUT / 'threshold_sensitivity_table.csv'}")


if __name__ == "__main__":
    main()
