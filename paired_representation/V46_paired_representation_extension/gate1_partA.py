"""Gate 1 (part A) — paired representation test: R1/R2 vs K_i.

R1 = myeloid aggregate SIGLEC15 (patient-level continuous readout, the closest
     analogue of a scalable clinical RNA readout)
R2 = SIGLEC15 RNA-detectable fraction of included tumour-myeloid cells

K_i = locked within-patient detectable-minus-non-detectable ACP5/MMP9 contrast.

All inputs are frozen archive products; no re-estimation of K_i is performed.
Statistics: Spearman rho, dataset-stratified patient bootstrap 95% CI (10,000
resamples, seed 15072026), and leave-one-dataset-out direction stability.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

V45 = Path(r".")
ANALYSIS = Path(r".\biological_validation"
                r"\analysis_F_cross_scale\outputs")
V04 = Path(r".\data\analysis\V04_discriminating_analysis")
OUT = Path(r".\data\analysis\V46_validation_extension\Gate1")
OUT.mkdir(parents=True, exist_ok=True)

SEED = 15072026
N_BOOT = 10000


def stratified_bootstrap_rho(df: pd.DataFrame, xcol: str, ycol: str,
                             seed: int = SEED, n_boot: int = N_BOOT) -> dict:
    """Resample patients within dataset; return rho CI."""
    rng = np.random.default_rng(seed)
    groups = [(g[xcol].to_numpy(float), g[ycol].to_numpy(float))
              for _, g in df.groupby("dataset")]
    vals = []
    for _ in range(n_boot):
        xs, ys = [], []
        for gx, gy in groups:
            idx = rng.integers(0, len(gx), len(gx))
            xs.append(gx[idx])
            ys.append(gy[idx])
        x = np.concatenate(xs)
        y = np.concatenate(ys)
        if len(np.unique(x)) < 3 or len(np.unique(y)) < 3:
            continue
        r = stats.spearmanr(x, y).statistic
        if np.isfinite(r):
            vals.append(r)
    v = np.asarray(vals)
    return {
        "n_boot_valid": int(v.size),
        "ci_low": float(np.percentile(v, 2.5)) if v.size else None,
        "ci_high": float(np.percentile(v, 97.5)) if v.size else None,
        "boot_sd": float(v.std(ddof=1)) if v.size > 1 else None,
    }


def leave_one_dataset_out(df: pd.DataFrame, xcol: str, ycol: str) -> dict:
    out = {}
    for ds in sorted(df["dataset"].unique()):
        sub = df.loc[df["dataset"] != ds]
        if sub[xcol].notna().sum() >= 4:
            out[f"omit_{ds}"] = {
                "n": int(sub[[xcol, ycol]].dropna().shape[0]),
                "rho": float(stats.spearmanr(sub[xcol], sub[ycol],
                                             nan_policy="omit").statistic),
            }
    full = float(stats.spearmanr(df[xcol], df[ycol], nan_policy="omit").statistic)
    signs = [np.sign(v["rho"]) for v in out.values()]
    out["_full_rho"] = full
    out["_direction_retained"] = bool(all(s == np.sign(full) for s in signs)) if signs else None
    out["_any_sign_reversal"] = bool(any(s != np.sign(full) for s in signs)) if signs else None
    return out


def main() -> None:
    metrics = pd.read_csv(V04 / "_v04_sample_metrics.csv")
    aggr = pd.read_csv(ANALYSIS / "myeloid_sample_aggregates_scored.csv")

    df = metrics.merge(
        aggr[["dataset", "sample_id", "SIGLEC15", "n_myeloid", "n_siglec15_positive",
              "Derosiers_ACP5_MMP9_core"]],
        on=["dataset", "sample_id"], how="left", validate="one_to_one")

    # readouts
    df["R1_myeloid_aggregate_SIGLEC15"] = df["SIGLEC15"]
    df["R2_positive_fraction"] = df["A_i_siglec15_pos_fraction"]
    df["R4_allcell_aggregate_SIGLEC15"] = np.nan   # filled by part B (h5ad)
    df["K_i"] = df["K_i_locked_unadjusted"]

    rep: dict = {
        "n_samples_total": int(len(df)),
        "n_with_K_i": int(df["K_i"].notna().sum()),
        "n_with_R1": int(df["R1_myeloid_aggregate_SIGLEC15"].notna().sum()),
        "n_with_R2": int(df["R2_positive_fraction"].notna().sum()),
        "merge_missing_R1": df.loc[df["R1_myeloid_aggregate_SIGLEC15"].isna(),
                                   ["dataset", "sample_id"]].to_dict("records"),
    }

    pairs = {
        "R1_myeloid_aggregate_SIGLEC15_vs_K_i": ("R1_myeloid_aggregate_SIGLEC15", "K_i"),
        "R2_positive_fraction_vs_K_i": ("R2_positive_fraction", "K_i"),
    }
    results = {}
    for name, (xc, yc) in pairs.items():
        sub = df[[ "dataset", "sample_id", xc, yc]].dropna()
        if len(sub) < 5:
            results[name] = {"error": "insufficient n", "n": int(len(sub))}
            continue
        rho, p = stats.spearmanr(sub[xc], sub[yc])
        boot = stratified_bootstrap_rho(sub, xc, yc)
        loo = leave_one_dataset_out(sub, xc, yc)
        # components for each dataset
        per_ds = {
            ds: {"n": int(len(g)), "rho": float(stats.spearmanr(g[xc], g[yc]).statistic)}
            for ds, g in sub.groupby("dataset") if len(g) >= 3
        }
        results[name] = {
            "n": int(len(sub)),
            "rho": float(rho), "p": float(p),
            **boot,
            "leave_one_dataset_out": loo,
            "per_dataset": per_ds,
        }

    rep["pairs"] = results
    df.to_csv(OUT / "GATE1_paired_readouts.csv", index=False)
    (OUT / "_gate1_partA_summary.json").write_text(
        json.dumps(rep, indent=2, ensure_ascii=False), encoding="utf-8")

    print(json.dumps(rep, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
