"""Gate 1 (part C) — final paired representation test: R1-R4 vs K_i.

R1 myeloid aggregate SIGLEC15        (closest analogue of a scalable RNA readout)
R2 SIGLEC15-detectable fraction
R3 mean SIGLEC15 intensity among detectable cells
R4 all-cell pseudo-whole-tumour aggregate

Statistics per readout: Spearman rho, dataset-stratified patient bootstrap 95% CI
(10,000 resamples, seed 15072026), leave-one-dataset-out direction stability,
and per-dataset rho.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

V04 = Path(r".\data\analysis\V04_discriminating_analysis")
G1 = Path(r".\data\analysis\V46_validation_extension\Gate1")
SEED, N_BOOT = 15072026, 10000


def boot_rho(df, xc, yc):
    rng = np.random.default_rng(SEED)
    groups = [(g[xc].to_numpy(float), g[yc].to_numpy(float))
              for _, g in df.groupby("dataset")]
    vals = []
    for _ in range(N_BOOT):
        xs, ys = [], []
        for gx, gy in groups:
            i = rng.integers(0, len(gx), len(gx))
            xs.append(gx[i]); ys.append(gy[i])
        x, y = np.concatenate(xs), np.concatenate(ys)
        if len(np.unique(x)) < 3 or len(np.unique(y)) < 3:
            continue
        r = stats.spearmanr(x, y).statistic
        if np.isfinite(r):
            vals.append(r)
    v = np.asarray(vals)
    return {"n_boot_valid": int(v.size),
            "ci_low": float(np.percentile(v, 2.5)) if v.size else None,
            "ci_high": float(np.percentile(v, 97.5)) if v.size else None}


def loo(df, xc, yc):
    full = float(stats.spearmanr(df[xc], df[yc], nan_policy="omit").statistic)
    out, signs = {}, []
    for ds in sorted(df["dataset"].unique()):
        s = df.loc[df["dataset"] != ds]
        if s[[xc, yc]].dropna().shape[0] >= 4:
            r = float(stats.spearmanr(s[xc], s[yc], nan_policy="omit").statistic)
            out[f"omit_{ds}"] = {"n": int(s[[xc, yc]].dropna().shape[0]), "rho": r}
            signs.append(np.sign(r))
    out["_full_rho"] = full
    out["_direction_retained"] = bool(signs) and all(s == np.sign(full) for s in signs)
    out["_any_sign_reversal"] = bool(signs) and any(s != np.sign(full) for s in signs)
    return out


def main() -> None:
    metrics = pd.read_csv(V04 / "_v04_sample_metrics.csv")
    pB = pd.read_csv(G1 / "GATE1_R3_R4_by_sample.csv")
    aggr = pd.read_csv(Path(r"."
                            r"\biological_validation\analysis_F_cross_scale\outputs"
                            r"\myeloid_sample_aggregates_scored.csv"))[
        ["dataset", "sample_id", "SIGLEC15"]]

    df = (metrics[["dataset", "sample_id", "n_myeloid", "n_positive",
                   "A_i_siglec15_pos_fraction", "K_i_locked_unadjusted"]]
          .merge(aggr.rename(columns={"SIGLEC15": "R1"}),
                 on=["dataset", "sample_id"], how="left")
          .merge(pB[["dataset", "sample_id", "R3_mean_positive_intensity",
                     "R4_allcell_aggregate", "n_all_cells_tumour"]],
                 on=["dataset", "sample_id"], how="left"))
    df["R2"] = df["A_i_siglec15_pos_fraction"]
    df["R3"] = df["R3_mean_positive_intensity"]
    df["R4"] = df["R4_allcell_aggregate"]
    df["K_i"] = df["K_i_locked_unadjusted"]
    # negative control: number of detectable cells (should not drive a positive relation)
    df["R_neg_npos"] = df["n_positive"]

    readouts = {
        "R1_myeloid_aggregate": ("R1", "myeloid aggregate SIGLEC15"),
        "R2_detectable_fraction": ("R2", "SIGLEC15-detectable cell fraction"),
        "R3_positive_intensity": ("R3", "mean SIGLEC15 among detectable cells"),
        "R4_allcell_aggregate": ("R4", "all-cell pseudo-whole-tumour aggregate"),
    }

    rep = {"n_samples": int(len(df)), "n_with_K_i": int(df.K_i.notna().sum()),
           "readouts": {}, "coverage": {}}
    for name, (col, desc) in readouts.items():
        sub = df[["dataset", "sample_id", col, "K_i"]].dropna()
        rep["coverage"][name] = {"n": int(len(sub)),
                                 "n_per_dataset": sub.groupby("dataset").size().to_dict()}
        if len(sub) < 5:
            rep["readouts"][name] = {"description": desc, "error": "insufficient n"}
            continue
        rho, p = stats.spearmanr(sub[col], sub["K_i"])
        b = boot_rho(sub, col, "K_i")
        rep["readouts"][name] = {
            "description": desc, "n": int(len(sub)),
            "rho": float(rho), "p": float(p), **b,
            "ci_excludes_moderate_positive": bool(b["ci_high"] is not None
                                                  and b["ci_high"] < 0.30),
            "leave_one_dataset_out": loo(sub, col, "K_i"),
            "per_dataset": {ds: {"n": int(len(g)),
                                 "rho": float(stats.spearmanr(g[col], g["K_i"]).statistic)}
                            for ds, g in sub.groupby("dataset") if len(g) >= 3},
        }

    df.to_csv(G1 / "GATE1_all_readouts_vs_Ki.csv", index=False)
    (G1 / "_gate1_final_summary.json").write_text(
        json.dumps(rep, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"samples={rep['n_samples']}  with K_i={rep['n_with_K_i']}\n")
    for name, r in rep["readouts"].items():
        if "error" in r:
            print(f"{name}: {r['error']}"); continue
        print(f"{name}  (n={r['n']})")
        print(f"   rho={r['rho']:+.3f}  CI [{r['ci_low']:+.3f}, {r['ci_high']:+.3f}]  p={r['p']:.3f}")
        print(f"   LOO retained={r['leave_one_dataset_out']['_direction_retained']}  "
              f"reversal={r['leave_one_dataset_out']['_any_sign_reversal']}")
        print(f"   per-dataset: " + ", ".join(
            f"{k}={v['rho']:+.2f}(n={v['n']})" for k, v in r["per_dataset"].items()))
        print()


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
