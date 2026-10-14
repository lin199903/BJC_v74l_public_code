"""S21b + L1 series heterogeneity: reviewer-driven sensitivity statistics.

1. Deep-vs-shallow formal contrast (S21b): equally weighted mean difference
   between deep and shallow stratum sample effects, with a dataset-stratified
   bootstrap CI (10,000 resamples, seed 15072026) and a paired comparison on
   the sample subset contributing both strata.
2. Level 1 series heterogeneity: dataset-specific means of the unadjusted
   core sample effects, with ranges, for reporting in Results.

Outputs: S21b_summary.json + L1_series_heterogeneity.json in the deposit
results directory; summary printed to stdout.
"""
import json

import numpy as np
import pandas as pd

BASE = r"./prior_analysis_v2\deposit_v24\results"
RNG = np.random.default_rng(15072026)
N_RES = 10_000


def main():
    df = pd.read_csv(BASE + r"\S21_depth_strata.csv")
    deep = df[df.stratum == "deep"]
    shal = df[df.stratum == "shallow"]

    d_mean = deep.difference_positive_minus_negative.mean()
    s_mean = shal.difference_positive_minus_negative.mean()
    diff_obs = d_mean - s_mean

    # stratified bootstrap (dataset-stratified, unpaired)
    diffs = []
    for _ in range(N_RES):
        ds = []
        for dataset, g in df.groupby("dataset"):
            n_deep = (g.stratum == "deep").sum()
            n_shal = (g.stratum == "shallow").sum()
            idx = RNG.choice(len(g), size=len(g), replace=True)
            gb = g.iloc[idx]
            dd = gb[gb.stratum == "deep"].difference_positive_minus_negative
            ss = gb[gb.stratum == "shallow"].difference_positive_minus_negative
            if len(dd) and len(ss):
                ds.append(dd.mean() - ss.mean())
        diffs.append(np.mean(ds))
    diffs = np.asarray(diffs)
    ci_lo, ci_hi = np.percentile(diffs, [2.5, 97.5])
    p_ci = float(np.mean(diffs <= 0) * 2 if diff_obs > 0 else np.mean(diffs >= 0) * 2)

    # paired subset (samples with both strata)
    piv = df.pivot_table(index=["dataset", "sample_id"], columns="stratum",
                         values="difference_positive_minus_negative").dropna()
    paired = piv["deep"] - piv["shallow"]
    paired_mean = float(paired.mean())
    pdiffs = []
    for _ in range(N_RES):
        rows = piv.iloc[RNG.choice(len(piv), size=len(piv), replace=True)]
        pdiffs.append((rows["deep"] - rows["shallow"]).mean())
    pdiffs = np.asarray(pdiffs)
    p_lo, p_hi = np.percentile(pdiffs, [2.5, 97.5])
    from scipy import stats as st
    wilc = st.wilcoxon(paired, zero_method="wilcox") if len(paired) >= 6 else None

    report = {
        "deep_mean": float(d_mean),
        "shallow_mean": float(s_mean),
        "difference_deep_minus_shallow": float(diff_obs),
        "unpaired_stratified_bootstrap_ci": [float(ci_lo), float(ci_hi)],
        "unpaired_n_deep": int(len(deep)),
        "unpaired_n_shallow": int(len(shal)),
        "paired_n_samples": int(len(paired)),
        "paired_mean_difference": paired_mean,
        "paired_bootstrap_ci": [float(p_lo), float(p_hi)],
        "paired_wilcoxon_p": float(wilc.pvalue) if wilc else None,
    }
    with open(BASE + r"\S21b_summary.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)

    # L1 series heterogeneity
    cells = pd.read_csv(BASE + r"\cell_level_sample_effects.csv")
    core = cells[cells.program == "Derosiers_ACP5_MMP9_core"]
    series = {}
    for dataset, g in core.groupby("dataset"):
        v = g.difference_positive_minus_negative.to_numpy(float)
        series[dataset] = {
            "n_samples": int(len(g)),
            "mean": float(v.mean()),
            "sd": float(v.std(ddof=1)) if len(g) > 1 else None,
            "min": float(v.min()),
            "max": float(v.max()),
        }
    with open(BASE + r"\L1_series_heterogeneity.json", "w", encoding="utf-8") as fh:
        json.dump(series, fh, indent=2)

    print("S21B_JSON_BEGIN")
    print(json.dumps(report, indent=2))
    print("S21B_JSON_END")
    print("L1_SERIES_BEGIN")
    print(json.dumps(series, indent=2))
    print("L1_SERIES_END")


if __name__ == "__main__":
    main()
