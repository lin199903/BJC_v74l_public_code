import json

import numpy as np
import pandas as pd
from scipy import stats as st

BASE = r"./prior_analysis_v2\deposit_v24\results"
RNG = np.random.default_rng(15072026)
N_RES = 10_000

df = pd.read_csv(BASE + r"\S21_depth_strata.csv")
piv = df.pivot_table(index=["dataset", "sample_id"], columns="stratum",
                     values="difference_positive_minus_negative").dropna()
piv["diff"] = piv["deep"] - piv["shallow"]

# dataset-stratified bootstrap of the paired difference
diffs = []
for _ in range(N_RES):
    parts = []
    for dataset, g in piv.groupby(level="dataset"):
        idx = RNG.choice(len(g), size=len(g), replace=True)
        parts.append(g["diff"].iloc[idx].mean())
    diffs.append(np.mean(parts))
d = np.asarray(diffs)
lo, hi = np.percentile(d, [2.5, 97.5])

p = st.wilcoxon(piv["diff"], zero_method="wilcox").pvalue
print("stratified paired: n=%d mean=%.4f ci=%.4f to %.4f p=%.4f"
      % (len(piv), piv["diff"].mean(), lo, hi, p))

# save canonical S21b summary (paired stratified + unpaired stratified)
deep = df[df.stratum == "deep"].difference_positive_minus_negative
shal = df[df.stratum == "shallow"].difference_positive_minus_negative
udiffs = []
for _ in range(N_RES):
    parts = []
    for dataset, g in df.groupby("dataset"):
        idx = RNG.choice(len(g), size=len(g), replace=True)
        gb = g.iloc[idx]
        dd = gb[gb.stratum == "deep"].difference_positive_minus_negative
        ss = gb[gb.stratum == "shallow"].difference_positive_minus_negative
        if len(dd) and len(ss):
            parts.append(dd.mean() - ss.mean())
    udiffs.append(np.mean(parts))
u = np.asarray(udiffs)
ulo, uhi = np.percentile(u, [2.5, 97.5])

report = {
    "deep_mean": float(deep.mean()),
    "shallow_mean": float(shal.mean()),
    "unpaired_difference_deep_minus_shallow": float(deep.mean() - shal.mean()),
    "unpaired_stratified_bootstrap_ci": [float(ulo), float(uhi)],
    "paired_n_samples": int(len(piv)),
    "paired_mean_difference": float(piv["diff"].mean()),
    "paired_stratified_bootstrap_ci": [float(lo), float(hi)],
    "paired_wilcoxon_p": float(p),
}
with open(BASE + r"\S21b_summary.json", "w", encoding="utf-8") as fh:
    json.dump(report, fh, indent=2)
print("S21b_summary.json updated")
