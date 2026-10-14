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
matched["sample_id"] = matched["sample_id"].astype(str).str.replace("b'", "").str.replace("'", "")

K = core.difference_positive_minus_negative.to_numpy(float)
N = 10_000

report = {}

loo_means = [float(np.mean(np.delete(K, i))) for i in range(len(K))]
report["locked_mean"] = float(np.mean(K))
report["n"] = int(len(K))
report["median"] = float(np.median(K))
report["iqr"] = [float(np.percentile(K, 25)), float(np.percentile(K, 75))]
report["mean_excluding_max"] = float(np.mean(K[K < K.max()]))
report["max_value"] = float(K.max())
report["max_z"] = float((K.max() - K.mean()) / K.std(ddof=1))
report["loo_mean_min"] = float(min(loo_means))
report["loo_mean_max"] = float(max(loo_means))
report["loo_mean_median"] = float(np.median(loo_means))
report["loo_mean_when_dropping_max"] = loo_means[int(np.argmax(K))]


def boot_ci(vals, datasets, seed):
    r = np.random.default_rng(seed)
    dset = pd.Series(datasets).to_numpy()
    out = []
    for _ in range(N):
        parts = []
        for ds in np.unique(dset):
            idx = np.where(dset == ds)[0]
            pick = r.choice(idx, size=len(idx), replace=True)
            parts.append(vals[pick].mean())
        out.append(np.mean(parts))
    out = np.asarray(out)
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


report["bootstrap_ci_full"] = boot_ci(K, core.dataset, 15072026)
mask = K < K.max()
report["bootstrap_ci_excluding_max"] = boot_ci(K[mask], core.dataset[mask], 15072026)

r_adj, p_adj = st.spearmanr(K, core.depth_adjusted_difference.to_numpy(float))
report["rank_raw_vs_depthadjusted"] = {"spearman": float(r_adj), "p": float(p_adj), "label": "consistency check"}
mm = core.merge(matched[["dataset", "sample_id", "difference_positive_minus_negative"]],
                on=["dataset", "sample_id"], suffixes=("", "_matched")).dropna(
    subset=["difference_positive_minus_negative", "difference_positive_minus_negative_matched"])
r_m, p_m = st.spearmanr(mm.difference_positive_minus_negative, mm.difference_positive_minus_negative_matched)
report["rank_raw_vs_matched"] = {"n": int(len(mm)), "spearman": float(r_m), "p": float(p_m),
                                 "label": "primary robustness evidence"}
mm2 = mm[mm.difference_positive_minus_negative < mm.difference_positive_minus_negative.max()]
r_m2, p_m2 = st.spearmanr(mm2.difference_positive_minus_negative, mm2.difference_positive_minus_negative_matched)
report["rank_raw_vs_matched_excluding_max"] = {"n": int(len(mm2)), "spearman": float(r_m2), "p": float(p_m2)}


def r_ci(r, n):
    z = np.arctanh(r)
    se = 1 / np.sqrt(n - 3)
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


core["prevalence"] = core.n_positive / core.n_myeloid
h5 = json.load(open(r"./temp\COUPLING_AUDIT_h5ad.json", encoding="utf-8"))
rows = []
for ds, recs in h5.items():
    if ds == "obs_columns":
        continue
    rows.extend(recs)
inten = pd.DataFrame(rows).dropna(subset=["mean_positive_intensity"])
m3 = core.merge(inten[["dataset", "sample_id", "mean_positive_intensity"]],
                on=["dataset", "sample_id"], how="left").dropna(subset=["mean_positive_intensity"])
r_p, p_p = st.spearmanr(core.prevalence, K)
r_i, p_i = st.spearmanr(m3.mean_positive_intensity, m3.difference_positive_minus_negative)
report["ki_vs_prevalence"] = {"n": int(len(core)), "spearman": float(r_p), "p": float(p_p),
                              "ci95": r_ci(r_p, len(core))}
report["ki_vs_intensity"] = {"n": int(len(m3)), "spearman": float(r_i), "p": float(p_i),
                             "ci95": r_ci(r_i, len(m3))}
report["variance_interpretation"] = ("Observed dispersion combines potential biological heterogeneity "
                                     "with heterogeneous measurement error; current data cannot separate them.")

with open(RES + r"\COUPLING_ROBUSTNESS_GATE.json", "w", encoding="utf-8") as fh:
    json.dump(report, fh, indent=2)
print(json.dumps(report, indent=2))
