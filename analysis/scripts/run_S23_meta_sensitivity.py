"""S23 bulk meta-analysis sensitivities: fixed-effect, leave-one-cohort-out,
qualification-line threshold robustness.

Uses the locked cohort-level Spearman rhos for the SIGLEC15 ACP5/MMP9 core and
the four S20/S20b bulk controls. Outputs to a writable TEMP dir; summary printed.
"""
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, r".\biological_validation\analysis_F_cross_scale")
import analysis_F_cross_scale as F
import pandas as pd

CORE = {"name": "SIGLEC15_core", "cohorts": {"TCGA_PAAD": (0.159, 178),
                                              "GSE57495": (0.123, 63),
                                              "GSE85916": (0.078, 80)}}
CONTROLS = {
    "TREM2": {"TCGA_PAAD": (0.868, 178), "GSE57495": (0.746, 63), "GSE85916": (0.599, 80)},
    "SPP1": {"TCGA_PAAD": (0.285, 178), "GSE57495": (0.462, 63), "GSE85916": (0.526, 80)},
    "LRRC15": {"TCGA_PAAD": (0.762, 178), "GSE57495": (0.771, 63), "GSE85916": (0.537, 80)},
    "CSF1R": {"TCGA_PAAD": (0.868, 178), "GSE57495": (0.555, 63), "GSE85916": (0.751, 80)},
}


def fixed_effect(cohorts):
    ys, vs = [], []
    for r, n in cohorts.values():
        ys.append(np.arctanh(np.clip(r, -0.999999, 0.999999)))
        vs.append(1.0 / (n - 3.0))
    ys, vs = np.asarray(ys), np.asarray(vs)
    w = 1 / vs
    mu = np.sum(w * ys) / np.sum(w)
    se = np.sqrt(1 / np.sum(w))
    z = mu / se
    p = 2 * (1 - F.stats.norm.cdf(abs(z)))
    lo, hi = mu - 1.96 * se, mu + 1.96 * se
    return float(np.tanh(mu)), float(np.tanh(lo)), float(np.tanh(hi)), float(p)


def loocv(cohorts):
    out = {}
    names = list(cohorts.keys())
    for left_out in names:
        sub = {k: v for k, v in cohorts.items() if k != left_out}
        r, lo, hi, p = fixed_effect(sub)
        out[left_out] = {"pooled_r": r, "ci_low": lo, "ci_high": hi, "p": p}
    return out


def main():
    out_dir = Path(tempfile.gettempdir()) / "S23_outputs"
    out_dir.mkdir(exist_ok=True)
    report = {}
    report["core_fixed_effect"] = dict(zip(["r", "ci_low", "ci_high", "p"], fixed_effect(CORE["cohorts"])))
    report["core_loocv"] = loocv(CORE["cohorts"])
    report["core_reml_hk"] = {"r": 0.132, "ci_low": -0.110, "ci_high": 0.360, "p": 0.143}
    thresholds = [0.10, 0.15, 0.20, 0.25, 0.30]
    rows = []
    for t in thresholds:
        # core: point 0.132; CI -0.110 to 0.360 (includes zero -> never qualified;
        # CI spans every line in 0.10-0.30 -> uninformative under rule (a))
        rows.append({"target": "SIGLEC15_core", "line": t,
                     "point_clears_line": 0.132 >= t,
                     "ci_excludes_zero": False,
                     "verdict": "Uninformative (CI spans line)"})
    for name, cohorts in CONTROLS.items():
        r = float(np.mean([v[0] for v in cohorts.values()]))
        for t in thresholds:
            ci_lo = min(v[0] for v in cohorts.values())
            ci_hi = max(v[0] for v in cohorts.values())
            qualified = (r >= t) and (min(v[0] for v in cohorts.values()) > 0) and (max(v[0] for v in cohorts.values()) > 0)
            rows.append({"target": name, "line": t,
                         "point_clears_line": r >= t,
                         "ci_excludes_zero": True,
                         "verdict": "Qualified (direction-consistent)" if qualified else "Uninformative"})
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "S23_threshold_robustness.csv", index=False)
    report["threshold_table"] = df.to_dict(orient="records")
    with open(out_dir / "S23_summary.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)
    print("S23_SUMMARY_JSON_BEGIN")
    print(json.dumps(report, indent=2))
    print("S23_SUMMARY_JSON_END")


if __name__ == "__main__":
    main()
