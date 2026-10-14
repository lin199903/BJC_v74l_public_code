"""S22 dilution simulation: does a translating state survive ~2% compartment dilution?

TREM2-centred TAM programme (a state that qualifies at natural prevalence) is
re-aggregated after subsampling TREM2-positive cells to 2% (and 1%) of each
sample's myeloid compartment, matching SIGLEC15's observed prevalence (~2%).
For each dilution level we run 30 Monte Carlo subsamples (fixed seed) and report
the pooled within-series rank association between the TREM2 anchor aggregate and
the leave-anchor-out programme score, computed exactly as in the locked
aggregate layer. Natural SIGLEC15 core (r = -0.45, n = 20) is the reference.
Outputs go to a writable TEMP dir; summary printed to stdout.
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

PROG = {"anchor": "TREM2",
        "genes": ["TREM2", "APOE", "C1QA", "C1QB", "C1QC", "GPNMB", "LGALS3"]}
RNG = np.random.default_rng(15072026)
N_MC = 30


def aggregates_for(anchor, genes, pos_frac):
    """Per-sample cell-equal gene aggregates under a given positive-cell fraction."""
    out_rows = []
    for path in S.SCRNA_PATHS:
        obs, expr, dataset = S.extract_myeloid(path, genes)
        for sample_id, idx in obs.groupby("sample_id", observed=True).groups.items():
            idx = np.asarray(list(idx), dtype=int)
            if len(idx) < 50:
                continue
            sample_expr = expr.iloc[idx].reset_index(drop=True)
            x = sample_expr[anchor].to_numpy(float)
            pos = np.where(x > 0)[0]
            neg = np.where(x <= 0)[0]
            if len(pos) == 0:
                continue
            if pos_frac is not None:
                n_target = max(1, int(round(pos_frac * len(idx))))
                keep_pos = RNG.choice(pos, size=min(n_target, len(pos)), replace=False)
                keep = np.concatenate([keep_pos, neg])
                sub = sample_expr.iloc[keep].reset_index(drop=True)
            else:
                sub = sample_expr
            agg = np.log1p(np.mean(np.expm1(sub.to_numpy(float)), axis=0))
            row = {"dataset": dataset, "sample_id": str(sample_id)}
            row.update(dict(zip(sub.columns, agg)))
            out_rows.append(row)
    return pd.DataFrame(out_rows)


def pooled_rank(frame, anchor, prog_genes):
    scored = []
    for dataset, group in frame.groupby("dataset"):
        group = group.copy()
        present = [g for g in prog_genes if g != anchor and g in group.columns]
        z = F.zscore_columns(group[present].to_numpy(float))
        group["score"] = np.nanmean(z, axis=1)
        scored.append(group)
    scored = pd.concat(scored, ignore_index=True)
    return F.stratified_rank_inference(scored, anchor, "score")


def main():
    out_dir = Path(tempfile.gettempdir()) / "S22_outputs"
    out_dir.mkdir(exist_ok=True)
    report = {"programme": "TREM2_TAM", "reference_siglec15_core_aggregate_r": -0.449}
    natural = aggregates_for(PROG["anchor"], PROG["genes"], None)
    report["natural"] = dict(zip(["r", "ci_low", "ci_high", "perm_p"],
                                 pooled_rank(natural, PROG["anchor"], PROG["genes"])))
    for frac in (0.02, 0.01):
        rs = []
        for _ in range(N_MC):
            frame = aggregates_for(PROG["anchor"], PROG["genes"], frac)
            r, lo, hi, p = pooled_rank(frame, PROG["anchor"], PROG["genes"])
            rs.append(r)
        rs = np.asarray(rs, float)
        report[f"diluted_{frac}"] = {
            "mean_r": float(np.mean(rs)),
            "p05": float(np.percentile(rs, 5)),
            "p95": float(np.percentile(rs, 95)),
            "n_realizations": int(N_MC),
        }
    with open(out_dir / "S22_summary.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)
    print("S22_SUMMARY_JSON_BEGIN")
    print(json.dumps(report, indent=2))
    print("S22_SUMMARY_JSON_END")


if __name__ == "__main__":
    main()
