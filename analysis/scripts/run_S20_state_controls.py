"""S20 state-level positive controls (post hoc demonstration).

Runs the TREM2-, SPP1- and LRRC15-centred state signatures through the same
cell / aggregate / bulk layers as the locked SIGLEC15 analysis, reusing the
locked functions (zscore, bootstrap, permutation, REML + modified
Hartung--Knapp) from analysis_F_cross_scale.py.

Outputs are written to a sandbox-writable TEMP directory and the summary JSON
is printed to stdout. No workspace writes.
"""
import json
import os
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy import stats

ANALYSIS = Path(r".\biological_validation\analysis_F_cross_scale")
sys.path.insert(0, str(ANALYSIS))
import analysis_F_cross_scale as F  # noqa: E402

SCRNA_PATHS = [
    Path(r".\biological_validation\data\scrna_pdac\scRNA_GSE212966.h5ad"),
    Path(r".\biological_validation\data\scrna_pdac\scRNA_GSE156405.h5ad"),
    Path(r".\biological_validation\data\scrna_pdac\scRNA_GSE155698_sub.h5ad"),
]
BULK_PATHS = {
    "TCGA_PAAD": r".\data\analysis\TCGA_PAAD_178_log2cpm.csv",
    "GSE57495": r".\data\GEO\GSE57495\GSE57495_expression.csv",
    "GSE85916": r".\data\GEO\GSE85916\GSE85916_expression.csv",
}
SOFT_PATHS = {
    "GSE57495": r".\data\GEO\GSE57495\GSE57495_family.soft.gz",
    "GSE85916": r".\data\GEO\GSE85916\GSE85916_family.soft.gz",
}

S20 = {
    "TREM2_TAM": {
        "anchor": "TREM2",
        "genes": ["TREM2", "APOE", "C1QA", "C1QB", "C1QC", "GPNMB", "LGALS3"],
        "cell": True,
        "aggregate": True,
    },
    "SPP1_TAM": {
        "anchor": "SPP1",
        "genes": ["SPP1", "FN1", "VEGFA", "ADM", "EGLN3", "NDRG1", "HIF1A"],
        "cell": True,
        "aggregate": True,
    },
    "LRRC15_CAF": {
        "anchor": "LRRC15",
        "genes": ["LRRC15", "MMP11", "COL10A1", "COL11A1", "TGFB3", "POSTN"],
        "cell": False,
        "aggregate": False,
    },
}


def s20_score(expr_frame: pd.DataFrame, prog: dict) -> tuple[np.ndarray, list[str]]:
    """Programme score with the anchor gene excluded, frozen-style coverage rule."""
    genes = [g for g in prog["genes"] if g != prog["anchor"]]
    present = [g for g in genes if g in expr_frame.columns]
    needed = max(3, int(np.ceil(0.60 * len(genes))))
    if len(present) < needed:
        return np.full(len(expr_frame), np.nan), present
    z = F.zscore_columns(expr_frame[present].to_numpy(float))
    return np.nanmean(z, axis=1), present


def extract_myeloid(path: Path, genes: list[str]):
    a = ad.read_h5ad(path, backed="r")
    obs = a.obs
    mask = obs["tissue"].astype(str).isin(F.TUMOUR_TISSUES) & obs["celltype"].astype(str).eq("Myeloid")
    present = [g for g in genes if g in a.var_names]
    subset = a[mask.to_numpy(), present].to_memory()
    matrix = subset.X.toarray() if scipy.sparse.issparse(subset.X) else np.asarray(subset.X)
    expression = pd.DataFrame(matrix, columns=present, index=subset.obs_names)
    obs_out = subset.obs.copy().reset_index(drop=True)
    expression = expression.reset_index(drop=True)
    dataset = str(subset.obs["dataset"].astype(str).iloc[0])
    a.file.close()
    return obs_out, expression, dataset


def cell_layer(name: str, prog: dict, out_dir: Path):
    genes = sorted(set(prog["genes"]))
    rows, coverage = [], []
    for path in SCRNA_PATHS:
        obs, expr, dataset = extract_myeloid(path, genes)
        coverage.append({"dataset": dataset, "program": name,
                         "genes_present": ";".join(g for g in prog["genes"] if g in expr)})
        for sample_id, idx in obs.groupby("sample_id", observed=True).groups.items():
            idx = np.asarray(list(idx), dtype=int)
            if len(idx) < 50:
                continue
            sample_expr = expr.iloc[idx].reset_index(drop=True)
            if prog["anchor"] not in sample_expr:
                continue
            positive = sample_expr[prog["anchor"]].to_numpy(float) > 0
            if positive.sum() == 0 or (~positive).sum() == 0:
                continue
            score, present = s20_score(sample_expr, prog)
            if not np.isfinite(score).any():
                continue
            rows.append({
                "dataset": dataset, "sample_id": str(sample_id),
                "n_myeloid": len(idx), "n_positive": int(positive.sum()),
                "n_negative": int((~positive).sum()),
                "n_genes_present": len(present),
                "difference_positive_minus_negative": float(np.nanmean(score[positive]) - np.nanmean(score[~positive])),
            })
    effects = pd.DataFrame(rows)
    if effects.empty:
        return {"program": name, "n_samples": 0, "error": "no estimable samples"}, None, coverage
    values = effects["difference_positive_minus_negative"].to_numpy(float)
    ci = F.stratified_bootstrap_mean(effects, "difference_positive_minus_negative")
    wilc = stats.wilcoxon(values, alternative="two-sided", zero_method="wilcox")
    summary = {
        "program": name, "n_samples": len(effects), "n_datasets": effects["dataset"].nunique(),
        "n_positive_cells": int(effects["n_positive"].sum()),
        "mean_sample_difference": float(values.mean()), "ci_low": float(ci[0]), "ci_high": float(ci[1]),
        "wilcoxon_p": float(wilc.pvalue),
    }
    effects.to_csv(out_dir / f"cell_level_{name}.csv", index=False)
    pd.DataFrame(coverage).to_csv(out_dir / f"coverage_{name}.csv", index=False)
    return summary, effects, coverage


def aggregate_layer(name: str, prog: dict, out_dir: Path):
    genes = sorted(set(prog["genes"]))
    agg_rows = []
    for path in SCRNA_PATHS:
        obs, expr, dataset = extract_myeloid(path, genes)
        for sample_id, idx in obs.groupby("sample_id", observed=True).groups.items():
            idx = np.asarray(list(idx), dtype=int)
            if len(idx) < 50:
                continue
            sample_expr = expr.iloc[idx].reset_index(drop=True)
            aggregate = np.log1p(np.mean(np.expm1(sample_expr.to_numpy(float)), axis=0))
            row = {"dataset": dataset, "sample_id": str(sample_id), "n_myeloid": len(idx)}
            row.update(dict(zip(sample_expr.columns, aggregate)))
            agg_rows.append(row)
    aggregates = pd.DataFrame(agg_rows)
    agg_scored = []
    for dataset, group in aggregates.groupby("dataset"):
        group = group.copy()
        score, _ = s20_score(group, prog)
        group[f"{name}_score"] = score
        agg_scored.append(group)
    agg_scored = pd.concat(agg_scored, ignore_index=True)
    rho, lo, hi, p = F.stratified_rank_inference(agg_scored, prog["anchor"], f"{name}_score")
    per_dataset = []
    for dataset, group in agg_scored.groupby("dataset"):
        x = group[prog["anchor"]].to_numpy(float)
        y = group[f"{name}_score"].to_numpy(float)
        valid = np.isfinite(x) & np.isfinite(y)
        r_d = float(stats.spearmanr(x[valid], y[valid])[0])
        per_dataset.append({"dataset": dataset, "n_samples": int(valid.sum()), "rho": r_d})
    summary = {
        "program": name, "n_samples": int(agg_scored[f"{name}_score"].notna().sum()),
        "pooled_rho": float(rho), "ci_low": float(lo), "ci_high": float(hi),
        "permutation_p": float(p), "per_dataset": per_dataset,
    }
    agg_scored.to_csv(out_dir / f"aggregate_{name}.csv", index=False)
    return summary


def load_geo_s20(cohort: str, genes: list[str]) -> pd.DataFrame:
    raw = pd.read_csv(BULK_PATHS[cohort], index_col=0)
    mapping = F.parse_soft_platform(Path(SOFT_PATHS[cohort]))
    genes_s = pd.Series(raw.index.map(mapping), index=raw.index)
    keep = genes_s.isin(genes)
    selected = raw.loc[keep].copy()
    selected.insert(0, "gene", genes_s.loc[keep].to_numpy())
    return selected.groupby("gene").mean(numeric_only=True)


def load_bulk_s20(genes: list[str]) -> dict[str, pd.DataFrame]:
    cohorts = {"TCGA_PAAD": pd.read_csv(BULK_PATHS["TCGA_PAAD"], index_col=0)}
    cohorts["GSE57495"] = load_geo_s20("GSE57495", genes)
    cohorts["GSE85916"] = load_geo_s20("GSE85916", genes)
    return cohorts


def bulk_layer(progs: dict, out_dir: Path):
    genes = sorted({g for p in progs.values() for g in p["genes"]})
    cohorts = load_bulk_s20(genes)
    results_rows, coverage = [], []
    for cohort, df in cohorts.items():
        frame = df.T
        for name, prog in progs.items():
            present = [g for g in prog["genes"] if g != prog["anchor"] and g in frame]
            coverage.append({"cohort": cohort, "program": name,
                             "genes_present": ";".join(present), "n_present": len(present)})
            if prog["anchor"] not in frame:
                continue
            score, pres = s20_score(frame, prog)
            if not np.isfinite(score).any():
                continue
            x = frame[prog["anchor"]].to_numpy(float)
            rho, p, lo, hi = F.spearman_bootstrap(x, score)
            pp = F.permutation_p(x, score, rho)
            results_rows.append({
                "cohort": cohort, "program": name, "n_samples": len(frame),
                "n_genes_present": len(pres), "rho": rho,
                "ci_low": lo, "ci_high": hi, "asymptotic_p": p, "permutation_p": pp,
            })
    results = pd.DataFrame(results_rows)
    metas = {}
    for name, group in results.groupby("program"):
        meta = F.random_effects_meta(group)
        cohort_rhos = group["rho"].to_numpy(float)
        direction_consistent = bool(np.all(cohort_rhos > 0) or np.all(cohort_rhos < 0))
        if abs(meta["r_random"]) >= 0.20 and direction_consistent:
            verdict = "Supportive (direction-consistent)"
        elif meta["ci_low"] > -0.36 and meta["ci_high"] < 0.36:
            verdict = "Non-supportive (failed qualification; excludes |r| >= 0.36)"
        else:
            verdict = "Uninformative"
        metas[name] = {
            "pooled_r": meta["r_random"], "ci_low": meta["ci_low"], "ci_high": meta["ci_high"],
            "p_modified_hartung_knapp": meta["p_modified_hartung_knapp"],
            "I2_pct": meta["I2_pct"], "n_cohorts": meta["n_cohorts"],
            "cohorts": meta["cohorts"], "verdict": verdict,
            "direction_consistent": direction_consistent,
            "cohort_rhos": {r["cohort"]: float(r["rho"]) for _, r in group.iterrows()},
            "cohort_perm_p": {r["cohort"]: float(r["permutation_p"]) for _, r in group.iterrows()},
        }
    results.to_csv(out_dir / "bulk_S20.csv", index=False)
    pd.DataFrame(coverage).to_csv(out_dir / "bulk_coverage_S20.csv", index=False)
    return metas


def tcga_median_check(out_dir: Path):
    frame = pd.read_csv(BULK_PATHS["TCGA_PAAD"], index_col=0).T
    out = {}
    for gene in ("SIGLEC15", "CD274"):
        if gene not in frame:
            continue
        v = frame[gene].to_numpy(float)
        out[gene] = {
            "n": int(len(v)),
            "median_log2cpm": float(np.median(v)),
            "median_cpm_if_2^x": float(np.median(2.0 ** v)),
            "median_cpm_if_2^x_minus_1": float(np.median(2.0 ** v - 1.0)),
        }
    with open(out_dir / "tcga_median_check.json", "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    return out


def main():
    out_dir = Path(tempfile.gettempdir()) / "S20_outputs"
    out_dir.mkdir(exist_ok=True)
    report = {}
    report["cell"] = {}
    report["aggregate"] = {}
    for name, prog in S20.items():
        if prog["cell"]:
            print(f"[S20] cell layer: {name}", flush=True)
            report["cell"][name], _, _ = cell_layer(name, prog, out_dir)
        if prog["aggregate"]:
            print(f"[S20] aggregate layer: {name}", flush=True)
            report["aggregate"][name] = aggregate_layer(name, prog, out_dir)
    print("[S20] bulk layer", flush=True)
    report["bulk"] = bulk_layer(S20, out_dir)
    print("[S20] TCGA median check", flush=True)
    report["tcga_median_check"] = tcga_median_check(out_dir)
    with open(out_dir / "S20_summary.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)
    print("S20_SUMMARY_JSON_BEGIN")
    print(json.dumps(report, indent=2))
    print("S20_SUMMARY_JSON_END")
    print("out_dir:", out_dir)


if __name__ == "__main__":
    main()
