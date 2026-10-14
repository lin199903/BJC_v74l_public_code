"""Locked cross-scale analysis of SIGLEC15 and osteoclast/ECM programs in PDAC.

The protocol and gene membership were frozen in LOCKED_PROTOCOL_2026-07-22.md
before association estimates were inspected. Biological samples/cases, never cells
or spots, are the inferential units.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import platform
import sys
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import optimize, stats


HERE = Path(__file__).resolve().parent
SUBMISSION = HERE.parents[1]
WORKSPACE = SUBMISSION.parents[1]
OUT = HERE / "outputs"
FIG = SUBMISSION / "figures"

SEED = 15072026
N_RESAMPLES = 10_000
RNG = np.random.default_rng(SEED)

TARGET = "SIGLEC15"
PROGRAMS = {
    "Derosiers_ACP5_MMP9_core": ["ACP5", "MMP9"],
    "canonical_osteoclast_panel": [
        "ACP5", "CTSK", "MMP9", "ATP6V0D2", "DCSTAMP", "OCSTAMP",
        "ITGB3", "CALCR", "NFATC1",
    ],
    "compact_ECM_degradation_panel": [
        "MMP2", "MMP9", "MMP12", "MMP14", "CTSB", "CTSK", "CTSL",
        "PLAU", "PLAUR",
    ],
}
DESCRIPTIVE = ["IL1B", "IL6"]
POSTHOC_MYELOID_COMPOSITION = [
    "LST1", "AIF1", "FCER1G", "TYROBP", "CTSS", "LILRB1", "CSF1R", "SPI1"
]
ALL_GENES = list(dict.fromkeys(
    [TARGET] + sum(PROGRAMS.values(), []) + DESCRIPTIVE + POSTHOC_MYELOID_COMPOSITION
))
PRIMARY = "Derosiers_ACP5_MMP9_core"
TUMOUR_TISSUES = {"PDAC", "PDAC_tumor", "PDAC_tissue"}

SCRNA_PATHS = [
    SUBMISSION / "biological_validation/data/scrna_pdac/scRNA_GSE212966.h5ad",
    SUBMISSION / "biological_validation/data/scrna_pdac/scRNA_GSE156405.h5ad",
    SUBMISSION / "biological_validation/data/scrna_pdac/scRNA_GSE155698_sub.h5ad",
]
SPATIAL_PATH = WORKSPACE / "data/analysis/Visium_GSE282302.h5ad"
BULK_PATHS = {
    "TCGA_PAAD": WORKSPACE / "data/analysis/TCGA_PAAD_178_log2cpm.csv",
    "GSE57495": WORKSPACE / "data/GEO/GSE57495/GSE57495_expression.csv",
    "GSE85916": WORKSPACE / "data/GEO/GSE85916/GSE85916_expression.csv",
}
SOFT_PATHS = {
    "GSE57495": WORKSPACE / "data/GEO/GSE57495/GSE57495_family.soft.gz",
    "GSE85916": WORKSPACE / "data/GEO/GSE85916/GSE85916_family.soft.gz",
}


def json_default(value):
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def zscore_columns(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    mean = np.nanmean(values, axis=0, keepdims=True)
    sd = np.nanstd(values, axis=0, ddof=0, keepdims=True)
    sd[(~np.isfinite(sd)) | (sd == 0)] = np.nan
    return (values - mean) / sd


def required_gene_count(program: str) -> int:
    if program == PRIMARY:
        return 2
    return max(3, int(np.ceil(0.60 * len(PROGRAMS[program]))))


def score_program(frame: pd.DataFrame, program: str) -> tuple[np.ndarray, list[str]]:
    present = [gene for gene in PROGRAMS[program] if gene in frame.columns]
    if len(present) < required_gene_count(program):
        return np.full(len(frame), np.nan), present
    z = zscore_columns(frame[present].to_numpy(float))
    return np.nanmean(z, axis=1), present


def residualize_depth(score: np.ndarray, obs: pd.DataFrame) -> np.ndarray:
    covariates = []
    for name in ("total_counts", "n_genes_by_counts"):
        if name in obs:
            covariates.append(np.log1p(pd.to_numeric(obs[name], errors="coerce").to_numpy(float)))
    if "pct_counts_mt" in obs:
        covariates.append(pd.to_numeric(obs["pct_counts_mt"], errors="coerce").to_numpy(float))
    if not covariates:
        return score.copy()
    x = np.column_stack(covariates)
    valid = np.isfinite(score) & np.all(np.isfinite(x), axis=1)
    residual = np.full(len(score), np.nan)
    if valid.sum() <= x.shape[1] + 2:
        return residual
    x_valid = zscore_columns(x[valid])
    x_valid = np.column_stack([np.ones(valid.sum()), np.nan_to_num(x_valid)])
    beta = np.linalg.lstsq(x_valid, score[valid], rcond=None)[0]
    residual[valid] = score[valid] - x_valid @ beta
    return residual


def holm_adjust(p_values: list[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.where(np.isfinite(p))[0]
    if not len(valid):
        return out.tolist()
    order = valid[np.argsort(p[valid])]
    running = 0.0
    m = len(order)
    for rank, idx in enumerate(order):
        adjusted = min(1.0, (m - rank) * p[idx])
        running = max(running, adjusted)
        out[idx] = running
    return out.tolist()


def stratified_bootstrap_mean(
    effects: pd.DataFrame, value_col: str, strata_col: str = "dataset"
) -> tuple[float, float]:
    strata = [g[value_col].dropna().to_numpy(float) for _, g in effects.groupby(strata_col)]
    boot = np.empty(N_RESAMPLES)
    for b in range(N_RESAMPLES):
        draw = np.concatenate([RNG.choice(x, size=len(x), replace=True) for x in strata if len(x)])
        boot[b] = np.mean(draw)
    return tuple(np.quantile(boot, [0.025, 0.975]))


def spearman_bootstrap(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, float]:
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[valid], np.asarray(y)[valid]
    rho, p = stats.spearmanr(x, y)
    idx = RNG.integers(0, len(x), size=(N_RESAMPLES, len(x)))
    xr = stats.rankdata(x[idx], axis=1)
    yr = stats.rankdata(y[idx], axis=1)
    xr -= xr.mean(axis=1, keepdims=True)
    yr -= yr.mean(axis=1, keepdims=True)
    denom = np.sqrt(np.sum(xr * xr, axis=1) * np.sum(yr * yr, axis=1))
    boot = np.divide(np.sum(xr * yr, axis=1), denom, out=np.full(N_RESAMPLES, np.nan), where=denom > 0)
    lo, hi = np.nanquantile(boot, [0.025, 0.975])
    return float(rho), float(p), float(lo), float(hi)


def permutation_p(x: np.ndarray, y: np.ndarray, observed: float) -> float:
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[valid], np.asarray(y)[valid]
    xr = stats.rankdata(x).astype(float)
    yr = stats.rankdata(y).astype(float)
    xr -= xr.mean(); yr -= yr.mean()
    denom = np.sqrt(np.sum(xr * xr) * np.sum(yr * yr))
    permutations = np.vstack([RNG.permutation(yr) for _ in range(N_RESAMPLES)])
    r = permutations @ xr / denom
    return float((np.sum(np.abs(r) >= abs(observed)) + 1) / (N_RESAMPLES + 1))


def pooled_within_dataset_rank(frame: pd.DataFrame, x_col: str, y_col: str) -> float:
    xs, ys = [], []
    for _, group in frame.groupby("dataset"):
        if len(group) < 4:
            continue
        x = stats.rankdata(group[x_col].to_numpy(float))
        y = stats.rankdata(group[y_col].to_numpy(float))
        x = (x - x.mean()) / (x.std(ddof=0) or 1)
        y = (y - y.mean()) / (y.std(ddof=0) or 1)
        xs.append(x)
        ys.append(y)
    if not xs:
        return np.nan
    return float(np.corrcoef(np.concatenate(xs), np.concatenate(ys))[0, 1])


def stratified_rank_inference(frame: pd.DataFrame, x_col: str, y_col: str):
    observed = pooled_within_dataset_rank(frame, x_col, y_col)
    boot, perm = [], []
    grouped = [
        (g[x_col].to_numpy(float), g[y_col].to_numpy(float))
        for _, g in frame.groupby("dataset")
    ]
    def combine(parts):
        xs, ys = [], []
        for x, y in parts:
            rx, ry = stats.rankdata(x), stats.rankdata(y)
            rx = (rx - rx.mean()) / (rx.std(ddof=0) or 1)
            ry = (ry - ry.mean()) / (ry.std(ddof=0) or 1)
            xs.append(rx); ys.append(ry)
        return float(np.corrcoef(np.concatenate(xs), np.concatenate(ys))[0, 1])
    for _ in range(N_RESAMPLES):
        draws, shuffled = [], []
        for x, y in grouped:
            idx = RNG.integers(0, len(x), len(x))
            draws.append((x[idx], y[idx]))
            shuffled.append((x, RNG.permutation(y)))
        boot.append(combine(draws))
        perm.append(combine(shuffled))
    boot = np.asarray(boot, float)
    perm = np.asarray(perm, float)
    ci = np.nanquantile(boot, [0.025, 0.975])
    p = (np.sum(np.abs(perm[np.isfinite(perm)]) >= abs(observed)) + 1) / (np.isfinite(perm).sum() + 1)
    return observed, float(ci[0]), float(ci[1]), float(p)


def parse_soft_platform(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    in_table = False
    header = None
    gene_col = None
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                for i, value in enumerate(header):
                    if "symbol" in value.casefold():
                        gene_col = i
                        break
                continue
            if gene_col is None or len(parts) <= gene_col:
                continue
            symbol = parts[gene_col].split("///")[0].strip()
            if symbol and symbol != "---":
                mapping[parts[0]] = symbol
    return mapping


def load_geo_expression(cohort: str, path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, index_col=0)
    mapping = parse_soft_platform(SOFT_PATHS[cohort])
    genes = pd.Series(raw.index.map(mapping), index=raw.index)
    keep = genes.isin(ALL_GENES)
    selected = raw.loc[keep].copy()
    selected.insert(0, "gene", genes.loc[keep].to_numpy())
    return selected.groupby("gene").mean(numeric_only=True)


def load_bulk() -> dict[str, pd.DataFrame]:
    cohorts = {"TCGA_PAAD": pd.read_csv(BULK_PATHS["TCGA_PAAD"], index_col=0)}
    cohorts["GSE57495"] = load_geo_expression("GSE57495", BULK_PATHS["GSE57495"])
    cohorts["GSE85916"] = load_geo_expression("GSE85916", BULK_PATHS["GSE85916"])
    return cohorts


def extract_h5ad(path: Path, spatial: bool = False) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    a = ad.read_h5ad(path, backed="r")
    obs = a.obs.copy()
    if spatial:
        mask = obs["spot_type"].astype(str).eq("Myeloid")
        dataset = "GSE282302"
    else:
        mask = obs["tissue"].astype(str).isin(TUMOUR_TISSUES) & obs["celltype"].astype(str).eq("Myeloid")
        dataset = str(obs.loc[mask, "dataset"].astype(str).iloc[0])
    present = [gene for gene in ALL_GENES if gene in a.var_names]
    subset = a[mask.to_numpy(), present].to_memory()
    matrix = subset.X.toarray() if scipy.sparse.issparse(subset.X) else np.asarray(subset.X)
    expression = pd.DataFrame(matrix, columns=present, index=subset.obs_names)
    obs = subset.obs.copy().reset_index(drop=True)
    expression = expression.reset_index(drop=True)
    a.file.close()
    return obs, expression, dataset


def analyse_scrna():
    sample_effect_rows = []
    aggregate_rows = []
    eligibility_rows = []
    coverage_rows = []
    for path in SCRNA_PATHS:
        obs, expr, dataset = extract_h5ad(path)
        for program, genes in PROGRAMS.items():
            coverage_rows.append({
                "layer": "scRNA", "dataset": dataset, "program": program,
                "n_present": sum(g in expr for g in genes), "n_total": len(genes),
                "genes_present": ";".join(g for g in genes if g in expr),
            })
        for sample_id, idx in obs.groupby("sample_id", observed=True).groups.items():
            idx = np.asarray(list(idx), dtype=int)
            if len(idx) < 50:
                continue
            sample_obs = obs.iloc[idx].reset_index(drop=True)
            sample_expr = expr.iloc[idx].reset_index(drop=True)
            positive = sample_expr[TARGET].to_numpy(float) > 0
            aggregate = np.log1p(np.mean(np.expm1(sample_expr.to_numpy(float)), axis=0))
            aggregate_row = {
                "dataset": dataset, "sample_id": str(sample_id), "n_myeloid": len(idx),
                "n_siglec15_positive": int(positive.sum()),
            }
            aggregate_row.update(dict(zip(sample_expr.columns, aggregate)))
            aggregate_rows.append(aggregate_row)
            eligibility_rows.append({
                "dataset": dataset, "sample_id": str(sample_id), "n_myeloid": len(idx),
                "n_positive": int(positive.sum()), "cell_contrast_estimable": bool(positive.sum() > 0 and (~positive).sum() > 0),
                "aggregate_included": True,
            })
            if positive.sum() == 0 or (~positive).sum() == 0:
                continue
            for program in PROGRAMS:
                score, present = score_program(sample_expr, program)
                if not np.isfinite(score).any():
                    continue
                residual = residualize_depth(score, sample_obs)
                sample_effect_rows.append({
                    "dataset": dataset, "sample_id": str(sample_id), "program": program,
                    "n_myeloid": len(idx), "n_positive": int(positive.sum()),
                    "n_negative": int((~positive).sum()), "n_genes_present": len(present),
                    "difference_positive_minus_negative": float(np.nanmean(score[positive]) - np.nanmean(score[~positive])),
                    "depth_adjusted_difference": float(np.nanmean(residual[positive]) - np.nanmean(residual[~positive])),
                })

    effects = pd.DataFrame(sample_effect_rows)
    aggregates = pd.DataFrame(aggregate_rows)
    summaries = []
    for program, group in effects.groupby("program"):
        values = group["difference_positive_minus_negative"].to_numpy(float)
        ci_lo, ci_hi = stratified_bootstrap_mean(group, "difference_positive_minus_negative")
        test = stats.wilcoxon(values, alternative="two-sided", zero_method="wilcox")
        adj_values = group["depth_adjusted_difference"].to_numpy(float)
        adj_lo, adj_hi = stratified_bootstrap_mean(group, "depth_adjusted_difference")
        adj_test = stats.wilcoxon(adj_values, alternative="two-sided", zero_method="wilcox")
        summaries.append({
            "program": program, "n_samples": len(group), "n_datasets": group["dataset"].nunique(),
            "mean_sample_difference": float(values.mean()), "ci_low": ci_lo, "ci_high": ci_hi,
            "wilcoxon_p": float(test.pvalue),
            "depth_adjusted_mean_difference": float(adj_values.mean()),
            "depth_adjusted_ci_low": adj_lo, "depth_adjusted_ci_high": adj_hi,
            "depth_adjusted_wilcoxon_p": float(adj_test.pvalue),
        })
    cell_summary = pd.DataFrame(summaries)
    secondary_mask = cell_summary["program"].ne(PRIMARY)
    cell_summary.loc[secondary_mask, "holm_p_secondary"] = holm_adjust(
        cell_summary.loc[secondary_mask, "wilcoxon_p"].tolist()
    )

    association_rows = []
    for dataset, group in aggregates.groupby("dataset"):
        for program in PROGRAMS:
            score, present = score_program(group, program)
            if not np.isfinite(score).any():
                continue
            rho, p, lo, hi = spearman_bootstrap(group[TARGET].to_numpy(float), score)
            pp = permutation_p(group[TARGET].to_numpy(float), score, rho)
            association_rows.append({
                "scope": "dataset", "dataset": dataset, "program": program,
                "n_samples": len(group), "n_genes_present": len(present), "rho": rho,
                "ci_low": lo, "ci_high": hi, "asymptotic_p": p, "permutation_p": pp,
            })

    aggregate_scored = []
    for dataset, group in aggregates.groupby("dataset"):
        group = group.copy()
        for program in PROGRAMS:
            group[program], _ = score_program(group, program)
        aggregate_scored.append(group)
    aggregate_scored = pd.concat(aggregate_scored, ignore_index=True)
    for program in PROGRAMS:
        rho, lo, hi, p = stratified_rank_inference(aggregate_scored, TARGET, program)
        association_rows.append({
            "scope": "pooled_within_dataset_ranks", "dataset": "ALL", "program": program,
            "n_samples": len(aggregate_scored), "n_genes_present": len(PROGRAMS[program]),
            "rho": rho, "ci_low": lo, "ci_high": hi, "asymptotic_p": np.nan,
            "permutation_p": p,
        })
    associations = pd.DataFrame(association_rows)
    for scope, idx in associations.groupby("scope").groups.items():
        secondary = [i for i in idx if associations.loc[i, "program"] != PRIMARY]
        associations.loc[secondary, "holm_p_secondary"] = holm_adjust(
            associations.loc[secondary, "permutation_p"].tolist()
        )
    return effects, cell_summary, aggregates, aggregate_scored, associations, pd.DataFrame(eligibility_rows), coverage_rows


def random_effects_meta(group: pd.DataFrame) -> dict:
    y = np.arctanh(np.clip(group["rho"].to_numpy(float), -0.999999, 0.999999))
    n = group["n_samples"].to_numpy(float)
    v = 1.0 / (n - 3.0)
    k = len(y)
    fixed_w = 1 / v
    fixed_mu = np.sum(fixed_w * y) / np.sum(fixed_w)
    q_fixed = np.sum(fixed_w * (y - fixed_mu) ** 2)
    i2 = max(0.0, (q_fixed - (k - 1)) / q_fixed) * 100 if q_fixed > 0 else 0.0

    def objective(tau2):
        w = 1.0 / (v + tau2)
        mu = np.sum(w * y) / np.sum(w)
        return np.sum(np.log(v + tau2)) + np.log(np.sum(w)) + np.sum(w * (y - mu) ** 2)

    result = optimize.minimize_scalar(objective, bounds=(0.0, 4.0), method="bounded")
    tau2 = max(0.0, float(result.x))
    if tau2 < 1e-8:
        tau2 = 0.0
    w = 1.0 / (v + tau2)
    mu = np.sum(w * y) / np.sum(w)
    q_hk_raw = np.sum(w * (y - mu) ** 2) / (k - 1)
    # Modified Hartung--Knapp safeguard: with very few homogeneous studies,
    # unmodified HK can make the interval paradoxically narrower than the
    # conventional random-effects interval. The scale factor is bounded below by 1.
    q_hk = max(1.0, float(q_hk_raw))
    se_hk = np.sqrt(q_hk / np.sum(w))
    critical = stats.t.ppf(0.975, df=k - 1)
    lo_z, hi_z = mu - critical * se_hk, mu + critical * se_hk
    pred_se = np.sqrt(tau2 + se_hk**2)
    pred_lo, pred_hi = mu - critical * pred_se, mu + critical * pred_se
    p = 2 * stats.t.sf(abs(mu / se_hk), df=k - 1) if se_hk > 0 else np.nan
    return {
        "program": group["program"].iloc[0], "n_cohorts": k,
        "cohorts": ";".join(group["cohort"]), "r_random": float(np.tanh(mu)),
        "ci_low": float(np.tanh(lo_z)), "ci_high": float(np.tanh(hi_z)),
        "p_modified_hartung_knapp": float(p), "hk_scale_raw": float(q_hk_raw),
        "hk_scale_used": float(q_hk), "tau2_fisher_z": tau2,
        "I2_pct": float(i2), "prediction_low": float(np.tanh(pred_lo)),
        "prediction_high": float(np.tanh(pred_hi)),
        "ci_wholly_inside_plus_minus_0_20": bool(np.tanh(lo_z) > -0.2 and np.tanh(hi_z) < 0.2),
    }


def partial_spearman(x: np.ndarray, y: np.ndarray, covariate: np.ndarray) -> float:
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(covariate)
    x, y, covariate = np.asarray(x)[valid], np.asarray(y)[valid], np.asarray(covariate)[valid]
    rx, ry, rc = stats.rankdata(x), stats.rankdata(y), stats.rankdata(covariate)
    design = np.column_stack([np.ones(len(rc)), rc])
    ex = rx - design @ np.linalg.lstsq(design, rx, rcond=None)[0]
    ey = ry - design @ np.linalg.lstsq(design, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def partial_spearman_bootstrap(x: np.ndarray, y: np.ndarray, covariate: np.ndarray):
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(covariate)
    x, y, covariate = np.asarray(x)[valid], np.asarray(y)[valid], np.asarray(covariate)[valid]
    observed = partial_spearman(x, y, covariate)
    boot = np.empty(N_RESAMPLES)
    for b in range(N_RESAMPLES):
        idx = RNG.integers(0, len(x), len(x))
        boot[b] = partial_spearman(x[idx], y[idx], covariate[idx])
    lo, hi = np.nanquantile(boot, [0.025, 0.975])
    df = len(x) - 3
    t_value = observed * np.sqrt(df / max(1e-12, 1 - observed**2))
    p = 2 * stats.t.sf(abs(t_value), df=df)
    return observed, float(lo), float(hi), float(p)


def analyse_bulk():
    cohorts = load_bulk()
    rows, coverage = [], []
    scored_for_plot = {}
    for cohort, genes_by_samples in cohorts.items():
        frame = genes_by_samples.T
        scored_for_plot[cohort] = frame
        for program, genes in PROGRAMS.items():
            present = [g for g in genes if g in frame]
            coverage.append({
                "layer": "bulk", "dataset": cohort, "program": program,
                "n_present": len(present), "n_total": len(genes),
                "genes_present": ";".join(present),
            })
            score, present = score_program(frame, program)
            if not np.isfinite(score).any() or TARGET not in frame:
                continue
            x = frame[TARGET].to_numpy(float)
            rho, p, lo, hi = spearman_bootstrap(x, score)
            pp = permutation_p(x, score, rho)
            rows.append({
                "cohort": cohort, "program": program, "n_samples": len(frame),
                "n_genes_present": len(present), "n_genes_total": len(genes),
                "rho": rho, "ci_low": lo, "ci_high": hi,
                "asymptotic_p": p, "permutation_p": pp,
                "ci_wholly_inside_plus_minus_0_20": bool(lo > -0.2 and hi < 0.2),
            })
    results = pd.DataFrame(rows)
    for cohort, idx in results.groupby("cohort").groups.items():
        secondary = [i for i in idx if results.loc[i, "program"] != PRIMARY]
        results.loc[secondary, "holm_p_secondary"] = holm_adjust(
            results.loc[secondary, "permutation_p"].tolist()
        )
    metas = pd.DataFrame([random_effects_meta(g) for _, g in results.groupby("program")])
    secondary = metas["program"].ne(PRIMARY)
    metas.loc[secondary, "holm_p_secondary"] = holm_adjust(
        metas.loc[secondary, "p_modified_hartung_knapp"].tolist()
    )
    posthoc_rows = []
    for cohort, frame in scored_for_plot.items():
        core_score, core_present = score_program(frame, PRIMARY)
        myeloid_present = [g for g in POSTHOC_MYELOID_COMPOSITION if g in frame]
        myeloid_score = np.nanmean(zscore_columns(frame[myeloid_present].to_numpy(float)), axis=1)
        rho, lo, hi, p = partial_spearman_bootstrap(
            frame[TARGET].to_numpy(float), core_score, myeloid_score
        )
        posthoc_rows.append({
            "cohort": cohort, "program": PRIMARY, "n_samples": len(frame),
            "n_core_genes": len(core_present), "n_myeloid_markers": len(myeloid_present),
            "myeloid_markers_present": ";".join(myeloid_present),
            "partial_rho": rho, "ci_low": lo, "ci_high": hi, "p_value": p,
        })
    posthoc = pd.DataFrame(posthoc_rows)
    meta_input = posthoc.rename(columns={"partial_rho": "rho"})
    posthoc_meta = pd.DataFrame([random_effects_meta(meta_input)])
    posthoc_meta["analysis"] = "posthoc_adjusted_for_pan_myeloid_score"
    return results, metas, posthoc, posthoc_meta, coverage, scored_for_plot


def analyse_spatial():
    obs, expr, dataset = extract_h5ad(SPATIAL_PATH, spatial=True)
    coverage = []
    eligibility_rows = []
    for program, genes in PROGRAMS.items():
        coverage.append({
            "layer": "spatial", "dataset": dataset, "program": program,
            "n_present": sum(g in expr for g in genes), "n_total": len(genes),
            "genes_present": ";".join(g for g in genes if g in expr),
        })
    capture_rows = []
    for sample_id, idx in obs.groupby("sample_id", observed=True).groups.items():
        idx = np.asarray(list(idx), dtype=int)
        capture_obs = obs.iloc[idx].reset_index(drop=True)
        capture_expr = expr.iloc[idx].reset_index(drop=True)
        positive = capture_expr[TARGET].to_numpy(float) > 0
        eligibility_rows.append({
            "dataset": dataset, "case": str(capture_obs["case"].iloc[0]),
            "sample_id": str(sample_id), "n_myeloid_spots": len(idx),
            "n_positive": int(positive.sum()), "n_negative": int((~positive).sum()),
            "eligible": bool(len(idx) >= 10 and positive.sum() >= 2 and (~positive).sum() >= 2),
        })
        if len(idx) < 10 or positive.sum() < 2 or (~positive).sum() < 2:
            continue
        case = str(capture_obs["case"].iloc[0])
        for program in PROGRAMS:
            score, present = score_program(capture_expr, program)
            if not np.isfinite(score).any():
                continue
            residual = residualize_depth(score, capture_obs)
            capture_rows.append({
                "dataset": dataset, "case": case, "sample_id": str(sample_id),
                "program": program, "n_myeloid_spots": len(idx),
                "n_positive": int(positive.sum()), "n_negative": int((~positive).sum()),
                "n_genes_present": len(present),
                "difference_positive_minus_negative": float(np.nanmean(score[positive]) - np.nanmean(score[~positive])),
                "depth_adjusted_difference": float(np.nanmean(residual[positive]) - np.nanmean(residual[~positive])),
            })
    captures = pd.DataFrame(capture_rows)
    cases = captures.groupby(["case", "program"], as_index=False).agg(
        n_captures=("sample_id", "nunique"),
        mean_capture_difference=("difference_positive_minus_negative", "mean"),
        depth_adjusted_mean_capture_difference=("depth_adjusted_difference", "mean"),
    )
    summaries = []
    for program, group in cases.groupby("program"):
        values = group["mean_capture_difference"].to_numpy(float)
        adjusted = group["depth_adjusted_mean_capture_difference"].to_numpy(float)
        crit = stats.t.ppf(0.975, df=len(values) - 1)
        mean = values.mean()
        se = values.std(ddof=1) / np.sqrt(len(values))
        adj_mean = adjusted.mean()
        adj_se = adjusted.std(ddof=1) / np.sqrt(len(adjusted))
        summaries.append({
            "program": program, "n_cases": len(values), "n_eligible_captures": int(group["n_captures"].sum()),
            "mean_case_difference": float(mean), "ci_low": float(mean - crit * se),
            "ci_high": float(mean + crit * se),
            "t_test_p": float(stats.ttest_1samp(values, 0).pvalue),
            "depth_adjusted_mean_case_difference": float(adj_mean),
            "depth_adjusted_ci_low": float(adj_mean - crit * adj_se),
            "depth_adjusted_ci_high": float(adj_mean + crit * adj_se),
            "depth_adjusted_t_test_p": float(stats.ttest_1samp(adjusted, 0).pvalue),
        })
    summary = pd.DataFrame(summaries)
    secondary = summary["program"].ne(PRIMARY)
    summary.loc[secondary, "holm_p_secondary"] = holm_adjust(
        summary.loc[secondary, "t_test_p"].tolist()
    )
    return captures, cases, summary, pd.DataFrame(eligibility_rows), coverage


def make_figure(
    cell_effects: pd.DataFrame,
    cell_summary: pd.DataFrame,
    aggregates_scored: pd.DataFrame,
    aggregate_assoc: pd.DataFrame,
    bulk: pd.DataFrame,
    meta: pd.DataFrame,
    bulk_posthoc: pd.DataFrame,
    posthoc_meta: pd.DataFrame,
    spatial_cases: pd.DataFrame,
    spatial_summary: pd.DataFrame,
    spatial_eligibility: pd.DataFrame,
):
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 8.5, "axes.linewidth": 0.8,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    palette = {"GSE155698": "#0072B2", "GSE156405": "#D55E00", "GSE212966": "#009E73"}
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.4), constrained_layout=True)

    ax = axes[0, 0]
    data = cell_effects[cell_effects["program"] == PRIMARY].copy()
    data = data.sort_values(["dataset", "difference_positive_minus_negative"])
    for i, row in enumerate(data.itertuples()):
        ax.scatter(row.difference_positive_minus_negative, i, s=23, color=palette[row.dataset], edgecolor="white", linewidth=0.4)
    summ = cell_summary[cell_summary["program"] == PRIMARY].iloc[0]
    ax.axvline(0, color="#777777", lw=0.8)
    ax.errorbar(summ.mean_sample_difference, len(data) + 0.7,
                xerr=[[summ.mean_sample_difference - summ.ci_low], [summ.ci_high - summ.mean_sample_difference]],
                fmt="D", color="black", capsize=3, markersize=4)
    ax.set_yticks(list(range(len(data))) + [len(data) + 0.7])
    ax.set_yticklabels(data["sample_id"].tolist() + ["Patient-level mean"])
    ax.set_xlabel("ACP5/MMP9 score difference\nSIGLEC15+ minus SIGLEC15− myeloid cells")
    ax.set_title("a  Tumour-myeloid cells", loc="left", fontweight="bold")

    ax = axes[0, 1]
    for dataset, group in aggregates_scored.groupby("dataset"):
        ax.scatter(group[TARGET], group[PRIMARY], s=32, alpha=0.85, label=dataset, color=palette[dataset], edgecolor="white", linewidth=0.4)
    pooled = aggregate_assoc[(aggregate_assoc["scope"] == "pooled_within_dataset_ranks") & (aggregate_assoc["program"] == PRIMARY)].iloc[0]
    ax.set_xlabel("Myeloid aggregate SIGLEC15")
    ax.set_ylabel("Myeloid aggregate ACP5/MMP9 score")
    ax.set_title(f"b  Sample-level myeloid aggregate\npooled within-series r={pooled.rho:.2f} ({pooled.ci_low:.2f}, {pooled.ci_high:.2f})", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7, loc="best")

    ax = axes[1, 0]
    b = bulk[bulk["program"] == PRIMARY].copy()
    labels = b["cohort"].tolist() + ["REML–mHK meta"]
    values = b["rho"].tolist() + [float(meta.loc[meta["program"] == PRIMARY, "r_random"].iloc[0])]
    lows = b["ci_low"].tolist() + [float(meta.loc[meta["program"] == PRIMARY, "ci_low"].iloc[0])]
    highs = b["ci_high"].tolist() + [float(meta.loc[meta["program"] == PRIMARY, "ci_high"].iloc[0])]
    y = np.arange(len(labels))[::-1]
    ax.axvspan(-0.2, 0.2, color="#EEEEEE", zorder=0)
    ax.axvline(0, color="#777777", lw=0.8)
    for yi, value, lo, hi, label in zip(y, values, lows, highs, labels):
        ax.errorbar(value, yi, xerr=[[value - lo], [hi - value]], fmt="D" if "meta" in label else "o", color="black", capsize=3, markersize=4)
    adjusted_values = bulk_posthoc["partial_rho"].tolist() + [float(posthoc_meta["r_random"].iloc[0])]
    adjusted_lows = bulk_posthoc["ci_low"].tolist() + [float(posthoc_meta["ci_low"].iloc[0])]
    adjusted_highs = bulk_posthoc["ci_high"].tolist() + [float(posthoc_meta["ci_high"].iloc[0])]
    for yi, value, lo, hi in zip(y - 0.16, adjusted_values, adjusted_lows, adjusted_highs):
        ax.errorbar(value, yi, xerr=[[value - lo], [hi - value]], fmt="s", color="#0072B2", capsize=2, markersize=3)
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.set_xlim(-0.6, 0.6)
    ax.set_xlabel("Spearman rho (95% CI)")
    ax.set_title("c  Bulk tumours", loc="left", fontweight="bold")

    ax = axes[1, 1]
    s = spatial_cases[spatial_cases["program"] == PRIMARY].copy().set_index("case")
    all_cases = sorted(spatial_eligibility["case"].astype(str).unique())
    ypos = {case: len(all_cases) - 1 - i for i, case in enumerate(all_cases)}
    for case, row in s.iterrows():
        ax.scatter(row["mean_capture_difference"], ypos[case], s=35, color="#CC79A7", edgecolor="white", linewidth=0.4)
    for case in all_cases:
        if case not in s.index:
            e = spatial_eligibility[spatial_eligibility["case"].astype(str).eq(case)]
            ax.text(-3.25, ypos[case], f"not estimable ({int(e.n_positive.sum())}/{int(e.n_myeloid_spots.sum())} positive)",
                    ha="left", va="center", fontsize=7, color="#666666")
    summ = spatial_summary[spatial_summary["program"] == PRIMARY].iloc[0]
    ax.axvline(0, color="#777777", lw=0.8)
    mean_y = len(all_cases) + 0.5
    ax.errorbar(summ.mean_case_difference, mean_y,
                xerr=[[summ.mean_case_difference - summ.ci_low], [summ.ci_high - summ.mean_case_difference]],
                fmt="D", color="black", capsize=3, markersize=4)
    ax.set_yticks([ypos[c] for c in all_cases] + [mean_y])
    ax.set_yticklabels(all_cases + ["Case-level mean"])
    ax.set_xlim(-3.6, 4.7)
    ax.set_xlabel("Within-capture ACP5/MMP9 difference\n(case-averaged; SIGLEC15+ minus negative)")
    ax.set_title("d  Spatial myeloid-enriched spots (2/3 cases estimable)", loc="left", fontweight="bold")

    for ax in axes.flat:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    pdf = FIG / "Figure_7_cross_scale_osteoclast_ECM.pdf"
    png = FIG / "Figure_7_cross_scale_osteoclast_ECM.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return pdf, png


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    cell_effects, cell_summary, aggregates, aggregate_scored, aggregate_assoc, scrna_eligibility, coverage_sc = analyse_scrna()
    bulk, meta, bulk_posthoc, posthoc_meta, coverage_bulk, _ = analyse_bulk()
    spatial_captures, spatial_cases, spatial_summary, spatial_eligibility, coverage_sp = analyse_spatial()

    tables = {
        "cell_level_sample_effects.csv": cell_effects,
        "cell_level_summary.csv": cell_summary,
        "myeloid_sample_aggregates.csv": aggregates,
        "myeloid_sample_aggregates_scored.csv": aggregate_scored,
        "myeloid_aggregate_associations.csv": aggregate_assoc,
        "scrna_sample_eligibility.csv": scrna_eligibility,
        "bulk_associations.csv": bulk,
        "bulk_random_effects_meta.csv": meta,
        "bulk_posthoc_myeloid_adjusted.csv": bulk_posthoc,
        "bulk_posthoc_myeloid_adjusted_meta.csv": posthoc_meta,
        "spatial_capture_effects.csv": spatial_captures,
        "spatial_case_effects.csv": spatial_cases,
        "spatial_summary.csv": spatial_summary,
        "spatial_capture_eligibility.csv": spatial_eligibility,
        "gene_coverage.csv": pd.DataFrame(coverage_sc + coverage_bulk + coverage_sp),
    }
    for name, table in tables.items():
        table.to_csv(OUT / name, index=False)

    pdf, png = make_figure(
        cell_effects, cell_summary, aggregate_scored, aggregate_assoc,
        bulk, meta, bulk_posthoc, posthoc_meta, spatial_cases, spatial_summary, spatial_eligibility,
    )

    inputs = SCRNA_PATHS + [SPATIAL_PATH] + list(BULK_PATHS.values()) + list(SOFT_PATHS.values())
    manifest_path = OUT / "input_manifest.csv"
    old_manifest = pd.read_csv(manifest_path) if manifest_path.exists() else pd.DataFrame()
    old_lookup = {
        (str(row.path), int(row.size_bytes)): str(row.sha256)
        for row in old_manifest.itertuples()
    } if not old_manifest.empty else {}
    manifest_rows = []
    for path in inputs:
        key = (str(path), path.stat().st_size)
        manifest_rows.append({"path": str(path), "size_bytes": key[1], "sha256": old_lookup.get(key) or sha256_file(path)})
    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(OUT / "input_manifest.csv", index=False)

    primary_cell = cell_summary[cell_summary["program"] == PRIMARY].iloc[0].to_dict()
    primary_aggregate = aggregate_assoc[(aggregate_assoc["scope"] == "pooled_within_dataset_ranks") & (aggregate_assoc["program"] == PRIMARY)].iloc[0].to_dict()
    primary_bulk = meta[meta["program"] == PRIMARY].iloc[0].to_dict()
    primary_spatial = spatial_summary[spatial_summary["program"] == PRIMARY].iloc[0].to_dict()
    summary = {
        "analysis": "F: locked cross-scale SIGLEC15 osteoclast/ECM analysis",
        "protocol": str(HERE / "LOCKED_PROTOCOL_2026-07-22.md"),
        "primary_program": PRIMARY,
        "primary_results": {
            "cell_level_patient_unit": primary_cell,
            "myeloid_aggregate": primary_aggregate,
            "bulk_random_effects_modified_hk": primary_bulk,
            "spatial_case_unit_exploratory": primary_spatial,
            "bulk_posthoc_pan_myeloid_adjusted": posthoc_meta.iloc[0].to_dict(),
        },
        "interpretation_rules": {
            "cell_state_recovered": bool(primary_cell["ci_low"] > 0),
            "bulk_equivalence_claimed": False,
            "spatial_is_exploratory_n_cases_3": True,
            "protein_or_treatment_response_inferred": False,
        },
        "figure_pdf": str(pdf), "figure_png": str(png),
        "software": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
            "scipy": scipy.__version__, "anndata": ad.__version__,
            "matplotlib": mpl.__version__,
        },
    }
    (OUT / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=json_default), encoding="utf-8"
    )
    print(json.dumps(summary["primary_results"], indent=2, default=json_default))


if __name__ == "__main__":
    main()
