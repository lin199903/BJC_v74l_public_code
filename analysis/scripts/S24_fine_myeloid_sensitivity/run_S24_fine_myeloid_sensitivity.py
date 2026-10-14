"""Post hoc S24: fine-myeloid-subtype sensitivity for the ACP5/MMP9 core.

The analysis follows PROTOCOL_LOCKED_BEFORE_RESULTS.md in this directory.
"""

from __future__ import annotations

import hashlib
import json
import platform
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy import sparse, stats


ROOT = Path(r".")
OUT = ROOT / "data" / "analysis" / "S24_fine_myeloid_sensitivity_v34"
ANNOTATION = ROOT / "data" / "analysis" / "Task5_myeloid_subtypes" / "myeloid_combined_obs.csv"
H5ADS = [
    ROOT / "data" / "analysis" / "scRNA_GSE155698_sub.h5ad",
    ROOT / "data" / "analysis" / "scRNA_GSE156405.h5ad",
    ROOT / "data" / "analysis" / "scRNA_GSE212966.h5ad",
]
TUMOUR_TISSUES = {"PDAC", "PDAC_tumor", "PDAC_tissue"}
SEED = 10092026
N_BOOT = 10_000


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x, ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return np.full(len(x), np.nan)
    return (x - np.nanmean(x)) / sd


def standardise_matrix(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    means = np.nanmean(x, axis=0)
    sds = np.nanstd(x, axis=0, ddof=0)
    sds[(~np.isfinite(sds)) | (sds == 0)] = 1.0
    return np.nan_to_num((x - means) / sds)


def exposure_coefficient(score: np.ndarray, positive: np.ndarray,
                         subtype: pd.Series | None,
                         covariates: np.ndarray | None = None) -> float:
    positive = np.asarray(positive, dtype=float)
    nuisance = [np.ones(len(score))]
    if subtype is not None:
        counts = subtype.astype(str).value_counts()
        ordered = list(counts.index)
        if len(ordered) > 1:
            dummies = pd.get_dummies(
                pd.Categorical(subtype.astype(str), categories=ordered),
                drop_first=True,
                dtype=float,
            ).to_numpy()
            if dummies.ndim == 1:
                dummies = dummies[:, None]
            nuisance.extend([dummies[:, i] for i in range(dummies.shape[1])])
    if covariates is not None:
        covariates = standardise_matrix(covariates)
        nuisance.extend([covariates[:, i] for i in range(covariates.shape[1])])
    nuisance_design = np.column_stack(nuisance)
    valid = np.isfinite(score) & np.isfinite(positive) & np.all(np.isfinite(nuisance_design), axis=1)
    if valid.sum() <= nuisance_design.shape[1] + 2 or len(np.unique(positive[valid])) < 2:
        return np.nan
    nuisance_valid = nuisance_design[valid]
    positive_valid = positive[valid]
    score_valid = score[valid]
    positive_residual = positive_valid - nuisance_valid @ np.linalg.lstsq(
        nuisance_valid, positive_valid, rcond=None
    )[0]
    score_residual = score_valid - nuisance_valid @ np.linalg.lstsq(
        nuisance_valid, score_valid, rcond=None
    )[0]
    denominator = float(positive_residual @ positive_residual)
    if denominator <= 1e-10:
        return np.nan
    return float((positive_residual @ score_residual) / denominator)


def stratified_bootstrap(frame: pd.DataFrame, column: str) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    strata = [g[column].dropna().to_numpy(float) for _, g in frame.groupby("dataset")]
    boot = np.empty(N_BOOT)
    for i in range(N_BOOT):
        draw = np.concatenate([rng.choice(v, size=len(v), replace=True) for v in strata if len(v)])
        boot[i] = draw.mean()
    return tuple(np.quantile(boot, [0.025, 0.975]))


def summarise(frame: pd.DataFrame, column: str) -> dict:
    values = frame[column].dropna().to_numpy(float)
    lo, hi = stratified_bootstrap(frame.dropna(subset=[column]), column)
    wilcox = stats.wilcoxon(values, alternative="two-sided", zero_method="wilcox")
    return {
        "n_samples": int(len(values)),
        "n_datasets": int(frame.loc[frame[column].notna(), "dataset"].nunique()),
        "mean": float(values.mean()),
        "median": float(np.median(values)),
        "ci_low": float(lo),
        "ci_high": float(hi),
        "wilcoxon_p": float(wilcox.pvalue),
        "n_positive_direction": int((values > 0).sum()),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    labels = pd.read_csv(ANNOTATION, index_col=0, low_memory=False)
    if not labels.index.is_unique:
        raise RuntimeError("Fine-subtype annotation contains duplicate cell barcodes")

    sample_rows: list[dict] = []
    match_rows: list[dict] = []
    subtype_counts: list[dict] = []

    for path in H5ADS:
        obj = ad.read_h5ad(path, backed="r")
        dataset_values = obj.obs["dataset"].astype(str)
        dataset = str(dataset_values.mode().iloc[0])
        mask = (
            obj.obs["tissue"].astype(str).isin(TUMOUR_TISSUES)
            & obj.obs["celltype"].astype(str).eq("Myeloid")
        )
        obs = obj.obs.loc[mask].copy()
        obs.index = obs.index.astype(str)
        lab = labels.loc[
            labels["dataset"].astype(str).eq(dataset)
            & labels["tissue"].astype(str).isin(TUMOUR_TISSUES),
            ["myeloid_subtype"],
        ].copy()
        lab.index = lab.index.astype(str)
        overlap = obs.index.intersection(lab.index)
        match_rows.append({
            "dataset": dataset,
            "n_parent_myeloid": int(len(obs)),
            "n_fine_annotation": int(len(lab)),
            "n_matched": int(len(overlap)),
            "n_unmatched_parent": int(len(obs) - len(overlap)),
            "n_unmatched_annotation": int(len(lab) - len(overlap)),
        })
        if len(overlap) != len(obs) or len(overlap) != len(lab):
            raise RuntimeError(f"Non-complete barcode match for {dataset}")
        obs = obs.loc[overlap].join(lab.loc[overlap])

        gene_idx = [obj.var_names.get_loc(g) for g in ["SIGLEC15", "ACP5", "MMP9"]]
        row_idx = obj.obs_names.get_indexer(overlap)
        x = obj.X[row_idx, :][:, gene_idx]
        if sparse.issparse(x):
            x = x.toarray()
        x = np.asarray(x, dtype=float)
        expr = pd.DataFrame(x, index=overlap, columns=["SIGLEC15", "ACP5", "MMP9"])
        obj.file.close()

        for sample_id, ids in obs.groupby("sample_id", observed=True).groups.items():
            ids = pd.Index(ids)
            if len(ids) < 50:
                continue
            o = obs.loc[ids]
            e = expr.loc[ids]
            positive = e["SIGLEC15"].to_numpy(float) > 0
            if positive.sum() == 0 or (~positive).sum() == 0:
                continue
            score = np.nanmean(
                np.column_stack([zscore(e["ACP5"]), zscore(e["MMP9"])]), axis=1
            )
            covars = np.column_stack([
                np.log1p(pd.to_numeric(o["total_counts"], errors="coerce")),
                np.log1p(pd.to_numeric(o["n_genes_by_counts"], errors="coerce")),
                pd.to_numeric(o["pct_counts_mt"], errors="coerce"),
            ])
            sub = o["myeloid_subtype"].astype(str)
            subtype_overlap = sum(
                int(positive[sub.eq(label).to_numpy()].sum() > 0
                    and (~positive[sub.eq(label).to_numpy()]).sum() > 0)
                for label in sub.unique()
            )
            unadjusted = float(np.nanmean(score[positive]) - np.nanmean(score[~positive]))
            subtype_adjusted = exposure_coefficient(score, positive, sub)
            subtype_depth_adjusted = exposure_coefficient(score, positive, sub, covars)

            trem2 = sub.eq("TREM2_TAM").to_numpy()
            trem2_pos = positive[trem2]
            trem2_score = score[trem2]
            trem2_covars = covars[trem2]
            trem2_diff = np.nan
            trem2_depth = np.nan
            if trem2.sum() > 0 and trem2_pos.sum() > 0 and (~trem2_pos).sum() > 0:
                trem2_diff = float(
                    np.nanmean(trem2_score[trem2_pos])
                    - np.nanmean(trem2_score[~trem2_pos])
                )
                trem2_depth = exposure_coefficient(
                    trem2_score, trem2_pos, None, trem2_covars
                )

            sample_rows.append({
                "dataset": dataset,
                "sample_id": str(sample_id),
                "n_myeloid": int(len(ids)),
                "n_positive": int(positive.sum()),
                "n_negative": int((~positive).sum()),
                "n_subtypes": int(sub.nunique()),
                "n_subtypes_with_exposure_overlap": int(subtype_overlap),
                "unadjusted_difference": unadjusted,
                "fine_subtype_adjusted_coefficient": subtype_adjusted,
                "fine_subtype_depth_adjusted_coefficient": subtype_depth_adjusted,
                "n_trem2_tam": int(trem2.sum()),
                "n_trem2_positive": int(trem2_pos.sum()),
                "n_trem2_negative": int((~trem2_pos).sum()) if trem2.sum() else 0,
                "trem2_restricted_difference": trem2_diff,
                "trem2_depth_adjusted_coefficient": trem2_depth,
            })
            for subtype, n in sub.value_counts().items():
                subtype_counts.append({
                    "dataset": dataset,
                    "sample_id": str(sample_id),
                    "myeloid_subtype": subtype,
                    "n_cells": int(n),
                    "n_positive": int((positive & sub.eq(subtype).to_numpy()).sum()),
                })

    samples = pd.DataFrame(sample_rows)
    matching = pd.DataFrame(match_rows)
    counts = pd.DataFrame(subtype_counts)
    metrics = [
        "unadjusted_difference",
        "fine_subtype_adjusted_coefficient",
        "fine_subtype_depth_adjusted_coefficient",
        "trem2_restricted_difference",
        "trem2_depth_adjusted_coefficient",
    ]
    summary = {metric: summarise(samples, metric) for metric in metrics}
    extreme_index = samples["unadjusted_difference"].abs().idxmax()
    extreme_sample = samples.loc[extreme_index, ["dataset", "sample_id"]].to_dict()
    without_extreme = samples.drop(index=extreme_index)
    posthoc_without_extreme = {
        metric: summarise(without_extreme, metric)
        for metric in [
            "fine_subtype_adjusted_coefficient",
            "fine_subtype_depth_adjusted_coefficient",
        ]
    }
    leave_one_dataset_out_rows = []
    for omitted in sorted(samples["dataset"].unique()):
        retained = samples.loc[~samples["dataset"].eq(omitted)].copy()
        for metric in [
            "fine_subtype_adjusted_coefficient",
            "fine_subtype_depth_adjusted_coefficient",
        ]:
            estimate = summarise(retained, metric)
            leave_one_dataset_out_rows.append({
                "omitted_dataset": omitted,
                "metric": metric,
                **estimate,
            })
    leave_one_dataset_out = pd.DataFrame(leave_one_dataset_out_rows)
    dataset_summary = (
        samples.groupby("dataset")[metrics]
        .agg(["count", "mean", "median"])
        .reset_index()
    )
    dataset_summary.columns = [
        "_".join([str(v) for v in col if str(v)]) if isinstance(col, tuple) else str(col)
        for col in dataset_summary.columns
    ]
    provenance = {
        "protocol": str(OUT / "PROTOCOL_LOCKED_BEFORE_RESULTS.md"),
        "seed": SEED,
        "bootstrap_resamples": N_BOOT,
        "software": {
            "python": platform.python_version(),
            "anndata": version("anndata"),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
        "input_sha256": {
            str(ANNOTATION): sha256(ANNOTATION),
            **{str(path): sha256(path) for path in H5ADS},
        },
        "barcode_matching": matching.to_dict(orient="records"),
        "summary": summary,
        "posthoc_influence_diagnostics": {
            "excluded_extreme_sample": extreme_sample,
            "without_extreme": posthoc_without_extreme,
            "leave_one_dataset_out": leave_one_dataset_out.to_dict(orient="records"),
        },
    }

    samples.to_csv(OUT / "S24_sample_level_effects.csv", index=False)
    counts.to_csv(OUT / "S24_sample_subtype_counts.csv", index=False)
    matching.to_csv(OUT / "S24_barcode_matching.csv", index=False)
    dataset_summary.to_csv(OUT / "S24_dataset_summary.csv", index=False)
    leave_one_dataset_out.to_csv(OUT / "S24_leave_one_dataset_out.csv", index=False)
    (OUT / "S24_summary.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(provenance, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
