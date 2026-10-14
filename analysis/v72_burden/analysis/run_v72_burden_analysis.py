"""Frozen v72 tumour-level burden and sparsity analysis.

The analysis implements V72_BURDEN_ESTIMAND_CONTRACT.md. Biological samples,
not cells, are inferential units. The v71 archive is read only and is used only
for numerical reproduction checks.
"""

from __future__ import annotations
import os

import hashlib
import json
import platform
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy import sparse, stats


ROOT = Path(os.environ.get("SIGLEC15_PROJECT_ROOT", r"E:\Tex\HMlong"))
BRANCH = ROOT / "BJC_v72_burden_estimand"
OUT = BRANCH / "results"
CONTRACT = BRANCH / "V72_BURDEN_ESTIMAND_CONTRACT.md"
ARCHIVE = (
    ROOT
    / "npjPrecOncol/FINAL_SUBMISSION_v6/BMC_MEDICINE_SUBMISSION_v47"
    / "_archives_translate/Additional_file_2_Analysis_archive_v39/results"
)
PATHS = [
    ROOT / "data/analysis/scRNA_GSE155698_sub.h5ad",
    ROOT / "data/analysis/scRNA_GSE156405.h5ad",
    ROOT / "data/analysis/scRNA_GSE212966.h5ad",
]
TUMOUR_TISSUES = {"PDAC", "PDAC_tumor", "PDAC_tissue"}
GENES = ["SIGLEC15", "ACP5", "MMP9"]
SEED = 20261001
N_RESAMPLES = 10_000
EXPECTED_CONTRACT_SHA256 = "fe15ecc8e2e235c59def41b657acdd6a14e90b751fa44b77a82d1a07ff0d0fc0"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def json_default(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    raise TypeError(type(value).__name__)


def dense(x) -> np.ndarray:
    return x.toarray() if sparse.issparse(x) else np.asarray(x)


def inverse_log_aggregate(values: np.ndarray) -> np.ndarray:
    return np.log1p(np.mean(np.expm1(np.asarray(values, float)), axis=0))


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, float)
    sd = np.nanstd(values, ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return np.full(values.shape, np.nan)
    return (values - np.nanmean(values)) / sd


def prepare_samples() -> tuple[pd.DataFrame, dict[str, dict[str, np.ndarray]]]:
    rows: list[dict] = []
    cell_payload: dict[str, dict[str, np.ndarray]] = {}
    for path in PATHS:
        a = ad.read_h5ad(path, backed="r")
        if any(gene not in a.var_names for gene in GENES):
            raise RuntimeError(f"Required genes missing from {path.name}")
        obs = a.obs.copy()
        tumour = obs["tissue"].astype(str).isin(TUMOUR_TISSUES).to_numpy()
        myeloid = tumour & obs["celltype"].astype(str).eq("Myeloid").to_numpy()
        dataset_values = obs.loc[myeloid, "dataset"].astype(str).unique()
        if len(dataset_values) != 1:
            raise RuntimeError(f"Unexpected dataset labels in {path.name}: {dataset_values}")
        dataset = str(dataset_values[0])

        my_obs = obs.loc[myeloid, ["sample_id"]].copy().reset_index(drop=True)
        my_expr = dense(a[myeloid, GENES].X).astype(float, copy=False)
        tumour_obs = obs.loc[tumour, ["sample_id"]].copy().reset_index(drop=True)
        tumour_target = dense(a[tumour, ["SIGLEC15"]].X).astype(float, copy=False).reshape(-1)
        a.file.close()

        sample_ids = []
        target_parts = []
        sample_codes = np.full(len(my_obs), -1, dtype=int)
        for sample_id, index in my_obs.groupby("sample_id", sort=True).groups.items():
            idx = np.asarray(list(index), dtype=int)
            if len(idx) < 50:
                continue
            sample_id = str(sample_id)
            sample_code = len(sample_ids)
            sample_ids.append(sample_id)
            sample_codes[idx] = sample_code
            values = my_expr[idx]
            aggregate = inverse_log_aggregate(values)
            target = values[:, 0]
            positive = target > 0
            all_idx = np.flatnonzero(tumour_obs["sample_id"].astype(str).to_numpy() == sample_id)
            all_target = tumour_target[all_idx]
            rows.append(
                {
                    "dataset": dataset,
                    "sample_id": sample_id,
                    "n_myeloid": len(idx),
                    "n_all_tumour": len(all_idx),
                    "n_siglec15_positive": int(positive.sum()),
                    "siglec15_positive_fraction": float(positive.mean()),
                    "siglec15_positive_mean": float(np.mean(target[positive])) if positive.any() else np.nan,
                    "siglec15_myeloid_aggregate": float(aggregate[0]),
                    "acp5_myeloid_aggregate": float(aggregate[1]),
                    "mmp9_myeloid_aggregate": float(aggregate[2]),
                    "siglec15_all_tumour_aggregate": float(inverse_log_aggregate(all_target[:, None])[0]),
                }
            )
            target_parts.append(target)

        keep = sample_codes >= 0
        cell_payload[dataset] = {
            "sample_ids": np.asarray(sample_ids, dtype=object),
            "sample_code": sample_codes[keep],
            "target": my_expr[keep, 0].astype(float, copy=True),
            "n_by_sample": np.bincount(sample_codes[keep], minlength=len(sample_ids)),
        }

    frame = pd.DataFrame(rows).sort_values(["dataset", "sample_id"]).reset_index(drop=True)
    frame["acp5_mmp9_burden"] = np.nan
    for _, idx in frame.groupby("dataset", sort=True).groups.items():
        idx = np.asarray(list(idx), dtype=int)
        acp5 = zscore(frame.loc[idx, "acp5_myeloid_aggregate"].to_numpy(float))
        mmp9 = zscore(frame.loc[idx, "mmp9_myeloid_aggregate"].to_numpy(float))
        frame.loc[idx, "acp5_mmp9_burden"] = np.nanmean(np.column_stack([acp5, mmp9]), axis=1)
    return frame, cell_payload


def pooled_rank(frame: pd.DataFrame, x_col: str, y_col: str, min_stratum_n: int = 2) -> tuple[float, int, str]:
    xs, ys, used = [], [], []
    for dataset, group in frame.groupby("dataset", sort=True):
        group = group[[x_col, y_col]].dropna()
        if len(group) < min_stratum_n:
            continue
        x = stats.rankdata(group[x_col].to_numpy(float)).astype(float)
        y = stats.rankdata(group[y_col].to_numpy(float)).astype(float)
        x_sd, y_sd = x.std(ddof=0), y.std(ddof=0)
        if x_sd == 0 or y_sd == 0:
            continue
        xs.append((x - x.mean()) / x_sd)
        ys.append((y - y.mean()) / y_sd)
        used.append(f"{dataset}:{len(group)}")
    if not xs:
        return np.nan, 0, ""
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    return float(np.corrcoef(x, y)[0, 1]), len(x), ";".join(used)


def combine_rank_arrays(parts: list[tuple[np.ndarray, np.ndarray]]) -> float:
    """Combine within-stratum standardised ranks without DataFrame overhead."""
    xs, ys = [], []
    for x_raw, y_raw in parts:
        if len(x_raw) < 2:
            continue
        x = stats.rankdata(x_raw).astype(float)
        y = stats.rankdata(y_raw).astype(float)
        x_sd, y_sd = x.std(ddof=0), y.std(ddof=0)
        if x_sd == 0 or y_sd == 0:
            continue
        xs.append((x - x.mean()) / x_sd)
        ys.append((y - y.mean()) / y_sd)
    if not xs:
        return np.nan
    return float(np.corrcoef(np.concatenate(xs), np.concatenate(ys))[0, 1])


def stratified_inference(frame: pd.DataFrame, x_col: str, y_col: str, seed: int) -> dict:
    frame = frame[["dataset", x_col, y_col]].dropna().copy()
    observed, n_used, composition = pooled_rank(frame, x_col, y_col)
    rng = np.random.default_rng(seed)
    boot = np.full(N_RESAMPLES, np.nan)
    perm = np.full(N_RESAMPLES, np.nan)
    groups = [
        (group[x_col].to_numpy(float), group[y_col].to_numpy(float))
        for _, group in frame.groupby("dataset", sort=True)
    ]
    for b in range(N_RESAMPLES):
        draw_parts, perm_parts = [], []
        for x, y in groups:
            idx = rng.integers(0, len(x), len(x))
            draw_parts.append((x[idx], y[idx]))
            perm_parts.append((x, rng.permutation(y)))
        boot[b] = combine_rank_arrays(draw_parts)
        perm[b] = combine_rank_arrays(perm_parts)
    finite_boot = boot[np.isfinite(boot)]
    finite_perm = perm[np.isfinite(perm)]
    ci = np.quantile(finite_boot, [0.025, 0.975])
    p = (np.sum(np.abs(finite_perm) >= abs(observed)) + 1) / (len(finite_perm) + 1)
    return {
        "rho": observed,
        "ci_low": float(ci[0]),
        "ci_high": float(ci[1]),
        "permutation_p": float(p),
        "n_used": n_used,
        "composition": composition,
        "n_boot_finite": len(finite_boot),
        "n_perm_finite": len(finite_perm),
    }


def leave_one_dataset_out(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for omitted in sorted(frame["dataset"].unique()):
        subset = frame.loc[frame["dataset"] != omitted]
        rho, n, composition = pooled_rank(subset, "siglec15_myeloid_aggregate", "acp5_mmp9_burden")
        rows.append({"omitted_dataset": omitted, "rho": rho, "n_used": n, "composition": composition})
    return pd.DataFrame(rows)


def cell_exchangeability_null(
    frame: pd.DataFrame, payload: dict[str, dict[str, np.ndarray]], observed: float
) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(SEED + 700)
    frame_lookup = frame.set_index(["dataset", "sample_id"])
    burden_by_dataset = {
        dataset: frame_lookup.loc[
            [(dataset, str(sample_id)) for sample_id in payload[dataset]["sample_ids"]],
            "acp5_mmp9_burden",
        ].to_numpy(float)
        for dataset in sorted(payload)
    }
    values = np.full(N_RESAMPLES, np.nan)
    for b in range(N_RESAMPLES):
        parts = []
        for dataset in sorted(payload):
            item = payload[dataset]
            shuffled = rng.permutation(item["target"])
            linear = np.expm1(shuffled)
            sums = np.bincount(item["sample_code"], weights=linear, minlength=len(item["sample_ids"]))
            aggregate = np.log1p(sums / item["n_by_sample"])
            parts.append((aggregate, burden_by_dataset[dataset]))
        values[b] = combine_rank_arrays(parts)
    finite = values[np.isfinite(values)]
    summary = {
        "observed_rho": observed,
        "n_simulations": len(finite),
        "null_mean": float(np.mean(finite)),
        "null_sd": float(np.std(finite, ddof=1)),
        "null_q025": float(np.quantile(finite, 0.025)),
        "null_q50": float(np.quantile(finite, 0.5)),
        "null_q975": float(np.quantile(finite, 0.975)),
        "two_sided_tail_probability": float((np.sum(np.abs(finite) >= abs(observed)) + 1) / (len(finite) + 1)),
    }
    return pd.DataFrame({"simulation": np.arange(1, len(values) + 1), "rho": values}), summary


def reproduction_check(frame: pd.DataFrame) -> dict:
    archived_samples = pd.read_csv(ARCHIVE / "myeloid_sample_aggregates_scored.csv")
    archived_assoc = pd.read_csv(ARCHIVE / "myeloid_aggregate_associations.csv")
    merged = frame.merge(
        archived_samples[["dataset", "sample_id", "SIGLEC15", "Derosiers_ACP5_MMP9_core"]],
        on=["dataset", "sample_id"],
        how="outer",
        indicator=True,
    )
    primary = archived_assoc.loc[
        (archived_assoc["scope"] == "pooled_within_dataset_ranks")
        & (archived_assoc["program"] == "Derosiers_ACP5_MMP9_core")
    ].iloc[0]
    rho, n, composition = pooled_rank(frame, "siglec15_myeloid_aggregate", "acp5_mmp9_burden")
    result = {
        "sample_keys_match": bool((merged["_merge"] == "both").all() and len(merged) == 20),
        "max_abs_r1_difference": float(np.nanmax(np.abs(merged["siglec15_myeloid_aggregate"] - merged["SIGLEC15"]))),
        "max_abs_burden_difference": float(np.nanmax(np.abs(merged["acp5_mmp9_burden"] - merged["Derosiers_ACP5_MMP9_core"]))),
        "computed_rho": rho,
        "archived_rho": float(primary["rho"]),
        "absolute_rho_difference": abs(rho - float(primary["rho"])),
        "n_used": n,
        "composition": composition,
    }
    if not result["sample_keys_match"] or result["max_abs_r1_difference"] > 1e-10 or result["max_abs_burden_difference"] > 1e-10 or result["absolute_rho_difference"] > 1e-12:
        raise RuntimeError(f"Archived result was not reproduced: {result}")
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if sha256(CONTRACT) != EXPECTED_CONTRACT_SHA256:
        raise RuntimeError("Frozen analysis contract hash mismatch")
    frame, payload = prepare_samples()
    reproduction = reproduction_check(frame)

    specifications = [
        ("primary_all_eligible", frame, SEED + 10),
        ("sparse_sensitivity_n_positive_ge_5", frame.loc[frame["n_siglec15_positive"] >= 5], SEED + 20),
        ("sparse_sensitivity_n_positive_ge_10", frame.loc[frame["n_siglec15_positive"] >= 10], SEED + 30),
    ]
    rows = []
    for name, subset, seed in specifications:
        result = stratified_inference(subset, "siglec15_myeloid_aggregate", "acp5_mmp9_burden", seed)
        rows.append({"specification": name, "n_available": len(subset), **result})
    primary = rows[0]["rho"]

    secondary_rows = []
    for index, (readout, label) in enumerate(
        [
            ("siglec15_positive_fraction", "R2_positive_fraction"),
            ("siglec15_positive_mean", "R3_positive_cell_mean"),
            ("siglec15_all_tumour_aggregate", "R4_all_tumour_aggregate"),
        ]
    ):
        result = stratified_inference(frame, readout, "acp5_mmp9_burden", SEED + 100 + index)
        secondary_rows.append({"readout": label, "column": readout, "n_available": int(frame[readout].notna().sum()), **result})

    null_draws, null_summary = cell_exchangeability_null(frame, payload, primary)
    loo = leave_one_dataset_out(frame)

    frame.to_csv(OUT / "v72_sample_level_burden_and_readouts.csv", index=False)
    pd.DataFrame(rows).to_csv(OUT / "v72_primary_and_sparsity_estimates.csv", index=False)
    pd.DataFrame(secondary_rows).to_csv(OUT / "v72_secondary_readout_estimates.csv", index=False)
    loo.to_csv(OUT / "v72_leave_one_dataset_out.csv", index=False)
    null_draws.to_csv(OUT / "v72_cell_exchangeability_null_draws.csv", index=False)
    pd.DataFrame([null_summary]).to_csv(OUT / "v72_cell_exchangeability_null_summary.csv", index=False)

    input_manifest = pd.DataFrame(
        [{"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)} for path in PATHS]
    )
    input_manifest.to_csv(OUT / "v72_input_manifest.csv", index=False)
    summary = {
        "contract": str(CONTRACT),
        "contract_sha256": sha256(CONTRACT),
        "seed": SEED,
        "n_resamples": N_RESAMPLES,
        "reproduction": reproduction,
        "primary_and_sparsity": rows,
        "secondary_readouts": secondary_rows,
        "cell_exchangeability_null": null_summary,
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "anndata": ad.__version__,
        },
        "claim_ceiling": (
            "Tumour-level SIGLEC15 RNA did not establish stable preservation of the between-tumour "
            "ordering of the myeloid ACP5/MMP9 aggregate burden across the available samples and "
            "sparsity specifications."
        ),
    }
    (OUT / "V72_ANALYSIS_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=json_default) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=json_default))


if __name__ == "__main__":
    main()
