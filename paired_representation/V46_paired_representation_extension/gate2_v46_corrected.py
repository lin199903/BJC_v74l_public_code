"""Corrected Gate 2 inference for the v46 manuscript.

The original composition sensitivity residualised only the protein-state
outcome.  This script implements a conventional partial Spearman correlation:
all variables are rank transformed, both the SIGLEC15 RNA exposure and the
protein-state outcome are residualised on the same composition proxies, and
the two residual series are correlated.  It also records bootstrap intervals
and construction sensitivities for the protein-state score.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(r".")
PG = ROOT / "data" / "analysis" / "proteogenomic"
OUT = (
    ROOT
    / "npjPrecOncol"
    / "FINAL_SUBMISSION_v6"
    / "BMC_MEDICINE_SUBMISSION_v46"
    / "Analysis"
    / "V46_validation_extension"
    / "Gate2"
)
OUT.mkdir(parents=True, exist_ok=True)

STATE = ["ACP5", "MMP9", "CTSK", "ITGB3", "NFATC1", "CTSB", "MMP14"]
CORE = ["ACP5", "MMP9"]
MYELOID_PROXY = ["AIF1", "CSF1R", "TYROBP"]
EPITHELIAL_PROXY = ["EPCAM", "KRT8"]
SEED = 15072026
N_BOOT = 10_000


def read_cct(path: Path) -> pd.DataFrame:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        lines = [line.rstrip("\n") for line in handle if not line.startswith("#")]
    header = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:] if line.strip()]
    return pd.DataFrame(
        [row[1:] for row in rows],
        index=[row[0] for row in rows],
        columns=header[1:],
    ).apply(pd.to_numeric, errors="coerce")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def zscore(series: pd.Series) -> pd.Series:
    return (series - series.mean()) / series.std(ddof=1)


def score(frame: pd.DataFrame, genes: list[str], complete: bool = False) -> pd.Series:
    matrix = pd.concat([zscore(frame.loc[gene]).rename(gene) for gene in genes], axis=1)
    if complete:
        return matrix.dropna().mean(axis=1)
    return matrix.mean(axis=1, skipna=True)


def spearman(x: pd.Series, y: pd.Series) -> dict[str, float | int]:
    paired = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
    result = stats.spearmanr(paired["x"], paired["y"])
    return {"n": int(len(paired)), "rho": float(result.statistic), "p": float(result.pvalue)}


def residuals_on_ranked_covariates(
    outcome: pd.Series, covariates: pd.DataFrame
) -> pd.Series:
    joined = pd.concat([outcome.rename("outcome"), covariates], axis=1).dropna()
    ranked = joined.rank(method="average")
    xmat = np.column_stack(
        [np.ones(len(ranked))]
        + [ranked[column].to_numpy(float) for column in covariates.columns]
    )
    beta, *_ = np.linalg.lstsq(xmat, ranked["outcome"].to_numpy(float), rcond=None)
    return pd.Series(
        ranked["outcome"].to_numpy(float) - xmat @ beta,
        index=ranked.index,
        name=outcome.name,
    )


def partial_spearman(
    x: pd.Series, y: pd.Series, covariates: pd.DataFrame
) -> dict[str, float | int]:
    joined = pd.concat([x.rename("x"), y.rename("y"), covariates], axis=1).dropna()
    x_resid = residuals_on_ranked_covariates(joined["x"], joined[covariates.columns])
    y_resid = residuals_on_ranked_covariates(joined["y"], joined[covariates.columns])
    rho = float(stats.pearsonr(x_resid, y_resid).statistic)
    k = len(covariates.columns)
    degrees_freedom = len(joined) - k - 2
    t_value = rho * np.sqrt(degrees_freedom / max(1.0 - rho**2, np.finfo(float).eps))
    p_value = float(2 * stats.t.sf(abs(t_value), degrees_freedom))
    return {
        "n": int(len(joined)),
        "rho": rho,
        "p": p_value,
        "degrees_freedom": int(degrees_freedom),
        "n_covariates": int(k),
    }


def bootstrap_interval(
    frame: pd.DataFrame,
    estimator,
    seed_offset: int,
) -> dict[str, float | int]:
    rng = np.random.default_rng(SEED + seed_offset)
    values: list[float] = []
    for _ in range(N_BOOT):
        indices = rng.integers(0, len(frame), len(frame))
        sample = frame.iloc[indices].reset_index(drop=True)
        try:
            value = float(estimator(sample))
        except (ValueError, np.linalg.LinAlgError):
            continue
        if np.isfinite(value):
            values.append(value)
    array = np.asarray(values)
    return {
        "bootstrap_seed": int(SEED + seed_offset),
        "n_boot_valid": int(len(array)),
        "ci_low": float(np.percentile(array, 2.5)),
        "ci_high": float(np.percentile(array, 97.5)),
    }


def main() -> None:
    protein_path = PG / "proteome_tumor.cct"
    rna_path = PG / "rna_tumor.cct"
    protein = read_cct(protein_path)
    rna = read_cct(rna_path)
    common = [sample for sample in rna.columns if sample in set(protein.columns)]
    protein = protein[common]
    rna = rna[common]

    state = score(protein, STATE)
    core = score(protein, CORE)
    complete_state = score(protein, STATE, complete=True)
    myeloid = score(protein, MYELOID_PROXY)
    epithelial = score(protein, EPITHELIAL_PROXY)

    frame = pd.DataFrame(
        {
            "rna_SIGLEC15": rna.loc["SIGLEC15"],
            "protein_state_7gene": state,
            "protein_state_core": core,
            "protein_myeloid_proxy": myeloid,
            "protein_epithelial_proxy": epithelial,
            "rna_EPCAM": rna.loc["EPCAM"],
            "protein_EPCAM": protein.loc["EPCAM"],
        }
    )

    unadjusted = spearman(frame["rna_SIGLEC15"], frame["protein_state_7gene"])
    unadjusted_frame = frame[["rna_SIGLEC15", "protein_state_7gene"]].dropna()
    unadjusted.update(
        bootstrap_interval(
            unadjusted_frame,
            lambda data: stats.spearmanr(
                data["rna_SIGLEC15"], data["protein_state_7gene"]
            ).statistic,
            0,
        )
    )

    covariate_sets = {
        "myeloid_proxy": ["protein_myeloid_proxy"],
        "myeloid_and_epithelial_proxies": [
            "protein_myeloid_proxy",
            "protein_epithelial_proxy",
        ],
    }
    adjusted: dict[str, dict[str, float | int]] = {}
    for offset, (name, columns) in enumerate(covariate_sets.items(), start=1):
        subset = frame[["rna_SIGLEC15", "protein_state_7gene", *columns]].dropna()
        result = partial_spearman(
            subset["rna_SIGLEC15"],
            subset["protein_state_7gene"],
            subset[columns],
        )
        result.update(
            bootstrap_interval(
                subset,
                lambda data, cols=columns: partial_spearman(
                    data["rna_SIGLEC15"],
                    data["protein_state_7gene"],
                    data[cols],
                )["rho"],
                offset,
            )
        )
        adjusted[name] = result

    construction_sensitivity: dict[str, dict[str, float | int]] = {}
    core_frame = frame[["rna_SIGLEC15", "protein_state_core"]].dropna()
    core_result = spearman(core_frame["rna_SIGLEC15"], core_frame["protein_state_core"])
    core_result.update(
        bootstrap_interval(
            core_frame,
            lambda data: stats.spearmanr(
                data["rna_SIGLEC15"], data["protein_state_core"]
            ).statistic,
            10,
        )
    )
    construction_sensitivity["ACP5_MMP9_core"] = core_result
    complete_frame = pd.DataFrame(
        {"rna_SIGLEC15": rna.loc["SIGLEC15"], "protein_state_complete": complete_state}
    ).dropna()
    complete_result = spearman(
        complete_frame["rna_SIGLEC15"], complete_frame["protein_state_complete"]
    )
    complete_result.update(
        bootstrap_interval(
            complete_frame,
            lambda data: stats.spearmanr(
                data["rna_SIGLEC15"], data["protein_state_complete"]
            ).statistic,
            11,
        )
    )
    construction_sensitivity["complete_case_7gene"] = complete_result
    construction_sensitivity["leave_one_protein_out"] = {
        omitted: spearman(rna.loc["SIGLEC15"], score(protein, [g for g in STATE if g != omitted]))
        for omitted in STATE
    }

    controls = {
        "myeloid_proxy_vs_state": spearman(
            frame["protein_myeloid_proxy"], frame["protein_state_7gene"]
        ),
        "EPCAM_RNA_vs_protein": spearman(frame["rna_EPCAM"], frame["protein_EPCAM"]),
    }
    myeloid_control_frame = frame[
        ["protein_myeloid_proxy", "protein_state_7gene"]
    ].dropna()
    controls["myeloid_proxy_vs_state"].update(
        bootstrap_interval(
            myeloid_control_frame,
            lambda data: stats.spearmanr(
                data["protein_myeloid_proxy"], data["protein_state_7gene"]
            ).statistic,
            20,
        )
    )
    epcam_control_frame = frame[["rna_EPCAM", "protein_EPCAM"]].dropna()
    controls["EPCAM_RNA_vs_protein"].update(
        bootstrap_interval(
            epcam_control_frame,
            lambda data: stats.spearmanr(
                data["rna_EPCAM"], data["protein_EPCAM"]
            ).statistic,
            21,
        )
    )

    result = {
        "analysis_identity": "paired orthogonal representation test",
        "method_correction": {
            "original_issue": (
                "The earlier sensitivity residualised only the protein-state outcome."
            ),
            "corrected_method": (
                "Partial Spearman correlation after rank transformation and "
                "residualisation of both exposure and outcome on the same proxies."
            ),
        },
        "inputs": {
            "proteome_tumor.cct": {
                "bytes": protein_path.stat().st_size,
                "sha256": sha256(protein_path),
            },
            "rna_tumor.cct": {
                "bytes": rna_path.stat().st_size,
                "sha256": sha256(rna_path),
            },
            "shared_samples": len(common),
        },
        "state_proteins": STATE,
        "unadjusted": unadjusted,
        "partial_spearman": adjusted,
        "construction_sensitivity": construction_sensitivity,
        "positive_controls": controls,
        "claim_boundary": (
            "Tumour-level SIGLEC15 RNA versus a downstream bulk-tissue protein-state "
            "score in the same patients; SIGLEC15 protein was not quantified and the "
            "protein score is not myeloid-cell specific."
        ),
    }

    frame.to_csv(OUT / "GATE2_v46_corrected_analysis.csv", index_label="sample_id")
    (OUT / "_gate2_v46_corrected_summary.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
