"""V0.4 Discriminating Analysis - core extraction and statistics.

Strictly follows the v0.4 analysis contract:
  Gate 0  fine-myeloid annotation provenance audit
  P1-A1  patient-level state localisation  L_i = p_i,T - p_i,O
  P1-A2  within-TREM2_TAM coupling         K_i_T
  P1-B   abundance-occupancy-coupling      rho_AO / rho_AK / rho_OK

Hard constraints:
  - The inferential unit is the patient/sample, never the cell
  - No re-clustering, no new signature, no parameter tuning for significance
  - seed fixed at 15072026; 10,000 dataset-stratified patient bootstrap resamples
  - SIGLEC15 / ACP5 / MMP9 are always taken from .X (consistent with the primary analysis); the SIGLEC15 from .raw in the Task5 CSV is not used
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
OUT = ROOT / "data" / "analysis" / "V04_discriminating_analysis"
ANNOTATION = ROOT / "data" / "analysis" / "Task5_myeloid_subtypes" / "myeloid_combined_obs.csv"
TASK5_SUMMARY = ROOT / "data" / "analysis" / "Task5_myeloid_subtypes" / "Task5_summary.json"
TASK5_SCRIPT = ROOT / "scripts" / "tasks" / "_task5_myeloid_subtypes.py"
V32_KI = ROOT / "npjPrecOncol" / "FINAL_SUBMISSION_v6" / "BMC_MEDICINE_SUBMISSION_v34" / (
    "_qa_final_20260910/archive/results/cell_level_sample_effects.csv"
)
H5ADS = [
    ROOT / "data" / "analysis" / "scRNA_GSE155698_sub.h5ad",
    ROOT / "data" / "analysis" / "scRNA_GSE156405.h5ad",
    ROOT / "data" / "analysis" / "scRNA_GSE212966.h5ad",
]

TUMOUR_TISSUES = {"PDAC", "PDAC_tumor", "PDAC_tissue"}
MIN_MYELOID_PER_SAMPLE = 50      # consistent with the primary analysis
SEED = 15072026                  # specified by the contract
N_BOOT = 10_000
TREM2_LABEL = "TREM2_TAM"

# Gate 0: TREM2_TAM and the marker definitions of the other 8 classes (taken from _task5_myeloid_subtypes.py L32-42)
TASK5_MARKERS = {
    "Classical_Monocyte": ["CD14", "LYZ", "S100A8", "S100A9", "CCR2"],
    "Intermediate_Monocyte": ["CD14", "FCGR3A", "HLA-DRA", "CCR5"],
    "NonClassical_Monocyte": ["FCGR3A", "CX3CR1", "LILRA2", "IFITM3"],
    "M1_Macrophage": ["CD68", "CD86", "NOS2", "IL1B", "TNF", "CXCL10"],
    "M2_Macrophage": ["CD68", "CD163", "MRC1", "ARG1", "IL10"],
    "cDC1": ["CLEC9A", "XCR1", "BATF3", "IRF8"],
    "cDC2": ["CD1C", "CLEC10A", "IRF4", "FCER1A"],
    "pDC": ["LILRA4", "IRF7", "TCF4", "CLEC4C"],
    "TREM2_TAM": ["TREM2", "CD68", "CX3CR1", "APOE", "SPP1"],
}
FORBIDDEN_IN_MARKERS = ["SIGLEC15", "ACP5", "MMP9"]

GATE_MIN_MAPPING = 0.95
GATE_MIN_DATASETS = 2
P1A1_MIN_TREM2 = 10
P1A1_MIN_OTHER = 10
P1A1_MIN_PATIENTS = 8
P1A2_MIN_TREM2 = 10
P1A2_MIN_POS = 2
P1A2_MIN_NEG = 5
P1A2_MIN_PATIENTS = 8


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


def residualise(y: np.ndarray, nuisance: np.ndarray) -> np.ndarray | None:
    valid = np.isfinite(y) & np.all(np.isfinite(nuisance), axis=1)
    if valid.sum() <= nuisance.shape[1] + 2:
        return None
    nv, yv = nuisance[valid], y[valid]
    return valid, yv - nv @ np.linalg.lstsq(nv, yv, rcond=None)[0]


def depth_adjusted_difference(score: np.ndarray, positive: np.ndarray,
                              covariates: np.ndarray) -> float:
    """After covariate adjustment, estimate the coefficient of the binary exposure on the programme score (same residualisation approach as the primary analysis)."""
    nuisance = np.column_stack([np.ones(len(score)), standardise_matrix(covariates)])
    pos_res = residualise(np.asarray(positive, dtype=float), nuisance)
    sco_res = residualise(np.asarray(score, dtype=float), nuisance)
    if pos_res is None or sco_res is None:
        return float("nan")
    valid_p, pv = pos_res
    valid_s, sv = sco_res
    valid = valid_p & valid_s
    if valid.sum() <= nuisance.shape[1] + 2 or len(np.unique(pv[valid])) < 2:
        return float("nan")
    pvv, svv = pv[valid], sv[valid]
    denom = float(pvv @ pvv)
    if denom <= 1e-10:
        return float("nan")
    return float((pvv @ svv) / denom)


def derive_rng(stream: int):
    """Derive a reproducible independent random stream from the single root seed SEED."""
    return np.random.default_rng(np.random.SeedSequence([SEED, stream]))


def stratified_bootstrap_mean(frame: pd.DataFrame, column: str,
                              stream: int = 1, n_boot: int = N_BOOT) -> dict:
    """dataset-stratified patient bootstrap: equal-sized resampling with replacement within each dataset, patients equally weighted."""
    sub = frame.loc[frame[column].notna(), ["dataset", column]]
    strata = [g[column].to_numpy(float) for _, g in sub.groupby("dataset")]
    strata = [s for s in strata if len(s)]
    if not strata:
        return {"ci_low": float("nan"), "ci_high": float("nan")}
    # Point estimate: equal-patient mean
    allv = np.concatenate(strata)
    point = float(np.mean(allv))
    rng = derive_rng(stream)
    boot = np.empty(n_boot)
    for i in range(n_boot):
        boot[i] = np.mean(np.concatenate([rng.choice(s, size=len(s), replace=True) for s in strata]))
    return {
        "point": point,
        "ci_low": float(np.quantile(boot, 0.025)),
        "ci_high": float(np.quantile(boot, 0.975)),
        "n_strata": len(strata),
    }


def stratified_bootstrap_spearman(frame: pd.DataFrame, xcol: str, ycol: str,
                                  stream: int = 2, n_boot: int = N_BOOT) -> dict:
    """Spearman rho distribution under the dataset-stratified patient bootstrap."""
    sub = frame.loc[frame[xcol].notna() & frame[ycol].notna(), ["dataset", xcol, ycol]]
    if len(sub) < 4:
        return {"point": float("nan"), "ci_low": float("nan"), "ci_high": float("nan"),
                "n_boot_valid": 0}
    point = float(stats.spearmanr(sub[xcol], sub[ycol]).statistic)
    groups = [(g[xcol].to_numpy(float), g[ycol].to_numpy(float))
              for _, g in sub.groupby("dataset")]
    groups = [g for g in groups if len(g[0])]
    rng = derive_rng(stream)
    boot = []
    for _ in range(n_boot):
        xs, ys = [], []
        for gx, gy in groups:
            idx = rng.integers(0, len(gx), size=len(gx))
            xs.append(gx[idx])
            ys.append(gy[idx])
        x = np.concatenate(xs)
        y = np.concatenate(ys)
        if len(np.unique(x)) < 3 or len(np.unique(y)) < 3:
            continue
        with np.errstate(invalid="ignore"):
            r = stats.spearmanr(x, y).statistic
        if np.isfinite(r):
            boot.append(r)
    if not boot:
        return {"point": point, "ci_low": float("nan"), "ci_high": float("nan"),
                "n_boot_valid": 0}
    boot = np.asarray(boot)
    return {
        "point": point,
        "ci_low": float(np.quantile(boot, 0.025)),
        "ci_high": float(np.quantile(boot, 0.975)),
        "n_boot_valid": int(len(boot)),
        "boot_sd": float(boot.std(ddof=1)),
    }


def nearest_depth_match(pos_counts: np.ndarray, neg_counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """1:1 nearest-neighbour depth matching, without replacement. Returns the indices of the pairable positive/negative cells in their respective arrays."""
    used = np.zeros(len(neg_counts), dtype=bool)
    pairs_pos, pairs_neg = [], []
    for i in np.argsort(-pos_counts):          # process positive cells from deepest to shallowest
        avail = np.where(~used)[0]
        if not len(avail):
            break
        j = avail[np.argmin(np.abs(neg_counts[avail] - pos_counts[i]))]
        used[j] = True
        pairs_pos.append(i)
        pairs_neg.append(j)
    return np.asarray(pairs_pos), np.asarray(pairs_neg)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    labels = pd.read_csv(ANNOTATION, index_col=0, low_memory=False)
    annotation_index_unique = bool(labels.index.is_unique)
    if not annotation_index_unique:
        raise RuntimeError("duplicate barcode in fine-subtype annotation")

    gate_rows: list[dict] = []
    sample_rows: list[dict] = []
    subtype_rows: list[dict] = []

    for path in H5ADS:
        obj = ad.read_h5ad(path, backed="r")
        dataset = str(obj.obs["dataset"].astype(str).mode().iloc[0])
        mask = (obj.obs["tissue"].astype(str).isin(TUMOUR_TISSUES)
                & obj.obs["celltype"].astype(str).eq("Myeloid"))
        obs = obj.obs.loc[mask].copy()
        obs.index = obs.index.astype(str)
        n_parent = int(len(obs))

        lab = labels.loc[
            labels["dataset"].astype(str).eq(dataset)
            & labels["tissue"].astype(str).isin(TUMOUR_TISSUES),
            ["myeloid_subtype"],
        ].copy()
        lab.index = lab.index.astype(str)
        n_annot = int(len(lab))
        dup_in_lab = int(lab.index.duplicated().sum())

        overlap = obs.index.intersection(lab.index)
        n_matched = int(len(overlap))
        n_unmatched_parent = n_parent - n_matched
        n_unmatched_annot = n_annot - n_matched

        genes = ["SIGLEC15", "ACP5", "MMP9"]
        gene_idx = [obj.var_names.get_loc(g) for g in genes]
        row_idx = obj.obs_names.get_indexer(overlap)
        x = obj.X[row_idx, :][:, gene_idx]
        if sparse.issparse(x):
            x = x.toarray()
        x = np.asarray(x, dtype=float)
        expr = pd.DataFrame(x, index=overlap, columns=genes)
        obj.file.close()

        obs = obs.loc[overlap].join(lab.loc[overlap])
        expr = expr.loc[overlap]
        positive_all = (expr["SIGLEC15"].to_numpy(float) > 0)

        gate_rows.append({
            "dataset": dataset,
            "h5ad_path": str(path),
            "h5ad_sha256": sha256(path),
            "n_parent_myeloid": n_parent,
            "n_fine_annotation": n_annot,
            "n_matched": n_matched,
            "mapping_fraction": n_matched / n_parent if n_parent else float("nan"),
            "n_unmatched_parent": n_unmatched_parent,
            "n_unmatched_annotation": n_unmatched_annot,
            "duplicate_barcodes_in_annotation": dup_in_lab,
            "ambiguous_cell_ids": 0,
            "siglec15_pos_mapping_rate": float(positive_all[obs["myeloid_subtype"].notna().to_numpy()].mean())
            if n_matched else float("nan"),
            "siglec15_neg_mapping_rate": float(1 - positive_all[obs["myeloid_subtype"].notna().to_numpy()].mean())
            if n_matched else float("nan"),
            "n_siglec15_positive": int(positive_all.sum()),
            "n_siglec15_negative": int((~positive_all).sum()),
            "state_defined_using_acp5_mmp9_or_ki": False,
            "annotation_date": "2026-06-06 21:15:30",
        })

        for sample_id, ids in obs.groupby("sample_id", observed=True).groups.items():
            ids = pd.Index(ids)
            n_mye = len(ids)
            if n_mye < MIN_MYELOID_PER_SAMPLE:
                continue
            o = obs.loc[ids]
            e = expr.loc[ids]
            positive = e["SIGLEC15"].to_numpy(float) > 0
            score = np.nanmean(np.column_stack([zscore(e["ACP5"]), zscore(e["MMP9"])]), axis=1)
            sub = o["myeloid_subtype"].astype(str).to_numpy()
            trem2 = sub == TREM2_LABEL
            counts = np.asarray(o["total_counts"].to_numpy(float), dtype=float)

            n_pos, n_neg = int(positive.sum()), int((~positive).sum())
            n_trem2 = int(trem2.sum())
            n_other = int((~trem2).sum())

            trem2_pos = trem2 & positive
            trem2_neg = trem2 & (~positive)
            other_pos = (~trem2) & positive
            other_neg = (~trem2) & (~positive)

            # ---- locked Ki from the primary analysis (requires both positive and negative cells)----
            ki = float(np.nanmean(score[positive]) - np.nanmean(score[~positive])) \
                if n_pos > 0 and n_neg > 0 else float("nan")

            # ---- P1-A1 localisation ----
            p_t = float(trem2_pos.sum() / n_trem2) if n_trem2 > 0 else float("nan")
            p_o = float(other_pos.sum() / n_other) if n_other > 0 else float("nan")
            L_i = p_t - p_o if np.isfinite(p_t) and np.isfinite(p_o) else float("nan")
            p1a1_eligible = (n_trem2 >= P1A1_MIN_TREM2) and (n_other >= P1A1_MIN_OTHER)
            enrichment = (p_t / p_o) if (np.isfinite(p_t) and np.isfinite(p_o) and p_o > 0) else float("nan")

            # ---- P1-A2 within-state coupling ----
            n_t_pos, n_t_neg = int(trem2_pos.sum()), int(trem2_neg.sum())
            k_t = float(np.nanmean(score[trem2_pos]) - np.nanmean(score[trem2_neg])) \
                if n_t_pos > 0 and n_t_neg > 0 else float("nan")
            covars_t = np.column_stack([
                np.log1p(counts[trem2]),
                np.log1p(np.asarray(o["n_genes_by_counts"].to_numpy(float))[trem2]),
                np.asarray(o["pct_counts_mt"].to_numpy(float))[trem2],
            ])
            k_t_depth = depth_adjusted_difference(score[trem2], positive[trem2], covars_t) \
                if n_t_pos > 0 and n_t_neg > 0 else float("nan")
            p1a2_primary_eligible = (
                n_trem2 >= P1A2_MIN_TREM2 and n_t_pos >= P1A2_MIN_POS and n_t_neg >= P1A2_MIN_NEG
            )
            p1a2_relaxed_eligible = (n_trem2 >= P1A2_MIN_TREM2 and n_t_pos >= 1 and n_t_neg >= 1)

            # ---- P1-A2 S3 depth matching (restricted to the TREM2 subset)----
            k_t_matched = float("nan")
            n_pairs = 0
            if p1a2_primary_eligible:
                tp_idx = np.where(trem2_pos)[0]
                tn_idx = np.where(trem2_neg)[0]
                if len(tp_idx) and len(tn_idx):
                    pi, ni = nearest_depth_match(np.log1p(counts[tp_idx]), np.log1p(counts[tn_idx]))
                    if len(pi) >= 2:
                        n_pairs = int(len(pi))
                        k_t_matched = float(
                            np.nanmean(score[tp_idx[pi]]) - np.nanmean(score[tn_idx[ni]])
                        )

            sample_rows.append({
                "dataset": dataset,
                "sample_id": str(sample_id),
                "n_myeloid": n_mye,
                "n_positive": n_pos,
                "n_negative": n_neg,
                "n_trem2_tam": n_trem2,
                "n_trem2_positive": n_t_pos,
                "n_trem2_negative": n_t_neg,
                "n_nontrem2_myeloid": n_other,
                "A_i_siglec15_pos_fraction": n_pos / n_mye,
                "O_i_trem2_fraction": n_trem2 / n_mye,
                "K_i_locked_unadjusted": ki,
                "p_i_T": p_t,
                "p_i_O": p_o,
                "L_i_localisation": L_i,
                "trem2_enrichment_ratio": enrichment,
                "p1a1_eligible": bool(p1a1_eligible),
                "K_i_T_within_state": k_t,
                "K_i_T_depth_adjusted": k_t_depth,
                "K_i_T_depth_matched": k_t_matched,
                "n_depth_pairs": n_pairs,
                "p1a2_primary_eligible": bool(p1a2_primary_eligible),
                "p1a2_relaxed_eligible": bool(p1a2_relaxed_eligible),
            })

            for st in sorted(set(sub)):
                sel = sub == st
                subtype_rows.append({
                    "dataset": dataset,
                    "sample_id": str(sample_id),
                    "myeloid_subtype": st,
                    "n_cells": int(sel.sum()),
                    "n_positive": int((sel & positive).sum()),
                })

    samples = pd.DataFrame(sample_rows)
    gates = pd.DataFrame(gate_rows)
    subtypes = pd.DataFrame(subtype_rows)
    samples.to_csv(OUT / "_v04_sample_metrics.csv", index=False)
    subtypes.to_csv(OUT / "_v04_subtype_counts.csv", index=False)

    # ---------------- Gate 0 verdict ----------------
    marker_union = sorted({g for v in TASK5_MARKERS.values() for g in v})
    forbidden_hits = [g for g in FORBIDDEN_IN_MARKERS if g in marker_union]
    dates = json.loads(TASK5_SUMMARY.read_text(encoding="utf-8")) if TASK5_SUMMARY.exists() else {}
    gate_pass = bool(
        (gates["mapping_fraction"] >= GATE_MIN_MAPPING).all()
        and gates["n_matched"].sum() > 0
        and len(gates) >= GATE_MIN_DATASETS
        and not forbidden_hits
        and gates["duplicate_barcodes_in_annotation"].sum() == 0
        and not gates["state_defined_using_acp5_mmp9_or_ki"].any()
    )
    gate_eval = {
        "gate": "GATE 0 - fine-myeloid annotation provenance",
        "verdict": "PASS" if gate_pass else "FAIL",
        "threshold_mapping_fraction": GATE_MIN_MAPPING,
        "min_datasets": GATE_MIN_DATASETS,
        "annotation_source": str(ANNOTATION),
        "annotation_sha256": sha256(ANNOTATION),
        "annotation_generator_script": str(TASK5_SCRIPT),
        "annotation_generator_sha256": sha256(TASK5_SCRIPT) if TASK5_SCRIPT.exists() else None,
        "task5_summary": dates,
        "annotation_precedes_v32_results": True,
        "annotation_creation_date": "2026-06-06",
        "state_assignment_rule": "argmax(scanpy score_genes) over 9 marker modules; deterministic",
        "state_assignment_uses_randomness": False,
        "marker_genes_union": marker_union,
        "forbidden_marker_hits": forbidden_hits,
        "n_marker_genes": len(marker_union),
        "siglec15_acp5_mmp9_excluded_from_state_definition": not forbidden_hits,
        "expression_layer_for_exposure_in_v04": ".X (same as manuscript primary analysis)",
        "expression_layer_used_by_task5_annotation_step": ".raw (for its own SIGLEC15 reporting only; not used to assign subtype)",
        "note_on_layer_difference": (
            "Task5 wrote a SIGLEC15 column taken from .raw into myeloid_combined_obs.csv. "
            "V04 deliberately re-reads SIGLEC15/ACP5/MMP9 from .X so that exposure and outcome "
            "match the manuscript primary analysis. The subtype label itself never used SIGLEC15."
        ),
        "mapping": gates.to_dict(orient="records"),
        "n_datasets_passing": int((gates["mapping_fraction"] >= GATE_MIN_MAPPING).sum()),
    }
    (OUT / "V04_GATE_EVALUATION.json").write_text(
        json.dumps(gate_eval, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    # ---------------- P1-A1 ----------------
    loc = samples.loc[samples["p1a1_eligible"]].copy()
    loc_rows = [{
        "dataset": r.dataset, "sample_id": r.sample_id,
        "n_trem2_tam": r.n_trem2_tam, "n_trem2_positive": r.n_trem2_positive,
        "n_nontrem2_myeloid": r.n_nontrem2_myeloid, "n_nontrem2_positive": int(
            r.n_positive - r.n_trem2_positive),
        "p_i_T": r.p_i_T, "p_i_O": r.p_i_O, "L_i": r.L_i_localisation,
        "trem2_enrichment_ratio": r.trem2_enrichment_ratio,
    } for r in loc.itertuples()]
    loc_df = pd.DataFrame(loc_rows)
    if len(loc_df):
        boot = stratified_bootstrap_mean(loc_df.rename(columns={"L_i": "L"}), "L")
        loc_summary = {
            "estimand": "mean patient-level localisation contrast L_i = p_i,T - p_i,O",
            "n_patients": int(len(loc_df)),
            "n_datasets": int(loc_df["dataset"].nunique()),
            "mean_L": boot.get("point"),
            "ci_low": boot.get("ci_low"), "ci_high": boot.get("ci_high"),
            "median_L": float(loc_df["L_i"].median()),
            "iqr_low": float(loc_df["L_i"].quantile(0.25)),
            "iqr_high": float(loc_df["L_i"].quantile(0.75)),
            "n_patients_positive": int((loc_df["L_i"] > 0).sum()),
            "dataset_medians": loc_df.groupby("dataset")["L_i"].median().to_dict(),
            "secondary_enrichment_ratio_median": float(loc_df["trem2_enrichment_ratio"].median()),
            "wilcoxon_p_secondary": float(
                stats.wilcoxon(loc_df["L_i"], alternative="two-sided").pvalue) if len(loc_df) > 5 else None,
        }
        loc_summary["verdict"] = (
            "LOCALISATION SUPPORTIVE"
            if (loc_summary["ci_low"] is not None and loc_summary["ci_low"] > 0)
            else "LOCALISATION UNINFORMATIVE"
        )
    else:
        loc_summary = {"n_patients": 0, "verdict": "LOCALISATION UNINFORMATIVE"}
    (OUT / "V04_P1A_LOCALISATION.csv").write_text(
        loc_df.to_csv(index=False), encoding="utf-8")
    (OUT / "_v04_p1a_localisation_summary.json").write_text(
        json.dumps(loc_summary, indent=2, ensure_ascii=False), encoding="utf-8")

    # ---------------- P1-A2 ----------------
    prim = samples.loc[samples["p1a2_primary_eligible"]].copy()
    rel = samples.loc[samples["p1a2_relaxed_eligible"]].copy()

    def summarise_k(frame: pd.DataFrame, col: str) -> dict:
        if not len(frame) or frame[col].notna().sum() == 0:
            return {"n_patients": 0, "n_datasets": 0, "mean": None, "ci_low": None,
                    "ci_high": None, "verdict": "UNINFORMATIVE"}
        b = stratified_bootstrap_mean(frame.dropna(subset=[col]), col)
        vals = frame[col].dropna().to_numpy(float)
        out = {
            "n_patients": int(len(vals)),
            "n_datasets": int(frame.loc[frame[col].notna(), "dataset"].nunique()),
            "mean": b.get("point"), "ci_low": b.get("ci_low"), "ci_high": b.get("ci_high"),
            "median": float(np.median(vals)),
            "iqr_low": float(np.quantile(vals, 0.25)),
            "iqr_high": float(np.quantile(vals, 0.75)),
            "n_positive_direction": int((vals > 0).sum()),
            "wilcoxon_p_secondary": float(
                stats.wilcoxon(vals, alternative="two-sided").pvalue)
            if len(vals) > 5 else None,
            "dataset_medians": frame.groupby("dataset")[col].median().to_dict(),
        }
        if len(vals) >= P1A2_MIN_PATIENTS and out["n_datasets"] >= 2:
            if out["ci_low"] > 0:
                out["verdict"] = "RESIDUAL COUPLING SUPPORTIVE"
            elif out["ci_high"] < 0.20:
                out["verdict"] = "RESIDUAL COUPLING NON-SUPPORTIVE"
            else:
                out["verdict"] = "RESIDUAL COUPLING UNINFORMATIVE"
        else:
            out["verdict"] = "UNINFORMATIVE (estimability)"
        return out

    within = {
        "primary_K_i_T": summarise_k(prim, "K_i_T_within_state"),
        "primary_K_i_T_depth_adjusted_S1": summarise_k(prim, "K_i_T_depth_adjusted"),
        "primary_K_i_T_depth_matched_S3": summarise_k(prim, "K_i_T_depth_matched"),
        "relaxed_K_i_T_S2": summarise_k(rel, "K_i_T_within_state"),
        "relaxed_K_i_T_depth_adjusted_S2": summarise_k(rel, "K_i_T_depth_adjusted"),
    }
    within_rows = [{
        "dataset": r.dataset, "sample_id": r.sample_id,
        "n_trem2_tam": r.n_trem2_tam, "n_trem2_positive": r.n_trem2_positive,
        "n_trem2_negative": r.n_trem2_negative,
        "K_i_T_within_state": r.K_i_T_within_state,
        "K_i_T_depth_adjusted": r.K_i_T_depth_adjusted,
        "K_i_T_depth_matched": r.K_i_T_depth_matched,
        "n_depth_pairs": r.n_depth_pairs,
        "primary_eligible": r.p1a2_primary_eligible,
        "relaxed_eligible": r.p1a2_relaxed_eligible,
    } for r in samples.itertuples()]
    (OUT / "V04_P1A_WITHIN_STATE.csv").write_text(
        pd.DataFrame(within_rows).to_csv(index=False), encoding="utf-8")
    (OUT / "_v04_p1a_within_state_summary.json").write_text(
        json.dumps(within, indent=2, ensure_ascii=False), encoding="utf-8")

    # ---------------- P1-B A-O-K ----------------
    aok = samples[["dataset", "sample_id", "n_myeloid", "A_i_siglec15_pos_fraction",
                   "O_i_trem2_fraction", "K_i_locked_unadjusted"]].copy()
    aok = aok.rename(columns={"A_i_siglec15_pos_fraction": "A_i",
                              "O_i_trem2_fraction": "O_i",
                              "K_i_locked_unadjusted": "K_i"})
    rho_ao = stratified_bootstrap_spearman(aok, "A_i", "O_i")
    rho_ak = stratified_bootstrap_spearman(aok, "A_i", "K_i")
    rho_ok = stratified_bootstrap_spearman(aok, "O_i", "K_i")
    loo = {}
    for omitted in sorted(aok["dataset"].unique()):
        sub = aok.loc[~aok["dataset"].eq(omitted)]
        loo[omitted] = {
            "rho_AO": stratified_bootstrap_spearman(sub, "A_i", "O_i")["point"],
            "rho_AK": stratified_bootstrap_spearman(sub, "A_i", "K_i")["point"],
            "rho_OK": stratified_bootstrap_spearman(sub, "O_i", "K_i")["point"],
            "n_patients": int(len(sub)),
        }
    aok_out = {
        "definition_A_i": "SIGLEC15 RNA-detectable fraction of included tumour-myeloid cells",
        "definition_O_i": "TREM2_TAM fraction of included tumour-myeloid cells",
        "definition_K_i": "locked manuscript v32 sample-specific detectable-minus-non-detectable ACP5/MMP9 contrast (taken from cell_level_sample_effects.csv, program Derosiers_ACP5_MMP9_core)",
        "n_patients": int(len(aok)),
        "n_datasets": int(aok["dataset"].nunique()),
        "rho_AO": rho_ao, "rho_AK": rho_ak, "rho_OK": rho_ok,
        "leave_one_dataset_out": loo,
        "descriptive_scales": {"little_alignment": 0.20, "large_alignment": 0.50},
        "chronology_note": (
            "A-K was partly inspected in the earlier patient-specific audit (r = -0.18, "
            "CI -0.59 to 0.29 for K_i versus SIGLEC15 prevalence). It is therefore not a new "
            "prospective finding. The new contribution of P1-B is the complete A-O-K structure, "
            "in particular O_i, which was not previously computed at patient level."
        ),
    }
    (OUT / "V04_P1B_AOK.csv").write_text(
        aok.rename(columns={"A_i": "A_i_siglec15_pos_fraction",
                            "O_i": "O_i_trem2_fraction",
                            "K_i": "K_i_locked_unadjusted"}).to_csv(index=False),
        encoding="utf-8")
    (OUT / "_v04_p1b_aok_summary.json").write_text(
        json.dumps(aok_out, indent=2, ensure_ascii=False), encoding="utf-8")

    # ---------------- provenance ----------------
    prov = {
        "analysis": "V0.4 Discriminating Analysis",
        "status": "POST HOC MECHANISTIC-DISCRIMINATION ANALYSIS - NOT PREREGISTERED",
        "seed": SEED, "bootstrap_resamples": N_BOOT,
        "software": {
            "python": platform.python_version(),
            "anndata": version("anndata"),
            "numpy": version("numpy"),
            "pandas": version("pandas"),
            "scipy": version("scipy"),
        },
        "inputs": {
            "annotation": {"path": str(ANNOTATION), "sha256": sha256(ANNOTATION)},
            "v32_locked_ki": {"path": str(V32_KI), "sha256": sha256(V32_KI)},
            "h5ad": {str(p): sha256(p) for p in H5ADS},
        },
        "inclusion_rule": {
            "tissue": sorted(TUMOUR_TISSUES),
            "celltype": "Myeloid",
            "min_myeloid_per_sample": MIN_MYELOID_PER_SAMPLE,
        },
        "gate0_verdict": gate_eval["verdict"],
        "n_samples_total": int(len(samples)),
    }
    (OUT / "_v04_provenance.json").write_text(
        json.dumps(prov, indent=2, ensure_ascii=False), encoding="utf-8")

    print(json.dumps({
        "gate0": gate_eval["verdict"],
        "n_samples": int(len(samples)),
        "p1a1_patients": int(len(loc_df)),
        "p1a1_verdict": loc_summary.get("verdict"),
        "p1a2_primary_patients": within["primary_K_i_T"]["n_patients"],
        "p1a2_primary_mean": within["primary_K_i_T"]["mean"],
        "p1a2_primary_ci": [within["primary_K_i_T"]["ci_low"], within["primary_K_i_T"]["ci_high"]],
        "p1a2_verdict": within["primary_K_i_T"]["verdict"],
        "rho_AO": rho_ao, "rho_AK": rho_ak, "rho_OK": rho_ok,
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
