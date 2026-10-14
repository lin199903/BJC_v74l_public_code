"""V0.4 influence diagnostics and Gate 0 provenance output.

The contract STOP RULES require that the inference be downgraded if extreme observations materially determine the result.
This script performs influence diagnostics only (leave-one-patient-out / leave-one-dataset-out / sparsity stratification) and introduces no new biological analysis.

It also writes the standalone file V04_STATE_ANNOTATION_PROVENANCE.csv required by the contract.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from run_V04_core import (
    OUT, SEED, N_BOOT, derive_rng, stratified_bootstrap_mean,
)

SPARSE_CUTOFF = 5  # a patient with <= 5 SIGLEC15+ TREM2 cells is treated as sparse


def boot_mean(vals: np.ndarray, datasets: np.ndarray, stream: int) -> dict:
    order = np.argsort(datasets)
    vals, datasets = vals[order], datasets[order]
    strata = [vals[datasets == d] for d in np.unique(datasets)]
    strata = [s for s in strata if len(s)]
    if not strata:
        return {"point": float("nan"), "ci_low": float("nan"), "ci_high": float("nan")}
    point = float(np.mean(np.concatenate(strata)))
    rng = derive_rng(stream)
    boot = np.empty(N_BOOT)
    for i in range(N_BOOT):
        boot[i] = np.mean(np.concatenate([rng.choice(s, len(s), replace=True) for s in strata]))
    return {"point": point,
            "ci_low": float(np.quantile(boot, 0.025)),
            "ci_high": float(np.quantile(boot, 0.975)),
            "n_patients": int(len(vals)),
            "n_datasets": int(len(strata))}


def verdict(d: dict, min_patients: int = 8) -> str:
    if d.get("n_patients", 0) < min_patients or d.get("n_datasets", 0) < 2:
        return "UNINFORMATIVE (estimability)"
    if not np.isfinite(d["ci_low"]):
        return "UNINFORMATIVE"
    if d["ci_low"] > 0:
        return "SUPPORTIVE"
    if d["ci_high"] < 0.20:
        return "NON-SUPPORTIVE"
    return "UNINFORMATIVE"


def main() -> None:
    # ---------- 1. provenance CSV ----------
    gate = json.loads((OUT / "V04_GATE_EVALUATION.json").read_text(encoding="utf-8"))
    prov_rows = []
    for r in gate["mapping"]:
        prov_rows.append({
            "dataset": r["dataset"],
            "annotation_source": "data/analysis/Task5_myeloid_subtypes/myeloid_combined_obs.csv",
            "annotation_sha256": gate["annotation_sha256"],
            "annotation_creation_date": gate["annotation_creation_date"],
            "annotation_generator": "_task5_myeloid_subtypes.py",
            "state_assignment_rule": gate["state_assignment_rule"],
            "n_parent_myeloid": r["n_parent_myeloid"],
            "n_fine_annotation": r["n_fine_annotation"],
            "n_matched": r["n_matched"],
            "mapping_fraction": r["mapping_fraction"],
            "mapping_percent": 100 * r["mapping_fraction"],
            "n_unmatched_parent": r["n_unmatched_parent"],
            "n_unmatched_annotation": r["n_unmatched_annotation"],
            "duplicate_barcodes": r["duplicate_barcodes_in_annotation"],
            "ambiguous_cell_ids": r["ambiguous_cell_ids"],
            "n_siglec15_positive": r["n_siglec15_positive"],
            "n_siglec15_negative": r["n_siglec15_negative"],
            "gate_min_mapping_fraction": gate["threshold_mapping_fraction"],
            "dataset_gate_pass": bool(r["mapping_fraction"] >= gate["threshold_mapping_fraction"]),
            "state_defined_using_acp5_mmp9_or_ki": r["state_defined_using_acp5_mmp9_or_ki"],
            "siglec15_acp5_mmp9_in_marker_union": False,
            "n_marker_genes_union": gate["n_marker_genes"],
            "exposure_layer": gate["expression_layer_for_exposure_in_v04"],
            "h5ad_sha256": r["h5ad_sha256"],
        })
    pd.DataFrame(prov_rows).to_csv(OUT / "V04_STATE_ANNOTATION_PROVENANCE.csv",
                                   index=False, encoding="utf-8")

    # ---------- 2. Influence diagnostics ----------
    samples = pd.read_csv(OUT / "_v04_sample_metrics.csv")
    rep: dict = {"note": (
        "Post hoc influence diagnostics executed to satisfy the contract STOP RULES. "
        "They are diagnostics, not new biological analyses, and were not locked before results."
    )}

    # --- P1-A1 ---
    loc = samples.loc[samples["p1a1_eligible"]].copy()
    ds = loc["dataset"].to_numpy()
    base_a1 = boot_mean(loc["L_i_localisation"].to_numpy(float), ds, stream=11)
    base_a1["verdict"] = verdict(base_a1)
    rep["p1a1_full"] = base_a1

    worst = loc["L_i_localisation"].abs().idxmax()
    ex = loc.drop(index=worst)
    d = boot_mean(ex["L_i_localisation"].to_numpy(float), ex["dataset"].to_numpy(), stream=12)
    d["verdict"] = verdict(d)
    d["excluded"] = f"{loc.loc[worst, 'dataset']}/{loc.loc[worst, 'sample_id']}"
    d["excluded_value"] = float(loc.loc[worst, "L_i_localisation"])
    rep["p1a1_without_extreme"] = d

    loo_ds = {}
    for omitted in sorted(loc["dataset"].unique()):
        sub = loc.loc[~loc["dataset"].eq(omitted)]
        dd = boot_mean(sub["L_i_localisation"].to_numpy(float),
                       sub["dataset"].to_numpy(), stream=13)
        dd["verdict"] = verdict(dd)
        loo_ds[omitted] = dd
    rep["p1a1_leave_one_dataset_out"] = loo_ds

    # Leave-one-patient-out: change in the mean and the lower CI bound
    loo_pt = []
    for idx in loc.index:
        sub = loc.drop(index=idx)
        dd = boot_mean(sub["L_i_localisation"].to_numpy(float),
                       sub["dataset"].to_numpy(), stream=14)
        loo_pt.append({
            "removed": f"{loc.loc[idx, 'dataset']}/{loc.loc[idx, 'sample_id']}",
            "removed_value": float(loc.loc[idx, "L_i_localisation"]),
            "mean_without": dd["point"], "ci_low": dd["ci_low"], "ci_high": dd["ci_high"],
            "verdict": verdict(dd),
        })
    loo_pt.sort(key=lambda r: r["mean_without"])
    rep["p1a1_leave_one_patient_out"] = {
        "min_mean_case": loo_pt[0], "max_mean_case": loo_pt[-1],
        "n_cases_flipping_verdict": int(sum(r["verdict"] != "SUPPORTIVE" for r in loo_pt)),
        "all_cases": loo_pt,
    }

    # --- P1-A2 ---
    ws = samples.loc[samples["p1a2_primary_eligible"]].copy()
    for col, key in [("K_i_T_within_state", "p1a2_primary"),
                     ("K_i_T_depth_adjusted", "p1a2_primary_depth_adjusted_S1"),
                     ("K_i_T_depth_matched", "p1a2_primary_depth_matched_S3")]:
        b = boot_mean(ws[col].to_numpy(float), ws["dataset"].to_numpy(), stream=21)
        b["verdict"] = verdict(b)
        rep[key + "_full"] = b

        w = ws[col].abs().idxmax()
        exd = ws.drop(index=w)
        dd = boot_mean(exd[col].to_numpy(float), exd["dataset"].to_numpy(), stream=22)
        dd["verdict"] = verdict(dd, min_patients=7)
        dd["excluded"] = f"{ws.loc[w, 'dataset']}/{ws.loc[w, 'sample_id']}"
        dd["excluded_value"] = float(ws.loc[w, col])
        rep[key + "_without_extreme"] = dd

        loo = {}
        for omitted in sorted(ws["dataset"].unique()):
            sub = ws.loc[~ws["dataset"].eq(omitted)]
            ddd = boot_mean(sub[col].to_numpy(float), sub["dataset"].to_numpy(), stream=23)
            ddd["verdict"] = verdict(ddd, min_patients=6)
            loo[omitted] = ddd
        rep[key + "_leave_one_dataset_out"] = loo

    # Sparsity stratification (primary K_i_T analysis only)
    sparse = ws.loc[ws["n_trem2_positive"] <= SPARSE_CUTOFF]
    dense = ws.loc[ws["n_trem2_positive"] > SPARSE_CUTOFF]
    rep["p1a2_sparsity_strata"] = {
        "cutoff_positive_cells": SPARSE_CUTOFF,
        "sparse": {**boot_mean(sparse["K_i_T_within_state"].to_numpy(float),
                               sparse["dataset"].to_numpy(), stream=31),
                   "samples": [f"{r.dataset}/{r.sample_id}" for r in sparse.itertuples()]},
        "dense": {**boot_mean(dense["K_i_T_within_state"].to_numpy(float),
                              dense["dataset"].to_numpy(), stream=32),
                  "samples": [f"{r.dataset}/{r.sample_id}" for r in dense.itertuples()]},
        "n_sparse": int(len(sparse)), "n_dense": int(len(dense)),
    }
    for k in ("sparse", "dense"):
        rep["p1a2_sparsity_strata"][k]["verdict"] = verdict(
            rep["p1a2_sparsity_strata"][k], min_patients=4)

    # --- Leave-one-patient-out (K_i_T) ---
    loo_k = []
    for idx in ws.index:
        sub = ws.drop(index=idx)
        dd = boot_mean(sub["K_i_T_within_state"].to_numpy(float),
                       sub["dataset"].to_numpy(), stream=24)
        loo_k.append({
            "removed": f"{ws.loc[idx, 'dataset']}/{ws.loc[idx, 'sample_id']}",
            "removed_value": float(ws.loc[idx, "K_i_T_within_state"]),
            "mean_without": dd["point"], "ci_low": dd["ci_low"], "ci_high": dd["ci_high"],
            "verdict": verdict(dd, min_patients=7),
        })
    loo_k.sort(key=lambda r: r["mean_without"])
    rep["p1a2_leave_one_patient_out"] = {
        "min_mean_case": loo_k[0], "max_mean_case": loo_k[-1],
        "n_cases_flipping_verdict": int(sum(r["verdict"] != "SUPPORTIVE" for r in loo_k)),
        "all_cases": loo_k,
    }

    # --- P1-B influence ---
    aok = pd.read_csv(OUT / "V04_P1B_AOK.csv")
    rep["p1b_leave_one_patient_out"] = {}
    for xc, yc, lbl in [("A_i_siglec15_pos_fraction", "O_i_trem2_fraction", "rho_AO"),
                        ("A_i_siglec15_pos_fraction", "K_i_locked_unadjusted", "rho_AK"),
                        ("O_i_trem2_fraction", "K_i_locked_unadjusted", "rho_OK")]:
        sub = aok.dropna(subset=[xc, yc])
        base = float(stats.spearmanr(sub[xc], sub[yc]).statistic)
        deltas = []
        for idx in sub.index:
            s2 = sub.drop(index=idx)
            if len(s2) > 3:
                r2 = float(stats.spearmanr(s2[xc], s2[yc]).statistic)
                deltas.append({"removed": f"{sub.loc[idx, 'dataset']}/{sub.loc[idx, 'sample_id']}",
                               "rho_without": r2, "delta": r2 - base})
        deltas.sort(key=lambda r: r["delta"])
        rep["p1b_leave_one_patient_out"][lbl] = {
            "point": base,
            "most_influential": deltas[0] if deltas else None,
            "least_influential": deltas[-1] if deltas else None,
            "max_abs_delta": float(max(abs(d["delta"]) for d in deltas)) if deltas else None,
            "n_flipping_sign": int(sum(np.sign(d["rho_without"]) != np.sign(base) for d in deltas)),
        }

    rep["seed"] = SEED
    rep["n_bootstrap"] = N_BOOT
    (OUT / "_v04_influence_diagnostics.json").write_text(
        json.dumps(rep, indent=2, ensure_ascii=False), encoding="utf-8")

    print(json.dumps({
        "provenance_csv": "written",
        "p1a1_full": rep["p1a1_full"],
        "p1a1_without_extreme": rep["p1a1_without_extreme"],
        "p1a1_flip": rep["p1a1_leave_one_patient_out"]["n_cases_flipping_verdict"],
        "p1a2_full": rep["p1a2_primary_full"],
        "p1a2_without_extreme": rep["p1a2_primary_without_extreme"],
        "p1a2_flip": rep["p1a2_leave_one_patient_out"]["n_cases_flipping_verdict"],
        "sparsity": {k: {kk: vv for kk, vv in v.items() if kk != "samples"}
                     for k, v in rep["p1a2_sparsity_strata"].items()
                     if isinstance(v, dict)},
        "p1b": rep["p1b_leave_one_patient_out"],
    }, indent=2, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
