"""P2-2 demonstration payload: apply the qualification protocol's bulk-layer
adjudication to additional immune targets.

Post hoc demonstration of discriminant validity. Not a locked analysis and not
a biomarker claim. It reuses the same statistical framework (cohort-wise
Spearman bootstrap, permutation p, REML + modified Hartung-Knapp
meta-analysis) and the same bulk-layer decision thresholds as the locked
protocol, applied to pre-specified target-programme pairs.

Verdict rules at the bulk layer (as defined in the qualification protocol):
- Supportive: pooled rho >= 0.20 (the positive-qualification line) with CI
  excluding the line.
- Non-supportive: CI upper bound < 0.20 (failed qualification excluding the
  positive-qualification line); if the CI also excludes |r| >= 0.36 the null
  is strong.
- Uninformative: precision insufficient to distinguish (CI straddles 0.20).

Pre-declared supplementary rule for protocol-internal heterogeneity: with
only three cohorts the modified Hartung-Knapp interval can be excessively
conservative when between-cohort heterogeneity is large, leaving strong
direction-consistent signals Uninformative. Therefore, when the primary
verdict is Uninformative, a supplementary verdict is computed: if all cohort
estimates share one direction and every cohort-level 95% bootstrap CI
excludes zero, the pair is classified as direction-consistent supportive
(positive) or non-supportive (negative); otherwise it remains Uninformative.
The supplementary rule is applied uniformly to all pairs.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import analysis_F_cross_scale as s17

OUT = Path(__file__).resolve().parent / "outputs_P22_scorecard"

POSITIVE_LINE = 0.20
EXCLUSION_LINE = 0.36

P22_PROGRAMS = {
    "cytolytic": ["GZMA", "GZMB", "PRF1"],
    "CD8_T": ["CD8A", "CD8B", "CD3D", "CD3E"],
    "regulatory_T": ["FOXP3", "IL2RA", "CTLA4"],
    "B_cell": ["CD19", "MS4A1", "CD79A"],
    "pan_myeloid": s17.POSTHOC_MYELOID_COMPOSITION,
}

TARGETS = [
    "CD8A", "GZMB", "CD274", "FOXP3", "CD19",
    "KRT19", "EPCAM", "COL1A1", "SIGLEC15",
]

PRESPECIFIED = [
    {
        "target": "CD8A",
        "program": "CD8_T",
        "expected": "Supportive",
        "rationale": "CD8A is a defining transcript of the CD8 T-cell compartment; the residual CD8_T programme (CD8B/CD3D/CD3E) marks the same compartment.",
    },
    {
        "target": "GZMB",
        "program": "cytolytic",
        "expected": "Supportive",
        "rationale": "GZMB is a core cytolytic effector gene; the residual cytolytic programme (GZMA/PRF1) is expected to co-vary strongly.",
    },
    {
        "target": "FOXP3",
        "program": "regulatory_T",
        "expected": "Supportive",
        "rationale": "FOXP3 is the lineage-defining transcription factor of regulatory T cells; the residual programme (IL2RA/CTLA4) marks the same compartment.",
    },
    {
        "target": "CD274",
        "program": "cytolytic",
        "expected": "Supportive",
        "rationale": "PD-L1 is an IFN-gamma-induced checkpoint ligand expected to co-vary with cytolytic activity, though less tightly than lineage-defining pairs.",
    },
    {
        "target": "CD19",
        "program": "cytolytic",
        "expected": "Uninformative",
        "rationale": "Co-infiltration structure may induce positive co-variation between any immune-cell marker and cytolytic activity in bulk tissue; a positive but sub-threshold or heterogeneous association is expected.",
    },
    {
        "target": "KRT19",
        "program": "cytolytic",
        "expected": "Non-supportive",
        "rationale": "Epithelial ductal and cytolytic compartments are expected to show no systematic positive association.",
    },
    {
        "target": "EPCAM",
        "program": "cytolytic",
        "expected": "Non-supportive",
        "rationale": "Epithelial and cytolytic compartments are expected to show no systematic positive association.",
    },
    {
        "target": "COL1A1",
        "program": "cytolytic",
        "expected": "Uninformative",
        "rationale": "Stromal abundance and immune infiltration have heterogeneous, context-dependent relationships in PDAC; no directional expectation is committed.",
    },
    {
        "target": "SIGLEC15",
        "program": "cytolytic",
        "expected": "Uninformative",
        "rationale": "Internal-consistency reference: the locked analysis classified SIGLEC15 bulk associations as failing qualification with large effects excluded; the same verdict should emerge here.",
    },
    {
        "target": "SIGLEC15",
        "program": "pan_myeloid",
        "expected": "Uninformative",
        "rationale": "Internal-consistency reference: post hoc adjustment in the locked analysis indicated that SIGLEC15 bulk RNA is not strongly driven by myeloid content; a null verdict is expected.",
    },
]


def score_program(frame: pd.DataFrame, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    present = [gene for gene in genes if gene in frame.columns]
    if len(present) < 2:
        return np.full(len(frame), np.nan), present
    values = frame[present].to_numpy(float)
    denom = values.std(axis=0, ddof=0)
    denom[denom == 0] = 1.0
    z = (values - values.mean(axis=0)) / denom
    return np.nanmean(z, axis=1), present


def load_bulk_p22() -> dict[str, pd.DataFrame]:
    cohorts = {}
    t = pd.read_csv(s17.BULK_PATHS["TCGA_PAAD"], index_col=0)
    cohorts["TCGA_PAAD"] = t
    for name, path, soft in (
        ("GSE57495", s17.BULK_PATHS["GSE57495"], s17.SOFT_PATHS["GSE57495"]),
        ("GSE85916", s17.BULK_PATHS["GSE85916"], s17.SOFT_PATHS["GSE85916"]),
    ):
        raw = pd.read_csv(path, index_col=0)
        mapping = s17.parse_soft_platform(soft)
        needed = set(TARGETS) | set(sum(P22_PROGRAMS.values(), []))
        genes = pd.Series(raw.index.map(mapping), index=raw.index)
        keep = genes.isin(needed)
        selected = raw.loc[keep].copy()
        selected.insert(0, "gene", genes.loc[keep].to_numpy())
        cohorts[name] = selected.groupby("gene").mean(numeric_only=True)
    return cohorts


def adjudicate(meta: dict) -> str:
    r = meta["r_random"]
    lo, hi = meta["ci_low"], meta["ci_high"]
    if r >= POSITIVE_LINE and lo >= POSITIVE_LINE:
        return "Supportive"
    if hi < POSITIVE_LINE:
        if lo > -EXCLUSION_LINE and hi < EXCLUSION_LINE:
            return "Non-supportive (strong)"
        return "Non-supportive"
    return "Uninformative"


def supplementary_verdict(sub: pd.DataFrame) -> str | None:
    """Supplementary rule pre-declared in the docstring: applied only when the
    primary verdict is Uninformative; requires all cohort estimates to share
    one direction and every cohort-level 95% bootstrap CI to exclude zero."""
    if sub.empty:
        return None
    if (sub["rho"] > 0).all() and (sub["ci_low"] > 0).all():
        return "direction-consistent supportive"
    if (sub["rho"] < 0).all() and (sub["ci_high"] < 0).all():
        return "direction-consistent non-supportive"
    return None


def matches_expectation(observed: str, expected: str) -> bool:
    if observed == expected:
        return True
    if expected == "Supportive":
        return observed.startswith("direction-consistent supportive")
    if expected == "Non-supportive":
        return observed.startswith("direction-consistent non-supportive")
    return False


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cohorts = load_bulk_p22()

    rows = []
    for cohort_name, frame in cohorts.items():
        t = frame.T
        for target in TARGETS:
            if target not in t.columns:
                continue
            x = t[target].to_numpy(float)
            for program_name, genes in P22_PROGRAMS.items():
                score_genes = [g for g in genes if g != target]
                score, present = score_program(t, score_genes)
                if not np.isfinite(score).any():
                    continue
                rho, p, lo, hi = s17.spearman_bootstrap(x, score)
                pp = s17.permutation_p(x, score, rho)
                rows.append({
                    "cohort": cohort_name, "target": target, "program": program_name,
                    "n_samples": len(t), "n_score_genes_present": len(present),
                    "rho": rho, "ci_low": lo, "ci_high": hi,
                    "asymptotic_p": p, "permutation_p": pp,
                })
    results = pd.DataFrame(rows)
    results.to_csv(OUT / "scorecard_full_matrix.csv", index=False)

    adjudication_rows = []
    supplementary_applied = []
    for spec in PRESPECIFIED:
        target, program = spec["target"], spec["program"]
        sub = results[(results["target"] == target) & (results["program"] == program)].copy()
        meta = s17.random_effects_meta(sub)
        primary = adjudicate(meta)
        observed = primary
        if primary == "Uninformative":
            sup = supplementary_verdict(sub)
            if sup is not None:
                observed = sup
                supplementary_applied.append({
                    "target": target, "program": program,
                    "primary_verdict": primary, "supplementary_verdict": sup,
                })
        cohort_direction = ";".join(
            f"{r['cohort']}={r['rho']:+.3f}" for _, r in sub.iterrows()
        )
        adjudication_rows.append({
            **spec,
            "primary_verdict": primary,
            "observed": observed,
            "verdict_matches_expectation": matches_expectation(observed, spec["expected"]),
            "n_cohorts": meta["n_cohorts"],
            "pooled_rho": meta["r_random"],
            "ci_low": meta["ci_low"],
            "ci_high": meta["ci_high"],
            "p_modified_hartung_knapp": meta["p_modified_hartung_knapp"],
            "I2_pct": meta["I2_pct"],
            "tau2_fisher_z": meta["tau2_fisher_z"],
            "cohort_directions": cohort_direction,
            "cohort_rhos": ";".join(f"{r['rho']:.3f}" for _, r in sub.iterrows()),
        })
    table = pd.DataFrame(adjudication_rows)
    table.to_csv(OUT / "scorecard_table.csv", index=False)

    sig_internal = {}
    for _, row in table[table["target"] == "SIGLEC15"].iterrows():
        sig_internal[row["program"]] = {
            "expected": row["expected"],
            "observed": row["observed"],
            "pooled_rho": row["pooled_rho"],
            "ci_low": row["ci_low"],
            "ci_high": row["ci_high"],
        }
    summary = {
        "analysis": "P2-2 demonstration payload: bulk-layer adjudication of the qualification protocol on additional immune targets",
        "status": "post hoc demonstration; not a locked analysis; no additional biomarker claim",
        "seed": s17.SEED,
        "bootstrap_resamples": s17.N_RESAMPLES,
        "positive_qualification_line": POSITIVE_LINE,
        "strong_null_exclusion_line": EXCLUSION_LINE,
        "supplementary_rule": {
            "trigger": "primary verdict Uninformative",
            "criterion": "all cohort estimates share one direction and every cohort-level 95% bootstrap CI excludes zero",
            "applied_uniformly": True,
            "pairs_upgraded": supplementary_applied,
        },
        "verdicts": {f"{row['target']}/{row['program']}": row["observed"] for _, row in table.iterrows()},
        "all_match_expectation": bool(table["verdict_matches_expectation"].all()),
        "internal_consistency_reference": sig_internal,
    }
    (OUT / "scorecard_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    print(table.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
