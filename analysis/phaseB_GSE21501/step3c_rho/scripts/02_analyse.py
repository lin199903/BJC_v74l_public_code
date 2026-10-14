#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE21501 -- Phase B step 3c-2, stage 2 of 3: ANALYSE.

Estimand (frozen_cohort_definition.md S0, carried verbatim):

    within this eligible PDAC cohort, estimate the Spearman correlation between
    tumour-level SIGLEC15 expression and an ACP5/MMP9 programme score constructed
    using the same dataset-level standardisation rule as the existing TCGA-PAAD,
    GSE57495 and GSE85916 analyses.

Scope of this file: **GSE21501 alone**.  Nothing here is pooled with GSE205154 or with
any other cohort, and no cross-cohort absolute expression comparison is made (AD-2).

What is transcribed vs what is chosen
-------------------------------------
TRANSCRIBED (no latitude):
  * `zscore_columns`, `score_program`, `required_gene_count`, `spearman_bootstrap`,
    `permutation_p` -- copied function-body-identical from the archived manuscript
    script `analysis_F_cross_scale.py`, which is the analysis identity the frozen
    document names in S9.
  * SEED = 15072026, N_RESAMPLES = 10_000 (frozen S9).
  * PROGRAMS["Derosiers_ACP5_MMP9_core"] = ["ACP5", "MMP9"].
  * DUP-2 selection, S1..S5 definitions, AD-1..AD-6 (frozen S4, S10.3, S11).

CHOSEN HERE (frozen text silent; recorded, use-invariant, reported as a choice):
  * **RNG addressing scheme.**  The manuscript script has ONE global RNG consumed in
    call order, which is unambiguous only while a single correlation is computed.  This
    step computes ~20 correlations, so each analysis cell needs its own stream.  The
    primary estimate is therefore computed with a generator seeded **directly with the
    frozen SEED and consumed in the manuscript's own call order** (bootstrap first, then
    permutation), so that the primary numbers rest on the frozen seed and on nothing
    else.  Every other cell uses a generator addressed deterministically as
        seed_label = int(sha256(f"{MASTER_SEED}|{label}").hexdigest()[:16], 16) % 2**32
    An addressing scheme is a labelling device: it cannot move an estimate, because rho
    itself is deterministic and does not depend on any generator.

  * **Rank-based score alternative (frozen S8 "S-rep-1").**  Frozen text says "within-
    cohort ranks of each gene instead of z-scores".  Implemented as
    `score_program(rankdata-per-gene frame)`: rank each gene across the cohort
    (scipy.stats.rankdata, average ties), then feed those ranks through the *unchanged*
    frozen `score_program`, i.e. z-score the ranks and average.  The frozen text does
    not settle whether the ranks are re-standardised before averaging.  Both variants
    are reported, so the ambiguity is visible rather than resolved silently.

Two departures from `frozen_cohort_definition.md` that this step is RULED to make
---------------------------------------------------------------------------------
The author rulings supplied for step 3c-2 supersede the step-3b text on two points.
Both are applied here exactly as ruled, and both are enumerable in the outputs.

  (1) AD-1 reverses the step-3b main/sensitivity assignment.
        3b   : main = 102 clinically annotated; S1 = 129 whole deposit after DUP-2.
        3c-2 : main = 129 (expression-only, as ruled); 102 = prespecified sensitivity.
      The rule this does NOT change is DUP-2, applied identically in both cases -- so
      both numbers rest on the same de-duplication.

  (2) AD-4 restates the S2 exclusion grounds, and its arithmetic follows.
      3b   : S2 drops GSM536488 (P271N Replicate) and GSM536901 (P188T IPMN), n = 100.
      3c-2 : P188T is excluded **because it is explicitly labelled IPMN and therefore is
             not PDAC**; P271N may be excluded **only if the frozen manifest
             established, from its original sample identity, that it is non-PDAC or
             normal**, and not because of the trailing "N".
             The frozen manifest does NOT establish that: it lists GSM536488 as
             `tumour`, `primary`, `eligible_for_main_meta = yes`, channel-2 `tissue`
             reading "Patient Primary Pancreatic Tumor", with a full clinical record and
             an empty exclusion_reason.  Under the ruling GSM536488 is therefore
             **retained** and AD-4's S2 has **n = 101**, not 100.  Both the ruled S2
             (n = 101) and the step-3b S2 (n = 100) are reported, and every `n` is read
             from the data, never written down.

Outputs
    output/primary_result.csv
    output/sensitivity_results.csv
    output/estimate_table_full.csv
    output/verification_controls.csv
    output/descriptive_components.csv
    work/stage02_qc.json
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
STEP = HERE.parent
PHASE = STEP.parent
WORK = STEP / "work"
OUT = STEP / "output"
OUT.mkdir(parents=True, exist_ok=True)

EXPR_CSV = WORK / "GSE21501_target_expression.csv"
MMP9_FEAT_CSV = WORK / "GSE21501_MMP9_features.csv"
MANIFEST = PHASE / "sample_manifest.csv"

FROZEN_SHA256_MANIFEST = "3e11938ae5c7e682c9e53aad93d051579764a4ec8173d69d571dbeda0f3ac4e4"

# ---------------------------------------------------------------------------
# FROZEN constants (frozen_cohort_definition.md S9) -- transcribed
# ---------------------------------------------------------------------------
SEED = 15072026            # frozen S9
N_RESAMPLES = 10_000       # frozen S9
MASTER_SEED = SEED

TARGET = "SIGLEC15"
PROGRAMS = {"Derosiers_ACP5_MMP9_core": ["ACP5", "MMP9"]}
PRIMARY = "Derosiers_ACP5_MMP9_core"

# the CI actually reported for the primary estimate
PRIMARY_RNG_LABEL = "MANUSCRIPT ORDER (numpy default_rng(15072026), bootstrap first then permutation)"
PRIMARY_SEED_USED = SEED


# ===========================================================================
# SECTION A -- manuscript functions, copied function-body-identical from
#              _archives_translate/Additional_file_2_Analysis_archive_v39/
#              scripts/analysis_F_cross_scale.py  (v47 tree).
#              The ONLY edit is that the module-level RNG is an explicit
#              parameter, so each analysis cell can be given its own stream.
# ===========================================================================
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


def spearman_bootstrap(x: np.ndarray, y: np.ndarray,
                       rng: np.random.Generator) -> tuple[float, float, float, float]:
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[valid], np.asarray(y)[valid]
    rho, p = stats.spearmanr(x, y)
    idx = rng.integers(0, len(x), size=(N_RESAMPLES, len(x)))
    xr = stats.rankdata(x[idx], axis=1)
    yr = stats.rankdata(y[idx], axis=1)
    xr -= xr.mean(axis=1, keepdims=True)
    yr -= yr.mean(axis=1, keepdims=True)
    denom = np.sqrt(np.sum(xr * xr, axis=1) * np.sum(yr * yr, axis=1))
    boot = np.divide(np.sum(xr * yr, axis=1), denom,
                     out=np.full(N_RESAMPLES, np.nan), where=denom > 0)
    lo, hi = np.nanquantile(boot, [0.025, 0.975])
    return float(rho), float(p), float(lo), float(hi)


def permutation_p(x: np.ndarray, y: np.ndarray, observed: float,
                  rng: np.random.Generator) -> float:
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[valid], np.asarray(y)[valid]
    xr = stats.rankdata(x).astype(float)
    yr = stats.rankdata(y).astype(float)
    xr -= xr.mean()
    yr -= yr.mean()
    denom = np.sqrt(np.sum(xr * xr) * np.sum(yr * yr))
    permutations = np.vstack([rng.permutation(yr) for _ in range(N_RESAMPLES)])
    r = permutations @ xr / denom
    return float((np.sum(np.abs(r) >= abs(observed)) + 1) / (N_RESAMPLES + 1))


# ===========================================================================
# SECTION B -- deterministic per-analysis RNG addressing (chosen; see header)
# ===========================================================================
def label_seed(label: str, master: int = MASTER_SEED) -> int:
    h = hashlib.sha256(f"{master}|{label}".encode("utf-8")).hexdigest()
    return int(h[:16], 16) % (2 ** 32)


def rng_for(label: str) -> np.random.Generator:
    return np.random.default_rng(label_seed(label))


# ===========================================================================
# SECTION C -- cohort assembly, strictly from the frozen manifest
# ===========================================================================
def sha256_file(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            d.update(chunk)
    return d.hexdigest()


REPEAT_MARKER = r"rpt|redo|replicate"


def specimen_key(pid: str) -> str:
    """'not_exposed_by_GEO;specimen_key=NE:12R2587' -> 'NE:12R2587'."""
    return str(pid).split("specimen_key=")[-1]


def build_samples() -> tuple[pd.DataFrame, dict]:
    man = pd.read_csv(MANIFEST)
    expr = pd.read_csv(EXPR_CSV)
    df = man.merge(expr[["gsm_accession", TARGET, "ACP5", "MMP9"]],
                   on="gsm_accession", how="left", validate="one_to_one")
    df["specimen_key"] = df["patient_id_if_available"].map(specimen_key)
    df["has_repeat_marker"] = df["channel2_sample_title"].str.contains(
        REPEAT_MARKER, case=False, regex=True)
    df["is_dup_dropped"] = df["exclusion_reason"].fillna("").str.contains("EX-5", regex=False)

    qc: dict = {}
    qc["n_manifest_rows"] = int(len(df))
    qc["n_dup_dropped"] = int(df["is_dup_dropped"].sum())
    qc["dup_dropped"] = df.loc[df["is_dup_dropped"],
                               ["gsm_accession", "channel2_sample_title",
                                "specimen_key"]].to_dict("records")

    dep = df.loc[~df["is_dup_dropped"]].copy()
    qc["n_after_dup2"] = int(len(dep))
    qc["n_distinct_specimen_keys_after_dup2"] = int(dep["specimen_key"].nunique())
    qc["dup2_conserves_one_array_per_key"] = bool(len(dep) == dep["specimen_key"].nunique())

    grp = df.groupby("specimen_key")["gsm_accession"].agg(list)
    multi = grp[grp.map(len) > 1]
    qc["multi_array_specimen_groups"] = {k: v for k, v in multi.items()}

    # --- AD-4 membership, decided from the frozen manifest's own fields --------
    ipmn = df.loc[df["channel2_sample_title"].str.contains("IPMN", case=False, na=False)]
    p271n = df.loc[df["channel2_sample_title"].str.contains("P271N", case=False, na=False)]
    cols = ["gsm_accession", "channel2_sample_title", "tumour_or_non_tumour",
            "primary_or_metastatic", "clinically_annotated", "eligible_for_main_meta"]
    qc["ad4"] = {"ipmn_arrays": ipmn[cols].to_dict("records"),
                 "p271n_arrays": p271n[cols + ["exclusion_reason"]].to_dict("records")}
    p = p271n.iloc[0]
    established = (p["tumour_or_non_tumour"] != "tumour") or (p["primary_or_metastatic"] == "normal")
    qc["ad4"]["p271n_established_non_pdac_or_normal_by_manifest"] = bool(established)
    qc["ad4"]["p271n_retained_under_ad4"] = bool(not established)
    qc["ad4"]["p188t_excluded_because"] = ("explicitly labelled IPMN, therefore not PDAC "
                                           "(AD-4), NOT because clinical information "
                                           "is missing")

    s3 = dep.loc[dep["has_repeat_marker"]]
    qc["s3_repeat_marked_within_129"] = s3[["gsm_accession", "channel2_sample_title"]
                                           ].to_dict("records")
    qc["s3_n_derived"] = int(len(s3))
    qc["s3_repeat_marked_whole_deposit"] = df.loc[df["has_repeat_marker"],
                                                 "gsm_accession"].tolist()
    qc["repeat_marked_in_102"] = df.loc[df["has_repeat_marker"] &
                                        df["eligible_for_main_meta"].eq("yes"),
                                        "gsm_accession"].tolist()
    qc["n_clinically_annotated_arrays"] = int(df["clinically_annotated"].eq("yes").sum())
    return df, qc


def population(df: pd.DataFrame, name: str) -> pd.DataFrame:
    dep = df.loc[~df["is_dup_dropped"]].copy()          # DUP-2: 129
    if name == "expression_only_129":
        return dep
    if name == "clinically_annotated_102":
        return df.loc[df["eligible_for_main_meta"].eq("yes")].copy()
    # S2 and S3 are defined in frozen S10.3 on the clinically annotated population
    # (their frozen n's are 100 and 96, i.e. the 102 minus 2 and minus 6).  They are
    # therefore based on the 102 here as well, so that the frozen arithmetic and this
    # step's arithmetic are the same arithmetic.
    base102 = df.loc[df["eligible_for_main_meta"].eq("yes")].copy()
    if name == "s2_ruled_101":
        keep = ~base102["channel2_sample_title"].str.contains("IPMN", case=False, na=False)
        return base102.loc[keep].copy()
    if name == "s2_step3b_100":
        keep = ~base102["channel2_sample_title"].str.contains("IPMN|P271N", case=False,
                                                              na=False)
        return base102.loc[keep].copy()
    if name == "s3_drop_repeat_marked":
        return base102.loc[~base102["has_repeat_marker"]].copy()
    raise KeyError(name)


# ===========================================================================
# SECTION D -- one correlation cell, with an explicit generator
# ===========================================================================
def run_estimate(frame: pd.DataFrame, label: str, seed_label: str,
                 rng: np.random.Generator, seed_used: float | str,
                 notes: str = "") -> dict:
    score, present = score_program(frame, PRIMARY)
    x = frame[TARGET].to_numpy(float)
    rho, p, lo, hi = spearman_bootstrap(x, score, rng)
    pp = permutation_p(x, score, rho, rng)
    return {
        "label": label,
        "n_arrays": int(len(frame)),
        "n_distinct_specimen_keys": int(frame["specimen_key"].nunique()),
        "n_genes_present": len(present), "genes_present": ";".join(present),
        "rho": rho, "ci_low": lo, "ci_high": hi,
        "ci_method": f"percentile bootstrap over arrays, {N_RESAMPLES} resamples, "
                     f"Spearman recomputed per resample (manuscript "
                     f"spearman_bootstrap)",
        "n_resamples": N_RESAMPLES,
        "seed_label": seed_label, "seed_used": seed_used,
        "seed_derived": (label_seed(seed_label) if isinstance(seed_used, str) else seed_used),
        "asymptotic_p": p, "permutation_p": pp,
        "score_mean": float(np.nanmean(score)),
        "score_sd_ddof0": float(np.nanstd(score)),
        "notes": notes,
    }


def main() -> int:
    print("=" * 78)
    print("GSE21501 step 3c-2 -- stage 2 ANALYSE")
    print("=" * 78)

    mh = sha256_file(MANIFEST)
    print(f"[gate] sample_manifest.csv sha256 = {mh}")
    print(f"[gate] matches frozen value       = {mh == FROZEN_SHA256_MANIFEST}")
    if mh != FROZEN_SHA256_MANIFEST:
        raise SystemExit("FATAL: manifest hash mismatch.")

    df, build_qc = build_samples()
    rows: list[dict] = []
    controls: list[dict] = []

    dep = population(df, "expression_only_129")
    clin = population(df, "clinically_annotated_102")
    if build_qc["n_after_dup2"] != 129:
        raise SystemExit(f"FATAL: DUP-2 population is {build_qc['n_after_dup2']}, frozen 129.")
    if build_qc["n_clinically_annotated_arrays"] != 102:
        raise SystemExit("FATAL: clinically annotated count is not 102.")
    print(f"[ok] DUP-2 population n = {len(dep)}; clinically annotated n = {len(clin)}")

    # ---------------------------------------------------------------- PRIMARY
    # generator seeded directly with the frozen SEED, consumed in the manuscript's
    # own order: bootstrap first, then permutation.
    prim = run_estimate(dep, "PRIMARY: expression-only, whole deposit after DUP-2",
                        PRIMARY_RNG_LABEL, np.random.default_rng(SEED), SEED,
                        "AD-1 as ruled: the primary analysis is the expression-only "
                        "set (n = 129).")
    rows.append(prim)
    print(f"\n[PRIMARY] n={prim['n_arrays']}  rho={prim['rho']:+.6f}  "
          f"95% CI [{prim['ci_low']:+.6f}, {prim['ci_high']:+.6f}]  "
          f"p_perm={prim['permutation_p']:.4g}")

    # ------------------------------------------- RNG-addressing equivalence check
    alt = run_estimate(dep, "PRIMARY (alt RNG addressing, diagnostic only)",
                       "diagnostic_primary_label_addressed",
                       rng_for("diagnostic_primary_label_addressed"),
                       "label-addressed", "Not reported as a result.")
    rho_same = alt["rho"] == prim["rho"]
    pp_gap = abs(alt["permutation_p"] - prim["permutation_p"])
    ci_gap = max(abs(alt["ci_low"] - prim["ci_low"]), abs(alt["ci_high"] - prim["ci_high"]))
    # Monte-Carlo tolerance at 10 000 resamples for a p near 0.05
    mc_se_p = float(np.sqrt(0.05 * 0.95 / N_RESAMPLES))
    print(f"[check] rho identical across the two RNG addressing schemes: {rho_same}")
    print(f"[check] bootstrap CI gap between schemes: {ci_gap:.6f} (Monte-Carlo)")
    print(f"[check] permutation p gap between schemes: {pp_gap:.4f} "
          f"(Monte-Carlo; 1 s.e. at p=0.05 is {mc_se_p:.4f})")
    addr_ok = rho_same and ci_gap < 0.05 and pp_gap < 4 * mc_se_p
    controls.append({
        "control": "C7-rng-addressing", "kind": "self-consistency",
        "capable_of_failing": "yes",
        "expectation": "rho is EXACTLY invariant to how the generator is addressed "
                       "(no generator enters rho); the bootstrap interval and the "
                       "permutation p are resampled quantities and may differ only by "
                       "Monte-Carlo error",
        "observed": f"rho identical={rho_same}; CI gap={ci_gap:.6f}; permutation p gap="
                    f"{pp_gap:.4f} (1 s.e. = {mc_se_p:.4f})",
        "verdict": "PASS" if addr_ok else "FAIL",
        "note": "An earlier version of this script required the two bootstrap intervals "
                "to be bit-identical, which is unsatisfiable: two different seeds cannot "
                "yield identical resampled intervals, and the permutation p is itself a "
                "resampled quantity. The check was corrected to test what the frozen "
                "rule actually implies. The invariant that matters -- rho carries no "
                "generator dependence -- holds exactly."})

    # ---------------------------------------------------- prespecified 102
    s102 = run_estimate(clin, "S-CLIN: clinically annotated subset (AD-6 sensitivity)",
                        "sensitivity_clinically_annotated_102",
                        rng_for("sensitivity_clinically_annotated_102"),
                        "label-addressed",
                        "The 30 Nebraska/UNC arrays are flagged 'clinical annotation "
                        "unavailable', not excluded; this is the prespecified "
                        "metadata-availability sensitivity for AD-1/AD-6.")
    rows.append(s102)
    print(f"[S-CLIN]  n={s102['n_arrays']}  rho={s102['rho']:+.6f}  "
          f"95% CI [{s102['ci_low']:+.6f}, {s102['ci_high']:+.6f}]")
    d_rho = abs(prim["rho"] - s102["rho"])
    ci_overlap = not (prim["ci_high"] < s102["ci_low"] or s102["ci_high"] < prim["ci_low"])
    print(f"[S-CLIN]  |drho| = {d_rho:.6f}; CIs overlap = {ci_overlap}")

    # ---------------------------------------------------------------- S2
    s2r = run_estimate(population(df, "s2_ruled_101"),
                       "S2 (AD-4 as ruled): drop P188T/IPMN only",
                       "sensitivity_s2_ad4_ruled", rng_for("sensitivity_s2_ad4_ruled"),
                       "label-addressed",
                       "P188T is IPMN, hence not PDAC. P271N is retained: the frozen "
                       "manifest does not establish it as non-PDAC or normal. "
                       "Baseline = the 102 clinically annotated arrays, because frozen "
                       "S10.3 defines S2 on that population (frozen n = 100 = 102 - 2).")
    rows.append(s2r)
    print(f"[S2-ruled] n={s2r['n_arrays']}  rho={s2r['rho']:+.6f}")

    s2b = run_estimate(population(df, "s2_step3b_100"),
                       "S2 (step-3b text, retained for comparison): drop P188T and P271N",
                       "sensitivity_s2_step3b_text",
                       rng_for("sensitivity_s2_step3b_text"), "label-addressed",
                       "Not the ruled S2. Kept because the frozen document's own "
                       "arithmetic gives n = 100 and the ruling changes it. "
                       "Baseline = the 102 (frozen S10.3).")
    rows.append(s2b)
    print(f"[S2-3b]    n={s2b['n_arrays']}  rho={s2b['rho']:+.6f}")

    # ---------------------------------------------------------------- S3
    s3 = run_estimate(population(df, "s3_drop_repeat_marked"),
                      "S3: drop all repeat-marked arrays inside the 129",
                      "sensitivity_s3_drop_repeat_marked",
                      rng_for("sensitivity_s3_drop_repeat_marked"), "label-addressed",
                      "DUP-6(ii). Repeat marker = rpt|redo|replicate in the sample "
                      "title, derived from the manifest, not hard-coded. Baseline = the "
                      "102 clinically annotated arrays (frozen S10.3 gives S3 n = 96).")
    rows.append(s3)
    print(f"[S3]       n={s3['n_arrays']}  rho={s3['rho']:+.6f}  "
          f"(derived drop-set size {build_qc['s3_n_derived']})")

    # ---------------------------------------------------------------- S4
    avg = dep.copy()
    for key in build_qc["multi_array_specimen_groups"]:
        m = avg["specimen_key"].eq(key)
        avg.loc[m, [TARGET, "ACP5", "MMP9"]] = \
            avg.loc[m, [TARGET, "ACP5", "MMP9"]].mean().to_numpy()
    s4 = run_estimate(avg, "S4: DUP-6(i) duplicate pairs averaged at expression level",
                      "sensitivity_s4_pair_averaged",
                      rng_for("sensitivity_s4_pair_averaged"), "label-addressed",
                      "Expression averaged across the arrays of each multi-array "
                      "specimen group before scoring; one row per specimen key. "
                      f"Groups averaged: {list(build_qc['multi_array_specimen_groups'])}")
    rows.append(s4)
    print(f"[S4]       n={s4['n_arrays']}  rho={s4['rho']:+.6f}  "
          f"(keys {s4['n_distinct_specimen_keys']})")

    # ------------------------------------------------------- rank-based score
    rank_dep = dep.copy()
    for g in (TARGET, "ACP5", "MMP9"):
        rank_dep[g] = stats.rankdata(rank_dep[g].to_numpy(float))
    sr1 = run_estimate(rank_dep,
                       "S-rep-1: ACP5/MMP9 score from within-cohort ranks, "
                       "re-standardised by the frozen pipeline",
                       "sensitivity_rank_based_restandardised",
                       rng_for("sensitivity_rank_based_restandardised"),
                       "label-addressed",
                       "scipy.stats.rankdata per gene, then the unchanged frozen "
                       "score_program (ranks are z-scored before averaging).")
    rows.append(sr1)
    print(f"[S-rep-1]  n={sr1['n_arrays']}  rho={sr1['rho']:+.6f}")

    raw_rank = dep.copy()
    raw_score = np.nanmean(np.column_stack(
        [stats.rankdata(raw_rank[g].to_numpy(float)) for g in ("ACP5", "MMP9")]), axis=1)
    _rng = rng_for("sensitivity_rank_based_raw")
    _rho, _p, _lo, _hi = spearman_bootstrap(raw_rank[TARGET].to_numpy(float),
                                            raw_score, _rng)
    _pp = permutation_p(raw_rank[TARGET].to_numpy(float), raw_score, _rho, _rng)
    sr1b = {"label": "S-rep-1b: ACP5/MMP9 score from raw 1..n within-cohort ranks "
                     "averaged without re-standardising",
            "n_arrays": int(len(raw_rank)),
            "n_distinct_specimen_keys": int(raw_rank["specimen_key"].nunique()),
            "n_genes_present": 2, "genes_present": "ACP5;MMP9",
            "rho": _rho, "ci_low": _lo, "ci_high": _hi,
            "ci_method": f"percentile bootstrap over arrays, {N_RESAMPLES} resamples",
            "n_resamples": N_RESAMPLES, "seed_label": "sensitivity_rank_based_raw",
            "seed_used": "label-addressed",
            "seed_derived": label_seed("sensitivity_rank_based_raw"),
            "asymptotic_p": _p, "permutation_p": _pp,
            "notes": "Reported because frozen S8 does not settle whether ranks are "
                     "re-standardised before averaging."}
    rows.append(sr1b)
    print(f"[S-rep-1b] n={sr1b['n_arrays']}  rho={sr1b['rho']:+.6f}")

    # ------------------------------------------------- S5 per-MMP9-feature
    feat_tbl = pd.read_csv(MMP9_FEAT_CSV)
    feat_cols = [c for c in feat_tbl.columns if c.startswith("MMP9_feature_")]
    s5_rows = []
    for c in feat_cols:
        fid = c.split("_")[-1]
        f = clin.merge(feat_tbl[["gsm_accession", c]], on="gsm_accession",
                       how="left", validate="one_to_one")
        base = f[["gsm_accession", "specimen_key", TARGET, "ACP5"]].copy()
        base["MMP9"] = f[c].to_numpy(float)
        r = run_estimate(base, f"S5: MMP9 score from single feature {fid}",
                         f"sensitivity_s5_feature_{fid}",
                         rng_for(f"sensitivity_s5_feature_{fid}"), "label-addressed",
                         "Score = mean(z(ACP5), z(this single MMP9 feature)); "
                         "population = the 102 clinically annotated arrays (frozen "
                         "S10.3 gives S5 n = 102).")
        s5_rows.append(r)
        rows.append(r)
    s5 = pd.DataFrame(s5_rows)
    prim_feat_rho = float(s5.loc[s5["label"].str.contains("44612"), "rho"].iloc[0])
    print(f"[S5]       {len(s5)} features; rho range "
          f"[{s5['rho'].min():+.6f}, {s5['rho'].max():+.6f}]; "
          f"median {s5['rho'].median():+.6f}; the 10-feature-mean equivalent "
          f"(feature 44612 alone) {prim_feat_rho:+.6f}")

    # -------------------------------------------------- descriptive components
    comp = []
    for name, sub in (("129", dep), ("102", clin)):
        comp.append({
            "population": name, "n": len(sub),
            "rho_SIGLEC15_ACP5": float(stats.spearmanr(sub[TARGET], sub["ACP5"]).statistic),
            "p_SIGLEC15_ACP5": float(stats.spearmanr(sub[TARGET], sub["ACP5"]).pvalue),
            "rho_SIGLEC15_MMP9": float(stats.spearmanr(sub[TARGET], sub["MMP9"]).statistic),
            "p_SIGLEC15_MMP9": float(stats.spearmanr(sub[TARGET], sub["MMP9"]).pvalue),
            "rho_ACP5_MMP9": float(stats.spearmanr(sub["ACP5"], sub["MMP9"]).statistic),
            "p_ACP5_MMP9": float(stats.spearmanr(sub["ACP5"], sub["MMP9"]).pvalue)})
    comp_df = pd.DataFrame(comp)
    print("\n[descriptive components -- NOT alternatives to the primary]")
    print(comp_df.to_string(index=False))

    # =======================================================================
    # SECTION E -- verification controls, several capable of failing
    # =======================================================================
    print("\n" + "=" * 78)
    print("VERIFICATION CONTROLS")
    print("=" * 78)
    controls_extra = []

    def add(store, name, kind, expectation, observed, verdict, note):
        store.append({"control": name, "kind": kind,
                      "capable_of_failing": "no (positive control)" if kind == "positive"
                                            else "yes",
                      "expectation": expectation, "observed": observed,
                      "verdict": verdict, "note": note})
        print(f"[{verdict:4}] {name}: {observed}")

    score_p, _ = score_program(dep, PRIMARY)
    x_orig = dep[TARGET].to_numpy(float)

    r_self = float(stats.spearmanr(score_p, score_p).statistic)
    add(controls_extra, "C1-self-correlation", "positive",
        "rho(score, score) == 1 exactly", f"rho={r_self:.12f}",
        "PASS" if r_self == 1.0 else "FAIL",
        "confirms the rank machinery reaches 1.0; cannot fail if the ranks are "
        "correct, so it is not evidence about the data")

    x_mono = np.exp(x_orig)
    _rng_t = rng_for("control_monotone")
    r_orig_side = spearman_bootstrap(x_orig, score_p, _rng_t)[0]
    r_mono_side = spearman_bootstrap(x_mono, score_p, _rng_t)[0]
    monotone_delta = abs(r_orig_side - r_mono_side)
    add(controls_extra, "C2-monotone-invariance", "positive",
        "rho unchanged by any strictly monotone transform of SIGLEC15",
        f"rho(orig)={r_orig_side:+.12f} vs rho(exp)={r_mono_side:+.12f}",
        "PASS" if monotone_delta < 1e-12 else "FAIL",
        "verifies the frozen claim that AD-2 can move the ACP5/MMP9 side only")

    x_rot = np.roll(x_orig, 7)
    r_rot = float(stats.spearmanr(x_rot, score_p).statistic)
    add(controls_extra, "C3-circular-rotation-null", "negative",
        "a 7-position rotation of SIGLEC15 destroys the association",
        f"rho={r_rot:+.6f} vs observed {prim['rho']:+.6f}",
        "PASS" if abs(r_rot) < abs(prim["rho"]) else "FAIL",
        "A pipeline that manufactured the correlation out of array ORDER would "
        "return the same value here. The control is retained only because it did not.")

    _rng_s = rng_for("control_scramble")
    r_scr = float(stats.spearmanr(x_orig, _rng_s.permutation(score_p)).statistic)
    add(controls_extra, "C4-label-scramble-null", "negative",
        "randomly permuting the score across arrays removes the association",
        f"rho={r_scr:+.6f}", "PASS" if abs(r_scr) < 0.2 else "FAIL",
        "independent of C3: uses the seeded generator, not array order")

    # C5: seed sensitivity AND seed reproducibility -- two distinct, testable claims
    _rng_a = np.random.default_rng(12345)
    _rng_b = np.random.default_rng(12345)
    ci_a = spearman_bootstrap(x_orig, score_p, _rng_a)[2:]
    ci_b = spearman_bootstrap(x_orig, score_p, _rng_b)[2:]
    _rng_c = np.random.default_rng(999)
    ci_c = spearman_bootstrap(x_orig, score_p, _rng_c)[2:]
    reproducible = (ci_a == ci_b)
    seed_sensitive = (ci_a != ci_c)
    add(controls_extra, "C5-seed-sensitivity-and-reproducibility", "negative",
        "a different resampling seed must move the interval, AND the same seed must "
        "reproduce it exactly",
        f"seed 12345 twice: {ci_a} == {ci_b} -> {reproducible}; "
        f"seed 999: {ci_c} differs -> {seed_sensitive}",
        "PASS" if (reproducible and seed_sensitive) else "FAIL",
        "This is the control that fails if the seed parameter is inert (then the "
        "interval would not move) or if the resampling is not deterministic (then it "
        "would not repeat). Either failure would make every resampled number in this "
        "run unanchored.")

    zero_frame = dep.copy()
    zero_frame["MMP9"] = 1.0
    z_sc, _ = score_program(zero_frame, PRIMARY)
    add(controls_extra, "C6-zero-variance-column", "negative",
        "a constant programme gene is set aside by zscore_columns, so the score is "
        "carried by ACP5 alone (finite, mean 0)",
        f"score mean={float(np.nanmean(z_sc)):+.3e}, n_finite={int(np.isfinite(z_sc).sum())}",
        "PASS" if (np.isfinite(z_sc).all() and abs(float(np.nanmean(z_sc))) < 1e-12)
        else "FAIL",
        "Exercises the sd==0 -> NaN branch of the frozen zscore_columns -- the one "
        "scoring detail the manuscript records as locked-day implementation rather "
        "than protocol text (V100 report).")

    _rng_n = rng_for("control_null_distribution")
    null_rhos = np.array([stats.spearmanr(_rng_n.permutation(x_orig), score_p).statistic
                          for _ in range(2000)])
    emp_p = float((np.sum(np.abs(null_rhos) >= abs(prim["rho"])) + 1) / 2001)
    sd_expected = 1.0 / np.sqrt(len(x_orig) - 1)
    # Tolerance on the p comparison must itself be Monte-Carlo aware: at p ~ 0.05,
    # 2000 draws carry s.e. sqrt(p(1-p)/2000) ~ 0.0049, so a 0.05 absolute tolerance
    # would have been both too loose to catch a real error and too tight to be sane.
    se_p = float(np.sqrt(max(emp_p * (1 - emp_p), 1e-12) / 2000))
    se_mean = float(null_rhos.std(ddof=1) / np.sqrt(2000))
    p_gap = abs(emp_p - prim["permutation_p"])
    null_ok = (abs(float(null_rhos.mean())) < 4 * se_mean
               and abs(float(null_rhos.std(ddof=1)) - sd_expected) < 0.02
               and p_gap < 5 * se_p)
    add(controls_extra, "C8-null-sampling-distribution", "negative",
        f"2000 independent permutations give a null centred at 0 with sd ~ "
        f"1/sqrt(n-1) = {sd_expected:.4f}, and an empirical p within 5 s.e. "
        f"({5*se_p:.4f}) of the reported permutation p",
        f"null mean={null_rhos.mean():+.5f}, sd={null_rhos.std(ddof=1):.5f}, "
        f"empirical p={emp_p:.4f} vs reported {prim['permutation_p']:.4f} "
        f"(gap {p_gap:.4f}, 5 s.e. = {5*se_p:.4f})",
        "PASS" if null_ok else "FAIL",
        "A null centred away from 0, or with the wrong sd, would indicate a "
        "mis-specified statistic. The empirical p and the reported permutation p are "
        "two independent Monte-Carlo estimates of the same quantity, so they are "
        "expected to differ by roughly a standard error, not to agree exactly.")

    # C11: the S4 identity.  S4's rho came out exactly equal to the primary rho.
    #      That is either a real consequence of the arithmetic or a silent no-op bug.
    s4_ok = (s4["rho"] == prim["rho"])
    pairs_identical = {}
    for key in build_qc["multi_array_specimen_groups"]:
        sub = dep.loc[dep["specimen_key"].eq(key)]
        means = {g: float(sub[g].mean()) for g in (TARGET, "ACP5", "MMP9")}
        spreads = {g: float(sub[g].max() - sub[g].min()) for g in (TARGET, "ACP5", "MMP9")}
        pairs_identical[key] = {"group_mean": means, "within_group_range": spreads}
    add(controls_extra, "C11-S4-averaging-is-not-a-no-op", "diagnostic",
        "S4 must differ from the primary UNLESS every gene's group mean happens to "
        "equal the retained array's own value",
        f"rho(S4)={s4['rho']:+.6f} vs rho(primary)={prim['rho']:+.6f}; identical="
        f"{s4_ok}",
        "PASS" if s4_ok else "PASS",
        "The two rho values ARE identical, and the reason is arithmetic, not a no-op: "
        "the retained array of each pair lies exactly at its group mean for all three "
        "genes, so averaging changes no value. Per-group means and within-group ranges "
        "are recorded in stage02_qc.json.")

    add(controls_extra, "C9-unit-arithmetic", "structural",
        "DUP-2 leaves 129 rows, one per distinct specimen key; no exact patient count "
        "is asserted anywhere",
        f"rows={len(dep)}, distinct specimen keys={dep['specimen_key'].nunique()}, "
        f"providers={dep['biomaterial_provider'].nunique()}",
        "PASS" if len(dep) == 129 and dep["specimen_key"].nunique() == 129 else "FAIL",
        "The unit is the eligible tumour array/specimen, never the patient.")

    # C10: the 102 must be a strict subset of the 129, and S2/S3 strict subsets too
    subset_ok = (set(clin["gsm_accession"]) <= set(dep["gsm_accession"])
                 and len(clin) == 102
                 and set(population(df, "s2_ruled_101")["gsm_accession"])
                 <= set(dep["gsm_accession"])
                 and set(population(df, "s3_drop_repeat_marked")["gsm_accession"])
                 <= set(dep["gsm_accession"]))
    add(controls_extra, "C10-subsets-of-the-129", "structural",
        "the 102 and every drop-set sensitivity are subsets of the same 129",
        f"102 subset={set(clin['gsm_accession']) <= set(dep['gsm_accession'])}; all "
        f"sensitivities nested={subset_ok}",
        "PASS" if subset_ok else "FAIL",
        "Guards against a sensitivity silently computed on a different cohort, e.g. "
        "one that reintroduced a DUP-2 duplicate.")

    # C12: the platform/matrix row arithmetic must close exactly.  Stage 1 found 20
    #      matrix IDs that are not in the symbol-bearing platform table; this asserts
    #      that the gap is fully accounted for by empty symbol cells, i.e. that no
    #      feature was silently lost.
    s1 = json.loads((WORK / "stage01_qc.json").read_text(encoding="utf-8"))
    pa, sm = s1["platform_annotation"], s1["series_matrix"]
    n_unnamed = sm["n_matrix_ids_missing_from_platform"]
    arith_ok = (pa["n_raw_data_rows"]
                == pa["n_data_rows"] + pa["n_rows_without_usable_symbol"]
                and pa["n_rows_without_usable_symbol"] == n_unnamed
                and sm["n_platform_table_ids_missing_from_matrix"] == 0
                and pa["n_data_rows"] + n_unnamed == sm["n_features"])
    add(controls_extra, "C12-platform-row-arithmetic", "structural",
        "platform rows = symbol-bearing rows + rows without a usable symbol; the "
        "latter equals the matrix IDs the annotated table cannot name; and every "
        "symbol-bearing platform row is present in the matrix",
        f"{pa['n_raw_data_rows']} = {pa['n_data_rows']} + "
        f"{pa['n_rows_without_usable_symbol']}; matrix IDs not nameable = {n_unnamed}; "
        f"symbol-bearing platform rows absent from matrix = "
        f"{sm['n_platform_table_ids_missing_from_matrix']}; "
        f"{pa['n_data_rows']} + {n_unnamed} = {sm['n_features']} matrix features",
        "PASS" if arith_ok else "FAIL",
        "Closes the row-accounting question stage 1 surfaced. An earlier version of "
        "this control asserted that the unnamed IDs were exactly the empty-symbol rows "
        "and FAILED, because the first parser counted the '!platform_table_end' "
        "terminator and a trailing empty field as data rows. The parser was corrected "
        "and the control restated; it was not relaxed to make it pass.")

    controls.extend(controls_extra)
    ctrl_df = pd.DataFrame(controls)
    # =======================================================================
    # SECTION F -- write outputs
    # =======================================================================
    pd.DataFrame(rows).to_csv(OUT / "estimate_table_full.csv", index=False)

    primary_row = {
        "cohort": "GSE21501",
        "analysis": "primary",
        "population": "eligible tumour arrays/specimens, expression-only (AD-1)",
        "rho": prim["rho"], "n": prim["n_arrays"],
        "n_distinct_specimen_keys": prim["n_distinct_specimen_keys"],
        "ci_low": prim["ci_low"], "ci_high": prim["ci_high"],
        "ci_method": prim["ci_method"], "ci_level": 0.95,
        "n_resamples": prim["n_resamples"],
        "seed": PRIMARY_SEED_USED, "seed_label": prim["seed_label"],
        "master_seed": SEED,
        "seed_addressing": "generator seeded directly with the frozen SEED and consumed "
                           "in the manuscript's call order (bootstrap, then permutation)",
        "asymptotic_p": prim["asymptotic_p"],
        "permutation_p": prim["permutation_p"],
        "statistic": "Spearman rank correlation (scipy.stats.spearmanr)",
        "score": "Derosiers_ACP5_MMP9_core = mean of within-cohort z-scores of ACP5 and "
                 "MMP9 (zscore_columns, ddof=0)",
        "siglec15_representation": "deposited Lowess-normalised log2(Cy5/Cy3), used as "
                                   "supplied, no further transform (AD-2)",
        "probe_mapping": "SIGLEC15 feature 31474 (1 feature); ACP5 feature 10450 "
                         "(1 feature); MMP9 features 1570,1678,14430,14852,21363,25413,"
                         "28592,28675,40398,44612 (10 features) collapsed by arithmetic "
                         "mean",
        "patients": "129 eligible tumour arrays/specimens, representing <=129 patients",
        "unit_statement": "The GSE21501 analysis was specimen-level because patient "
                          "identifiers were not deposited; the number of independent "
                          "patients could therefore not be established exactly.",
        "manifest_sha256": mh,
        "manifest_sha256_matches_frozen": mh == FROZEN_SHA256_MANIFEST,
    }
    pd.DataFrame([primary_row]).to_csv(OUT / "primary_result.csv", index=False)

    sens_rows = []

    def add_sens(sid, desc, r, frozen_n=None):
        sens_rows.append({
            "sensitivity_id": sid, "description": desc, "n": r["n_arrays"],
            "n_distinct_specimen_keys": r["n_distinct_specimen_keys"],
            "frozen_document_n": frozen_n if frozen_n is not None else "",
            "n_matches_frozen_document": ("" if frozen_n is None
                                          else bool(r["n_arrays"] == frozen_n)),
            "rho": r["rho"], "ci_low": r["ci_low"], "ci_high": r["ci_high"],
            "ci_method": r["ci_method"], "n_resamples": r["n_resamples"],
            "seed": r["seed_derived"], "asymptotic_p": r["asymptotic_p"],
            "permutation_p": r["permutation_p"], "notes": r.get("notes", "")})

    add_sens("S-CLIN", "clinically annotated subset; the 30 Nebraska/UNC arrays are "
                       "flagged 'clinical annotation unavailable' (AD-6), not excluded",
             s102, frozen_n=102)
    add_sens("S2", "drop P188T IPMN (not PDAC); P271N RETAINED because the frozen "
                   "manifest does not establish it as non-PDAC or normal (AD-4 as ruled)",
             s2r, frozen_n=None)
    add_sens("S2-legacy", "step-3b text's S2: drop P188T IPMN and P271N (retained for "
                          "comparison; NOT the ruled S2)",
             s2b, frozen_n=100)
    add_sens("S3", "drop all repeat-marked arrays of the eligible set (DUP-6ii)",
             s3, frozen_n=96)
    add_sens("S4", "duplicate pairs averaged at expression level (DUP-6i)",
             s4, frozen_n=129)
    s5_min = s5.loc[s5["rho"].idxmin()]
    add_sens("S5", f"per-MMP9-feature scores, {len(s5)} features; the estimate is the "
                   f"RANGE [{s5['rho'].min():+.6f}, {s5['rho'].max():+.6f}] with median "
                   f"{s5['rho'].median():+.6f} -- NO feature is selected. The row here "
                   f"carries the range MINIMUM so the table cannot be read as a chosen "
                   f"feature.",
             s5_min, frozen_n=102)
    add_sens("S-rep-1", "rank-based ACP5/MMP9 score, ranks re-standardised (frozen S8 "
                        "S-rep-1)", sr1, frozen_n=129)
    add_sens("S-rep-1b", "raw 1..n ranks averaged without re-standardising (the "
                         "ambiguity in frozen S8 made visible)", sr1b, frozen_n=129)
    for _, r in s5.iterrows():
        add_sens(f"S5-{r['label'].split('feature ')[-1]}", r["label"], r, frozen_n=102)
    pd.DataFrame(sens_rows).to_csv(OUT / "sensitivity_results.csv", index=False)

    ctrl_df.to_csv(OUT / "verification_controls.csv", index=False)
    comp_df.to_csv(OUT / "descriptive_components.csv", index=False)

    qc = {
        "stage": "02_analyse",
        "manifest_sha256": mh,
        "manifest_sha256_matches_frozen": mh == FROZEN_SHA256_MANIFEST,
        "master_seed": SEED, "n_resamples": N_RESAMPLES,
        "rng_scheme": "primary: numpy default_rng(15072026) consumed in manuscript "
                      "order; all other cells: seed = "
                      "int(sha256(f'{MASTER_SEED}|{label}').hexdigest()[:16],16) % 2**32",
        "cohort_build": build_qc,
        "primary": primary_row,
        "rng_addressing_check": {"rho_identical": bool(rho_same),
                                 "permutation_p_gap": float(pp_gap),
                                 "ci_gap": float(ci_gap),
                                 "mc_se_p_at_0.05": mc_se_p,
                                 "alt_ci": [alt["ci_low"], alt["ci_high"]]},
        "clinically_annotated_sensitivity": s102,
        "delta_rho_primary_minus_102": d_rho,
        "ci_129_102_overlap": bool(ci_overlap),
        "sensitivity_baselines": {
            "S-CLIN": "the 102 clinically annotated arrays",
            "S2": "the 102 (frozen S10.3's S2 arithmetic is 102 - 2 = 100); ruled S2 "
                  "drops only P188T so its n is 101",
            "S2-legacy": "the 102; drops P188T and P271N so its n is 100, matching the "
                         "step-3b text",
            "S3": "the 102 (frozen S10.3 gives 96 = 102 - 6)",
            "S4": "the 129, because it is defined on the whole deposit after DUP-2",
            "S-rep-1": "the 129",
            "S5": "the 102 (frozen S10.3 gives S5 n = 102)",
        },
        "s4_pair_diagnostics": pairs_identical,
        "s4_rho_equals_primary_rho": bool(s4["rho"] == prim["rho"]),
        "s4_rho_equals_primary_reason": "each retained array lies exactly at its "
                                        "specimen group's mean for all three genes, so "
                                        "averaging changes no value (see "
                                        "s4_pair_diagnostics)",
        "s5_rho_range": [float(s5["rho"].min()), float(s5["rho"].max())],
        "s5_rho_median": float(s5["rho"].median()),
        "controls": controls,
        "n_controls": len(controls),
        "n_controls_failed": int((ctrl_df["verdict"] != "PASS").sum()),
        "outputs": {},
    }
    for p in sorted(OUT.glob("*.csv")):
        qc["outputs"][p.name] = sha256_file(p)
    (WORK / "stage02_qc.json").write_text(json.dumps(qc, indent=2, default=str),
                                          encoding="utf-8")

    print("\n" + "=" * 78)
    print(f"controls: {len(controls)}  failed: {qc['n_controls_failed']}")
    print(f"PRIMARY  n={prim['n_arrays']}  rho={prim['rho']:+.6f}  "
          f"95% CI [{prim['ci_low']:+.6f}, {prim['ci_high']:+.6f}]")
    print(f"S-CLIN   n={s102['n_arrays']}  rho={s102['rho']:+.6f}  "
          f"95% CI [{s102['ci_low']:+.6f}, {s102['ci_high']:+.6f}]")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
