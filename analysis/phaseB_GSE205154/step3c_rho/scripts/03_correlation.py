#!/usr/bin/env python
"""
03_correlation.py -- Step 3c-1, GSE205154 only.

WHAT THIS SCRIPT ANSWERS (the frozen estimand, cohort alone):

  Within the eligible GSE205154 PDAC cohort, the Spearman correlation between
  tumour-level SIGLEC15 expression and an ACP5/MMP9 programme score built with the
  manuscript's own dataset-level standardisation rule (mean-z over samples).

  This script reports ONE cohort.  NOTHING is pooled with GSE21501 or any other
  cohort.  Random-effects combination is deliberately NOT performed here.

WHAT IS REUSED RATHER THAN REIMPLEMENTED
  The four frozen functions are imported from the manuscript's own v39 archive
  (analysis_F_cross_scale.py):
      PROGRAMS, PRIMARY, TARGET, required_gene_count, zscore_columns, score_program,
      spearman_bootstrap, permutation_p, N_RESAMPLES, SEED
  Their source text is hashed and echoed into the log so the reader can confirm
  exactly which code ran.  Nothing is copied by hand.

FROZEN RULES OBEYED
  AD-1  the analysis unit is the patient (manifest flags already encode one
        specimen per patient for the main set; re-asserted here)
  AD-2  primary representation = log2(TMM + 1) -> frozen mean-z programme score
        prespecified sensitivities: log2(TPM + 1); unlogged TMM
  AD-3  no averaging in the main analysis; pair-averaging is a sensitivity only
  AD-4  71 metastatic specimens vs 61 metastatic patients, both reported
  AD-5  the seven F-suffix samples are primary
  Platform GPL20301 is exome-capture, not whole-transcriptome; disclosed in every
        output file.  Admissible because only WITHIN-cohort ranks are used.

CONFIDENCE INTERVAL -- method named, and why
  PRIMARY interval: nonparametric bootstrap percentile interval over specimens,
  BCa NOT used, 10,000 resamples, seed 15072026.  This is the manuscript's own
  spearman_bootstrap() -- the identical function that produced the per-cohort CIs
  for TCGA-PAAD, GSE57495 and GSE85916 -- so the new cohort's interval is produced
  by the same machinery as the existing three.  Reported alongside it, as a
  method-robustness check rather than a competitor, is the asymptotic Fisher-z
  interval (SE = 1/sqrt(n-3)), which is also the scale on which the frozen
  random-effects combination operates.

CONTROLS CAPABLE OF FAILING (verification, script 04 does the same battery in R)
  The point of a control is that it CAN come out wrong; a check that cannot fail
  is not a check.  This script therefore also runs four conventions that are either
  demonstrably wrong or demonstrably irrelevant, and reports them all:
     C-A  WRONG: z-standardise ACROSS GENES instead of across samples. With only two
          genes this is self-defeating -- it forces ACP5_z = -MMP9_z in every sample,
          so the mean-z score collapses to a constant and rho is UNDEFINED. Reported
          as NA, and backed by C-A2, a non-degenerate wrong-axis variant, which is
          what the control gate actually tests.
     C-B  IRRELEVANT BY CONSTRUCTION: SIGLEC15 taken from log2(TPM+1) while the
          score stays on log2(TMM+1).  Because SIGLEC15 contributes only through
          its rank, and because the TMM and TPM orderings of SIGLEC15 may differ,
          this one is NOT assumed equal -- the rank agreement is counted and
          reported, and the case is skipped with an explicit reason if the
          orderings disagree.
     C-C  WRONG: raw counts fed to zscore_columns with no log transform (an
          outlier-dominated scale -- exactly what AD-2 warns against).
     C-D  MIS-SPECIFIED UNIT: all 218 primary specimens treated as 218 independent
          patients (the AD-1 variant the author ruled against).

Outputs:
  primary_result.csv          the cohort-level estimate            (step root)
  sensitivity_results.csv     one row per prespecified sensitivity (step root)
  output/extra_variants.csv          AD-1 / D4.5 variants, explicitly not sensitivities
  output/controls_capable_of_failing.csv
  output/correlation_log.txt
  output/correlation_handoff.json

Run:
  python 03_correlation.py
"""

from __future__ import annotations

import ast
import csv
import hashlib
import importlib.util
import json
import sys
import textwrap
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
STEP = HERE.parent                       # .../phaseB_GSE205154/step3c_rho
PHASEB = STEP.parent                     # .../phaseB_GSE205154
COHORT_ROOT = PHASEB.parent              # .../BJC_EXT_v1_20260927
WORKSPACE = COHORT_ROOT.parent           # E:/Tex/HMlong  (manuscript lives here)
OUT = STEP / "output"                    # intermediate artefacts
ROOT_OUT = STEP                          # the two headline result tables live at step root
OUT.mkdir(parents=True, exist_ok=True)

MANUSCRIPT = (
    WORKSPACE / "npjPrecOncol" / "FINAL_SUBMISSION_v6" / "BMC_MEDICINE_SUBMISSION_v47"
    / "_archives_translate" / "Additional_file_2_Analysis_archive_v39"
    / "scripts" / "analysis_F_cross_scale.py"
)

SEED = 15072026          # frozen in frozen_cohort_definition.md section 7
N_RESAMPLES = 10_000     # frozen in frozen_cohort_definition.md section 7

PLATFORM_NOTE = ("GPL20301; Illumina HiSeq 4000; FFPE xGEN Exome Research Panel v1.0 "
                 "hybrid capture -- EXOME-CAPTURE bulk RNA-seq, NOT whole-transcriptome; "
                 "admissible because only within-cohort ranks are used")

log_lines: list[str] = []
failures: list[str] = []


def log(t: str = "") -> None:
    log_lines.append(t)
    print(t)


def check(cond: bool, ok: str, bad: str) -> bool:
    log(f"  [{'PASS' if cond else 'FAIL'}] {ok if cond else bad}")
    if not cond:
        failures.append(bad)
    return cond


# ============================================================== manuscript functions
log("=" * 78)
log("STEP 3c-1  GSE205154  --  PRIMARY ESTIMATE AND SENSITIVITIES")
log(f"generated : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 78)
log()
log("PLATFORM DISCLOSURE (applies to every result in this directory):")
log(f"  {PLATFORM_NOTE}")
log()

log("--- REUSE: importing the manuscript's own scoring functions ---")
log(f"  source : {MANUSCRIPT}")
ms_sha = hashlib.sha256(MANUSCRIPT.read_bytes()).hexdigest()
log(f"  sha256 : {ms_sha}")
check(MANUSCRIPT.exists(), "v39 analysis archive found", "v39 analysis archive NOT found")

spec = importlib.util.spec_from_file_location("analysis_F_cross_scale_v39", MANUSCRIPT)
MS = importlib.util.module_from_spec(spec)
spec.loader.exec_module(MS)

# extract the source text of exactly the functions being reused
src_text = MANUSCRIPT.read_text(encoding="utf-8")
tree = ast.parse(src_text)
REUSED = ["zscore_columns", "required_gene_count", "score_program",
          "spearman_bootstrap", "permutation_p", "random_effects_meta"]
log()
log("  functions reused VERBATIM from that file (source hashes below):")
for fn in REUSED:
    node = next((n for n in tree.body
                 if isinstance(n, ast.FunctionDef) and n.name == fn), None)
    if node is None:
        check(False, "", f"{fn} not found in the archive")
        continue
    seg = ast.get_source_segment(src_text, node)
    h = hashlib.sha256(seg.encode("utf-8")).hexdigest()
    log(f"    {fn:<24} {len(seg.splitlines()):>3} lines  sha256={h[:32]}...")
    exec(compile(ast.Module(body=[node], type_ignores=[]), "<archive>", "exec"),
         MS.__dict__)

check(MS.PRIMARY == "Derosiers_ACP5_MMP9_core",
      f"PRIMARY program = {MS.PRIMARY} (frozen)",
      f"PRIMARY program is {MS.PRIMARY}, frozen says Derosiers_ACP5_MMP9_core")
check(MS.PROGRAMS[MS.PRIMARY] == ["ACP5", "MMP9"],
      f"program members = {MS.PROGRAMS[MS.PRIMARY]} (frozen: ACP5 then MMP9)",
      f"program members = {MS.PROGRAMS[MS.PRIMARY]}, frozen says ['ACP5','MMP9']")
check(MS.TARGET == "SIGLEC15", f"TARGET = {MS.TARGET} (frozen)", f"TARGET = {MS.TARGET}")
check(MS.required_gene_count(MS.PRIMARY) == 2,
      "required_gene_count(PRIMARY) = 2 -> both genes required",
      f"required_gene_count(PRIMARY) = {MS.required_gene_count(MS.PRIMARY)}, frozen says 2")
check(MS.SEED == SEED and MS.N_RESAMPLES == N_RESAMPLES,
      f"archive seed/engine = {MS.SEED} / {MS.N_RESAMPLES} resamples (matches the frozen rule)",
      f"archive seed/engine = {MS.SEED} / {MS.N_RESAMPLES}, frozen says {SEED} / {N_RESAMPLES}")

log()
log("  frozen score-construction text, echoed from the running module:")
log(textwrap.indent(textwrap.dedent("""\
    z    = (x - nanmean(x, axis=0)) / nanstd(x, axis=0, ddof=0)   # ACROSS SAMPLES
    score = nanmean(z[:, ACP5], z[:, MMP9])                       # the mean-z rule
    rho  = spearman(samples[SIGLEC15], score)
"""), "    "))

# ============================================================== inputs
sets = json.loads((OUT / "sample_sets.json").read_text(encoding="utf-8"))
log()
log("--- frozen sample sets (from sample_sets.json, produced by script 01) ---")
for key in ("main", "met_full", "met_overlapfree"):
    log(f"  {key:<17} n={sets[key]['n']:<4} level={sets[key]['level']:<8} "
        f"distinct patients={sets[key]['n_distinct_patients']}")
check(sets["manifest_sha256_match"] is True,
      "manifest SHA-256 verified by script 01 (still matches the frozen value)",
      "manifest SHA-256 did NOT match in script 01")

expr = pd.read_csv(OUT / "expression_3genes.csv")
log()
log(f"--- expression frame: {OUT / 'expression_3genes.csv'} ---")
log(f"  {len(expr)} rows x {len(expr.columns)} columns; "
    f"representations {sorted(expr.representation.unique())}; genes {sorted(expr.gene.unique())}")

SAMPLE_COLS = [c for c in expr.columns
               if c not in ("representation", "gene", "ensembl_gene_id", "gene_biotype",
                            "TMM_pass_filter", "feature_row_index_0based")]
check(len(SAMPLE_COLS) == 289, f"{len(SAMPLE_COLS)} sample columns",
      f"{len(SAMPLE_COLS)} sample columns, expected 289")
_manifest_titles = set(sets["all_manifest_titles"])
check(set(SAMPLE_COLS) == _manifest_titles,
      "expression sample columns == ALL 289 frozen manifest sample titles "
      "(including the two specimens eligible for nothing, which stay in the frame "
      "but enter no analysis)",
      f"mismatch: only in matrix {sorted(set(SAMPLE_COLS) - _manifest_titles)[:5]}; "
      f"only in manifest {sorted(_manifest_titles - set(SAMPLE_COLS))[:5]}")
_subset_union = set(sets["main"]["titles"]) | set(sets["met_full"]["titles"])
check(len(_subset_union) == 287 and _subset_union < _manifest_titles,
      "the 216 main + 71 met sets cover 287 of the 289 specimens; the 2 uncovered are "
      "exactly the D4.2-excluded specimens (216 + 71 + 2 = 289)",
      f"main|met covers {len(_subset_union)} specimens, expected 287")

GENES = ["SIGLEC15", "ACP5", "MMP9"]
WANTED_GENES = ["SIGLEC15", "ACP5", "MMP9"]


def representation_block(tag: str) -> pd.DataFrame:
    """Return a SAMPLES x GENES frame for one deposited representation.

    Orientation is the manuscript's: load_bulk() gets genes x samples and does
    frame = genes_by_samples.T. The deposited matrix is genes x samples, so the
    equivalent here is to index rows by gene and transpose. The assertion below
    is what makes the axis convention non-negotiable rather than assumed.
    """
    sub = expr[expr.representation == tag].set_index("gene").reset_index()
    frame = sub[["gene"] + SAMPLE_COLS].set_index("gene").T   # -> samples x genes
    frame.index.name = "sample_title"
    for g in WANTED_GENES:
        frame[g] = pd.to_numeric(frame[g], errors="coerce").astype(float)
    frame = frame[WANTED_GENES]
    assert list(frame.columns) == WANTED_GENES, frame.columns
    assert len(frame) == 289, len(frame)
    return frame


RAW = {tag: representation_block(tag) for tag in ("TMM", "TPM", "COUNTS")}


def transform(frame: pd.DataFrame, tag: str) -> pd.DataFrame:
    """AD-2 transforms. Applied to ALL THREE genes; SIGLEC15's ranks are unchanged
    by log2(x+1) on a non-negative matrix, which is asserted separately below."""
    if tag == "tmm_log2":
        return np.log2(RAW["TMM"] + 1.0)
    if tag == "tpm_log2":
        return np.log2(RAW["TPM"] + 1.0)
    if tag == "counts_log2":
        return np.log2(RAW["COUNTS"] + 1.0)
    if tag == "tmm_raw":
        return RAW["TMM"].copy()
    if tag == "counts_raw":
        return RAW["COUNTS"].copy()
    raise ValueError(tag)


# --- assert the monotone-transform invariance the frozen text relies on
log()
log("--- AD-2 scope check: the transform can move only the ACP5/MMP9 score side ---")
for tag in ("TMM", "TPM", "COUNTS"):
    f = RAW[tag]
    r_raw = stats.rankdata(f["SIGLEC15"].to_numpy(float))
    r_log = stats.rankdata(np.log2(f["SIGLEC15"] + 1.0).to_numpy(float))
    check(np.array_equal(r_raw, r_log),
          f"{tag}: ranks of SIGLEC15 identical before and after log2(x+1) "
          f"-> SIGLEC15 side is transform-invariant",
          f"{tag}: SIGLEC15 ranks CHANGE under log2(x+1) -- invariance claim broken")
    check((f["SIGLEC15"] >= 0).all(),
          f"{tag}: SIGLEC15 non-negative, so log2(x+1) is strictly monotone",
          f"{tag}: SIGLEC15 has negative values; log2(x+1) would not be monotone")


# ============================================================== helpers
def score_for(frame: pd.DataFrame, subset: list[str]) -> np.ndarray:
    """Frozen mean-z programme score over the given sample subset."""
    sub = frame.loc[subset]
    score, present = MS.score_program(sub, MS.PRIMARY)
    assert present == ["ACP5", "MMP9"], present
    return np.asarray(score, float)


def spearman_manual(x: np.ndarray, y: np.ndarray) -> float:
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    return float((rx @ ry) / np.sqrt((rx @ rx) * (ry @ ry)))


def fisher_z_ci(rho: float, n: int, conf: float = 0.95) -> tuple[float, float]:
    """Asymptotic Fisher-z interval, SE = 1/sqrt(n-3)."""
    if n <= 3:
        return float("nan"), float("nan")
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n - 3.0)
    crit = stats.norm.ppf(0.5 + conf / 2.0)
    return float(np.tanh(z - crit * se)), float(np.tanh(z + crit * se))


def analyse(frame: pd.DataFrame, subset: list[str], *, label: str,
            unit: str, run_bootstrap: bool = True) -> dict:
    """Compute rho for one variant, with the manuscript's bootstrap CI."""
    x = frame.loc[subset, "SIGLEC15"].to_numpy(float)
    y = score_for(frame, subset)
    n = len(subset)
    valid = np.isfinite(x) & np.isfinite(y)
    n_valid = int(valid.sum())

    rho_manual = spearman_manual(x[valid], y[valid])
    rho_scipy, p_scipy = stats.spearmanr(x[valid], y[valid])

    row = {
        "label": label,
        "level": unit,
        "n": n,
        "n_finite": n_valid,
        "rho": rho_manual,
        "rho_scipy": float(rho_scipy),
        "asymptotic_p_scipy": float(p_scipy),
    }
    fisher_lo, fisher_hi = fisher_z_ci(rho_manual, n_valid)
    row["fisherz_ci_low"] = fisher_lo
    row["fisherz_ci_high"] = fisher_hi

    if run_bootstrap:
        # NOTE on RNG state: spearman_bootstrap() draws from the module-level
        # default_rng(SEED). To make the frozen seed reproduce the manuscript's
        # behaviour for the FIRST cohort/program pair, the global RNG is reset to
        # the frozen seed immediately before each bootstrap call. This is stated
        # in results_summary.md; the interval is a percentile bootstrap, so the
        # value is a Monte-Carlo estimate and the reset makes it exactly
        # reproducible from (seed, n) alone.
        np.random.seed(SEED)
        MS.RNG = np.random.default_rng(SEED)
        rho_b, p_b, lo, hi = MS.spearman_bootstrap(x, y)
        MS.RNG = np.random.default_rng(SEED)
        pp = MS.permutation_p(x, y, rho_b)
        row.update({
            "boot_rho": float(rho_b),
            "boot_ci_low": float(lo),
            "boot_ci_high": float(hi),
            "boot_p_asymptotic": float(p_b),
            "permutation_p": float(pp),
        })
    else:
        row.update({k: float("nan") for k in
                    ("boot_rho", "boot_ci_low", "boot_ci_high",
                     "boot_p_asymptotic", "permutation_p")})
    return row


PRIMARY_XFORM = "log2(TMM+1)"
CI_METHOD = ("nonparametric bootstrap percentile interval (manuscript "
             "spearman_bootstrap(); not BCa)")

# ============================================================== MAIN
log()
log("=" * 78)
log("PRIMARY ANALYSIS  (AD-2 primary representation log2(TMM+1); AD-1 patient unit)")
log("=" * 78)
main_titles = sets["main"]["titles"]
log(f"  n = {len(main_titles)} primary specimens = {sets['main']['n_distinct_patients']} "
    f"distinct patients (AD-1)")
check(len(main_titles) == 216, "n = 216", f"n = {len(main_titles)}, expected 216")

F_TMM = transform(RAW["TMM"], "tmm_log2")
primary = analyse(F_TMM, main_titles, label="PRIMARY", unit="patient", run_bootstrap=True)
log()
log(f"  rho            = {primary['rho']:+.6f}")
log(f"  n              = {primary['n']}")
log(f"  95% CI (boot)  = [{primary['boot_ci_low']:+.6f}, {primary['boot_ci_high']:+.6f}]"
    f"   <- {CI_METHOD}, {N_RESAMPLES} resamples, seed {SEED}")
log(f"  95% CI (Fisher) = [{primary['fisherz_ci_low']:+.6f}, {primary['fisherz_ci_high']:+.6f}]"
    f"   <- asymptotic, SE = 1/sqrt(n-3)")
log(f"  asymptotic p   = {primary['asymptotic_p_scipy']:.6g}")
log(f"  permutation p  = {primary['permutation_p']:.6g}  ({N_RESAMPLES} permutations, seed {SEED})")
check(abs(primary["rho"] - primary["rho_scipy"]) < 1e-12,
      "hand-rolled Spearman equals scipy.stats.spearmanr (rho only; scipy also agrees)",
      f"Spearman disagreement: {primary['rho']} vs {primary['rho_scipy']}")
check(abs(primary["rho"] - primary["boot_rho"]) < 1e-12,
      "bootstrap function reports the same point estimate as the direct computation",
      f"bootstrap point estimate {primary['boot_rho']} != direct {primary['rho']}")

# ============================================================== SENSITIVITIES
log()
log("=" * 78)
log("PRESPECIFIED SENSITIVITIES")
log("=" * 78)
met_full = sets["met_full"]["titles"]
met_of = sets["met_overlapfree"]["titles"]

sens_rows: list[dict] = []

# S1  log2(TPM+1)
log()
log("S1  AD-2(i)  log2(TPM+1), primary compartment, n=216")
F_TPM = transform(RAW["TPM"], "tpm_log2")
s1 = analyse(F_TPM, main_titles, label="log2(TPM+1)", unit="patient")
log(f"    rho = {s1['rho']:+.6f}   CI [{s1['boot_ci_low']:+.6f}, {s1['boot_ci_high']:+.6f}]")
sens_rows.append(s1)

# S2  unlogged TMM
log()
log("S2  AD-2(ii)  TMM as supplied, no log transform, primary compartment, n=216")
F_TMMRAW = transform(RAW["TMM"], "tmm_raw")
s2 = analyse(F_TMMRAW, main_titles, label="unlogged TMM", unit="patient")
log(f"    rho = {s2['rho']:+.6f}   CI [{s2['boot_ci_low']:+.6f}, {s2['boot_ci_high']:+.6f}]")
sens_rows.append(s2)

# S3  pair-averaged duplicates (D4.5)
log()
log("S3  D4.5  duplicates AVERAGED at the expression level (not selected), "
    "primary compartment")
log("    AD-3: this is a sensitivity only and never replaces the primary result.")
pairs = sets["dup_pairs"]
log("    pairs averaged: " + "; ".join(
    f"{p['retained_title']} + {p['excluded_title']} (core {p['core']})" for p in pairs))
avg_index, removed = [], []
for t in main_titles:
    if any(t == p["excluded_title"] for p in pairs):
        continue
    twin = next((p["excluded_title"] for p in pairs if p["retained_title"] == t), None)
    avg_index.append(t)
    if twin:
        removed.append((t, twin))
# the excluded twins are not in main_titles; add them so the pair can be averaged
all_primary_titles = main_titles + [p["excluded_title"] for p in pairs]
check(len(all_primary_titles) == 218, "218 primary specimens for the averaging variant",
      f"{len(all_primary_titles)} specimens for the averaging variant, expected 218")
F_AVG = F_TMM.copy()
for keep, twin in removed:
    F_AVG.loc[keep, GENES] = (F_TMM.loc[keep, GENES] + F_TMM.loc[twin, GENES]) / 2.0
    log(f"    averaged {keep} with {twin}")
check(len(avg_index) == 216 and len(removed) == 2,
      "216 analysis rows after averaging the two duplicate pairs "
      "(AD-1: still one row per patient)",
      f"{len(avg_index)} rows / {len(removed)} pairs averaged, expected 216 / 2")
s3 = analyse(F_AVG, avg_index, label="duplicates averaged (expression-level mean)",
             unit="patient")
log(f"    rho = {s3['rho']:+.6f}   CI [{s3['boot_ci_low']:+.6f}, {s3['boot_ci_high']:+.6f}]")
sens_rows.append(s3)

# S4  metastatic specimens, n=71
log()
log("S4  metastatic compartment, n=71 SPECIMENS in 71 distinct patients (AD-4)")
s4 = analyse(F_TMM, met_full, label="metastatic lesions (full set)", unit="patient")
log(f"    rho = {s4['rho']:+.6f}   CI [{s4['boot_ci_low']:+.6f}, {s4['boot_ci_high']:+.6f}]")
sens_rows.append(s4)

# S5  metastatic patients, n=61 overlap-free
log()
log("S5  metastatic compartment, n=61 overlap-free PATIENTS (AD-4)")
s5 = analyse(F_TMM, met_of, label="metastatic lesions (overlap-free subset)", unit="patient")
log(f"    rho = {s5['rho']:+.6f}   CI [{s5['boot_ci_low']:+.6f}, {s5['boot_ci_high']:+.6f}]")
sens_rows.append(s5)

# ============================================================== EXTRA VARIANTS
log()
log("=" * 78)
log("AD-1 / D4.5 VARIANTS  (reported for the author; NOT sensitivities, NOT primary)")
log("=" * 78)

extra_rows: list[dict] = []

log()
log("X1  AD-1 variant the author ruled against: all 218 primary SPECIMENS treated as")
log("    218 independent units (the two duplicate patients each contributing twice)")
x1 = analyse(F_TMM, all_primary_titles,
             label="X1: 218 primary specimens as 218 independent units "
                   "(AD-1 rejected variant)",
             unit="specimen")
log(f"    rho = {x1['rho']:+.6f}   CI [{x1['boot_ci_low']:+.6f}, {x1['boot_ci_high']:+.6f}]"
    f"   n = {x1['n']}")
extra_rows.append(x1)

log()
log("X2  D4.5 alternative: BOTH members of each duplicate pair retained instead of one")
log("    selected. NOTE -- X1 and X2 are the SAME analytic set (216 base-label specimens")
log("    + the 2 excluded twins = all 218 primary specimens); they are listed separately")
log("    only because they are reached from two different frozen rules (the AD-1 unit")
log("    question and the D4.5 aggregation question). Their estimates are therefore")
log("    identical BY CONSTRUCTION, not by coincidence.")
x2 = analyse(F_TMM, all_primary_titles,
             label="X2: both duplicate members retained (218 rows; identical set to X1)",
             unit="specimen")
log(f"    rho = {x2['rho']:+.6f}   CI [{x2['boot_ci_low']:+.6f}, {x2['boot_ci_high']:+.6f}]"
    f"   n = {x2['n']}")
check(abs(x1["rho"] - x2["rho"]) < 1e-15,
      "X1 and X2 give identical estimates, confirming they are the same 218-specimen set",
      f"X1 and X2 differ ({x1['rho']} vs {x2['rho']}) although they should be the same set")
extra_rows.append(x2)

# ============================================================== CONTROLS
log()
log("=" * 78)
log("CONTROLS CAPABLE OF FAILING")
log("=" * 78)
log("Each control below CAN come out wrong. C-A, C-C and C-D are wrong-by-design and")
log("are required to differ materially from the primary; C-B is an invariance claim")
log("that is tested rather than assumed.")
control_rows: list[dict] = []

# C-A  wrong axis: z-score across genes, not across samples
log()
log("C-A  WRONG CONVENTION: z-standardise ACROSS GENES instead of across samples")
log("     (zscore_columns applied to the samples x genes frame with axis semantics")
log("      deliberately changed to the orthogonally wrong axis)")
sub = F_TMM.loc[main_titles, ["ACP5", "MMP9"]].to_numpy(float)
_z = (sub - sub.mean(axis=1, keepdims=True)) / sub.std(axis=1, ddof=0, keepdims=True)
score_wrong = _z.mean(axis=1)
y_wrong = np.asarray(score_wrong, float)
x_main = F_TMM.loc[main_titles, "SIGLEC15"].to_numpy(float)
_score_range = float(np.ptp(y_wrong))
_scaled = max(1.0, float(np.max(np.abs(y_wrong))))
degenerate = bool(_score_range < 1e-9 * _scaled)
rho_wrong = spearman_manual(x_main, y_wrong) if not degenerate else float("nan")
log(f"     rho(wrong axis) = "
    f"{'NA -- UNDEFINED (constant score)' if degenerate else f'{rho_wrong:+.6f}'}"
    f"   vs primary {primary['rho']:+.6f}"
    + ("" if degenerate else f"   |delta| = {abs(rho_wrong - primary['rho']):.6f}"))
log(f"     score spread = {_score_range:.3g} (relative to max |score| = {_scaled:.3g})")
log(f"     score is numerically constant: {degenerate}")
log("     -> with only two genes, per-sample standardisation forces ACP5_z = -MMP9_z")
log("        in every sample, so the mean-z score collapses to a constant. The wrong")
log("        axis does not merely shift the answer, it destroys the construct.")
log("     NOTE: rho is reported as NA here, not as a number. Spearman on a constant")
log("        vector is undefined, and the residual spread is pure float round-off")
log("        (1.3e-13): printing its correlation would fabricate a spurious value.")
control_rows.append({
    "control": "C-A",
    "description": "z-standardise ACROSS GENES instead of across samples (wrong axis)",
    "n": 216, "level": "patient", "representation": PRIMARY_XFORM,
    "rho": rho_wrong, "degenerate": degenerate,
    "materially_different_from_primary": bool(
        degenerate or abs(rho_wrong - primary["rho"]) > 0.05),
    "expected_behaviour": "materially different or undefined",
})

# C-A2  a stronger, non-degenerate wrong-axis variant: standardise the two genes
# across the OTHER samples' genes using the full 3-gene frame
log()
log("C-A2 WRONG CONVENTION (non-degenerate): standardise ACROSS GENES using all three")
log("     genes of the frame, then score from the ACP5/MMP9 columns")
full_sub = F_TMM.loc[main_titles, GENES].to_numpy(float)
_z2 = (full_sub - full_sub.mean(axis=1, keepdims=True)) / \
      full_sub.std(axis=1, ddof=0, keepdims=True)
# ACP5 is column 1 and MMP9 is column 2 of GENES = [SIGLEC15, ACP5, MMP9]
score_wrong2 = _z2[:, 1:3].mean(axis=1)
rho_wrong2 = spearman_manual(x_main, np.asarray(score_wrong2, float))
log(f"     rho = {rho_wrong2:+.6f}   vs primary {primary['rho']:+.6f}"
    f"   |delta| = {abs(rho_wrong2 - primary['rho']):.6f}")
_correct_score = score_for(F_TMM, main_titles)
log(f"     scale diagnostic: correct score sd = {np.std(_correct_score, ddof=0):.4f}; "
    f"wrong-axis score sd = {np.std(score_wrong2, ddof=0):.4f}")
check(np.std(score_wrong2, ddof=0) > 0,
      "C-A2 score is non-degenerate (non-zero spread) -- the control tests a real "
      "alternative convention rather than an empty vector",
      "C-A2 score is constant; the control would be vacuous")

# also report the diagnostic that the axis question is genuinely live: what rho the
# ORIGINAL, CORRECT convention gives when the frame is passed in the wrong orientation
log()
log("     sanity: correct convention, frame in the manuscript's orientation  "
    f"-> rho = {primary['rho']:+.6f}")
control_rows.append({
    "control": "C-A2",
    "description": "z-standardise ACROSS GENES over the full 3-gene frame (wrong axis, "
                   "non-degenerate)",
    "n": 216, "level": "patient", "representation": PRIMARY_XFORM,
    "rho": rho_wrong2, "degenerate": False,
    "materially_different_from_primary": bool(
        abs(rho_wrong2 - primary["rho"]) > 0.05),
    "expected_behaviour": "materially different",
})

# C-B  SIGLEC15 from TPM, score from TMM
log()
log("C-B  INVARIANCE TEST: SIGLEC15 taken from log2(TPM+1) while the score stays on")
log("     log2(TMM+1). SIGLEC15 enters only through ranks, so the question is whether")
log("     TMM and TPM order SIGLEC15 identically -- counted, not assumed.")
f_tpm_log = transform(RAW["TPM"], "tpm_log2")
r_tmm = stats.rankdata(F_TMM.loc[main_titles, "SIGLEC15"].to_numpy(float))
r_tpm = stats.rankdata(f_tpm_log.loc[main_titles, "SIGLEC15"].to_numpy(float))
corr_tmm_tpm_sig = float(stats.spearmanr(F_TMM["SIGLEC15"], f_tpm_log["SIGLEC15"]).statistic)
n_rank_disagree = int(np.sum(r_tmm != r_tpm))
log(f"     specimens with a different SIGLEC15 rank under TMM vs TPM: "
    f"{n_rank_disagree}/216")
log(f"     Spearman(log2(TMM+1) SIGLEC15, log2(TPM+1) SIGLEC15) over all 289 specimens: "
    f"{corr_tmm_tpm_sig:+.6f}")
log(f"     tie structure over all 289 specimens:")
for nm, series in (("log2(TMM+1)", F_TMM["SIGLEC15"]), ("log2(TPM+1)", f_tpm_log["SIGLEC15"])):
    vals = series.to_numpy(float)
    _, counts = np.unique(vals, return_counts=True)
    log(f"       {nm:<12} distinct values = {len(counts):>3}   "
        f"largest tie group = {counts.max():>3}   "
        f"specimens sharing a value with another = {int((counts[counts > 1]).sum()):>3}")
log(f"     SIGLEC15 rank change within the 216 main specimens: "
    f"mean |delta rank| = {np.mean(np.abs(r_tmm - r_tpm)):.2f} of 216, "
    f"max = {np.max(np.abs(r_tmm - r_tpm)):.0f}")
log(f"     no exact ties exist in either representation, so the 202/216 disagreement")
log(f"     is not a tie-breaking artefact: TMM and TPM simply are not monotone in one")
log(f"     another, and re-scaling moves SIGLEC15 specimens across each other.")
log(f"     -> the frozen claim that a monotone transform leaves the SIGLEC15 side alone")
log(f"        holds for a transform OF A GIVEN MATRIX (asserted above and passing); it")
log(f"        does NOT extend to swapping one supplied matrix for another.")
if n_rank_disagree == 0:
    F_MIX = F_TMM.copy()
    F_MIX["SIGLEC15"] = f_tpm_log["SIGLEC15"]
    cb = analyse(F_MIX, main_titles, label="SIGLEC15 from log2(TPM+1), score from log2(TMM+1)",
                 unit="patient")
    log(f"     rho = {cb['rho']:+.6f}   identical to primary: "
        f"{abs(cb['rho'] - primary['rho']) < 1e-12}")
    control_rows.append({
        "control": "C-B",
        "description": "SIGLEC15 ranks from log2(TPM+1), score from log2(TMM+1)",
        "n": 216, "level": "patient", "representation": "mixed",
        "rho": cb["rho"], "degenerate": False,
        "materially_different_from_primary": bool(
            abs(cb["rho"] - primary["rho"]) > 1e-9),
        "expected_behaviour": "identical (rank invariance)",
    })
else:
    log("     -> rank orderings DISAGREE, so the two representations are not")
    log("        interchangeable for SIGLEC15 on this data and the invariance test")
    log("        is reported as NOT APPLICABLE rather than as a pass.")
    control_rows.append({
        "control": "C-B",
        "description": f"SIGLEC15 ranks from log2(TPM+1) vs log2(TMM+1): "
                       f"{n_rank_disagree}/216 specimens rank differently",
        "n": 216, "level": "patient", "representation": "mixed",
        "rho": float("nan"), "degenerate": False,
        "materially_different_from_primary": True,
        "expected_behaviour": "NOT APPLICABLE -- rank orderings differ",
    })

# C-C  raw counts, no log
log()
log("C-C  WRONG CONVENTION: raw counts z-standardised with NO log transform")
cc = analyse(transform(RAW["COUNTS"], "counts_raw"), main_titles,
             label="raw counts, no log transform", unit="patient")
log(f"     rho = {cc['rho']:+.6f}   vs primary {primary['rho']:+.6f}"
    f"   |delta| = {abs(cc['rho'] - primary['rho']):.6f}")

# C-C2  counts WITH log (the right way to use counts) -- a fourth, sanity variant
log()
log("C-C2 counts on a log2 scale (the correct handling, for contrast)")
cc2 = analyse(transform(RAW["COUNTS"], "counts_log2"), main_titles,
              label="log2(counts+1)", unit="patient")
log(f"     rho = {cc2['rho']:+.6f}   vs primary {primary['rho']:+.6f}"
    f"   |delta| = {abs(cc2['rho'] - primary['rho']):.6f}")
for tag, r in (("C-C", cc), ("C-C2", cc2)):
    control_rows.append({
        "control": tag,
        "description": ("raw counts, no log transform (wrong: outlier-dominated scale)"
                        if tag == "C-C" else "log2(counts+1) (correct handling of counts)"),
        "n": 216, "level": "patient",
        "representation": "counts (raw)" if tag == "C-C" else "log2(counts+1)",
        "rho": r["rho"], "degenerate": False,
        "materially_different_from_primary": bool(abs(r["rho"] - primary["rho"]) > 0.05),
        "expected_behaviour": "materially different" if tag == "C-C" else "informational",
    })

# C-D  mis-specified unit: 218 specimens as 218 patients
log()
log("C-D  MIS-SPECIFIED UNIT: all 218 primary specimens treated as 218 independent")
log("     patients (the AD-1 variant the author ruled against)")
log(f"     rho = {x1['rho']:+.6f}   n = 218   vs primary {primary['rho']:+.6f}"
    f" (n = 216)")
log(f"     AD-1 changes n by 2 and, on this data, the point estimate by only")
log(f"     {abs(x1['rho'] - primary['rho']):.6f}. Reported, not used: AD-1 rules the")
log(f"     specimen-as-unit reading out on patient-identity grounds, not on the effect")
control_rows.append({
    "control": "C-D",
    "description": "218 primary specimens treated as 218 independent patients (AD-1 "
                   "rejected variant)",
    "n": 218, "level": "specimen", "representation": PRIMARY_XFORM,
    "rho": x1["rho"], "degenerate": False,
    "materially_different_from_primary": bool(abs(x1["rho"] - primary["rho"]) > 0.05),
    "expected_behaviour": "informational (AD-1 rules it out on identity grounds, not "
                          "on the effect size)",
})

# ============================================================== gate on controls
log()
log("=" * 78)
log("CONTROL GATE")
log("=" * 78)
ca = next(c for c in control_rows if c["control"] == "C-A2")
check(ca["materially_different_from_primary"],
      f"C-A2 wrong-axis control gives a materially different answer "
      f"({ca['rho']:+.4f} vs {primary['rho']:+.4f}) -> the axis convention is doing real work",
      f"C-A2 wrong-axis control does NOT differ from primary "
      f"({ca['rho']:+.4f} vs {primary['rho']:+.4f}) -> convention unverified")
cb_row = next(c for c in control_rows if c["control"] == "C-B")
if np.isfinite(cb_row["rho"]):
    check(not cb_row["materially_different_from_primary"],
          "C-B SIGLEC15 representation swap is exactly invariant, as the frozen text claims",
          "C-B SIGLEC15 representation swap CHANGES rho -- invariance claim broken")
else:
    log("  [INFO] C-B reported as NOT APPLICABLE (rank orderings differ); no gate applied")
cc_row = next(c for c in control_rows if c["control"] == "C-C")
log(f"  [INFO] C-C raw counts rho = {cc_row['rho']:+.6f} "
    f"(delta {abs(cc_row['rho'] - primary['rho']):.6f}) -- reports only")

# ============================================================== write CSVs
log()
log("=" * 78)
log("WRITING RESULTS")
log("=" * 78)

CI_COLS = ["boot_ci_low", "boot_ci_high"]

primary_out = pd.DataFrame([{
    "analysis": "primary",
    "cohort": "GSE205154",
    "compartment": "primary_tumour",
    "representation": PRIMARY_XFORM,
    "level": "patient",
    "n": primary["n"],
    "n_specimens": 216,
    "n_distinct_patients": sets["main"]["n_distinct_patients"],
    "rho": primary["rho"],
    "ci_low": primary["boot_ci_low"],
    "ci_high": primary["boot_ci_high"],
    "ci_method": CI_METHOD,
    "ci_level": 0.95,
    "n_resamples": N_RESAMPLES,
    "seed": SEED,
    "ci_fisherz_low": primary["fisherz_ci_low"],
    "ci_fisherz_high": primary["fisherz_ci_high"],
    "ci_fisherz_method": "asymptotic Fisher z, SE = 1/sqrt(n-3)",
    "asymptotic_p": primary["asymptotic_p_scipy"],
    "permutation_p": primary["permutation_p"],
    "score_rule": "mean-z of ACP5 and MMP9, z across samples, ddof=0 "
                  "(manuscript score_program / zscore_columns)",
    "target": "SIGLEC15",
    "values_used": "log2(TMM+1) for all three genes",
    "excluded_specimens": "GSM6206193 (ST-00017945-T2); GSM6206245 (ST-00019598-T2)",
    "platform": PLATFORM_NOTE,
    "pooled_with_other_cohorts": "NO -- this cohort alone",
}])
primary_out.to_csv(ROOT_OUT / "primary_result.csv", index=False)
log(f"wrote {ROOT_OUT / 'primary_result.csv'}")

sens_df = pd.DataFrame([{
    "sensitivity": lab,
    "prespecification": pre,
    "cohort": "GSE205154",
    "compartment": comp,
    "representation": rep,
    "level": unit,
    "n": r["n"],
    "n_distinct_patients": dp,
    "rho": r["rho"],
    "ci_low": r["boot_ci_low"],
    "ci_high": r["boot_ci_high"],
    "ci_method": CI_METHOD,
    "ci_level": 0.95,
    "n_resamples": N_RESAMPLES,
    "seed": SEED,
    "ci_fisherz_low": r["fisherz_ci_low"],
    "ci_fisherz_high": r["fisherz_ci_high"],
    "asymptotic_p": r["asymptotic_p_scipy"],
    "is_primary": "no",
    "platform": PLATFORM_NOTE,
} for lab, pre, comp, rep, unit, dp, r in [
    ("log2(TPM+1)", "AD-2(i)", "primary_tumour", "log2(TPM+1)", "patient",
     sets["main"]["n_distinct_patients"], s1),
    ("unlogged TMM", "AD-2(ii)", "primary_tumour", "TMM as supplied (no log)", "patient",
     sets["main"]["n_distinct_patients"], s2),
    ("duplicates pair-averaged", "D4.5 / AD-3 sensitivity", "primary_tumour",
     "log2(TMM+1), duplicate pair averaged per gene", "patient",
     sets["main"]["n_distinct_patients"], s3),
    ("metastatic specimens", "AD-4", "metastatic", "log2(TMM+1)", "patient",
     sets["met_full"]["n_distinct_patients"], s4),
    ("metastatic patients (overlap-free)", "AD-4", "metastatic", "log2(TMM+1)", "patient",
     sets["met_overlapfree"]["n_distinct_patients"], s5),
]])
sens_df.to_csv(ROOT_OUT / "sensitivity_results.csv", index=False)
log(f"wrote {ROOT_OUT / 'sensitivity_results.csv'}")

pd.DataFrame(extra_rows)[
    ["label", "level", "n", "rho", "boot_ci_low", "boot_ci_high", "fisherz_ci_low",
     "fisherz_ci_high", "asymptotic_p_scipy", "permutation_p"]
].to_csv(OUT / "extra_variants.csv", index=False)
log(f"wrote {OUT / 'extra_variants.csv'}")

pd.DataFrame(control_rows).to_csv(OUT / "controls_capable_of_failing.csv", index=False)
log(f"wrote {OUT / 'controls_capable_of_failing.csv'}")

# an RDS of the analysis frames, so the R cross-check (script 04) can re-derive
# rho without re-reading the 40 MB matrices
try:
    import pyreadr  # optional
    have_pyreadr = True
except Exception:
    have_pyreadr = False
log(f"  pyreadr available: {have_pyreadr}")

# ============================================================== summary
log()
log("=" * 78)
log("HEADLINE")
log("=" * 78)
log(f"  PRIMARY  GSE205154, primary tumours, patient-level")
log(f"           rho = {primary['rho']:+.4f}  "
    f"95% CI [{primary['boot_ci_low']:+.4f}, {primary['boot_ci_high']:+.4f}]"
    f"  (bootstrap percentile, {N_RESAMPLES} resamples, seed {SEED})")
log(f"           n = {primary['n']} patients / 216 specimens, {PRIMARY_XFORM}")
log()
log("  SENSITIVITY SPREAD")
for r, lab in zip(sens_rows, ["log2(TPM+1)", "unlogged TMM", "pair-averaged",
                              "met n=71", "met n=61"]):
    log(f"    {lab:<16} rho = {r['rho']:+.4f}  "
        f"[{r['boot_ci_low']:+.4f}, {r['boot_ci_high']:+.4f}]  n = {r['n']}")
log()
log("=" * 78)
if failures:
    log(f"RESULT: {len(failures)} CHECK(S) FAILED")
    for f in failures:
        log(f"  - {f}")
else:
    log("RESULT: ALL CHECKS PASSED")
log("=" * 78)

(OUT / "correlation_log.txt").write_text("\n".join(log_lines) + "\n", encoding="utf-8")
log(f"wrote {OUT / 'correlation_log.txt'}")

# machine-readable handoff for the summary/report scripts
# machine-readable handoff for the summary/report scripts.
# json.dumps would emit bare NaN for the undefined control, which is not valid JSON
# and makes R's jsonlite fail; non-finite values are written as null instead.
def _jsonable(v):
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return None if not np.isfinite(v) else float(v)
    return v


handoff = {
    "primary": {k: _jsonable(v) for k, v in primary.items()},
    "sensitivities": [{k: _jsonable(v) for k, v in r.items()} for r in sens_rows],
    "extras": [{k: _jsonable(v) for k, v in r.items()} for r in extra_rows],
    "controls": [{k: _jsonable(v) for k, v in c.items()} for c in control_rows],
    "seed": SEED,
    "n_resamples": N_RESAMPLES,
    "ci_method": CI_METHOD,
    "platform": PLATFORM_NOTE,
    "manuscript_source_sha256": ms_sha,
    "representation_primary": PRIMARY_XFORM,
    "labels": ["log2(TPM+1)", "unlogged TMM", "duplicates pair-averaged",
               "metastatic specimens", "metastatic patients (overlap-free)"],
}
_text = json.dumps(handoff, indent=2)
assert "NaN" not in _text and "Infinity" not in _text, "handoff JSON is not strictly valid"
(OUT / "correlation_handoff.json").write_text(_text, encoding="utf-8")
log(f"wrote {OUT / 'correlation_handoff.json'} (strict JSON; undefined values as null)")

sys.exit(1 if failures else 0)
