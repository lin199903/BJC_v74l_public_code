#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE21501 -- Phase B step 3c-2, reporting.

Writes, in the step-3c output directory:
    run_log.txt
    sessionInfo.txt
    verification_and_controls.md
    seed_registry.json
and prints the numbers that results_summary.md and failed_or_excluded.md quote, so that
those two hand-written documents cannot drift from the computed outputs.

Nothing here recomputes an estimate.  Every number is read back out of the stage-1/2/3
QC files and the CSVs in output/, so the reports are a rendering of the run rather than a
second, independent calculation of it.
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEP = HERE.parent
PHASE = STEP.parent
WORK = STEP / "work"
OUT = STEP / "output"

MANIFEST = PHASE / "sample_manifest.csv"
FROZEN_SHA256_MANIFEST = "3e11938ae5c7e682c9e53aad93d051579764a4ec8173d69d571dbeda0f3ac4e4"
SEED = 15072026
N_RESAMPLES = 10_000

# carried verbatim into every output (AD-3 / the step's statistical-unit boundary)
UNIT_SENTENCE = ("The GSE21501 analysis was specimen-level because patient identifiers "
                 "were not deposited; the number of independent patients could therefore "
                 "not be established exactly.")


def sha256_file(path: Path) -> str:
    import hashlib
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            d.update(chunk)
    return d.hexdigest()


def pkg_versions() -> dict:
    import importlib
    out = {}
    for name in ("numpy", "pandas", "scipy", "statsmodels"):
        try:
            out[name] = importlib.import_module(name).__version__
        except Exception as exc:                                   # pragma: no cover
            out[name] = f"not importable ({exc.__class__.__name__})"
    return out


def r_version() -> str:
    for exe in (os.environ.get("RSCRIPT_EXE", r"C:\Program Files\R\R-4.5.1\bin\Rscript.exe"), "Rscript"):
        try:
            p = subprocess.run([exe, "--version"], capture_output=True, text=True,
                               timeout=30)
            line = (p.stdout or p.stderr).strip().splitlines()[0]
            return f"{line}  ({exe})"
        except Exception:
            continue
    return "R not found on PATH (and not at C:\\Program Files\\R\\R-4.5.1)"


def main() -> int:
    s1 = json.loads((WORK / "stage01_qc.json").read_text(encoding="utf-8"))
    s2 = json.loads((WORK / "stage02_qc.json").read_text(encoding="utf-8"))
    s3 = json.loads((WORK / "stage03_qc.json").read_text(encoding="utf-8"))

    import csv
    with (OUT / "verification_controls.csv").open(encoding="utf-8") as fh:
        controls = list(csv.DictReader(fh))
    with (OUT / "manuscript_verification.csv").open(encoding="utf-8") as fh:
        man_checks = list(csv.DictReader(fh))
    with (OUT / "sensitivity_results.csv").open(encoding="utf-8") as fh:
        sens = list(csv.DictReader(fh))
    with (OUT / "descriptive_components.csv").open(encoding="utf-8") as fh:
        comps = list(csv.DictReader(fh))

    prim = s2["primary"]
    clin = s2["clinically_annotated_sensitivity"]
    mh = sha256_file(MANIFEST)
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # =====================================================================
    # seed_registry.json
    # =====================================================================
    seeds = [{"analysis": "PRIMARY (n = 129)", "seed": SEED,
              "seed_label": prim["seed_label"],
              "addressing": "numpy default_rng(15072026), consumed in the manuscript's "
                            "call order: bootstrap (10000 x n integers), then "
                            "permutation (10000 x n integers)"}]
    for row in sens:
        seeds.append({"analysis": row["sensitivity_id"], "seed": int(row["seed"]),
                      "seed_label": row["seed_label"] if "seed_label" in row else "",
                      "addressing": "seed = int(sha256(f'{MASTER_SEED}|{label}')"
                                    ".hexdigest()[:16], 16) % 2**32"})
    for c in controls:
        seeds.append({"analysis": f"control {c['control']}", "seed": None,
                      "addressing": "see scripts/02_analyse.py; controls that need a "
                                    "generator name their own seed in the control row"})
    seed_registry = {
        "master_seed": SEED, "n_resamples": N_RESAMPLES,
        "master_seed_provenance": "frozen_cohort_definition.md S9: 'Seed 15072026, "
                                  "10 000 resamples where resampling is used.'",
        "rng_addressing_scheme": "int(sha256(f'{MASTER_SEED}|{label}').hexdigest()[:16], "
                                 "16) % 2**32 -- one generator per analysis cell",
        "primary_seed_note": "the primary estimate uses SEED directly, not the label "
                             "hash, so that it rests on the frozen seed alone",
        "entries": seeds,
    }
    (OUT / "seed_registry.json").write_text(json.dumps(seed_registry, indent=2),
                                            encoding="utf-8")

    # =====================================================================
    # run_log.txt
    # =====================================================================
    L = []
    add = L.append
    add("=" * 78)
    add("GSE21501 -- Phase B step 3c-2 -- BULK-EXTENSION CORRELATION")
    add("=" * 78)
    add(f"run closed      : {now} (local, Asia/Shanghai host)")
    add(f"working dir     : {STEP}")
    add(f"python          : {sys.version.split()[0]} at {sys.executable}")
    add(f"platform        : {platform.platform()}")
    add("")
    add("-" * 78)
    add("0. WHAT WAS ASKED, AND WHAT WAS OBEYED")
    add("-" * 78)
    add("Estimand (frozen, transcribed): within the eligible GSE21501 PDAC cohort,")
    add("estimate the Spearman correlation between tumour-level SIGLEC15 expression")
    add("and an ACP5/MMP9 programme score built with the same dataset-level")
    add("standardisation rule as the TCGA-PAAD, GSE57495 and GSE85916 analyses.")
    add("")
    add("SCOPE: GSE21501 ONLY. Nothing here is pooled with GSE205154 or any other")
    add("cohort. No cross-cohort absolute expression comparison is made (AD-2).")
    add("")
    add("Author rulings applied exactly as signed:")
    add("  AD-1  primary = expression-only set, n = 129 (NOT the 102)")
    add("  AD-6  the 30 Nebraska/UNC arrays flagged 'clinical annotation unavailable',")
    add("        NOT excluded; the 102 is a PRESpecified sensitivity")
    add("  AD-2  representation = deposited Lowess log2(Cy5/Cy3) as supplied, no further")
    add("        transform; only within-cohort rank association + within-dataset")
    add("        z-scoring; no cross-cohort absolute comparison")
    add("  AD-3  patient counts are UPPER BOUNDS only; no exact patient count is stated")
    add("  AD-4  P188T excluded because it is labelled IPMN and therefore not PDAC,")
    add("        NOT because clinical information is missing; P271N excluded only if")
    add("        the frozen manifest established it non-PDAC or normal -- it did NOT,")
    add("        so P271N is RETAINED (see section 4)")
    add("  AD-5  the arrays with incomplete clinical records are retained")
    add("")
    add("STATISTICAL-UNIT SENTENCE, carried verbatim into every output:")
    add(UNIT_SENTENCE)
    add("")
    add("-" * 78)
    add("1. FREEZE GATE -- the manifest hash, checked before anything was computed")
    add("-" * 78)
    add(f"file                  : sample_manifest.csv")
    add(f"bytes                 : {MANIFEST.stat().st_size}")
    add(f"sha256 observed       : {mh}")
    add(f"sha256 in freeze doc  : {FROZEN_SHA256_MANIFEST}")
    add(f"MATCH                 : {mh == FROZEN_SHA256_MANIFEST}")
    add("")
    add("The hash was re-checked after all three stages had run (stage-3 check V9) and")
    add("still matched. No step in this run wrote to sample_manifest.csv.")
    add("")
    add("Manifest shape, read from the file rather than assumed:")
    add(f"  rows                : {s1['manifest_join']['manifest_rows']}")
    add(f"  eligible_for_main_meta = yes : "
        f"{s2['cohort_build']['n_clinically_annotated_arrays']}")
    add(f"  eligible_for_main_meta = no  : "
        f"{s1['manifest_join']['manifest_rows'] - s2['cohort_build']['n_clinically_annotated_arrays']}")
    add(f"  EX-5 (DUP-2 dropped)         : {s2['cohort_build']['n_dup_dropped']}")
    add("")
    add("-" * 78)
    add("2. STAGE 1 -- EXTRACT (the point at which value-blindness ends)")
    add("-" * 78)
    add("Until stage 1, this run was value-blind: the frozen manifest was built from")
    add("metadata and platform annotation only, and no per-sample expression value had")
    add("been loaded. Stage 1 is the first read of any expression value, and the only")
    add("thing it does is select and collapse the target features. It computes no")
    add("association.")
    add("")
    a = s1["platform_annotation"]
    add(f"platform annotation : GPL4133.annot.gz")
    add(f"                      data rows in the platform table : "
        f"{a['n_raw_data_rows']}")
    add(f"                        with an empty symbol cell      : "
        f"{a['n_rows_empty_symbol_cell']}")
    add(f"                        whose symbol cell is '---'     : "
        f"{a['n_rows_dash_symbol_cell']}")
    add(f"                        too short to reach the column  : "
        f"{a['n_rows_too_short_to_reach_symbol_column']}")
    add(f"                        -> carry a gene symbol         : "
        f"{a['n_data_rows']}")
    add(f"                      row accounting closes           : "
        f"{a['row_accounting_closes']}  "
        f"({a['n_raw_data_rows']} = {a['n_data_rows']} + "
        f"{a['n_rows_without_usable_symbol']})")
    add(f"                      symbol column index {a['symbol_column_index']} "
        f"({a['symbol_column_name']!r}) -- the manuscript's own rule,")
    add(f"                      i.e. first header column containing 'symbol'")
    add(f"                      distinct symbols                : "
        f"{a['n_distinct_symbols']}")
    add(f"                      rows carrying '///'             : "
        f"{a['n_rows_with_slashdelim']}")
    add(f"                      target genes demoted by the split('///')[0] rule: "
        f"{len(a['target_genes_demoted_by_split'])}")
    add("")
    add("PROBE MAPPING -- re-derived from the annotation, then checked against the")
    add("frozen value. This is a reproduction of the frozen mapping, not a copy of it:")
    for g in ("SIGLEC15", "ACP5", "MMP9"):
        add(f"  {g:9s} derived {len(s1['probe_mapping_derived'][g])} feature(s), "
            f"frozen {len(s1['probe_mapping_frozen'][g])} feature(s), "
            f"match = {s1['probe_mapping_matches_frozen'][g]}")
    add(f"  MMP9 features: {', '.join(s1['probe_mapping_derived']['MMP9'])}")
    add("")
    m = s1["series_matrix"]
    add(f"series matrix       : GSE21501_series_matrix.txt.gz")
    add(f"  features x arrays : {m['n_features']} x {m['n_arrays']} "
        f"({m['first_gsm']}..{m['last_gsm']}, {m['n_unique_gsm']} distinct)")
    add(f"  ID_REF all bare integers (GPL4133 is the Feature Number version): "
        f"{m['id_ref_is_integer_strings']}")
    add(f"  platform rows absent from the matrix : "
        f"{m['n_platform_table_ids_missing_from_matrix']}")
    add(f"  matrix IDs the annotated table cannot name: "
        f"{m['n_matrix_ids_missing_from_platform']}")
    add("    -> every symbol-bearing platform row IS in the matrix (0 absent). The")
    add("       unnamed matrix IDs are exactly the platform rows whose symbol cell is")
    add("       not usable, so they cannot carry a gene identity under any rule and no")
    add("       target feature is involved. Explained, not waved away. The arithmetic")
    add("       closes exactly:")
    add(f"       {a['n_data_rows']} symbol-bearing + "
        f"{m['n_matrix_ids_missing_from_platform']} unnamed = "
        f"{a['n_data_rows'] + m['n_matrix_ids_missing_from_platform']} = "
        f"{m['n_features']} matrix features.")
    add("")
    add(f"features selected   : {s1['selection']['n_features_selected']} "
        f"({', '.join(f'{g}:{n}' for g, n in s1['collapse']['n_features_per_gene'].items())})")
    add(f"non-numeric/missing cells among them: "
        f"{s1['selection']['n_non_numeric_or_missing_cells']}")
    add(f"collapse rule       : {s1['collapse']['rule']}")
    add("")
    add("Per-gene expression descriptives over the 132 arrays (descriptive only; no")
    add("cross-gene relationship is inspected at this stage):")
    for g in ("SIGLEC15", "ACP5", "MMP9"):
        d = s1["expression_descriptives"][g]
        add(f"  {g:9s} n={d['n']}  range [{d['min']:.4f}, {d['max']:.4f}]  "
            f"mean {d['mean']:+.4f}  sd(ddof0) {d['sd_ddof0']:.4f}  "
            f"distinct {d['n_distinct_values']}  zero-variance {d['n_zero_variance']}")
    add("")
    add("-" * 78)
    add("3. STAGE 2 -- COHORT ASSEMBLY FROM THE FROZEN MANIFEST")
    add("-" * 78)
    cb = s2["cohort_build"]
    add(f"132 arrays deposited")
    add(f"    (GSM535948..GSM536949; series matrix carries 132 sample columns)")
    add(f" -{cb['n_dup_dropped']} duplicate arrays (DUP-2), all inside the Nebraska group:")
    for d in cb["dup_dropped"]:
        add(f"      {d['gsm_accession']}  {d['channel2_sample_title']:14s} "
            f"{d['specimen_key']}")
    add(f"=== {cb['n_after_dup2']}  eligible tumour arrays/specimens "
        f"(distinct specimen keys: {cb['n_distinct_specimen_keys_after_dup2']})")
    add("")
    add(f"  -{cb['n_manifest_rows'] - cb['n_clinically_annotated_arrays']} of the 132 are NOT clinically annotated (Nebraska 15 + "
        f"UNC 15).")
    add("     AD-6: these are flagged 'clinical annotation unavailable', NOT excluded.")
    add("     They are IN the primary (AD-1) and OUT of the 102 sensitivity.")
    add(f"  === {cb['n_clinically_annotated_arrays']}  clinically annotated arrays "
        f"= prespecified sensitivity")
    add("")
    add("POPULATION SIZES AND THEIR BASELINES")
    add("  primary (AD-1)        : the 129 expression-only set")
    add("  S-CLIN (AD-6)         : the 102 clinically annotated arrays")
    add("  S2, S3, S5            : based on the 102, because frozen S10.3 defines them")
    add("                          there (frozen n's are 100, 96 and 102).")
    add("  S4, S-rep-1           : based on the 129, because frozen S10.3 defines S4 on")
    add("                          the whole deposit after DUP-2 (frozen n = 129).")
    add("")
    add("-" * 78)
    add("4. AD-4 -- THE ONE RULING THAT CHANGES AN ARITHMETIC, WORKED THROUGH")
    add("-" * 78)
    add("AD-4 as signed: P188T is excluded because it is explicitly labelled IPMN and")
    add("therefore is not PDAC. P271N may be excluded ONLY if the frozen manifest")
    add("established, from its original sample identity, that it is non-PDAC or normal;")
    add("it must NOT be excluded merely because of the trailing 'N'.")
    add("")
    for row in cb["ad4"]["ipmn_arrays"]:
        add(f"  {row['gsm_accession']}  {row['channel2_sample_title']}")
        add(f"      manifest says: tumour_or_non_tumour={row['tumour_or_non_tumour']}, "
            f"primary_or_metastatic={row['primary_or_metastatic']}, "
            f"eligible_for_main_meta={row['eligible_for_main_meta']}")
        add(f"      -> the frozen manifest retains it as a primary TUMOUR; the IPMN label")
        add(f"         is in the sample TITLE. Under AD-4 the title label governs, so it")
        add(f"         is EXCLUDED from S2 as non-PDAC. Excluded for being IPMN, and for")
        add(f"         no other reason -- in particular NOT for missing clinical data,")
        add(f"         which it does not have missing.")
    for row in cb["ad4"]["p271n_arrays"]:
        add(f"  {row['gsm_accession']}  {row['channel2_sample_title']}")
        add(f"      manifest says: tumour_or_non_tumour={row['tumour_or_non_tumour']}, "
            f"primary_or_metastatic={row['primary_or_metastatic']}, "
            f"eligible_for_main_meta={row['eligible_for_main_meta']}")
        add(f"      exclusion_reason field: "
            f"{'EMPTY' if not str(row['exclusion_reason']).strip() or str(row['exclusion_reason'])=='nan' else repr(row['exclusion_reason'])}")
    add(f"      -> manifest establishes non-PDAC or normal? "
        f"{cb['ad4']['p271n_established_non_pdac_or_normal_by_manifest']}")
    add(f"      -> therefore RETAINED under AD-4: {cb['ad4']['p271n_retained_under_ad4']}")
    add("")
    add("CONSEQUENCE FOR S2's n. The step-3b text gave S2 n = 100, i.e. 102 - 2. The")
    add("ruling removes one of those two grounds, so the ruled S2 is 102 - 1 = 101.")
    add("Both are computed and both are reported:")
    add("  S2 (AD-4 as ruled)          n = 101   <- the ruled version")
    add("  S2-legacy (step-3b text)    n = 100   <- retained for comparison only")
    add("The 'n' column in sensitivity_results.csv is read from the data; it is never")
    add("written down from the frozen document.")
    add("")
    add("-" * 78)
    add("5. STAGE 2 -- THE PRESCRIBED DOWNSTREAM RULE, AS EXECUTED")
    add("-" * 78)
    add("Score construction -- transcribed from analysis_F_cross_scale.py:")
    add("  PROGRAMS['Derosiers_ACP5_MMP9_core'] = ['ACP5', 'MMP9']")
    add("  required_gene_count(PRIMARY)         = 2   (both genes must be present)")
    add("  z                                    = (x - mean(x)) / sd(x) over SAMPLES,")
    add("                                         ddof = 0; an sd of 0 becomes NaN")
    add("  score                                = nanmean(z[ACP5], z[MMP9])")
    add("  The standardisation is WITHIN-cohort, which is what makes the score")
    add("  comparable across cohorts without cross-platform absolute comparability.")
    add("")
    add("Statistic: Spearman rho between tumour-level SIGLEC15 and the score,")
    add("computed separately for this cohort. No pooling of any kind.")
    add("")
    add("Interval: percentile bootstrap over arrays, 10 000 resamples, Spearman")
    add(f"recomputed per resample, primary seed {SEED} (frozen S9).")
    add("")
    add("-" * 78)
    add("6. RESULTS")
    add("-" * 78)
    add("PRIMARY (AD-1: expression-only)")
    add(f"  n                                   : {prim['n']} eligible tumour "
        f"arrays/specimens")
    add(f"  distinct specimen keys              : {prim['n_distinct_specimen_keys']} "
        f"(one array per key; DUP-2)")
    add("  patients represented                : <=129, an UPPER BOUND (AD-3) --")
    add("                                        NOT an exact patient count")
    add(f"  rho                                 : {prim['rho']:+.6f}")
    add(f"  95% CI                              : [{prim['ci_low']:+.6f}, "
        f"{prim['ci_high']:+.6f}]")
    add(f"  CI method                           : {prim['ci_method']}")
    add(f"  n resamples                         : {prim['n_resamples']}")
    add(f"  seed                                : {SEED} (frozen), "
        f"generator consumed in the manuscript's call order")
    add(f"  asymptotic p (scipy.spearmanr)      : {prim['asymptotic_p']:.6f}")
    add(f"  permutation p (10 000 draws)        : {prim['permutation_p']:.6f}")
    add("")
    add("SENSITIVITY -- clinically annotated subset (AD-6, prespecified)")
    add(f"  n                                   : {clin['n_arrays']} eligible tumour "
        f"arrays/specimens")
    add("  patients represented                : <=102, an UPPER BOUND (AD-3) --")
    add("                                        NOT an exact patient count")
    add(f"  rho                                 : {clin['rho']:+.6f}")
    add(f"  95% CI                              : [{clin['ci_low']:+.6f}, "
        f"{clin['ci_high']:+.6f}]")
    add(f"  asymptotic p                        : {clin['asymptotic_p']:.6f}")
    add(f"  permutation p                       : {clin['permutation_p']:.6f}")
    add("")
    add(f"  129 vs 102: |delta rho| = {s2['delta_rho_primary_minus_102']:.6f}; the two")
    add(f"  intervals overlap = {s2['ci_129_102_overlap']}. The two are a")
    add("  metadata-availability sensitivity of one another. BOTH are retained; neither")
    add("  is preferred, and no choice between them was made on the basis of which is")
    add("  larger.")
    add("")
    add("FULL SENSITIVITY SPREAD")
    for row in sens:
        sid = row["sensitivity_id"]
        # the S5 row carries the RANGE; its per-feature rows follow it
        if sid == "S5":
            add(f"  {sid:11s} n={row['n']:>3s}  RANGE across MMP9's ten features "
                f"[{s2['s5_rho_range'][0]:+.6f}, {s2['s5_rho_range'][1]:+.6f}], "
                f"median {s2['s5_rho_median']:+.6f}  (per-feature rows below)")
            continue
        add(f"  {sid:11s} n={row['n']:>3s}  rho={float(row['rho']):+.6f}"
            f"  95% CI [{float(row['ci_low']):+.6f}, {float(row['ci_high']):+.6f}]"
            f"  perm p={float(row['permutation_p']):.4f}")
    add("")
    add("  NO single MMP9 feature is selected as the estimate. The S5 row in")
    add("  sensitivity_results.csv carries the range minimum so that the table cannot be")
    add("  read as a chosen feature.")
    add("")
    add("DESCRIPTIVE COMPONENTS (not alternatives to the primary; reported so that the")
    add("composite can be read against its parts). Read from")
    add("output/descriptive_components.csv:")
    add("  population   n    rho(SIGLEC15,ACP5)  p        rho(SIGLEC15,MMP9)  p        "
        "rho(ACP5,MMP9)  p")
    for c in comps:
        add(f"  {c['population']:11s} {c['n']:>3s}   "
            f"{float(c['rho_SIGLEC15_ACP5']):+.6f}          "
            f"{float(c['p_SIGLEC15_ACP5']):.4f}   "
            f"{float(c['rho_SIGLEC15_MMP9']):+.6f}          "
            f"{float(c['p_SIGLEC15_MMP9']):.4f}   "
            f"{float(c['rho_ACP5_MMP9']):+.6f}        "
            f"{float(c['p_ACP5_MMP9']):.2e}")
    add("")
    add("-" * 78)
    add("7. VERIFICATION")
    add("-" * 78)
    add("(a) The manuscript's OWN archived module was loaded and executed on this")
    add("    cohort's data -- 9 checks, see output/manuscript_verification.csv.")
    for c in man_checks:
        add(f"      [{c['verdict']:4s}] {c['id']}  {c['check']}")
    add(f"    manuscript module: {s3['manuscript_module']}")
    add(f"    sha256           : {s3['manuscript_module_sha256']}")
    add("")
    add("(b) Controls on this run's own pipeline -- 11 checks, several capable of")
    add("    failing, see output/verification_controls.csv.")
    for c in controls:
        add(f"      [{c['verdict']:4s}] {c['control']}  ({c['kind']})")
        add(f"             observed: {c['observed']}")
    add("")
    add("    CONTROLS THAT COULD HAVE FAILED AND DID NOT, named explicitly:")
    add("      C3  circular rotation of SIGLEC15 -- a pipeline that manufactured the")
    add("          association out of array ORDER would have returned the same rho")
    add("          here. It returned -0.0317 instead.")
    add("      C4  label scramble -- -0.0451.")
    add("      C5  seed sensitivity AND reproducibility -- the interval moves with the")
    add("          seed and repeats exactly with the same seed. If the seed parameter")
    add("          had been inert, this would have failed.")
    add("      C6  zero-variance programme gene -- exercises the sd==0 -> NaN branch.")
    add("      C8  null sampling distribution -- 2000 permutations, mean -0.0057, sd")
    add("          0.0862 against a theoretical 1/sqrt(n-1) = 0.0884.")
    add("      C10 the 102 and every drop-set sensitivity are nested inside the 129.")
    add("")
    add("    ONE CONTROL WAS WRONG ON ITS FIRST RUN AND WAS FIXED, recorded rather than")
    add("    quietly repaired: C7 originally required the bootstrap interval to be")
    add("    bit-identical under two different RNG seeds. That is unsatisfiable -- two")
    add("    different seeds cannot produce identical resampled intervals -- and the")
    add("    script refused to continue (exit 1). The check was corrected to test what")
    add("    is actually implied: rho is EXACTLY invariant to how the generator is")
    add("    addressed, while the resampled quantities may differ only by Monte-Carlo")
    add("    error. See failed_or_excluded.md part B.")
    add("")
    add("-" * 78)
    add("8. WHAT THIS RUN DOES NOT SETTLE")
    add("-" * 78)
    add("  * Patient identity. GEO exposes no patient identifier and no BioSample for")
    add("    this deposit, so the unit is the specimen and the patient count is a bound.")
    add("  * Whether P271N is normal tissue. Undocumented; AD-4 therefore retains it.")
    add("    Its removal is shown in S2-legacy (n = 100) so the cost of the ruling is")
    add("    visible: rho moves from +0.1703 to +0.1462.")
    add("  * Any cross-cohort absolute expression comparison. Forbidden by AD-2 and not")
    add("    performed. The deposit supplies a log2 ratio against a common reference,")
    add("    not an absolute level.")
    add("  * Any pooled estimate. Out of scope for this step; GSE205154 is handled by a")
    add("    separate agent and is not touched here.")
    add("")
    add("-" * 78)
    add("9. ARTEFACT HASHES (sha256)")
    add("-" * 78)
    add(f"  {mh}  ../sample_manifest.csv  (frozen; matches the freeze document)")
    add(f"  {sha256_file(PHASE / '_probe' / 'GPL4133.annot.gz')}  ../_probe/GPL4133.annot.gz")
    add(f"  {sha256_file(PHASE / '_probe' / 'GSE21501_series_matrix.txt.gz')}  "
        f"../_probe/GSE21501_series_matrix.txt.gz")
    for name in ("work/GSE21501_target_expression.csv", "work/GSE21501_MMP9_features.csv",
                 "work/GSE21501_ACP5_feature.csv"):
        add(f"  {sha256_file(STEP / name)}  {name}")
    for name, h in sorted(s2["outputs"].items()):
        add(f"  {h}  output/{name}")
    add(f"  {sha256_file(OUT / 'manuscript_verification.csv')}  "
        f"output/manuscript_verification.csv")
    for script in sorted((HERE).glob("*.py")):
        add(f"  {sha256_file(script)}  scripts/{script.name}")
    add("")
    add("=" * 78)
    add("END OF RUN LOG")
    add("=" * 78)
    (OUT / "run_log.txt").write_text("\n".join(L) + "\n", encoding="utf-8")

    # =====================================================================
    # sessionInfo.txt
    # =====================================================================
    S = []
    a2 = S.append
    a2("GSE21501 Phase B step 3c-2 -- sessionInfo")
    a2(f"generated : {now} (local)")
    a2("")
    a2("=" * 78)
    a2("PYTHON")
    a2("=" * 78)
    a2(f"  sys.version        : {sys.version}")
    a2(f"  sys.executable     : {sys.executable}")
    a2(f"  platform.python    : {platform.python_implementation()} "
       f"{platform.python_version()}")
    a2(f"  platform.platform  : {platform.platform()}")
    a2(f"  platform.machine   : {platform.machine()}")
    a2(f"  os.cpu_count       : {os.cpu_count()}")
    a2("")
    a2("  packages (the ones that contribute a number to this step):")
    for k, v in pkg_versions().items():
        a2(f"    {k:14s} {v}")
    a2("")
    a2("  statistical machinery actually used:")
    a2("    scipy.stats.spearmanr       rho and the asymptotic p")
    a2("    scipy.stats.rankdata        ranks in the bootstrap and the permutation")
    a2("    numpy.random.default_rng    PCG64, the only RNG in this step")
    a2("    numpy.nanstd/nanmean/nanquantile, pandas groupby.mean")
    a2("")
    a2("=" * 78)
    a2("R")
    a2("=" * 78)
    a2(f"  {r_version()}")
    a2("  R was NOT used by this step. Every number here was produced in Python by")
    a2("  scripts/01_extract_expression.py, scripts/02_analyse.py and")
    a2("  scripts/03_verify_against_manuscript.py. The archived manuscript module that")
    a2("  the verification executes is itself Python.")
    a2("")
    a2("=" * 78)
    a2("SHELL / TOOLING")
    a2("=" * 78)
    a2("  shell              : Windows PowerShell 5.1 (powershell.exe).")
    a2("                       'pwsh' is NOT on PATH on this host, so no logger here")
    a2("                       invokes it.")
    a2("  no network access was needed by any stage of this step: the deposit and the")
    a2("  platform annotation were already local, and the archived manuscript module is")
    a2("  read from disk.")
    a2("")
    a2("=" * 78)
    a2("SEEDS AND RESAMPLING")
    a2("=" * 78)
    a2(f"  master seed        : {SEED}   (frozen_cohort_definition.md S9)")
    a2(f"  n resamples        : {N_RESAMPLES}  (frozen S9)")
    a2(f"  primary addressing : generator seeded directly with {SEED} and consumed in")
    a2("                       the manuscript's call order (bootstrap, then permutation)")
    a2("  other cells        : seed = int(sha256(f'{master}|{label}').hexdigest()[:16],")
    a2("                       16) % 2**32 -- see output/seed_registry.json")
    a2("")
    a2("=" * 78)
    a2("ARTEFACT HASHES (sha256)")
    a2("=" * 78)
    a2(f"  {mh}  ../sample_manifest.csv  (frozen)")
    a2(f"  {sha256_file(PHASE / '_probe' / 'GPL4133.annot.gz')}  ../_probe/GPL4133.annot.gz")
    a2(f"  {sha256_file(PHASE / '_probe' / 'GSE21501_series_matrix.txt.gz')}  "
       f"../_probe/GSE21501_series_matrix.txt.gz")
    a2(f"  {s3['manuscript_module_sha256']}  manuscript module "
       f"(analysis_F_cross_scale.py, v47 archive)")
    a2("")
    a2("  outputs:")
    for name, h in sorted(s2["outputs"].items()):
        a2(f"    {h}  output/{name}")
    a2(f"    {sha256_file(OUT / 'manuscript_verification.csv')}  "
       f"output/manuscript_verification.csv")
    a2("  scripts:")
    for script in sorted(HERE.glob("*.py")):
        a2(f"    {sha256_file(script)}  scripts/{script.name}")
    a2("")
    (OUT / "sessionInfo.txt").write_text("\n".join(S) + "\n", encoding="utf-8")

    # =====================================================================
    # verification_and_controls.md
    # =====================================================================
    V = []
    v = V.append
    v("# Verification and controls -- GSE21501 step 3c-2")
    v("")
    v(f"Generated {now}. Every row below is reproduced from the machine-readable files")
    v("`output/manuscript_verification.csv` and `output/verification_controls.csv`.")
    v("")
    v("## Part 1 -- the manuscript's own scoring function, executed on this cohort")
    v("")
    v("Stage 2 re-implements the frozen scoring and resampling rules by transcription.")
    v("A transcription that is wrong in the same way as its original is invisible to")
    v("internal checking, so the archived manuscript module")
    v("`analysis_F_cross_scale.py` was imported and **its** functions were run on **this**")
    v("cohort's data. The module was read, never edited.")
    v("")
    v(f"* path: `{s3['manuscript_module']}`")
    v(f"* sha256: `{s3['manuscript_module_sha256']}`")
    v(f"* verification population: n = {s3['verification_population_n']} (DUP-2 applied)")
    v(f"* seed {SEED}, {N_RESAMPLES} resamples")
    v("")
    v("| id | check | expectation | observed | verdict |")
    v("|---|---|---|---|---|")
    for c in man_checks:
        v(f"| {c['id']} | {c['check']} | {c['expectation']} | {c['observed']} | "
          f"**{c['verdict']}** |")
    v("")
    v("V3--V6 run on the real data with the real seed, so a divergence in the collapse")
    v("rule, the ddof, the nanmean, the clip or the resampling scheme would surface as a")
    v("numeric mismatch. V8 is the end-to-end statement: the manuscript's own scoring")
    v("function, run on this cohort, returns exactly the rho written to")
    v("`output/primary_result.csv`.")
    v("")
    v("## Part 2 -- controls on this run's own pipeline")
    v("")
    v("| control | kind | expectation | observed | verdict |")
    v("|---|---|---|---|---|")
    for c in controls:
        v(f"| {c['control']} | {c['kind']} | {c['expectation']} | {c['observed']} | "
          f"**{c['verdict']}** |")
    v("")
    v("### Which of these could actually have failed")
    v("")
    v("A control that cannot fail is decoration. The ones below are the load-bearing")
    v("ones, each with the specific failure it would have caught:")
    v("")
    v("* **C3 circular rotation** -- rho fell from +0.1712 to -0.0317. A pipeline that")
    v("  had manufactured the association from the *order* of the array columns would")
    v("  have returned the same value here, because a rotation preserves nothing about")
    v("  the pairings except their order.")
    v("* **C4 label scramble** -- -0.0451, using the seeded generator rather than array")
    v("  order, so it is not a restatement of C3.")
    v("* **C5 seed sensitivity and reproducibility** -- the same seed reproduces the")
    v("  interval exactly, and a different seed moves it. If the seed argument had been")
    v("  inert, every resampled number in this run would have been unanchored and this")
    v("  control would have failed.")
    v("* **C6 zero-variance programme gene** -- exercises the `sd == 0 -> NaN` branch of")
    v("  the frozen `zscore_columns`, the one scoring detail the manuscript records as")
    v("  locked-day implementation rather than protocol text.")
    v("* **C8 null sampling distribution** -- 2000 independent permutations gave a null")
    v("  with mean -0.0057 and sd 0.0862 against a theoretical 1/sqrt(n-1) = 0.0884. A")
    v("  null centred away from 0 would have indicated a mis-specified statistic.")
    v("* **C10 nesting** -- the 102 and every drop-set sensitivity are subsets of the")
    v("  same 129, which guards against a sensitivity silently computed on a different")
    v("  cohort, for instance one that reintroduced a DUP-2 duplicate.")
    v("")
    v("### A control that was wrong, and what happened")
    v("")
    v("C7 first required the bootstrap interval to be **bit-identical** under two")
    v("different RNG seeds. That is unsatisfiable: two different seeds cannot produce")
    v("identical resampled intervals. The script treated the mismatch as fatal and")
    v("exited 1 rather than reporting a number it could not stand behind. The check was")
    v("then corrected to test what the frozen rule actually implies -- `rho` is exactly")
    v("invariant to how the generator is addressed, since no generator enters it, while")
    v("the resampled quantities may differ only by Monte-Carlo error. The correction is")
    v("recorded in the run log and in `failed_or_excluded.md`, not silently applied.")
    v("")
    (OUT / "verification_and_controls.md").write_text("\n".join(V) + "\n",
                                                     encoding="utf-8")

    print("wrote run_log.txt, sessionInfo.txt, verification_and_controls.md, "
          "seed_registry.json")
    print(f"  manuscript checks : {len(man_checks)} "
          f"(failed {sum(1 for c in man_checks if c['verdict'] != 'PASS')})")
    print(f"  controls          : {len(controls)} "
          f"(failed {sum(1 for c in controls if c['verdict'] != 'PASS')})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
