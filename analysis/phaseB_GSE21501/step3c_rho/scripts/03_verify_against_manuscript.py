#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE21501 -- Phase B step 3c-2, stage 3 of 3: VERIFY AGAINST THE MANUSCRIPT ITSELF.

The point of this file
----------------------
Stage 2 re-implements the frozen scoring and resampling rules by transcription.  A
transcription that is wrong in the SAME way as the thing it transcribes is invisible to
any amount of internal checking.  The one check that is not self-confirming is to load
the manuscript's own archived module and run **its** functions on **this** cohort, then
require agreement with the stage-2 numbers.

    manuscript module:
      npjPrecOncol/FINAL_SUBMISSION_v6/BMC_MEDICINE_SUBMISSION_v47/
        _archives_translate/Additional_file_2_Analysis_archive_v39/
        scripts/analysis_F_cross_scale.py

That module is imported, never edited.  It has no side effects on import beyond
defining paths and a module-level RNG, and this script never calls its `main()`.

What is checked
---------------
  V1  the archived module loads, and its SHA-256 is recorded
  V2  the four frozen functions exist in it with the expected signatures
  V3  `zscore_columns` agrees cell-for-cell on the real 129 x 3 matrix
  V4  `score_program` produces the identical score vector
  V5  `spearman_bootstrap` returns the identical (rho, p, lo, hi) given the same seed
  V6  `permutation_p` returns the identical p given the same seed
  V7  `required_gene_count(PRIMARY) == 2`, i.e. the two-gene rule the frozen text names

A control that could have failed: V3-V6 are run on the real data with the real seeds, so
any divergence in the collapse rule, the ddof, the nanmean, the clip, or the resampling
scheme would surface as a numeric mismatch rather than as a passing check.

Output
    output/manuscript_verification.csv
    work/stage03_qc.json
"""

from __future__ import annotations
import os

import hashlib
import importlib.util
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

CHOSEN_HOME = Path(
    os.path.join(os.environ.get("SIGLEC15_SUBMISSION_ROOT", r"E:\Tex\HMlong\npjPrecOncol\FINAL_SUBMISSION_v6"), "BMC_MEDICINE_SUBMISSION_v47"))
ARCHIVE_REL = Path("_archives_translate/Additional_file_2_Analysis_archive_v39/"
                   "scripts/analysis_F_cross_scale.py")

SEED = 15072026
N_RESAMPLES = 10_000
TARGET = "SIGLEC15"
MANIFEST = PHASE / "sample_manifest.csv"
FROZEN_SHA256_MANIFEST = "3e11938ae5c7e682c9e53aad93d051579764a4ec8173d69d571dbeda0f3ac4e4"


def sha256_file(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            d.update(chunk)
    return d.hexdigest()


def load_module(path: Path, name: str):
    """Import a file by path.  Neither module runs main() on import (both guard it)."""
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_manuscript_module(path: Path):
    return load_module(path, "manuscript_analysis_F")


def main() -> int:
    print("=" * 78)
    print("GSE21501 step 3c-2 -- stage 3 VERIFY AGAINST THE MANUSCRIPT MODULE")
    print("=" * 78)

    checks = []

    def record(vid, what, expectation, observed, verdict, note=""):
        checks.append({"id": vid, "check": what, "expectation": expectation,
                       "observed": observed, "verdict": verdict, "note": note})
        print(f"[{verdict:4}] {vid} {what}: {observed}")

    # -- locate the module, preferring the frozen document's own home tree ----
    cand = CHOSEN_HOME / ARCHIVE_REL
    found_from = "frozen document's home tree (BMC_MEDICINE_SUBMISSION_v47)"
    if not cand.exists():
        cand = None
        for p in sorted(PHASE.parents[0].rglob("analysis_F_cross_scale.py")):
            if "Additional_file_2_Analysis_archive" in str(p):
                cand = p
                found_from = "fallback rglob"
                break
    if cand is None or not cand.exists():
        record("V1", "archived manuscript module located", "a file exists",
               "NOT FOUND", "FAIL",
               "the manuscript's own scoring function is unavailable, so this run's "
               "verification would rest entirely on transcription")
        pd.DataFrame(checks).to_csv(OUT / "manuscript_verification.csv", index=False)
        (WORK / "stage03_qc.json").write_text(json.dumps(
            {"stage": "03_verify_manuscript", "checks": checks,
             "n_failed": 1}, indent=2), encoding="utf-8")
        return 1

    man_hash = sha256_file(cand)
    mod = load_manuscript_module(cand)
    record("V1", "archived manuscript module loaded",
           "imports without calling main()",
           f"ok, sha256={man_hash[:16]}... from {found_from}", "PASS",
           f"full path: {cand}")

    # -- V2 signatures --------------------------------------------------------
    have = {n: hasattr(mod, n) for n in
            ("zscore_columns", "score_program", "required_gene_count",
             "spearman_bootstrap", "permutation_p")}
    record("V2", "the four frozen functions exist in the archived module",
           "all present", str(have),
           "PASS" if all(have.values()) else "FAIL",
           "these are the functions the frozen document names in S9")

    # -- data -----------------------------------------------------------------
    man = pd.read_csv(MANIFEST)
    expr = pd.read_csv(WORK / "GSE21501_target_expression.csv")
    df = man.merge(expr[["gsm_accession", TARGET, "ACP5", "MMP9"]],
                   on="gsm_accession", how="left", validate="one_to_one")
    df["is_dup_dropped"] = df["exclusion_reason"].fillna("").str.contains("EX-5",
                                                                         regex=False)
    dep = df.loc[~df["is_dup_dropped"]].reset_index(drop=True)
    print(f"[info] verification population: n = {len(dep)} (DUP-2 applied)")
    if len(dep) != 129:
        raise SystemExit(f"FATAL: expected 129, got {len(dep)}")

    frame = dep[[TARGET, "ACP5", "MMP9"]].copy()
    raw = frame.to_numpy(float)

    # -- V3 zscore_columns cell-for-cell -------------------------------------
    stage2 = load_module(HERE / "02_analyse.py", "stage2_analyse")
    z_mine = stage2.zscore_columns
    z_ms = mod.zscore_columns(raw)
    z_me = z_mine(raw)
    d_z = float(np.nanmax(np.abs(z_ms - z_me)))
    record("V3", "zscore_columns agrees cell-for-cell on the real 129 x 3 matrix",
           "max abs difference == 0", f"max|diff| = {d_z:.3e}",
           "PASS" if d_z == 0.0 else "FAIL",
           "covers ddof, the sd==0 -> NaN branch, and the centring")

    # -- V4 score_program identical ------------------------------------------
    sp_mine = stage2.score_program
    s_ms, pres_ms = mod.score_program(frame, "Derosiers_ACP5_MMP9_core")
    s_me, pres_me = sp_mine(frame, "Derosiers_ACP5_MMP9_core")
    d_s = float(np.nanmax(np.abs(s_ms - s_me)))
    same_pres = list(pres_ms) == list(pres_me)
    record("V4", "score_program produces the identical score vector",
           "max abs difference == 0 and identical gene list",
           f"max|diff| = {d_s:.3e}; genes {pres_ms} vs {pres_me}",
           "PASS" if (d_s == 0.0 and same_pres) else "FAIL",
           "this is the ACP5/MMP9 score the primary rho is computed against")

    # -- V5 spearman_bootstrap identical under one shared seed ---------------
    # The archived function has signature (x, y) and draws from the MODULE-LEVEL
    # `RNG`.  Rather than rewrite it, this script reaches into the module and seats a
    # generator it controls.  The archived code path is therefore executed verbatim --
    # only its entropy source is supplied from outside.
    mod.N_RESAMPLES = N_RESAMPLES
    stage2.N_RESAMPLES = N_RESAMPLES
    x = dep[TARGET].to_numpy(float)

    mod.RNG = np.random.default_rng(SEED)
    ms_rho, ms_p, ms_lo, ms_hi = mod.spearman_bootstrap(x, s_ms)
    me_rho, me_p, me_lo, me_hi = stage2.spearman_bootstrap(
        x, s_me, np.random.default_rng(SEED))
    ident = (ms_rho == me_rho and ms_p == me_p and ms_lo == me_lo and ms_hi == me_hi)
    record("V5", "spearman_bootstrap identical given the same seed",
           "rho, asymptotic p and both CI limits agree exactly",
           f"rho {ms_rho:+.12f} vs {me_rho:+.12f}; CI [{ms_lo:+.9f},{ms_hi:+.9f}] vs "
           f"[{me_lo:+.9f},{me_hi:+.9f}]",
           "PASS" if ident else "FAIL",
           "the archived function was executed as written, with its module-level RNG "
           "seated to the frozen SEED; the stage-2 transcription was given an "
           "independently constructed generator with the same seed")

    # -- V6 permutation_p identical ------------------------------------------
    mod.RNG = np.random.default_rng(SEED)
    ms_pp = mod.permutation_p(x, s_ms, ms_rho)
    me_pp = stage2.permutation_p(x, s_me, me_rho, np.random.default_rng(SEED))
    record("V6", "permutation_p identical given the same seed",
           "p agrees exactly", f"{ms_pp:.10f} vs {me_pp:.10f}",
           "PASS" if ms_pp == me_pp else "FAIL",
           "independent of the bootstrap: a separate 10 000-draw generator")

    # -- V7 the two-gene rule -------------------------------------------------
    rgc = int(mod.required_gene_count("Derosiers_ACP5_MMP9_core"))
    record("V7", "required_gene_count(PRIMARY) == 2", "== 2", str(rgc),
           "PASS" if rgc == 2 else "FAIL",
           "the frozen text's 'both genes must be present' rule")

    # -- V8 does the archived module's own PRIMARY value reproduce stage 2? ----
    stage2_qc = json.loads((WORK / "stage02_qc.json").read_text(encoding="utf-8"))
    prim_rho = stage2_qc["primary"]["rho"]
    record("V8", "the archived module reproduces the stage-2 primary rho",
           "identical to the value written to primary_result.csv",
           f"manuscript-module rho = {ms_rho:+.12f} vs stage-2 rho = {prim_rho:+.12f}",
           "PASS" if ms_rho == prim_rho else "FAIL",
           "This is the end-to-end check: the manuscript's own scoring function, run "
           "on this cohort's data, returns exactly the number in primary_result.csv.")

    # -- manifest hash, re-asserted at the end of the run ---------------------
    mh = sha256_file(MANIFEST)
    record("V9", "frozen manifest hash still matches after the whole run",
           "== the frozen value", mh,
           "PASS" if mh == FROZEN_SHA256_MANIFEST else "FAIL",
           "no step in this run wrote to sample_manifest.csv")

    dfc = pd.DataFrame(checks)
    dfc.to_csv(OUT / "manuscript_verification.csv", index=False)
    qc = {"stage": "03_verify_manuscript",
          "manuscript_module": str(cand),
          "manuscript_module_sha256": man_hash,
          "manuscript_module_found_from": found_from,
          "verification_population_n": int(len(dep)),
          "seeds_used": {"bootstrap": SEED, "permutation": SEED},
          "n_resamples": N_RESAMPLES,
          "checks": checks,
          "n_checks": len(checks),
          "n_failed": int((dfc["verdict"] != "PASS").sum())}
    (WORK / "stage03_qc.json").write_text(json.dumps(qc, indent=2, default=str),
                                          encoding="utf-8")
    print("\n" + "=" * 78)
    print(f"manuscript-module checks: {len(checks)}  failed: {qc['n_failed']}")
    print("=" * 78)
    return 0 if qc["n_failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
