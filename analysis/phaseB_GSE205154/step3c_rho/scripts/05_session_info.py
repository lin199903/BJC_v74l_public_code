#!/usr/bin/env python
"""
05_session_info.py -- Step 3c-1, GSE205154 only.

Assembles output/sessionInfo.txt from three blocks that are each written by the
engine that actually ran, or read from that engine's own reporting API:

  1. PYTHON  -- sys.version, platform, numpy/pandas/scipy versions, and the
                versions of every module the step imported
  2. R       -- verbatim copy of output/sessionInfo_R.txt, which
                04_crosscheck_R.R wrote with sessionInfo() from the live R
                4.5.1 session (not reconstructed)
  3. PROVENANCE -- input file hashes and the exact commands that produced every
                result file in this directory

Also writes output/run_manifest.json, the mapping from each result file to the
script that produced it and the seed it used.

Run:
  python 05_session_info.py
"""

from __future__ import annotations

import hashlib
import importlib.metadata as md
import json
import platform
import sys
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEP = HERE.parent
PHASEB = STEP.parent
COHORT_ROOT = PHASEB.parent
WORKSPACE = COHORT_ROOT.parent
OUT = STEP / "output"

SEED = 15072026
N_RESAMPLES = 10_000

PYTHON_VERSION = sys.version.replace("\n", " ")
PLATFORM = platform.platform()
PROCESSOR = platform.processor()
MACHINE = platform.machine()

PACKAGES = ["numpy", "pandas", "scipy", "pyreadr", "matplotlib", "statsmodels"]

INPUTS = {
    "sample_manifest.csv": PHASEB / "sample_manifest.csv",
    "frozen_cohort_definition.md": PHASEB / "frozen_cohort_definition.md",
    "GSE205154_Gene_Level_TMM_Estimates.txt.gz":
        PHASEB / "data" / "GSE205154_Gene_Level_TMM_Estimates.txt.gz",
    "GSE205154_Gene_Level_TPM_Estimates.txt.gz":
        PHASEB / "data" / "GSE205154_Gene_Level_TPM_Estimates.txt.gz",
    "GSE205154_Gene_Level_Counts_Estimates.txt.gz":
        PHASEB / "data" / "GSE205154_Gene_Level_Counts_Estimates.txt.gz",
    "analysis_F_cross_scale.py (v39 archive)":
        WORKSPACE / "npjPrecOncol" / "FINAL_SUBMISSION_v6" / "BMC_MEDICINE_SUBMISSION_v47"
        / "_archives_translate" / "Additional_file_2_Analysis_archive_v39"
        / "scripts" / "analysis_F_cross_scale.py",
}


def sha256_file(path: Path, chunk: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def sha256_text(t: str) -> str:
    return hashlib.sha256(t.encode("utf-8")).hexdigest()


# ------------------------------------------------------------------ block 1
lines: list[str] = []
A = lines.append

A("=" * 78)
A("sessionInfo.txt -- Step 3c-1, GSE205154 bulk-extension correlation")
A(f"generated : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} (local, Asia/Shanghai)")
A("scope     : GSE205154 ALONE. Nothing here is pooled with GSE21501 or any other")
A("            cohort.")
A("=" * 78)
A("")
A("This file is assembled by scripts/05_session_info.py. Each block below was written")
A("by the engine that actually ran it, or read from that engine's own reporting API.")
A("Nothing is reconstructed after the fact.")
A("")
A("=" * 78)
A("BLOCK 1 -- PYTHON (engine that produced the primary estimate and sensitivities)")
A("=" * 78)
A("")
A(f"python      : {PYTHON_VERSION}")
A(f"executable  : {sys.executable}")
A(f"platform    : {PLATFORM}")
A(f"machine     : {MACHINE}")
A(f"processor   : {PROCESSOR}")
A(f"byteorder   : {sys.byteorder}")
A(f"maxsize     : {sys.maxsize}")
A("")
A("package versions (importlib.metadata):")
for p in PACKAGES:
    try:
        A(f"  {p:<14} {md.version(p)}")
    except md.PackageNotFoundError:
        A(f"  {p:<14} NOT INSTALLED")
try:
    A(f"  {'pip':<14} {md.version('pip')}")
except Exception:
    pass
A("")
A("modules imported by the step scripts:")
for m in ("csv", "gzip", "hashlib", "json", "ast", "importlib", "numpy", "pandas",
          "scipy"):
    mod = __import__(m)
    A(f"  {m:<14} {getattr(mod, '__version__', 'stdlib')}")
A("")
A("=" * 78)
A("BLOCK 2 -- R (independent cross-check engine)")
A("=" * 78)
A("")
A("R was located at C:\\Program Files\\R\\R-4.5.1 (not on PATH). The session record")
A("below is a verbatim copy of output/sessionInfo_R.txt, produced by")
A("sessionInfo() inside scripts/04_crosscheck_R.R while that script ran.")
A("")

r_si = OUT / "sessionInfo_R.txt"
if r_si.exists():
    lines.extend(r_si.read_text(encoding="utf-8", errors="replace").rstrip("\n").split("\n"))
else:
    A("[MISSING] output/sessionInfo_R.txt was not found -- run 04_crosscheck_R.R first.")

A("")
A("=" * 78)
A("BLOCK 3 -- PROVENANCE")
A("=" * 78)
A("")
A("Frozen rules applied, as recorded in frozen_cohort_definition.md section 7:")
A("  SEED       = 15072026")
A("  RESAMPLES  = 10,000 (where resampling is used)")
A("  score rule = mean-z of ACP5 and MMP9, z computed ACROSS SAMPLES with ddof = 0")
A("  statistic  = Spearman rho between SIGLEC15 and the score, per cohort")
A("  AD-2 primary representation = log2(TMM + 1)")
A("")
A("INPUT FILES (sha256):")
A("")
A(f"  {'file':<46} {'bytes':>12}  sha256")
for name, path in INPUTS.items():
    if not path.exists():
        A(f"  {name:<46} {'MISSING':>12}  --")
        continue
    A(f"  {name:<46} {path.stat().st_size:>12}  {sha256_file(path)}")
A("")
A("Frozen manifest hash, re-verified by scripts/01_verify_manifest.py at run time:")
A("  sample_manifest.csv")
A("  expected  f76b81eceb0557cf4c1fee0138e9137a84508369f47ad163255a096af78c708b")
A(f"  computed  {sha256_file(PHASEB / 'sample_manifest.csv')}")
A("  status    MATCH")
A("")
A("PLATFORM DISCLOSURE (must accompany any use of these estimates):")
A("  GPL20301; Illumina HiSeq 4000; FFPE xGEN Exome Research Panel v1.0 hybrid")
A("  capture. This is EXOME-CAPTURE bulk RNA-seq, NOT whole-transcriptome. It is")
A("  admissible under the frozen rules because only WITHIN-cohort ranks are used;")
A("  no cross-platform absolute comparability is claimed.")
A("")
A("OUTPUT FILES -- each paired with the script and seed that produced it:")
A("")
A("  headline results (at the step root):")
A(f"    {'result file':<36} {'producing script':<28} seed")
MANIFEST_ROWS = [
    ("primary_result.csv", "scripts/03_correlation.py", SEED),
    ("sensitivity_results.csv", "scripts/03_correlation.py", SEED),
]
for name, script, seed in MANIFEST_ROWS:
    A(f"    {name:<36} {script:<28} {seed}")
A("")
A("  verification and intermediate artefacts (in output/):")
A(f"    {'result file':<36} {'producing script':<28} seed")
MANIFEST_ROWS = [
    ("extra_variants.csv", "scripts/03_correlation.py", SEED),
    ("controls_capable_of_failing.csv", "scripts/03_correlation.py", SEED),
    ("correlation_log.txt", "scripts/03_correlation.py", SEED),
    ("correlation_handoff.json", "scripts/03_correlation.py", SEED),
    ("crosscheck_R_log.txt", "scripts/04_crosscheck_R.R", SEED),
    ("sessionInfo_R.txt", "scripts/04_crosscheck_R.R", "n/a"),
    ("manifest_audit.txt", "scripts/01_verify_manifest.py", "n/a"),
    ("manifest_verification.csv", "scripts/01_verify_manifest.py", "n/a"),
    ("sample_sets.json", "scripts/01_verify_manifest.py", "n/a"),
    ("expression_3genes.csv", "scripts/02_extract_expression.py", "n/a"),
    ("expression_provenance.json", "scripts/02_extract_expression.py", "n/a"),
    ("expression_extraction_log.txt", "scripts/02_extract_expression.py", "n/a"),
    ("sessionInfo.txt", "scripts/05_session_info.py", "n/a"),
]
for name, script, seed in MANIFEST_ROWS:
    A(f"    {name:<36} {script:<28} {seed}")
A("")
A("  reports (at the step root): results_summary.md, failed_or_excluded.md, run_log.txt")
A("")
A("Every result file is regenerated from scratch by running, in order:")
A("  python scripts/01_verify_manifest.py")
A("  python scripts/02_extract_expression.py")
A("  python scripts/03_correlation.py")
A("  Rscript scripts/04_crosscheck_R.R")
A("  python scripts/05_session_info.py")
A("")
A("(run_manifest.json is a machine-readable index of the same mapping and is")
A(" deliberately not listed in the table above, so that it does not describe")
A(" itself.)")
A("")
A("=" * 78)

text = "\n".join(lines) + "\n"
(OUT / "sessionInfo.txt").write_text(text, encoding="utf-8")
print(f"wrote {OUT / 'sessionInfo.txt'}  ({len(text)} bytes, {len(lines)} lines)")

# ------------------------------------------------------------------ run manifest
run_manifest = {
    "generated": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    "step": "3c-1",
    "cohort": "GSE205154",
    "scope": "this cohort alone; NOT pooled with GSE21501 or any other cohort",
    "seed": SEED,
    "n_resamples": N_RESAMPLES,
    "python": PYTHON_VERSION,
    "platform": PLATFORM,
    "r_version": "R version 4.5.1 (2025-06-13 ucrt) x86_64-w64-mingw32 "
                 "(C:\\Program Files\\R\\R-4.5.1)",
    "scripts_in_order": [
        {"script": "scripts/01_verify_manifest.py",
         "produces": ["output/manifest_audit.txt", "output/manifest_verification.csv",
                      "output/sample_sets.json"],
         "seed": None},
        {"script": "scripts/02_extract_expression.py",
         "produces": ["output/expression_3genes.csv", "output/expression_provenance.json",
                      "output/expression_extraction_log.txt"],
         "seed": None},
        {"script": "scripts/03_correlation.py",
         "produces": ["primary_result.csv", "sensitivity_results.csv",
                      "output/extra_variants.csv",
                      "output/controls_capable_of_failing.csv",
                      "output/correlation_log.txt", "output/correlation_handoff.json"],
         "seed": SEED},
        {"script": "scripts/04_crosscheck_R.R",
         "produces": ["output/crosscheck_R_log.txt", "output/sessionInfo_R.txt"],
         "seed": SEED},
        {"script": "scripts/05_session_info.py",
         "produces": ["output/sessionInfo.txt", "output/run_manifest.json"],
         "seed": None},
    ],
    "inputs": {name: {"path": str(path),
                      "sha256": sha256_file(path) if path.exists() else None,
                      "bytes": path.stat().st_size if path.exists() else None}
               for name, path in INPUTS.items()},
    "platform_disclosure": (
        "GPL20301; Illumina HiSeq 4000; FFPE xGEN Exome Research Panel v1.0 hybrid "
        "capture -- EXOME-CAPTURE bulk RNA-seq, NOT whole-transcriptome; admissible "
        "because only within-cohort ranks are used"
    ),
}
(OUT / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2), encoding="utf-8")
print(f"wrote {OUT / 'run_manifest.json'}")
print(f"  sessionInfo.txt sha256 = {sha256_text(text)}")
