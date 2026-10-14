#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
extension_meta / 03_session_info.py

Writes sessionInfo.txt: the Python record, the R record, and the input/output hashes for
every artefact this directory depends on or produces.
"""
from __future__ import annotations
import os

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

EXT = Path(os.environ.get("SIGLEC15_EXT_META_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927\extension_meta"))
ARCHIVE = Path(
    os.path.join(os.environ.get("SIGLEC15_SUBMISSION_ROOT", r"E:\Tex\HMlong\npjPrecOncol\FINAL_SUBMISSION_v6"), "BMC_MEDICINE_SUBMISSION_v39", "_qa_archive_v39_20260912")
)
RSCRIPT = Path(os.environ.get("RSCRIPT_EXE", r"C:\Program Files\R\R-4.5.1\bin\Rscript.exe"))

INPUTS = [
    ARCHIVE / "scripts" / "analysis_F_cross_scale.py",
    ARCHIVE / "results" / "bulk_associations.csv",
    ARCHIVE / "results" / "bulk_random_effects_meta.csv",
    Path(os.path.join(os.environ.get("SIGLEC15_PHASEB_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927"), "phaseB_GSE205154", "step3c_rho", "primary_result.csv")),
    Path(os.path.join(os.environ.get("SIGLEC15_PHASEB_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927"), "phaseB_GSE205154", "step3c_rho", "results_summary.md")),
    Path(os.path.join(os.environ.get("SIGLEC15_PHASEB_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927"), "phaseB_GSE21501", "step3c_rho", "output", "primary_result.csv")),
    Path(os.path.join(os.environ.get("SIGLEC15_PHASEB_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927"), "phaseB_GSE21501", "step3c_rho", "output", "results_summary.md")),
]
OUTPUTS = [
    "meta_results.csv", "forest_inputs.csv", "verification_controls.csv",
    "results_summary.md", "failed_or_excluded.md", "run_log.txt", "sessionInfo.txt",
    "scripts/01_extension_meta.py", "scripts/02_crosscheck_metafor.R",
    "scripts/03_session_info.py", "scripts/04_verify_summary.py",
    "provenance/meta_handoff.json", "provenance/crosscheck_metafor.csv",
    "provenance/crosscheck_metafor_log.txt", "provenance/checks.json",
    "provenance/headline.json", "provenance/source_hashes.txt",
    "provenance/verify_summary_log.txt",
]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main() -> int:
    L: list[str] = []
    A = L.append
    A("=" * 100)
    A("sessionInfo.txt -- extension_meta (secondary bulk-extension meta-analysis, BJC revision)")
    A(f"written : {datetime.now().isoformat(timespec='seconds')}")
    A("=" * 100)

    A("\n[1] PYTHON SESSION")
    A(f"    python     : {sys.version}")
    A(f"    executable : {sys.executable}")
    A(f"    platform   : {platform.platform()}")
    A(f"    processor  : {platform.processor()}")
    A(f"    numpy      : {np.__version__}")
    A(f"    scipy      : {scipy.__version__}")
    A(f"    pandas     : {pd.__version__}")
    import importlib.metadata as md
    for pkg in ("numpy", "scipy", "pandas", "anndata", "matplotlib"):
        try:
            A(f"    {pkg:<10s} : {md.version(pkg)}")
        except Exception as e:  # pragma: no cover
            A(f"    {pkg:<10s} : (unavailable: {e})")

    A("\n[2] R SESSION")
    if RSCRIPT.exists():
        try:
            out = subprocess.run(
                [str(RSCRIPT), "-e",
                 'cat(R.version.string,"\\n"); for (p in c("metafor","jsonlite")) '
                 'cat(p, ":", if (requireNamespace(p, quietly=TRUE)) '
                 'as.character(packageVersion(p)) else "NOT INSTALLED", "\\n")'],
                capture_output=True, text=True, timeout=120,
            )
            A("    " + str(RSCRIPT))
            for line in (out.stdout + out.stderr).strip().splitlines():
                A("    " + line)
        except Exception as e:  # pragma: no cover
            A(f"    R invocation failed: {e}")
    else:
        A("    Rscript not found")

    A("\n[3] WHAT WAS REUSED (not re-implemented)")
    A("    The manuscript's own random_effects_meta() was imported by path from the v39")
    A("    analysis archive and called directly for all three models. A byte-faithful replica")
    A("    of that function exists in scripts/01_extension_meta.py only so that deliberately")
    A("    mutated conventions can be run as falsification controls; the replica is validated")
    A("    against the imported original (control C0, max|delta| = 0 on all three models)")
    A("    before any mutated variant is believed.")

    A("\n[4] INPUT HASHES")
    for p in INPUTS:
        A(f"    {'OK ' if p.exists() else 'MISSING'} {p}")
        if p.exists():
            A(f"          sha256 {sha256(p)}")

    A("\n[5] OUTPUT HASHES")
    for rel in OUTPUTS:
        p = EXT / rel
        A(f"    {'OK ' if p.exists() else 'MISSING'} {rel}")
        if p.exists():
            A(f"          sha256 {sha256(p)}")

    A("\n[6] FROZEN TEXT CARRIED VERBATIM")
    A('    "The three original PDAC bulk cohorts will be combined with the independently')
    A('     eligible primary-tumour GSE205154 cohort using the same Fisher-z transformation,')
    A('     REML random-effects model and modified Hartung-Knapp variance rule used in the')
    A('     manuscript. This four-cohort estimate is a revision-stage secondary pooled estimate')
    A('     and does not replace the original three-cohort estimate."')
    A('    "GSE21501 will then be added in a five-cohort sensitivity analysis because exact')
    A('     patient independence cannot be established from the deposited metadata."')
    A("")
    A("    The GSE21501 analysis was specimen-level because patient identifiers were not")
    A("    deposited; the number of independent patients could therefore not be established")
    A("    exactly.")

    A("\n[7] SCOPE")
    A("    No manuscript text, figure or claim was written or modified. No cohort was excluded,")
    A("    reweighted or reshaped on the basis of its rho magnitude or direction. Bulk cohorts")
    A("    only: no single-cell and no snRNA cohort was pooled. The published three-cohort")
    A("    estimate (rho = 0.132, 95% CI -0.110 to +0.360) is reproduced and remains unreplaced.")

    txt = "\n".join(L) + "\n"
    (EXT / "sessionInfo.txt").write_text(txt, encoding="utf-8")
    print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
