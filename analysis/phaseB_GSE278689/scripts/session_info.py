#!/usr/bin/env python
# =============================================================================
# session_info.py -- Phase B step 2, GSE278689
# Records the exact software environment that produced every result file, and
# the GSM <-> barcode-prefix <-> patient map taken from the GEO series matrix.
# =============================================================================
from __future__ import annotations

import os
import importlib.metadata as md
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(os.environ.get("SIGLEC15_GSE278689_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927\phaseB_GSE278689"))
PROV = ROOT / "provenance"
RSCRIPT = os.environ.get("RSCRIPT_EXE", r"C:\Program Files\R\R-4.5.1\bin\Rscript.exe")

PY_PKGS = ["numpy", "pandas", "scipy", "anndata", "scanpy", "h5py", "statsmodels"]
R_PKGS = ["data.table", "Matrix", "boot", "stats", "Seurat", "SeuratObject", "scran", "scater"]

# from !Sample_title / !Sample_geo_accession of GSE278689_series_matrix.txt.gz,
# joined to the barcode prefix by patient ID (00_audit_inputs.R, verified 1:1)
ACCESSION = [
    ("GSM8561922", "PA40", "2048540_FFPE"),
    ("GSM8561923", "PA21", "2107809_FFPE"),
    ("GSM8561924", "PA85", "1869485_FFPE"),
    ("GSM8561925", "PA98", "2151198_FFPE"),
    ("GSM8561926", "PA72", "2051072_FFPE"),
    ("GSM8561927", "PA91", "1942191_FFPE"),
    ("GSM8561928", "PA29", "2049529_FFPE"),
    ("GSM8561929", "PA24", "2021124_FFPE"),
]


def main() -> None:
    PROV.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    A = lines.append
    A("=" * 78)
    A("sessionInfo -- Phase B step 2, GSE278689 (snRNA) -- cross-modality sensitivity analysis")
    A("=" * 78)
    A(f"generated            : {time.strftime('%Y-%m-%d %H:%M:%S %Z')}")
    A(f"working directory    : {Path.cwd()}")
    A("")
    A("--- operating system ---")
    A(f"platform             : {platform.platform()}")
    A(f"machine / processor  : {platform.machine()} / {platform.processor()}")
    A(f"python               : {sys.version.splitlines()[0]}")
    A(f"python executable    : {sys.executable}")
    A("")
    A("--- python packages ---")
    for p in PY_PKGS:
        try:
            A(f"{p:<21}: {md.version(p)}")
        except md.PackageNotFoundError:
            A(f"{p:<21}: (not installed)")
    A("")
    A("--- R ---")
    try:
        rv = subprocess.run([RSCRIPT, "--version"], capture_output=True, text=True, timeout=60)
        A((rv.stdout or rv.stderr).strip().splitlines()[0])
    except Exception as e:  # noqa: BLE001
        A(f"Rscript not reachable: {e}")
    try:
        code = ("ip <- installed.packages()[, 'Version']; "
                "for (p in c(" + ",".join(f'"{p}"' for p in R_PKGS) + ")) "
                "cat(sprintf('%-21s: %s\\n', p, if (p %in% names(ip)) ip[[p]] else '(not installed)'))")
        rv = subprocess.run([RSCRIPT, "-e", code], capture_output=True, text=True, timeout=180)
        A(rv.stdout.strip() or rv.stderr.strip())
    except Exception as e:  # noqa: BLE001
        A(f"R package query did not run: {e}")
    A("")
    A("--- frozen constants ---")
    A("seed                 : 15072026")
    A("bootstrap resamples  : 10000")
    A("z-score ddof         : 0 (population SD), matching zscore_columns()")
    A("detection            : non-zero RAW deposited count")
    A("normalisation        : CP10K then log1p (deposit verified as raw integer counts)")
    A("programme score      : ACP5 + MMP9 only; SIGLEC15 never enters it")
    A("contrast             : mean(score | detectable) - mean(score | non-detectable)")
    A("sample weighting     : equal")
    A("")
    A("--- sample accession map (GEO series matrix -> barcode prefix) ---")
    A(f"{'GSM':<12}{'patient':<10}{'barcode prefix':<18}{'title':<26}")
    for gsm, pat, pref in ACCESSION:
        A(f"{gsm:<12}{pat:<10}{pref:<18}{pat + ', Tumor, snRNA':<26}")
    A("")
    A("Series title         : Integrated Single-Cell and Spatial Transcriptomics Uncover")
    A("                       Distinct Cellular Subtypes Involved in Neural Invasion in")
    A("                       PDAC [snRNA-Seq]")
    A("SubSeries of         : GSE278694 (SuperSeries; shares patients with GSE278688/GSE278687)")
    A("Platform             : 10x Genomics, Cell Ranger 7.1.0, GRCh38, Illumina NovaSeq 6000")
    A("Modality             : snRNA (single NUCLEI) -- NOT scRNA. This run is a")
    A("                       modality-specific sensitivity analysis and is never pooled.")
    A("=" * 78)

    (PROV / "sessionInfo.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    # the brief also names sessionInfo.txt at the run root
    (ROOT / "sessionInfo.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    acc = ["gsm,patient,sample_id_barcode_prefix,title,tissue,modality"]
    for gsm, pat, pref in ACCESSION:
        acc.append(f"{gsm},{pat},{pref},\"{pat}, Tumor, snRNA\",Tumor,snRNA")
    (PROV / "sample_accession_map.csv").write_text("\n".join(acc) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print("\nwrote provenance/sessionInfo.txt and sessionInfo.txt")
    print("wrote provenance/sample_accession_map.csv")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
