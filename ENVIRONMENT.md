# Software environment

This note separates two things that are easy to conflate: the environment that **produced** the reported results, and the environment a third party needs to **re-run** the code. They are not necessarily the same, and the original versions below are transcribed from the `sessionInfo.txt` records in this deposit rather than restated from a current machine.

## Original environment

As recorded at run time by the scripts that wrote each `sessionInfo.txt`:

| Component | Version | Recorded in |
|---|---|---|
| Python | 3.13.0 (MSC v.1941 64 bit, AMD64) | `analysis/phaseB_GSE205154/sessionInfo.txt`, `analysis/phaseB_GSE278689/sessionInfo.txt`, `analysis/phaseB_GSE278689/provenance/sessionInfo.txt` |
| Platform | Windows-11-10.0.26200-SP0 | as above |
| numpy | 2.4.6 | as above |
| pandas | 2.3.3 | as above |
| scipy | 1.18.0 | as above |
| R | 4.5.1 (2025-06-13 ucrt), x86_64-w64-mingw32/x64 | `analysis/phaseB_GSE205154/step3c_rho/output/sessionInfo_R.txt`, `analysis/phaseB_GSE278688/sessionInfo.txt` |
| Rscript path used | `C:\Program Files\R\R-4.5.1\bin\Rscript.exe` | overridable via `RSCRIPT_EXE` |

The full package-level listings — including the Bioconductor packages pulled in for the R cross-checks — are in those `sessionInfo.txt` files, which are reproduced unedited apart from the user-name substitution noted in the root README.

Third-party Python packages imported by the code, with the number of files that import each. The list was produced by walking the import statements in the deposit rather than from memory, so it is complete for this code:

| Package | Files | Used for |
|---|---:|---|
| `pandas` | 44 | tabular handling throughout |
| `numpy` | 40 | numerical work throughout |
| `scipy` | 30 | correlation, distributions, statistics |
| `anndata` | 8 | single-cell object handling |
| `matplotlib` | 7 | figure generation |
| `statsmodels` | 1 | Benjamini–Hochberg correction in the cross-target calibration panel (`analysis/scripts/cross_target_panel/run_cross_target.py`) |
| `h5py` | 1 | HDF5 access |
| `pyreadr` | 1 | reading R data files |

R packages called via `library()` / `require()` across the R scripts:

`ComplexHeatmap`, `data.table`, `dplyr`, `ggplot2`, `ggrepel`, `grDevices`, `grid`, `jsonlite`, `metafor`, `patchwork`, `ragg`, `readr`, `scales`, `svglite`, `tidyr`

with the Bioconductor packages pulled in for the cross-checks (`Matrix`, `MatrixGenerics`, `BiocManager`) recorded per run in the R session records.

## Reproduction environment

There is no lockfile. The original run predates this deposit, and reconstructing exact transitive pins for a Python 3.13 / R 4.5 mixed environment after the fact would produce a file that looks authoritative without being verifiable. Rather than publish invented pins, the guidance is:

1. Match the **major/minor** versions of Python and R above. Patch-level differences were not tested by re-running the analysis, so no claim is made here that results are invariant across patch releases. The one place where a version change is known to be able to move a reported interval is the `metafor` cross-check: install the version recorded in `analysis/extension_meta/provenance/crosscheck_metafor_log.txt`, and treat any interval computed under a different `metafor` version as needing its own verification.
2. Install the Python packages listed above and let the resolver pick compatible versions.
3. Install the R packages named in the relevant `sessionInfo_R.txt` / `sessionInfo.txt`.
4. Verify your environment reproduces the recorded one by re-running the relevant `session_info.py` and diffing against the committed record.

If you need a hard pin set, generate it from your own resolved environment and treat it as *your* reproduction environment, not as a reconstruction of the original.

## Expected input layout

The scripts locate their inputs through the environment variables documented in the root README, each defaulting to the path used in the original run. The layout those defaults describe is:

```
<project root>/                          # SIGLEC15_PROJECT_ROOT
└── BJC_EXT_v1_20260927/                 # SIGLEC15_PHASEB_ROOT
    ├── annotations/                     # per-dataset metadata files
    │   ├── GSE278688_sc_metadata.csv.gz
    │   └── GSE278689_sn_metadata.csv.gz
    ├── phaseB_GSE205154/
    │   └── step3c_rho/
    │       ├── primary_result.csv
    │       └── results_summary.md
    ├── phaseB_GSE21501/
    │   └── step3c_rho/output/
    │       ├── primary_result.csv
    │       └── results_summary.md
    ├── phaseB_GSE278688/                # SIGLEC15_GSE278688_ROOT
    ├── phaseB_GSE278689/                # SIGLEC15_GSE278689_ROOT
    └── extension_meta/                  # SIGLEC15_EXT_META_ROOT

<submission root>/                       # SIGLEC15_SUBMISSION_ROOT
├── biological_validation/analysis_F_cross_scale/
└── BMC_MEDICINE_SUBMISSION_v47/         # external single-cell gate outputs
```

Set the variables to your own locations and the layout is free-form; the tree above is only what the defaults assume.

### Cross-target calibration panel

`analysis/scripts/cross_target_panel/run_cross_target.py` does not use the variables above. It walks four directories up from its own location to find a project root, then reads expression matrices from `<root>/data/analysis/` and clinical tables from `<root>/data/TCGA/`. In the original tree the script lived at `<root>/data/analysis/cross_target_panel/`, so that walk landed on `<root>`. It does not live there in this deposit, and those two input directories are not redistributed.

To run it, set `SIGLEC15_PROJECT_ROOT` to a tree containing `data/analysis/` and `data/TCGA/`:

```
<SIGLEC15_PROJECT_ROOT>/
├── data/
│   ├── analysis/     # expression matrices (TCGA-PAAD and the panel inputs)
│   └── TCGA/         # clinical tables
```

When the variable is unset the original four-directory walk is used unchanged, so the script behaves in the original tree exactly as it did. Outputs are written beside the script itself.

## What has not been verified

The code compiles, resolves its paths, and every path expression was evaluated with the variables unset to confirm it reproduces the original location. The analysis has **not** been re-executed end to end from public source data in this deposit, because the expression matrices are not redistributed here. Text-and-code agreement is not the same as independent numerical reproduction, and this deposit should not be read as claiming the latter.
