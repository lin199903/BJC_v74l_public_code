# SIGLEC15 cross-scale measurement — analysis code

Analysis code, frozen analysis rules and software-environment records supporting:

> **Cross-scale evaluation of SIGLEC15 RNA for measurement of an ACP5/MMP9-high myeloid state in pancreatic ductal adenocarcinoma**
> Junyi Lin, Yan Huang, Molong Hou

This repository is the public companion to the manuscript's Code availability statement. It holds the **code**; the machine-readable result tables, provenance dumps and integrity inventories are published with the paper as Additional files 2–4 and are not duplicated here.

## Contents

| Path | What it holds |
|---|---|
| `analysis/scripts/` | The primary cross-scale analysis and its prespecified sensitivity analyses (`analysis_F_cross_scale.py`, `run_S18_S19.py`, `run_S20*`, `run_S21*`, `run_S22_dilution.py`, `run_S23_meta_sensitivity.py`, `run_all_validation.py`, `run_P22_scorecard.py`), the coupling audits, and the figure-generation code (`generate_figures_v24.R`, `generate_figures_v34.R`, `fig1_*`, `fig2_v34.py`). |
| `analysis/v72_burden/` | The v72 tumour-level burden and sparsity extension: its frozen estimand contract and the analysis script. |
| `analysis/V04_V05_state_restricted/` | The state-restricted analyses and their claim-ceiling audit, with the analysis specification. |
| `analysis/phaseB_GSE205154/`, `analysis/phaseB_GSE21501/` | Bulk extension cohorts: per-dataset scripts, and the frozen cohort definition that pins the eligibility rules for each dataset. |
| `analysis/phaseB_GSE278688/`, `analysis/phaseB_GSE278689/` | Cellular sensitivity datasets, analysed separately and not pooled with the bulk cohorts: GSE278688 is whole-cell single-cell RNA-seq (`_sc_*` matrices) and GSE278689 is single-nucleus RNA-seq (`_sn_*` matrices). Input audits, eligibility checks, the final analysis and the independent verification scripts. |
| `analysis/extension_meta/scripts/` | The four-cohort and five-cohort random-effects meta-analysis, its `metafor` cross-check in R, and the summary verifier. |
| `paired_representation/` | The V46 paired-representation extension: gate scripts and the gate report. |
| `external_single_cell/` | The external single-cell assessment: the frozen analysis rules, the eligibility and independent-check scripts, and the input inventory. |

Documents that fix the rules are kept alongside the code that implements them — `analysis/scripts/LOCKED_PROTOCOL_2026-07-22.md`, `analysis/scripts/ANALYSIS_SPECIFICATION_AMENDMENTS.md`, `analysis/v72_burden/V72_BURDEN_ESTIMAND_CONTRACT.md`, the per-dataset `frozen_cohort_definition.md` files, `external_single_cell/FROZEN_ANALYSIS_RULES.md` and `paired_representation/V46_paired_representation_extension/README_V46_EXTENSION.md`. Analysis rules and results were separated in time by design; this repository preserves the documents that date the rules.

Each analysis directory that produced a run record carries its `sessionInfo.txt`, `SOFTWARE_VERSIONS.txt` or `ENVIRONMENT.txt` capturing the software versions used; the cross-check runs under `analysis/phaseB_GSE205154/` also carry the R-side record. Three directories — `analysis/scripts/`, `analysis/v72_burden/` and `analysis/V04_V05_state_restricted/` — have no per-directory record of their own; the environment for those runs is documented in `ENVIRONMENT.md` and in the session records of the Phase B runs that shared the same interpreter.

## Notes on the code as run

The scripts are reproduced as they were run, rather than refactored, so that the code that produced the reported results is the code you can read. Two things follow.

**Running the scripts.** The original run placed its inputs and outputs at absolute paths on one machine. Those paths are now resolved through environment variables, each defaulting to the location used in the original run, so the scripts run unchanged on that layout and can be pointed anywhere else by setting the variable. The variables are:

| Variable | Resolves |
|---|---|
| `SIGLEC15_GSE278688_ROOT` | Root of the `phaseB_GSE278688` analysis directory |
| `SIGLEC15_GSE278689_ROOT` | Root of the `phaseB_GSE278689` analysis directory |
| `SIGLEC15_EXT_META_ROOT` | Root of the `extension_meta` analysis directory |
| `SIGLEC15_PHASEB_ROOT` | Root of the Phase B extension tree, which holds the `annotations/` metadata files |
| `SIGLEC15_SUBMISSION_ROOT` | Root of the submission tree, where the external single-cell gate outputs are written |
| `SIGLEC15_PROJECT_ROOT` | Project root, used by the v72 burden extension |
| `RSCRIPT_EXE` | Full path to `Rscript.exe`; most scripts fall back to `Rscript` on `PATH` |

Only these path literals were changed. Each substituted literal was restored and the resulting source compared to the original as an abstract syntax tree; all 65 Python files compile and are logically identical to the code that was run. The eight PowerShell scripts were checked the same way, line by line. As a runtime check, every path expression was evaluated with the variables unset and confirmed to reproduce the original location exactly, and a third-party static analyser reported no undefined names.

`ENVIRONMENT.md` records the software versions captured at run time, sets out what a third party needs to reproduce the analysis, and gives the input layout these defaults assume. It keeps the original environment and the reproduction environment separate, because they need not be the same.

**The `sessionInfo.txt` records.** These capture the software environment at the time of the run. The Windows user name appearing in those interpreter paths has been replaced with `C:\Users\<user>` before publication. That substitution is the only edit made to these records, and it changes no analysis output.

## Verifying this release

```bash
python verify_manifest.py
```

`SHA256_MANIFEST.txt` locks every file in this release. Verify with:

```bash
python verify_manifest.py
```

The script re-hashes every listed file and fails — exit status 1 — on any of four conditions: a listed file is missing, a listed file's hash differs, a manifest line is malformed, or a file is present that the manifest does not account for. The last condition is what makes it usable as a release gate: a file added after the manifest was generated stops the release rather than shipping an archive that disagrees with its own inventory. Nothing is written.

Archive files, editor temporaries and Python bytecode caches are excluded from the inventory by construction, in the generator and the verifier alike, and are reported as `[ignored]` rather than failing the run. If those two exclusion lists ever drift apart, a clean tree would start failing or a dirty one would start passing; `release.ps1` runs this same check and refuses to tag if it fails, so a drift surfaces at the next release attempt rather than in the published archive.

## Data sources

No primary data are redistributed here. All are public and remain under the terms of their own repositories:

- TCGA-PAAD — NCI Genomic Data Commons
- GTEx — GTEx Portal
- CPTAC-PDAC — LinkedOmics
- PXD048644 — PRIDE Archive
- NCBI GEO, main cross-scale analysis and comparators — GSE57495, GSE85916, GSE155698, GSE156405, GSE212966, GSE282302
- NCBI GEO, bulk extension cohorts — GSE205154, GSE21501
- NCBI GEO, additional cellular sensitivity datasets (not pooled with the bulk cohorts) — GSE278688 (whole-cell single-cell RNA-seq), GSE278689 (single-nucleus RNA-seq)
- NCBI GEO, external single-cell assessment — GSE217845 (selected); GSE237846, GSE277782, GSE310352 (screened but not used)

Accessions, preprocessing roles and inclusion rules are listed in the manuscript and its Supplementary Table S14. The ICGC Data Portal has retired; the PACA-CA detectability result derives from expression files retained from the ICGC 25K open object store, and the retained table, extraction script and retrieval log are in Additional file 2.

## Where each analysis lives

The code deposit and the paper's supplementary files are complementary, and the split is deliberate: code, frozen rules and run records here; machine-readable results there. This table pairs them.

| Analysis | Code in this deposit | Reported results |
|---|---|---|
| Primary single-cell and bulk cross-scale analysis | `analysis/scripts/` | Additional file 2 |
| Tumour-level burden and sparsity extension | `analysis/v72_burden/` | Additional file 2, §23 and Table S40 |
| State-restricted V04/V05 analyses | `analysis/V04_V05_state_restricted/` | Additional file 2 |
| Bulk extension cohorts | `analysis/phaseB_GSE205154/`, `phaseB_GSE21501/` | Additional file 2, §§20–22 |
| Additional cellular sensitivity datasets (GSE278688 whole-cell, GSE278689 single-nucleus) | `analysis/phaseB_GSE278688/`, `phaseB_GSE278689/` | Additional file 2, §§20–21 |
| Four- and five-cohort meta-analysis | `analysis/extension_meta/scripts/` | Additional file 2, §§20–22 |
| Paired RNA–protein representation | `paired_representation/` | Additional file 3 |
| External single-cell assessment (GSE217845) | `external_single_cell/` | Additional file 4 |

Two exceptions to the split, both small and both included here because the reading guide's warnings are about them: `analysis/extension_meta/meta_results.csv`, whose historical notes are superseded by the errata, and the two errata tables under `analysis/v74h_reporting_corrections/`. A reader checking the reporting corrections needs those files next to the guide that describes them. Everything else in the result category stays in Additional files 2–4.

## Analysis rules and pre-specification

The cross-scale question, its estimands and its scoring were frozen on 22 July 2026 (Asia/Shanghai), before any cross-scale estimate was inspected; the decision rules were written into the protocol document on 10 August 2026. `analysis/PROVENANCE_TIMELINE.md`, the locked protocol, the per-dataset frozen cohort definitions and the manuscript's Supplementary §1 record that chronology.

Where a result proved sensitive to an analysis choice, the code reproduces the sensitivity rather than a single selected value — for example the four-cohort and five-cohort modified Hartung–Knapp intervals, and the paired-representation results under `paired_representation/`.

## Licence

MIT — see `LICENSE`. Third-party data are not covered by this licence.

## Citation

See `CITATION.cff`. If you use this code, please cite the article and this repository.
