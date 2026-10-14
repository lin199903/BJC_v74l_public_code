# Cross-target calibration panel: implementation amendments

> **Note on file locations.** This record cites paths from the project's working tree (`data/analysis/...`, `npjPrecOncol/...`). Two of them are included in this code deposit, at the locations below; the rest are result files and are published with the paper as Additional file 2.
>
> | Cited as | In this deposit |
> |---|---|
> | `data/analysis/cross_target_panel/PRE_SPECIFICATION.md` | `analysis/scripts/cross_target_panel/PRE_SPECIFICATION.md` |
> | `data/analysis/cross_target_panel/run_cross_target.py` | `analysis/scripts/cross_target_panel/run_cross_target.py` |
> | `data/analysis/cross_target_panel_verified/`, `npjPrecOncol/*.csv` | Additional file 2 |

## Purpose

This record documents differences between the initial cross-target analysis plan and the executable analysis that generated the reported results. It is provided to prevent the calibration panel from being misrepresented as an immutable preregistration.

## File chronology

- `data/analysis/cross_target_panel/PRE_SPECIFICATION.md`: modified 2026-07-14 01:39:35.
- `data/analysis/cross_target_panel/run_cross_target.py`: modified 2026-07-14 01:48:50.
- Original script output `cross_target_results_flat.csv`: modified 2026-07-14 01:50:18.

The target--context pairs and broad program assignments were recorded before the output file was generated. The short interval between the plan and executable script, together with the changes below, means that this panel is reported as a same-system calibration rather than as a preregistered validation study.

## Amendments

1. **MUC-state module.** The plan listed four genes. The executable script used the 12-gene MUC1/CD24 module used in the main analysis. The four-gene listing was treated as a clerical mismatch, not as the implemented module.
2. **SIGLEC15 expected program.** The plan referred to Hallmark, KEGG and Reactome PI3K definitions. To retain one expected-program test per target in the calibration panel, the executable script used the leave-target-out Hallmark PI3K module. The broader definitions remained in the main 16-specification PI3K analysis.
3. **Covariate adjustment.** The plan proposed tumor-purity and deconvolution adjustment. The executable calibration panel used unadjusted target--program Spearman correlations. Scaden was not used after leave-one-donor-out pseudobulk validation produced near-constant estimates and failed to recover known epithelial and myeloid proportions. The marker-adjusted analyses in the main study are separate from this calibration panel.
4. **Multiplicity family.** The executable script applied Benjamini--Hochberg correction to all 30 tissue-association tests. EPCAM--Epithelial and CSF1R--Myeloid are structurally identical to their expected-program tests. They are marked with daggers to prevent double interpretation but remain in the 30-test Family B correction.

## Reproducibility check

The original executable script was rerun in `data/analysis/cross_target_panel_verified/`. The rerun reproduced all 36 Spearman correlations and bootstrap confidence intervals exactly; p-values and FDR values differed only at floating-point precision (maximum absolute difference approximately `1e-16`). The canonical flat output for this manuscript version is:

`npjPrecOncol/cross_target_results_verified.csv`

The later manually assembled `npjPrecOncol/cross_target_results_flat.csv` is not used in version 5 because its expected-program values and multiplicity labels were not traceable to the executable script.

## Interpretation boundary

The panel asks whether the same analysis can recover distinct RNA-associated programs for reference targets. It does not establish equal sensitivity for low-abundance targets, validate a universal target-classification method, or replace matched protein and perturbation evidence.
