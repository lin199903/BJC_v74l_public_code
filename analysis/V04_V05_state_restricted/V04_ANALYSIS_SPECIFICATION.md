# V04 Analysis Specification — SIGLEC15–PDAC Birdseye Reconstruction, v0.4 Discriminating Analysis

**Status: POST HOC MECHANISTIC-DISCRIMINATION ANALYSIS. NOT PREREGISTERED.**
This specification must never be described as prespecified relative to the original manuscript (v32/v34) results.
The original cross-scale question was locked on 22 July 2026; the v0.4 discriminative question was raised after those results were known.

---

## 1. Background

The submitted manuscript establishes, on public data:

1. SIGLEC15 RNA is detected in 99.4% of bulk TCGA-PAAD tumours but in only 1.9% of tumour-myeloid cells.
2. Within patients, SIGLEC15 RNA-detectable myeloid cells carry higher ACP5/MMP9 programme scores than non-detectable cells.
3. That contrast survives depth adjustment, within-sample depth matching, and post hoc fine-subtype adjustment.
4. Once cellular conditioning is removed, the same relation is not stably reproduced by myeloid sample aggregates or whole-tumour bulk RNA.
5. SIGLEC15-positive cells are enriched in the pre-existing TREM2_TAM annotation.
6. Diluting a TREM2-centred state to ~2% prevalence attenuates the sample-level correlation from 0.67 to 0.16.
7. Patient-specific coupling estimates vary substantially and show no clear monotonic relation to SIGLEC15 prevalence or positive-cell intensity.

Birdseye reconstruction produced four competing explanations:

| Model | Statement |
|---|---|
| **M1** rarity / dilution | The cellular coupling is real but too rare to survive aggregation. |
| **M2** state localisation / composition | SIGLEC15-positive cells preferentially occupy a TREM2-like remodelling TAM state, and that composition alone generates the overall contrast. |
| **M3** residual / context-dependent coupling | Even within the same broad TAM state, SIGLEC15-positive cells retain higher ACP5/MMP9, with patient-varying magnitude. |
| **M4** multicompartment mixing / measurement substitution | Whole-tumour marginal RNA is a different estimand from the conditional within-myeloid quantity. |

The purpose of v0.4 is **not** to discover another mechanism. It is to discriminate M2 from residual M3 and to test whether target abundance, state occupancy and target–state coupling are interchangeable patient-level quantities.

---

## 2. Data and fixed analysis contract

| Item | Value |
|---|---|
| Datasets | GSE155698, GSE156405, GSE212966 (only) |
| Inclusion | `tissue ∈ {PDAC, PDAC_tumor, PDAC_tissue}` and `celltype == "Myeloid"`; sample retained if ≥50 included myeloid cells |
| Exposure | SIGLEC15 stored normalised expression > 0, read from `.X` |
| Outcome | ACP5/MMP9 programme score = mean of within-sample z-scores of ACP5 and MMP9, computed over all included myeloid cells of that sample |
| Inferential unit | patient / sample. Never the cell. |
| Seed | 15072026 (single root; independent streams derived as `SeedSequence([15072026, k])`) |
| Bootstrap | 10,000 dataset-stratified patient resamples; patients equally weighted |
| Programme | `Derosiers_ACP5_MMP9_core` only |

**Explicitly out of scope for v0.4:** new pathway analyses, survival analyses, machine learning, deconvolution, new TAM signatures (SPP1/GPNMB/other), new clustering solutions, new detectability thresholds, ambient-RNA branch.

---

## 3. GATE 0 — fine-myeloid annotation provenance audit

**Question.** Can the pre-existing TREM2_TAM fine-myeloid annotation be mapped deterministically onto the exact cells used by the manuscript expression objects?

**Acceptable sources.** Published cell-level metadata; previously generated deterministic cell-ID annotation tables; fixed cell-ID mappings created before v0.4.

**Forbidden.** Generating a new clustering solution; selecting a favourable cluster; tuning the state definition after seeing outcomes.

**Reported per dataset.** Annotation source and date; number of included myeloid cells; number uniquely mapped; mapping percentage; unmatched counts; duplicate/ambiguous cell IDs; SIGLEC15-positive and -negative mapping rates; whether ACP5, MMP9 or current coupling results were used to define the state.

**PASS criterion.** ≥95% of included myeloid cells uniquely and deterministically mapped, in ≥2 datasets, with no forbidden gene in the marker definitions.

**Stop rule.** If one dataset fails, it is excluded from fine-state analyses. If fewer than two datasets remain, P1-A and the occupancy part of P1-B are not inferred. If the state definition used the current SIGLEC15/ACP5/MMP9 result, STOP.

---

## 4. P1-A1 — patient-level state localisation

**Estimand.** For each eligible patient *i*:

```
p_i,T = (# SIGLEC15+ and TREM2_TAM) / (# TREM2_TAM)
p_i,O = (# SIGLEC15+ and non-TREM2 myeloid) / (# non-TREM2 myeloid)
L_i   = p_i,T − p_i,O
```

Population estimand: the equal-patient mean of `L_i`.

**Eligibility.** ≥10 TREM2_TAM cells **and** ≥10 other myeloid cells. Zero-positive compartments are valid data and must not trigger exclusion.

**Inference.** Requires ≥8 eligible patients and ≥2 datasets. Primary interval: dataset-stratified patient bootstrap.

**Verdict rule.** `LOCALISATION SUPPORTIVE` if the 95% CI for mean `L` lies entirely above 0. Otherwise `LOCALISATION UNINFORMATIVE`. **No formal "absence of localisation" verdict is permitted**, because a zero or negative estimate with a wide interval cannot distinguish absent localisation from insufficient power.

**Secondary descriptor (not a verdict input).** TREM2_TAM enrichment ratio.

---

## 5. P1-A2 — within-TREM2_TAM coupling

**Estimand.** For each eligible patient *i*:

```
K_i_T = mean ACP5/MMP9 among SIGLEC15+ TREM2_TAM cells
      − mean ACP5/MMP9 among SIGLEC15− TREM2_TAM cells
```

**Primary eligibility.** ≥10 TREM2_TAM cells, ≥2 SIGLEC15+ TREM2_TAM cells, ≥5 SIGLEC15− TREM2_TAM cells.

**Inference.** Requires ≥8 estimable patients and ≥2 datasets; otherwise descriptive only and `MECHANISTIC DISCRIMINATION = UNINFORMATIVE`.

**Pre-specified sensitivities.**

- **S1** depth-adjusted programme score, covariates `log1p(total counts)`, `log1p(detected genes)`, mitochondrial fraction (same residualisation as the manuscript).
- **S2** relaxed estimability (≥1 SIGLEC15+ cell). Descriptive/estimability sensitivity only; must never replace the primary result.
- **S3** within-TREM2_TAM 1:1 nearest-neighbour depth matching without replacement on `log1p(total counts)`, run only if the primary analysis has ≥8 estimable patients.

**Verdict rule.**

| Condition | Verdict |
|---|---|
| 95% CI for mean `K_i_T` entirely > 0 | residual coupling SUPPORTIVE |
| upper 95% CI < 0.20 | residual coupling NON-SUPPORTIVE |
| otherwise | residual coupling UNINFORMATIVE |

The 0.20 line is an operational descriptive boundary inherited from the manuscript. It is not a biological or clinical cutoff.

**Interpretation matrix.**

| Localisation | Residual coupling | Reading |
|---|---|---|
| supportive | supportive | Broad localisation exists but cannot explain the entire cellular relation; M2-only is insufficient and an M3-like architecture becomes observationally supported. |
| supportive | non-supportive | Localisation-dominant explanation; M2 receives primary support. |
| supportive | uninformative | Localisation established; M2 vs residual M3 not discriminable. |
| uninformative | supportive | TREM2_TAM localisation is not necessary to explain the association; a within-state association remains. |
| uninformative | non-supportive | The current TREM2_TAM construct does not explain the original cellular contrast. Do **not** search for a new TAM state in this cycle. |
| uninformative | uninformative | No mechanistic upgrade. |

---

## 6. P1-B — abundance–occupancy–coupling decomposition

Three patient-level quantities, defined once:

| Symbol | Definition |
|---|---|
| `A_i` | SIGLEC15 RNA-detectable fraction of included tumour-myeloid cells (target abundance) |
| `O_i` | TREM2_TAM fraction of included tumour-myeloid cells (state occupancy) |
| `K_i` | the existing locked manuscript patient-specific detectable-minus-non-detectable ACP5/MMP9 contrast (target–state coupling) — **not redefined** |

**Primary quantities.** Exactly three rank relationships: `rho_AO`, `rho_AK`, `rho_OK`, each with a dataset-stratified patient-bootstrap 95% CI, a scatterplot carrying the dataset identifier, and a leave-one-dataset-out sensitivity.

**Forbidden in P1-B.** Dichotomising patients; mediation analysis; SEM; PCA; patient clustering; predictive models; latent variables.

**Descriptive scales.** |rho| < 0.20 little observed rank alignment; 0.20–0.50 intermediate; ≥0.50 large. Descriptive only; the p value is not the decision quantity.

**Interpretation branches.**

1. Strong A–O with weak A–K → abundance behaves more like occupancy than coupling.
2. Strong O–K with weak A–O → state burden relates to coupling while target abundance is a poor state proxy.
3. Strong A–K with other relationships weak → abundance may contain coupling information; abundance/coupling separation claims must be reduced.
4. All three strongly positive → a single-axis approximation remains plausible; cross-scale loss should be attributed more heavily to rarity, mixing and measurement scale.
5. All three weak/mixed and CIs exclude large concordance → empirical support for a multidimensional patient architecture.
6. All three weak/mixed but CIs wide → **UNINFORMATIVE**; multidimensionality must not be claimed.

**Chronology constraint.** A–K was already partly inspected in the earlier patient-specific audit. It must not be presented as a new prospective discovery. The new contribution is the complete A–O–K structure, especially `O_i`.

---

## 7. Technical null

The sequencing-depth / ambient-RNA branch is not reopened. The manuscript already reports depth adjustment, depth matching, stricter detectability thresholds, deep/shallow strata and patient-level inference. Raw counts are unavailable, so formal SoupX/CellBender decontamination cannot be performed. This is recorded as a limitation, not addressed with further arbitrary threshold analyses.

---

## 8. Stop rules

Inference is stopped or downgraded if: fine-state mapping is not reliably recoverable; mapping <95%; fewer than two datasets remain; the state definition was outcome-driven; the primary within-state analysis has <8 estimable patients; only one dataset contributes; SIGLEC15-positive cells are too sparse for stable patient-level estimates; TREM2_TAM occupancy shows insufficient variation; or extreme observations materially determine the result.

A failed analysis may **not** be rescued by lowering eligibility thresholds after seeing results, generating a second TAM definition, switching to SPP1/GPNMB, changing the ACP5/MMP9 outcome, creating a composite signature, or adding pathway analysis.

A negative, null or uninformative result is an acceptable endpoint.

---

## 9. Final mechanism adjudication

After all analyses, only M1–M4 are updated. The adjudication report must separate **directly observed**, **supported interpretation**, **biologically plausible but untested**, and **not answerable with current data**.

It must never claim from public RNA data: specific glycan binding in patient tissue; DAP12 recruitment or phosphorylation; SYK/PI3K activation caused by SIGLEC15; causal induction of ACP5/MMP9 in human PDAC TAM; treatment-response prediction; or protein-level biomarker validity.

---

## 10. Required outputs

1. `V04_ANALYSIS_SPECIFICATION.md` (this file)
2. `V04_STATE_ANNOTATION_PROVENANCE.csv`
3. `V04_GATE_EVALUATION.json`
4. `V04_P1A_LOCALISATION.csv`
5. `V04_P1A_WITHIN_STATE.csv`
6. `V04_P1B_AOK.csv`
7. `V04_MECHANISM_ADJUDICATION.md`
8. `V04_FIG_MECHANISM_DISCRIMINATION.pdf`

The manuscript is not modified during analysis, and the Discussion is not rewritten until `V04_MECHANISM_ADJUDICATION.md` is complete.

**Success criterion.** Not a significant result, but that the same frozen analysis remains scientifically interpretable whether M2 wins, residual M3 is supported, both coexist, or the data remain uninformative.
