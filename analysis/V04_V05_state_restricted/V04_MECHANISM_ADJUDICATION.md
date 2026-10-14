# V04 Mechanism Adjudication — SIGLEC15–PDAC v0.4 Discriminating Analysis

**Date of analysis:** 10 September 2026
**Status:** POST HOC mechanistic-discrimination analysis. **Not preregistered.** The original cross-scale question was locked on 22 July 2026; the v0.4 discriminative question was raised only after the manuscript results were known.
**Seed:** 15072026 (single root; independent streams via `SeedSequence([15072026, k])`). **Bootstrap:** 10,000 dataset-stratified patient resamples.
**Inferential unit:** patient/sample throughout. Cells never entered an inferential test.

All evidence partitions below use the four mandated categories. Nothing is promoted across partitions without a stated reason.

---

## 1. Completion summary

| Item | Result |
|---|---|
| **Gate 0 verdict** | **PASS** — 3/3 datasets, 100% deterministic mapping, target genes absent from marker definitions |
| P1-A1 eligible | **16 patients / 3 datasets** |
| P1-A2 eligible (primary) | **13 patients / 3 datasets** |
| P1-A2 eligible (relaxed S2) | 15 patients / 3 datasets |
| P1-A1 localisation | **L = 0.0445, 95% CI 0.0176–0.0859 → SUPPORTIVE** |
| P1-A2 within-state | **K<sub>i,T</sub> = 0.996, 95% CI 0.560–1.558 → SUPPORTIVE** |
| P1-B A–O–K | ρ<sub>AO</sub> = −0.450 [−0.790, 0.021]; ρ<sub>AK</sub> = −0.184 [−0.568, 0.253]; ρ<sub>OK</sub> = 0.221 [−0.329, 0.600] |
| P1-B robustness | **Only ρ<sub>AO</sub> is stable** (negative under all leave-one-dataset-out configurations). ρ<sub>AK</sub> and ρ<sub>OK</sub> both flip sign when GSE156405 is omitted. |
| **Interpretation matrix cell** | **A** — localisation supportive + residual coupling supportive |
| **M2-only sufficient?** | **No.** |
| **M3 status** | **Observationally supported. Not causal.** |
| Single-patient influence | **0 of 16 (P1-A1) and 0 of 13 (P1-A2) leave-one-out cases flipped the verdict** |

---

## 2. Gate 0 — annotation provenance (DIRECTLY OBSERVED)

| Dataset | Included myeloid | Annotated | Matched | Mapping | Unmatched | Duplicates |
|---|---|---|---|---|---|---|
| GSE155698 | 1,715 | 1,715 | 1,715 | 100.00% | 0 | 0 |
| GSE156405 | 4,635 | 4,635 | 4,635 | 100.00% | 0 | 0 |
| GSE212966 | 4,081 | 4,081 | 4,081 | 100.00% | 0 | 0 |
| **Total** | **10,431** | **10,431** | **10,431** | **100.00%** | **0** | **0** |

Supporting facts, each independently verified:

- **Source.** `data/analysis/Task5_myeloid_subtypes/myeloid_combined_obs.csv`, generated 2026-06-06 21:15:30 by `scripts/tasks/_task5_myeloid_subtypes.py`. This **precedes** all SIGLEC15–ACP5/MMP9 cross-scale results (July 2026 onward).
- **Assignment rule.** Each cell takes the highest of nine `scanpy.tl.score_genes` marker-module scores (`idxmax`). Deterministic; no random step, no clustering, no outcome information.
- **Marker union.** 41 module entries collapsing to **36 unique genes** across nine subtypes. **SIGLEC15, ACP5 and MMP9 are absent from every module.** The TREM2_TAM module is `TREM2, CD68, CX3CR1, APOE, SPP1`.
- **Expression layer.** Task5 scored modules with `use_raw=True` and wrote a SIGLEC15 column taken from `.raw`. V04 deliberately re-read SIGLEC15/ACP5/MMP9 from `.X` to match the manuscript. **Binarisation is identical across layers**: 56/1,715, 69/4,635 and 72/4,081 positives, i.e. 197/10,431 = 1.889% — exactly the count reported in the manuscript. The layer difference is therefore immaterial here, and is documented rather than hidden.
- **Not outcome-driven.** The state definition used neither ACP5, MMP9, nor any coupling estimate. Recorded as an explicit field in `V04_GATE_EVALUATION.json`.

**Verdict: PASS.** Fine-state analyses were permitted to proceed.

---

## 3. DIRECTLY OBSERVED

1. **SIGLEC15-positive tumour-myeloid cells are enriched in the pre-existing TREM2_TAM label at patient level.** Mean patient-level contrast L = p<sub>i,T</sub> − p<sub>i,O</sub> = **0.0445** (95% CI 0.0176–0.0859), 16 patients, 13/16 with positive direction, median 0.0226, dataset medians 0.0332 / 0.0160 / 0.0306. Median enrichment ratio 3.32.
2. **The localisation contrast is small in absolute terms.** Four percentage points of detection prevalence. Its CI excludes zero but does not describe a large separation.
3. **A within-TREM2_TAM SIGLEC15 contrast persists.** K<sub>i,T</sub> = **0.996** (95% CI 0.560–1.558), 13 patients, **13/13 positive**. With depth adjustment 0.634 [0.272, 1.036] (11/13); with within-state depth matching 0.738 [0.252, 1.355] (10/13); relaxed estimability 0.904 [0.523, 1.401] (15/15).
4. **The positive within-state contrast was not confined to the samples with the fewest SIGLEC15-positive cells, and it remained after exclusion of the largest patient-specific estimate.** Stratifying the primary eligible set at ≤5 versus >5 SIGLEC15-positive TREM2_TAM cells gives 0.733 [0.242, 1.224] (n=6) and 1.222 [0.646, 2.177] (n=7); excluding GSE212966/PDAC5 gives 0.754 [0.476, 1.027]. Cell sparsity remains a component of measurement uncertainty and is not claimed to be absent — what is shown is that the result does not depend on the sparsest samples or on the most extreme one.
5. **Neither verdict depends on the dominating patient.** Excluding GSE212966/PDAC5 (the extreme sample already identified in S24) gives L = 0.0247 [0.0143, 0.0352] and K<sub>i,T</sub> = 0.754 [0.476, 1.027]. All 16 and all 13 leave-one-patient-out cases retain their verdict.
6. **The three patient-level quantities are not one rank axis, but only one of the three separations is robust.**

   Full sample: ρ<sub>AO</sub> = −0.450 [−0.790, **0.021**]; ρ<sub>AK</sub> = −0.184 [−0.568, **0.253**]; ρ<sub>OK</sub> = 0.221 [−0.329, 0.600].

   Leave-one-dataset-out:

   | Omitted | n | ρ<sub>AO</sub> | ρ<sub>AK</sub> | ρ<sub>OK</sub> |
   |---|---|---|---|---|
   | — (full) | 20 / 19 | **−0.450** | **−0.184** | **+0.221** |
   | GSE155698 | 15 | −0.596 | −0.169 | +0.191 |
   | GSE156405 | 11 | −0.373 | **+0.173** | **−0.155** |
   | GSE212966 | 14 | −0.363 | −0.363 | **+0.505** |

   **ρ<sub>AO</sub> is negative under every configuration and never flips sign.** ρ<sub>AK</sub> and ρ<sub>OK</sub> both flip sign when GSE156405 is omitted, and ρ<sub>OK</sub> reaches 0.505 when GSE212966 is omitted. No leave-one-**patient**-out case flips any sign, but single-dataset influence on the two weak relationships is decisive. Any statement about ρ<sub>AK</sub> or ρ<sub>OK</sub> is therefore unreliable.

---

## 4. SUPPORTED INTERPRETATION

1. **Broad TREM2_TAM localisation is real but modest.** The interval lies wholly above zero on 16 patients across three datasets, survives removal of the extreme patient, and is consistent in direction within each dataset. It supports M2 as a genuine phenomenon rather than a thresholding accident.
2. **M2 alone is insufficient.** A residual within-state contrast is present in every one of the 13 primary-eligible patients, survives depth adjustment, within-state depth matching, relaxed estimability and stratification by positive-cell density. Whatever localisation contributes, it does not account for the cellular contrast on its own. Under the pre-specified matrix this is cell **A**, which is the outcome the contract anticipated as "M2-only is insufficient; M3-like architecture becomes observationally supported".
3. **The rarity explanation and the localisation explanation coexist rather than compete.** M1 remains necessary: at 1.9% prevalence an aggregate cannot track a conditional relation. M2 and residual M3 describe *where* the relation sits, M1 describes why aggregation destroys it.
4. **Target abundance is a poor proxy for state occupancy in these data, and this is the one A–O–K finding that survives every perturbation tested.** ρ<sub>AO</sub> is negative in the full sample and under all three leave-one-dataset-out configurations (−0.596, −0.373, −0.363), and its interval upper bound is 0.021. A shared-axis reading would require a strong positive alignment; that reading is effectively excluded. Abundance versus coupling points the same way (upper bound 0.253) but does not survive leave-one-dataset-out.
5. **The claim that abundance, occupancy and coupling are three interchangeable patient-level quantities is not supported.** One pairwise alignment (abundance–occupancy) is robustly excluded from being strong. The other two are either uninformative or sign-unstable across datasets, so they neither support nor refute interchangeability.

---

## 5. BIOLOGICALLY PLAUSIBLE BUT UNTESTED

1. **Context-gated coupling.** K<sub>i,T</sub> varies markedly between patients (dataset medians 0.31, 1.19, 0.65 for the primary set), which is compatible with coupling being modulated by niche, ligand or adaptor environment. Nothing in these data identifies such a modulator.
2. **Why the within-state contrast exceeds the whole-myeloid contrast** (0.996 versus 0.639) is unresolved. A state-restricted effect, a shift in the ACP5/MMP9 distribution inside TREM2_TAM cells, and residual depth or sparsity structure are all compatible. This is an observation awaiting an explanation, not evidence of a stronger mechanism.
3. **A ligand- or glycan-defined subset.** The reconstruction's historical reading treats SIGLEC15 as a receptor embedded in a myeloid state rather than a tumour-wide scalar. Whether the detectability contrast indexes a specific ligand context cannot be assessed here.
4. **Any causal SIGLEC15 → DAP12/SYK → ACP5/MMP9 chain in human PDAC.** Consistent with the source literature, but not tested, and not testable with these measurements.

---

## 6. NOT ANSWERABLE WITH CURRENT DATA

1. Protein-level co-expression, compartment-resolved SIGLEC15, or any statement about protein abundance. Public discovery proteomics provides no quantifiable SIGLEC15 row.
2. DAP12 recruitment or phosphorylation; SYK/PI3K activation attributable to SIGLEC15.
3. Whether the TREM2_TAM labels correspond to an externally validated lineage. SPP1 sits inside the TREM2 module and is biologically adjacent to osteoclast-like programmes, so the restriction analysis tests robustness to a defined annotation construct, not a lineage.
4. Whether patient-level K<sub>i</sub> heterogeneity reflects niche biology or estimation noise. With 2–24 SIGLEC15-positive TREM2_TAM cells per patient these two cannot be separated.
5. Treatment-response relevance, clinical qualification, or biomarker validity.
6. Whether a paired same-tumour cross-scale design would reproduce the localisation gradient. No such material was analysed.

---

## 7. Model adjudication (M1–M4 updated)

| Model | Prior status | **Updated status** | Basis |
|---|---|---|---|
| **M1** rarity / dilution | Supported, insufficient alone | **Unchanged.** Necessary, not sufficient. | Not re-tested in v0.4; the S22 dilution result stands. |
| **M2** state localisation / composition | Supported, sufficiency unresolved | **Supported but demonstrably insufficient.** | L = 0.0445 [0.0176, 0.0859] > 0 in 13/16 patients, robust to extremes; but K<sub>i,T</sub> remains positive in 13/13. |
| **M3** residual / context-dependent coupling | Plausible and data-compatible, not demonstrated | **Observationally supported. Explicitly not causal.** | Positive in every eligible patient; survives depth adjustment, matching, relaxed estimability and density stratification. |
| **M4** multicompartment mixing / measurement substitution | Contributing; off-myeloid mixing insufficient | **Unchanged.** | v0.4 adds nothing new; the myeloid-only aggregate still fails. |

**Layered reading, now with patient-level evidence at two of its levels:**

```
state localisation            (M2, supported: L > 0, 13/16 patients)
   ↓
residual within-state coupling (M3-like, observationally supported: K_i,T > 0, 13/13)
   ↓
rarity-driven attenuation      (M1, unchanged)
   ↓
loss of conditioning during aggregation (M4, unchanged)
   ↓
poor patient-level marginal representation
```

This is more specific than "SIGLEC15 is too rare for bulk RNA to see it", and it still does not license "SIGLEC15 activates DAP12/PI3K, therefore ACP5/MMP9 is high".

---

## 8. Manuscript implications

| Manuscript element | Direction | Reason |
|---|---|---|
| Primary cross-scale conclusion | **Unchanged** | v0.4 does not touch the aggregate or bulk layers. |
| Characterisation of the TREM2_TAM enrichment as a compositional nuisance | **Tighten** | It is now partly a real, replicable localisation gradient and partly a source of confounding. Presenting it purely as confounding to be adjusted away understates the data. |
| Scope of the "rarity" explanation | **Increase, modestly and explicitly** | Rarity remains necessary; it is no longer the only describable level. Two further levels now carry patient-level evidence. |
| Any statement about mechanism | **Must not increase** | Everything above is post hoc, RNA-only, and based on 2–24 positive cells per patient. |
| Claims about protein, causality, treatment response, clinical use | **Must not change** | Not answerable with these data. |
| Post hoc labelling | **Must not change** | v0.4 is post hoc and must be labelled as such wherever it is cited. |

If v0.4 is written into the manuscript at all, it belongs as an explicitly post hoc mechanistic-discrimination paragraph with the localisation and within-state estimates, their intervals, the sparsity strata and the extreme-patient sensitivity reported together. An unqualified "we found residual within-state coupling" would overstate 2–24 cells per patient.

---

## 9. Deviations from the protocol

All deviations are recorded with their reason. None was made after seeing an outcome in order to change a verdict.

1. **Added leave-one-patient-out diagnostics** (the contract specified leave-one-dataset-out only). Reason: the contract's STOP RULES require judging whether extreme observations materially determine the result, and with n = 13–16 a single-patient leave-out is the only way to test that. No verdict changed; the diagnostics are reported in full.
2. **Added a positive-cell-density stratification** (≤5 vs >5 SIGLEC15-positive TREM2_TAM cells). Reason: required by the same STOP RULE on sparse positive-cell support. Both strata were pre-specified as descriptive before the values were read, and neither replaces the primary analysis.
3. **Random-stream derivation.** Independent streams were produced as `SeedSequence([15072026, k])` instead of reusing one generator, to prevent accidental correlation between analyses. The root seed is the contract seed and the derivation is recorded in `_v04_provenance.json`.
4. **Implementation cross-validation.** The relaxed (S2) estimate reproduced the earlier S24 value (0.904 [0.523, 1.401] versus 0.904 [0.525, 1.392]) from an independently written script and seed, confirming that the extraction path matches the manuscript pipeline.
5. **No deviations in scope.** No new pathway, survival, machine-learning, deconvolution, threshold, signature or clustering analysis was run. The ambient-RNA branch was not reopened. No eligibility threshold was changed after results were seen, and no alternative TAM definition was generated.
6. **Manuscript untouched.** No manuscript file was modified during this analysis.

---

## 10. Outputs

| File | Content |
|---|---|
| `V04_ANALYSIS_SPECIFICATION.md` | Frozen design, estimands, eligibility, verdict rules, stop rules |
| `V04_STATE_ANNOTATION_PROVENANCE.csv` | Per-dataset annotation provenance and mapping audit |
| `V04_GATE_EVALUATION.json` | Gate 0 verdict, marker union, layer note, per-dataset mapping |
| `V04_P1A_LOCALISATION.csv` | Per-patient p<sub>i,T</sub>, p<sub>i,O</sub>, L<sub>i</sub>, enrichment ratio |
| `V04_P1A_WITHIN_STATE.csv` | Per-patient K<sub>i,T</sub>, depth-adjusted, depth-matched, eligibility flags |
| `V04_P1B_AOK.csv` | Per-patient A<sub>i</sub>, O<sub>i</sub>, K<sub>i</sub> |
| `V04_MECHANISM_ADJUDICATION.md` | This report |
| `V04_FIG_MECHANISM_DISCRIMINATION.pdf` | Gate 0, localisation, within-state and A–O–K panels |

Supporting internal files (not part of the eight deliverables): `_v04_sample_metrics.csv`, `_v04_subtype_counts.csv`, `_v04_influence_diagnostics.json`, `_v04_layer_crosscheck.json`, `_v04_provenance.json`, `_v04_p1a_*_summary.json`, `_v04_p1b_aok_summary.json`, and the three analysis scripts.

**Success criterion met.** The frozen analysis is interpretable under the cell that actually occurred (matrix cell A), and the same script would have returned a readable verdict for any other cell, including a null or uninformative one.
