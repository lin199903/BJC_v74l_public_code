# Cross-target Reference Panel — Pre-Registered Specification

> **Reading note — this is an internal freeze record, not a formal pre-registration.**
>
> "Pre-registered" here means only that the specification below was written down and frozen
> inside this project before the cross-target association results were inspected. It was
> **not** deposited with, or registered by, any third-party registry, and the panel is
> **not** reported as a preregistered validation study.
>
> The plan and the executable script that was actually run differ in four respects, and the
> two were written nine minutes apart. Those differences are recorded in
> `../ANALYSIS_SPECIFICATION_AMENDMENTS.md`, which governs how this panel is to be read:
> the panel is reported as a **same-system calibration**. Where this file and the amendments
> record disagree, the amendments record is authoritative.
>
> The frozen text below is reproduced unedited, because it is part of the provenance chain.
> Its title is not evidence of external registration and should not be cited as such.

**Frozen date:** 2026-07-14
**Status:** IMMUTABLE. All targets, contexts, expected identities, programs, ecology modules, and test families below were locked BEFORE inspecting any cross-target association results. No target will be added, removed, or reclassified based on observed outcomes.

## 1. Purpose

To assess whether the patient-tissue identity framework that separates *target presence*, *state association*, and *evidence of functional control* recovers distinct expected tissue patterns across pre-specified target-context pairs, rather than mechanically producing a uniform negative or epithelial-state assignment. This provides initial cross-target discriminative face validity, not formal validation of a universal target-classification method.

## 2. Pre-Specified Target-Context Panel

| # | Target | Context | Expected Identity | Reference Class |
|---|--------|---------|-------------------|-----------------|
| 1 | CD274 | TCGA-PAAD (n=178) | Immune / cytolytic ecology, PI3K-associated | Checkpoint-centred (same-cancer positive calibration) |
| 2 | TIGIT | TCGA-NSCLC (LUAD n=161 + LUSC n=501) | T-cell exhaustion / checkpoint co-expression | Checkpoint-centred (cross-cancer calibration) |
| 3 | ERBB2 | TCGA-BRCA (n=847) | HER2/ERBB pathway, copy-number/protein anchoring available | Pathway-coupled |
| 4 | EPCAM | TCGA-PAAD (n=178) | Epithelial differentiation state, no immune-centred signal expected | Epithelial-state |
| 5 | CSF1R | TCGA-PAAD (n=178) | Macrophage / myeloid ecology | Myeloid-state |
| 6 | SIGLEC15 | TCGA-PAAD (n=178) | Prior proposed: PI3K-coupled checkpoint (NOT expected to recover) | Identity-mismatch (core case) |

**Selection rationale:** Pairs were chosen based on established clinical or tissue identities and availability of matched patient-tissue data, BEFORE inspecting the corresponding association results.

## 3. Pre-Specified Programs (Layer 2)

Each target has ONE pre-specified expected program. Program scoring uses mean-z of the program gene set. **Leave-target-out scoring is mandatory**: the target gene itself is removed from its expected program gene set before scoring. No other gene is removed based on observed correlations.

| Target | Expected Program | Program Gene Source |
|--------|------------------|---------------------|
| CD274 | Cytolytic immune (GZMA, GZMB, PRF1, IFNG, NKG7, GNLY, KLRD1, KLRB1) | Literature-defined cytolytic signature |
| TIGIT | T-cell exhaustion / checkpoint (PDCD1, CTLA4, HAVCR2, LAG3, TIGIT, CD28, ICOS) | Canonical immune checkpoint module |
| ERBB2 | HER2/ERBB signaling (ERBB3, ERBB4, NRG1, NDRG1, GRB7, MED1) | KEGG/Reactome ERBB signaling (ERBB2 removed) |
| EPCAM | Epithelial differentiation (KRT18, KRT19, KRT8, EPCAM, CDH1, CLDN4, CLDN7) | Hallmark epithelial |
| CSF1R | Macrophage / myeloid (CD68, CSF1R, CD14, ITGAM, AIF1, LYZ, FCGR3A) | Literature-defined myeloid signature |
| SIGLEC15 | PI3K/AKT/mTOR (Hallmark PI3K, KEGG PI3K-Akt, Reactome PI3K) | Main manuscript PI3K specifications |

## 4. Pre-Specified Competitive Ecology (Layer 3)

All six targets use the SAME five ecology modules, scored uniformly:

1. **T-cell / cytolytic score** (same as CD274 expected program)
2. **Myeloid / macrophage score** (same as CSF1R expected program)
3. **Epithelial differentiation score** (same as EPCAM expected program)
4. **Stromal / fibroblast score** (DCN, LUM, COL1A1, COL3A1, PDGFRB, ACTA2)
5. **MUC1/CD24 mucin-state score** (MUC1, CD24, MUC3A, MUC4 from main manuscript)

## 5. Pre-Specified Test Families and Multiplicity Correction

### Family A — Expected program tests
- 6 tests (one per target)
- Effect size threshold: |ρ| ≥ 0.20
- BH-FDR computed across the 6 tests in this family

### Family B — Competitive ecology tests
- 30 tests (6 targets × 5 ecology modules)
- Effect size threshold: |ρ| ≥ 0.20
- BH-FDR computed across the 30 tests in this family

### Family C — Orthogonal anchoring (descriptive only)
- Per-target: protein detection, copy-number (ERBB2 only), spatial (SIGLEC15 only)
- Reported as descriptive results with 95% CIs; no FDR correction applied
- Not mixed with Families A and B

### Reporting
For every test:
- Spearman ρ
- Patient-level bootstrap 95% CI (B = 5,000)
- Raw p-value
- BH-FDR within the pre-specified family

All six targets' results will be reported in full, regardless of whether they match the expected pattern.

## 6. Statistical Method

- Spearman correlation between target expression and each program/ecology score
- Patient-level non-parametric bootstrap (B = 5,000) for 95% CI
- BH procedure for FDR within each pre-specified family
- No TOST equivalence testing on the reference panel (TOST is reserved for the SIGLEC15 core PI3K null result in the main text)

## 7. Adjustments

For all targets, minimum adjustments:
- Tumor purity (via ESTIMATE or TCGA-provided)
- Cell composition (epithelial, myeloid, stromal fractions via Scaden or similar)

For ERBB2-BRCA, additional pre-specified sensitivity analyses:
- Full TCGA-BRCA cohort
- HER2-enriched (PAM50) subgroup
- ERBB2-amplified subgroup (copy-number > 2)

## 8. Output

- Supplementary Figure S2: Cross-target heatmap
- Supplementary Table S21: Pre-specified panel (this document, tabular form)
- Supplementary Table S22: Full results table

## 9. Interpretation Boundary

This panel provides **initial cross-target discriminative face validity**. It does NOT:
- Prove the framework is universally applicable
- Establish sensitivity/specificity of a target-classification tool
- Replace direct protein or perturbation evidence
- Constitute external validation on independent cohorts

The panel tests ONE question: does the framework recover distinct expected tissue patterns across different target identities, rather than assigning all targets to a common category?
