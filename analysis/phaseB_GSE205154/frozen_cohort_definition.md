# frozen_cohort_definition.md — GSE205154

**Status:** FROZEN. The later ρ step must obey this document unconditionally.
**Step:** Phase B step 3a — manifest construction and cohort-definition freeze only.
**Date frozen:** 2026-09-28 (Asia/Shanghai)
**Frozen artefact:** `sample_manifest.csv` (289 rows × 13 columns)
**SHA-256 of `sample_manifest.csv`:**

```
f76b81eceb0557cf4c1fee0138e9137a84508369f47ad163255a096af78c708b
```

> A later run proves it ran after this freeze by re-hashing `sample_manifest.csv` and
> checking it against the value above. If the hash differs, the manifest was edited after
> the freeze and the estimate that used it is not covered by this definition.

---

## 0. Governing rule, and what was NOT done in this step

> **Lock the cohort BEFORE looking at the effect.**

**No correlation between SIGLEC15 and any ACP5/MMP9 programme score was computed,
inspected or reported in this step.** No per-sample expression value of SIGLEC15, ACP5 or
MMP9 was loaded. Matrix access was restricted to (i) header lines — column labels — and
(ii) the leading feature-annotation fields of each row, which are feature metadata
identical for every sample. The guards that enforce this are written into
`scripts/01_fetch_metadata.py` and `scripts/03_scan_features.py` and are described in
`run_log.txt`.

**Frozen estimand (verbatim, for context only — not computed here):**

> within each eligible PDAC cohort, estimate the Spearman correlation between
> tumour-level SIGLEC15 expression and an ACP5/MMP9 programme score constructed using the
> same dataset-level standardisation rule as the existing TCGA-PAAD, GSE57495 and GSE85916
> analyses. Cohort-specific estimates will be combined only after cohort eligibility is
> established independently of the observed correlation.

**Exclusion rule carried over unchanged:** a cohort is **never** excluded because of the
observed ρ magnitude or direction.

---

## 1. Cohort identity and deposit

| item | value |
|---|---|
| accession | GSE205154 |
| title | *Bulk RNA-Seq gene expression data for 289 primary and metastatic pancreatic ductal adenocarcinoma (PDAC) tumor samples from FFPE blocks* |
| SuperSeries | GSE281129 (children: GSE205154, GSE281005) |
| BioProject / SRA | PRJNA844022 / SRP377758 |
| publication | PMID 39789181 — *Ongoing replication stress tolerance and clonal T cell responses distinguish liver and lung recurrence and outcomes in pancreatic cancer*, **Nature Cancer** 2025, doi 10.1038/s43018-024-00881-3 |
| institution | Oregon Health & Science University (pathology review OHSU; extraction and sequencing Tempus) |
| platform | GPL20301, Illumina HiSeq 4000, `library_strategy = RNA-Seq`, `library_selection = cDNA`, `molecule = total RNA` |
| library construction | FFPE sections; pathologist-marked regions of **>20 % tumour nuclei** scraped; 100 ng RNA fragmented to ~200 bp; **xGEN Exome Research Panel v1.0 hybrid capture** (this is exome-capture RNA-seq, not whole-transcriptome) |
| quantification | kallisto 0.44.0 pseudoalignment → tximport 1.12.1; GRCh38.p5 / gencode v24 |
| deposited samples | 289 GSM, 289 BioSample, 289 distinct sample titles (all 1:1) |
| public since | 2024-09-30 (submitted 2022-05-31; series last updated 2025-02-19) |

**Raw sequence data were not used.** The deposit states that raw reads are in dbGaP under
`phs003597.v1.p1` (controlled access). Only the processed gene-level matrices were used.

### 1.1 What representation the deposit supplies, and how it was handled

The deposit supplies **four** things:

1. `GSE205154_RAW.tar` (1.77 GB) — 289 per-sample kallisto `_RNA_abundance.h5` files.
2. `GSE205154_Gene_Level_Counts_Estimates.txt.gz` — gene-level estimated counts
   (deposited 2022-05-31), sha256 `4c7004cfe15a94d868e24cf651f0e5d05779709bf233a74b48b084efdf238fdc`
3. `GSE205154_Gene_Level_TPM_Estimates.txt.gz` — gene-level TPM (2022-05-31),
   sha256 `686582f0c00467de844d9ee93356e192da2c03c047f804f2b9d5270c4b62f152`
4. `GSE205154_Gene_Level_TMM_Estimates.txt.gz` — gene-level TMM (edgeR method), **not log
   scaled**, plus a per-gene `TMM_pass_filter` column (added/updated 2024-09-27),
   sha256 `9f362bb06995e28d3f0f2297f89ddf0f1190235c76a1a8cfb5c199d8211895c3`

All three matrices have **60,554 feature rows × 289 sample columns**, with the sample
columns in **identical order** across the three, and identical gene row indices — verified
in `_probe/GSE205154_feature_scan.txt`.

**The CP10K + log1p single-cell rule was not applied.** It governs raw single-cell counts,
and this is a bulk deposit. Nothing was re-quantified from raw reads.

**There is no ambiguity about what the deposit supplies.** What the frozen text does *not*
settle is **which** of the three supplied matrices should feed the programme score. That
is an author decision and is flagged as **AD-2** in §8.

---

## 2. How primary vs metastatic was determined

**This is determinable, and was determined. There is no fallback in force.**

**Primary metadata field:**

```
!Sample_characteristics_ch1 = "tumor type: <value>"
```

**Value mapping (exhaustive over this deposit):**

| value | maps to | n |
|---|---|---|
| `Primary` | primary PDAC tumour | **218** |
| `Met` | PDAC metastasis | **71** |

**Two independent corroborating fields were checked against it, and all 289 samples agree
on all three:**

| field | primary samples | metastatic samples |
|---|---|---|
| `!Sample_description` | `Primary pancreas tumor` (218/218) | `Metastisis to the <site>` (71/71) |
| `!Sample_title` suffix | `T` (207), `T2` (2), `A-T` (2), `F` (7) | `M` (70), `M2` (1) |
| `!Sample_characteristics_ch1 = tissue:` | `Pancreas` (218/218) | 16 distinct sites (71/71) |

`Sample_description` vs `tumor type` conflicts: **0 of 289.**
`tissue` vs `source_name` mismatches: **0 of 289.**

**Tumour vs non-tumour:** all 289 samples are tumour. `cell type: scraped tumor` is set for
289/289, and the extraction protocol states that a clinical pathologist marked regions of
high tumour content on H&E slides before scraping. **The series contains no
adjacent-normal or other non-tumour sample**, so the "tumour tissue" condition excludes
nothing.

**Two edge cases, resolved explicitly (neither is a defect):**

* **`GSM6206096` / `ST-00020064-M`** — a metastasis *located in the pancreas*
  (`tissue: Pancreas`, `description: Metastisis to the Pancreas`). A pancreatic metastasis
  is **not** a primary tumour. It is classified **metastatic**, because the frozen rule
  keys on `tumor type` and not on `tissue`. Had the rule keyed on `tissue`, this sample
  would have been misclassified as primary.
* **`GSM6206144` / `ST-00020092-M`** — `tissue: NA`, `description: Metastisis to the NA`.
  The metastatic **site** is unstated; the **compartment** is unambiguous
  (`tumor type: Met`). Retained as metastatic; the unstated site is recorded in its
  `exclusion_reason`.

**Disease condition:** all 289 samples are PDAC — stated in the series title, summary and
overall design, and consistent with the specimen annotations and the linked publication.
All 218 primary specimens are pancreatic. **Condition met for all 289 samples.**

---

## 3. Rule for patient identity

### 3.1 The prior-scan belief was wrong, and this is the correction

The prior scan recorded **"no patient identifier exposed"**. Verified in this step:

* the GSE205154 series record exposes **no patient field** — its sample-level
  characteristics are only `tumor type`, `tissue` and `cell type`;
* **all 289 BioSample records** (`SAMN2877…`, retrieved by efetch) expose only four
  attributes — `source_name`, `tumor type`, `tissue`, `cell type` — and **no**
  patient/subject/donor attribute; all 867 of their `<Id>` values are distinct.

**But the deposit does expose a patient identifier, through its sibling SubSeries.**
GSE281129 is a SuperSeries with two children, GSE205154 and **GSE281005** (Adaptive
Biotechnologies TCRβ sequencing of the same OHSU patient series). GSE281005 publishes a
characteristic literally named:

```
!Sample_characteristics_ch1 = "public patient_id: ST-########"
```

for **all 505 of its samples**, and it shares the `ST-########` specimen-title namespace
with GSE205154.

### 3.2 The frozen patient rule

> **FROZEN — patient identity rule.**
> The patient key for GSE205154 is the **`ST-########` specimen core** of the sample
> title — the digit string that precedes the specimen suffix (`-T`, `-T2`, `-A-T`, `-M`,
> `-M2`, `-F`) — written as `ST-########`. It is recorded for every sample in the
> manifest column `patient_id_if_available`.
>
> This key is **not** a mere convention. In GSE281005, where it is published as
> `public patient_id`, the map between specimen core and patient identifier is **strictly
> 1:1 across all 290 cores**: no core carries more than one patient identifier and no
> patient identifier owns more than one core. The same 1:1 relation therefore holds a
> fortiori over GSE205154.
>
> **The effective unit of analysis is the patient, not the sample.**

**Coverage and its limit.** 202 of the 289 GSE205154 samples sit in one of the 194 cores
that also appear in GSE281005, so their key is the *published* `public patient_id`.
The other 87 samples sit in one of the 83 cores with no GSE281005 counterpart; their key
is the same `ST-########` string read off the sample title, without a published
counterpart. For those 83 cores we cannot exclude that two *differently numbered* cores
belong to one patient — the deposit exposes no way to detect it. Across the 290 cores
where both series overlap there is no instance of core-splitting, which is evidence
against it but not proof for the unmatched 83.

**Interpretation limits that follow, and must be carried into any write-up:**

1. The unit is the patient, but for 83 of 277 cores that is a *within-deposit* key rather
   than an externally published identifier. A reader cannot independently re-derive the
   87 samples' patient grouping from GSE205154 alone.
2. `Sex` and `age at_diagnosis` exist only for the 202 samples whose core appears in
   GSE281005; they are absent for the other 87. Neither field enters any frozen
   eligibility rule, so nothing in the arithmetic depends on them.
3. Patient-level overlap between GSE205154 and **GSE282302** (Visium, same institution)
   **cannot be assessed**: GSE282302 titles are ROI identifiers with no shared namespace
   (GSM intersection 0, specimen-core overlap 0, title overlap 0, but that is not proof of
   patient disjointness). This does not affect the arithmetic — GSE282302 contributes no
   bulk expression and is not a bulk meta-analysis cohort.
4. Patient-level overlap with **TCGA-PAAD** is not assessable either: no shared accession
   namespace and a different institution.
5. GSE205154 and GSE281005 **are** companion SubSeries over the same patients and share
   170 verbatim sample titles. GSE281005 is TCRβ sequencing and contributes no expression
   data, so it is not a cohort and cannot duplicate this one — it is the source of the
   patient key.

---

## 4. Rule for duplicate / replicate specimens from the same patient

> **FROZEN — duplicate/replicate rule, written down before any effect was examined.**
>
> **D4.1** Eligibility is assessed **within a compartment**. A patient may legitimately
> appear once in the primary compartment and once in the metastatic compartment; those are
> two different analyses and are never pooled, so this is not a duplicate.
>
> **D4.2** If a patient contributes **more than one specimen inside the same compartment**,
> exactly **one** of them enters the analysis. The retained specimen is the **base-label
> specimen** — the specimen whose title carries no numeric replicate suffix, i.e. `-T`
> rather than `-T2`, `-M` rather than `-M2`. If a core contains no base-label specimen,
> the lexicographically smallest `sample_title` in that core is retained.
>
> **D4.3** Expression values are **not** averaged across duplicate specimens. Selection was
> chosen over aggregation because aggregation would introduce a specimen-combination rule
> that no existing cohort in this manuscript uses, and the frozen text requires only that
> *a* rule be fixed in advance.
>
> **D4.4** The choice is outcome-blind by construction: D4.2 depends only on the sample
> title, so it cannot be influenced by any expression value.
>
> **D4.5 Pre-specified sensitivity (must be reported whatever it shows):** the main
> analysis is repeated with the duplicate pair **averaged at the expression level**
> (per-gene arithmetic mean over the two specimens of that core), and with both members
> retained. Neither variant may replace the primary result.

**What D4 actually removes in this deposit.** Only two cores contain more than one specimen
in the same compartment, and both are primary:

| core | specimens | compartment | retained | excluded |
|---|---|---|---|---|
| `ST-00017945` | `ST-00017945-T`, `ST-00017945-T2` | primary | `ST-00017945-T` (GSM6206194) | `ST-00017945-T2` (GSM6206193) |
| `ST-00019598` | `ST-00019598-T`, `ST-00019598-T2` | primary | `ST-00019598-T` (GSM6206244) | `ST-00019598-T2` (GSM6206245) |

**Cross-compartment pairs are not duplicates.** Ten cores contribute one primary and one
metastasis each (`ST-00004898`, `ST-00006291`, `ST-00007307`, `ST-00010984`,
`ST-00014524`, `ST-00015839`, `ST-00017804`, `ST-00017838`, `ST-00019429`, `ST-00019601`).
The primary member enters the main analysis and the metastatic member enters the
sensitivity analysis. Because the two analyses are never pooled, no within-analysis
duplication arises. The consequence is that **10 patients appear in both analyses**; an
**overlap-free sensitivity subset of n = 61** is pre-specified for reporting alongside the
full 71.

**Metastatic compartment duplication: none.** All 71 metastatic specimens sit in 71
distinct cores, so D4.2 removes nothing there.

---

## 5. Frozen eligibility checks, condition by condition

| condition | frozen rule | result for GSE205154 |
|---|---|---|
| **Disease** | histologically identified PDAC | **PASS** — series title/summary/design and PMID 39789181 all PDAC; pathologist-marked tumour regions; all 218 primary specimens pancreatic, all 71 metastases PDAC metastases |
| **Specimen** | tumour tissue; primary preferred for the main meta-analysis | **PASS** — 289/289 are tumour (`cell type: scraped tumor`, >20 % tumour nuclei by pathology review); 218 primary → main, 71 metastatic → prespecified sensitivity |
| **Unit** | one independent tumour/patient; duplicate specimens need a pre-defined aggregation/selection rule | **PASS after D4** — patient key established (§3); 277 cores; 2 within-compartment duplicate cores resolved by D4.2 → 216 primary specimens |
| **SIGLEC15** | interpretable gene-level expression | **PASS** — `hgnc_symbol = SIGLEC15`, ENSG00000197046, `protein_coding`, exactly one row (index 17332); `TMM_pass_filter = TRUE` |
| **ACP5 / MMP9** | both measurable, score constructible by the established rule | **PASS** — ACP5 ENSG00000102575 (row 2686) and MMP9 ENSG00000100985 (row 2418), both `protein_coding`, both `TMM_pass_filter = TRUE`; the frozen rule requires exactly these two genes and both are present |
| **Metadata** | sufficient to distinguish primary / metastatic / non-tumour | **PASS** — `tumor type` populated for 289/289, corroborated by `Sample_description` and the title suffix; 0 conflicts; no non-tumour samples exist in the series |
| **Platform** | within-cohort rank information obtainable; cross-platform absolute comparability NOT required | **PASS** — single platform GPL20301, all 289 samples in one gene-level matrix with identical sample-column order across the Counts, TMM and TPM matrices; rank information is available within the cohort. **Cross-platform absolute comparability is neither required nor claimed** |
| **Analysis** | cohort-specific Spearman → random-effects meta-analysis | **DEFERRED** — not performed in this step |
| **Exclusion** | never exclude a cohort because of the observed ρ magnitude or direction | **ACKNOWLEDGED AND BINDING** — no ρ was computed in this step, so this condition could not have operated |

### 5.1 Overlap with the eight named datasets — verified, not assumed

Checked against the 289 GSMs, the 277 specimen cores and the 289 titles of GSE205154:

| dataset | series n | GSM overlap | specimen-core overlap | title overlap |
|---|---|---|---|---|
| GSE155698 | 41 | 0 | 0 | 0 |
| GSE156405 | 17 | 0 | 0 | 0 |
| GSE212966 | 12 | 0 | 0 | 0 |
| GSE217845 | 16 | 0 | 0 | 0 |
| GSE282302 | 108 | 0 | 0 | 0 |
| GSE57495 | 63 | 0 | 0 | 0 |
| GSE85916 | 80 | 0 | 0 | 0 |
| TCGA-PAAD | n/a | not assessable | not assessable | not assessable |

No overlap is identified with any comparator. The TCGA-PAAD row is not assessable because
GSE205154 samples carry no TCGA barcode and no shared accession namespace exists; the
cohorts are from different institutions. See §3.2 limits 3 and 4 for the two overlap
questions that public metadata cannot settle.

---

## 6. FINAL n — the numbers this manifest fixes

### 6.1 Arithmetic

```
289  deposited samples  (289 GSM = 289 BioSample = 289 distinct titles)
  0  non-tumour samples           all 289 are pathologist-marked tumour tissue
289  tumour samples
       split on the frozen compartment field 'tumor type':
218  primary tumours              -> main meta-analysis compartment
 71  metastatic lesions           -> prespecified sensitivity compartment

MAIN META-ANALYSIS (primary tumours, one specimen per patient)
218  primary specimens
 -2  second primary specimen from an already-represented patient (D4.2)
       ST-00017945-T2  (GSM6206193)  -> retained ST-00017945-T  (GSM6206194)
       ST-00019598-T2  (GSM6206245)  -> retained ST-00019598-T  (GSM6206244)
 -0  gene-usability exclusions
===
216  FINAL n for the main meta-analysis

PRESPECIFIED SENSITIVITY ANALYSIS (metastatic lesions, one specimen per patient)
 71  metastatic specimens
 -0  within-compartment duplicate patients
===
 71  FINAL n for the prespecified sensitivity analysis
       of which 61 are from patients that do not also appear in the main analysis
       (the overlap-free subset, pre-specified for reporting)
```

### 6.2 Frozen totals

| quantity | frozen value |
|---|---|
| **final n, main meta-analysis (primary tumours)** | **216** |
| **final n, prespecified sensitivity analysis (metastatic lesions)** | **71** |
| overlap-free sensitivity subset | 61 |
| distinct patients (whole deposit) | 277 |
| samples eligible for nothing | 2 (`GSM6206193`, `GSM6206245`) |
| samples eligible for both analyses | 0 |

**No pooled or meta-analytic sample size is stated here.** What the other cohorts
contribute, and what any combined estimate would rest on, is outside what this manifest
establishes and is deliberately not pre-written.

---

## 7. Frozen downstream rule that the ρ step must obey

These are **not** choices made here; they are transcribed from the existing analysis
identity so that the later step cannot silently re-specify them.

**Score construction — `analysis_F_cross_scale.py` (v39 analysis archive), functions
`zscore_columns()` and `score_program()`:**

```
PROGRAMS["Derosiers_ACP5_MMP9_core"] = ["ACP5", "MMP9"]
required_gene_count(PRIMARY)         = 2          # both genes must be present
z                                    = (x - mean(x)) / sd(x)   over SAMPLES, ddof = 0
score                                = nanmean(z[ACP5], z[MMP9])   # the "mean-z" method
```

The standardisation is **within-cohort** (each gene is centred and scaled across the
cohort's samples), which is what makes the score comparable across cohorts without
requiring cross-platform absolute comparability.

**Statistic:** Spearman ρ between tumour-level SIGLEC15 and the score, computed
**separately per cohort**.

**Combination — `random_effects_meta()`:** Fisher z = `arctanh(clip(ρ, ±0.999999))`;
`v = 1/(n − 3)`; τ² by REML; modified Hartung–Knapp interval with the scale factor floored
at 1 (`q_HK = max(1, Σw(y−μ)²/(k−1))`, `se = sqrt(q_HK/Σw)`, t critical with `df = k − 1`);
back-transform by `tanh`. Seed 15072026, 10 000 resamples where resampling is used.

**Compartments:** the primary-tumour ρ and the metastatic ρ are **separate cohort-level
estimates and are never pooled into one dataset-level ρ.**

---

## 8. AUTHOR DECISIONS — flagged, not buried

The frozen text did not settle the following. Each was decided explicitly here, before any
effect was examined, and each must be reported as an author decision rather than presented
as a pre-registered rule.

### AD-1 — Patient key = `ST-########` specimen core (adopted)

*Decision:* use the specimen core as the patient key for all 289 samples, as set out in §3.2.
*Reason:* in GSE281005 the same key is published as `public patient_id` and maps 1:1 onto
290 patient identifiers with no exception, so it is an identity key and not a naming
convention.
*Consequence if the author disagrees:* treating the sample as the unit would admit the two
`-T2` specimens as independent tumours and give **n = 218** instead of 216 for the main
analysis. That variant is recorded here so the author can rule on it, but it is **not** the
frozen rule.

### AD-2 — Which deposited matrix feeds the programme score (PRIMARY: log2(TMM + 1))

*Decision:* the primary representation is the deposited **TMM** matrix (edgeR method,
library-composition normalised, not length normalised), transformed as **log2(x + 1)**
before the frozen within-cohort z-standardisation.
*Reason:* the three existing cohorts all entered `score_program()` on a **log2 scale** —
TCGA-PAAD as `log2(CPM + 1)`, GSE57495 and GSE85916 as log2 microarray intensities — and
the mean-z rule was developed and applied in that input domain. Of the three representations
this deposit supplies, TMM is the closest structural analogue of CPM: library-composition
normalised but **not** length normalised, so `log2(TMM + 1)` reproduces the input identity
of `log2(CPM + 1)`. Standardising untransformed counts would put the z-scores on a
different, outlier-dominated scale from every comparator cohort and would silently change
what the programme score measures.
*Why this is a flagged decision:* the frozen text says gene expression is "used as
supplied" for bulk and is silent about which of the three supplied matrices to choose. This
is that choice.
*Scope of the decision:* it can move the **ACP5/MMP9 score** side only. SIGLEC15 enters the
Spearman correlation through its **ranks**, which are invariant under any strictly monotone
transform of the same matrix, so no part of this decision touches the SIGLEC15 side.
*Pre-specified sensitivities (must be reported whatever they show, and must never replace
the primary):* (i) **log2(TPM + 1)**; (ii) **TMM as supplied, no log transform**. All three
matrices carry SIGLEC15, ACP5 and MMP9 with identical row indices, so all three variants
are runnable.

### AD-3 — Duplicate rule: selection by base specimen label (adopted)

*Decision:* D4.2 in §4 — retain the base-label specimen, do not average.
*Reason:* selection needs no new aggregation identity; averaging would introduce a
specimen-combination rule no comparator cohort uses. The rule is outcome-blind because it
reads only the sample title.
*Pre-specified sensitivity:* expression-level averaging of the duplicate pair (§4, D4.5).

### AD-4 — Overlap-free sensitivity subset reported alongside the full 71

*Decision:* report the metastatic ρ at n = 71 **and** at n = 61 (the ten patients that also
contribute a primary to the main analysis removed).
*Reason:* it costs nothing and it pre-empts the question of whether the two compartment
analyses are patient-disjoint. The frozen text does not require it, so it is flagged.

### AD-5 — `F`-suffix samples counted as primary

*Decision:* the seven `F`-suffix samples are treated as primary tumours.
*Reason:* all three independent metadata fields agree (`tumor type: Primary`,
`tissue: Pancreas`, `description: Primary pancreas tumor`), and the frozen rule reads
`tumor type` rather than the title suffix. The deposit never defines what `F` denotes; the
meaning remains unknown (§9). No eligibility decision depends on resolving it.

---

## 9. What could not be determined

1. **The meaning of the `F` title suffix** (7 samples). All three metadata fields are
   concordant, so this affects nothing, but the code is undocumented in the deposit.
2. **Whether two differently-numbered specimen cores belong to one patient**, for the 83
   cores with no GSE281005 counterpart (87 samples). The 290-core 1:1 mapping is evidence
   against core-splitting but is not proof for the unmatched cores.
3. **Patient-level overlap with GSE282302** (same institution, Visium). Not assessable from
   public metadata; does not affect this cohort's arithmetic.
4. **Patient-level overlap with TCGA-PAAD.** No shared accession namespace.
5. **Raw sequencing data.** Controlled access via dbGaP `phs003597.v1.p1`; not obtained and
   not needed for the frozen protocol.
6. **Sex and age for 87 samples.** Published only through GSE281005; neither field enters
   any frozen eligibility rule.

---

## 10. Attestation

* The manifest `sample_manifest.csv` was built by `scripts/10_build_manifest.py` from
  metadata and feature-identifier information only.
* Every one of the 289 rows carries an explicit `yes`/`no` on `eligible_for_main_meta` and
  on `eligible_for_sensitivity`, and every `no` carries a written `exclusion_reason`.
  No sample is silently skipped.
* No Spearman correlation between SIGLEC15 and any ACP5/MMP9 programme score was computed,
  inspected or reported. No per-sample expression value of SIGLEC15, ACP5 or MMP9 was
  loaded.
* The SHA-256 recorded at the top of this document is the hash of the frozen manifest. Any
  later edit to the manifest invalidates this definition and must be recorded as a new
  version with its own hash.
