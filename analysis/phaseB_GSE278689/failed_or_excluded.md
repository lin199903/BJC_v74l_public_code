# failed_or_excluded.md — GSE278689, Phase B step 2

> **On this filename.** `failed_or_excluded.md` is the fixed deliverable name set by the run
> brief. Nothing in this run is characterised as a breakdown of any analysis. Every sample
> either **passes** the frozen eligibility rules or is **excluded with a stated reason**; every
> convention question is either resolved by a frozen constant or **flagged as a fork**.
>
> Per the author's pre-registered interpretation criteria, the four words named in those
> criteria are reserved for nothing: they are not applied to any sample, estimate, interval or
> contrast in this document. Where a verification control is recorded as meeting or not meeting
> its criterion, that is a technical verdict on the control, never a characterisation of a
> study result.

Scope: GSE278689 only, snRNA (nuclei), cross-modality sensitivity analysis, single stratum.
Never pooled with the primary 19-sample scRNA estimate or with GSE278688.

---

## 1. Sample-by-sample outcome (no silent exclusions)

Source: `output/sample_manifest.csv`, produced by `scripts/07_eligibility_check.R`.

| sample_id | patient | tissue | total cells | myeloid | detectable | non-detectable | ≥50 | both groups | final | reason |
|---|---|---|---|---|---|---|---|---|---|---|
| 2107809_FFPE | PA21 | Tumor | 4,684 | 820 | 13 | 807 | pass | pass | **included** | — |
| 2021124_FFPE | PA24 | Tumor | 6,816 | 1,346 | 11 | 1,335 | pass | pass | **included** | — |
| 2049529_FFPE | PA29 | Tumor | 1,604 | 235 | 3 | 232 | pass | pass | **included** | — |
| 1869485_FFPE | PA85 | Tumor | 6,109 | 194 | 1 | 193 | pass | pass | **included** | — |
| 1942191_FFPE | PA91 | Tumor | 3,113 | 520 | 6 | 514 | pass | pass | **included** | — |
| 2151198_FFPE | PA98 | Tumor | 6,336 | 1,426 | 10 | 1,416 | pass | pass | **included** | — |
| 2048540_FFPE | PA40 | Tumor | 1,693 | 252 | **0** | 252 | pass | **excluded** | excluded | 0 of 252 tumour-myeloid cells carry a non-zero SIGLEC15 count; the detectable group is empty and the primary contrast requires both groups |
| 2051072_FFPE | PA72 | Tumor | 1,650 | 598 | **0** | 598 | pass | **excluded** | excluded | 0 of 598 tumour-myeloid cells carry a non-zero SIGLEC15 count; the detectable group is empty |

Every one of the 8 deposited samples appears in the manifest with an explicit verdict. No
sample is dropped without a reason, and no sample that meets the rules is dropped.

**Patient overlap / duplicates / integrated objects:** 8 samples, 8 distinct patients, one
prefix per patient and one patient per prefix — no overlap to exclude. No duplicate
`sample_id`. No integrated or batch-corrected object was read; only the deposited count matrix.

---

## 2. Forks that were flagged rather than resolved silently

### 2.1 Normalisation (flagged; manuscript convention adopted)
The deposit is **raw integer counts**, so the brief's first branch applies and CP10K + log1p is
used. The alternatives are reported because the choice is not immaterial:

| convention | mean contrast | Δ vs adopted | max single-sample Δ |
|---|---|---|---|
| **CP10K + log1p (adopted)** | **+0.9175** | — | — |
| CP10K without log1p | +1.2463 | +0.3287 | 1.5866 |
| log1p(counts), no library scaling | +1.2813 | +0.3638 | 1.2163 |
| raw counts, z-scored | +1.3022 | +0.3847 | 1.0635 |

The adopted convention is the one the manuscript's own pipeline uses
(`np.log1p(Mm / lib * 1e4)` in `Scripts/e2b_independent_check.py`). Source:
`output/_normalisation_fork_test.json`, `scripts/05_normalisation_fork_test.py`.

### 2.2 Aggregation of the specimen-level myeloid summary (flagged; both reported)
`log1p(mean(cpm))` (the only rule the manuscript writes down) versus `mean(log1p(cpm))` are
**not** monotone-linked, so the myeloid Spearman moves **+0.4286 → +0.3714 (Δ 0.0571)**. The
all-cell aggregate is exactly unaffected (ρ +0.9429 under both). Both values are in
`output/specimen_level_correlations.csv`. Source: `scripts/11_aggregation_convention.R`.

### 2.3 Interval rule (flagged)
The percentile bootstrap CI (+0.0392 to +2.3055) clears zero; the **normal-approximation
(−0.4030 to +2.2380)** and **basic bootstrap (−0.4736 to +1.7959)** intervals of the same six
values both include zero. The frozen brief specifies the bootstrap interval, so the percentile
interval is primary, but the disagreement is a real property of this estimate at n = 6.
Source: `scripts/10_ci_sanity.R`.

---

## 3. A defect found in this run and corrected

`scripts/05_normalisation_fork_test.py`, first revision, passed **all three** captured columns
(`SIGLEC15`, `ACP5`, `MMP9`) into the programme score. That is precisely the target-in-program
error the frozen brief forbids, and it inflated the mean contrast from **+0.9175 to +4.0071**.

It was caught by the `WRONG 6` control in `scripts/09_independent_verification.py`, which
deliberately commits the same error and is *expected* to move the answer. The script was
corrected to assert that the programme matrix has exactly two columns, and the reference branch
now reproduces the R primary result to 8.9 × 10⁻¹⁵. The superseded first-revision log is
retained in `logs/05_normalisation_fork_20260928_072101.log`; the corrected run is
`logs/05_normalisation_fork_FIXED_20260928_072355.log`.

Similarly, first revisions of `10_ci_sanity.R` and `11_aggregation_convention.R` contained two
printed claims that their own checks contradicted (that every interval rule agrees on the sign,
and that the aggregation rules are rank-equivalent). Both claims were wrong and were replaced
with the measured result.

---

## 4. Rules deliberately NOT applied

* **The restricted TREM2_TAM rule** (`≥10 cells with ≥2 detectable and ≥5 non-detectable`). The
  frozen brief assigns it to a separate analysis. It was **not** applied to the primary
  contrast, and no restricted result is reported. This matters: **2 of the 6 estimable samples**
  (PA29 with 3 detectable, PA85 with 1) would be dropped by that rule. Applying it would give
  n = 5 and mean +1.1704 — recorded only as a deliberately-wrong control in
  `output/_verification_independent.json`, never as a result.
* **Re-annotation of any kind.** The deposited `all_celltype` column is used verbatim. Marker
  genes appear only in verification (step 06), never to call a cell type.
* **Any second modality or dataset.** Only GSE278689 files were read for analysis.
* **Protein inference.**

---

## 5. Not determinable from this run

1. Anything about whole-cell (scRNA) measurement. The detectability comparison is a
   between-dataset, between-modality statement; no pooling and no combined estimate is formed.
   The denominators were not equated (published fraction = 197/9,918; here = 44/4,541 over the
   6 estimable samples).
2. Whether these six samples are comparable to the manuscript's nineteen — they are not
   like-for-like, being tumour-only by the frozen rule and about half the detectable density.
3. Any restricted TREM2_TAM result (§4).
4. Whether an invariant programme gene would change anything — **no eligible sample had one**;
   all six used ACP5 and MMP9. `output/programme_gene_variance_by_sample.csv` records the
   within-sample SDs that establish this.
