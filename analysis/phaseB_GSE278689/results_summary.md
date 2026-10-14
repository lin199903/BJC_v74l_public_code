# results_summary.md — GSE278689, Phase B step 2

**Scope.** GSE278689 only. **Modality: snRNA (single nuclei).** This is a **cross-modality
sensitivity analysis** and a **single stratum**. It is **not pooled** with the manuscript's
19-sample scRNA estimate and **not pooled** with GSE278688; no combined estimate is formed
anywhere in this run. It answers one question only: *does the target–state relation remain
estimable, and in what direction, under a different transcript-capture process?*

Frozen constants used verbatim: eligible label `tissue == "Tumor"`; cell type = deposited
`all_celltype == "Myeloid"` (no re-annotation); ≥50 tumour-myeloid cells per sample; **both**
detectable and non-detectable groups required with **no minimum group size**; programme score
= **ACP5 + MMP9, both required**, each z-standardised within the biological sample (ddof = 0)
then averaged; **SIGLEC15 never enters the programme score**; samples equally weighted;
detection = non-zero **raw** count; CP10K + log1p; seed **15072026**; **10,000** resamples.

---

## Detectability — the number that governs the reading

The deposit is **raw integer UMI counts** (35,060,262 entries, 0 non-integer, 0 negative), so
the frozen first branch applies: CP10K + log1p. Detectability is read from the raw count and is
invariant to that transform (all four candidate normalisations reproduce the identical split).

| Measurement process | SIGLEC15-detectable / tumour-myeloid | Fraction |
|---|---|---|
| Published scRNA series (primary, 19 samples) | 197 / 9,918 | **2.00%** |
| GSE278688 — scRNA, whole cells | 52 / 9,102 | **0.57%** |
| **GSE278689 — snRNA, nuclei (this run), 6 estimable samples** | **44 / 4,541** | **0.97%** |
| GSE278689 — all 8 deposited samples | 44 / 5,391 | 0.82% |

Nuclear detectability is **≈2.1-fold lower** than the published scRNA series (ratio 0.48) and
≈1.7-fold higher than GSE278688. Across every cell type in this series, 157 / 32,005 cells
(0.49%) carry a non-zero SIGLEC15 count.

**This is a modality/dataset difference, not a quality ranking.** It is the reason the
permitted statement for this run is that **the relation was difficult to estimate under a
lower-detectability nuclear RNA measurement process**. It is *not* evidence about cellular
biology.

---

## Result

**Exclusion arithmetic (8 → final estimable n).**

```
 8  samples deposited
-0  ineligible tissue label (all 8 are tissue == "Tumor")
 8  tumour-labelled samples
-0  patient overlap between samples     (8 distinct patients, 1:1 with the 8 samples)
-0  duplicate samples
-0  integrated objects
-0  fewer than 50 tumour-myeloid cells  (minimum 194)
-2  SIGLEC15-detectable group empty     (PA40/2048540_FFPE: 0 of 252; PA72/2051072_FFPE: 0 of 598)
 6  = FINAL ESTIMABLE n
```

**Final estimable n = 6.** No invariant programme gene occurred in any eligible sample, so all
six used **ACP5 and MMP9**.

| sample_id | patient | detectable | non-detectable | contrast |
|---|---|---|---|---|
| 2107809_FFPE | PA21 | 13 | 807 | **+4.2017** |
| 2151198_FFPE | PA98 | 10 | 1,416 | +0.7915 |
| 1942191_FFPE | PA91 | 6 | 514 | +0.4036 |
| 2049529_FFPE | PA29 | 3 | 232 | +0.2368 |
| 2021124_FFPE | PA24 | 11 | 1,335 | +0.2183 |
| 1869485_FFPE | PA85 | 1 | 193 | −0.3467 |

**Mean contrast +0.9175**, 95% bootstrap CI **+0.0392 to +2.3055** (10,000 resamples, seed
15072026, sample-level within the single stratum); signed-rank p = 0.1159; sign test
p = 0.2188; **positive direction in 5 of 6 samples**; median +0.3202.

> A positive contrast was observed in a modality-specific sensitivity analysis.

**Four same-sample specimen-level SIGLEC15 summaries** (manuscript inverse-log aggregation),
each against the sample-specific contrast, Spearman with a sample bootstrap interval:

| Summary | ρ | 95% CI | p |
|---|---|---|---|
| Aggregated myeloid SIGLEC15 RNA | **+0.4286** | −0.8000 to +1.0000 | 0.397 |
| Detectable-cell fraction | **+0.6000** | −0.8000 to +1.0000 | 0.208 |
| Mean among detectable cells | **+0.0286** | −1.0000 to +1.0000 | 0.957 |
| All-cell tumour aggregate | **+0.9429** | +0.5148 to +1.0000 | 0.005 |

At n = 6 Spearman can only take k/15, so these are +6/15, +9/15, +0/15 and +14/15. The first
three intervals span most of the range and include zero.

---

## Stable in this dataset

1. **The eligibility arithmetic and the manifest.** 32,005 barcodes, 32,005 metadata rows, zero
   duplicates, barcode multisets identical, sample prefixes 1:1 with patients. All 8 samples are
   `Tumor`; the minimum tumour-myeloid count is 194. Every exclusion is named with a reason.
2. **The deposit's representation and the detectability split.** Two independent parsers
   (chunked pandas reader, `scipy.io.mmread`) agree entry by entry. Values are integral, so the
   CP10K+log1p branch is forced, and the detectable/non-detectable split is identical under all
   four candidate normalisations — eligibility and every group size are unambiguous.
3. **The gene→row mapping.** Six scale-free positive controls on `all_celltype`, no margin
   constant: KRT19 top in Ductal, EPCAM Ductal, COL1A1 CAF (5.47 vs runner-up 2.51), PTPRC in
   T cells, CD68 Myeloid (2.07 vs 0.53), CPA3 Myeloid. **6/6 meet the criterion**; a
   deliberately rotated gene→column assignment meets 2/6, so the checks discriminate.
4. **The implemented contrast.** An independent Python re-implementation that imports the
   manuscript's own `zscore_columns()` and `score_program()` reproduces all six sample contrasts
   to **1.8 × 10⁻¹⁵**; the CI and p reproduce across languages and RNG streams.
5. **The direction of the contrast under the frozen rule.** 5/6 positive; the percentile
   bootstrap interval lies above zero under every seed tried (lower bound +0.034 to +0.039).

## Unstable in this dataset

1. **How much uncertainty the interval really carries.** The percentile bootstrap CI clears zero
   by +0.039, but the **normal-approximation and basic-bootstrap intervals of the same six
   values both include zero** (−0.40 and −0.47). The exact bootstrap (46,656 resamples, no RNG)
   puts P(resampled mean ≤ 0) at **1.5%**. The margin is narrow, not comfortable.
2. **The point estimate is driven by one sample.** PA21/2107809_FFPE contributes +4.2017 against
   a median of +0.3202; without it the mean falls to **+0.2607**. Dropping the single negative
   sample (PA85/1869485_FFPE) moves the mean to **+1.1704**.
3. **All four specimen-level correlations.** Only the all-cell aggregate has an interval
   excluding zero, and at n = 6 these estimates carry very little information (135/10,000
   bootstrap resamples land on exactly ρ = −1).
4. **The magnitude is normalisation-dependent.** Under the deposit-consistent alternatives the
   mean moves from +0.9175 to +1.2463 (CP10K without log1p), +1.2813 (log1p of counts, no
   library scaling) and +1.3022 (raw counts z-scored) — up to **+0.3847 on the mean** and up to
   1.59 for a single sample. The manuscript convention is used and every alternative is reported.
5. **The aggregation rule for the myeloid summary.** `log1p(mean(cpm))` versus `mean(log1p(cpm))`
   are not monotone-linked, and the myeloid ρ moves from **+0.4286 to +0.3714 (Δ 0.0571)**. Both
   values are reported. The all-cell aggregate is exactly unaffected.

## Not determinable

1. **Anything about whole cells.** The comparison above is a between-dataset, between-modality
   statement. No cross-modality pooling, no combined estimate, and no inference about cellular
   biology is made or supported by this run.
2. **Whether these six samples are comparable to the manuscript's nineteen.** They are not
   like-for-like: the eligible set is tumour-only by the frozen rule and the detectable-cell
   density is about half. The denominators were not equated.
3. **Any restricted TREM2_TAM sensitivity result.** The `≥10 cells with ≥2 detectable and ≥5
   non-detectable` rule belongs to a separate analysis and was **not** applied; two of the six
   estimable samples have fewer than five detectable cells and are retained by the primary rule
   precisely because it sets no minimum.
4. **Protein status of any kind.** Never inferred.

---

*Seed 15072026 · 10,000 bootstrap resamples · single stratum · snRNA (nuclei).*
*Every result file is paired with the script that produced it in `run_log.txt`; environment in
`sessionInfo.txt`; exclusions, forks and defects in `failed_or_excluded.md`.*
