# results_summary.md — GSE278688, Phase B step 1

**Scope:** GSE278688 only. Frozen rules applied verbatim. Seed **15072026**, **10,000** bootstrap
resamples. Full command/version trail in `run_log.txt` and `sessionInfo.txt`.

**Final estimable n = 8.**
29 deposited → −16 non-tumour tissue labels (`Adjacent_normal`, `PBMC`) → 13 tumour samples
→ −0 (<50 tumour-myeloid: none failed; minimum was 80) → −5 with zero SIGLEC15-detectable
tumour-myeloid cells → **8**. No patient contributed more than one eligible sample; no duplicated
sample; no invariant programme gene in any eligible sample (so all 8 scores used **ACP5 and MMP9**).

---

## Stable in this dataset

1. **The eligibility arithmetic itself.** The 29-sample manifest reconciles exactly with the
   barcode and metadata files (199,184 cells; barcode multiset identical, 0 duplicates,
   `setequal` TRUE). Tissue and patient labels parsed from barcodes agree with the metadata
   columns for all 29 samples. Verified in `07_eligibility_check.R`.
2. **SIGLEC15 detectability.** 123/199,184 cells are SIGLEC15-positive; 52 of these are myeloid
   (0.57% of 9,102 myeloid cells). Detectability is **invariant to the normalisation choice** —
   every candidate normalisation is monotone with f(0)=0 — so the eligible set and all
   detectable/non-detectable group sizes are unambiguous.
3. **The gene→row mapping.** Positive controls behave: KRT19 Ductal 33.73 (next 0.29),
   COL1A1 CAF 81.59, CD68 Myeloid 7.16, CPA3 Myeloid 7.67, ACP5 highest in Myeloid.
4. **The implemented contrast.** An independent Python re-implementation that imports the
   manuscript's own `zscore_columns()` reproduces all 8 sample contrasts to **6.7 × 10⁻¹⁶**
   (positive control: the two rejected conventions differ by up to 0.76, so the agreement is
   not vacuous).
5. **Direction of the programme contrast.** Mean contrast is positive, **+0.252**.

## Unstable in this dataset

1. **The mean contrast and its interval.** +0.252 (95% CI **−0.201 to +0.812**); signed-rank
   p = 0.484; sign test p = 1.00; **4/8 samples positive**. The CI crosses zero and the
   point estimate is driven by one sample (PDAC06-Tumor, +1.987 versus a median of +0.083).
   Excluding that sample alone would move the mean to about −0.00.
2. **All four specimen-level correlations.** ρ = +0.048, +0.095, −0.071, +0.262 for aggregated
   myeloid RNA, detectable-cell fraction, mean among detectable cells and all-cell tumour
   aggregate respectively. Every CI spans nearly the full range (e.g. −1.000 to +0.848) and no
   interval excludes zero. **These estimates carry essentially no information at n=8**: the
   attainable Spearman grid is k/28, and 319/10,000 bootstrap resamples land on exactly −1.0.
3. **The normalisation-dependence of the magnitude.** The deposit is **raw integer UMI counts**,
   not stored log1p-normalised values. Under the manuscript's own documented pipeline
   (normalise to 10,000 CPM, then log1p) the mean contrast is **+0.252**; without library-size
   scaling (`log1p(counts)`) it is **+0.388**. The gap (0.136; up to 0.76 for a single sample)
   is a real fork. The manuscript's convention was adopted and is recorded in every output file,
   but it was **not** specified by the frozen brief — see `failed_or_excluded.md` §3.
4. **Effect of the z-score degrees-of-freedom.** ddof=0 (manuscript) vs ddof=1 changes the mean
   contrast by 0.0004 — immaterial here, but recorded because it is not zero.

## Could not be determined

1. **Whether the contrast replicates.** Not answerable from one dataset; no claim is made. The
   comparison is in `output/primary_summary.csv` and is stated numerically in the message to the
   requester.
2. **The all-cell tumour aggregate aggregation rule — but it does not matter here.** The
   manuscript's SI §15 defines the measure as "every cell in the retained tumour samples" without
   stating the aggregation, whereas the myeloid aggregate uses inverse-log. Both candidates were
   computed: a plain CP10K mean and inverse-log give the **identical** ρ = +0.2619 (p = 0.53),
   because Spearman depends only on ranks and the two are monotone transforms of each other.
   Only the bootstrap interval differs slightly (−0.924 to +0.974 vs −0.803 to +0.918), since
   resampling can break the monotone relation. The plain CP10K mean is reported in the table for
   unit-consistency with the other three summaries. See `scripts/11_allcell_convention.R`.
3. **Any TREM2_TAM restricted-sensitivity result.** The frozen brief explicitly separates the
   `≥10 cells with ≥2 detectable and ≥5 non-detectable` rule into a different analysis. It was
   **not** applied here and no such result is reported.
4. **Whether these 8 samples are comparable to the manuscript's 19.** They are not like-for-like
   in composition: the eligible set is tumour-only by the frozen rule, and detectable-cell
   density is lower here (52/9,102 = 0.57% of myeloid cells) than in the published series
   (197/9,918 = 2.0%). The relevant denominators differ and were not equated.
5. **Protein status.** Never inferred, per the frozen rule.
