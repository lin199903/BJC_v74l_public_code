# failed_or_excluded.md — GSE278688, Phase B step 1

Two sections: **(A)** every excluded sample with its concrete reason, and **(B)** everything this
run attempted that failed or that could not be resolved as specified.

---

## A. Excluded samples — full accounting (29 → 8)

Machine-readable in `output/sample_manifest.csv` (all 29 rows, with `exclusion_reason`).

### A1. Excluded on tissue label — 16 samples

Frozen rule: eligible = `tissue == "Tumor"`; `Adjacent_normal` and `PBMC` are not eligible.
None of these was evaluated for cell count; the label excludes them first.

| sample_id | patient | tissue | myeloid cells | reason |
|---|---|---|---|---|
| PDAC04-Normal | PA04 | Adjacent_normal | 61 | not primary PDAC tumour |
| PDAC05_PBMC | PA05 | PBMC | 41 | not primary PDAC tumour |
| PDAC06-Normal | PA06 | Adjacent_normal | 489 | not primary PDAC tumour |
| PDAC06-PBMC | PA06 | PBMC | 262 | not primary PDAC tumour |
| PDAC07-PBMC | PA07 | PBMC | 178 | not primary PDAC tumour |
| PDAC08-Normal | PA08 | Adjacent_normal | 146 | not primary PDAC tumour |
| PDAC08-PBMC | PA08 | PBMC | 589 | not primary PDAC tumour |
| PDAC09-Normal | PA09 | Adjacent_normal | 231 | not primary PDAC tumour |
| PDAC09-PBMC | PA09 | PBMC | 20 | not primary PDAC tumour |
| PDAC10-Normal | PA10 | Adjacent_normal | 772 | not primary PDAC tumour |
| PDAC10-PBMC | PA10 | PBMC | 105 | not primary PDAC tumour |
| PDAC11-PBMC | PA11 | PBMC | 103 | not primary PDAC tumour |
| PDAC12-PBMC | PA12 | PBMC | 94 | not primary PDAC tumour |
| PDAC13-PBMC | PA13 | PBMC | 174 | not primary PDAC tumour |
| PDAC14-PBMC | PA14 | PBMC | 86 | not primary PDAC tumour |
| PDAC15-PBMC | PA15 | PBMC | 97 | not primary PDAC tumour |

*(Note: 11 PBMC + 5 Adjacent_normal = 16. An earlier scan had reported 9 PBMC; the deposit
contains 11 PBMC samples. The manifest reconciles exactly with the barcode file.)*

### A2. Excluded on <50 tumour-myeloid cells — 0 samples

No tumour sample failed this test. The smallest tumour-myeloid compartment was PDAC14-Tumor
with **80** cells; the next smallest was PDAC07-Tumor with 93. The ≥50 rule therefore removed
nothing in this dataset — it was applied, and it bound on no sample.

### A3. Excluded because the SIGLEC15-detectable group was empty — 5 samples

Frozen rule: the primary contrast requires **both** detectable and non-detectable myeloid cells.
These five tumour samples passed the ≥50 test but had **zero** SIGLEC15-positive myeloid cells.

| sample_id | patient | tumour-myeloid cells | detectable | non-detectable |
|---|---|---|---|---|
| PDAC05-Tumor | PA05 | 191 | 0 | 191 |
| PDAC07-Tumor | PA07 | 93 | 0 | 93 |
| PDAC08-Tumor | PA08 | 308 | 0 | 308 |
| PDAC09-Tumor | PA09 | 277 | 0 | 277 |
| PDAC14-Tumor | PA14 | 80 | 0 | 80 |

### A4. Excluded on patient overlap / duplication / integration — 0 samples

- **Patient overlap:** 11 of 14 patients contributed more than one sample to the deposit, but
  **no patient contributed more than one eligible tumour sample** (each of PA01, PA02, PA04,
  PA06, PA10, PA11, PA12, PA13 appears exactly once in the final set). The exclusion rule
  therefore bound on no sample.
- **Duplicate samples:** 0. The 29 `sample_id` values are distinct.
- **Integrated objects:** 0 exclusions on this ground — see §B4; the deposit is a concatenation
  with sample-prefixed barcodes, not an integrated object with shared barcodes.

**Arithmetic: 29 − 16 − 0 − 5 = 8.**

---

## B. What failed, and what was not resolvable as specified

### B1. No `RAW.tar` and no `filelist.txt` (informational, resolved)

The task brief described a manifest at `.../suppl/filelist.txt` listing the contents of a
`RAW.tar`. Neither exists for GSE278688: `filelist.txt` returns **HTTP 404** and
`GSE278688_RAW.tar` returns **HTTP 404**. The deposit is instead four flat supplementary files
containing **one combined matrix** for all 29 samples, not 29 per-sample matrices. This was
resolved by using the directory listing and HTTP `HEAD` requests to obtain names and byte sizes.

### B2. Reading GEO SOFT from PowerShell returned byte codes (resolved)

`Invoke-WebRequest ... $r.Content` returned a space-separated list of decimal byte values rather
than text (e.g. `94 83 69 82 73 69 83 ...`), because the response body was being coerced to a
string one byte-code at a time. Writing that string to disk produced an unreadable file on the
first attempt. Resolved by reassembling with
`[System.Text.Encoding]::UTF8.GetString([byte[]]$codes)`. Recorded in
`provenance/RUN_PROVENANCE.md` so the next run does not repeat it.

### B3. Download throughput collapsed repeatedly (resolved, cost ~35 minutes)

The 1.04 GB matrix could not be fetched at a usable sustained rate by any simple method:

| method | observed |
|---|---|
| `Invoke-WebRequest -OutFile` | buffers the whole body in memory; **no bytes on disk after 2 min** |
| single stream, fresh connection | up to 242 KB/s for ~4 MB, then decayed to ~35–50 KB/s |
| 16 concurrent connections, small ranges | 176 KB/s aggregate — **slower** than one stream |
| 8 concurrent connections, 4 MB ranges | 365 KB/s aggregate |
| **32 concurrent connections, 34 MB ranges** | **1.8 MB/s sustained** (13× single stream) |

A further failure: the resumed run stalled on its last two parts (30/32 complete, 34 MB
remaining) at ~50 KB/min, because NCBI throttles long-lived connections. Relaunching with fresh
connections cleared it. Net effect: three kill/relaunch cycles. The downloader is resumable and
validates every part length before concatenating; the final file matched the expected
**1,093,755,688** bytes exactly. One relaunch attempt was **self-terminated** by this run's own
process-cleanup filter, which matched the wrapper's command line; the filter was narrowed.

### B4. The deposit is a single combined object, not per-sample files

The brief's exclusion rule "integrated objects (avoid double counting)" was checked rather than
assumed. The matrix is 30,543 genes × 199,184 cells with 358,929,111 non-zero entries; barcodes
are sample-prefixed (`PDAC01-Tumor-GEX_out_...`) and the barcode list matches the metadata
row-for-row as a multiset with **0 duplicates**. Summing each sample's cells reproduces the
dataset total exactly (4,807 + 4,824 + … + 7,899 = 199,184). No cell is counted twice, so no
sample was excluded on this ground. **Reported as resolved, not as an exclusion.**

### B5. The deposited matrix is RAW COUNTS, and the frozen rule says "stored log1p-normalised"

**This is the one substantive specification gap, and it was resolved by reference, not by choice.**

- The brief defines detectability as *stored log1p-normalised expression > 0*. The deposit
  contains **raw integer UMI counts** (MatrixMarket flag `integer`; all 358,929,111 values
  integral; min 1, max 16,277). There is no stored log1p layer.
- **Detectability is unaffected.** Any candidate normalisation is monotone with f(0) = 0, so
  `>0` selects the identical cells either way. Verified explicitly: for SIGLEC15, ACP5 and MMP9
  the detectable sets are bit-identical under both candidates. The eligible set, the group
  sizes and the exclusion arithmetic are therefore **not** in question.
- **The programme-score magnitude IS affected**, because the z-standardisation is applied to
  normalised values in the original pipeline:

  | reconstruction | mean contrast |
  |---|---|
  | log1p(counts), no library-size scaling | **+0.388** |
  | CP10K then log1p (manuscript pipeline) | **+0.252** |

  Per-sample differences reach 0.76.
- **How it was resolved:** the manuscript's own code was located and read. It documents the
  pipeline as "Normalize to 10,000 CPM, log1p transform"
  (`biological_validation/download_data.py:370`), and the manuscript's scoring module performs
  `zscore_columns()` with **ddof = 0** on those stored values
  (`analysis_F_cross_scale.py:95-100`), with `N_RESAMPLES = 10_000` and `SEED = 15072026`
  (`analysis_F_cross_scale.py:32-34`). The published anchor was confirmed in the pipeline's own
  output: `cell_level_summary.csv` gives n = 19, mean **0.6393547**, CI 0.3275–0.9890. The
  reconstruction, the seed, the resample count and the ddof were therefore all taken from the
  existing frozen pipeline, not chosen here. Recorded in `scripts/08_final_analysis.R` header and
  in `output/primary_summary.csv`.
- **Residual uncertainty:** if the intended reading was literally `log1p(counts)` with no
  library-size scaling, the mean contrast is +0.388 rather than +0.252. Both are reported;
  the adopted one is the manuscript-consistent one.

### B6. One convention was genuinely ambiguous and was reported both ways

The "all-cell tumour aggregate" (SI §15) is defined as "every cell in the retained tumour
samples" but its aggregation is not stated (unlike the myeloid aggregate, which is inverse-log).
Both candidates were computed: plain CP10K mean and inverse-log give the **identical**
Spearman ρ = **+0.2619**, because Spearman is rank-based and the two are monotone transforms of
one another; only the bootstrap interval differs slightly. The ambiguity therefore does **not**
change any reported correlation. See `scripts/11_allcell_convention.R`.

### B7. No invariant programme gene occurred — checked, not assumed

The brief asks for any sample where a programme gene is invariant to be reported. Checked across
all 29 samples: MMP9 is invariant within sample in **6 samples, all PBMC**
(PDAC05_PBMC, PDAC07-PBMC, PDAC09-PBMC, PDAC11-PBMC, PDAC12-PBMC, PDAC14-PBMC), where no cell
has any MMP9 transcript. **Zero tumour samples and zero of the 8 eligible samples had an
invariant programme gene**, so all 8 scores were formed from **both ACP5 and MMP9** and no
sample was carried by a single gene. (Contrast with the published series, where MMP9 was
invariant in 1 of 19.) PDAC05-Tumor, which has the smallest MMP9 signal, is **not** invariant:
13 of its 191 myeloid cells have MMP9 transcripts.

### B8. Not attempted

- **No re-annotation of any cell.** The deposited `all_celltype` was used as given, as required.
- **No TREM2_TAM restricted-sensitivity analysis.** The `≥10 cells with ≥2 detectable and
  ≥5 non-detectable` rule belongs to that separate analysis and was deliberately **not** applied;
  conflating it with the primary rule was explicitly prohibited.
- **No protein inference.** Detectability is transcript detection only.
- **No other dataset was accessed.** GSE278688 only.
