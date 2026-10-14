# failed_or_excluded.md — GSE21501, Phase B step 3c-2

Two ledgers, kept apart on purpose:

* **Part A — EXCLUDED / RETAINED-BY-RULING SAMPLES.** What is not in which analysis, and on
  what authority. Nothing is silently dropped.
* **Part B — FAILURES.** Every attempt in this step that did not work, what was tried, what
  was substituted, and — where a check of mine was itself wrong — that the check, not the
  data, was what failed.

Scope reminder: this step reports **GSE21501 alone**. Nothing here is pooled with GSE205154
or any other cohort.

---

# PART A — EXCLUDED AND RETAINED-BY-RULING SAMPLES

## A.1 The head-line arithmetic

```
132  arrays deposited (GSM535948..GSM536949)
 -3  duplicate arrays, frozen rule DUP-2 (all three inside the Nebraska group)
===
129  eligible tumour arrays/specimens  ->  THE PRIMARY ANALYSIS (AD-1)
```

Of the 129, **102 are clinically annotated** (Johns Hopkins 34, NorthShore 19, Northwestern
49) and **30 are not** (Nebraska 15 + UNC 15). Under **AD-6 the 30 are flagged `clinical
annotation unavailable`, NOT excluded**: they are *in* the primary and *out* of the
clinically annotated subset, which is the prespecified sensitivity of §2 of the summary.

## A.2 The three arrays dropped by DUP-2 — the only genuine exclusions

| GSM | title | specimen key | dropped under |
|---|---|---|---|
| `GSM535950` | `12R2587.2` | `NE:12R2587` | DUP-2b (both members unmarked; lexicographically larger title) |
| `GSM536033` | `13R2778Rpt` | `NE:13R2778` | DUP-2a (the repeat-marked member) |
| `GSM536041` | `9R1701B` | `NE:9R1701` | DUP-2b (both members unmarked) |

All three carry `EX-5 duplicate array` in the frozen manifest's `exclusion_reason`, and all
three sit inside the **non**-clinically-annotated Nebraska group. This step did not choose
them: it read `EX-5` from the frozen manifest. **DUP-2 removes nothing from the 102**, so
102 is also the pre-de-duplication count for that subset.

**Note on S3.** Because the whole three-array duplicate block lies in Nebraska, the derived
S3 drop-set — repeat-marked arrays inside the **102** — is six arrays
(`GSM536488 P271N Replicate`, `GSM536894 P164T redo`, `GSM536897 P180T Replicate`,
`GSM536922 P275T Replicate`, `GSM536923 P31T Replicate`, `GSM536926 P41T Replicate`),
giving n = 96 as the frozen document states. This was derived by matching
`rpt|redo|replicate` against the manifest titles, not copied from the frozen table, and it
reproduced the frozen set exactly.

## A.3 `GSM536901` / `P188T IPMN` — excluded from S2, on the AD-4 ground

The array is excluded from S2 **because it is explicitly labelled IPMN and therefore is not
PDAC** — a precursor lesion rather than invasive PDAC.

It is **not** excluded for missing information, and it does not have missing information:
the frozen manifest records `tumour_or_non_tumour = tumour`, `primary_or_metastatic =
primary`, `eligible_for_main_meta = yes`, an **empty** `exclusion_reason`, and a full
clinical record (`os time 6`, `os event 1`, `t stage 1`, `n stage 1`, High Risk). The IPMN
label lives in the sample **title**, and AD-4 rules that the title label governs. Excluding
it for absent annotation would be factually wrong about this array.

**It remains in the primary (n = 129).** AD-4 governs S2 only.

## A.4 `GSM536488` / `P271N Replicate` — RETAINED, and this departs from the step-3b arithmetic

AD-4 permits excluding `P271N` **only if the frozen manifest established, from its original
sample identity, that it is non-PDAC or normal**, and forbids excluding it merely because of
the trailing `N`.

**The frozen manifest does not establish that.** For `GSM536488` it records:

| field | value |
|---|---|
| `tumour_or_non_tumour` | `tumour` |
| `primary_or_metastatic` | `primary` |
| `channel2_sample_title` | `P271N Replicate` |
| `tissue_type` (channel 2) | `Patient Primary Pancreatic Tumor` |
| `clinically_annotated` | `yes` |
| `eligible_for_main_meta` | `yes` |
| `exclusion_reason` | *(empty)* |

The manifest's own test — `established_non_pdac_or_normal = (tumour_or_non_tumour !=
"tumour") or (primary_or_metastatic == "normal")` — evaluates to **False**. The array is
therefore **retained**, including in S2.

**Consequence, stated plainly because it changes a frozen number.** The step-3b text gave
S2 n = 100 (102 − 2). The ruling removes one of those two grounds, so the ruled S2 is
**102 − 1 = 101**. Both are computed and both are reported:

| id | n | ρ | status |
|---|---|---|---|
| **S2** (AD-4 as ruled) | **101** | +0.1703 | the ruled version |
| S2-legacy (step-3b text) | 100 | +0.1462 | retained for comparison only |

This is a ruling-driven arithmetic change, so it is surfaced rather than absorbed. The `n`
column in `output/sensitivity_results.csv` is **read from the data**; it is never written
down from the frozen document, which is why the departure shows up as a value rather than as
a silent inconsistency.

## A.5 Everything else is retained

* **`GSM536940 PE621T`** (no N stage) and **`GSM536892 P152T`**, **`GSM536924 P35T`**,
  **`GSM536927 P42T`** (no T stage) — retained under AD-5. No clinical covariate enters this
  estimand, so an incomplete clinical record is not a ground for exclusion.
* **The five arrays the 3b record called "excluded for nothing"** — retained identically.
* **No array was excluded, dropped, re-weighted, or re-shaped because of the observed ρ
  magnitude or direction.** The exclusion rule was applied before ρ was computed, and ρ was
  not recomputed under any alternative cohort that the frozen rules do not prespecify.

---

# PART B — FAILURES, WITH WHAT WAS TRIED

## B.1 A verification check of mine was wrong and halted the run (C7, first form)

**What I tried.** A self-consistency check requiring the primary bootstrap interval to be
**bit-identical** under two different RNG-addressing schemes.

**What happened.** `FATAL: the primary estimate is not stable across the two RNG addressing
schemes` — the script exited 1 rather than emit a number it could not stand behind.

**Why it failed, and why the check was wrong rather than the data.** Two different random
seeds cannot produce identical resampled intervals; the assertion was **unsatisfiable by
construction**. The permutation *p* is also a resampled quantity and is legitimately
generator-dependent, which I had wrongly assumed it was not.

**What was corrected.** The invariant that the frozen rule actually implies: **ρ is exactly
invariant to how the generator is addressed, because no generator enters ρ.** That holds
exactly (identical to the last digit). The resampled quantities are compared against a
Monte-Carlo tolerance derived from the resample count (CI gap 0.0032; permutation-*p* gap
0.0052 against one standard error of 0.0022). The correction is recorded here and in
`output/run_log.txt` §7 rather than applied silently.

**What this cost.** One aborted run. No number changed: the primary ρ, its interval and its
seeds are the same before and after, because the check never touched the estimate.

## B.2 A structural control caught a real parser defect (C12, first form)

**What I tried.** A row-accounting control asserting that the platform table and the matrix
reconcile: symbol-bearing platform rows + symbol-less platform rows = matrix features, and
the matrix IDs the annotation cannot name are exactly the symbol-less rows.

**What happened.** `FAIL: 45221 - 12226 = 32995; matrix features = 45220; matrix-only IDs =
20`. The numbers disagreed with each other, so at least one was wrong.

**What was actually wrong — a defect in my own parser, plus a second defect in my own
diagnostic of it.**

1. **The parser (real defect).** The platform SOFT file ends with a
   **`!platform_table_end`** terminator. My first parser screened for comment lines (`#`,
   `!`, `^`) only *before* it had seen the `ID` header, so once the table had started the
   terminator was counted as a **data row**: 45 221 rows for a 45 220-row table. It failed
   the symbol-column length test, so it was filed under "empty symbol" and inflated that
   tally too — 12 226 where the true value is 12 225. Both off-by-one errors landed in the
   same subtraction, and `45221 - 12226 = 32995 = 45220 - 12225` was the coincidence that
   made the middle number look right while the outer two were wrong.

2. **My first diagnostic script (second defect, and the one that actually misled me).** Its
   job was to find which matrix IDs the annotation could not name. It had *its own* parsing
   difference and so reported "matrix-only: 0, platform-only: 1, that one being
   `!platform_table_end`" — i.e. it claimed there was **no** discrepancy at all. That was
   wrong, and because it disagreed with the control I had just seen fail, I went looking for
   a third explanation.

3. **What that misreading cost.** The honest third fact was sitting in plain sight and I
   misread it as a discrepancy when it was a **display truncation**. The QC field
   `matrix_ids_missing_from_platform` was written as `unnamed[:20]` — a *sample* of the
   first twenty entries — and I compared its **length (20)** against the empty-symbol tally
   (12 226) and concluded the two sets were unrelated. They were not unrelated. The real
   count was 12 225 and those twenty IDs were simply the first twenty members of it.

**What is true now, established by re-running the original parser logic verbatim and
comparing row by row** (`work/_diag_orig_parser.py`, retained):

```
matrix rows                                 45220
OLD parser: raw rows counted                45221   <- the terminator
OLD parser: rows failing the symbol test        1   <- ['!platform_table_end']
OLD parser: rows with an empty symbol cell  12225
OLD parser: rows carrying a gene symbol     32995
OLD parser: matrix IDs not in the symbol map 12225
of the first 20 of those, empty-symbol rows     20
annotation rows too short                       0
TRUE matrix IDs not nameable                12225
```

So the twenty IDs were **correctly the empty-symbol rows from the beginning**; they were
only ever truncated in the report. The genuinely broken quantity was the raw row count, and
that is exactly what the control was built to catch.

**What was corrected.** The parser now classifies every row explicitly and reports the
ledger instead of a pair of totals, refuses to continue unless the accounting closes, and
writes the **full** unnamed ID list with its count rather than a truncated sample. Stage 1
now prints:

```
45220 platform rows = 32995 with a gene symbol + 12225 without a usable symbol
32995 + 12225       = 45220 matrix features
symbol-bearing platform rows absent from the matrix = 0
```

C12 was restated against those identities and passes. **It was not relaxed to make it pass**
— the defect was fixed and the check re-derived from the corrected ledger.

**What this cost.** Two diagnostic scripts (one obsolete, both retained under `work/`), two
extractor revisions, and re-runs of all four stages. **The target genes were never at risk.**
SIGLEC15 (31474), ACP5 (10450) and all ten MMP9 features carry valid symbols, were present
in the matrix under every version of the parser, and matched the frozen probe mapping in
every run. ρ is identical to the last digit across the defective and corrected versions. The
defect was entirely in the row *bookkeeping around* the target features.

## B.3 A control of mine was mis-toleranced (C8, first form)

**What I tried.** Requiring a 2000-permutation empirical null *p* to agree with the reported
10 000-draw permutation *p* to within 0.05.

**What happened.** `PASS` — but for the wrong reason, and it was flagged as marginal
(0.0450 vs 0.0527, gap 0.0077 against a 0.05 threshold).

**Why it was wrong.** The two quantities are **independent Monte-Carlo estimates of the same
*p***, so they are expected to differ by roughly a standard error, not to agree exactly. At
*p* ≈ 0.05 with 2000 draws, one standard error is 0.0049; a 0.05 tolerance is ten of them
and would have passed almost any disagreement, including a real defect. The tolerance was
nonsense in both directions at once.

**What was corrected.** The comparison is now against **5 standard errors** computed from
the observed *p* (threshold 0.0232 for the observed gap of 0.0077), and the *null
distribution's own* moments are checked separately against theory (mean −0.0057 against 0,
sd 0.0862 against 1/√(n−1) = 0.0884). Both parts pass, and now they mean something.

## B.4 `pwsh` is not available

`pwsh` is not on `PATH` on this host. All console work was done with Windows PowerShell 5.1
via `powershell.exe`. No logger in this step shells out to a shell at all; every stage is a
plain `python <script>` invocation, so no stage depends on `pwsh` existing.

## B.5 Things deliberately NOT attempted

* **No HTTP request was made by any stage.** The deposit (`GSE21501_series_matrix.txt.gz`),
  the platform annotation (`GPL4133.annot.gz`) and the manuscript module were all already
  local, and each was hash-checked against the value recorded in the step-3b `sessionInfo.txt`.
  The Europe PMC `fullTextXML` 406 and `pmc.ncbi.nlm.nih.gov` HTML workaround noted in the
  host brief were therefore not needed here, and no workaround was substituted for them.
* **`GSE21501_family.soft.gz` (187 687 633 bytes) was never obtained** — the step-3b transfer
  stalled, and it is recorded as unobtained there. Nothing in this step depends on it: the
  gene-identity rule was applied directly to `GPL4133.annot.gz` plus the series matrix, which
  is the same pair of inputs the frozen probe mapping was built from. This remains an
  **unobtained** file, not a failure of this step.
* **No pooling with GSE205154.** Out of scope by instruction and not performed; GSE205154 is
  handled by a separate agent.
* **No manuscript text, figure or claim was written or modified** by this step. Outputs are
  confined to `step3c_rho/`.

## B.6 Not a failure, but recorded because it looks like one

**S4 returns exactly the primary ρ (+0.171188), which invites the suspicion that averaging
was a no-op.** It is not a no-op; it is arithmetic. The three duplicate groups are
**byte-identical technical replicates** — the within-group range is **0.0000 for all three
genes in all three groups** — so every retained array already lies exactly at its group mean
and averaging changes no value. Control C11 records the per-group means and ranges so the
identity can be checked rather than taken on trust.

This is itself a finding about the deposit worth one line: GSE21501's duplicate pairs are
exact technical replicates, which is why the frozen choice between *selection* (DUP-2) and
*aggregation* (S4) cannot move this cohort's estimate at all. It could move a cohort whose
replicates disagree, so the frozen rule's outcome-blindness still matters in general.
