# failed_or_excluded.md — Step 3c‑1, GSE205154

Everything that failed, was excluded, or could not be determined, with what was tried.

Two distinct categories are kept separate, because they have different consequences:

* **Part A — specimens excluded from the analysis.** These are *by design*, frozen in
  advance, and are the reason n = 216 rather than 218.
* **Part B — things that went wrong during this step.** Script bugs, a mis‑specified
  guard, controls that came back inconclusive. None of these touched a reported number;
  each was caught by a check and is recorded with what was tried.

---

## Part A — excluded specimens (by design, frozen)

### A1. Two specimens excluded from every analysis

| GSM | title | patient | reason | recorded in |
|---|---|---|---|---|
| `GSM6206193` | `ST-00017945-T2` | ST‑00017945 | second primary specimen from a patient already represented; frozen D4.2 one‑specimen‑per‑patient rule retained `ST-00017945-T` (`GSM6206194`) | `sample_manifest.csv` `exclusion_reason`, reproduced verbatim in `output/manifest_audit.txt` |
| `GSM6206245` | `ST-00019598-T2` | ST‑00019598 | same, retained `ST-00019598-T` (`GSM6206244`) | as above |

Both specimens remain present in `output/expression_3genes.csv` (all 289 columns are
extracted) and enter no analysis. They are used **only** in the D4.5 pair‑averaging
sensitivity S3 and in the D4.5 alternative X2, both explicitly labelled as such.

Verified independently in two engines: the Python guard G6 and its R equivalent both
confirm the base‑label specimen is the one flagged into the main set for each core.

### A2. A replicate‑looking title that is NOT a duplicate

`GSM6206205` / `ST-00020168-M2` carries an `M2` suffix but its core contains **no** `-M`
specimen, so it is an `M2`‑labelled metastasis from a distinct patient, not a
within‑compartment duplicate. The metastatic compartment therefore has **zero** duplicates
(71 specimens in 71 cores). Checked explicitly (guard G5c) rather than assumed from the
suffix, per the AD‑5 principle of never inferring identity from a filename suffix.

### A3. Nothing was excluded on the basis of ρ

No specimen, patient, compartment or representation was dropped, reweighted or reshaped
because of the magnitude or direction of the observed correlation. The cohort was fixed
before ρ was computed, and the frozen manifest was used unchanged (SHA‑256 verified at run
time). The metastatic compartment returned ρ ≈ 0 and remains in the report at both n = 71
and n = 61.

### A4. Eligibility conditions that excluded nothing

Recorded for completeness, because a condition that excludes nothing is still a condition
that was applied:

| condition | samples excluded |
|---|---|
| non‑tumour tissue | **0** — all 289 samples are pathologist‑marked tumour |
| missing `tumor type` | **0** — populated for 289/289 |
| gene unusable (SIGLEC15 / ACP5 / MMP9) | **0** — all three `protein_coding`, one row each, `TMM_pass_filter = TRUE` |
| metadata conflict between the three compartment fields | **0** of 289 |

---

## Part B — failures during this step

### B1. A mis‑specified guard in `01_verify_manifest.py` (fixed)

**What failed.** Two guards were written from the frozen narrative rather than from the
manifest itself, and both fired on the first run:

* G5 asserted that the 216 main specimens sit in **214** distinct patient cores.
* G6 looked for the two duplicate cores inside the **216‑row** main set.
* G12 asserted the platform string contains the lowercase phrase `exome capture`.

**What was tried / found.** The manifest itself was inspected directly. The 216 main rows
already have the D4.2 exclusion **baked into their eligibility flags**: the 218 primary
*rows* sit in 216 distinct patients, and subtracting the two `-T2` rows leaves the 216
`eligible_for_main_meta = yes` rows mapping 1:1 onto 216 patients. G6's stated intent —
"does the manifest's flag agree with the frozen base‑label rule?" — therefore had to be
checked against the **218 primary specimens**, not the 216 eligible ones. G12 failed on
letter case only (`xGEN Exome Research Panel` vs the lowercase literal I searched for).

**Resolution.** The guards were corrected to the values the frozen document actually
fixes (216 main specimens → 216 distinct patients; 218 primary specimens → 216 patients of
which 2 contribute two; duplicate cores checked within the 218) and **given a stronger
form**: the main set is now asserted to contain no patient twice at all, and the D4.2
base‑label retention is asserted to agree with the manifest flag for both cores. All
guards pass. The superseded first‑run output is preserved in `run_log.txt` with its
`[FAIL]` lines, so the correction is auditable rather than silent.

### B2. Wrong axis convention for `zscore_columns` in the C‑A control (fixed)

**What failed.** The first version of the C‑A control crashed:
`ValueError: the 'keepdims' parameter is not supported in the pandas implementation of mean()`.

**What was tried.** The control was rewritten on raw NumPy arrays rather than
`pandas.DataFrame.mean()`. While fixing it, a second, more interesting problem surfaced:
with only **two** programme genes, per‑sample standardisation forces `ACP5_z = −MMP9_z` in
every sample, so the mean‑z score collapses to a numerical constant (spread 1.3e‑13) and
Spearman on it is **undefined**.

**Resolution.** (i) C‑A is reported as **ρ = NA**, not as the small spurious value that
floating‑point round‑off would have produced (an earlier draft printed `+0.014754` from
exactly that artefact — it was caught and removed, not explained away). (ii) A
**non‑degenerate** wrong‑axis control, **C‑A2** (standardise across all three genes of the
frame), was added; that is the control the gate actually tests, and it returns
ρ = −0.2016 vs the primary +0.1317, Δ = 0.333. (iii) A spread assertion confirms C‑A2 is
testing a real alternative convention rather than an empty vector.

### B3. The invariance control C‑B came back NOT APPLICABLE

**What was tried.** The frozen text states that a monotone transform can move only the
ACP5/MMP9 score side, because SIGLEC15 enters through ranks. The intended control was to
take SIGLEC15 from `log2(TPM+1)` while the score stayed on `log2(TMM+1)`, and confirm the
result was bit‑identical.

**What was found.** The control **failed its precondition**, and is reported as
**NOT APPLICABLE** rather than converted into a pass:

* SIGLEC15 ranks under `log2(TMM+1)` and `log2(TPM+1)` differ for **202 of 216** main
  specimens (only 14 keep the same rank).
* Mean absolute rank shift 10.4 of 216; maximum 52.
* Spearman between the two SIGLEC15 vectors across all 289 specimens is only **+0.9707**.
* **Neither matrix contains a single exact tie** among the 289 SIGLEC15 values, so this is
  not a tie‑breaking artefact: TMM and TPM are simply not monotone in one another.

**What this does and does not mean.** The frozen claim is about a monotone transform *of a
given matrix*, and that claim was tested separately and **passes**: SIGLEC15 ranks are
bit‑identical before and after `log2(x+1)` for TMM, TPM and counts, and all values are
non‑negative so the transform is strictly monotone. The claim does **not** extend to
swapping one supplied matrix for another, and the control demonstrates that it does not.
This is a genuine finding about representation sensitivity and is carried into
`results_summary.md` §3.

### B4. R cross‑check errors (all fixed)

The independent R script failed three times before running clean. Each failure is recorded
because each was a real defect, not a formatting nuisance:

| # | error | cause | fix |
|---|---|---|---|
| 1 | `Error in z_by_sample[, c("ACP5","MMP9")] : subscript out of bounds` | base R's `apply()` drops column names | restore `dimnames()` after `apply()` |
| 2 | `Error in parse_con(...) : invalid char in json text … "rho": NaN` | Python's `json.dumps` emits bare `NaN`, which is not valid JSON | the Python handoff now serialises non‑finite values as `null`, with an assertion that the emitted text contains neither `NaN` nor `Infinity` |
| 3 | `$ operator is invalid for atomic vectors` | `jsonlite::fromJSON` simplified the array of sensitivities into a data frame | `simplifyVector = FALSE` |

### B5. Two reporting defects caught in review (fixed)

* `extra_variants.csv` listed "218 specimens as independent units" and "both duplicate
  members retained" as two variants. They are the **same analytic set** reached from two
  different frozen rules, so their estimates are identical *by construction*. This is now
  stated explicitly, and an assertion confirms the two agree.
* `run_manifest.json` initially listed itself among the files it describes. The
  self‑reference was removed.
* The two headline result tables were first written to `output/` alongside the
  intermediates, but the step's output list names them at the step root. They were moved to
  the step root, the scripts were updated to write them there, and the run was repeated.
  Byte‑identical copies left behind in `output/` were removed after confirming the hashes
  matched, so there is exactly one copy of each result.

---

## Part C — what could not be determined

Carried in full in `results_summary.md` §5. Summarised here so that nothing is only in one
place:

1. **Meaning of the `F` title suffix** (7 samples). Undocumented in the deposit. All three
   metadata fields agree they are primary; the decision rests on `tumor type`, not on the
   suffix (AD‑5), so nothing depends on resolving it.
2. **Whether two differently numbered specimen cores belong to one patient**, for the 83
   cores with no GSE281005 counterpart (87 samples). If core‑splitting exists there, n = 216
   is an overcount. Not assessable from public metadata.
3. **Patient‑level overlap with GSE282302 and with TCGA‑PAAD.** Not assessable; no shared
   namespace. Does not affect this cohort's arithmetic.
4. **Any pooled or meta‑analytic estimate.** Deliberately not computed. Out of scope for
   this step.
5. **Raw sequence data.** Controlled access (dbGaP `phs003597.v1.p1`); not obtained and not
   needed.
6. **Sex and age for 87 samples.** Available only through GSE281005; neither enters any
   frozen eligibility rule.
7. **Whether a third representation would resolve the TMM/TPM divergence.** The deposit
   supplies no further representation, so the C‑B question cannot be closed with this data.

---

## Part D — limitations of the verification itself

Stated so the verification is not read as stronger than it is:

* The **bootstrap percentile endpoints cannot be reproduced across engines** exactly.
  Python's `default_rng(15072026)` (PCG64) and R's Mersenne‑Twister draw different
  resamples by design. Agreement is therefore bounded by a Monte‑Carlo tolerance; observed
  maximum disagreement across all six estimates was **0.016** (metastatic n = 71), and
  0.0027 on the primary. **Point estimates are not Monte‑Carlo and agree to < 1e‑9.**
* The R cross‑check is independent in *implementation* — it re‑derives the sample sets
  from the manifest and re‑implements the frozen scoring rule in base R — but it is **not
  independent in data**. Both engines read the same three extracted gene rows from the same
  deposited matrices. A defect in the deposited matrices themselves, or in the
  extraction script, would survive both.
* The extraction script is guarded by 40 provenance assertions (file hashes, row counts,
  gene identity, row indices, column identity and order, value scale), but those checks
  verify the extracted rows against the frozen record, not against an orthogonal source. No
  independent re‑quantification from raw reads was possible: raw data are controlled‑access.
* Across the four logs there are 102 assertions and 0 failures, but **an assertion is only
  as good as its specification.** The two specification errors found in this step (B1) were
  both cases where a guard passed its own intent while testing the wrong object; they were
  caught by a human reading the first‑run failures, not by the guards themselves.

---

## Part E — nothing was excluded for convenience

No sample, gene, representation, compartment or analysis was dropped, replaced or
re‑specified at any point in this step in response to a result. The five prespecified
sensitivities are reported **all five**, whatever they showed; the one control that was
expected to pass did not, and is reported as not applicable; and the one estimate that
could have been quietly dropped for being null (the metastatic compartment) is reported at
both frozen n.

**Choices this step had to make, where the frozen text did not settle them.** Each is
disclosed and given an explicit reason rather than resolved silently:

| choice | resolution | reason |
|---|---|---|
| which per‑cohort CI method | the manuscript's own `spearman_bootstrap()` (percentile, not BCa) | it is the identical function that produced the existing TCGA‑PAAD / GSE57495 / GSE85916 per‑cohort intervals, so this cohort's interval comes off the same machinery |
| RNG state for the bootstrap | the module RNG is reset to seed 15072026 immediately before each bootstrap call | makes every interval exactly reproducible from (seed, n) alone, and matches the manuscript's behaviour for the first cohort/program pair |
| whether to pool | **not pooled**; no random‑effects combination performed | the task fixes this cohort as the whole scope |
| whether to report an interval for the null compartment | reported, at both n = 71 and n = 61 | AD‑4 requires both; the frozen exclusion rule forbids dropping a cohort for a ρ of any magnitude or direction |

Nothing was escalated as a blocking question, because each of the above is either fixed by
the frozen text or resolved by reusing the existing analysis identity unchanged. The one
point genuinely not resolvable from the public deposit — patient identity for the 83
unmatched cores, and the meaning of the `F` suffix — is recorded in Part C as
indeterminate rather than decided by assumption.
