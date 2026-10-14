# results_summary.md — Step 3c‑1: GSE205154 bulk‑extension correlation

**Scope: GSE205154 alone.** Nothing in this directory is pooled with GSE21501 or with any
other cohort. No random‑effects combination was performed. No manuscript text, figure or
claim was written or modified.

---

## 0. Integrity of the frozen cohort

| item | value |
|---|---|
| `sample_manifest.csv` SHA‑256 recorded in the freeze | `f76b81eceb0557cf4c1fee0138e9137a84508369f47ad163255a096af78c708b` |
| SHA‑256 recomputed at run time | `f76b81eceb0557cf4c1fee0138e9137a84508369f47ad163255a096af78c708b` |
| **verdict** | **MATCH — the manifest was not edited after the freeze** |

The three deposited matrices were also re‑hashed against §1.1 of the frozen definition:
all three match. Assertion tally across the four logs, all passing:
**29** frozen‑cohort guards (`output/manifest_audit.txt`),
**40** provenance guards (`output/expression_extraction_log.txt`),
**24** analysis and control guards (`output/correlation_log.txt`),
**9** independent R cross‑checks (`output/crosscheck_R_log.txt`) — **102 assertions, 0 failures**.

**Platform disclosure.** GPL20301 is **exome‑capture** bulk RNA‑seq (FFPE sections, xGEN
Exome Research Panel v1.0 hybrid capture), **not whole‑transcriptome**. It is admissible
under the frozen rules because only **within‑cohort** ranks are used and no cross‑platform
absolute comparability is claimed. This disclosure is carried in every result file.

---

## 1. Primary estimate

| | |
|---|---|
| cohort | GSE205154 |
| compartment | primary tumours |
| **unit** | **patient** (AD‑1) |
| **n** | **216** patients / 216 specimens |
| representation | `log2(TMM+1)`, then the frozen mean‑z programme score (AD‑2) |
| score | mean of the across‑sample z‑scores of ACP5 and MMP9, `ddof = 0` |
| target | SIGLEC15 (enters through ranks) |
| **ρ** | **+0.1317** |
| **95% CI** | **[−0.0024, +0.2648]** |
| **CI method** | **nonparametric bootstrap percentile interval, 10 000 resamples, seed 15072026** — the manuscript's own `spearman_bootstrap()`, not BCa |
| CI, secondary method | asymptotic Fisher z, SE = 1/√(n−3): [−0.0018, +0.2606] |
| asymptotic p | 0.0533 |
| permutation p | 0.0532 (10 000 permutations, seed 15072026) |

Both interval methods place the interval **across zero**: the bootstrap lower limit is
−0.0024 and the Fisher lower limit is −0.0018. The point estimate is positive and small.

**One choice this step had to make, disclosed rather than buried.** The frozen text fixes
the combination machinery (Fisher z, REML τ², modified Hartung–Knapp, seed 15072026,
10 000 resamples) but does not name a *per‑cohort* interval method. It was resolved by
reusing the manuscript's **own** `spearman_bootstrap()` — the identical function that
produced the existing TCGA‑PAAD, GSE57495 and GSE85916 per‑cohort intervals — so this new
cohort's interval comes off the same machinery as the three already reported. The
asymptotic Fisher‑z interval is reported alongside it because Fisher z is the scale the
frozen combination step operates on. Bootstrap percentile endpoints are Monte‑Carlo
estimates (the R cross‑check, using a different generator with the same seed, reproduces
them to ≤ 0.016); the point estimate, which is not Monte‑Carlo, is reproduced exactly.

---

## 2. Sensitivity spread

All five prespecified sensitivities are in `output/sensitivity_results.csv`; each row is
labelled specimen‑level or patient‑level.

| # | sensitivity | level | n | ρ | 95% CI (bootstrap) |
|---|---|---|---|---|---|
| S1 | `log2(TPM+1)` (AD‑2 i) | patient | 216 | **+0.1547** | [+0.0175, +0.2904] |
| S2 | unlogged TMM (AD‑2 ii) | patient | 216 | **+0.1456** | [+0.0103, +0.2776] |
| S3 | duplicate pair averaged (D4.5 / AD‑3) | patient | 216 | **+0.1352** | [+0.0012, +0.2694] |
| S4 | metastatic specimens (AD‑4) | patient | **71** | **−0.0351** | [−0.2664, +0.2002] |
| S5 | metastatic patients, overlap‑free (AD‑4) | patient | **61** | **+0.0230** | [−0.2327, +0.2849] |
| — | **primary** | patient | 216 | **+0.1317** | [−0.0024, +0.2648] |

Spread within the **primary compartment**: ρ ∈ [+0.1317, +0.1547], a range of 0.023.
Spread in the **metastatic compartment**: ρ ∈ [−0.0351, +0.0230], a range of 0.058, both
essentially null and both straddling zero.

The 71 vs 61 distinction (AD‑4): the full metastatic set is **71 specimens from 71
patients**; 10 of those patients also contribute a primary specimen to the main analysis,
so the overlap‑free subset is **61 patient‑level units**. Neither metastatic estimate is
distinguishable from zero, and the two differ by 0.058.

`output/extra_variants.csv` additionally reports the AD‑1 variant the author ruled against
(218 primary specimens as 218 independent units): ρ = +0.1321, n = 218. On this data that
choice moves the point estimate by 0.0004. **AD‑1 is nonetheless applied**, because it
rules the specimen‑as‑unit reading out on patient‑identity grounds, not on the effect size.
The rejection‑of‑ρ‑based‑exclusion rule is respected: the frozen cohort was used as frozen
regardless of the magnitude or direction of ρ.

---

## 3. What is stable, and what is not

**Stable**

* The **sign and rough magnitude** of the primary‑compartment association. All three
  primary‑compartment variants (S1, S2, S3) and the primary itself lie in
  ρ ∈ [+0.13, +0.16] with the same sign.
* The **conclusion that the primary‑compartment interval crosses zero** is stable across
  both interval methods and across every primary‑compartment representation.
* The **metastatic compartment is null** under both n = 71 and n = 61, and under the
  primary representation.
* **SIGLEC15's side of the correlation is invariant to the log transform of a given
  matrix.** Asserted directly: SIGLEC15 ranks are identical before and after `log2(x+1)`
  for TMM, TPM and counts (all non‑negative, so the transform is strictly monotone).
  This is exactly the scope AD‑2 claims, and it holds.
* **The result is reproducible in two independent engines.** Base‑R code that
  re‑implements the frozen rule from the frozen text, and re‑derives the sample sets from
  the manifest rather than reading the Python pipeline's output, reproduces **every** point
  estimate to < 1e‑9 (mostly < 1e‑16). See `output/crosscheck_R_log.txt`.

**Not stable / not established**

* **Whether the primary association is non‑null.** It is not: both intervals include zero,
  and p = 0.053 is above 0.05 by a margin small enough that it should not be leaned on.
  Nothing here supports a claim of a positive association in this cohort.
* **The exact value of ρ is representation‑dependent to a modest degree.** Moving from
  `log2(TMM+1)` to `log2(TPM+1)` moves ρ by +0.023, which is about the same size as the
  distance from the interval to zero. The sensitivity analysis was prespecified precisely
  because AD‑2 flagged this, and it shows the representation choice is not inert.
* **Whether TMM and TPM are interchangeable for SIGLEC15.** They are **not**. Only 14 of
  216 specimens keep the same SIGLEC15 rank under both representations (202 differ; mean
  absolute rank shift 10.4 of 216; max 52), even though Spearman(log2(TMM+1),
  log2(TPM+1)) = +0.9707 across all 289 specimens and **neither matrix contains a single
  exact tie** among the 289 SIGLEC15 values. TMM and TPM are simply not monotone in one
  another. Consequence: the frozen statement that the transform moves only the
  ACP5/MMP9 side is true for a monotone transform **of one matrix**, but it does **not**
  extend to swapping one supplied matrix for another. The control that was supposed to
  demonstrate invariance instead returned **NOT APPLICABLE**, and it is reported that way
  rather than massaged into a pass.
* **The metastatic estimates should not be read as "no association".** At n = 71 and
  n = 61 the intervals are wide (±0.25), so they are compatible with effects as large as
  the primary compartment's and larger. They are uninformative about direction, not
  evidence of absence.

---

## 4. Verification — controls capable of failing

Verification is not self‑confirming: the scoring rule was executed by importing the
manuscript's own functions verbatim from the v39 archive
(`analysis_F_cross_scale.py`, sha256 `7b6a311a…`, with per‑function source hashes recorded
in `output/correlation_log.txt`), **and** independently re‑implemented in base R, **and**
run against conventions that are supposed to fail. Full table:
`output/controls_capable_of_failing.csv`.

| control | convention | ρ | outcome |
|---|---|---|---|
| **C‑A2** | z‑standardise **across genes** instead of across samples (wrong axis) | **−0.2016** | **materially different** (Δ = 0.333, sign flips). Score sd falls from 0.826 to 0.061. **Gate PASSES** — the axis convention is doing real work. |
| **C‑A** | wrong axis with only the two programme genes | **NA** | Score collapses to a numerical constant (spread 1.3e‑13) because per‑sample z forces ACP5_z = −MMP9_z. Spearman on a constant is undefined, so it is reported as **NA**, not as a spurious near‑zero number. |
| **C‑C** | raw counts, no log transform | **+0.4274** | materially different (Δ = 0.296) — exactly the outlier‑dominated scale AD‑2 warns against. |
| **C‑C2** | `log2(counts+1)` | **+0.3873** | Δ = 0.256. Informational: this cohort's counts do not reproduce the TMM/TPM answer, which is itself a caution about representation. |
| **C‑B** | SIGLEC15 from `log2(TPM+1)`, score from `log2(TMM+1)` | **NOT APPLICABLE** | The invariance claim was tested and **failed its precondition** (202/216 rank disagreements). Reported as not applicable, with the rank counts, rather than as a pass. |
| **C‑D** | 218 specimens treated as 218 patients (AD‑1 rejected variant) | +0.1321, n = 218 | Informational. Δ = 0.0004. |

The C‑A2 gate is the one that could have failed and did not; C‑B is the one that was
expected to confirm an invariance claim and instead came back inconclusive. Both outcomes
are reported as obtained.

---

## 5. What could not be determined

1. **Whether the `F` title suffix means anything biological.** All 289 samples were
   classified from the `tumor type` field, never from the filename suffix (AD‑5); the
   seven `F` samples are primary because `tumor type = Primary`, and the deposit never
   defines `F`. No decision depends on it.
2. **Whether two differently numbered specimen cores belong to one patient**, for the 83
   cores with no GSE281005 counterpart (87 samples). If any such split exists, n = 216 is
   an overcount. Public metadata cannot settle it.
3. **Patient‑level overlap with GSE282302 and with TCGA‑PAAD.** Not assessable from public
   metadata. Does not affect this cohort's arithmetic, but does bear on any later
   combination.
4. **Whether the C‑B discrepancy would survive with a different representation choice.**
   The 202/216 rank disagreement between TMM and TPM for SIGLEC15 is a property of this
   deposit; no third representation exists in the deposit to break the tie.
5. **Any pooled or meta‑analytic estimate.** Deliberately not computed here. Combining
   this cohort's ρ with others requires the frozen REML / modified Hartung–Knapp
   machinery, and is outside this step.
6. **Raw sequence data.** Controlled access via dbGaP `phs003597.v1.p1`; not obtained and
   not needed for the frozen protocol.

---

## 6. Files

Every result file is paired with the script and seed that produced it; see
`output/sessionInfo.txt` (block 3) and `output/run_manifest.json`. The two headline result
tables sit at the step root; every verification and intermediate artefact is under
`output/`.

| file | contents |
|---|---|
| `primary_result.csv` | the cohort‑level estimate, ρ, n, CI + named method, resamples, seed |
| `sensitivity_results.csv` | one row per prespecified sensitivity, level labelled |
| `results_summary.md` | this file |
| `failed_or_excluded.md` | everything that failed, with what was tried |
| `run_log.txt` | append‑only log of every command run, including superseded failing runs |
| `output/extra_variants.csv` | AD‑1 / D4.5 variants, explicitly not sensitivities |
| `output/controls_capable_of_failing.csv` | the control battery |
| `output/manifest_audit.txt` | 29 frozen‑cohort guards |
| `output/expression_extraction_log.txt` | 40 provenance guards on the three matrices |
| `output/correlation_log.txt` | the primary run, verbatim, with per‑function source hashes |
| `output/crosscheck_R_log.txt` | the independent R re‑derivation and its agreement table |
| `output/sessionInfo.txt` | Python and R session records; input and output hashes |
| `output/sessionInfo_R.txt` | the live R `sessionInfo()` output |
| `output/manifest_verification.csv` | all 289 samples with their set membership |
| `output/expression_3genes.csv` | the three extracted genes × 289 samples, 3 representations |
| `output/sample_sets.json` | the frozen sample sets, with the manifest hash |
| `output/run_manifest.json` | machine‑readable script → output → seed mapping |

Reproduce from scratch, in this order:

```
python  scripts/01_verify_manifest.py
python  scripts/02_extract_expression.py
python  scripts/03_correlation.py
Rscript scripts/04_crosscheck_R.R
python  scripts/05_session_info.py
```
