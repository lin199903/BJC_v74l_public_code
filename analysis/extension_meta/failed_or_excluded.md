# failed_or_excluded.md — secondary bulk-extension meta-analysis

Everything that failed, was excluded, or could not be determined, with what was tried.

Two ledgers, kept apart on purpose:

* **Part A — excluded / not analysed by design.** What is not in which estimate, and on what
  authority. Nothing is silently dropped from a cohort that was eligible.
* **Part B — failures during this step.** Every attempt that did not work, every check of mine
  that was itself wrong, and — separately — **Part C: controls that came back NULL**. A control
  that fails to detect a wrong convention is a result about the data, not a licence to relabel
  it as a pass.

---

# PART A — WHAT IS NOT IN THE ESTIMATES

## A.1 Nothing was excluded on the basis of ρ

No cohort, specimen, patient or representation was dropped, reweighted or reshaped because of
the magnitude or direction of its correlation.

* All **three original PDAC bulk cohorts** enter every model (TCGA_PAAD, GSE57495, GSE85916).
* **GSE205154** enters the four-cohort estimate at its frozen n = 216, ρ = +0.1317 — a positive
  but zero-crossing cohort estimate — and was not dropped for being non-significant.
* **GSE21501** enters the five-cohort sensitivity at its frozen n = 129, ρ = +0.1712, which is
  the **highest** ρ of the five. It was not promoted to the four-cohort estimate for being the
  largest, and it was not demoted for being specimen-level: the frozen text places it in the
  five-cohort sensitivity, and that is where it is.
* GSE57495 (the lowest-ρ original cohort, ρ = +0.1231) and GSE85916 (ρ = +0.0780) are in every
  model. No cohort was dropped to narrow an interval.

## A.2 Bulk only

No single-cell and no snRNA cohort was pooled. GSE278688, GSE278689, GSE212966, GSE156405 and
GSE155698 appear nowhere in `meta_results.csv`; this is asserted, not merely intended
(`scripts/04_verify_summary.py`, scope check).

## A.3 The GSE21501 patient count is not determinable — carried as a bound

GEO exposes no patient identifier and no BioSample accession for GSE21501, so the 129 eligible
tumour arrays/specimens cannot be converted into a patient count. The estimate is therefore
labelled **specimen-level, n = 129, representing ≤ 129 patients** and is confined to the
five-cohort **sensitivity**. No exact patient count is asserted anywhere in this directory.

## A.4 The published three-cohort estimate is reproduced, not replaced

The three-cohort row in `meta_results.csv` is the same estimate as the manuscript's, rerun
through the manuscript's own function so that its machinery could be verified. It is **not** a
new primary and does not replace the published ρ = 0.132 (95% CI −0.110 to +0.360).

---

# PART B — FAILURES DURING THIS STEP

## B.1 The manuscript function does not return `Q` (crash on run 1)

**What failed.** `KeyError: 'Q'`. `random_effects_meta()` returns `I2_pct` but not the Q it was
computed from, and `meta_results.csv` requires a Q column.

**What was tried.** Reading Q out of the function was impossible without editing it, and editing
the archived module was not acceptable — the module hash `7b6a311a…` is the whole point of the
reuse.

**Resolution.** Q is computed by a helper using the **identical formula that sits inside the
function** (lines 418–421 of the module: `w = 1/v`, `v = 1/(n−3)`, fixed-effect μ, then
`Q = Σw(y−μ)²`). This is not an independent estimate that could disagree: it is checked against
the function's own returned `I2_pct`, which is a function of Q. That check (run log §4b) passes
for all three models to < 1e-9. Q is therefore the manuscript's Q, read out rather than assumed.

## B.2 A double-join bug in my own logging (run 1)

**What failed.** The cohort list printed as `T;C;G;A;_;P;A;A;D;…`. The function already returns
`cohorts` as a joined string; my log line joined it a second time, character by character.

**Resolution.** Fixed in the log line.

## B.3 The SAME double-join bug reached `meta_results.csv` — and a clean run did not catch it

**This is the most important failure recorded here.** The `cohorts_included` column of
`meta_results.csv` was written as `T;C;G;A;…` for every row. It was written twice: first by the
run that went on to crash later at `mu_fisher_z` (exit 1), and then again — unchanged — by the
next run, which **completed with exit code 0 and 49 checks passing**. No assertion looked at the
shape of that column. It was caught only by opening the written CSV and reading it.

**What this cost.** Two bad artefacts written to disk. The numerical columns were never affected:
the pooled ρ, CI, τ², I², Q and p in those same files were correct throughout, and the 49
passing checks were all genuinely passing — they simply did not cover the one column that was
wrong.

**Resolution.** The join was fixed at source, and the verifier (`04_verify_summary.py`) now
asserts (i) that every row's cohort list has exactly `k` entries and (ii) that the pooled set is
exactly the five frozen bulk cohorts. The lesson is recorded rather than absorbed: **a zero exit
code and a green check ledger are not evidence that an artefact is well formed**, and the only
reason this was found is that the output file was opened and read instead of being trusted.

## B.4 `mu` and `se` are also not returned (crash on run 2)

**What failed.** `KeyError: 'mu_fisher_z'`. `forest_inputs.csv` needs the z-scale mean and the
standard error actually used, neither of which the function returns.

**Resolution.** A helper reproduces lines 432–439 of the module (`w = 1/(v+τ²)`,
`μ = Σwy/Σw`, `se = √(max(1, Q_hk_raw)/Σw)`), validated by checking that `tanh(μ)` reproduces
the function's returned `r_random` (1e-15) for all three models, and independently by the
verifier's reconciliation of the forest weights against the se (exact to 1e-12).

## B.5 The anchor-reproduction tolerance was too tight, and was NOT simply loosened

**What failed.** With a 1e-12 tolerance the reproduction of the published interval **failed** on
both CI limits, `|Δ| = 2.983e-12` (ci_low) and `2.629e-12` (ci_high). The point estimate matched
to 2.8e-17 and τ² matched bit-for-bit, so this was not an input or model difference — but "it is
probably floating point" is not a finding.

**What was tried.** The residual was **attributed**. The per-cohort ρ values read from the
archive are bit-identical to those used here, so the difference had to be in the numerics. The
archive's implied critical value was recovered from its own stored interval,
`t = (μ − arctanh(ci_low_archive)) / se`, giving **4.302652729696140**, against this build's
`t₀.₉₇₅,₂ = 4.302652729749462` — a relative difference of **1.24e-11**, which is exactly the
scale of the observed residual.

**Resolution.** The tolerance is now 1e-9 **and** the attribution is itself asserted (the
implied-vs-actual t relative difference must be < 1e-9). The tolerance was widened only after
the cause was identified and pinned to a scipy-build difference in the t quantile. At every
published digit the reproduction is exact: 0.132 / −0.110 / +0.360.

## B.6 The hand-arithmetic re-derivation was mis-specified, and it was mine that was wrong

**What failed.** An independent hand re-derivation of the three-cohort interval returned
`max|Δ| = 9.568e-05` in z (8.5e-06 in ρ) against the function.

**What was actually wrong.** My re-derivation set τ² = 0, on the reasoning that the raw HK scale
is < 1 and Q < k−1. But the function carries its returned τ² = 5.628673589008228e-06 into the
weights, and at n = 60–178 that is enough to move the pooled z by ~1e-4. **The assumption was
wrong, not the function.**

**Resolution.** Re-derived carrying the function's own τ² → `max|Δ| = 0.000e+00`, and the hand
computation of the HK raw scale reproduces the function's to 0.0. The mis-specified first form is
retained as an observation in `run_log.txt` §5d and in `provenance/checks.json`.

## B.7 The materiality criterion for the controls was under-specified

**What failed.** The first criterion scored a convention as "materially different" only if it
moved ρ by > 0.01 or changed the zero-crossing verdict. Under it, the unmodified-HK control on
the **five-cohort** model scored "no material change" — even though it narrowed the interval from
0.213 to 0.078 while leaving ρ and the crossing verdict untouched.

**Resolution.** A third clause was added: the interval **width** changing by more than 25 %.
Re-run. The criterion was fixed; **the control was not relabelled to make anything pass**, and
the affected control (C3, five-cohort) still passes on its own merits (width × 0.366).

## B.8 The first independent-engine cross-check reported 6 failures out of 12

`metafor` 4.8.0 disagreed with the manuscript on τ² in all three models: it returns **exactly 0**,
the manuscript returns **5.628673589008228e-06**. Propagating that difference made the
reconstructed modified-HK interval disagree by 5.0e-05 to 9.6e-05. Nothing was relaxed; both
were attributed:

1. **τ².** The REML criterion is minimised at the lower boundary τ² = 0 in all three models
   (Q < k−1 in every case), which is confirmed by evaluating the criterion at 0, at the returned
   value, and at 1e-3 in each model (run log §4c). The manuscript's value is the bounded
   optimiser's tolerance-limited approach to that boundary (scipy `xatol = 1e-5`), which the
   module's own `τ² < 1e-8 → 0` snap does not catch. `metafor` reaches the boundary exactly.
   The check was restated at the tolerance that is meaningful (≤ `xatol`) **with the boundary
   fact asserted**, not deleted — and the fact is disclosed in `meta_results.csv` and in
   `results_summary.md` §2 because the stored τ² column would otherwise be read at face value.
2. **The CI gap.** Re-running the closed-form modified-HK rule in base R **with the manuscript's
   τ² carried** matches the manuscript's CI to **≤ 5.6e-17**. The full metafor-vs-manuscript gap
   is then asserted to equal the gap induced by the τ² difference alone (9.568e-05, 6.872e-05 and
   5.042e-05 for k = 3, 4, 5 respectively — equal to the last digit in each case).

Re-run: **15 checks, 0 failed**, including that both engines agree on every zero-crossing verdict
and that the unmodified-HK **control** values are metafor's own `knha` answer (read from
`verification_controls.csv`, not hard-coded).

## B.9 R emitted a package-build warning

`package 'Matrix' was built under R version 4.5.2` while running under R 4.5.1. Benign — it
affects no computed quantity, and `metafor` loaded and ran. Recorded only so that the warning in
the log is not mistaken for a failed step.

## B.10 Raw console captures preserved for the failures above

`run_log.txt` is append-only and holds every run that reached its final write step, so the
earlier **crashing** runs would otherwise exist only as prose in this file. Both raw console
captures are therefore kept:

| file | what it is |
|---|---|
| `provenance/superseded_run_20260928T124713_all_checks_green_but_csv_column_malformed.txt` | the full stdout of the run that reported **49 checks, 0 failed** and still wrote the malformed `cohorts_included` column (B.3). This is the evidence that a green ledger did not catch it. |
| `provenance/superseded_run_crash_KeyError_mu_fisher_z.txt` | the tail of the run that wrote `meta_results.csv` and then crashed at `mu_fisher_z` (B.4), which is how the malformed column first reached disk. |

Neither is a current artefact. Both are retained so the failure record can be checked against
raw output rather than taken on trust.

---

# PART C — CONTROLS THAT CAME BACK NULL (reported as null)

These were **expected to be materially different and were not**. They are reported as null
results, not converted into passes, and not quietly dropped from the control table.

| control | convention | Δρ vs frozen | width ratio | verdict |
|---|---|---|---|---|
| C1 | within-study variance `1/(n−1)` instead of `1/(n−3)` | −0.000135 | ×0.993 | **no material change** |
| C2 | within-study variance `1/n` instead of `1/(n−3)` | −0.000201 | ×0.989 | **no material change** |
| C7 | pool raw ρ on the ρ scale, no Fisher-z transform | −0.000086 | ×1.024 | **no material change** |

**What this means.** On these data the frozen estimate is insensitive to the within-study
variance offset and to the Fisher-z transform. Both are immaterial for the same two reasons: all
cohorts have large n (63–216, so `n−3`, `n−1` and `n` differ by ≤ 5 % in v), and every |ρ| < 0.18,
where `arctanh` is nearly the identity. **The convention that is load-bearing here is the
Hartung–Knapp rule, not the variance offset and not the transform** — and that is what the
controls that did fire (C3, C3b, C4, C5, C6, C8) demonstrate.

A control battery in which every control fires would have been the suspicious outcome. Three
nulls are reported as three nulls.

---

# PART D — WHAT COULD NOT BE DETERMINED

Carried in full in `results_summary.md` §6. Summarised so that nothing lives in only one place:

1. **Whether GSE205154's 216 specimens are 216 distinct patients.** 83 cores have no counterpart
   in GSE281005; if any of those split, n = 216 is an overcount. Not assessable from public
   metadata.
2. **How many independent patients GSE21501's 129 specimens represent.** Not determinable by
   construction — GEO exposes no patient identifier. This is why the estimate is specimen-level
   and confined to a sensitivity.
3. **Patient-level overlap between GSE205154 and TCGA-PAAD.** No shared identifier namespace, so
   independence cannot be verified. It bears on the four-cohort estimate's independence
   assumption and is carried as a limitation, not modelled.
4. **Any single-cell or snRNA contribution.** Out of scope by instruction; not pooled.
5. **Raw sequence data** for GSE205154 (controlled access, dbGaP `phs003597.v1.p1`). Not obtained
   and not needed for this step.

---

# PART E — CHOICES MADE, AND WHY NONE REQUIRED A STOP-AND-ASK

The frozen text did not need to be amended and no blocking question was raised. Each judgement
call below was resolved by **reusing an existing analysis identity** rather than by choosing,
and each is stated here so it is auditable.

| choice | resolution | why it is not a silent choice |
|---|---|---|
| which interval is "the" pooled CI | the function's own `ci_low`/`ci_high` (modified-HK t interval on Fisher z, back-transformed) | the archive's three-cohort row stores exactly those two columns and they equal the published −0.110/+0.360, so the archive settles it |
| which p-value | `p_modified_hartung_knapp` from the same function | it is the only p the function returns, and it reproduces the archived 0.143148 |
| which gene program | `Derosiers_ACP5_MMP9_core` | it is the module's own `PRIMARY` constant, and its three-cohort row is the published 0.132 |
| where Q comes from, the function not returning it | the identical in-function formula, validated against the function's own I² | the value is not chosen, only its read-out path; the check is in run log §4b |
| which τ² to store | the value the function returns, with the boundary fact disclosed in the notes and in the summary | changing it to 0 would be altering the frozen machinery's output; both readings are given |
| GSE21501's n | 129, as the frozen input fixes | the frozen text fixes the cohort as specimen-level; substituting an "effective n" would be reweighting beyond the frozen rule |

**The one judgement that carries real weight, stated plainly.** GSE21501's `v = 1/(n−3)` assumes
129 independent units. If the deposit in fact contains fewer than 129 patients, that variance is
understated and the five-cohort interval is too narrow — which is precisely why its lower limit
(+0.0317) should not be read as a clean exclusion of zero. **No effective-n correction was
applied**, because the frozen text does not authorise one and applying one would be reweighting a
cohort beyond the frozen rule. If the author wants that variant it is a separate, non-frozen
analysis, and it is flagged here rather than performed.

---

# PART F — LIMITATIONS OF THE VERIFICATION ITSELF

* The reproduction of the published anchor is exact at every published digit and attributed to
  ≤ 1.3e-11 relative at full precision. It is **not** an independent re-derivation of the
  manuscript's estimate from raw expression: it reruns the manuscript's function on the
  manuscript's own archived per-cohort values. A defect in `bulk_associations.csv` would survive
  this check. What it does establish is that the extension models use the same machinery as the
  published estimate, which is the property the frozen rule requires.
* The R cross-check is independent in implementation (different language, different package,
  different optimiser) but **not independent in data**: both engines read the same per-cohort ρ
  and n. It verifies the meta-analytic arithmetic, not the cohort estimates themselves. Those
  were verified upstream, in `phaseB_GSE205154/step3c_rho/` and `phaseB_GSE21501/step3c_rho/`.
* The control battery establishes that specific wrong conventions would have changed the answer.
  It does not establish that no other wrong convention would change it, and the three null
  controls show that some wrong conventions would not.
* `verification_controls.csv` is written by the same script that writes the estimates, so a
  systematic error in the replica would affect both. That risk is mitigated — not eliminated — by
  control C0 (replica ≡ imported original, max|Δ| = 0 on all three models) and by the R
  cross-check, which reimplements the modified-HK rule from its description and agrees to
  ≤ 5.6e-17.

---

# PART G — NOTHING WAS EXCLUDED FOR CONVENIENCE

No cohort, specimen, representation or analysis was dropped, replaced or re-specified at any
point in response to a result. The three models are the three the frozen text names, in the
order it names them. Both of the estimates that could have been collapsed into one — the
four-cohort and the five-cohort — are reported side by side, and **neither is selected as the
correct answer**, because the frozen text forbids exactly that selection when they disagree.

**No manuscript text, figure or claim was written or modified. All output is confined to
`BJC_EXT_v1_20260927/extension_meta/`.**
