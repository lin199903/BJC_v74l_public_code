# results_summary.md — secondary bulk-extension meta-analysis

**Status: revision-stage SECONDARY estimates.** The published three-cohort bulk estimate
(ρ = 0.132, 95% CI −0.110 to +0.360) is **reproduced here, unchanged and unreplaced**. The
four- and five-cohort estimates below are additional and do not supersede it.

**Scope.** Bulk cohorts only — no single-cell and no snRNA cohort is pooled. No cohort was
excluded, reweighted or reshaped on the basis of its ρ magnitude or direction. No manuscript
text, figure or claim was written or modified by this step.

**Machinery.** The manuscript's own `random_effects_meta()` was imported by path from the v39
analysis archive (`analysis_F_cross_scale.py`, sha256 `7b6a311a…`, function at lines 413–454,
function source sha256 `1b501f10…`) and called directly: Fisher z → REML τ² → modified
Hartung–Knapp. The three original cohorts' ρ and n were read from the archive's own output
(`results/bulk_associations.csv`, sha256 `91761ef0…`), not re-derived.

---

## 1. The three estimates

| estimate | k | cohorts | pooled ρ | 95% CI | crosses zero? |
|---|---|---|---|---|---|
| 3-cohort published — **reproduced** | 3 | TCGA_PAAD, GSE57495, GSE85916 | 0.132332 | −0.110117 to +0.359917 | yes |
| **4-cohort secondary extension** | 4 | + GSE205154 | **0.132074** | **−0.006107 to +0.265305** | **yes** — lower limit −0.006107 |
| **5-cohort sensitivity** | 5 | + GSE21501 | **0.139682** | **+0.031725 to +0.244417** | **no** — lower limit +0.031725 |

**The four- and five-cohort intervals disagree about excluding zero.** The four-cohort interval
still crosses zero, by 0.0061; the five-cohort interval excludes zero, by 0.0317. Per the frozen
rule, **neither is selected as the correct answer** — the disagreement is itself the
inferential-unit / cohort-composition sensitivity, and both are reported side by side.

## 2. Heterogeneity

| estimate | τ² (Fisher z) | I² | Q | p (modified HK) |
|---|---|---|---|---|
| 3-cohort (reproduced) | 5.628673589008228e-06 | 0.0 % | 0.369367 | 0.143148 |
| 4-cohort | 5.628673589008228e-06 | 0.0 % | 0.369422 | 0.055754 |
| 5-cohort | 5.628673589008228e-06 | 0.0 % | 0.532291 | 0.023048 |

**Read τ² as 0 at the REML lower boundary, not as 5.63e-06.** In all three models Q < k−1, so
the REML criterion is minimised at τ² = 0; the stored value is the bounded optimiser's
tolerance-limited approach to that boundary, and the manuscript's own `τ² < 1e-8 → 0` snap does
not catch it. `metafor` 4.8.0 returns exactly 0 for all three models (run log §4c;
`provenance/crosscheck_metafor.csv`). Carrying 0 instead of 5.63e-06 moves each interval by
≤ 7e-5 and reproduces every zero-crossing verdict below. **I² = 0 %: no detectable
between-cohort heterogeneity in any model.** Q rises from 0.369367 (k=3) to 0.369422 (k=4) to
0.532291 (k=5), still far below k−1 in every case.

## 3. Why the two intervals disagree

Three facts, all from the run log, none of them a preference:

1. **GSE205154 adds weight but no new signal.** It carries **40.56 %** of the four-cohort weight
   (212.74 of 524.52) and its ρ (+0.1317) is almost exactly the three-cohort pooled value, so the
   four-cohort point estimate barely moves (0.132332 → 0.132074). What changes is precision: k
   rises 3 → 4 and the wider t quantile is replaced by a narrower one, so the interval narrows
   from 0.470 wide to 0.271 wide — but its lower limit lands at −0.0061, just inside zero.
2. **GSE21501 is the marginal inclusion.** It carries **19.36 %** of the five-cohort weight
   (125.91 of 650.43) and is the highest-ρ cohort (+0.1712), so adding it moves the point
   estimate up to 0.139682 and lifts the lower limit clear of zero. This is a cohort-composition
   effect, not a precision effect, and it is the effect the frozen text anticipated.
3. **The modified Hartung–Knapp floor is binding in both new models.** The raw HK scale is
   0.123063 (four-cohort) and 0.132994 (five-cohort), both below 1, so the rule substitutes 1 and
   both intervals are drawn at t₀.₉₇₅,k₋₁ width. Control C3 (run log §5b) shows this is exactly
   what keeps the four-cohort interval crossing zero: with the *unmodified* HK variance the same
   four cohorts give +0.083906 to +0.179626, which excludes zero. **No estimate in
   `meta_results.csv` uses the unmodified rule** — the frozen modified rule is what is reported —
   but the fact that it is decision-relevant here is recorded rather than left implicit. The same
   control applied to the published three-cohort model (C3b) gives +0.028418 to +0.233417, so the
   published anchor's own zero-crossing is likewise a property of this frozen convention.

## 4. GSE21501 carries less inferential weight, and this is stated

> The GSE21501 analysis was specimen-level because patient identifiers were not deposited; the number of independent patients could therefore not be established exactly.

Its unit is the **specimen**, not the patient: n = 129 eligible tumour arrays/specimens
representing **≤ 129 patients**, and no exact patient count is asserted. Its variance
v = 1/(n−3) therefore assumes 129 independent units, which cannot be verified from the deposited
metadata. **It is not a patient-level cohort and does not carry the same inferential weight as
the other four**; the five-cohort estimate is consequently a sensitivity analysis, and the
four-cohort estimate is the one built only from patient-level cohorts. GSE21501 carries 19.36 %
of the five-cohort weight on an assumption of independence that the deposit cannot confirm.

## 5. Verification (detail in `run_log.txt`, `verification_controls.csv`)

* **Anchor reproduces exactly**: ρ 0.13233228554949808 vs archive 0.1323322855494981
  (|Δ| 2.8e-17); CI to |Δ| ≤ 3.0e-12. The residual is attributed, not waved away — the archive's
  implied t₀.₉₇₅,₂ is 4.302652729696140 against this build's 4.302652729749462 (relative
  1.24e-11), i.e. a scipy-build difference in the t quantile. At every published digit the
  reproduction is exact: 0.132, −0.110, +0.360.
* **Independent engine**: `metafor` 4.8.0 reproduces Q to 2.2e-16 and I² exactly, and the
  closed-form modified-HK rule reimplemented in base R matches the manuscript's CI to ≤ 5.6e-17
  once the manuscript's τ² is carried (15/15 checks; `provenance/crosscheck_metafor_log.txt`).
* **Controls capable of failing** (`verification_controls.csv`, 13 rows): the **unmodified HK**
  control changes the four-cohort zero-crossing verdict (C3); the **conventional normal interval**
  does the same (C4); a **wrong sign convention** flips the pooled estimate (C5) or moves it by
  0.086 when applied to one cohort only (C6); a **perturbed input** fails the anchor reproduction
  (C8, |Δ| 5.8e-03). The replica used for these controls is validated against the imported
  original first (C0, max|Δ| = 0 on all three models).
* **Controls that came back NULL** — reported as null, not hidden: the within-study variance rule
  (1/(n−1) or 1/n; |Δρ| ≤ 2.0e-04) and dropping the Fisher-z transform (|Δρ| 8.6e-05) make **no
  material difference** here, because all four/five cohorts have large n, τ² is at the boundary,
  and all |ρ| < 0.18 where arctanh is nearly the identity. The load-bearing convention is the
  Hartung–Knapp rule, not the variance offset or the transform.

## 6. Not determined

* Whether GSE205154's 216 specimens are 216 distinct patients (83 cores lack a counterpart in
  GSE281005; if cores split, 216 is an overcount). Public metadata cannot settle it.
* Whether GSE21501's 129 specimens come from fewer than 129 patients — **not determinable by
  construction**, which is why the estimate is labelled specimen-level and reported as a
  sensitivity only.
* Patient-level overlap between GSE205154 and TCGA-PAAD; not assessable, no shared identifier
  namespace. It bears on independence and is carried as a limitation, not modelled.
* Anything about the single-cell or snRNA cohorts: out of scope, not pooled.

*Provenance: `meta_results.csv`, `forest_inputs.csv`, `verification_controls.csv`, `run_log.txt`,
`sessionInfo.txt`, `provenance/`.*
