# Verification and controls -- GSE21501 step 3c-2

Generated 2026-09-28 12:27:21. Every row below is reproduced from the machine-readable files
`output/manuscript_verification.csv` and `output/verification_controls.csv`.

## Part 1 -- the manuscript's own scoring function, executed on this cohort

Stage 2 re-implements the frozen scoring and resampling rules by transcription.
A transcription that is wrong in the same way as its original is invisible to
internal checking, so the archived manuscript module
`analysis_F_cross_scale.py` was imported and **its** functions were run on **this**
cohort's data. The module was read, never edited.

* path: `E:\Tex\HMlong\npjPrecOncol\FINAL_SUBMISSION_v6\BMC_MEDICINE_SUBMISSION_v47\_archives_translate\Additional_file_2_Analysis_archive_v39\scripts\analysis_F_cross_scale.py`
* sha256: `7b6a311a0a73835f66596b7ea21070a6355c1fe6de44b518394dfe89f7c4334e`
* verification population: n = 129 (DUP-2 applied)
* seed 15072026, 10000 resamples

| id | check | expectation | observed | verdict |
|---|---|---|---|---|
| V1 | archived manuscript module loaded | imports without calling main() | ok, sha256=7b6a311a0a73835f... from frozen document's home tree (BMC_MEDICINE_SUBMISSION_v47) | **PASS** |
| V2 | the four frozen functions exist in the archived module | all present | {'zscore_columns': True, 'score_program': True, 'required_gene_count': True, 'spearman_bootstrap': True, 'permutation_p': True} | **PASS** |
| V3 | zscore_columns agrees cell-for-cell on the real 129 x 3 matrix | max abs difference == 0 | max|diff| = 0.000e+00 | **PASS** |
| V4 | score_program produces the identical score vector | max abs difference == 0 and identical gene list | max|diff| = 0.000e+00; genes ['ACP5', 'MMP9'] vs ['ACP5', 'MMP9'] | **PASS** |
| V5 | spearman_bootstrap identical given the same seed | rho, asymptotic p and both CI limits agree exactly | rho +0.171188105944 vs +0.171188105944; CI [-0.012823752,+0.340530829] vs [-0.012823752,+0.340530829] | **PASS** |
| V6 | permutation_p identical given the same seed | p agrees exactly | 0.0506949305 vs 0.0506949305 | **PASS** |
| V7 | required_gene_count(PRIMARY) == 2 | == 2 | 2 | **PASS** |
| V8 | the archived module reproduces the stage-2 primary rho | identical to the value written to primary_result.csv | manuscript-module rho = +0.171188105944 vs stage-2 rho = +0.171188105944 | **PASS** |
| V9 | frozen manifest hash still matches after the whole run | == the frozen value | 3e11938ae5c7e682c9e53aad93d051579764a4ec8173d69d571dbeda0f3ac4e4 | **PASS** |

V3--V6 run on the real data with the real seed, so a divergence in the collapse
rule, the ddof, the nanmean, the clip or the resampling scheme would surface as a
numeric mismatch. V8 is the end-to-end statement: the manuscript's own scoring
function, run on this cohort, returns exactly the rho written to
`output/primary_result.csv`.

## Part 2 -- controls on this run's own pipeline

| control | kind | expectation | observed | verdict |
|---|---|---|---|---|
| C7-rng-addressing | self-consistency | rho is EXACTLY invariant to how the generator is addressed (no generator enters rho); the bootstrap interval and the permutation p are resampled quantities and may differ only by Monte-Carlo error | rho identical=True; CI gap=0.003213; permutation p gap=0.0052 (1 s.e. = 0.0022) | **PASS** |
| C1-self-correlation | positive | rho(score, score) == 1 exactly | rho=1.000000000000 | **PASS** |
| C2-monotone-invariance | positive | rho unchanged by any strictly monotone transform of SIGLEC15 | rho(orig)=+0.171188105944 vs rho(exp)=+0.171188105944 | **PASS** |
| C3-circular-rotation-null | negative | a 7-position rotation of SIGLEC15 destroys the association | rho=-0.031678 vs observed +0.171188 | **PASS** |
| C4-label-scramble-null | negative | randomly permuting the score across arrays removes the association | rho=-0.045095 | **PASS** |
| C5-seed-sensitivity-and-reproducibility | negative | a different resampling seed must move the interval, AND the same seed must reproduce it exactly | seed 12345 twice: (-0.007953377687870697, 0.342038162880203) == (-0.007953377687870697, 0.342038162880203) -> True; seed 999: (-0.010195331375882594, 0.3427851247040422) differs -> True | **PASS** |
| C6-zero-variance-column | negative | a constant programme gene is set aside by zscore_columns, so the score is carried by ACP5 alone (finite, mean 0) | score mean=-9.639e-17, n_finite=129 | **PASS** |
| C8-null-sampling-distribution | negative | 2000 independent permutations give a null centred at 0 with sd ~ 1/sqrt(n-1) = 0.0884, and an empirical p within 5 s.e. (0.0232) of the reported permutation p | null mean=-0.00568, sd=0.08615, empirical p=0.0450 vs reported 0.0527 (gap 0.0077, 5 s.e. = 0.0232) | **PASS** |
| C11-S4-averaging-is-not-a-no-op | diagnostic | S4 must differ from the primary UNLESS every gene's group mean happens to equal the retained array's own value | rho(S4)=+0.171188 vs rho(primary)=+0.171188; identical=True | **PASS** |
| C9-unit-arithmetic | structural | DUP-2 leaves 129 rows, one per distinct specimen key; no exact patient count is asserted anywhere | rows=129, distinct specimen keys=129, providers=5 | **PASS** |
| C10-subsets-of-the-129 | structural | the 102 and every drop-set sensitivity are subsets of the same 129 | 102 subset=True; all sensitivities nested=True | **PASS** |
| C12-platform-row-arithmetic | structural | platform rows = symbol-bearing rows + rows without a usable symbol; the latter equals the matrix IDs the annotated table cannot name; and every symbol-bearing platform row is present in the matrix | 45220 = 32995 + 12225; matrix IDs not nameable = 12225; symbol-bearing platform rows absent from matrix = 0; 32995 + 12225 = 45220 matrix features | **PASS** |

### Which of these could actually have failed

A control that cannot fail is decoration. The ones below are the load-bearing
ones, each with the specific failure it would have caught:

* **C3 circular rotation** -- rho fell from +0.1712 to -0.0317. A pipeline that
  had manufactured the association from the *order* of the array columns would
  have returned the same value here, because a rotation preserves nothing about
  the pairings except their order.
* **C4 label scramble** -- -0.0451, using the seeded generator rather than array
  order, so it is not a restatement of C3.
* **C5 seed sensitivity and reproducibility** -- the same seed reproduces the
  interval exactly, and a different seed moves it. If the seed argument had been
  inert, every resampled number in this run would have been unanchored and this
  control would have failed.
* **C6 zero-variance programme gene** -- exercises the `sd == 0 -> NaN` branch of
  the frozen `zscore_columns`, the one scoring detail the manuscript records as
  locked-day implementation rather than protocol text.
* **C8 null sampling distribution** -- 2000 independent permutations gave a null
  with mean -0.0057 and sd 0.0862 against a theoretical 1/sqrt(n-1) = 0.0884. A
  null centred away from 0 would have indicated a mis-specified statistic.
* **C10 nesting** -- the 102 and every drop-set sensitivity are subsets of the
  same 129, which guards against a sensitivity silently computed on a different
  cohort, for instance one that reintroduced a DUP-2 duplicate.

### A control that was wrong, and what happened

C7 first required the bootstrap interval to be **bit-identical** under two
different RNG seeds. That is unsatisfiable: two different seeds cannot produce
identical resampled intervals. The script treated the mismatch as fatal and
exited 1 rather than reporting a number it could not stand behind. The check was
then corrected to test what the frozen rule actually implies -- `rho` is exactly
invariant to how the generator is addressed, since no generator enters it, while
the resampled quantities may differ only by Monte-Carlo error. The correction is
recorded in the run log and in `failed_or_excluded.md`, not silently applied.

