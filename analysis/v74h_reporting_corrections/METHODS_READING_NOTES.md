# Executable reading notes for frozen methods

## TREM2 subsampling

In `scripts/run_S22_dilution.py`, retained detectable-cell count is `min(n_positive, max(1, round(f*n_original)))`, for nominal f=0.02 or 0.01. All non-detectable cells are retained. The achieved fraction is therefore `k/(n_negative+k)`. Both target and leave-target-out programme are reaggregated from the retained set. Across archived samples, retained fractions are 1.37–6.79% for the 2% original-cell quota and 0.98–3.82% for the 1% quota.

The companion denominator table is derived deterministically from archived sample counts. It is not a new simulation. The archived 0.165/0.118 means and 5th–95th simulation percentiles remain unchanged. The experiment changes composition and reference together; it does not isolate target capture loss or identify a biological mechanism.

## Within-dataset rank pooling and bootstrap

`v72_burden/analysis/run_v72_burden_analysis.py` uses average ranks for ties and population standard deviations (ddof=0) before concatenation. For estimable strata this gives sample-count weights, 5/20, 9/20 and 6/20 in the primary burden analysis. It is not Fisher-z meta-analysis or n−1 pooling, and the separate K_i correlation does not use this B_i estimand.

The full eligible dataset supplies gene-standardisation parameters and sample burdens. Pair bootstrap holds that mapping fixed, recomputes ranks within resamples and omits strata with fewer than two samples or constant ranked variables. Draws with no estimable stratum are discarded. The ten-detectable-cell sensitivity contains 9,979 finite estimates from 10,000 draws and a percentile interval of −1 to 1. These notes clarify the archived conditional procedure; they do not replace it.

## Two-cohort leave-one-out sensitivity

The archived leave-one-cohort-out bulk sensitivity uses fixed-effect inverse-variance pooling for the remaining two cohorts. It changes both cohort composition and pooling model. Two-cohort random-effects models are not structurally impossible; no new random-effects sensitivity is introduced here.
