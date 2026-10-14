# Implementation audit

## Changes made after the first complete run

1. A sample with no detected SIGLEC15-positive myeloid cell (GSE156405, LuM) was
   initially omitted from both the within-sample cell contrast and the sample-level
   aggregate. This was incorrect for the latter estimand. The code was corrected so
   that LuM remains non-estimable only for the positive-versus-negative cell contrast
   and is included in the aggregate analysis. Final counts are 19 estimable samples
   for the cell contrast and 20 for the aggregate analysis.
2. The locked protocol anticipated three estimable spatial cases. Under the frozen
   capture threshold, case C1 had only one SIGLEC15-positive myeloid-labelled spot
   across 379 such spots and no eligible capture. The spatial contrast is therefore
   estimable in two cases (four captures), not three. C1 is retained in the
   eligibility table and its non-estimability is reported; it is not silently
   discarded.
3. The first implementation used the unmodified Hartung--Knapp variance scale. With
   three homogeneous cohorts its scale was below one, producing an implausibly narrow
   interval. The final implementation uses the standard modified Hartung--Knapp
   safeguard `max(1, q)`. This is a conservative inferential correction and does not
   change cohort estimates, gene sets or included samples. Both the raw and used
   scale factors are exported.
4. After observing a weak positive bulk association, a pan-myeloid-adjusted analysis
   was added and permanently labelled post hoc. Its gene set and rationale are frozen
   in `POSTHOC_SENSITIVITY_ADDENDUM_2026-07-22.md`.

No association-driven sample deletion, gene replacement or endpoint switching was
performed.
