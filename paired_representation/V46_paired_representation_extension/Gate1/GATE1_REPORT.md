# Gate 1 report — paired representation test

**Status: COMPLETE. Correctness gate EXACT.**

## Correctness gate

The single-cell aggregation was recomputed from the three files listed in the
archived `input_manifest.csv` (SHA-256 verified against the manifest) using the
frozen selection rule (`tissue in [PDAC, PDAC_tumor, PDAC_tissue]` and
`celltype == Myeloid`, >=50 myeloid cells per sample). The recomputed myeloid
aggregate reproduces the archived values with a maximum absolute difference of
9.8e-17 across 20 samples
(myeloid counts 20/20 and detectable-cell counts 20/20 identical).

Cell totals: 4,081 (GSE212966) + 4,635 (GSE156405) + 1,715 (GSE155698) =
**10,431 tumour-myeloid cells across 20 samples**, matching the manuscript.

## Readouts tested

| Readout | Definition |
|---|---|
| R1 | myeloid aggregate SIGLEC15 — closest analogue of a scalable RNA readout |
| R2 | SIGLEC15-detectable fraction of included tumour-myeloid cells |
| R3 | mean SIGLEC15 among detectable tumour-myeloid cells |
| R4 | all-cell pseudo-whole-tumour aggregate over the tumour compartment |

K_i is the locked within-patient detectable-minus-non-detectable ACP5/MMP9
contrast (19 estimable samples).

## Results

| Readout | n | Spearman rho | 95% CI | CI upper bound | LO direction |
|---|---:|---:|---|---:|---|
| myeloid aggregate SIGLEC15 | 19 | -0.363 | [-0.765, +0.105] | +0.105 | retained |
| SIGLEC15-detectable cell fraction | 19 | -0.184 | [-0.560, +0.254] | +0.254 | sign reversal |
| mean SIGLEC15 among detectable cells | 19 | -0.358 | [-0.709, +0.139] | +0.139 | retained |
| all-cell pseudo-whole-tumour aggregate | 19 | +0.075 | [-0.401, +0.520] | +0.520 | sign reversal |

Per-dataset Spearman rho:

| Readout | GSE155698 (n=5) | GSE156405 (n=8) | GSE212966 (n=6) |
|---|---:|---:|---:|
| myeloid aggregate SIGLEC15 | +0.30 | -0.55 | -0.14 |
| SIGLEC15-detectable cell fraction | +0.30 | -0.40 | +0.37 |
| mean SIGLEC15 among detectable cells | -0.90 | -0.71 | -0.26 |
| all-cell pseudo-whole-tumour aggregate | +0.90 | +0.10 | -0.77 |

## Interpretation

1. **Both continuous target-quantity readouts exclude moderate-to-strong
   positive concordance.** R1 and R3 have 95% CI upper bounds of
   +0.105 and +0.139. On the same
   patient samples, a higher target quantity did not accompany a larger
   cellular contrast, and the data are inconsistent with a positive
   association above roughly 0.15.
2. **R3 is the most internally consistent readout**, being negative in all
   three series separately (-0.90, -0.71, -0.26) with direction retained
   under leave-one-dataset-out.
3. **R1 retains direction** under leave-one-dataset-out (all three omissions
   negative), so its negative point estimate is not driven by one series.
4. **The detectable fraction (R2) is unstable** (one sign reversal), so no
   directional statement is supported for prevalence alone.
5. **The all-cell aggregate (R4) is uninformative** (rho +0.075,
   CI [-0.401, +0.520]) and does not support strong positive concordance
   either.

## Claim ceiling

This is a paired comparison **within the same public single-cell samples**; it
is not an independent clinical assay pairing, which would require a scalable
readout and the cellular reference measured in the same specimen on different
platforms. The supported statement is:

> In the same patient samples, neither of the two continuous target-quantity
> readouts showed positive concordance with the cellular state, and
> moderate-to-strong positive concordance is excluded; prevalence and the
> all-cell aggregate were uninformative.

Not supported: that any readout has been validated, qualified, or shown to be
interchangeable; and not supported: a biological negative relation.

## What this changes in the manuscript

v45 Discussion states that the unresolved quantity is same-patient concordance
and that paired specimens are required. Gate 1 supplies a within-sample paired
test on the existing public data, which partially closes that gap at the RNA
level while leaving the clinical-assay pairing open. Recommended placement:

- **Results**: 2-3 sentences reporting R1 and R3 with CIs and the exclusion
  bound; R2/R4 as uninformative; pointing to a new supplementary table.
- **Discussion**: revise the 'unresolved quantity' paragraph so that the open
  item becomes independent clinical-assay pairing rather than any same-patient
  comparison at all.
- **Conclusions**: unchanged claim ceiling; optionally note that the
  within-sample RNA comparison did not support interchangeability.
- **Supplementary**: new table with all four readouts, per-dataset estimates
  and leave-one-dataset-out behaviour.
