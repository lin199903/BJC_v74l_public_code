# BJC v72 burden-estimand and sparsity analysis contract

Status: **FROZEN BEFORE OUTCOME INSPECTION**  
Freeze date: 2026-10-01 (Asia/Shanghai)  
Parent release: `SUBMISSION_PACKAGE_BJC_v71`  
Branch: `BJC_v72_burden_estimand`  

## 1. Scientific question

Does a tumour-level SIGLEC15 RNA summary preserve the between-tumour ordering of an ACP5/MMP9-high myeloid-state burden measured in the same cell-resolved tumour samples?

This question is distinct from the previously reported within-sample cellular coupling, \(K_i\), which compares ACP5/MMP9 scores between SIGLEC15 RNA-detectable and non-detectable myeloid cells. The new reference quantity is a tumour-sample burden. It is not \(K_i\), target prevalence, or treatment response.

## 2. Population, unit and eligibility

- Population: tumour-derived myeloid cells in GSE155698, GSE156405 and GSE212966.
- Inferential unit: biological sample, because a complete sample-to-patient map is unavailable.
- Primary eligibility: at least 50 tumour-myeloid cells, matching the locked cross-scale pipeline.
- Primary analysis set: all eligible samples, including samples with no SIGLEC15 RNA-detectable myeloid cell.
- Dataset is a design stratum. No cell is treated as an independent population-level replicate.

## 3. Frozen primary reference estimand

For sample \(i\) and gene \(g\), aggregate expression over all eligible tumour-myeloid cells is

\[
A_{ig}=\log\left[1+\operatorname{mean}_{c\in i}\{\exp(x_{icg})-1\}\right].
\]

Within each dataset, \(A_{i,ACP5}\) and \(A_{i,MMP9}\) are standardised across samples. The burden reference is

\[
B_i=\operatorname{mean}\{z(A_{i,ACP5}),z(A_{i,MMP9})\}.
\]

`B_i` is the **sample-level myeloid ACP5/MMP9 aggregate burden**. It is the primary tumour-ordering reference for v72.

The mean of the existing within-sample cell-standardised programme score is prohibited as a burden estimand because its sample mean is approximately zero by construction. No threshold-defined “high-state fraction” is primary because the stored cell scores do not provide a prespecified cross-sample common scale.

## 4. Frozen target readouts

Primary proxy:

- `R1`: myeloid SIGLEC15 aggregate, constructed with the same inverse-log aggregation over all eligible tumour-myeloid cells.

Secondary decompositions:

- `R2`: fraction of eligible tumour-myeloid cells with stored SIGLEC15 expression greater than zero.
- `R3`: mean stored SIGLEC15 expression among detectable cells; undefined when no cell is detectable.
- `R4`: all-tumour-cell SIGLEC15 aggregate using the same inverse-log aggregation.

These readouts are not interchangeable. `R2` and `R3` are decomposition variables; `R4` changes compartment.

## 5. Primary association and inference

- Estimand: pooled within-dataset Spearman rank correlation between `R1` and `B_i`.
- Pooling: rank and standardise both variables within each dataset, concatenate the standardised ranks, and calculate their correlation. This equals a sample-size-weighted combination of the dataset-specific rank correlations when strata are complete.
- Uncertainty: 10,000 dataset-stratified sample bootstraps, percentile 95% interval.
- Null test: 10,000 within-dataset permutations of the burden ordering relative to the proxy; two-sided Monte Carlo p value.
- Dataset-specific correlations are descriptive because the strata contain only 5–9 samples.

## 6. Prespecified sparsity analyses

The following analyses are sensitivities, not replacement primary estimands:

1. Re-estimate the `R1`–`B_i` pooled rank correlation after requiring at least 5 detectable SIGLEC15-positive myeloid cells per sample.
2. Repeat after requiring at least 10 detectable cells.
3. Report retained sample counts and dataset composition at every threshold.
4. Run a cell-exchangeability null simulation within each dataset: permute observed SIGLEC15 cell-level values across all eligible tumour-myeloid cells, retain each sample's cell count and the dataset-wide sparsity/intensity distribution, reconstruct `R1`, and recalculate the pooled `R1`–`B_i` correlation for 10,000 draws.

The simulation tests whether the observed correlation is unusual under exchangeable sparse detection. It does not identify a biological mechanism and does not prove that sparsity caused or did not cause the observed association.

## 7. Secondary analyses

- Estimate `R2`, `R3` and `R4` against `B_i` using the same stratified rank framework where estimable.
- Report leave-one-dataset-out estimates for `R1`–`B_i`.
- Verify exact numerical reproduction of the archived 20-sample `R1`–`B_i` result before accepting any new output.
- Preserve the existing `K_i` analyses unchanged and explicitly separate them from the new burden analysis.

No reliability-weighted correlation is added. Weighting by positive-cell count would redefine the sample-level estimand and give larger target-positive samples more influence without a validated measurement-error model.

## 8. Claim hierarchy and ceiling

Permitted primary conclusion:

> Tumour-level SIGLEC15 RNA did not establish stable preservation of the between-tumour ordering of the myeloid ACP5/MMP9 aggregate burden across the available samples and sparsity specifications.

Permitted stronger wording only if all prespecified estimates and intervals justify it: none by default.

Prohibited inferences:

- no association;
- formal non-observability or a non-injective measurement function;
- causation by sparsity, composition or spatial heterogeneity;
- SIGLEC15 protein abundance, target engagement or treatment response;
- clinical utility or predictive value.

## 9. Source lock and production rules

- All v72 numbers must be generated by the frozen v72 script and carried into tables, prose and gates from machine-readable outputs.
- The v71 package is read-only. v72 is assembled in separate work, build and package directories.
- Methods state the estimand before Results report it. Results distinguish `B_i` from `K_i`. Discussion interprets the measurement consequence without re-litigating the analysis identity.
- New analyses enter the supplementary archive and SI before or at the same build as the revised main text.
- Manuscript, SI, cover letter, manifest, release gate and rendered QA must all refer to the same final bytes.

## 10. Acceptance gates

The branch cannot be released unless:

1. archived `R1`–`B_i` result is reproduced within numerical tolerance;
2. all prespecified outputs exist and are internally consistent;
3. every reported sample count is traceable to the sample-level output;
4. no sentence conflates `K_i` with `B_i`;
5. v71 upload bytes remain unchanged;
6. the v72 upload set contains exactly 11 files and the package manifest matches their final hashes;
7. DOCX and PDF rendering checks pass on the final package files.
