# Locked protocol: cross-scale SIGLEC15--osteoclast/ECM analysis

Lock date: 2026-07-22 (Asia/Shanghai)

Status: prospectively frozen for the new cross-scale analysis before inspecting any
SIGLEC15--program association estimates. This is not a registered protocol and must
not be described as preregistration.

## Scientific question

Derosiers et al. (Glycobiology, 2026; DOI: 10.1093/glycob/cwag035) reported that
SIGLEC15 marks an osteoclast-like, matrix-remodelling tumour-associated myeloid
state. The public abstract specifically names DAP12-dependent induction of ACP5 and
MMP9 and secretion of IL-1beta and IL-6. It does not expose a portable fixed gene
signature. We therefore test whether this article-anchored state is detectable in
public PDAC data, and whether its association with SIGLEC15 is preserved as data are
aggregated from myeloid cells to samples, bulk tumours and spatial regions.

## Frozen programs

SIGLEC15 is excluded from every program score.

1. Primary, article-anchored functional core: `ACP5`, `MMP9`.
2. Secondary canonical osteoclast panel: `ACP5`, `CTSK`, `MMP9`, `ATP6V0D2`,
   `DCSTAMP`, `OCSTAMP`, `ITGB3`, `CALCR`, `NFATC1`.
3. Secondary compact ECM-degradation effector panel: `MMP2`, `MMP9`, `MMP12`,
   `MMP14`, `CTSB`, `CTSK`, `CTSL`, `PLAU`, `PLAUR`.
4. Descriptive inflammatory outputs only: `IL1B`, `IL6`; these are not combined
   with the primary score and are not confirmatory endpoints.

The two expanded panels are literature-informed sensitivity panels, not signatures
claimed to have been published by Derosiers et al. Their exact membership is frozen
here before association estimation.

## Data contract

- Independent scRNA-seq series: GSE212966, GSE156405 and GSE155698. The integrated
  GSE156405/GSE155698 object is excluded to prevent duplicate counting.
- scRNA inclusion: exact tumour labels `PDAC`, `PDAC_tumor` or `PDAC_tissue`; exact
  top-level cell label `Myeloid`; sample with at least 50 included myeloid cells.
- scRNA expression: stored natural-log `log1p` normalized matrix. A cell is
  SIGLEC15-positive when stored SIGLEC15 expression is greater than zero.
- Myeloid sample aggregate: for each gene, `log1p(mean(expm1(log_expression)))`
  across included cells. Because raw molecule counts are unavailable in these
  processed objects, this is explicitly described as a cell-equal pseudobulk-like
  aggregate rather than a count-summed pseudobulk.
- Bulk cohorts: TCGA-PAAD (178 tumours), GSE57495 (63) and GSE85916 (80). GEO probe
  identifiers are mapped with their deposited SOFT platform annotations; multiple
  probes per gene are averaged exactly as in the existing external-validation code.
- Spatial dataset: GSE282302, three post-neoadjuvant resected human PDAC cases,
  processed public object with case, capture (`sample_id`) and spot annotations.
  The primary spatial analysis is restricted to spots labelled `Myeloid`.

No sample is excluded on the basis of an observed association or outcome.

## Frozen scoring

- At cell/spot level, each included gene is z-scored within its biological sample
  (scRNA sample or spatial capture), then averaged. At least two genes must be
  present for a multi-gene program; the two-gene primary score requires both genes.
- At sample/bulk level, each gene is z-scored across samples within its dataset,
  then averaged. The primary score requires both genes. Secondary scores require at
  least 60% of their frozen genes and at least three genes.
- Program scores are leave-SIGLEC15-out by construction.

## Estimands and inference

### Layer 1: SIGLEC15-positive myeloid cells

For each eligible tumour sample, estimate the difference in mean program score
between SIGLEC15-positive and SIGLEC15-negative myeloid cells. The primary summary
is the equally weighted mean of these sample-level differences, with a 95% CI from
a dataset-stratified patient bootstrap (10,000 replicates; seed 15072026). A
two-sided Wilcoxon signed-rank test across sample differences is reported. Cells are
never treated as independent inferential units.

Detection-depth sensitivity: within each sample, residualize the program score on
`log1p(total_counts)`, `log1p(n_genes_by_counts)` and `pct_counts_mt` using ordinary
least squares, then repeat the sample-level contrast. This is a sensitivity analysis,
not a replacement endpoint.

### Layer 2: sample-level myeloid aggregate

Within each scRNA dataset, compute Spearman correlation between aggregated SIGLEC15
and each program score. Also report a pooled within-dataset rank association after
removing dataset-specific rank means. The pooled 95% CI and two-sided P value use a
dataset-stratified patient bootstrap/permutation. Dataset-specific estimates remain
descriptive because series sizes are 5--9 samples.

### Layer 3: bulk tumours

Within each cohort, compute Spearman rho with a 95% percentile bootstrap CI (10,000
replicates) and a two-sided permutation P value. Combine Fisher-z transformed cohort
estimates using REML random effects and Hartung--Knapp confidence limits; report
heterogeneity and prediction intervals. With only three cohorts, meta-analytic
intervals are expected to be imprecise. No equivalence or non-inferiority conclusion
will be made from a point estimate or from the +/-0.20 descriptive region. An
equivalence statement would require the relevant 95% interval to lie wholly within
that region, which is not assumed.

### Layer 4: hierarchical spatial validation

For every spatial capture with at least 10 myeloid-labelled spots, at least two
SIGLEC15-positive spots and at least two negative spots, estimate the within-capture
program-score difference. Average eligible capture effects within each case, then
average the three case effects. Report the case effects and a t-based 95% interval
with 2 degrees of freedom. This is exploratory because only three cases are
available; spots and captures are not treated as independent patients. A sensitivity
analysis repeats the comparison after depth residualization within capture.

## Multiplicity and interpretation

- The ACP5/MMP9 core is the sole primary program at each layer.
- Within each layer, the two secondary program P values are Holm-adjusted.
- IL1B and IL6 are descriptive gene-level outputs.
- Evidence for a cell-restricted state requires a positive primary cell-level
  effect with a 95% CI excluding zero and directionally concordant secondary panels.
- Failure to reproduce at sample or bulk level requires imprecise/null association
  estimates to be described as such; it does not prove biological absence or
  equivalence.
- Spatial evidence is supportive/exploratory regardless of P value because n=3
  cases.
- RNA results cannot establish SIGLEC15 protein abundance, target engagement or
  treatment response.

## Outputs

All tables are CSV, the complete machine-readable summary is JSON, and the principal
figure is exported as vector PDF and 300-dpi PNG. The analysis script records input
paths, SHA-256 hashes, software versions, frozen gene coverage and all exclusions.
