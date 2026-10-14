# results_summary.md — GSE21501, Phase B step 3c-2

**Cohort:** GSE21501 (GPL4133, Agilent-014850 Feature Number version), 132 deposited arrays.
**Frozen manifest:** `sample_manifest.csv`, sha256 `3e11938a…f3ac4e4` — **re-hashed at the
start and again at the end of this run; it still matched.**

> **The GSE21501 analysis was specimen-level because patient identifiers were not deposited; the number of independent patients could therefore not be established exactly.**

The unit throughout is the **eligible tumour array/specimen**, never "patients". The
figures below are therefore **129 eligible tumour arrays/specimens, representing ≤129
patients** (and, for the subset, **≤102 patients**). No exact patient count is asserted
anywhere.

---

## 1. Primary estimate (AD-1: expression-only set, n = 129)

| | |
|---|---|
| **ρ** | **+0.1712** |
| **n** | **129** eligible tumour arrays/specimens (129 distinct specimen keys; ≤129 patients) |
| **95 % CI** | **−0.0128 (−0.0128237521) to +0.3405 (+0.3405308288)** |
| **CI method** | percentile bootstrap over arrays, **10 000 resamples**, Spearman recomputed per resample (the manuscript's own `spearman_bootstrap`) |
| **seed** | **15072026** (frozen), consumed in the manuscript's call order |
| asymptotic *p* | 0.0524 |
| permutation *p* | 0.0527 (10 000 draws) |

The permutation *p* is a simulation quantity, not an exact statistic: this 0.0527 is drawn
from the generator left after the bootstrap has consumed its draws in the manuscript's call
order, and a fresh generator at the same seed 15072026 gives 0.0507. Both are reported in
`output/run_log.txt` and in the verification section. Neither crosses 0.05, so nothing in
this summary depends on which draw is quoted. ρ and the interval are unaffected — ρ carries
no generator dependence at all.

The estimate is the Spearman correlation between the deposited SIGLEC15 value and
`Derosiers_ACP5_MMP9_core`, the mean of the within-cohort z-scores of ACP5 and MMP9,
built by the frozen pipeline's own `zscore_columns`/`score_program`. SIGLEC15 is
`log2(Cy5/Cy3)` as supplied, with no further transform (AD-2); because ρ uses ranks, AD-2
can move only the score side.

## 2. Prespecified sensitivity: clinically annotated subset (AD-6, n = 102)

| | |
|---|---|
| **ρ** | **+0.1478** |
| **n** | **102** eligible tumour arrays/specimens (≤102 patients) |
| **95 % CI** | −0.0577 to +0.3441 |
| permutation *p* | 0.1445 |

**The 30 Nebraska/UNC arrays are flagged `clinical annotation unavailable`, not excluded.**
They are in the primary and out of this subset. The two estimates differ by **|Δρ| =
0.0234** and their intervals overlap across almost their whole length. This is a
**metadata-availability sensitivity**: both numbers are retained and reported, and neither
was chosen on the basis of being larger.

## 3. Sensitivity spread — every prespecified version

| id | what it varies | n | ρ | 95 % CI |
|---|---|---|---|---|
| primary (AD-1) | expression-only, whole deposit after DUP-2 | 129 | **+0.1712** | −0.0128, +0.3405 |
| S-CLIN (AD-6) | clinically annotated subset | 102 | +0.1478 | −0.0577, +0.3441 |
| S2 (AD-4 as ruled) | drop `P188T IPMN` (not PDAC); `P271N` retained | 101 | +0.1703 | −0.0302, +0.3604 |
| S2-legacy (step-3b text) | additionally drop `P271N` — comparison only | 100 | +0.1462 | −0.0514, +0.3397 |
| S3 | drop all six repeat-marked arrays | 96 | +0.1189 | −0.0868, +0.3137 |
| S4 | duplicate pairs averaged at expression level | 129 | +0.1712 | −0.0094, +0.3439 |
| S5 | per-MMP9-feature, 10 features (**range, no feature selected**) | 102 | +0.1403 to +0.1602 (median +0.1480) | −0.0619, +0.3549 across the ten features |
| S-rep-1 | rank-based ACP5/MMP9 score, re-standardised (frozen S8) | 129 | +0.1683 | −0.0104, +0.3387 |
| S-rep-1b | rank-based, raw 1..n ranks averaged | 129 | +0.1698 | −0.0107, +0.3368 |

**No "no-log-transform" alternative exists for this deposit.** GSE21501 supplies a single
representation — the Lowess-normalised `log2(Cy5/Cy3)` ratio — so there is no untransformed
quantity to re-run on without raw feature intensities, which are not deposited. This
absence is a property of the deposit, not a choice made here.

## 4. Stable versus unstable

**Stable.** All **18** estimates computed for this cohort lie between **+0.1189 and
+0.1712** — a total spread of 0.052 — and **every one of their 95 % intervals contains
zero**. De-duplication is inert in the primary and in S4: the three duplicate groups are
byte-identical technical replicates, so DUP-2 and pair-averaging return the same ρ. The
rank-based score alternatives move ρ by at most 0.003, and S5's ten single-feature scores
span only 0.020. The estimate does not depend on the standardisation device, the duplicate
rule, the representation, or the MMP9 collapse.

**Also stable, and the sharper fact.** No version reaches significance at the 0.05 level:
the smallest permutation *p* across all 18 is 0.0527. The primary interval clears zero by
only 0.013 at its lower limit. The two *n* values are therefore not a choice between a
positive and a null result; they are two renderings of the same indeterminate one.

**Unstable / worth flagging.** S3 is the one sensitivity that moves ρ appreciably, from
+0.1712 to +0.1189, i.e. by 0.052 — the six submitter-flagged re-runs are not a random
subset and their removal weakens the association by about a third of its own size. It is
reported, not preferred.

**Not part of the estimate but relevant to reading it.** The composite's component
associations are heterogeneous: over the 129, ρ(SIGLEC15, ACP5) = −0.0390 (*p* = 0.66)
while ρ(SIGLEC15, MMP9) = +0.3404 (*p* = 7.9 × 10⁻⁵). The composite sits where an unweighted
average of a null and a modest positive association would put it. These are descriptive
component values recorded in `output/descriptive_components.csv`; they are **not**
alternatives to the primary and were not used to select it.

## 5. Not determinable from this deposit

- **Any exact patient count.** GEO exposes no patient identifier, no `!Sample_relation` and
  no BioSample accession for GSE21501, so nothing can convert 129 specimens into a patient
  number. All patient figures are upper bounds (AD-3).
- **Whether `P271N` (`GSM536488`) is normal tissue.** Its trailing `N` is undocumented and
  the frozen manifest does not establish it as non-PDAC or normal, so under AD-4 it is
  **retained**; S2-legacy shows what excluding it costs (ρ +0.1703 → +0.1462). Which of the
  two readings is right is not decidable from public data.
- **Whether `P188T` (`GSM536901`) is invasive PDAC or a pure precursor lesion.** It is
  excluded from S2 as IPMN, and for that reason alone — **not** for missing clinical data.
- **Any cross-cohort absolute expression comparison.** Forbidden by AD-2 and not performed:
  the deposit supplies a ratio against a common reference, not an absolute level.
- **Any pooled or meta-analytic estimate.** Out of scope for this step. **GSE205154 is not
  touched here and nothing is pooled with it.** The TCGA-PAAD, GSE57495 and GSE85916
  estimates are unaffected by this run.

## 6. Verification in one line

The manuscript's **own archived scoring module** was imported and executed on this cohort's
data: `zscore_columns`, `score_program`, `spearman_bootstrap` and `permutation_p` all agree
cell-for-cell and to the last digit (**9/9 checks**), four of them capable of failing.
Twelve further controls were run on this run's own pipeline (**12/12 pass**), including a
circular-rotation null and a label-scramble null that would both have failed had the
association been manufactured from array order, and a seed control that fails if the
resampling is not reproducible. Full detail in `output/verification_and_controls.md`.

---

*Artefacts: `output/primary_result.csv`, `output/sensitivity_results.csv`,
`output/estimate_table_full.csv`, `output/verification_controls.csv`,
`output/manuscript_verification.csv`, `output/descriptive_components.csv`,
`output/run_log.txt`, `output/sessionInfo.txt`, `output/seed_registry.json`. Scripts in
`scripts/`, each paired with its seed in `output/seed_registry.json`.*
