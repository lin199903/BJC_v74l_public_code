# v46 paired-representation analysis extension

This archive contains the analyses added after the v45 production freeze. It
does not redistribute the public expression matrices.

## Same-sample RNA summaries

The cellular reference is the within-sample difference in the ACP5/MMP9 score
between SIGLEC15 RNA-detectable and non-detectable tumour-myeloid cells. Four
SIGLEC15 RNA summaries from the same samples are compared with that reference:
the myeloid aggregate, detectable-cell fraction, mean expression among
detectable cells and an all-cell tumour aggregate.

Run order:

1. `gate1_partA.py`
2. `gate1_partB.py`
3. `gate1_partB2.py`
4. `gate1_partC.py`

The final machine-readable result is
`Gate1/_gate1_final_summary.json`. Input files are public GEO deposits or
files documented in Additional file 2; their archived hashes are checked by the
scripts.

## Paired CPTAC-PDAC RNA and protein state

`gate2_p2.py` records protein availability in the paired LinkedOmics
CPTAC-PDAC matrices. `gate2_v46_corrected.py` calculates the final v46
associations, bootstrap intervals, score-construction sensitivities and positive
controls. For composition adjustment, the final script rank-transforms all
variables and residualises both SIGLEC15 RNA and the protein-state score on the
same covariates before correlating the residuals.

Required local inputs:

- `data/analysis/proteogenomic/rna_tumor.cct`
- `data/analysis/proteogenomic/proteome_tumor.cct`

Their byte counts and SHA-256 hashes are recorded in
`Gate2/_gate2_v46_corrected_summary.json`. The matrices can be obtained from
the LinkedOmics CPTAC-PDAC release cited in the manuscript.

## Reproducibility

Python package versions are stored in `SOFTWARE_VERSIONS.txt`. All bootstrap
seeds are fixed in the scripts. `SHA256SUMS.txt` inventories every other
payload file in this archive.
