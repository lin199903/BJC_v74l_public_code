# Additional file 4. External single-cell assessment

This archive supports the external single-cell assessment reported in the manuscript
(Results 3.2, Discussion; Additional file 1, §17 and Table S29).

## Contents

    FROZEN_ANALYSIS_RULES.md     analysis rules and reading criteria, written before any
                                 expression matrix was downloaded or any effect computed
    EXECUTION_REPORT.md          what was run, what was found, and what the result does not
                                 establish
    results/E2B_RESULT.json      per-sample myeloid cell counts, detectable-cell counts,
                                 contrasts, stratum summaries and the bootstrap interval
    eligibility/                 eligibility check of the candidate external datasets,
                                 including the two that failed on gene-panel coverage
    scripts/                     the analysis scripts, unmodified
    logs/                        execution log
    INPUT_INVENTORY.tsv          accession, size and SHA-256 of every downloaded input file
    ENVIRONMENT.txt              software versions and platform
    SHA256SUMS.txt               integrity records for every file in this archive

## Data

The input expression matrices are public and are not redistributed here. They are the
filtered 10x matrices for GSE217845, retrieved per sample from the NCBI Gene Expression
Omnibus. `INPUT_INVENTORY.tsv` lists every file that was downloaded, its size and its
SHA-256, so the same inputs can be retrieved and verified.

## Scope

The assessment estimates the cellular ACP5/MMP9 contrast in a separate public PDAC
single-cell dataset. It is an independent data-source check of the cellular association. It
is not a patient-level proxy qualification, it does not include ambient-RNA correction, and
it does not provide spatial or protein localisation.
