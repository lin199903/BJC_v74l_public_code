# v1.0 — analysis code

Analysis code, frozen analysis rules and software-environment records for:

> Cross-scale evaluation of SIGLEC15 RNA for measurement of an ACP5/MMP9-high myeloid state in pancreatic ductal adenocarcinoma
> Junyi Lin, Yan Huang, Molong Hou

## What is in this release

| Path | Contents |
|---|---|
| `analysis/scripts/` | Primary cross-scale analysis and prespecified sensitivity analyses; coupling audits; figure-generation code. |
| `analysis/v72_burden/` | v72 tumour-level burden and sparsity extension, with its frozen estimand contract. |
| `analysis/V04_V05_state_restricted/` | State-restricted analyses and their claim-ceiling audit. |
| `analysis/phaseB_GSE205154/`, `phaseB_GSE21501/` | Bulk extension cohorts, with each dataset's frozen cohort definition. |
| `analysis/phaseB_GSE278688/`, `phaseB_GSE278689/` | Additional cellular sensitivity datasets, analysed separately: GSE278688 whole-cell single-cell RNA-seq, GSE278689 single-nucleus RNA-seq. Independent verification scripts included. |
| `analysis/extension_meta/scripts/` | Four-cohort and five-cohort random-effects meta-analysis and its `metafor` cross-check. |
| `paired_representation/` | V46 paired-representation extension: gate scripts and gate report. |
| `external_single_cell/` | External single-cell assessment: frozen analysis rules, eligibility and independent-check scripts, input inventory. |

## What is deliberately not in this release

The machine-readable result tables are **not** duplicated here — they are published with the paper as Additional files 2–4. This deposit is the code, together with the documents that fix the analysis rules, the run and provenance logs, the environment records, and the two small errata tables that the reading guide refers to.

No primary data are redistributed. TCGA, GTEx, CPTAC/LinkedOmics, PRIDE PXD048644 and the listed GEO series remain under the terms of their own repositories; accessions and inclusion rules are in the manuscript and its Supplementary Table S14.

## Where to start

- `README.md` — what is here, the artifact map, and the environment variables
- `ENVIRONMENT.md` — recorded original environment versus reproduction requirements
- `analysis/CURRENT_READING_GUIDE.md` — reporting corrections to apply when reading historical status headings

## Before publishing

Two fields describe the release itself and are only correct once it exists:

1. `CITATION.cff` → `date-released`. It currently reads `2026-10-14`, the date this candidate was assembled. If the GitHub Release is tagged on a later date, update this field to the tag date before releasing, so that the citation metadata and the archive agree.
2. `.zenodo.json` → `version`. It reads `v1.0` and must match the Git tag exactly. Zenodo's GitHub integration prefers `.zenodo.json` over `CITATION.cff` when both are present, so this is the file that decides the deposited version string.

`release.ps1` refuses to tag if `verify_manifest.py` fails. It does not check these two fields; reconcile them by hand immediately before running it.

## Integrity

`SHA256_MANIFEST.txt` locks every file in this release. Verify with:

```bash
python verify_manifest.py
```

## Licence

MIT — see `LICENSE`. Third-party data are not covered.

## Archived copy

This release is intended for archival by Zenodo, which mints the DOI cited in the manuscript's Code availability statement. The DOI is issued at publication; until then, cite the repository URL.
