# RUN PROVENANCE -- Phase B step 1, GSE278688 (frozen production run)

**Run started (local, Asia/Shanghai):** 2026-09-28 00:45
**Operator:** delegated subagent (parent `session-bbd7f29d-e275-4015-ae46-70635913f222`)
**Scope:** GSE278688 ONLY. No other dataset touched.
**Output root:** `E:\Tex\HMlong\BJC_EXT_v1_20260927\phaseB_GSE278688\`
**Frozen input dir (read-only by discipline):** `E:\Tex\HMlong\BJC_EXT_v1_20260927\annotations\`

## Environment

| item | value |
|---|---|
| OS | Windows 11 x64 build 26200 |
| shell | Windows PowerShell 5.1.26100.9444 (NOTE: `pwsh` 7 is NOT installed; `ForEach-Object -Parallel` unavailable) |
| R | 4.5.1 (2025-06-13 ucrt), x86_64-w64-mingw32 |
| Python | 3.13.0 |
| CPU logical cores | 24 |
| RAM | 31.36 GB |
| Free disk E: | 398.6 GB |

## Source files (GEO FTP over HTTPS)

Root: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn/GSE278688/suppl/`

| file | bytes | sha256 |
|---|---|---|
| GSE278688_sc_barcodes.tsv.gz | 911,362 | (see hashes.txt) |
| GSE278688_sc_features.tsv.gz | 221,703 | (see hashes.txt) |
| GSE278688_sc_metadata.csv.gz | 1,166,884 | (see hashes.txt) |
| GSE278688_sc_matrix.mtx.gz | 1,093,755,688 | (see hashes.txt) |

**There is NO `RAW.tar` for GSE278688** (confirmed 404). `filelist.txt` also does not exist (404).
The deposit is a **single combined matrix** + barcodes + features + metadata, i.e. one object
containing all 29 samples, not 29 per-sample matrices.

## Authoritative SOFT metadata (retrieved 2026-09-28)

Endpoint: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278688&targ=self&form=text&view=brief`
Saved at `_probe/GSE278688_soft_brief.txt`.

- `!Series_title` = Integrated Single-Cell and Spatial Transcriptomics Uncover Distinct Cellular Subtypes Involved in Neural Invasion in PDAC [scRNA-Seq]
- `!Series_platform_id` = GPL24676
- `!Series_pubmed_id` = 40680743 ; DOI 10.1016/j.ccell.2025.06.020 (Cancer Cell)
- `!Series_overall_design` = "We generated scRNA-seq and snRNA-seq data from **25 patients** ... including tumor and normal tissue."
- `!Series_sample_id` = **29 GSMs** (GSM8582554–GSM8582582)
- `!Series_relation` = **SubSeries of: GSE278694**
- Raw sequencing data deposited at GSA **HRA008568** (not at GEO)

### Access note (reproducible recipe)
`web_fetch` cannot read GEO SOFT, and GEO HTML pages are reCAPTCHA-protected.
Working path: PowerShell `Invoke-WebRequest`. **Trap:** in this shell `$r.Content` comes back as a
space-separated list of decimal byte codes; it must be reassembled with
`[System.Text.Encoding]::UTF8.GetString([byte[]]$codes)`.

## Files written by this run

| dir | contents |
|---|---|
| `data/` | downloaded deposit files (unmodified) |
| `scripts/` | every script, one per result file |
| `logs/` | one log per script run, with timestamps |
| `output/` | result tables |
| `provenance/` | hashes, seeds, this file |
