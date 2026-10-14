# failed_or_excluded.md — GSE205154

**Step:** Phase B step 3a — manifest construction and cohort-definition freeze (bulk RNA-seq)
**Dataset:** GSE205154, *Bulk RNA-Seq gene expression data for 289 primary and metastatic
pancreatic ductal adenocarcinoma (PDAC) tumor samples from FFPE blocks* (OHSU;
SuperSeries GSE281129; PMID 39789181, *Nature Cancer* 2025, doi 10.1038/s43018-024-00881-3)
**Deposited samples:** 289
**Companion files:** `sample_manifest.csv`, `frozen_cohort_definition.md`, `run_log.txt`,
`sessionInfo.txt`
**Manifest SHA-256:** `f76b81eceb0557cf4c1fee0138e9137a84508369f47ad163255a096af78c708b`

> **Governing rule observed in this step: lock the cohort BEFORE looking at the effect.**
> No correlation between SIGLEC15 and any ACP5/MMP9 programme score was computed,
> inspected or reported in this run. No per-sample expression value of SIGLEC15, ACP5 or
> MMP9 was loaded.

---

## 0. How to read this file

Every one of the 289 deposited samples carries an explicit pass/fail on both eligibility
columns in `sample_manifest.csv`, each with a written `exclusion_reason`. There are no
silent skips, and there is no "partially processed" state: a sample is either `yes` or
`no` on `eligible_for_main_meta` and on `eligible_for_sensitivity`, and every `no` has a
concrete reason recorded in the same row.

This file lists the exclusions grouped by rule, then records everything attempted during
this step that failed, and everything that could not be determined.

---

## 1. Excluded from the MAIN meta-analysis — 73 of 289

### 1.1 Metastatic lesions — 71 samples

**Reason (uniform, applied by the frozen specimen rule).** The main meta-analysis is
restricted to primary tumours; the frozen protocol places metastatic lesions in a
separate prespecified sensitivity analysis and forbids pooling the two compartments into
one dataset-level estimate. Each of these 71 samples is therefore `eligible_for_main_meta = no`
and `eligible_for_sensitivity = yes`.

| GSM | title | site | patient key |
|---|---|---|---|
| GSM6206008 | ST-00006650-M | Lung | ST-00006650 |
| GSM6206018 | ST-00007377-M | Peritoneum | ST-00007377 |
| GSM6206021 | ST-00007402-M | Mesocolon | ST-00007402 |
| GSM6206026 | ST-00007369-M | Liver | ST-00007369 |
| GSM6206030 | ST-00007223-M | Liver | ST-00007223 |
| GSM6206033 | ST-00007381-M | Peritoneum | ST-00007381 |
| GSM6206036 | ST-00007380-M | Gallbladder | ST-00007380 |
| GSM6206043 | ST-00021353-M | Liver | ST-00021353 |
| GSM6206044 | ST-00007382-M | Peritoneum | ST-00007382 |
| GSM6206048 | ST-00007389-M | Liver | ST-00007389 |
| GSM6206049 | ST-00007288-M | Peritoneum | ST-00007288 |
| GSM6206053 | ST-00007595-M | Liver | ST-00007595 |
| GSM6206055 | ST-00006597-M | Falciform Ligament | ST-00006597 |
| GSM6206057 | ST-00007594-M | Liver | ST-00007594 |
| GSM6206058 | ST-00007590-M | Liver | ST-00007590 |
| GSM6206059 | ST-00007307-M | Jejunum | ST-00007307 |
| GSM6206063 | ST-00009788-M | Liver | ST-00009788 |
| GSM6206065 | ST-00009790-M | Peritoneum | ST-00009790 |
| GSM6206067 | ST-00006291-M | Retroperitoneum | ST-00006291 |
| GSM6206074 | ST-00011356-M | Liver | ST-00011356 |
| GSM6206078 | ST-00012967-M | Liver | ST-00012967 |
| GSM6206080 | ST-00010984-M | Peritoneum | ST-00010984 |
| GSM6206082 | ST-00012975-M | Liver | ST-00012975 |
| GSM6206086 | ST-00020045-M | Liver | ST-00020045 |
| GSM6206089 | ST-00020057-M | Peritoneum | ST-00020057 |
| GSM6206096 | ST-00020064-M | **Pancreas** | ST-00020064 |
| GSM6206104 | ST-00020070-M | Liver | ST-00020070 |
| GSM6206106 | ST-00020069-M | Peritoneum | ST-00020069 |
| GSM6206108 | ST-00020072-M | Omentum | ST-00020072 |
| GSM6206110 | ST-00013526-M | Lung | ST-00013526 |
| GSM6206120 | ST-00014524-M | Liver | ST-00014524 |
| GSM6206127 | ST-00020104-M | Peritoneum | ST-00020104 |
| GSM6206132 | ST-00020105-M | Liver | ST-00020105 |
| GSM6206136 | ST-00020106-M | Liver | ST-00020106 |
| GSM6206138 | ST-00020112-M | Liver | ST-00020112 |
| GSM6206141 | ST-00015839-M | Ovary | ST-00015839 |
| GSM6206144 | ST-00020092-M | **NA** | ST-00020092 |
| GSM6206145 | ST-00020119-M | Lymph Node | ST-00020119 |
| GSM6206147 | ST-00020118-M | Liver | ST-00020118 |
| GSM6206156 | ST-00020126-M | Peritoneum | ST-00020126 |
| GSM6206157 | ST-00020131-M | Peritoneum | ST-00020131 |
| GSM6206162 | ST-00017003-M | Liver | ST-00017003 |
| GSM6206173 | ST-00020089-M | Liver | ST-00020089 |
| GSM6206175 | ST-00017804-M | Lung | ST-00017804 |
| GSM6206180 | ST-00020149-M | Peritoneum | ST-00020149 |
| GSM6206182 | ST-00020144-M | Lung | ST-00020144 |
| GSM6206188 | ST-00017838-M | Liver | ST-00017838 |
| GSM6206190 | ST-00020196-M | Falciform Ligament | ST-00020196 |
| GSM6206192 | ST-00020154-M | Liver | ST-00020154 |
| GSM6206200 | ST-00017915-M | Lung | ST-00017915 |
| GSM6206205 | ST-00020168-M2 | Liver | ST-00020168 |
| GSM6206208 | ST-00018232-M | Lung | ST-00018232 |
| GSM6206209 | ST-00020164-M | Liver | ST-00020164 |
| GSM6206216 | ST-00020213-M | Liver | ST-00020213 |
| GSM6206225 | ST-00020190-M | Liver | ST-00020190 |
| GSM6206226 | ST-00020185-M | Omentum | ST-00020185 |
| GSM6206232 | ST-00004898-M | Lung | ST-00004898 |
| GSM6206234 | ST-00019607-M | Peritoneum | ST-00019607 |
| GSM6206237 | ST-00021359-M | Perigastric Nodule | ST-00021359 |
| GSM6206239 | ST-00020203-M | Liver | ST-00020203 |
| GSM6206248 | ST-00020308-M | Peritoneum | ST-00020308 |
| GSM6206255 | ST-00019043-M | LymphNode | ST-00019043 |
| GSM6206262 | ST-00020232-M | Retroperitoneum | ST-00020232 |
| GSM6206265 | ST-00019171-M | Duodenum | ST-00019171 |
| GSM6206266 | ST-00020251-M | Peritoneum | ST-00020251 |
| GSM6206267 | ST-00020323-M | Omentum | ST-00020323 |
| GSM6206268 | ST-00020281-M | Liver | ST-00020281 |
| GSM6206273 | ST-00019429-M | Duodenum | ST-00019429 |
| GSM6206274 | ST-00019601-M | Duodenum | ST-00019601 |
| GSM6206289 | ST-00020944-M | Liver | ST-00020944 |
| GSM6206290 | ST-00020824-M | Liver | ST-00020824 |

**Two rows in this table need a note, and neither is a defect:**

* **GSM6206096 `ST-00020064-M`** — a metastasis *located in the pancreas*. Three
  independent fields agree that it is metastatic (`tumor type: Met`,
  `description: Metastisis to the Pancreas`, title suffix `M`). A pancreatic
  metastasis is not a primary tumour, so it stays in the metastatic compartment.
  This is exactly the case where `tissue:` alone would have mis-classified the sample;
  the frozen rule keys on `tumor type`, not on `tissue`.
* **GSM6206144 `ST-00020092-M`** — `tissue: NA`, `description: Metastisis to the NA`.
  The metastatic *site* is not stated, but the *compartment* is unambiguous from
  `tumor type: Met`. The sample is therefore retained in the sensitivity analysis, and
  the unstated site is recorded in its `exclusion_reason`. Only the site annotation is
  missing; nothing about eligibility depends on it.

### 1.2 Second primary specimen from an already-represented patient — 2 samples

**Reason.** The frozen specimen rule is "one independent tumour/patient" and requires a
pre-defined aggregation-or-selection rule for duplicate specimens. Both rows below are a
*second* primary specimen from a patient who already contributes a first one; the frozen
rule retains the base-label specimen.

| GSM | title | patient key | reason |
|---|---|---|---|
| GSM6206193 | ST-00017945-T2 | ST-00017945 | second primary specimen from this patient; the frozen one-specimen-per-patient rule retained `ST-00017945-T` (GSM6206194) |
| GSM6206245 | ST-00019598-T2 | ST-00019598 | second primary specimen from this patient; the frozen one-specimen-per-patient rule retained `ST-00019598-T` (GSM6206244) |

---

## 2. Excluded from the prespecified SENSITIVITY analysis — 218 of 289

**Reason (uniform).** The sensitivity analysis is the metastatic compartment only, by the
frozen rule that the two compartments are never pooled into one dataset-level estimate.
All 218 primary specimens are `eligible_for_sensitivity = no`, including the two
duplicate specimens in §1.2 (which are excluded from both analyses).

The 218 GSMs, for the avoidance of any silent omission:

```
GSM6206002 GSM6206003 GSM6206004 GSM6206005 GSM6206006 GSM6206007 GSM6206009 GSM6206010
GSM6206011 GSM6206012 GSM6206013 GSM6206014 GSM6206015 GSM6206016 GSM6206017 GSM6206019
GSM6206020 GSM6206022 GSM6206023 GSM6206024 GSM6206025 GSM6206027 GSM6206028 GSM6206029
GSM6206031 GSM6206032 GSM6206034 GSM6206035 GSM6206037 GSM6206038 GSM6206039 GSM6206040
GSM6206041 GSM6206042 GSM6206045 GSM6206046 GSM6206047 GSM6206050 GSM6206051 GSM6206052
GSM6206054 GSM6206056 GSM6206060 GSM6206061 GSM6206062 GSM6206064 GSM6206066 GSM6206068
GSM6206069 GSM6206070 GSM6206071 GSM6206072 GSM6206073 GSM6206075 GSM6206076 GSM6206077
GSM6206079 GSM6206081 GSM6206083 GSM6206084 GSM6206085 GSM6206087 GSM6206088 GSM6206090
GSM6206091 GSM6206092 GSM6206093 GSM6206094 GSM6206095 GSM6206097 GSM6206098 GSM6206099
GSM6206100 GSM6206101 GSM6206102 GSM6206103 GSM6206105 GSM6206107 GSM6206109 GSM6206111
GSM6206112 GSM6206113 GSM6206114 GSM6206115 GSM6206116 GSM6206117 GSM6206118 GSM6206119
GSM6206121 GSM6206122 GSM6206123 GSM6206124 GSM6206125 GSM6206126 GSM6206128 GSM6206129
GSM6206130 GSM6206131 GSM6206133 GSM6206134 GSM6206135 GSM6206137 GSM6206139 GSM6206140
GSM6206142 GSM6206143 GSM6206146 GSM6206148 GSM6206149 GSM6206150 GSM6206151 GSM6206152
GSM6206153 GSM6206154 GSM6206155 GSM6206158 GSM6206159 GSM6206160 GSM6206161 GSM6206163
GSM6206164 GSM6206165 GSM6206166 GSM6206167 GSM6206168 GSM6206169 GSM6206170 GSM6206171
GSM6206172 GSM6206174 GSM6206176 GSM6206177 GSM6206178 GSM6206179 GSM6206181 GSM6206183
GSM6206184 GSM6206185 GSM6206186 GSM6206187 GSM6206189 GSM6206191 GSM6206193 GSM6206194
GSM6206195 GSM6206196 GSM6206197 GSM6206198 GSM6206199 GSM6206201 GSM6206202 GSM6206203
GSM6206204 GSM6206206 GSM6206207 GSM6206210 GSM6206211 GSM6206212 GSM6206213 GSM6206214
GSM6206215 GSM6206217 GSM6206218 GSM6206219 GSM6206220 GSM6206221 GSM6206222 GSM6206223
GSM6206224 GSM6206227 GSM6206228 GSM6206229 GSM6206230 GSM6206231 GSM6206233 GSM6206235
GSM6206236 GSM6206238 GSM6206240 GSM6206241 GSM6206242 GSM6206243 GSM6206244 GSM6206245
GSM6206246 GSM6206247 GSM6206249 GSM6206250 GSM6206251 GSM6206252 GSM6206253 GSM6206254
GSM6206256 GSM6206257 GSM6206258 GSM6206259 GSM6206260 GSM6206261 GSM6206263 GSM6206264
GSM6206269 GSM6206270 GSM6206271 GSM6206272 GSM6206275 GSM6206276 GSM6206277 GSM6206278
GSM6206279 GSM6206280 GSM6206281 GSM6206282 GSM6206283 GSM6206284 GSM6206285 GSM6206286
GSM6206287 GSM6206288
```

---

## 3. Excluded from BOTH analyses — 2 samples

Exactly two samples are ineligible for everything: `GSM6206193` (`ST-00017945-T2`) and
`GSM6206245` (`ST-00019598-T2`), both by the duplicate-specimen rule (§1.2).

## 4. No sample failed a gene-usability condition

All three required genes are present exactly once in the deposited gene-level matrices,
as `protein_coding` rows with the expected Ensembl identifiers, and all three carry
`TMM_pass_filter = TRUE` in the authors' own deposited filter column:

| gene | Ensembl | biotype | row index (0-based) | TMM_pass_filter |
|---|---|---|---|---|
| SIGLEC15 | ENSG00000197046 | protein_coding | 17332 | TRUE |
| ACP5 | ENSG00000102575 | protein_coding | 2686 | TRUE |
| MMP9 | ENSG00000100985 | protein_coding | 2418 | TRUE |

Consequently **zero** samples were excluded for gene usability, and the arithmetic is
`289 − 71 (metastatic) − 2 (duplicate primary) = 216` for the main analysis and
`289 − 218 (primary) = 71` for the sensitivity analysis.

---

## 5. Attempts that failed, and what was tried

Everything below is also recorded with timestamps and exit codes in `run_log.txt`.

### 5.1 GEO web access

* **GEO HTML pages are reCAPTCHA-protected.** Attempted route abandoned; all metadata
  was taken from the **NCBI FTP over HTTPS** mirror instead
  (`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205154/...`), which returns the
  same series-matrix and family-SOFT records without a challenge.
* **GEO's SOFT CGI endpoint is rejected by the `web_fetch` tool.** Not retried; the FTP
  `soft/GSE205154_family.soft.gz` copy was used.

### 5.2 Download stalls — two, both recovered without restarting from zero

* **Attempt 1** (`scripts/03_scan_features.py`, single-connection read loop). The read
  stalled at exactly **8,388,608 bytes** with the socket open and no exception: the
  process had consumed 0.33 s of CPU over roughly five minutes. There was no per-read
  timeout, so the loop could not detect the stall.
  **Fix:** a dedicated resumable downloader, `scripts/03a_download_matrices.py`, with a
  45 s socket timeout, range-based resume from the last good byte, and verification of
  the final size against the server's `Content-Length`.
* **Attempt 2** (`scripts/03a_download_matrices.py`, still single-connection). It resumed
  from the existing 9,437,184 bytes and ran at roughly **18 KB/s**, and stalled a second
  time around 27 MB before recovering. Two concurrent single-connection downloads shared
  roughly 30 KB/s, leaving an estimated 33 minutes for the TMM and TPM matrices.
  **Fix:** `scripts/03b_parallel_download.py` — each file split into 8 byte ranges
  fetched on 8 concurrent connections, each range written to its own resumable part file,
  every part length-checked, parts concatenated and the assembled size checked against
  `Content-Length`. Throughput rose to roughly 550 KB/s. The partially downloaded
  single-connection copies were discarded, not appended to.
  Results: TMM `40,598,676` bytes, sha256 `9f362bb06995e28d3f0f2297f89ddf0f1190235c76a1a8cfb5c199d8211895c3`;
  Counts `30,091,494` bytes, sha256 `4c7004cfe15a94d868e24cf651f0e5d05779709bf233a74b48b084efdf238fdc`.
* **Guard added as a result.** `scripts/03_scan_features.py` now refuses to scan any file
  whose size differs from the server's `Content-Length`, because a truncated gzip stream
  would silently under-report the feature set and could have produced a false
  "gene absent" verdict.

### 5.3 Script defects found and corrected during this step

These are recorded because a silent, unrecorded correction is exactly the failure mode
this project's audit discipline exists to prevent.

1. **`04_parse_soft_samples.py`, run 1 — BioSample identifier stored as a URL.**
   The `BioSample:` relation was stored whole
   (`https://www.ncbi.nlm.nih.gov/biosample/SAMN28776628`) instead of the bare accession,
   so the per-core join in `06_patient_identity_probe.py` matched nothing and printed
   `attrs={}` for every sample. The *aggregate* attribute census in that run was computed
   straight from the efetch XML and was still correct, so the run-1 headline conclusion
   ("no patient-key attribute") happened to be right — but the supporting detail was
   empty and could not be audited.
   **Fix:** store the bare accession; re-ran 04, 05 and 06. Run 2 shows the real
   attributes. *Note on run 2's output:* the line
   `duplicate cores sharing an attribute value` lists all 12 duplicate cores, and the
   per-core comparison shows shared values such as `Pancreas`, `Primary`,
   `scraped tumor`. **These are shared because both specimens in a core are primary
   pancreatic tumours — they are not evidence of a patient key.** The only two
   discriminating fields the record exposes are `tumor type` and `tissue`, and both are
   *supposed* to agree within a patient.
2. **`06_patient_identity_probe.py`, run 1 — missing table guard on the SuperSeries
   matrix.** The GSE281129 series matrix was decompressed and split into lines without
   the data-block guard used everywhere else in this step.
   **Fix:** run 2 applies the same guard and reports how many data-block lines were
   dropped (1).
3. **`09_overlap_check.py`, run 1 — wrong separator for series-matrix records.**
   Series matrices are TAB-separated; only family-SOFT files use ` = `. The run-1 parser
   split on ` = `, so every series parsed as **0 samples** and the "no overlap" result was
   vacuously true.
   **Fix:** `split_row()` now handles both. Re-run.
4. **`09_overlap_check.py`, run 1 — wrong series-matrix filename for GSE155698.**
   GSE155698 has two platform-specific matrices
   (`GSE155698-GPL20301_...` and `GSE155698-GPL24676_...`), not `GSE155698_series_matrix.txt.gz`;
   the single-URL fetch returned HTTP 404.
   **Fix:** the series table now carries a list of candidate filenames, all of which are
   fetched and merged.
5. **`09_overlap_check.py`, run 2 — `!Series_sample_id` is space-packed.**
   In a GEO series matrix, all sample accessions sit in **one quoted, space-separated
   field**, so splitting on tab alone yields a single token: run 2 reported
   `samples in series : 1` for every dataset and the GSM-overlap test compared against a
   one-element set. The specimen-core and title tests were unaffected and already showed
   zero overlap, but the accession test was meaningless.
   **Fix:** `unquote_ids()` splits on tabs *and* whitespace. Run 3 gives the real series
   sizes (GSE155698 41, GSE156405 17, GSE212966 12, GSE217845 16, GSE282302 108,
   GSE57495 63, GSE85916 80) and confirms zero overlap on all three tests.
6. **`03_scan_features.py` — dangling call to a removed function.**
   After the download step was moved out to `03a`, the call site still invoked
   `download(...)`, raising `name 'download' is not defined`.
   **Fix:** the call site now invokes `ensure_complete(...)`, which verifies size against
   `Content-Length` and refuses partial files.

### 5.4 Tool-level failures

* The **first `09_overlap_check.py` invocation exceeded the 120 s command timeout**
  because the script buffered all output until completion and its 8 MB range fetches were
  competing with the matrix download. **Fix:** moved to a background job, reduced the
  fetch chunk to 512 KiB, and made the fetcher stop as soon as the `!Series_sample_id`
  and `!Sample_title` rows have been seen.
* **`Get-ChildItem` reported a stale size** for files another process had open for
  writing (it showed 9 MB while the file was actually 23 MB), which briefly looked like a
  download shrinking. **Fix:** true sizes were read with
  `[System.IO.File]::Open($path,'Open','Read','ReadWrite').Length`.

### 5.5 A conclusion carried in from the prior scan was WRONG, and is corrected here

The prior scan recorded **"no patient identifier exposed"**. That is true of the
GSE205154 series record and of all 289 of its BioSample records — but it is **not true of
the deposit as a whole**. The SuperSeries GSE281129 has a second child, **GSE281005**,
which publishes a characteristic literally named **`public patient_id` for all 505 of its
samples**, and shares the `ST-########` specimen-title namespace with GSE205154. The
mapping between `public patient_id` and specimen core is strictly 1:1 across 290 cores (no
core with two ids, no id with two cores).

**Consequence:** patient identity *can* be established for this cohort, so the effective
unit is the patient, not the sample — see `frozen_cohort_definition.md` §3. Had the prior
scan's belief been carried forward unverified, the manifest would have been built on the
sample as the unit and the two duplicate primary specimens (§1.2) would have entered the
main analysis as if they were independent tumours.

---

## 6. Could not be determined

1. **The meaning of the `F` title suffix.** Seven samples carry a suffix `F`
   (`ST-00021339-F`, `ST-00021340-F`, `ST-00020206-F`, `ST-00020208-F`, `ST-00020215-F`,
   `ST-00020926-F`, `ST-00021060-F`). All seven are concordantly annotated
   `tumor type: Primary`, `tissue: Pancreas`, `description: Primary pancreas tumor`. The
   deposit does not define what `F` denotes. **This does not affect eligibility**, because
   the frozen compartment rule reads the `tumor type` field and not the title suffix, and
   all three independent fields agree. Recorded so the vocabulary is not mistaken for
   something it may not be.
2. **Whether two *different* specimen cores belong to one patient.** For 83 of the 277
   cores, no GSE281005 counterpart exists, so no published `public patient_id` is
   available and no external attribute can corroborate the core. Across the 290 cores
   where both series overlap, core ↔ patient is exactly 1:1 with no exceptions, which is
   evidence against core-splitting but is not proof for the 83 unmatched cores. If a
   patient did contribute two differently-numbered cores, the deposit exposes no way to
   detect it.
3. **Patient-level overlap with GSE282302.** GSE282302 is a Visium series from the same
   institution (OHSU) and is used descriptively in the manuscript. Its sample titles are
   ROI identifiers (`C1_D10_ROI1_s1`) that share no namespace with `ST-########`, so
   `GSM` intersection is 0, specimen-core overlap is 0 and title overlap is 0 — but a
   patient-level overlap **cannot be excluded** from public metadata. It does not affect
   this cohort's arithmetic: GSE282302 contributes no bulk expression and is not one of
   the bulk meta-analysis cohorts.
4. **Patient-level overlap with TCGA-PAAD.** Not assessable: GSE205154 samples carry no
   TCGA barcode and there is no shared accession namespace. The cohorts are from different
   institutions (OHSU vs. the TCGA network).
5. **Raw sequence data.** Not obtained. The series states that raw data are deposited in
   dbGaP under `phs003597.v1.p1` (controlled access). Only the processed gene-level
   matrices were used, which is what the frozen protocol requires.
6. **Sex and age for 87 samples.** `Sex` and `age at_diagnosis` are exposed only through
   GSE281005, so they exist for the 202 samples whose core appears there and are absent
   for the other 87. Neither field enters any frozen eligibility rule.

---

## 7. Verification status of this step's artefacts

| check | result |
|---|---|
| deposited samples reconciled to GSMs | 289 / 289, all unique |
| BioSample accessions | 289 / 289, all unique |
| matrix sample columns vs SOFT titles | identical as a set for Counts, TMM (289 each) |
| compartment field populated for every sample | 289 / 289 |
| `Sample_description` vs `tumor type` conflicts | 0 |
| `tissue` vs `source_name` mismatches | 0 |
| excluded samples with a written reason | 289 / 289 rows carry an explicit reason |
| gene presence (SIGLEC15 / ACP5 / MMP9) | present once each, protein_coding, TMM filter TRUE |
| overlap with GSE155698, GSE156405, GSE212966, GSE217845, GSE282302, GSE57495, GSE85916 | 0 GSM, 0 specimen core, 0 title |
| correlations computed in this step | **0** |
| per-sample SIGLEC15 / ACP5 / MMP9 values loaded in this step | **0** |
