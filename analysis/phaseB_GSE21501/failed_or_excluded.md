# failed_or_excluded.md — GSE21501

Phase B step 3b. Two separate ledgers, kept apart on purpose:

* **Part A — EXCLUDED SAMPLES.** Every array that is not in the main meta-analysis, with its
  reason. There are 30, and none is silently skipped.
* **Part B — FAILED / ABANDONED ATTEMPTS.** Every command, endpoint and approach that did
  not work, what was substituted, and what it cost.

No ρ was computed and no per-sample expression value of SIGLEC15, ACP5 or MMP9 was loaded
at any point in either ledger.

---

# PART A — EXCLUDED SAMPLES

**30 of 132 arrays are excluded from the main meta-analysis. 0 arrays are excluded from
everything** — all 30 remain runnable as prespecified sensitivity **S1 (n = 129)**.

## A.1 Exclusion reason codes

| code | meaning |
|---|---|
| **EX-1** | not clinically annotated: the series overall design states the 30 Nebraska+UNC arrays "were not analyzed with clinical data so it is not provided", and no clinical characteristic is present for the array |
| **EX-2** | tumour-cell fraction not established: the channel-2 free-text `tumor type` declares a **mixed tissue composition with explicit percentages** |
| **EX-3** | tissue compartment contradicted: the channel-2 free-text reads `KV:Tu in LN` (tumour in lymph node), contradicting the deposited `tissue: Patient Primary Pancreatic Tumor` label |
| **EX-4** | per-sample eligibility not established from the deposit: only the single channel-2 `tissue` string; the publication (PMC2903589) corroborates these arrays only at **group** level and publishes no GSM-level mapping |
| **EX-5** | duplicate array: another array in the deposit represents the same specimen (frozen rule DUP-2) |

## A.2 The 30, one row each

### Nebraska — Rapid Autopsy Pancreatic Program (15 arrays, 12 distinct specimens)

| # | GSM | title | reasons | note |
|---|---|---|---|---|
| 1 | GSM535948 | `12R2585B` | EX-1, EX-2, EX-4 | `Mix Pancreas (20/80/0)` |
| 2 | GSM535949 | `12R2587.1` | EX-1, EX-2, EX-4 | `Mix(50%N) P`; **retained** member of `NE:12R2587` |
| 3 | GSM535950 | `12R2587.2` | EX-1, EX-2, EX-4, **EX-5** | `Mix(80%N) P`; dropped by DUP-2b |
| 4 | GSM535951 | `12R2685T` | EX-1, EX-4 | `Panc`; no composition declared |
| 5 | GSM536032 | `13R2778` | EX-1, EX-2, EX-4 | `Mix(90%N) Panc`; **retained** member of `NE:13R2778` |
| 6 | GSM536033 | `13R2778Rpt` | EX-1, EX-4, **EX-5** | submitter-marked `Rpt`; dropped by DUP-2a |
| 7 | GSM536034 | `2R14B` | EX-1, EX-2, **EX-3**, EX-4 | `Mix Pancreas (KV:Tu in LN) 30/60/10` — nodal |
| 8 | GSM536035 | `3R261B` | EX-1, EX-4 | `T Panc` |
| 9 | GSM536036 | `4R678B2` | EX-1, EX-2, EX-4 | `Mix Panc` |
| 10 | GSM536037 | `5R913` | EX-1, EX-2, EX-4 | `Mix Pancreas 50/30/20` |
| 11 | GSM536038 | `7R1353` | EX-1, EX-2, EX-4 | `Brpt Panc 70/30/0` |
| 12 | GSM536039 | `8R1405` | EX-1, EX-2, EX-4 | `Mix 20/60/20` |
| 13 | GSM536040 | `9R1701.2` | EX-1, EX-2, EX-4 | `P Mix`; **retained** member of `NE:9R1701` |
| 14 | GSM536041 | `9R1701B` | EX-1, EX-2, EX-4, **EX-5** | `Mix Panc (10/90/0)`; dropped by DUP-2b |
| 15 | GSM536042 | `9R1703` | EX-1, EX-2, EX-4 | `Mix Pancreas` |

### UNC — Yeh Lab (15 arrays, 15 distinct specimens)

| # | GSM | title | reasons |
|---|---|---|---|
| 16 | GSM536043 | `P01-0505B` | EX-1, EX-4 |
| 17 | GSM536044 | `P03-0175AB` | EX-1, EX-4 |
| 18 | GSM536045 | `P04-0234B` | EX-1, EX-4 |
| 19 | GSM536046 | `P04-0316B` | EX-1, EX-4 |
| 20 | GSM536443 | `P04-0379B` | EX-1, EX-4 |
| 21 | GSM536444 | `P05-0175B` | EX-1, EX-4 |
| 22 | GSM536445 | `P05-0268B` | EX-1, EX-4 |
| 23 | GSM536446 | `P06-0925B` | EX-1, EX-4 |
| 24 | GSM536447 | `P06-1118B` | EX-1, EX-4 |
| 25 | GSM536448 | `P06-1366B` | EX-1, EX-4 |
| 26 | GSM536449 | `P06-1552B` | EX-1, EX-4 |
| 27 | GSM536450 | `P070710B` | EX-1, EX-4 |
| 28 | GSM536451 | `P071118B` | EX-1, EX-4 |
| 29 | GSM536452 | `P071155B` | EX-1, EX-4 |
| 30 | GSM536453 | `P99-0320B` | EX-1, EX-4 |

**Reason-group tally (exactly reproducible from `sample_manifest.csv`):**

| reason combination | n |
|---|---|
| EX-1 + EX-2 + EX-3 + EX-4 | 1 |
| EX-1 + EX-2 + EX-4 | 7 |
| EX-1 + EX-2 + EX-4 + EX-5 | 1 |
| EX-1 + EX-4 | 19 |
| EX-1 + EX-4 + EX-5 | 2 |
| **total** | **30** |

## A.3 The counter-argument, recorded

**14 of the 30 are not excluded on any positive evidence of ineligibility.** The 15 UNC1
arrays are, per the publication's own Table 1, resected primary PDAC with **M0 = 15/15**,
pathologist-reviewed ("All available samples were reviewed by a single pathologist (KAV)"),
with T/N stage, grade, margin and survival reported in aggregate. They are excluded because
that record **is not in the deposit** and **cannot be attached to an individual GSM** — not
because anything says they are not tumours.

An author who reads the frozen estimand as expression-only may prefer **n = 129** as the
main analysis. Prespecified sensitivity **S1** exists so that this is measured rather than
argued. See `frozen_cohort_definition.md` §5.3 and AD-1.

## A.4 Arrays excluded from everything

**None.** The three arrays removed by the frozen duplicate rule (`GSM535950`,
`GSM536033`, `GSM536041`) are inside the excluded Nebraska group, so they are absent from
both the main analysis and S1 — but they are recorded in
`_probe/GSE21501_duplicate_resolution.csv` with the rule that removed them and the title
that was retained in their place. Nothing disappears without a row.

---

# PART B — FAILED / ABANDONED ATTEMPTS

## B.1 `pwsh` is not available on this host

* **Attempt:** `pwsh -NoProfile -File scripts/run_logged.ps1 ...` — the pattern inherited
  from the sibling GSE205154 run.
* **Result:** `The term 'pwsh' is not recognized as the name of a cmdlet, function, script
  file, or operable program.` **exit 1.**
* **Diagnosis:** the shell on this host is **Windows PowerShell 5.1**
  (`$PSVersionTable.PSVersion` = 5.1.26100.9444, `$PSHOME` =
  `C:\Windows\System32\WindowsPowerShell\v1.0`), and there is no PowerShell 7 anywhere on
  `PATH`.
* **Substitution:** invoke the logger as
  `powershell.exe -NoProfile -ExecutionPolicy Bypass -File scripts\run_logged.ps1 ...`.
* **Second-order fix:** PS 5.1's `Add-Content -Encoding UTF8` writes a **BOM**, which the
  sibling's log does not have. `run_logged.ps1` was changed to append through
  `[System.IO.File]::AppendAllText(...)` with `UTF8Encoding($false)` so the log stays
  BOM-free and byte-comparable in style with the sibling run.

## B.2 R was not on `PATH`

* **Attempt:** `Rscript --version`.
* **Result:** `The term 'Rscript' is not recognized...` **exit 1.**
* **Diagnosis:** R is installed at `C:\Program Files\R\R-4.5.1` but is not on `PATH`.
* **Resolution:** **no substitution was needed and R was not used.** Every operation in this
  step is a metadata operation — parsing a series-matrix header, parsing a platform
  annotation, deriving title keys and intersecting GSM sets. None of it needs R, and no R
  package contributes to any number in this directory. This is recorded rather than
  silently dropped because `sessionInfo.txt` is asked for "exact R/Python package versions"
  and the honest answer for R is "not used", not a version string.

## B.3 `GSE21501_family.soft.gz` — 187 MB, stalled, abandoned

* **Attempt:** download `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE21nnn/GSE21501/soft/GSE21501_family.soft.gz`
  in the same script as the series matrix, via a single `urllib` connection.
* **Result:** the series matrix (12 575 998 bytes) completed in ~4 minutes; the family SOFT
  transfer then ran **past the 600 s tool timeout without completing**, and was killed.
  Two orphaned `python` processes were left holding the connection and had to be
  `Stop-Process -Force`d.
* **Diagnosis:** a `HEAD` request put the file at **187 687 633 bytes**. For a two-channel
  family SOFT the bulk of that is 132 `!sample_table` blocks of 45 220 rows each — i.e. the
  file is ~99 % expression data, which this step must not read anyway. This is the same
  single-connection stall class the sibling run hit at exactly 8 388 608 bytes.
* **Root cause of the failure, precisely:** this was **not** purely a size problem. It was a
  **wrong-source problem** — the family SOFT was being fetched for its *metadata*, which is a
  tiny fraction of it, and every byte of the rest was data the guard would have had to
  discard.
* **Substitution, two independent sources, both cheap and both metadata-only:**
  1. `GSE21501_series_matrix.txt.gz` (12 575 998 bytes) — its **header** holds all 132
     samples' channel-1 and channel-2 fields. Download completed; sha256 recorded in
     `_probe/GSE21501_download_manifest.txt`.
  2. `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21501&targ=gsm&form=text&view=brief`
     — 351 980 bytes, all 132 `^SAMPLE` brief records, no sample tables.
* **Outcome:** the two sources **agree field-for-field** on every eligibility-relevant
  value (all 132 channel-2 `tissue`, all five channel-2 providers, all `tumor type`
  strings, all clinical characteristics). Nothing in this directory depends on the family
  SOFT, and `sessionInfo.txt` says so explicitly.
* **Not retried** with a parallel range downloader. The sibling's 8-way downloader exists
  and works, but the file it would fetch is 99 % data this step is forbidden to read, so
  retrying it would have added a large amount of prohibited content to the workspace for no
  metadata gain. This is a deliberate decision, not an omission.

## B.4 GEO's SOFT endpoint via `web_fetch` — rejected

* **Attempt (inherited note, re-verified):** fetch GEO's `geo/text` SOFT endpoint with
  `web_fetch`.
* **Result:** rejected by the fetch layer. GEO HTML accession pages are reCAPTCHA-protected.
* **Substitution:** PowerShell `Invoke-WebRequest`, which works against both the FTP SOFT
  paths and the `acc.cgi` query endpoint (used for all of B.3's source 2, and for the eleven
  comparator series).
* **Side effect worth recording:** `Invoke-WebRequest` returned the `acc.cgi` brief payload
  as a **byte array** rather than text on the first try, producing a file of decimal byte
  values (`94 83 65 77 80 76 69 ...`). Fixed by using `-OutFile`, which writes the response
  body verbatim. The first, corrupt file was overwritten; the delivered
  `_probe/GSE21501_gsm_brief.txt` is 351 980 bytes with sha256
  `7b911201628a1626b7c504b598c4ab84242379c15b6b73c9722db9c0e708b2d1` and 132 `^SAMPLE`
  blocks.

## B.5 Europe PMC full-text XML — HTTP 406

* **Attempt:** `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2903589/fullTextXML`
  to read the publication's clinical tables.
* **Result:** **HTTP 406 Not Acceptable** (`Whitelabel Error Page`; the endpoint has no
  mapping for the request as sent).
* **Also tried:** `https://www.ebi.ac.uk/europepmc/webservices/rest/MED/20644708/textMinedTerms/GENE_PROTEIN`
  → **HTTP 404**.
* **Substitution:** `https://pmc.ncbi.nlm.nih.gov/articles/PMC2903589/` (**HTTP 200**) and
  `https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000307`
  (**HTTP 200**). The Methods, Table 1 and Table 2 were retrieved and are quoted verbatim in
  `_probe/GSE21501_publication_evidence.txt`.
* **Residual loss:** both retrievals truncated in the middle of Table 2, and the Discussion
  and Supplementary sections were not read. Mitigation: Table 2's remaining columns carry
  treatment and margin breakdowns that no frozen eligibility condition depends on, and the
  one quantity this step needed from the publication — the cohort sizes — appears in the
  Methods text and is fully captured. The truncation is recorded rather than papered over.

## B.6 `GPL4133_family.soft.gz` — 6.4 GB, never attempted

* A `HEAD` request put the platform family SOFT at **6 407 144 478 bytes**. It was not
  downloaded. The 9 258 564-byte `GPL4133.annot.gz` supplies the same probe→symbol mapping
  used by the frozen analysis, and was downloaded instead.
* Note for anyone repeating this: `GPL4133.annot.gz` is the **Feature Number version** of
  the platform, so its `ID` column is a bare integer, not an `A_23_P…` probe name. That
  looks like a row index and is not one — it matches the series matrix `ID_REF` values
  exactly, 45 220 distinct values in the same order.

## B.7 Two scripts produced wrong output on the first run and were corrected

Recorded because a manifest is only worth what its error history is worth.

1. **`06_patient_identity.py`, section 3 narration.** The first version printed
   "for every repeat-marked array the OTHER array is absent from the deposit". That is false
   for `13R2778Rpt`, whose base array `13R2778` **is** present. Corrected to state that
   6 of 7 repeat-marked arrays have no sibling, that the 7th does, and that DUP-2 therefore
   drops the repeat for that one and retains it for the other six. No number changed; the
   claim was wrong and is now right.

2. **`07_build_manifest.py`, EX-5 scope.** The first version attached reason **EX-5**
   ("duplicate array") to **both** members of every duplicate group, including the member the
   frozen rule **retained**. EX-5 must apply only to the array that was actually dropped.
   Corrected; EX-5 now appears on exactly 3 rows (`GSM535950`, `GSM536033`, `GSM536041`)
   instead of 6. `eligible_for_main_meta` was unaffected (all six are in the excluded 30),
   but the reason text was wrong. **This correction changed the manifest SHA-256** from
   `1a6593c6285691a93f3d7bfcef18edeaea2db5b8a72fb86e84a3aac304c21105` to the frozen
   `3e11938ae5c7e682c9e53aad93d051579764a4ec8173d69d571dbeda0f3ac4e4`; only the latter is
   recorded in `frozen_cohort_definition.md`.

3. **`07_build_manifest.py`, EX-2 and EX-4 wording.** The first version asserted that the
   Nebraska free text declared a tumour/normal "admixture with an explicit mixing ratio".
   The publication does not define the three percentages, so the stronger word "admixture"
   was not supportable. Reworded to "declares a MIXED tissue composition with explicit
   percentages, so the tumour-cell fraction of this specimen is unknown", and EX-4 was
   reworded to cite the publication's group-level corroboration explicitly. This also
   changed the SHA-256 (to the frozen value).

## B.8 A stray per-sample table was downloaded and immediately destroyed

* **Attempt:** to test for a hidden patient/relation field, `acc.cgi?acc=GSM536894&…&view=full`
  and `acc.cgi?acc=GSM536488&…&view=full` were fetched.
* **Result:** **HTTP 200**, 3 567 297 and 3 563 804 bytes. `view=full` embeds the complete
  45 220-row sample table, i.e. the per-sample expression values of SIGLEC15, ACP5 and MMP9
  for those two arrays.
* **What was and was not read:** the files were read **only** through a
  `Select-String` filter for the patterns `relation|BioSample|SRA|title|characteristics|
  description|data_row|supplementary`. The `#CORR = correlation coefficient among pixels`
  line matched that filter by accident (the substring `relation` inside `correlation`), which
  is how it is known that the match was accidental rather than a `!Sample_relation` field:
  **no `!Sample_relation` tag exists anywhere in the 132 brief records** (verified separately,
  count 0). **No expression value was displayed, extracted or bound.**
* **Action taken:** both files were **deleted immediately** after the check
  (`Remove-Item`), so that no per-sample expression table is retained anywhere in this
  directory. Verified deleted.
* **Full-disclosure note:** the bytes passed through the process, as they unavoidably must
  for any HTTP read. The guarantee this run makes — in `frozen_cohort_definition.md` §0 and
  in the `05_probe_mapping.py` docstring — is about **retention and interpretation**, not
  about I/O, and this is the one place where a broader payload than intended was received
  before being discarded. It is recorded here rather than omitted.

## B.9 A subagent delegation was refused

* **Attempt:** delegate the publication reading (PMC2903589, PMC11779630) to a background
  subagent.
* **Result:** `Error: subagent depth 2 exceeds maxDepth 1`. This session is itself a
  delegated subagent, so it cannot spawn children.
* **Substitution:** the publications were read directly with `web_fetch` (see B.5), which
  worked. No information was lost, but it means the publication evidence carries no
  independent second reading.

---

## PART C — WHAT WAS NEVER ATTEMPTED, AND WHY

| not attempted | why |
|---|---|
| `GSE21501_RAW.tar` | per-sample raw Agilent feature files; the series matrix already supplies the deposited quantity and the raw files are forbidden content for this step |
| dbGaP `phs003597.v1.p1` | controlled access, and it belongs to the **GSE205154** sibling cohort, not to GSE21501. It could not settle GSE21501's patient identity in any case |
| any re-quantification from raw intensities | the deposit supplies no raw intensities; the series matrix is the only representation (AD-2) |
| any expression summary, detection rate, variance filter or QC on gene values | forbidden by the governing rule; gene usability was established structurally (presence + non-empty token only) |
| any computation involving SIGLEC15, ACP5 or MMP9 values | forbidden by the governing rule |
