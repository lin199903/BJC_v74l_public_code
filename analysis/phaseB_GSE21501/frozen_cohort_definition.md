# frozen_cohort_definition.md — GSE21501

**Status:** FROZEN. The later ρ step must obey this document unconditionally.
**Step:** Phase B step 3b — manifest construction and cohort-definition freeze only.
**Date frozen:** 2026-09-28 (Asia/Shanghai)
**Frozen artefact:** `sample_manifest.csv` (132 rows × 14 columns)
**SHA-256 of `sample_manifest.csv`:**

```
3e11938ae5c7e682c9e53aad93d051579764a4ec8173d69d571dbeda0f3ac4e4
```

> A later run proves it ran after this freeze by re-hashing `sample_manifest.csv` and
> checking it against the value above. If the hash differs, the manifest was edited after
> the freeze and the estimate that used it is not covered by this definition.

---

## 0. Governing rule, and what was NOT done in this step

> **Lock the cohort BEFORE looking at the effect.**

**No correlation between SIGLEC15 and any ACP5/MMP9 programme score was computed,
inspected or reported in this step. No per-sample expression value of SIGLEC15, ACP5 or
MMP9 was loaded.**

The guard that enforces this is **value-blind guard v2**, implemented in
`scripts/05_probe_mapping.py`:

1. only the data-block **header** row is read, to obtain column labels;
2. for every data row, **only the bytes before the first tab** are sliced out — that is the
   probe identifier, which is platform feature metadata identical for every sample;
3. for the rows whose identifier is one of the target probes, a cell is recorded as empty
   **iff two consecutive tab positions are adjacent** (`cell_end <= cell_start`) and as null
   iff its 4-byte slice equals `b"null"`. The cell's bytes are never bound to a variable,
   never decoded, never converted to a number. Only booleans and counts leave the loop;
4. the emptiness detector is **self-tested** against seven synthetic lines with known empty
   and null cells before it is pointed at real data, so that "0 empty cells" is evidence
   about the data and not about a broken detector. The self-test passes
   (`_probe/GSE21501_probe_mapping.txt` §0).

There is no code path in this run that turns a cell into a number, and no output in this
directory contains one.

**Frozen estimand (verbatim, for context only — not computed here):**

> within each eligible PDAC cohort, estimate the Spearman correlation between
> tumour-level SIGLEC15 expression and an ACP5/MMP9 programme score constructed using the
> same dataset-level standardisation rule as the existing TCGA-PAAD, GSE57495 and GSE85916
> analyses. Cohort-specific estimates will be combined only after cohort eligibility is
> established independently of the observed correlation.

**Exclusion rule carried over unchanged:** a cohort is **never** excluded because of the
observed ρ magnitude or direction.

---

## 1. Cohort identity and deposit

| item | value |
|---|---|
| accession | GSE21501 |
| title | *A six-gene signature predicts survival of patients with localized pancreatic ductal adenocarcinoma* |
| BioProject | PRJNA125999 |
| platform | **GPL4133** — Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (**Feature Number version**) |
| design | `!Series_overall_design = "reference x sample"`; `!Sample_channel_count = 2` for all 132 |
| arrays | 132 GSM (GSM535948–GSM536949 range), one series-matrix column each |
| features per array | 45 220, identical to the platform row count |
| publication | PMID 20644708 — Stratford JK *et al.*, *A six-gene signature predicts survival of patients with localized pancreatic ductal adenocarcinoma*, **PLoS Medicine** 2010;7(7):e1000307, doi 10.1371/journal.pmed.1000307 (PMC2903589) |
| second PubMed id | 28380042 |
| public since | 2010-07-19 |
| suppliers | Northwestern Memorial Hospital 49; Johns Hopkins Medical Institutions 34; NorthShore University HealthSystem 19; University of Nebraska Medical Center Rapid Autopsy Pancreatic Program 15; Yeh Lab-UNC-CH 15 |

Evidence for every row above is in `_probe/GSE21501_series_level.txt`,
`_probe/GSE21501_metadata_field_audit.txt` and
`_probe/GSE21501_publication_evidence.txt`.

### 1.1 What representation the deposit supplies

The deposit supplies **one** matrix: `GSE21501_series_matrix.txt.gz`. Its `VALUE` column is
the Agilent Feature Extraction output described identically for all 132 samples:

> `!Sample_data_processing = "Data for both channels were Lowess-normalized and then log(2)
> ratio was taken using Agilent Feature Extraction (v10.5) software."`

and the publication states the same thing from the other side:

> "cDNA was labeled with Cy5-dUTP and a reference control (Stratagene) was labeled with
> Cy3-dUTP … The log2 ratio of the mean red intensity over mean green intensity was
> calculated for each gene and went through LOWESS normalization."

So the deposited quantity is **log2(tumour / reference)**, Lowess-normalised. It is already
on a log2 scale. **The CP10K + log1p single-cell rule was not applied** — it governs raw
single-cell counts, and this is a two-colour bulk array deposit. Nothing was re-quantified.

### 1.2 ⚠️ GSE21501 IS A TWO-CHANNEL ARRAY — and this is the trap

| | channel 1 | channel 2 |
|---|---|---|
| label | **Cy3** | **Cy5** |
| source name | `Human Reference` (all 132) | `Pancreas Tissue` (all 132) |
| `characteristics` | `tissue: Stratagene Human reference RNA` (**131**) and `reference: Stratagene …` (**1**) | `tissue: Patient Primary Pancreatic Tumor` (**132**) |
| provider | `Stratagene (Part # 74000_41, lot 1130623)` (all 132) | five institutions (see §1) |

**All sample-level metadata used for eligibility is taken from a channel-2 (`_ch2`) field.
Nothing is taken from channel 1.** The channel-1 fields were parsed, but only to quantify
how misleading they are:

* **Channel 1 is the same reference in every array.** Its `tissue` value is
  `Stratagene Human reference RNA` for 131 arrays and `reference: Stratagene …` for one,
  and its provider is the same Stratagene lot with the same part number for all 132.
* **A reader who took the series-level or channel-1 characteristics would conclude that
  all 132 samples are reference RNA**, and would either empty the cohort or, worse, treat a
  constant as a specimen label. The `tissue` field that actually discriminates specimens
  exists only at `_ch2`.
* The channel-1 key is not even spelled consistently: one array carries `reference:` where
  131 carry `tissue:`. A pipeline keyed on `tissue:` alone would silently drop that array
  from channel 1 — and would never notice, because the array's real specimen data is on
  channel 2.
* The array title (`!Sample_title`) is **per array, not per channel**: GEO emits one
  `!Sample_title` per GSM and there is no `!Sample_title_ch1`/`_ch2` tag anywhere in the
  deposit. The manifest column `channel2_sample_title` therefore carries that single title,
  which names the channel-2 specimen, i.e. the tumour. This is stated so that the column
  name cannot be mistaken for a claim that GEO publishes two titles per array.

### 1.3 What `!Sample_title` does and does not encode

Titles carry an institution-specific specimen code (`PH10162T`, `PE133T`, `P109T`,
`12R2585B`, `P01-0505B`). They do **not** carry a patient identifier, and six of them carry
the submitter's own word that the array is a repeat (`P164T redo`, `P180T Replicate`,
`P271N Replicate`, `P275T Replicate`, `P31T Replicate`, `P41T Replicate`) plus one more
(`13R2778Rpt`). See §4.

---

## 2. How it was established that a sample is a primary PDAC tumour

**Primary metadata field:**

```
!Sample_characteristics_ch2 = "tissue: <value>"
```

**This deposit has exactly one value:**

| value | maps to | n |
|---|---|---|
| `Patient Primary Pancreatic Tumor` | primary pancreatic tumour | **132** |

**No channel-2 field declares a metastasis.** The only other channel-2 characteristic that
could describe a specimen is the free-text `tumor type`, present for 32 arrays, whose value
is a specimen code (`"tumor type: 12R2585B Mix Pancreas (20/80/0)"`), not a
primary/metastatic classification. It is **not** the analogue of GSE205154's `tumor type`
field and must not be read as one.

**Tumour vs non-tumour.** All 132 arrays are tumour by the deposited label, and the series
contains no adjacent-normal and no non-tumour array. **The "tumour tissue" condition
therefore excludes nothing on this axis.** The publication confirms the specimens are
tumour: "tissues were macrodissected in order to preserve the normal adjacent tissue and
stroma of the tumors", and "All available samples were reviewed by a single pathologist
(KAV)."

**Disease.** The publication identifies the cohort as PDAC throughout and reports TNM,
grade and margin for the majority of patients. T stage is populated for 98 of the 132
arrays, N stage for 101, risk group for 102, overall-survival time and event for 102.

**Two specimens within the 132 are flagged here rather than silently absorbed:**

* **`GSM536488` / `P271N Replicate`** — the only Northwestern title ending in `N`; all 48
  others end in `T`. The `N` is undocumented in both the deposit and the publication. The
  deposited channel-2 `tissue` calls it `Patient Primary Pancreatic Tumor`, it carries the
  full clinical record (`os time 11`, `os event 1`, `t stage 3`, `n stage 1`, `Low Risk`)
  and its `description` is character-for-character the same form as the tumour arrays, so
  the frozen rule (key on the deposited `tissue` field, not on the title suffix) **retains**
  it. If `N` denotes adjacent normal tissue this is wrong. Flagged as **AD-4**; removed in
  prespecified sensitivity S2.
* **`GSM536901` / `P188T IPMN`** — IPMN (intraductal papillary mucinous neoplasm) is a
  precursor lesion, not invasive PDAC. The deposited `tissue` field nevertheless calls it
  `Patient Primary Pancreatic Tumor`, and it carries a full clinical record
  (`os time 6`, `os event 1`, `t stage 1`, `n stage 1`, `High Risk`). §6 shows this array
  is also the most likely explanation of a one-array discrepancy against the publication's
  own cohort count. Retained by the frozen rule; flagged as **AD-4**; removed in S2.

Neither the deposit nor the publication mentions IPMN anywhere else.

---

## 3. Rule for patient identity

### 3.1 What GEO exposes

**Nothing.** Verified over all 132 brief GSM records (`_probe/GSE21501_gsm_brief.txt`):

* the distinct `!Sample_*` tags in the deposit are 42, and **none** of them is
  `!Sample_relation`, `!Sample_subject`, `!Sample_donor`, `!Sample_individual` or any
  other patient field;
* there is **no BioSample accession** in the deposit at all — the strings `BioSample`,
  `SRA` and `bioproject` occur zero times across the 132 records — so there is no
  secondary record through which a patient identifier could be requested;
* the only patient-linked handle is the **specimen code inside `!Sample_title`**.

This is a genuine absence, not a lookup that failed. The sibling GSE205154 run was able to
reach a published `public patient_id` through a companion SubSeries; **GSE21501 has no
companion series, no SuperSeries and no BioSample linkage**, so that route does not exist
here.

### 3.2 The frozen specimen-key rule

> **FROZEN — specimen-key rule.**
> The key is the **provider-scoped specimen code derived from `!Sample_title`**, prefixed by
> a two-letter provider namespace so that the five unrelated naming namespaces can never be
> accidentally joined:
>
> | provider | title pattern | key |
> |---|---|---|
> | Johns Hopkins | `PH[\d-]+T` | `JH:PH10162` |
> | NorthShore | `PE\d+T` | `NS:PE133` |
> | Northwestern | `P\d+T…` or `P\d+N…` | `NW:P109` |
> | Nebraska | `\d+R\d+…` | `NE:12R2587` |
> | UNC | `P\d{2}-?\d{4}` | `UN:P01-0505` |
>
> The key is recorded for every sample in the manifest column `patient_id_if_available`,
> which is prefixed `not_exposed_by_GEO;` precisely so that it cannot be mistaken for a
> patient identifier issued by the archive.
>
> **The key is an identity handle, and it is also an upper bound, never a count.** Two
> differently coded specimens could belong to one patient and the deposit exposes nothing
> that would reveal it.

### 3.3 Arrays ≠ patients — proved from the deposit's own text

Seven arrays carry the submitter's word that they are a repeat:

| GSM | title | marker | sibling array in the deposit? |
|---|---|---|---|
| GSM536894 | `P164T redo` | redo | **no** |
| GSM536897 | `P180T Replicate` | replicate | **no** |
| GSM536488 | `P271N Replicate` | replicate | **no** |
| GSM536922 | `P275T Replicate` | replicate | **no** |
| GSM536923 | `P31T Replicate` | replicate | **no** |
| GSM536926 | `P41T Replicate` | replicate | **no** |
| GSM536033 | `13R2778Rpt` | rpt | **yes** — `13R2778` |

So the deposit states in its own words that six of its arrays are a second assay of a
specimen, and one of those pairs is present in full. Three further same-specimen pairs
carry no marker at all and are detected only from the code (`12R2587.1`/`12R2587.2`;
`9R1701.2`/`9R1701B`).

**Independent cross-check against the publication.** PMC2903589 reports its clinical
cohorts as:

> "The training cohort included 34 patients with resected PDAC from Johns Hopkins Medical
> Institutions (JHMI). The testing or validation cohort included patients from two
> institutions: 48 from Northwestern Memorial Hospital (NW) and 19 from NorthShore
> University HealthSystem (NSU)."

34 + 48 + 19 = **101 published patients**. The deposit carries 34 + **49** + 19 = **102**
clinically annotated arrays. Johns Hopkins agrees at 34 and NorthShore agrees at 19; the
Northwestern arm is one array larger than the published patient count. **Therefore at least
one of the 49 Northwestern arrays does not correspond to a distinct patient in the
published cohort.** The two candidates consistent with everything else in the record are
`GSM536901 P188T IPMN` and `GSM536488 P271N Replicate`; the deposit cannot decide between
them, and neither can the publication, which publishes no GSM-level table.

### 3.4 Interpretation limits that follow, and must be carried into any write-up

1. The effective unit is the **specimen**, and the specimen is at best a proxy for the
   patient. It is a **within-deposit** key, not an externally published identifier. A
   reader **cannot** re-derive the patient grouping from public metadata alone — the
   publications would have to be obtained, and even then they publish only aggregate
   tables.
2. **No exact patient count is claimed anywhere in this record.** The main analysis rests
   on 102 arrays and the number of distinct patients behind them is **≤ 102**; the whole
   deposit rests on 132 arrays and **≤ 129** specimens (and, by the publication's own
   cohort totals, ≤ 101 of the 102 clinically annotated arrays are published patients —
   see AD-4).
3. **Patient-level overlap with TCGA-PAAD is not assessable** — no shared accession
   namespace, different institutions.
4. **Cross-provider identity is not assessable.** Five unrelated naming namespaces; no
   field licenses joining them. The provider prefix in the key makes this explicit rather
   than silently merging same-looking codes.
5. **The meaning of the trailing `N` in `P271N Replicate` is unknown** (§2, AD-4).

---

## 4. Rule for duplicate / replicate / redo arrays

> **FROZEN — duplicate rule, written down before any effect was examined.**
>
> **DUP-1** The specimen key is the provider-scoped title-derived code of §3.2.
>
> **DUP-2** If one specimen key is represented by **more than one array**, exactly **one**
> enters the analysis:
> **(a)** prefer the array whose title carries **no** repeat marker (`rpt`, `redo`,
> `Replicate`, case-insensitive);
> **(b)** if that does not single out one array, retain the **lexicographically smallest
> `sample_title`** in the group.
>
> **DUP-3** A repeat-marked array whose specimen key has **no other array** in the deposit
> is **retained**: it is the sole representation of that specimen, and dropping it would
> delete a specimen rather than de-duplicate one.
>
> **DUP-4** Expression values are **never averaged** across the group. Selection was chosen
> over aggregation because aggregation would introduce a specimen-combination rule that no
> existing cohort in this manuscript uses, and the frozen text requires only that *a* rule
> be fixed in advance.
>
> **DUP-5** The rule is outcome-blind by construction: DUP-2 reads only the sample title
> and the specimen key, so no expression value can influence it.
>
> **DUP-6 Pre-specified sensitivity (must be reported whatever it shows):** the analysis is
> repeated with (i) the duplicate pair **averaged at the expression level** where a pair
> exists, and (ii) **all** repeat-marked arrays **dropped**.

**What DUP-2 actually removes.** Exactly three arrays, all inside the Nebraska group:

| specimen key | arrays | retained | dropped | rule |
|---|---|---|---|---|
| `NE:12R2587` | `12R2587.1`, `12R2587.2` | `12R2587.1` (GSM535949) | `12R2587.2` (GSM535950) | DUP-2b (both unmarked) |
| `NE:13R2778` | `13R2778`, `13R2778Rpt` | `13R2778` (GSM536032) | `13R2778Rpt` (GSM536033) | DUP-2a |
| `NE:9R1701` | `9R1701.2`, `9R1701B` | `9R1701.2` (GSM536040) | `9R1701B` (GSM536041) | DUP-2b (both unmarked) |

**DUP-2 removes nothing from the main analysis.** Each of the 102 clinically annotated
arrays has a specimen key of its own: 34 Johns Hopkins + 19 NorthShore + 49 Northwestern =
102 keys for 102 arrays. The six `redo`/`Replicate` arrays of the 102 are retained under
DUP-3 because their originals are not in the deposit.

---

## 5. Which arrays carry clinical annotation, and the treatment of the 30

### 5.1 How the annotation split was derived (not hard-coded)

Three independent derivations, all from the deposit, all agreeing:

| derivation | result |
|---|---|
| **(a) provider split of channel 2** | Nebraska 15 + UNC 15 = **30** |
| **(b) arrays carrying ≥ 1 of `t stage`, `n stage`, `risk group`, `os time`, `os event`** | **102** |
| **(c) arrays carrying none of those** | **30** |
| **(d) do (a) and (c) coincide?** | **yes, exactly** — every array with no clinical field is Nebraska or UNC, and every Nebraska/UNC array has no clinical field |

The deposit says so itself. `!Series_overall_design` carries two lines:

> `"reference x sample"`
> `"The 30 Nebraska and UNC samples were not analyzed with clinical data so it is not
> provided. One of the PE samples is missing some clinical data."`

**The second sentence of that design line is resolvable, and was resolved.**
"One of the PE samples is missing some clinical data" is **`GSM536940` / `PE621T`**
(NorthShore): it carries `os time 17`, `os event 0`, `t stage 2` and `Risk Group: High Risk`
but **no `n stage`**. It is the only PE array with an incomplete record. Four further
arrays lack a T stage (`GSM536892 P152T`, `GSM536924 P35T`, `GSM536927 P42T`,
`GSM536946 PE734T`) and carry the description `"Pancreatic Tumor. , N Stage/…"`.

**All five incomplete-record arrays stay in the main analysis.** The frozen eligibility
conditions do not require a complete clinical record; they require metadata sufficient to
place the specimen in a compartment, which `tissue` supplies for every one of the 132.

### 5.2 What the publication adds about the 30

| | UNC1, n = 15 | NEB, n = 15 |
|---|---|---|
| what they are (publication, verbatim) | "15 patients with resected primary PDAC from the University of North Carolina at Chapel Hill (UNC)" | "15 patients with metastatic PDAC from the University of Nebraska Medical Center Rapid Autopsy Pancreatic Program (NEB)" |
| specimen handling | flash frozen at operation, 1999–2007 | "organs were harvested within 3 h post mortem and the specimens flash frozen in liquid nitrogen" |
| M stage (Table 1) | **M0 = 15 (100 %)** | **M1 = 15 (100 %)** |
| T / N / grade / margin (Table 1) | given, as aggregates | `NA` throughout |
| survival | median 9 (1–35) mo | `NA` |
| published at GSM level? | **no** | **no** |

The publication therefore corroborates the 30 **only at group level**, and publishes no
mapping from any GSM to a row of Table 1. In the deposit itself the 30 carry one
specimen field and nothing else.

Two per-sample facts inside the Nebraska branch **contradict** the deposited label:

* **12 of the 15 Nebraska arrays** declare a **mixed tissue composition with explicit
  percentages** in the channel-2 free text — `12R2585B Mix Pancreas (20/80/0)`,
  `12R2587.1 Mix(50%N) P`, `12R2587.2 Mix(80%N) P`, `13R2778 Mix(90%N) Panc`,
  `2R14B Mix Pancreas (KV:Tu in LN) 30/60/10`, `4R678B2 Mix Panc`,
  `5R913 Mix Pancreas 50/30/20`, `7R1353Brpt Panc 70/30/0`, `8R1405 Mix 20/60/20`,
  `9R1701.2 P Mix`, `9R1701B Mix Panc (10/90/0)`, `9R1703 Mix Pancreas`. **The tumour-cell
  fraction of those specimens is therefore not established.** The publication does not
  define the three numbers, so they cannot be converted into a purity estimate — but the
  specimen is on its face not a pure tumour, and the remaining 3 Nebraska arrays
  (`12R2685T Panc`, `13R2778Rpt`, `3R261B T Panc`) carry no such declaration only because
  they carry almost no annotation at all.
* **`2R14B`** reads `Mix Pancreas (KV:Tu in LN)` — **tumour in lymph node**. That directly
  contradicts `tissue: Patient Primary Pancreatic Tumor` for the same array. With the
  publication reporting M1 = 15/15 for the Nebraska branch, a nodal specimen is entirely
  plausible. For this array primary and metastatic **cannot be told apart from the
  deposit**, which is exactly what the frozen "Metadata" condition forbids.

### 5.3 The decision, and its reasoning

**The frozen text does not settle whether clinical annotation is required.** The estimand
is a within-cohort correlation between two expression quantities; no clinical covariate
enters it, and the eligibility table's Disease/Specimen/Unit/SIGLEC15/ACP5/MMP9/Platform
conditions say nothing about annotation. Requiring it is therefore a **choice**, and it is
recorded as **AD-1**, not presented as a pre-registered rule.

**Chosen treatment:**

| | population | n |
|---|---|---|
| **Main meta-analysis** | clinically annotated arrays (Johns Hopkins 34, NorthShore 19, Northwestern 49) | **102** |
| **Excluded from the main analysis, retained as prespecified sensitivity S1** | Nebraska 15 + UNC 15 | **30** |
| **Prespecified sensitivity S1** | whole deposit after DUP-2 | **129** |

**Reasoning, in the order that it weighs:**

1. The 30 are the deposit's own partition: **it** says they were not analysed with clinical
   data. Taking the deposit at its word is the conservative reading and it is checkable
   (§5.1(d)).
2. For the 30, **eligibility is only ever a group statement**. No field, no table and no
   accession ties any individual one of the 30 GSMs to a pathologist-reviewed, staged case.
   The frozen table requires "metadata sufficient to distinguish primary / metastatic /
   non-tumour" — a group-level assertion plus a single label does not satisfy that for an
   *individual array*, which is the unit the manifest grades.
3. For the Nebraska branch there is worse than absence of evidence: **12 of 15 declare a
   mixed tissue composition** and **one declares nodal tumour**, contradicting the label
   that the whole cohort's "primary" classification rests on.
4. The Nebraska branch is **post-mortem rapid-autopsy tissue** with a pre-analytical
   history (up to 3 h post mortem) different from every other array in the deposit, which
   is surgical resection flash-frozen at operation. The publication itself separates
   "decedent and resected PDAC samples" throughout.
5. The 30 contain **all three** of the deposit's same-specimen duplicate groups, i.e. they
   are the only part of the deposit where the "one independent tumour/patient" condition
   actually binds.

**And the counter-argument, recorded because it is real:** the 15 UNC1 arrays are resected
primary PDAC, M0 by the publication's own Table 1 (15/15), pathologist-reviewed, and look
eligible on every axis except that their clinical record is absent from GEO. An author who
reads the frozen estimand as expression-only may reasonably prefer **n = 129** as the main
analysis and 102 as the sensitivity. **S1 exists precisely so that this can be seen rather
than argued**, and it must be reported whatever it shows.

---

## 6. Frozen eligibility checks, condition by condition

| condition | frozen rule | result for GSE21501 |
|---|---|---|
| **Disease** | histologically identified PDAC | **PASS for the 102** — PDAC throughout the deposit and the publication; "All available samples were reviewed by a single pathologist (KAV)"; TNM/grade/margin reported. **NOT ESTABLISHED per array for the 30** — group-level only (§5.2, EX-4) |
| **Specimen** | tumour tissue; primary tumour preferred for the main meta-analysis | **PASS for the 102** — all resected primary pancreatic tumours, flash frozen at operation 1999–2007. **FAILS for 12 of the 30** — declared mixed tissue compositions, one of them nodal (EX-2, EX-3) |
| **Unit** | one independent tumour/patient; duplicate specimens require a pre-defined aggregation/selection rule | **PASS after DUP-2** — rule frozen in §4 before any effect was examined; it removes 3 arrays, all inside the 30, and 0 from the 102 |
| **SIGLEC15** | interpretable gene-level expression | **PASS** — exactly one GPL4133 feature maps to SIGLEC15 (feature 31474, `Gene title = "sialic acid binding Ig like lectin 15"`); present in the matrix and complete for all 132 arrays (§7) |
| **ACP5 / MMP9** | both measurable, score constructible by the established rule | **PASS** — ACP5: exactly one feature (10450, `"acid phosphatase 5, tartrate resistant"`); MMP9: ten features (1570, 1678, 14430, 14852, 21363, 25413, 28592, 28675, 40398, 44612), all `"matrix metallopeptidase 9"`, collapsed to one gene value by the existing rule (§7) |
| **Metadata** | sufficient to distinguish primary / metastatic / non-tumour | **PASS for the 102** — channel-2 `tissue` plus a full clinical record. **FAILS for the 30** — one specimen field only; for `2R14B` the free text contradicts it (EX-3) |
| **Platform** | within-cohort rank information obtainable; cross-platform absolute comparability NOT required | **PASS** — a single platform (GPL4133) and a single 45 220-feature matrix covering all 132 arrays with an identical common reference. Rank information is available within the cohort. **Cross-platform absolute comparability is neither required nor claimed**, which matters here because a two-colour log2 ratio is a *relative* quantity (§8) |
| **Analysis** | cohort-specific Spearman → random-effects meta-analysis | **DEFERRED** — not performed in this step |
| **Exclusion** | never exclude a cohort because of the observed ρ magnitude or direction | **ACKNOWLEDGED AND BINDING** — no ρ was computed in this step, so this condition could not have operated |

### 6.1 Overlap with the manuscript's comparator series — verified, not assumed

Checked against the 132 GSMs, the 132 titles and the 129 specimen keys:

| dataset | series n | GSM overlap |
|---|---|---|
| GSE155698 | 41 | 0 |
| GSE156405 | 17 | 0 |
| GSE205154 | 289 | 0 |
| GSE212966 | 12 | 0 |
| GSE217845 | 16 | 0 |
| GSE278688 | 29 | 0 |
| GSE278689 | 8 | 0 |
| GSE281005 | 505 | 0 |
| GSE282302 | 108 | 0 |
| GSE57495 | 63 | 0 |
| GSE85916 | 80 | 0 |
| TCGA-PAAD | n/a | not assessable |

No overlap is identified with any comparator. The TCGA-PAAD row is not assessable: GSE21501
samples carry no TCGA barcode and no shared accession namespace exists, and the cohorts are
from different institutions.

---

## 7. Probe-level gene mapping (this is an array platform)

**Platform annotation:** `GPL4133.annot.gz`, sha256
`0b82bf430eb95b985d0ca1043bb9ff98678d3ead613f74cb65ec5ae0416c3f85`, 45 220 rows, symbol
column `Gene symbol` (index 2) — the **same column-selection rule the frozen analysis uses**
(`parse_soft_platform`, first header column containing "symbol", `split("///")[0]`).

GPL4133 is the **Feature Number version** of the Agilent 4x44K design, so the feature
identifier is a bare integer and not an `A_23_P…` probe name. This is stated because it
looks like a row index and is not one: the matrix `ID_REF` column carries the same integers,
all 45 220 of them distinct, in the same order as the platform table.

| gene | features | rule |
|---|---|---|
| **SIGLEC15** | **1** — feature `31474` | single feature; no collapse |
| **ACP5** | **1** — feature `10450` | single feature; no collapse |
| **MMP9** | **10** — features `1570`, `1678`, `14430`, `14852`, `21363`, `25413`, `28592`, `28675`, `40398`, `44612` | **arithmetic mean** over the ten |

**How the multi-probe situation was resolved, without looking at the correlation.** The
existing cohorts entered the analysis through
`load_geo_expression()` in `analysis_F_cross_scale.py`:

```python
symbol = annotation_symbol.split("///")[0]
selected.groupby("gene").mean(numeric_only=True)
```

i.e. **the arithmetic mean over every probe carrying that symbol**. The identical rule is
applied to GSE21501, so the cohort reaches the programme score on the same gene-level
footing as GSE57495 and GSE85916. The rule was fixed from the existing analysis identity
and from the platform annotation; **no expression value entered the decision**, and the
choice of MMP9's ten features was not informed by how they behave.

**Two robustness checks, both expression-free:**

* **`///` multi-symbol check.** 116 platform rows carry `///` in the symbol cell. **Zero**
  target genes appear as a non-first symbol, so no target feature is lost by the
  `split("///")[0]` rule.
* **alias check by `Gene title`.** Every feature whose title matches the target gene's
  Entrez name also carries the target symbol, and no extra feature does —
  SIGLEC15 1/1, ACP5 1/1, MMP9 10/10.

**Usability.** All 45 220 data rows × 132 sample cells contain a token: **0 rows contain an
empty cell and 0 rows contain `null`** (whole-matrix value-blind scan; detector self-tested).
Every target feature is therefore present and complete for every array, and
`usable_siglec15_probe`, `usable_acp5_probe` and `usable_mmp9_probe` are `yes` for all 132
rows. **Gene usability excludes nothing.**

---

## 8. Representation — AUTHOR DECISION, FLAGGED

> **AD-2 — frozen primary representation: the deposited series-matrix `VALUE` as supplied,
> i.e. Lowess-normalised `log2(Cy5 / Cy3)`, with NO further transform.**

**What the platform supplies.** One quantity per feature per array:
`log2(tumour intensity / Stratagene reference intensity)`, Lowess-normalised (§1.1). It is
already on a **log2** scale, and it is a **ratio**, not an intensity.

**Why this is the representation chosen, on analysis-identity grounds.** Every existing
cohort entered `score_program()` on a **log2 scale** — TCGA-PAAD as `log2(CPM + 1)`,
GSE57495 and GSE85916 as log2 microarray intensities — and the mean-z rule was developed and
applied in that input domain (§9). The deposit supplies no second representation, so the
choice is not *which matrix* (as it was for GSE205154, which supplied Counts, TPM and TMM)
but **whether to leave the supplied log2 quantity alone or transform it**. Leaving it alone
keeps the input domain identical to every comparator. Applying an inverse log (i.e. using
the raw ratio) would move the score onto a linear, outlier-dominated scale that no
comparator cohort uses.

**Why this is flagged rather than settled.** The frozen text says bulk expression is "used
as supplied" and is silent about the two-colour ratio case. Two consequences deserve the
author's explicit attention, and neither is resolvable from the frozen text:

1. **A log2 ratio is a *relative* measurement.** SIGLEC15's value for an array is its
   abundance relative to the *same* Stratagene lot, not an absolute level. For a
   within-cohort Spearman ρ and a within-cohort z-score this is admissible — both are
   invariant to any strictly monotone per-array rescaling, and the reference is identical
   across all 132 arrays, so no array is measured against a different yardstick. But the
   resulting SIGLEC15 value is **not on the same absolute scale** as `log2(CPM + 1)` in
   TCGA-PAAD, and no cross-cohort absolute comparison may be made from it. The frozen
   Platform condition explicitly does not require such comparability, which is why this is
   acceptable — but it must be stated, not assumed.
2. **The scope of the decision is asymmetric, and this bounds its risk.** It can move the
   **ACP5/MMP9 score** side only. SIGLEC15 enters the Spearman correlation through its
   **ranks**, which are invariant under any strictly monotone transform of the same matrix,
   so no part of AD-2 touches the SIGLEC15 side.

**Pre-specified sensitivity alternatives (must be reported whatever they show, and must
never replace the primary):**

* **S-rep-1** — build the ACP5/MMP9 score from **within-cohort ranks** of each gene instead
  of z-scores, testing the score's sensitivity to the standardisation's outlier behaviour.
* **S-rep-2** — see **S5** in §10: rebuild the score from each of MMP9's ten features
  individually. SIGLEC15 and ACP5 have one feature each, so no probe-choice sensitivity
  exists for them.
* There is **no** "no-log-transform" alternative to offer, because the deposit supplies no
  untransformed quantity. That absence is itself worth recording: unlike GSE205154, this
  cohort cannot be re-run on a second representation without going back to raw feature
  intensities, which are not deposited.

---

## 9. Frozen downstream rule that the ρ step must obey

These are **not** choices made here; they are transcribed from the existing analysis
identity so that the later step cannot silently re-specify them.

**Score construction — `analysis_F_cross_scale.py`, functions `zscore_columns()` and
`score_program()`:**

```
PROGRAMS["Derosiers_ACP5_MMP9_core"] = ["ACP5", "MMP9"]
required_gene_count(PRIMARY)         = 2          # both genes must be present
z                                    = (x - mean(x)) / sd(x)   over SAMPLES, ddof = 0
score                                = nanmean(z[ACP5], z[MMP9])   # the "mean-z" method
```

The standardisation is **within-cohort** (each gene is centred and scaled across the
cohort's samples), which is what makes the score comparable across cohorts without
requiring cross-platform absolute comparability.

**Statistic:** Spearman ρ between tumour-level SIGLEC15 and the score, computed
**separately per cohort**.

**Combination — `random_effects_meta()`:** Fisher z = `arctanh(clip(ρ, ±0.999999))`;
`v = 1/(n − 3)`; τ² by REML; modified Hartung–Knapp interval with the scale factor floored
at 1 (`q_HK = max(1, Σw(y−μ)²/(k−1))`, `se = sqrt(q_HK/Σw)`, t critical with `df = k − 1`);
back-transform by `tanh`. Seed 15072026, 10 000 resamples where resampling is used.

**Compartments:** GSE21501 contributes **one** cohort-level estimate (primary tumours).
There is no metastatic compartment in this deposit, so no second estimate is generated here
and none may be manufactured from the Nebraska arrays.

---

## 10. FINAL n — the numbers this manifest fixes

### 10.1 Arithmetic

```
132  arrays deposited
       (GSM535948..GSM536949; !Series_sample_id lists 132; the series matrix has
        132 sample columns; every row has !Sample_data_row_count = 45220)
 -0  arrays lost to two-channel resolution   all 132 resolved to the channel-2 specimen
 -0  non-tumour arrays                       channel-2 'tissue' = Patient Primary
                                             Pancreatic Tumor for 132/132
 -0  metastatic arrays                       no channel-2 field declares a metastasis
===
132  primary pancreatic tumour arrays
 -30  Nebraska 15 + UNC 15                   no clinical annotation (deposit: "were not
                                             analyzed with clinical data so it is not
                                             provided"); 12 of 15 Nebraska arrays declare
                                             a mixed tissue composition and one declares
                                             nodal tumour; eligibility for these 30 is
                                             group-level only, never per array
===
102  clinically annotated arrays
 -0  duplicate arrays (DUP-2)                each of the 102 has its own specimen key:
                                             34 JH + 19 NS + 49 NW = 102 keys
 -0  arrays lost to gene usability           SIGLEC15, ACP5 and MMP9 all mapped and
                                             complete for all 132 arrays
===
102  FINAL n for the main meta-analysis

  For completeness, the same duplicate rule applied to the whole deposit:
132  arrays  ->  129 distinct specimen keys
 -3  duplicate arrays: 12R2587.2 (GSM535950), 13R2778Rpt (GSM536033),
     9R1701B (GSM536041) -- all three inside the excluded Nebraska group
===
129  final n for prespecified sensitivity S1
```

### 10.2 Frozen totals

| quantity | frozen value |
|---|---|
| **final n, main meta-analysis** | **102** |
| **final n, prespecified sensitivity S1** (whole deposit after DUP-2) | **129** |
| S2 (drop the two undocumented-suffix specimens, `GSM536488`, `GSM536901`) | 100 |
| S3 (drop the six `redo`/`Replicate` arrays of the 102) | 96 |
| arrays eligible for nothing in the deposit | 3 (`GSM535950`, `GSM536033`, `GSM536041`) |
| distinct specimen keys, main analysis | 102 |
| distinct specimen keys, whole deposit | 129 |
| **patients, main analysis** | **≤ 102 — NOT an exact count** |
| **patients, whole deposit** | **≤ 129 — NOT an exact count** |
| specimens excluded from every analysis | 0 — every one of the 132 is either in the main analysis or in S1 |

**No pooled or meta-analytic sample size is stated here.** What the other cohorts
contribute, and what any combined estimate would rest on, is outside what this manifest
establishes and is deliberately not pre-written.

### 10.3 Prespecified sensitivity analyses

| id | n | what it varies | why |
|---|---|---|---|
| **S1** | 129 | admits the 30 Nebraska/UNC arrays after DUP-2 | the frozen text does not require clinical annotation; this shows what the choice in AD-1 costs |
| **S2** | 100 | drops `GSM536488` (`P271N Replicate`) and `GSM536901` (`P188T IPMN`) | the two specimens whose histology/purity is not independently corroborated (AD-4) |
| **S3** | 96 | drops all six `redo`/`Replicate` arrays of the 102 | tests whether submitter-flagged re-runs behave differently (DUP-6) |
| **S4** | 129 | the duplicate pair **averaged** at expression level instead of selected | DUP-6; no pair exists inside the main 102, so S4 acts on S1's population |
| **S5** | 102 | score rebuilt from each of MMP9's ten features individually | the only multi-feature gene; report the range, not a chosen feature |

None of S1–S5 may replace the primary result.

---

## 11. AUTHOR DECISIONS — flagged, not buried

The frozen text did not settle the following. Each was decided explicitly here, before any
effect was examined, and each must be reported as an author decision rather than presented
as a pre-registered rule.

### AD-1 — Main-analysis eligibility requires clinical annotation (adopted)

*Decision:* the main meta-analysis is the **102** clinically annotated arrays; the 30
Nebraska/UNC arrays are excluded from it and carried as prespecified sensitivity **S1
(n = 129)**.
*Reason:* §5.3 — the deposit's own partition, per-array rather than group-level eligibility,
and direct per-sample contradiction inside the Nebraska branch.
*Consequence if the author disagrees:* the main analysis becomes **n = 129** and the 102
becomes the sensitivity. Both numbers are computed from the same frozen duplicate rule, so
the switch is a one-line change of `eligible_for_main_meta` — but it must be made by the
author, and the frozen manifest hash above covers **only** the 102 choice.
*Why this needed deciding:* the frozen eligibility table has no annotation condition; the
estimand needs none. This is a judgement about what "tumour tissue" and "sufficient
metadata" mean per array.

### AD-2 — Frozen representation is the supplied log2 ratio, unchanged (adopted)

Set out in full in **§8**. Flagged because the frozen text is silent on two-colour ratios.

### AD-3 — Specimen key is the provider-scoped title-derived code, and is reported as a bound (adopted)

*Decision:* §3.2. GEO exposes no patient identifier and no BioSample, so the title-derived
provider-scoped specimen code is the key, recorded in `patient_id_if_available` with a
`not_exposed_by_GEO;` prefix.
*Reason:* it is the only patient-linked handle that exists, and the provider prefix prevents
five unrelated namespaces from being silently merged.
*Consequence:* the unit is the **specimen**, which is at best a proxy for the patient.
**102 arrays, ≤ 102 patients.** No exact patient count is claimed.

### AD-4 — `P271N Replicate` and `P188T IPMN` are retained by the frozen rule, and removed in S2 (adopted)

*Decision:* both stay in the main analysis, on the strength of the deposited channel-2
`tissue` field and their full clinical records; both are removed in **S2 (n = 100)**.
*Reason:* the frozen rule keys on the deposited compartment field, not on the title suffix;
but the trailing `N` is undocumented and could denote normal tissue, and IPMN is a precursor
lesion rather than invasive PDAC. The publication's own cohort totals make one of these two
the most likely explanation of the 49-array / 48-patient Northwestern discrepancy (§3.3).
*Consequence if the author disagrees:* drop either or both from the primary, giving
n = 101 or n = 100.

### AD-5 — Four arrays with an incomplete clinical record are retained (adopted)

*Decision:* `GSM536940 PE621T` (no N stage — the array the deposit's own design note refers
to), `GSM536892 P152T`, `GSM536924 P35T`, `GSM536927 P42T` (no T stage) stay in the main
analysis.
*Reason:* the frozen conditions require metadata sufficient to place the specimen in a
compartment, which channel-2 `tissue` supplies; they do not require a complete clinical
record, and no clinical covariate enters the estimand.
*Consequence:* none for the arithmetic; recorded so the author can see that five of the 102
have at least one missing clinical field.

### AD-6 — The 30 are excluded from the main analysis but not from the record (adopted)

*Decision:* every one of the 30 carries an explicit `eligible_for_main_meta = no` **and** a
written `exclusion_reason`, and all 30 remain runnable as S1.
*Reason:* "excluded" must not mean "silently dropped"; the frozen text allows prespecified
sensitivity analyses and S1 is exactly the test of AD-1.

---

## 12. What could not be determined

1. **Any patient identifier.** GEO exposes none, and there is no BioSample accession to
   query. The publications publish only aggregate tables with no GSM mapping. The patient
   count is therefore reported as a bound throughout, never as a number.
2. **Whether the `N` in `P271N Replicate` denotes normal tissue.** Undocumented in both the
   deposit and the publication (AD-4).
3. **Whether `P188T IPMN` is invasive PDAC or a pure precursor lesion.** Undocumented; the
   publication never mentions IPMN (AD-4).
4. **Which of the 49 Northwestern arrays is the extra one** relative to the publication's 48
   Northwestern patients (§3.3). Two candidates, not decidable from public data.
5. **What the percentages in the Nebraska `Mix` titles mean** (`20/80/0`, `30/60/10`,
   `50/30/20`, `70/30/0`, `10/90/0`, `…Mix(90%N)…`). The publication does not define them,
   and they cannot be converted into a tumour-purity estimate.
6. **Why the 15 UNC1 and 15 NEB arrays have no clinical annotation in GEO** although
   Table 1 of the publication reports the UNC1 clinical characteristics in aggregate. The
   deposit's design line is the only statement, and it does not explain itself.
7. **How many probes per gene the original analysis used.** The publication reports neither
   a probe count nor a collapse rule; the rule applied here was taken from the manuscript's
   own frozen analysis identity (§7).
8. **Patient-level overlap with TCGA-PAAD and across the five providers.** No shared
   accession namespace or naming namespace exists; not assessable from public metadata.
9. **`GSE21501_family.soft.gz`** (187 687 633 bytes) was not obtained — the single-connection
   transfer stalled. Nothing in this directory depends on it; see `failed_or_excluded.md`.
10. **A representation-free re-analysis.** Unlike GSE205154, this deposit supplies only one
    quantity (the log2 ratio), so there is no second matrix to re-run on without raw feature
    intensities, which are not deposited (§8).

---

## 13. Attestation

* The manifest `sample_manifest.csv` was built by `scripts/07_build_manifest.py` from
  metadata, identifiers and platform annotation only.
* Every one of the **132** rows carries an explicit `yes`/`no` on `eligible_for_main_meta`;
  the build asserts this, and asserts that every `no` carries a non-empty
  `exclusion_reason`. **No sample is silently skipped.**
* The eligibility tally is exactly `yes` = 102 and `no` = 30, and the 30 are exactly the
  Nebraska 15 + UNC 15.
* All sample-level metadata was taken from **channel-2** fields. No eligibility decision
  reads a channel-1 field.
* **No Spearman correlation between SIGLEC15 and any ACP5/MMP9 programme score was
  computed, inspected or reported. No per-sample expression value of SIGLEC15, ACP5 or
  MMP9 was loaded.** The value-blind guard and its self-test are described in §0 and in
  `scripts/05_probe_mapping.py`.
* The SHA-256 recorded at the top of this document is the hash of the frozen manifest. Any
  later edit to the manifest invalidates this definition and must be recorded as a new
  version with its own hash.
