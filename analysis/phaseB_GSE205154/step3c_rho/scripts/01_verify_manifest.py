#!/usr/bin/env python
"""
01_verify_manifest.py -- Step 3c-1, GSE205154 only.

Purpose
-------
Re-verify the FROZEN cohort manifest before any expression value is touched, and
derive the four sample sets this step needs, using ONLY the frozen columns:

  * MAIN            n = 216   primary tumours, one specimen per patient (D4.2 selection)
  * MET_FULL        n =  71   metastatic lesions, one specimen per patient
  * MET_OVERLAPFREE n =  61   the 71 minus patients that also contribute a primary
  * DUP_PAIRS       2 pairs   the two within-compartment duplicate primary cores (D4.5)

This script reads NO expression value. It is manifest-only by construction: the only
input is sample_manifest.csv, and the only columns it reads are the accession /
title / compartment / eligibility columns.

Guardrails asserted here (fail loudly, never silently):
  G1  SHA-256 of sample_manifest.csv equals the value frozen in frozen_cohort_definition.md
  G2  exactly 289 data rows
  G3  eligible_for_main_meta == yes  ->  exactly 216
  G4  eligible_for_sensitivity == yes ->  exactly  71
  G5  216 main specimens -> EXACTLY 216 distinct patients (AD-1 enforced by construction;
      this is the guard that makes the main analysis patient-level, and it is the one that
      would fire if the manifest's flags stopped encoding the D4.2 selection)
  G5b the 218 primary specimens sit in 216 distinct patients (214 patients with exactly
      one primary specimen + 2 patients with two), and exactly two patients contribute
      more than one primary specimen -- the D4.2 duplicate cores
  G6  those two duplicate primary cores are exactly ST-00017945 and ST-00019598, and in
      each the base-label specimen is the one flagged into the main set
  G7  the 71 metastatic specimens sit in 71 distinct cores (no metastatic duplication)
  G8  the overlap between main patients and met patients is exactly 10
        -> MET_OVERLAPFREE n = 61
  G9  the two specimens excluded from everything are GSM6206193 and GSM6206245
  G10 every 'no' carries a non-empty exclusion_reason (nothing silently skipped)
  G11 the three gene-usability columns are 'yes' for all 289 rows
  G12 the platform string is GPL20301 / exome capture, for all 289 rows (AD/platform disclosure)
  G13 no expression-looking column is read (schema assertion: the manifest has 13 columns
      and none of them is a gene symbol)

Outputs (to step3c_rho/output/):
  manifest_audit.txt
  sample_sets.json          <- the frozen sample sets, consumed by every later script
  manifest_verification.csv <- one row per sample with its set membership

Run:
  python 01_verify_manifest.py
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from datetime import datetime
from pathlib import Path

# ----------------------------------------------------------------------------- paths
HERE = Path(__file__).resolve().parent
STEP = HERE.parent                       # .../phaseB_GSE205154/step3c_rho
PHASEB = STEP.parent                     # .../phaseB_GSE205154
ROOT = PHASEB.parent                     # .../BJC_EXT_v1_20260927

MANIFEST = PHASEB / "sample_manifest.csv"
FROZEN_DEF = PHASEB / "frozen_cohort_definition.md"
OUT = STEP / "output"
OUT.mkdir(parents=True, exist_ok=True)

EXPECTED_SHA256 = "f76b81eceb0557cf4c1fee0138e9137a84508369f47ad163255a096af78c708b"

# The frozen cohort definition fixes these and nothing in this step may change them.
EXPECT_MAIN_N = 216
EXPECT_MET_N = 71
EXPECT_MET_OVERLAPFREE_N = 61
EXPECT_ROWS = 289
EXPECT_PRIMARY_SPECIMENS = 218
EXPECT_PRIMARY_DISTINCT_PATIENTS = 216
EXPECT_DUP_CORES = ["ST-00017945", "ST-00019598"]
EXPECT_EXCLUDED_ALL = ["GSM6206193", "GSM6206245"]
EXPECT_CROSS_COMPARTMENT_PATIENTS = 10

MANIFEST_COLUMNS = [
    "gsm_accession",
    "sample_title",
    "tissue_type",
    "tumour_or_non_tumour",
    "primary_or_metastatic",
    "patient_id_if_available",
    "platform_or_library",
    "usable_siglec15_probe_or_gene",
    "usable_acp5_probe_or_gene",
    "usable_mmp9_probe_or_gene",
    "eligible_for_main_meta",
    "eligible_for_sensitivity",
    "exclusion_reason",
]

log_lines: list[str] = []
failures: list[str] = []


def log(text: str = "") -> None:
    log_lines.append(text)
    print(text)


def check(cond: bool, ok_msg: str, fail_msg: str, *, hard: bool = True) -> bool:
    if cond:
        log(f"  [PASS] {ok_msg}")
        return True
    log(f"  [FAIL] {fail_msg}")
    failures.append(fail_msg)
    if hard:
        return False
    return False


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def core_of(title: str) -> str:
    """Frozen patient key: the ST-######## core preceding the specimen suffix."""
    m = re.match(r"^(ST-\d+)-", title.strip())
    if not m:
        raise ValueError(f"title does not match ^(ST-\\d+)- : {title!r}")
    return m.group(1)


def suffix_of(title: str) -> str:
    m = re.match(r"^ST-\d+-(.+)$", title.strip())
    return m.group(1) if m else ""


# ============================================================== begin
log("=" * 78)
log("STEP 3c-1  GSE205154  --  MANIFEST RE-VERIFICATION (value-blind)")
log(f"generated : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 78)
log()
log("Governing rule obeyed by this script: the cohort is re-verified from the frozen")
log("manifest BEFORE any SIGLEC15 / ACP5 / MMP9 expression value is loaded. This script")
log("opens exactly one data file -- sample_manifest.csv -- and it contains only sample")
log("metadata and yes/no eligibility flags, never an expression measurement.")
log()

# ---------------------------------------------------------------- G1 hash
log("--- G1  SHA-256 of the frozen manifest ---")
log(f"  file     : {MANIFEST}")
actual_sha = sha256_file(MANIFEST)
log(f"  computed : {actual_sha}")
log(f"  expected : {EXPECTED_SHA256}")
log(f"  source   : frozen_cohort_definition.md, top of file")
hash_ok = actual_sha == EXPECTED_SHA256
check(hash_ok, "SHA-256 MATCHES the frozen value",
      f"SHA-256 MISMATCH -- manifest edited after freeze; estimate not covered by the definition")

# the definition text must actually contain that hash (no silent redefinition)
def_text = FROZEN_DEF.read_text(encoding="utf-8")
check(EXPECTED_SHA256 in def_text,
      "the same hash literal appears in frozen_cohort_definition.md",
      "hash literal absent from frozen_cohort_definition.md")

log()
log("--- G13  schema assertion: no expression column was read ---")
with MANIFEST.open("r", encoding="utf-8", newline="") as fh:
    reader = csv.DictReader(fh)
    fieldnames = list(reader.fieldnames or [])
    rows = list(reader)
check(fieldnames == MANIFEST_COLUMNS,
      f"columns are exactly the 13 frozen metadata columns",
      f"unexpected columns: {fieldnames}")
gene_like = [c for c in fieldnames if re.search(r"siglec|acp5|mmp9", c, re.I)
             and c not in ("usable_siglec15_probe_or_gene", "usable_acp5_probe_or_gene",
                           "usable_mmp9_probe_or_gene")]
check(not gene_like,
      "no column holds a SIGLEC15/ACP5/MMP9 measurement (only usability flags)",
      f"expression-like columns present: {gene_like}")

log()
log("--- G2  row count ---")
check(len(rows) == EXPECT_ROWS, f"{len(rows)} data rows", f"{len(rows)} rows, expected {EXPECT_ROWS}")
log(f"  gsm distinct      : {len({r['gsm_accession'] for r in rows})}")
log(f"  title distinct    : {len({r['sample_title'] for r in rows})}")

# ---------------------------------------------------------------- core / eligibility
for r in rows:
    r["_core"] = core_of(r["sample_title"])
    r["_suffix"] = suffix_of(r["sample_title"])
    r["_main"] = r["eligible_for_main_meta"].strip().lower() == "yes"
    r["_sens"] = r["eligible_for_sensitivity"].strip().lower() == "yes"

log()
log("--- G11  gene usability flags (all 289) ---")
for col, gname in (("usable_siglec15_probe_or_gene", "SIGLEC15"),
                   ("usable_acp5_probe_or_gene", "ACP5"),
                   ("usable_mmp9_probe_or_gene", "MMP9")):
    n_yes = sum(1 for r in rows if r[col].strip().lower() == "yes")
    check(n_yes == EXPECT_ROWS, f"{gname}: usable = yes for {n_yes}/{EXPECT_ROWS}",
          f"{gname}: only {n_yes}/{EXPECT_ROWS} usable")

log()
log("--- G12  platform disclosure string (all 289 must be identical) ---")
plats = sorted({r["platform_or_library"] for r in rows})
for p in plats:
    log(f"  {p}")
_plat_l = plats[0].lower() if len(plats) == 1 else ""
check(len(plats) == 1 and "gpl20301" in _plat_l and "exome" in _plat_l
      and "hybrid capture" in _plat_l,
      "single platform string: GPL20301, explicitly exome-capture "
      "(FFPE xGEN Exome Research Panel v1.0 hybrid capture), not whole-transcriptome",
      f"platform string not the expected GPL20301 exome-capture string: {plats}")

log()
log("--- G10  every 'no' carries an exclusion_reason ---")
bad_reason = [r["gsm_accession"] for r in rows
              if (not r["_main"] or not r["_sens"]) and not r["exclusion_reason"].strip()]
check(not bad_reason, "every non-eligible row has a written exclusion_reason",
      f"{len(bad_reason)} non-eligible rows without a reason: {bad_reason[:5]}")

# ---------------------------------------------------------------- the four sets
main = [r for r in rows if r["_main"]]
met = [r for r in rows if r["_sens"]]
neither = [r for r in rows if not r["_main"] and not r["_sens"]]
both = [r for r in rows if r["_main"] and r["_sens"]]

log()
log("--- G3/G4  eligibility counts ---")
check(len(main) == EXPECT_MAIN_N, f"eligible_for_main_meta : {len(main)}",
      f"eligible_for_main_meta : {len(main)}, expected {EXPECT_MAIN_N}")
check(len(met) == EXPECT_MET_N, f"eligible_for_sensitivity: {len(met)}",
      f"eligible_for_sensitivity: {len(met)}, expected {EXPECT_MET_N}")
check(len(both) == 0, "no sample is eligible for both analyses (frozen: 0)",
      f"{len(both)} samples eligible for both")
check(len(neither) == 2, f"eligible for nothing       : {len(neither)}",
      f"eligible for nothing: {len(neither)}, expected 2")

from collections import Counter

log()
log("--- G5  main set is patient-level: 216 specimens -> 216 distinct patients ---")
main_cores = sorted({r["_core"] for r in main})
check(len(main_cores) == EXPECT_MAIN_N,
      f"main: {len(main)} specimens in {len(main_cores)} distinct patients "
      f"-> one specimen per patient, AD-1 enforced",
      f"main: {len(main)} specimens in {len(main_cores)} distinct patients "
      f"(expected {EXPECT_MAIN_N}); the manifest flags do not encode the D4.2 selection")
# the strongest form of the same guard: no patient appears twice in the main set
main_core_counts_chk = Counter(r["_core"] for r in main)
_repeat = sorted(c for c, k in main_core_counts_chk.items() if k > 1)
check(not _repeat,
      "no patient contributes more than one specimen to the main analysis",
      f"patients with >1 main specimen: {_repeat}")

# G5b -- the D4.2 duplicate cores are visible in the 218 PRIMARY SPECIMENS, not in the
# 216-row eligible set (which has already had them removed by the manifest's own flags).
prim = [r for r in rows if r["primary_or_metastatic"].strip().lower() == "primary"]
prim_cores = sorted({r["_core"] for r in prim})
log()
log("--- G5b  the 218 primary specimens sit in 216 distinct patients (D4.2 material) ---")
check(len(prim) == EXPECT_PRIMARY_SPECIMENS,
      f"{len(prim)} primary specimens (AD-4: specimens)",
      f"{len(prim)} primary specimens, expected {EXPECT_PRIMARY_SPECIMENS}")
check(len(prim_cores) == EXPECT_PRIMARY_DISTINCT_PATIENTS,
      f"{len(prim)} primary specimens in {len(prim_cores)} distinct patients "
      f"({len(prim_cores) - 2} with one primary specimen + 2 with two)",
      f"{len(prim)} primary specimens in {len(prim_cores)} patients, "
      f"expected {EXPECT_PRIMARY_DISTINCT_PATIENTS}")

# which primary patients contribute >1 primary specimen
from collections import Counter
prim_core_counts = Counter(r["_core"] for r in prim)
dup_cores = sorted(c for c, k in prim_core_counts.items() if k > 1)
log()
log("--- G6  duplicate primary cores (D4.2): frozen rule reproduced ---")
check(dup_cores == EXPECT_DUP_CORES,
      f"exactly the two frozen duplicate cores: {dup_cores}",
      f"duplicate cores {dup_cores}, expected {EXPECT_DUP_CORES}")

dup_pairs = []
for c in dup_cores:
    members = sorted([r for r in prim if r["_core"] == c], key=lambda r: r["sample_title"])
    retained = [r for r in members if r["_suffix"] == "T"]
    log(f"  core {c}:")
    for r in members:
        mark = "RETAINED (base label -T)" if r["_suffix"] == "T" else "excluded from main by D4.2"
        log(f"      {r['gsm_accession']}  {r['sample_title']:<20} suffix={r['_suffix']:<4} {mark}")
    if len(retained) != 1:
        failures.append(f"core {c}: {len(retained)} base-label specimens, expected exactly 1")
        log(f"  [FAIL] core {c}: {len(retained)} base-label specimens")
        continue
    excluded = [r for r in members if r is not retained[0]]
    # the frozen rule says the BASE LABEL specimen is retained AND is the one in main
    check(retained[0]["_main"] and not excluded[0]["_main"],
          f"core {c}: base-label {retained[0]['sample_title']} is the one flagged eligible_for_main_meta",
          f"core {c}: D4.2 base-label rule and the manifest flag disagree")
    dup_pairs.append({
        "core": c,
        "retained_title": retained[0]["sample_title"],
        "retained_gsm": retained[0]["gsm_accession"],
        "excluded_title": excluded[0]["sample_title"],
        "excluded_gsm": excluded[0]["gsm_accession"],
    })

# G5c -- an M2 label whose core differs from the M specimen's core is NOT a duplicate
log()
log("--- G5c  replicate-looking titles that are NOT within-compartment duplicates ---")
for r in rows:
    if r["_suffix"] in ("T2", "M2"):
        same_core = [x for x in rows
                     if x["_core"] == r["_core"] and x["_suffix"] == r["_suffix"][0]]
        log(f"  {r['gsm_accession']}  {r['sample_title']:<20} compartment={r['primary_or_metastatic']:<11} "
            f"base-label specimen in SAME core: "
            f"{same_core[0]['sample_title'] if same_core else 'NONE'}")
        if r["_suffix"] == "M2" and not same_core:
            check(True,
                  f"{r['sample_title']}: no -M specimen in core {r['_core']} -> an M2-labelled "
                  f"metastasis from a distinct patient, NOT a within-compartment duplicate",
                  "")

log()
log("--- G9  specimens eligible for nothing ---")
neither_gsm = sorted(r["gsm_accession"] for r in neither)
check(neither_gsm == EXPECT_EXCLUDED_ALL,
      f"exactly {EXPECT_EXCLUDED_ALL}",
      f"unexpected: {neither_gsm}")
for r in neither:
    log(f"  {r['gsm_accession']}  {r['sample_title']:<20} reason: {r['exclusion_reason']}")

log()
log("--- G7  metastatic set: 71 specimens in 71 distinct cores ---")
met_cores = sorted({r["_core"] for r in met})
met_core_counts = Counter(r["_core"] for r in met)
met_dup = sorted(c for c, k in met_core_counts.items() if k > 1)
check(len(met_cores) == EXPECT_MET_N and not met_dup,
      f"met: {len(met)} specimens in {len(met_cores)} distinct cores; duplication = none",
      f"met: {len(met_cores)} cores, duplicates {met_dup}")

log()
log("--- G8  cross-compartment patient overlap (AD-4) ---")
overlap = sorted(set(r["_core"] for r in main) & set(r["_core"] for r in met))
check(len(overlap) == EXPECT_CROSS_COMPARTMENT_PATIENTS,
      f"{len(overlap)} patients contribute both a primary and a metastasis",
      f"{len(overlap)} overlapping patients, expected {EXPECT_CROSS_COMPARTMENT_PATIENTS}")
for c in overlap:
    m = [r for r in main if r["_core"] == c][0]
    s = [r for r in met if r["_core"] == c][0]
    log(f"  {c}   primary {m['gsm_accession']} {m['sample_title']:<22} "
        f"met {s['gsm_accession']} {s['sample_title']}")

met_overlapfree = [r for r in met if r["_core"] not in set(overlap)]
check(len(met_overlapfree) == EXPECT_MET_OVERLAPFREE_N,
      f"overlap-free metastatic subset: {len(met_overlapfree)}",
      f"overlap-free subset {len(met_overlapfree)}, expected {EXPECT_MET_OVERLAPFREE_N}")

# ---------------------------------------------------------------- AD-5 disclosure
log()
log("--- AD-5  F-suffix samples are primary (counted, not inferred biologically) ---")
f_samples = [r for r in rows if r["_suffix"].startswith("F")]
f_primary = [r for r in f_samples if r["_main"]]
log(f"  {len(f_samples)} F-suffix samples; {len(f_primary)} of them eligible_for_main_meta")
for r in f_samples:
    log(f"      {r['gsm_accession']}  {r['sample_title']:<20} tumor_type_field={r['primary_or_metastatic']:<8} "
        f"tissue={r['tissue_type']:<10} main={r['eligible_for_main_meta']}")
check(len(f_samples) == 7, f"7 F-suffix samples", f"{len(f_samples)} F-suffix samples, expected 7")
check(all(r["_main"] for r in f_samples), "all F-suffix samples classified primary by tumor type field",
      "some F-suffix sample is not classified primary")

log()
log("--- AD-1 / compartment tally ---")
comp = Counter(r["primary_or_metastatic"].strip().lower() for r in rows)
log(f"  {dict(comp)}")
check(comp.get("primary") == 218 and comp.get("metastatic") == 71,
      "218 primary specimens / 71 metastatic specimens (AD-4: specimens, not patients)",
      f"compartment tally unexpected: {dict(comp)}")

# ---------------------------------------------------------------- write
sample_sets = {
    "generated": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    "manifest": str(MANIFEST),
    "manifest_sha256": actual_sha,
    "manifest_sha256_expected": EXPECTED_SHA256,
    "manifest_sha256_match": bool(hash_ok),
    "main": {
        "label": "primary tumours, one specimen per patient (D4.2 selection)",
        "level": "patient",
        "n": len(main),
        "n_distinct_patients": len(main_cores),
        "gsm": [r["gsm_accession"] for r in main],
        "titles": [r["sample_title"] for r in main],
        "patients": sorted({r["_core"] for r in main}),
    },
    "met_full": {
        "label": "metastatic lesions, one specimen per patient",
        "level": "patient",
        "n": len(met),
        "n_distinct_patients": len(met_cores),
        "gsm": [r["gsm_accession"] for r in met],
        "titles": [r["sample_title"] for r in met],
        "patients": sorted({r["_core"] for r in met}),
    },
    "met_overlapfree": {
        "label": "metastatic lesions from patients with no primary in the main analysis",
        "level": "patient",
        "n": len(met_overlapfree),
        "n_distinct_patients": len({r["_core"] for r in met_overlapfree}),
        "gsm": [r["gsm_accession"] for r in met_overlapfree],
        "titles": [r["sample_title"] for r in met_overlapfree],
        "patients": sorted({r["_core"] for r in met_overlapfree}),
    },
    "primary_specimens_all_218": {
        "label": "all primary specimens BEFORE the D4.2 selection (the n=218 variant of AD-1)",
        "level": "specimen",
        "n": len(prim),
        "n_distinct_patients": len(prim_cores),
        "gsm": [r["gsm_accession"] for r in prim],
    },
    "dup_pairs": dup_pairs,
    "cross_compartment_cores": overlap,
    "excluded_from_everything": [
        {"gsm": r["gsm_accession"], "title": r["sample_title"], "reason": r["exclusion_reason"]}
        for r in neither
    ],
    "f_suffix_samples": [r["gsm_accession"] for r in f_samples],
    "all_manifest_titles": [r["sample_title"] for r in rows],
    "all_manifest_gsm": [r["gsm_accession"] for r in rows],
    "platform_or_library": plats[0] if len(plats) == 1 else plats,
    "distinct_patients_whole_deposit": len({r["_core"] for r in rows}),
}

(OUT / "sample_sets.json").write_text(json.dumps(sample_sets, indent=2), encoding="utf-8")
log()
log(f"wrote {OUT / 'sample_sets.json'}")

with (OUT / "manifest_verification.csv").open("w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["gsm_accession", "sample_title", "specimen_core", "specimen_suffix",
                "primary_or_metastatic", "eligible_for_main_meta", "eligible_for_sensitivity",
                "in_main_216", "in_met_71", "in_met_overlapfree_61",
                "excluded_by_D4.2", "tumour_or_non_tumour", "exclusion_reason"])
    excl_gsm = {d["excluded_gsm"] for d in dup_pairs}
    of_gsm = {r["gsm_accession"] for r in met_overlapfree}
    for r in sorted(rows, key=lambda r: r["gsm_accession"]):
        w.writerow([
            r["gsm_accession"], r["sample_title"], r["_core"], r["_suffix"],
            r["primary_or_metastatic"], r["eligible_for_main_meta"], r["eligible_for_sensitivity"],
            "yes" if r["_main"] else "no",
            "yes" if r["_sens"] else "no",
            "yes" if r["gsm_accession"] in of_gsm else "no",
            "yes" if r["gsm_accession"] in excl_gsm else "no",
            r["tumour_or_non_tumour"], r["exclusion_reason"],
        ])
log(f"wrote {OUT / 'manifest_verification.csv'}")

log()
log("=" * 78)
if failures:
    log(f"RESULT: {len(failures)} CHECK(S) FAILED")
    for f in failures:
        log(f"  - {f}")
else:
    log("RESULT: ALL CHECKS PASSED -- frozen cohort reproduced exactly")
    log(f"  G1  manifest SHA-256 MATCHES the frozen value")
    log(f"  G3  main               n = {len(main)}   patient-level "
        f"({len(main)} specimens -> {len(main_cores)} distinct patients)")
    log(f"  G4  met                n = {len(met)}    patient-level "
        f"({len(met)} specimens -> {len(met_cores)} distinct patients)")
    log(f"  G8  met (overlap-free) n = {len(met_overlapfree)}")
    log(f"  G9  excluded from everything: {neither_gsm}")
    log(f"  G12 platform           : {plats[0]}")
    log(f"  distinct patients, whole deposit: {len({r['_core'] for r in rows})}")
log("=" * 78)

(OUT / "manifest_audit.txt").write_text("\n".join(log_lines) + "\n", encoding="utf-8")
log(f"wrote {OUT / 'manifest_audit.txt'}")

sys.exit(1 if failures else 0)
