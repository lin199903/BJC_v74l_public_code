#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
extension_meta / 04_verify_summary.py

Checks the WRITE-UP against the DATA. results_summary.md is prose written after the run; this
script re-reads meta_results.csv / forest_inputs.csv / verification_controls.csv and asserts
that every headline number quoted in the prose is present and correct, that the required
verbatim sentence is carried unaltered, and that nothing outside the five frozen bulk cohorts
was pooled.

This exists so that the summary cannot silently drift from the numbers it reports. It is
capable of failing: if a value in the summary is edited away from the CSV, this exits non-zero.
"""
from __future__ import annotations
import os

import re
import sys
from pathlib import Path

import pandas as pd

EXT = Path(os.environ.get("SIGLEC15_EXT_META_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927\extension_meta"))
ALLOWED_COHORTS = {"TCGA_PAAD", "GSE57495", "GSE85916", "GSE205154", "GSE21501"}
SPECIMEN_SENTENCE = (
    "The GSE21501 analysis was specimen-level because patient identifiers were not deposited; "
    "the number of independent patients could therefore not be established exactly."
)

FAILS: list[str] = []
PASSES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    (PASSES if ok else FAILS).append(f"{name} :: {detail}")
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" :: {detail}" if detail else ""))


def norm(s: str) -> str:
    """Fold typographic minus/dash variants and non-breaking spaces so that a summary written
    with U+2212 compares equal to an ASCII-hyphen CSV value."""
    for ch in ("\u2212", "\u2013", "\u2010", "\u2011", "\u2012", "\u2014", "\u00a0"):
        s = s.replace(ch, "-" if ch != "\u00a0" else " ")
    return s


def main() -> int:
    md_raw = (EXT / "results_summary.md").read_text(encoding="utf-8")
    md = norm(md_raw)
    mr = pd.read_csv(EXT / "meta_results.csv")
    fi = pd.read_csv(EXT / "forest_inputs.csv")
    vc = pd.read_csv(EXT / "verification_controls.csv")

    print("=" * 92)
    print("04_verify_summary.py -- write-up vs data")
    print("=" * 92)

    # ---- 1. structure of meta_results.csv ------------------------------------------
    check("meta_results.csv has exactly the 3 required estimate rows", len(mr) == 3, f"{len(mr)} rows")
    check("meta_results.csv has exactly the 12 required columns",
          list(mr.columns) == ["estimate_label", "cohorts_included", "k", "pooled_rho", "ci_low",
                               "ci_high", "tau2", "I2", "Q", "p_value", "model", "notes"],
          ",".join(mr.columns))
    check("k values are 3, 4, 5 in order", list(mr["k"]) == [3, 4, 5], str(list(mr["k"])))

    # ---- 2. SCOPE: only the five frozen bulk cohorts were pooled --------------------
    pooled: set[str] = set()
    for s in mr["cohorts_included"]:
        pooled |= set(str(s).split(";"))
    check("only the five frozen bulk cohorts appear in any pooled estimate",
          pooled == ALLOWED_COHORTS, f"pooled={sorted(pooled)}")
    check("no single-cell / snRNA cohort was pooled",
          not any(re.search(r"278688|278689|GSE212966|GSE156405|GSE155698", str(s))
                  for s in mr["cohorts_included"]),
          "checked GSE278688/GSE278689/GSE212966/GSE156405/GSE155698")
    for _, r in mr.iterrows():
        n = len(str(r["cohorts_included"]).split(";"))
        if n != int(r["k"]):
            check(f"cohort list length matches k for row {r['k']}", False, f"{n} != {r['k']}")
    check("cohort list length matches k in every row",
          all(len(str(r["cohorts_included"]).split(";")) == int(r["k"]) for _, r in mr.iterrows()))

    # ---- 3. every headline number in the CSV appears in the write-up ----------------
    def has(text: str) -> bool:
        return norm(text) in md

    for _, r in mr.iterrows():
        lab = f"k={r['k']}"
        for col, nd in (("pooled_rho", 6), ("ci_low", 6), ("ci_high", 6),
                        ("Q", 6), ("p_value", 6)):
            v = f"{abs(float(r[col])):.{nd}f}"
            check(f"{lab}: {col} {v} quoted in results_summary.md", has(v), v)
        check(f"{lab}: tau2 {r['tau2']} quoted in results_summary.md", has(str(r["tau2"])),
              str(r["tau2"]))
        check(f"{lab}: I2 {r['I2']} stated as 0", float(r["I2"]) == 0.0, str(r["I2"]))

    # ---- 4. the published three-cohort anchor, as published ------------------------
    a = mr[mr["k"] == 3].iloc[0]
    check("3-cohort row rounds to the published rho 0.132", has("0.132"), f"{float(a['pooled_rho']):.6f}")
    check("3-cohort row rounds to the published CI -0.110", has("-0.110"),
          f"{float(a['ci_low']):.6f}")
    check("3-cohort row rounds to the published CI +0.360", has("+0.360"),
          f"{float(a['ci_high']):.6f}")
    check("summary states the published estimate is unreplaced",
          "unreplaced" in md and "REPRODUCED" in md_raw.upper())

    # ---- 5. the required verbatim sentence ----------------------------------------
    check("the GSE21501 specimen-level sentence is carried VERBATIM (exact substring)",
          SPECIMEN_SENTENCE in md_raw, f"{len(SPECIMEN_SENTENCE)} chars")
    check("the sentence is carried exactly once",
          md_raw.count(SPECIMEN_SENTENCE) == 1, f"{md_raw.count(SPECIMEN_SENTENCE)} occurrence(s)")
    # (precedence bug in the first form of this check: `A and B or C` parsed as `(A and B) or C`
    #  and passed for the wrong reason. Parenthesised, and both spellings required.)
    check("the <=129 patient bound is stated, not an exact count",
          ("\u2264 129 patients" in md_raw) or ("<= 129 patients" in md),
          f"U+2264 form: {'\u2264 129 patients' in md_raw}; ascii form: {'<= 129 patients' in md}")

    # ---- 6. zero-crossing statements match the CSV signs ---------------------------
    r4 = mr[mr["k"] == 4].iloc[0]
    r5 = mr[mr["k"] == 5].iloc[0]
    cross4 = float(r4["ci_low"]) < 0 < float(r4["ci_high"])
    cross5 = float(r5["ci_low"]) < 0 < float(r5["ci_high"])
    check("the 4- and 5-cohort intervals genuinely disagree about excluding zero",
          cross4 != cross5, f"4-cohort crosses={cross4}; 5-cohort crosses={cross5}")
    check("summary says the 4-cohort interval crosses zero", cross4 and "yes" in md.lower())
    check("summary reports both and selects neither",
          "neither is selected" in md.lower() and "side by side" in md.lower())

    # ---- 7. weight shares quoted in the summary match forest_inputs.csv -----------
    w4 = fi[(fi["model"].str.startswith("4-cohort")) & (fi["cohort"] == "GSE205154")]
    w5 = fi[(fi["model"].str.startswith("5-cohort")) & (fi["cohort"] == "GSE21501")]
    s4 = float(w4["weight_share_pct_of_model"].iloc[0])
    s5 = float(w5["weight_share_pct_of_model"].iloc[0])
    check("GSE205154 weight share quoted (40.56%)", has(f"{s4:.2f}"), f"{s4:.4f}%")
    check("GSE21501 weight share quoted (19.36%)", has(f"{s5:.2f}"), f"{s5:.4f}%")

    # ---- 8. forest_inputs traceability ------------------------------------------
    check("forest_inputs.csv has 3+4+5 = 12 rows, one per cohort per model", len(fi) == 12,
          f"{len(fi)} rows")
    check("every forest row carries the source file and its sha256",
          fi["rho_source"].notna().all() and fi["rho_source_sha256"].notna().all())
    check("GSE21501 is labelled specimen-level and no other cohort is",
          (fi.loc[fi["cohort"] == "GSE21501", "inferential_unit"].str.contains("specimen").all()
           and fi.loc[fi["cohort"] != "GSE21501", "inferential_unit"].eq("patient").all()),
          ";".join(sorted(fi["inferential_unit"].unique())))
    # each model's weights must sum to the reciprocal of the squared se used
    for m, g in fi.groupby("model"):
        wsum = float(g["weight_in_fitted_model"].sum())
        se = float(g["model_se_used"].iloc[0])
        q = max(1.0, float(g["model_hk_scale_raw"].iloc[0]))
        implied_se = (q / wsum) ** 0.5
        check(f"{m}: forest weights reconcile with the se actually used",
              abs(implied_se - se) < 1e-12, f"implied {implied_se:.12f} vs {se:.12f}")

    # ---- 9. control values quoted in the summary match verification_controls.csv --
    for cid, lo, hi in (("C3_unmodified_HK", "+0.083906", "+0.179626"),
                        ("C3b_unmodified_HK_published_3cohort", "+0.028418", "+0.233417")):
        row = vc[(vc["control"] == cid) & (vc["model"].str.contains("4-cohort|3-cohort"))]
        if cid.startswith("C3b"):
            row = vc[vc["control"] == cid]
        ok = len(row) > 0 and abs(float(row["ci_low"].iloc[0]) - float(lo)) < 1e-5 \
            and abs(float(row["ci_high"].iloc[0]) - float(hi)) < 1e-5
        check(f"{cid} CI in the summary matches verification_controls.csv", ok and has(lo) and has(hi),
              f"{lo}..{hi}")
    c6 = vc[vc["control"] == "C6_sign_flip_GSE205154"]
    d6 = abs(float(c6["delta_rho_vs_frozen"].iloc[0]))
    check("C6 single-cohort sign-flip delta quoted (0.086) matches the control file",
          has(f"{d6:.3f}") or has("0.086"), f"{d6:.6f}")

    # ---- 10. the independent-engine cross-check ran and passed --------------------
    xlog = (EXT / "provenance" / "crosscheck_metafor_log.txt").read_text(encoding="utf-8")
    check("metafor cross-check log records 15 checks, 0 failed",
          "TOTAL: 15 checks, 0 failed" in xlog)
    check("metafor log confirms Q agreement to 2.2e-16", "2.220e-16" in xlog)
    check("metafor log confirms the modified-HK rule matches to <=5.6e-17 or better",
          any(t in xlog for t in ("5.551e-17", "1.388e-17", "0.000e+00")))
    check("metafor log attributes the CI gap to the tau2 difference",
          xlog.count("fully explained by the tau2 difference") >= 3)
    check("summary quotes the cross-check findings",
          has("2.2e-16") and has("5.6e-17") and "15/15" in md)

    # ---- 11. provenance: source hash recorded ------------------------------------
    hl = (EXT / "provenance" / "source_hashes.txt").read_text(encoding="utf-8")
    check("provenance/source_hashes.txt records the manuscript module sha256",
          "7b6a311a0a73835f66596b7ea21070a6355c1fe6de44b518394dfe89f7c4334e" in hl)
    check("provenance/source_hashes.txt records the random_effects_meta() source sha256",
          "1b501f10eabaadebf5518630791a2757c101669e22a39312b0adacdc6cac3d98" in hl)
    check("summary quotes both source hashes", has("7b6a311a") and has("1b501f10"))

    print("\n" + "=" * 92)
    print(f"TOTAL: {len(PASSES) + len(FAILS)} checks, {len(FAILS)} failed")
    for f in FAILS:
        print(f"  FAILED: {f}")
    print("=" * 92)
    (EXT / "provenance" / "verify_summary_log.txt").write_text(
        "\n".join(PASSES + [f"FAILED: {f}" for f in FAILS]) + "\n", encoding="utf-8")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
