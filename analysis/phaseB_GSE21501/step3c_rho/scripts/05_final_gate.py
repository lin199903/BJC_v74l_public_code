#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE21501 -- Phase B step 3c-2, final coherence gate.

The two hand-written documents (output/results_summary.md, output/failed_or_excluded.md)
quote numbers that were computed elsewhere.  This script re-reads them and checks every
number it can mechanically extract against the machine-readable outputs, so a typo in
prose cannot survive to the deliverable.

It also asserts the task's explicit output contract:
    primary_result.csv       rho, n, CI, CI method named, n resamples, seed
    sensitivity_results.csv  one row per prespecified version, including the n=102
                             subset, S2, S3, S4 pair-averaged, S5 per-MMP9-feature,
                             and a rank-based score alternative
    run_log.txt, sessionInfo.txt
    results_summary.md       contains the statistical-unit sentence VERBATIM
    failed_or_excluded.md    exists

Exit code 0 = every check passed.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
STEP = HERE.parent
OUT = STEP / "output"

UNIT_SENTENCE = ("The GSE21501 analysis was specimen-level because patient identifiers "
                 "were not deposited; the number of independent patients could therefore "
                 "not be established exactly.")

problems: list[str] = []
passed = 0


def check(name: str, ok: bool, detail: str = "") -> None:
    global passed
    if ok:
        passed += 1
        print(f"[PASS] {name}" + (f" -- {detail}" if detail else ""))
    else:
        problems.append(name)
        print(f"[FAIL] {name} -- {detail}")


def main() -> int:
    prim = pd.read_csv(OUT / "primary_result.csv")
    sens = pd.read_csv(OUT / "sensitivity_results.csv")
    est = pd.read_csv(OUT / "estimate_table_full.csv")
    ctrl = pd.read_csv(OUT / "verification_controls.csv")
    manv = pd.read_csv(OUT / "manuscript_verification.csv")
    summary = (OUT / "results_summary.md").read_text(encoding="utf-8")
    failed = (OUT / "failed_or_excluded.md").read_text(encoding="utf-8")
    runlog = (OUT / "run_log.txt").read_text(encoding="utf-8")
    sess = (OUT / "sessionInfo.txt").read_text(encoding="utf-8")

    print("=" * 74)
    print("A. OUTPUT CONTRACT")
    print("=" * 74)

    # ---- primary_result.csv -------------------------------------------------
    p = prim.iloc[0]
    check("primary_result.csv has exactly one row", len(prim) == 1, f"{len(prim)} rows")
    for col in ("rho", "n", "ci_low", "ci_high", "ci_method", "ci_level",
                "n_resamples", "seed"):
        check(f"primary_result.csv carries {col!r}", col in prim.columns)
    check("primary n == 129", int(p["n"]) == 129, f"n = {int(p['n'])}")
    check("primary CI method is NAMED (not blank/placeholder)",
          isinstance(p["ci_method"], str) and len(str(p["ci_method"])) > 20,
          str(p["ci_method"])[:70])
    check("primary n_resamples == 10000", int(p["n_resamples"]) == 10000)
    check("primary seed == 15072026 (frozen S9)", int(p["seed"]) == 15072026)
    check("primary CI is ordered and brackets rho",
          p["ci_low"] < p["rho"] < p["ci_high"],
          f"[{p['ci_low']:+.6f}, {p['ci_high']:+.6f}] around {p['rho']:+.6f}")
    check("primary_result.csv states the patient figure as a bound",
          "<=129" in str(p["patients"]).replace("≤", "<=")
          or "129 eligible" in str(p["patients"]),
          str(p["patients"])[:80])
    check("primary_result.csv carries the unit sentence verbatim",
          UNIT_SENTENCE in str(p["unit_statement"]))

    # ---- sensitivity_results.csv -------------------------------------------
    ids = set(sens["sensitivity_id"])
    required_sens = {
        "S-CLIN": "the clinically annotated subset (n = 102)",
        "S2": "S2",
        "S3": "S3",
        "S4": "S4 (duplicate pairs averaged)",
        "S5": "S5 (per-MMP9-feature)",
        "S-rep-1": "a rank-based score alternative",
    }
    for sid, what in required_sens.items():
        check(f"sensitivity_results.csv contains {sid} ({what})", sid in ids)
    check("S-CLIN n == 102", int(sens.loc[sens.sensitivity_id.eq("S-CLIN"), "n"].iloc[0]) == 102)
    check("S3 n == 96 (frozen S10.3)", int(sens.loc[sens.sensitivity_id.eq("S3"), "n"].iloc[0]) == 96)
    check("S4 n == 129 (frozen S10.3)", int(sens.loc[sens.sensitivity_id.eq("S4"), "n"].iloc[0]) == 129)
    check("S5 n == 102 (frozen S10.3)", int(sens.loc[sens.sensitivity_id.eq("S5"), "n"].iloc[0]) == 102)
    check("S-rep-1 n == 129", int(sens.loc[sens.sensitivity_id.eq("S-rep-1"), "n"].iloc[0]) == 129)
    check("every sensitivity row reports an rho, a CI, resamples and a seed",
          sens[["rho", "ci_low", "ci_high", "n_resamples", "seed"]].notna().all().all())
    check("the S5 row carries the RANGE and names no chosen feature",
          "RANGE" in str(sens.loc[sens.sensitivity_id.eq("S5"), "description"].iloc[0]))
    check("a per-MMP9-feature row exists for all ten features",
          sum(1 for i in ids if re.fullmatch(r"S5-\d+", i)) == 10,
          f"{sum(1 for i in ids if re.fullmatch(r'S5-[0-9]+', i))} per-feature rows")
    check("the 102 subset is inside the 129 (a sensitivity is not a different cohort)",
          int(sens.loc[sens.sensitivity_id.eq("S-CLIN"), "n"].iloc[0])
          <= int(p["n"]))

    # ---- files that must exist ---------------------------------------------
    for fname in ("run_log.txt", "sessionInfo.txt", "results_summary.md",
                  "failed_or_excluded.md", "seed_registry.json",
                  "verification_controls.csv", "manuscript_verification.csv"):
        check(f"{fname} exists and is non-empty",
              (OUT / fname).exists() and (OUT / fname).stat().st_size > 0)

    print()
    print("=" * 74)
    print("B. THE STATISTICAL-UNIT SENTENCE AND THE AD-3 BOUND")
    print("=" * 74)
    check("results_summary.md contains the unit sentence VERBATIM",
          UNIT_SENTENCE in summary)
    check("run_log.txt contains the unit sentence VERBATIM", UNIT_SENTENCE in runlog)
    check("the unit sentence appears as a single unbroken line in the summary",
          any(line.strip() == "> **" + UNIT_SENTENCE + "**"
              or line.strip() == UNIT_SENTENCE
              or UNIT_SENTENCE in line for line in summary.splitlines()))
    for doc_name, doc in (("results_summary.md", summary),
                          ("failed_or_excluded.md", failed),
                          ("run_log.txt", runlog)):
        # An exact patient count would read "129 patients" with no bound marker next to
        # it.  The marker must be ADJACENT (within a few characters), because a long
        # labelled field can otherwise carry "<=102 (UPPER BOUND, AD-3)" on a wrapped
        # continuation line tens of characters away from the digits, and a distant-marker
        # heuristic would wave that through.  Whitespace is collapsed first so that a
        # count straddling a line break is read as one token rather than as its last
        # digit -- that mistake produced a false positive earlier.
        flat = re.sub(r"\s+", " ", doc)
        bad = []
        for m in re.finditer(r"(\d+)\s+patients", flat):
            prefix = flat[max(0, m.start() - 12):m.start()]
            if not re.search(r"(<=|≤)\s*$", prefix):
                bad.append(flat[max(0, m.start() - 45):m.end() + 10])
        check(f"{doc_name}: every '<N> patients' is immediately preceded by a bound "
              f"marker", not bad, str(bad[:3]))

    print()
    print("=" * 74)
    print("C. THE SUMMARY'S NUMBERS MATCH THE COMPUTED OUTPUTS")
    print("=" * 74)
    rho_prim, rho_clin = float(p["rho"]), float(
        sens.loc[sens.sensitivity_id.eq("S-CLIN"), "rho"].iloc[0])

    # Compare NUMBERS, not substrings.  A 4-decimal rendering in the prose is the
    # right precision for a summary; the gate's job is to confirm that the number the
    # prose states really is the computed one, so it rounds the computed value to the
    # prose's precision and looks for that token.  Matching raw substrings instead
    # produced a false failure when the prose wrote the full-precision value in
    # parentheses next to the rounded one.
    def states(doc: str, value: float, decimals: int = 4) -> bool:
        token = f"{value:+.{decimals}f}".replace("-", "\u2212")       # U+2212 minus
        ascii_token = f"{value:+.{decimals}f}"
        return (token in doc) or (ascii_token in doc)

    check(f"summary quotes the primary rho {rho_prim:+.4f}", states(summary, rho_prim))
    check(f"summary quotes the 102 rho {rho_clin:+.4f}", states(summary, rho_clin))
    check("summary quotes the primary CI lower limit",
          states(summary, float(p["ci_low"])), f"{float(p['ci_low']):+.4f}")
    check("summary quotes the primary CI upper limit",
          states(summary, float(p["ci_high"])), f"{float(p['ci_high']):+.4f}")
    check("summary states the seed", "15072026" in summary)
    check("summary states the resample count", "10 000" in summary or "10000" in summary)
    check("summary names the CI method",
          "percentile bootstrap" in summary.lower())

    # the sensitivity table in the summary vs the CSV
    for sid in ("S2", "S3", "S4", "S-rep-1"):
        r = float(sens.loc[sens.sensitivity_id.eq(sid), "rho"].iloc[0])
        check(f"summary quotes {sid} rho {r:+.4f}", states(summary, r), f"{r:+.4f}")

    check(f"summary quotes the number of estimates ({len(est)})",
          f"**{len(est)}**" in summary, f"{len(est)} estimates in estimate_table_full.csv")
    check("summary's 'every CI contains zero' claim is verified",
          ((est.ci_low > 0) | (est.ci_high < 0)).sum() == 0,
          f"{int(((est.ci_low > 0) | (est.ci_high < 0)).sum())} intervals exclude 0")
    check("summary's rho range matches the estimate table",
          states(summary, float(est.rho.min())) and states(summary, float(est.rho.max())),
          f"range [{est.rho.min():+.4f}, {est.rho.max():+.4f}]")
    check("summary's 'smallest permutation p' claim matches",
          f"{est.permutation_p.min():.4f}" in summary,
          f"min permutation p = {est.permutation_p.min():.4f}")
    check("summary states no no-log-transform alternative exists",
          "no-log-transform" in summary.lower())

    print()
    print("=" * 74)
    print("D. VERIFICATION ACTUALLY RAN")
    print("=" * 74)
    check("all manuscript-module checks passed",
          (manv["verdict"] == "PASS").all(),
          f"{len(manv)} checks, {(manv['verdict'] != 'PASS').sum()} failed")
    check("all pipeline controls passed",
          (ctrl["verdict"] == "PASS").all(),
          f"{len(ctrl)} controls, {(ctrl['verdict'] != 'PASS').sum()} failed")
    check("at least one control is capable of failing",
          (ctrl["capable_of_failing"] == "yes").sum() >= 1,
          f"{int((ctrl['capable_of_failing'] == 'yes').sum())} capable of failing")
    check("the manuscript module was the archive the frozen doc names",
          "Additional_file_2_Analysis_archive" in str(manv.iloc[0].get("note", ""))
          or "analysis_F_cross_scale" in (OUT / "verification_and_controls.md")
          .read_text(encoding="utf-8"))
    check("summary reports the manuscript-module check tally",
          f"{len(manv)}/{len(manv)}" in summary, f"{len(manv)}/{len(manv)}")
    check("summary reports the control tally",
          f"{len(ctrl)}/{len(ctrl)}" in summary, f"{len(ctrl)}/{len(ctrl)}")

    print()
    print("=" * 74)
    print("E. SCOPE AND PROHIBITIONS")
    print("=" * 74)
    # no pooling language that could be read as pooling GSE215154
    for term in ("pooled with GSE205154", "meta-analys", "random-effects"):
        hits = [d for d, doc in (("summary", summary), ("failed", failed))
                if re.search(term, doc, re.I)
                and "not" not in doc[max(0, doc.lower().find(term.lower()) - 120):
                                     doc.lower().find(term.lower())].lower()]
        check(f"no unpooled-context use of {term!r} in the prose", not hits, str(hits))
    check("summary states GSE205154 is not touched",
          "GSE205154" in summary and "not" in summary.lower())
    check("summary forbids cross-cohort absolute comparison (AD-2)",
          "absolute" in summary.lower() and "AD-2" in summary)
    check("no manuscript prose marker in the summary (no Results/Discussion headings)",
          not re.search(r"^#+\s*(results|discussion|introduction|conclusion)\b",
                        summary, re.I | re.M))
    for doc_name, doc in (("results_summary.md", summary),
                          ("failed_or_excluded.md", failed)):
        check(f"{doc_name} stays out of the manuscript directories",
              "FINAL_SUBMISSION" not in doc or "npjPrecOncol" in doc,
              "references only point at the archive it read")

    print()
    print("=" * 74)
    print(f"RESULT: {passed} passed, {len(problems)} failed")
    if problems:
        for x in problems:
            print(f"  FAILED: {x}")
    print("=" * 74)
    return 0 if not problems else 1


if __name__ == "__main__":
    sys.exit(main())
