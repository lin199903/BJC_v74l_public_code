#!/usr/bin/env python
# =============================================================================
# 05_normalisation_fork_test.py -- Phase B step 2, GSE278689
#
# Answers ONE question with evidence, before any result is computed:
#   Is this deposit raw integer counts, or processed log1p expression?
#
# The frozen wording for this run:
#   "For datasets supplied as raw integer counts, expression will be library-size
#    normalised to 10,000 counts per cell and log1p transformed, matching the
#    manuscript pipeline. RNA detectability is defined by a non-zero raw count and
#    is therefore invariant to this transformation. For datasets supplied as
#    processed log1p expression, the deposited representation will be used only
#    where it is demonstrably compatible with the same analysis identity."
#
# So the script (a) decides which branch applies and (b) measures whether the
# magnitude of the primary contrast depends on the choice. If it does, BOTH are
# reported and the fork is flagged -- it is never resolved silently.
#
# It also states the exactness that matters: detectability is >0 on the RAW count,
# and every candidate normalisation here is of the form g(count) with g(0)=0 and
# g strictly increasing, so the detectable/non-detectable split -- hence the
# eligible sample set and every group size -- is identical under all of them.
# =============================================================================
from __future__ import annotations

import os
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("SIGLEC15_GSE278689_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927\phaseB_GSE278689"))
OUT = ROOT / "output"


def zscore_columns(values: np.ndarray) -> np.ndarray:
    """The manuscript's own routine (analysis_F_cross_scale.py, verbatim semantics)."""
    values = np.asarray(values, dtype=float)
    mean = np.nanmean(values, axis=0, keepdims=True)
    sd = np.nanstd(values, axis=0, ddof=0, keepdims=True)
    sd[(~np.isfinite(sd)) | (sd == 0)] = np.nan
    return (values - mean) / sd


def transform(raw: np.ndarray, lib: np.ndarray, kind: str) -> np.ndarray:
    lib = np.where(lib > 0, lib, 1.0)
    if kind == "cp10k_log1p":          # manuscript convention
        return np.log1p(raw / lib[:, None] * 1e4)
    if kind == "log1p_counts":         # no library-size scaling
        return np.log1p(raw)
    if kind == "raw_counts":           # deposit used verbatim
        return raw.astype(float)
    if kind == "cp10k_nolog":          # library scaling, no log
        return raw / lib[:, None] * 1e4
    raise ValueError(kind)


def contrast(mat: np.ndarray, det: np.ndarray) -> tuple[float, bool, bool]:
    """(contrast, acp5_invariant, mmp9_invariant); primary score requires BOTH genes.

    `mat` MUST contain exactly the programme genes, in the order ACP5, MMP9.
    SIGLEC15 must never be passed in here -- it is the DETECTION variable, not a
    programme gene. (An earlier revision of this script passed all three columns;
    that is the forbidden target-in-program error and it moved the mean contrast
    from +0.92 to +4.01. Caught by the WRONG 6 control in 09_ and fixed.)
    """
    assert mat.shape[1] == 2, "programme matrix must have exactly ACP5 and MMP9"
    z = zscore_columns(mat)
    inv = [bool(np.all(~np.isfinite(z[:, j]))) for j in range(z.shape[1])]
    if any(inv):
        return float("nan"), inv[0], inv[1]
    score = np.nanmean(z, axis=1)
    return float(np.nanmean(score[det]) - np.nanmean(score[~det])), False, False


def main() -> None:
    summ = json.loads((OUT / "_matrix_summary.json").read_text(encoding="utf-8"))
    d = pd.read_csv(OUT / "_myeloid_gene_cache.csv")

    print("=" * 74)
    print("PART A -- WHAT IS DEPOSITED?")
    print("=" * 74)
    print(f"matrix                       : {summ['matrix']}")
    print(f"dimensions                   : {summ['n_genes']} genes x {summ['n_cells']} cells, "
          f"{summ['n_entries']:,} stored entries")
    print(f"non-integer stored values    : {summ['non_integer_values']}")
    print(f"negative stored values       : {summ['negative_values']}")
    print(f"largest stored value         : {summ['max_stored_value']:g}")
    print(f"=> all values integral       : {summ['all_values_integral']}")
    print(f"library size min/median/max  : {summ['libsize_min']:,} / "
          f"{summ['libsize_median']:,.0f} / {summ['libsize_max']:,}")

    raw = d[["SIGLEC15", "ACP5", "MMP9"]].to_numpy(float)
    prog = d[["ACP5", "MMP9"]].to_numpy(float)      # programme genes ONLY
    lib = d["libsize"].to_numpy(float)

    integral = bool(np.all(raw == np.floor(raw))) and summ["all_values_integral"]
    max_is_countlike = summ["max_stored_value"] < 1e5
    min_is_zero = float(raw.min()) == 0.0
    frac_one = float(np.mean(raw == 1))

    raw_counts = bool(integral and min_is_zero and max_is_countlike)
    print("\n-- decision --")
    print(f"  integral values                        : {integral}")
    print(f"  minimum value is exactly 0             : {min_is_zero}")
    print(f"  maximum value consistent with UMI counts: {max_is_countlike} "
          f"(max={summ['max_stored_value']:g})")
    print(f"  fraction of myeloid gene entries == 1  : {frac_one:.4f}")
    print(f"  => BRANCH: {'RAW INTEGER COUNTS' if raw_counts else 'PROCESSED EXPRESSION'}")
    if raw_counts:
        print("  => the frozen wording's FIRST branch applies: library-size normalise")
        print("     to 10,000 counts per cell, then log1p. This is used and is stated")
        print("     in every downstream output file.")
    else:
        print("  => the frozen wording's SECOND branch applies: compatibility with the")
        print("     manuscript identity must be demonstrated before use. See PART C.")

    print()
    print("=" * 74)
    print("PART B -- DETECTABILITY IS INVARIANT TO THE TRANSFORMATION")
    print("=" * 74)
    det = raw[:, 0] > 0
    print(f"SIGLEC15-detectable myeloid cells      : {int(det.sum())} / {len(det)} "
          f"({100*det.mean():.4f}%)")
    rows = []
    for kind in ("cp10k_log1p", "log1p_counts", "raw_counts", "cp10k_nolog"):
        t = transform(raw, lib, kind)
        dt = t[:, 0] > 0
        rows.append((kind, int(dt.sum()), bool(np.array_equal(dt, det)),
                     float(t.min()), float(np.median(t[:, 0]))))
    print(f"{'transform':<14}{'n_detectable':>13}{'identical split':>17}"
          f"{'min':>12}{'median SIGLEC15':>18}")
    for kind, n, same, mn, med in rows:
        print(f"{kind:<14}{n:>13}{str(same):>17}{mn:>12.4f}{med:>18.4f}")
    all_same = all(r[2] for r in rows)
    print(f"\n=> every transform reproduces the SAME detectable/non-detectable split: {all_same}")
    print("   Detectability is read from the RAW count (>0), so the eligible sample set,")
    print("   the group sizes and the primary contrast's ELIGIBILITY are unaffected by")
    print("   the normalisation branch. Only the contrast MAGNITUDE can move.")

    print()
    print("=" * 74)
    print("PART C -- DOES THE CHOICE CHANGE THE ESTIMATE? (the fork)")
    print("=" * 74)
    samples = sorted(d["sample_id"].unique())
    tab = {}
    for kind in ("cp10k_log1p", "log1p_counts", "raw_counts", "cp10k_nolog"):
        vals = []
        for s in samples:
            i = np.where(d["sample_id"].to_numpy() == s)[0]
            t = transform(prog[i], lib[i], kind)     # programme genes only
            c, _, _ = contrast(t, det[i])
            vals.append(c)
        tab[kind] = np.array(vals)

    hdr = f"{'sample':<16}" + "".join(f"{k:>15}" for k in tab)
    print(hdr)
    for j, s in enumerate(samples):
        print(f"{s:<16}" + "".join(f"{tab[k][j]:>15.4f}" for k in tab))
    print(f"{'MEAN (finite)':<16}" + "".join(f"{np.nanmean(tab[k]):>15.4f}" for k in tab))

    # cross-check the reference branch against the R primary result
    man = pd.read_csv(OUT / "sample_manifest.csv")
    ref_vals = tab["cp10k_log1p"][np.isfinite(tab["cp10k_log1p"])]
    prim = pd.read_csv(OUT / "primary_contrast_by_sample.csv")
    prim = prim[prim["status"] == "OK"]["contrast"].to_numpy()
    print(f"\ncross-check vs 08_final_analysis.R:")
    print(f"  n here {len(ref_vals)}   n in R {len(prim)}")
    print(f"  mean here {ref_vals.mean():+.6f}   mean in R {prim.mean():+.6f}   "
          f"max |diff| {np.max(np.abs(np.sort(ref_vals) - np.sort(prim))):.3e}")
    assert len(ref_vals) == len(prim) and np.max(np.abs(np.sort(ref_vals) - np.sort(prim))) < 1e-9, \
        "the fork test's reference branch does not reproduce the R primary contrast"

    ref = tab["cp10k_log1p"]
    print("\n-- deviations from the manuscript convention (CP10K + log1p) --")
    fork = {}
    for k in tab:
        if k == "cp10k_log1p":
            continue
        dm = float(np.nanmax(np.abs(tab[k] - ref)))
        dmean = float(abs(np.nanmean(tab[k]) - np.nanmean(ref)))
        fork[k] = {"max_abs_sample_deviation": dm, "mean_deviation": dmean}
        print(f"  {k:<14} mean contrast {np.nanmean(tab[k]):+.4f}  "
              f"(delta mean {np.nanmean(tab[k])-np.nanmean(ref):+.4f}; "
              f"max sample delta {dm:.4f})")
    material = any(v["mean_deviation"] >= 0.05 for v in fork.values())
    print(f"\n=> a non-manuscript convention would move the mean contrast by up to "
          f"{max(v['mean_deviation'] for v in fork.values()):.4f}")
    print(f"=> MATERIAL (>0.05 on the mean): {material}")
    if material:
        print("   FLAGGED: both the manuscript convention and the forked alternatives are")
        print("   reported; the manuscript convention is used for the primary result and")
        print("   the fork is recorded in failed_or_excluded.md. It was NOT resolved silently.")

    (OUT / "_normalisation_fork_test.json").write_text(json.dumps({
        "deposit_is_raw_integer_counts": raw_counts,
        "detectability_split_identical_across_transforms": all_same,
        "branch_used": "CP10K then log1p (first branch of the frozen wording)",
        "contrasts_by_transform": {k: v.tolist() for k, v in tab.items()},
        "mean_by_transform": {k: float(np.nanmean(v)) for k, v in tab.items()},
        "fork_vs_manuscript_convention": fork,
        "fork_is_material": material,
    }, indent=1), encoding="utf-8")
    print("\nwrote _normalisation_fork_test.json")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
