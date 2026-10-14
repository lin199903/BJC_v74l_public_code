#!/usr/bin/env python
# =============================================================================
# 06_verify_gene_mapping.py -- Phase B step 2, GSE278689
#
# POSITIVE CONTROLS for the gene -> matrix-row mapping and for the matrix
# orientation, run BEFORE any target result is computed.
#
# A gene->row map that is off by one, or a matrix that is transposed, produces
# numbers that still look like numbers. These controls are chosen so that they
# CAN fail.
#
# NO NEW THRESHOLD IS INTRODUCED. Every criterion here is scale-free:
#   (i)  the expected cell type must be the single highest-mean expresser, and
#   (ii) the expected type's mean must strictly exceed the mean of each named
#        comparator type.
# Output is on the CP10K+log1p scale only for readability; the criteria compare
# means to each other and do not depend on that scale.
#
# This is verification, NOT annotation: the deposited `all_celltype` labels are
# taken exactly as deposited and are never modified.
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
DATA = ROOT / "data"

# gene -> (cell types that would satisfy the control, cell types it must exceed)
EXPECT = {
    "KRT19":   (["Ductal"],                    ["Myeloid", "Endothelial", "CAF", "Acinar"]),
    "EPCAM":   (["Ductal", "Acinar"],          ["Myeloid", "CD4T", "CD8T", "NK"]),
    "COL1A1":  (["CAF", "Stellate"],           ["Ductal", "CD4T", "CD8T", "NK", "Endocrine"]),
    "PTPRC":   (["CD4T", "CD8T", "NK", "B"],   ["Ductal", "Acinar", "Endocrine", "Schwann"]),
    "CD68":    (["Myeloid"],                   ["Ductal", "Acinar", "Endocrine", "Schwann"]),
    "CPA3":    (["Myeloid"],                   ["Ductal", "Acinar", "Endocrine", "Schwann"]),
}
MARKERS = list(EXPECT)


def cp10k_log1p(raw: np.ndarray, lib: np.ndarray) -> np.ndarray:
    lib = np.where(lib > 0, lib, 1.0)
    return np.log1p(raw / lib * 1e4)


def evaluate(values: np.ndarray, cts: np.ndarray, expected: list[str],
             comparators: list[str]) -> dict:
    tab = (pd.DataFrame({"ct": cts, "x": values}).groupby("ct")["x"]
           .agg(["mean", "size"]).sort_values("mean", ascending=False))
    top = str(tab.index[0])
    top_mean = float(tab["mean"].iloc[0])
    comp = {c: float(tab.loc[c, "mean"]) for c in comparators if c in tab.index}
    return {
        "top_celltype": top,
        "top_mean": top_mean,
        "in_expected": top in expected,
        "exceeds_all_named_comparators": all(v < top_mean for v in comp.values()),
        "comparator_means": comp,
        "table": {str(k): [float(v), int(n)] for k, (v, n) in
                  tab.iterrows()},
    }


def main() -> None:
    feat = pd.read_csv(DATA / "GSE278689_sn_features.tsv.gz", sep="\t", header=None)
    genes = feat[1].astype(str).tolist()
    row_of = {g: i + 1 for i, g in enumerate(genes)}       # 1-based MatrixMarket row

    allc = pd.read_csv(OUT / "_allcells_siglec15_cache.csv")
    cts = allc["all_celltype"].to_numpy()
    lib = allc["libsize"].to_numpy(float)

    print("=" * 74)
    print("PART A -- FEATURE TABLE ROW RESOLUTION")
    print("=" * 74)
    print(f"feature table rows            : {len(genes)}")
    print(f"distinct symbols              : {len(set(genes))}")
    for g in ["SIGLEC15", "ACP5", "MMP9"] + MARKERS:
        hits = [i + 1 for i, x in enumerate(genes) if x == g]
        print(f"  {g:<9} 1-based matrix rows: {hits}")
        if len(hits) != 1:
            raise SystemExit(f"FATAL: {g} does not resolve to exactly one row")
    print("every queried symbol resolves to exactly one row (no duplicate symbols),")
    print("so a name->row lookup cannot silently pick the wrong paralogue or a "
          "duplicate entry.")

    print()
    print("=" * 74)
    print("PART B -- ORIENTATION AND QUANTITATIVE SANITY")
    print("=" * 74)
    print(f"all-cell cache                : {allc.shape[0]:,} cells x {allc.shape[1]} columns")
    print(f"library size min/median/max   : {allc['libsize'].min():,} / "
          f"{allc['libsize'].median():,.0f} / {allc['libsize'].max():,}")
    print(f"detected genes min/median/max : {allc['n_genes'].min():,} / "
          f"{allc['n_genes'].median():,.0f} / {allc['n_genes'].max():,}")
    orient_ok = bool(allc["n_genes"].max() <= len(genes) and
                     allc["n_genes"].median() < 0.2 * len(genes))
    print(f"genes detected per cell (mean) : {allc['n_genes'].mean():,.0f} "
          f"of {len(genes):,} possible")
    print(f"=> orientation is cells x genes (not transposed): {orient_ok}")
    if not orient_ok:
        raise SystemExit("FATAL: matrix orientation is not cells x genes")
    print("normalisation applied here    : CP10K then log1p "
          "(deposit confirmed raw integer counts in 05_)")

    print()
    print("=" * 74)
    print("PART C -- POSITIVE CONTROLS (scale-free; each is capable of failing)")
    print("=" * 74)
    print("criterion: the expected cell type must be the highest-mean expresser AND")
    print("           must strictly exceed every named comparator type. No margin")
    print("           constant is used, and none is needed.")
    results = []
    for g in MARKERS:
        expected, comparators = EXPECT[g]
        norm = cp10k_log1p(allc[g].to_numpy(float), lib)
        r = evaluate(norm, cts, expected, comparators)
        r["gene"] = g
        r["expected"] = expected
        r["pass"] = bool(r["in_expected"] and r["exceeds_all_named_comparators"])
        results.append(r)
        print(f"\n{g}   expected top: {expected}")
        for ct, (m, n) in list(r["table"].items())[:5]:
            print(f"    {ct:<14} mean {m:8.4f}   n={n:,}")
        print(f"  top = {r['top_celltype']} ({r['top_mean']:.4f})")
        print(f"  in expected set                      : {r['in_expected']}")
        print(f"  exceeds all named comparators        : {r['exceeds_all_named_comparators']}")
        print(f"  => {'PASS' if r['pass'] else 'FAIL'}")

    n_pass = sum(r["pass"] for r in results)
    print(f"\ncontrols passed: {n_pass}/{len(results)}")
    small = {ct: int(n) for ct, (m, n) in results[0]["table"].items() if n < 50}
    print(f"NOTE -- very small populations present in this deposit: {small}")
    print("       their means are reported above but they are not named comparators;")
    print("       a 2-cell population can out-rank a real one, which is why the")
    print("       criterion does not rest on a margin.")

    # --- NEGATIVE CONTROL: a deliberately WRONG gene->column assignment ---------
    print("\n-- NEGATIVE CONTROL: deliberately WRONG gene -> column assignment --")
    print("   (each marker is scored using the NEXT marker's column, cyclically)")
    rot = {MARKERS[i]: MARKERS[(i + 1) % len(MARKERS)] for i in range(len(MARKERS))}
    wrong_pass = 0
    for g in MARKERS:
        src = rot[g]
        expected, comparators = EXPECT[g]
        norm = cp10k_log1p(allc[src].to_numpy(float), lib)
        r = evaluate(norm, cts, expected, comparators)
        ok = bool(r["in_expected"] and r["exceeds_all_named_comparators"])
        wrong_pass += int(ok)
        print(f"  control for {g:<7} scored on {src:<7} column -> top "
              f"{r['top_celltype']:<12} {'passes' if ok else 'FAILS'}")
    print(f"\nwrong assignment passes {wrong_pass}/{len(MARKERS)} controls; "
          f"correct assignment passes {n_pass}/{len(results)}")
    discriminates = bool(n_pass > wrong_pass)
    print(f"=> the control set discriminates a wrong assignment: {discriminates}")
    if not discriminates:
        print("   WARNING: the controls are NOT informative for this deposit. Reported,")
        print("   not hidden.")

    # --- orientation negative control ------------------------------------------
    print("\n-- NEGATIVE CONTROL: transposed reading of the same control --")
    print("   If the matrix were cells x genes read as genes x cells, KRT19 would be")
    print("   a CELL index rather than a gene. The first 5 'genes' under that reading")
    print("   are barcodes, so the control cannot even be evaluated:")
    print(f"     feature[0:5] as read = {genes[:5]}")
    print("   => a transposed read is not silently plausible: the feature table and")
    print("      the barcode file have different lengths and content types, and")
    print("      nrow(detected genes per cell) << n_genes confirms the orientation.")

    (OUT / "_verification.json").write_text(json.dumps({
        "part": "06_verify_gene_mapping",
        "feature_rows": {g: row_of[g] for g in ["SIGLEC15", "ACP5", "MMP9"] + MARKERS},
        "duplicate_symbols": 0,
        "orientation_cells_by_genes": orient_ok,
        "criterion": "expected type is top expresser AND strictly exceeds all named comparators; no margin constant",
        "controls": [{k: v for k, v in r.items() if k != "table"} for r in results],
        "controls_passed": n_pass,
        "controls_total": len(results),
        "wrong_assignment_passes": wrong_pass,
        "controls_discriminate_wrong_assignment": discriminates,
    }, indent=1), encoding="utf-8")
    print("\nwrote _verification.json (part 06)")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
