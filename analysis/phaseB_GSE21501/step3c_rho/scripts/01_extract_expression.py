#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE21501 -- Phase B step 3c-2, stage 1 of 3: EXTRACT.

Purpose
-------
Turn the frozen deposit into a per-array table of the three frozen target genes,
and nothing else.  This stage computes no correlation and inspects no relationship
between genes: it only produces
    work/GSE21501_target_expression.csv   (132 rows x 4 cols + provenance cols)
and a QC record for it.

Why this stage exists separately
--------------------------------
Until this step the run was **value-blind** (frozen_cohort_definition.md S0): the
manifest was built from metadata and platform annotation only, and no per-sample
expression value of SIGLEC15, ACP5 or MMP9 had ever been loaded.  This script is the
first place in the GSE21501 run where an expression value is read.  Keeping it as its
own stage means the moment the blindness ends is a single, auditable file.

Gene-identity rule (frozen; transcribed, not re-chosen)
-------------------------------------------------------
frozen_cohort_definition.md S7 fixes the probe->gene rule as the **manuscript's own**
identity, i.e. `load_geo_expression()` of the archived
`analysis_F_cross_scale.py`:

    symbol = annotation_symbol.split("///")[0]
    selected.groupby("gene").mean(numeric_only=True)

That is: take the annotation row's symbol cell, keep the text before the first "///",
select every platform feature whose resulting symbol is a target gene, and collapse
each gene to **one** value per array by **arithmetic mean** over its features.

GPL4133 is the Feature Number version of Agilent-014850, so the feature identifier is
a bare integer and the series-matrix `ID_REF` column carries those same integers.  The
target features are **re-derived here from the platform annotation** rather than
hard-coded from the frozen document, so that the frozen values can be checked against
an independent derivation instead of being assumed.

Inputs (read-only; hashes recorded by 03_report.py)
    ../_probe/GPL4133.annot.gz                  platform annotation
    ../_probe/GSE21501_series_matrix.txt.gz     the deposited matrix
    ../sample_manifest.csv                      frozen cohort manifest
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
STEP = HERE.parent
PHASE = STEP.parent
WORK = STEP / "work"
WORK.mkdir(parents=True, exist_ok=True)

ANNOT = PHASE / "_probe" / "GPL4133.annot.gz"
MATRIX = PHASE / "_probe" / "GSE21501_series_matrix.txt.gz"
MANIFEST = PHASE / "sample_manifest.csv"

# ---------------------------------------------------------------------------
# FROZEN constants -- copied from frozen_cohort_definition.md, never chosen here
# ---------------------------------------------------------------------------
FROZEN_SHA256_MANIFEST = "3e11938ae5c7e682c9e53aad93d051579764a4ec8173d69d571dbeda0f3ac4e4"

# S7: the frozen probe mapping this script must *reproduce*, not assume.
FROZEN_FEATURES = {
    "SIGLEC15": ["31474"],
    "ACP5": ["10450"],
    "MMP9": ["1570", "1678", "14430", "14852", "21363",
             "25413", "28592", "28675", "40398", "44612"],
}
TARGET_GENES = ["SIGLEC15", "ACP5", "MMP9"]

# S9: the frozen programme definition.
PROGRAMS = {"Derosiers_ACP5_MMP9_core": ["ACP5", "MMP9"]}
PRIMARY_PROGRAM = "Derosiers_ACP5_MMP9_core"
TARGET = "SIGLEC15"


def sha256_file(path: Path, chunk: int = 8 * 1024 * 1024) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            d.update(b)
    return d.hexdigest()


# ---------------------------------------------------------------------------
# Step 1. Platform annotation -> symbol map, by the manuscript's own column rule
# ---------------------------------------------------------------------------
def parse_soft_platform_annot(path: Path) -> tuple[dict[str, str], dict]:
    """Manuscript rule: first header column whose name contains 'symbol'.

    Three counts are kept apart on purpose, because collapsing them is how a
    row-accounting error hides:

      n_raw_data_rows           every row of the platform data table, INCLUDING rows
                                that carry no usable symbol
      n_rows_without_usable_symbol
                                empty symbol cell, '---', or a row with too few fields
                                to reach the symbol column
      n_data_rows               rows that carry a gene symbol and therefore can reach
                                the gene-identity rule

    The matrix carries one column-value per row of `n_data_rows` plus the
    symbol-less rows, so `n_raw_data_rows - n_data_rows` must equal the count of
    matrix IDs that the annotated table cannot name.  That identity is asserted by
    stage 2's control C12 rather than assumed here.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    header = None
    rows = {}
    n_data_rows = 0
    n_raw_data_rows = 0
    n_empty_symbol_rows = 0
    n_dash_symbol_rows = 0
    n_short_rows = 0
    n_header_rows = 0
    symbol_idx = None
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            if header is None and parts[0] == "ID":
                header = parts
                n_header_rows += 1
                hits = [i for i, h in enumerate(header) if "symbol" in h.lower()]
                if not hits:
                    raise SystemExit("no header column containing 'symbol'")
                symbol_idx = hits[0]
                continue
            # GEO SOFT markers, not data: '^Annotation', '!Annotation_*',
            # '#ID = ...', and the '!platform_table_end' terminator.
            if parts[0].startswith(("#", "!", "^")):
                n_header_rows += 1
                continue
            if header is None:
                continue
            n_raw_data_rows += 1
            if len(parts) <= symbol_idx:
                n_short_rows += 1
                continue
            raw_cell = parts[symbol_idx].strip()
            if raw_cell == "":
                n_empty_symbol_rows += 1
                continue
            if raw_cell == "---":
                n_dash_symbol_rows += 1
                continue
            n_data_rows += 1
            symbol = raw_cell.split("///")[0].strip()
            if symbol:
                rows[parts[0].strip()] = symbol
    n_without_usable = n_empty_symbol_rows + n_dash_symbol_rows + n_short_rows
    info = {
        "header_row": header,
        "symbol_column_index": symbol_idx,
        "symbol_column_name": header[symbol_idx] if header else None,
        "n_header_rows_skipped": n_header_rows,
        "n_raw_data_rows": n_raw_data_rows,
        "n_rows_empty_symbol_cell": n_empty_symbol_rows,
        "n_rows_dash_symbol_cell": n_dash_symbol_rows,
        "n_rows_too_short_to_reach_symbol_column": n_short_rows,
        "n_rows_without_usable_symbol": n_without_usable,
        "n_data_rows": n_data_rows,
        "row_accounting_closes": bool(n_raw_data_rows == n_data_rows + n_without_usable),
        "n_features_with_symbol": len(rows),
        "n_distinct_symbols": len(set(rows.values())),
        "n_multi_symbol_cells": None,  # filled below
    }
    return rows, info


def audit_multi_symbol(path: Path, symbol_idx: int) -> dict:
    """How many annotation rows carry the '///' delimiter, and are any targets
    demoted out of first position by the split rule?"""
    opener = gzip.open if path.suffix == ".gz" else open
    started = False
    n_multi = 0
    target_demoted = []
    targets = set(TARGET_GENES)
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not started:
                if parts[0] == "ID":
                    started = True
                continue
            if len(parts) <= symbol_idx:
                continue
            cell = parts[symbol_idx]
            if "///" not in cell:
                continue
            n_multi += 1
            first = cell.split("///")[0].strip()
            rest = [s.strip() for s in cell.split("///")[1:]]
            for t in targets:
                if t in rest and t != first:
                    target_demoted.append({"feature": parts[0], "cell": cell})
    return {"n_rows_with_slashdelim": n_multi,
            "target_genes_demoted_by_split": target_demoted}


# ---------------------------------------------------------------------------
# Step 2. Series matrix -> expression frame
# ---------------------------------------------------------------------------
def read_series_matrix(path: Path) -> tuple[pd.DataFrame, list[str]]:
    """Return (frame indexed by ID_REF with one column per GSM, list of GSMs)."""
    opener = gzip.open if path.suffix == ".gz" else open
    in_table = False
    header = None
    records = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not in_table:
                if line.startswith("!series_matrix_table_begin"):
                    in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            parts = line.split("\t")
            if header is None:
                header = [p.strip('"') for p in parts]
                continue
            records.append(parts)
    index = [r[0].strip('"') for r in records]
    gsms = header[1:]
    data = {gsm: [r[i + 1].strip('"') for r in records] for i, gsm in enumerate(gsms)}
    frame = pd.DataFrame(data, index=index)
    frame.index.name = "ID_REF"
    return frame, gsms


def main() -> int:
    qc: dict = {"stage": "01_extract_expression"}

    # -- hash gate ---------------------------------------------------------
    manifest_hash = sha256_file(MANIFEST)
    qc["manifest_sha256_observed"] = manifest_hash
    qc["manifest_sha256_frozen"] = FROZEN_SHA256_MANIFEST
    qc["manifest_sha256_matches"] = manifest_hash == FROZEN_SHA256_MANIFEST
    if not qc["manifest_sha256_matches"]:
        raise SystemExit(
            "FATAL: sample_manifest.csv hash does not match the frozen value; "
            "the frozen cohort definition does not cover this file."
        )
    print(f"[ok] manifest sha256 matches frozen value: {manifest_hash}")

    qc["inputs"] = {}
    for p in (ANNOT, MATRIX, MANIFEST):
        qc["inputs"][str(p.relative_to(PHASE))] = {
            "sha256": sha256_file(p), "bytes": p.stat().st_size}

    # -- platform annotation ----------------------------------------------
    symbols, an_info = parse_soft_platform_annot(ANNOT)
    multi = audit_multi_symbol(ANNOT, an_info["symbol_column_index"])
    an_info.update(multi)
    qc["platform_annotation"] = an_info
    print(f"[ok] platform: {an_info['n_raw_data_rows']} data rows "
          f"({an_info['n_rows_empty_symbol_cell']} empty symbol cell, "
          f"{an_info['n_rows_dash_symbol_cell']} '---', "
          f"{an_info['n_rows_too_short_to_reach_symbol_column']} too short); "
          f"{an_info['n_data_rows']} rows carry a gene symbol; "
          f"row accounting closes = {an_info['row_accounting_closes']}; "
          f"symbol col {an_info['symbol_column_index']!r} "
          f"({an_info['symbol_column_name']!r}), "
          f"{an_info['n_distinct_symbols']} distinct symbols")
    if not an_info["row_accounting_closes"]:
        raise SystemExit("FATAL: platform row accounting does not close.")

    # -- re-derive the probe mapping --------------------------------------
    derived = {}
    for gene in TARGET_GENES:
        feats = [f for f, s in symbols.items() if s == gene]
        derived[gene] = sorted(feats, key=lambda x: int(x))
    qc["probe_mapping_derived"] = derived
    qc["probe_mapping_frozen"] = FROZEN_FEATURES
    qc["probe_mapping_matches_frozen"] = {
        g: derived[g] == FROZEN_FEATURES[g] for g in TARGET_GENES}
    for g in TARGET_GENES:
        ok = qc["probe_mapping_matches_frozen"][g]
        print(f"[{'ok' if ok else 'FAIL'}] {g}: derived {len(derived[g])} feature(s), "
              f"frozen {len(FROZEN_FEATURES[g])} feature(s), match={ok}")
    if not all(qc["probe_mapping_matches_frozen"].values()):
        raise SystemExit("FATAL: re-derived probe mapping disagrees with the frozen "
                         "mapping; refusing to continue.")

    # -- series matrix ------------------------------------------------------
    frame, gsms = read_series_matrix(MATRIX)
    unnamed = sorted(set(frame.index) - set(symbols), key=lambda x: int(x))
    qc["series_matrix"] = {
        "n_features": int(frame.shape[0]),
        "n_arrays": int(frame.shape[1]),
        "first_gsm": gsms[0], "last_gsm": gsms[-1],
        "n_unique_gsm": len(set(gsms)),
        "id_ref_is_integer_strings": bool(
            all(re.fullmatch(r"\d+", i) for i in frame.index)),
        "platform_table_ids_missing_from_matrix": sorted(
            set(symbols) - set(frame.index))[:20],
        "n_platform_table_ids_missing_from_matrix": len(set(symbols) - set(frame.index)),
        "matrix_ids_missing_from_platform": unnamed[:40],
        "n_matrix_ids_missing_from_platform": len(unnamed),
    }
    print(f"[ok] matrix: {frame.shape[0]} features x {frame.shape[1]} arrays; "
          f"ID_REF all bare integers = {qc['series_matrix']['id_ref_is_integer_strings']}")
    print(f"[ok] matrix IDs the annotated table cannot name: {len(unnamed)} "
          f"(platform rows whose symbol cell is not usable)")

    # -- select target features (manuscript identity rule) -----------------
    gene_of = pd.Series(frame.index.map(symbols), index=frame.index)
    keep = gene_of.isin(TARGET_GENES).to_numpy()
    qc["selection"] = {
        "n_features_selected": int(keep.sum()),
        "n_rows_in_matrix_with_no_symbol": int(gene_of.isna().sum()),
        "selected_features_by_gene": {
            g: sorted([f for f, k in zip(frame.index, keep)
                       if k and symbols.get(f) == g], key=lambda x: int(x))
            for g in TARGET_GENES},
    }
    print(f"[ok] selected {int(keep.sum())} features for {TARGET_GENES}")

    numeric = frame.loc[keep].apply(pd.to_numeric, errors="coerce")
    n_na = int(numeric.isna().to_numpy().sum())
    qc["selection"]["n_cells_selected"] = int(numeric.size)
    qc["selection"]["n_non_numeric_or_missing_cells"] = n_na
    print(f"[ok] non-numeric/missing cells among the selected features: {n_na}")

    # gene-level collapse: mean over features, exactly as load_geo_expression does
    genes = pd.Series(gene_of.loc[keep].to_numpy(), index=numeric.index)
    selected = numeric.copy()
    selected.insert(0, "gene", genes.to_numpy())
    collapsed = selected.groupby("gene").mean(numeric_only=True)          # genes x arrays
    expr = collapsed.T                                                     # arrays x genes
    expr.index.name = "gsm_accession"
    expr = expr[TARGET_GENES]

    # per-gene collapse arity, so the 10-feature MMP9 mean is visible, not assumed
    qc["collapse"] = {
        "rule": "selected.groupby('gene').mean(numeric_only=True)  "
                "[load_geo_expression identity]",
        "n_features_per_gene": {g: len(qc["selection"]["selected_features_by_gene"][g])
                                for g in TARGET_GENES},
        "n_arrays_per_gene": {g: int(collapsed.loc[g].notna().sum()) for g in TARGET_GENES},
    }
    for g in TARGET_GENES:
        print(f"     {g}: {qc['collapse']['n_features_per_gene'][g]} feature(s) "
              f"-> mean over {qc['collapse']['n_arrays_per_gene'][g]} arrays")

    # -- attach frozen manifest columns (metadata only; no expression joined) --
    man = pd.read_csv(MANIFEST)
    merged = man.merge(expr.reset_index(), on="gsm_accession", how="left",
                       validate="one_to_one")
    qc["manifest_join"] = {
        "manifest_rows": int(len(man)),
        "expression_rows": int(len(expr)),
        "merged_rows": int(len(merged)),
        "gsm_missing_expression": merged.loc[merged["SIGLEC15"].isna(),
                                             "gsm_accession"].tolist(),
        "gsm_extra_in_expression": sorted(set(expr.index) - set(man["gsm_accession"])),
    }

    out_cols = ["gsm_accession", "channel2_sample_title", "biomaterial_provider",
                "tissue_type", "tumour_or_non_tumour", "primary_or_metastatic",
                "patient_id_if_available", "platform", "clinically_annotated",
                "usable_siglec15_probe", "usable_acp5_probe", "usable_mmp9_probe",
                "eligible_for_main_meta"] + TARGET_GENES
    merged[out_cols].to_csv(WORK / "GSE21501_target_expression.csv", index=False)
    print(f"[ok] wrote {WORK / 'GSE21501_target_expression.csv'} "
          f"({len(merged)} rows)")

    # -- per-array expression QC (descriptive only; no cross-gene relationship) --
    desc = {}
    for g in TARGET_GENES:
        v = merged[g].to_numpy(float)
        desc[g] = {
            "n": int(np.isfinite(v).sum()),
            "min": float(np.nanmin(v)), "max": float(np.nanmax(v)),
            "mean": float(np.nanmean(v)), "sd_ddof0": float(np.nanstd(v)),
            "median": float(np.nanmedian(v)),
            "n_distinct_values": int(len(np.unique(v[np.isfinite(v)]))),
            "n_zero_variance": int(np.nanstd(v) == 0),
        }
    qc["expression_descriptives"] = desc
    for g in TARGET_GENES:
        d = desc[g]
        print(f"     {g}: n={d['n']} range [{d['min']:.4f}, {d['max']:.4f}] "
              f"sd={d['sd_ddof0']:.4f} distinct={d['n_distinct_values']}")

    # MMP9's ten features are collapsed by mean; record their individual values too,
    # because S5 (per-MMP9-feature) needs them and re-reading the 12 MB deposit for
    # that would be wasteful.
    mmp9_feats = FROZEN_FEATURES["MMP9"]
    feat_tbl = numeric.loc[mmp9_feats].T
    feat_tbl.index.name = "gsm_accession"
    feat_tbl.columns = [f"MMP9_feature_{c}" for c in feat_tbl.columns]
    feat_tbl.reset_index().to_csv(WORK / "GSE21501_MMP9_features.csv", index=False)
    acp5_tbl = numeric.loc[FROZEN_FEATURES["ACP5"]].T
    acp5_tbl.index.name = "gsm_accession"
    acp5_tbl.columns = ["ACP5_feature_10450"]
    acp5_tbl.reset_index().to_csv(WORK / "GSE21501_ACP5_feature.csv", index=False)
    qc["side_tables"] = {
        "mmp9_features": mmp9_feats, "acp5_feature": FROZEN_FEATURES["ACP5"],
        "note": "side tables carry the SAME numeric cells as the collapsed frame; "
                "they exist so S5 does not require a second pass over the deposit",
    }
    qc["outputs"] = {
        "work/GSE21501_target_expression.csv":
            sha256_file(WORK / "GSE21501_target_expression.csv"),
        "work/GSE21501_MMP9_features.csv":
            sha256_file(WORK / "GSE21501_MMP9_features.csv"),
        "work/GSE21501_ACP5_feature.csv":
            sha256_file(WORK / "GSE21501_ACP5_feature.csv"),
    }

    (WORK / "stage01_qc.json").write_text(
        json.dumps(qc, indent=2, default=str), encoding="utf-8")
    print(f"[ok] wrote {WORK / 'stage01_qc.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
