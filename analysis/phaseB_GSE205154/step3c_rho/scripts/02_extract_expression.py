#!/usr/bin/env python
"""
02_extract_expression.py -- Step 3c-1, GSE205154 only.

Purpose
-------
Extract the three frozen genes (SIGLEC15, ACP5, MMP9) from the three deposited
gene-level matrices, and nothing else. This is the ONLY script in this step that
opens an expression matrix.

Provenance re-verification, done here rather than assumed:

  P1  SHA-256 of all three .txt.gz files equals the values recorded in
      frozen_cohort_definition.md section 1.1
  P2  each matrix has exactly 60,554 feature rows
  P3  each matrix has exactly 289 sample columns, in the SAME order across all three
  P4  each of the three genes occurs on exactly ONE row in each matrix, is
      protein_coding, carries the frozen Ensembl id, and (TMM matrix) has
      TMM_pass_filter == TRUE
  P5  the row index of each gene is identical across the three matrices
      (frozen values: SIGLEC15 = 17332, ACP5 = 2686, MMP9 = 2418, 0-based)
      -- checked independently here, not copied from the freeze
  P6  the 289 matrix sample columns equal the 289 manifest sample_titles as a set
  P7  the value scale of each matrix is recorded, so that the log2 transform applied
      downstream is demonstrably the right one (TMM/TPM are not log scaled; counts are
      large integers)

Outputs (to step3c_rho/output/):
  expression_3genes.csv        genes x samples, one column block per representation
  expression_extraction_log.txt

Run:
  python 02_extract_expression.py
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEP = HERE.parent
PHASEB = STEP.parent
DATA = PHASEB / "data"
OUT = STEP / "output"
OUT.mkdir(parents=True, exist_ok=True)

# frozen in frozen_cohort_definition.md section 1.1
MATRICES = {
    "TMM": {
        "path": DATA / "GSE205154_Gene_Level_TMM_Estimates.txt.gz",
        "sha256": "9f362bb06995e28d3f0f2297f89ddf0f1190235c76a1a8cfb5c199d8211895c3",
        "value_note": "gene-level TMM (edgeR method), not log scaled; carries TMM_pass_filter",
    },
    "TPM": {
        "path": DATA / "GSE205154_Gene_Level_TPM_Estimates.txt.gz",
        "sha256": "686582f0c00467de844d9ee93356e192da2c03c047f804f2b9d5270c4b62f152",
        "value_note": "gene-level TPM, not log scaled",
    },
    "COUNTS": {
        "path": DATA / "GSE205154_Gene_Level_Counts_Estimates.txt.gz",
        "sha256": "4c7004cfe15a94d868e24cf651f0e5d05779709bf233a74b48b084efdf238fdc",
        "value_note": "gene-level estimated counts (kallisto/tximport)",
    },
}

# frozen identifiers, INDEPENDENTLY re-derived below
GENE_ENSEMBL = {
    "SIGLEC15": "ENSG00000197046",
    "ACP5": "ENSG00000102575",
    "MMP9": "ENSG00000100985",
}
FROZEN_ROW_INDEX = {"SIGLEC15": 17332, "ACP5": 2686, "MMP9": 2418}
WANTED = ["SIGLEC15", "ACP5", "MMP9"]

EXPECT_FEATURE_ROWS = 60554
EXPECT_SAMPLES = 289

log_lines: list[str] = []
failures: list[str] = []


def log(t: str = "") -> None:
    log_lines.append(t)
    print(t)


def check(cond: bool, ok: str, bad: str) -> bool:
    log(f"  [{'PASS' if cond else 'FAIL'}] {ok if cond else bad}")
    if not cond:
        failures.append(bad)
    return cond


def sha256_file(path: Path, chunk: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def norm(s: str) -> str:
    return s.strip().lstrip("\ufeff")


# ============================================================== begin
log("=" * 78)
log("STEP 3c-1  GSE205154  --  EXPRESSION EXTRACTION (3 genes only)")
log(f"generated : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 78)
log()
log("Frame of reference for the extraction:")
log("  The manuscript's own loader (analysis_F_cross_scale.py, v39 archive,")
log("  load_bulk()) reads each cohort as a GENES x SAMPLES table and then transposes")
log("  it:  frame = genes_by_samples.T   ->  SAMPLES as rows, GENES as columns.")
log("  zscore_columns() standardises axis=0, i.e. down the ROWS, i.e. ACROSS SAMPLES.")
log("  The deposited GSE205154 matrices are already GENES x SAMPLES, so the extraction")
log("  reproduces the same orientation directly. This orientation is re-asserted")
log("  explicitly in script 03.")
log()

# ---------------------------------------------------------------- read manifest titles
manifest_titles: list[str] = []
manifest_gsm: list[str] = []
with (PHASEB / "sample_manifest.csv").open("r", encoding="utf-8", newline="") as fh:
    for row in csv.DictReader(fh):
        manifest_titles.append(row["sample_title"].strip())
        manifest_gsm.append(row["gsm_accession"].strip())
log(f"manifest: {len(manifest_titles)} sample titles loaded (titles only; no expression)")

# ---------------------------------------------------------------- extract
extracted: dict[str, dict] = {}
sample_order_ref: list[str] | None = None

for tag, spec in MATRICES.items():
    path: Path = spec["path"]
    log()
    log("=" * 78)
    log(f"MATRIX : {tag}")
    log(f"file   : {path.name}")
    log(f"note   : {spec['value_note']}")
    log("=" * 78)

    log("--- P1  SHA-256 re-verified against the frozen record ---")
    actual = sha256_file(path)
    check(actual == spec["sha256"],
          f"SHA-256 matches frozen value ({actual[:16]}...)",
          f"SHA-256 MISMATCH: computed {actual}, frozen {spec['sha256']}")

    with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as fh:
        header = norm(fh.readline()).split("\t")
        # locate annotation columns
        def colidx(name: str) -> int:
            return header.index(name)

        i_symbol = colidx("hgnc_symbol")
        i_ensembl = colidx("ensembl_gene_id")
        i_biotype = colidx("gene_biotype")
        i_filter = header.index("TMM_pass_filter") if "TMM_pass_filter" in header else None
        i_gene_name = header.index("gene_name") if "gene_name" in header else None

        # The sample columns are labelled by SAMPLE TITLE (ST-########-X), not by GSM.
        # Detect the annotation block as the leading fixed set of named fields and treat
        # everything after it as samples -- this is robust to either labelling.
        FIXED = {"ID", "ensembl_gene_id", "hgnc_symbol", "gene_biotype", "gene_name",
                 "TMM_pass_filter"}
        sample_pos = [i for i, c in enumerate(header) if c not in FIXED]
        annotation_cols = [c for c in header if c in FIXED]
        sample_cols = [norm(header[i]) for i in sample_pos]
        log()
        log(f"  annotation columns : {annotation_cols}")
        log(f"  sample columns     : {len(sample_cols)}  "
            f"(positions {sample_pos[0]}..{sample_pos[-1]})")
        log(f"  sample labels are  : SAMPLE TITLES (ST-########-X), not GSM accessions")
        log(f"  first 3 / last 3   : {sample_cols[:3]} ... {sample_cols[-3:]}")

        check(len(sample_cols) == EXPECT_SAMPLES,
              f"289 sample columns", f"{len(sample_cols)} sample columns, expected {EXPECT_SAMPLES}")

        # P6  matrix columns == manifest titles as a set
        missing = sorted(set(manifest_titles) - set(sample_cols))
        extra = sorted(set(sample_cols) - set(manifest_titles))
        check(not missing and not extra,
              "matrix sample columns equal the 289 manifest sample titles (as a set)",
              f"column/title mismatch: missing={missing[:5]} extra={extra[:5]}")
        # and the ORDER must be the deposited order; record it for the identity check
        check(len(set(sample_cols)) == EXPECT_SAMPLES,
              "289 distinct sample titles in the matrix header (no duplicate column)",
              f"{len(set(sample_cols))} distinct titles, expected {EXPECT_SAMPLES}")

        # P3  identical order across matrices
        if sample_order_ref is None:
            sample_order_ref = sample_cols
            check(True, "column order fixed from this matrix (reference order)", "")
        else:
            check(sample_cols == sample_order_ref,
                  "column order identical to the first matrix",
                  "column order DIFFERS between matrices")

        # stream rows
        hits: dict[str, dict] = {}
        n_rows = 0
        for line in fh:
            if not line.strip():
                continue
            row_index = n_rows
            n_rows += 1
            # cheap pre-filter on the symbol field, then split fully only on a hit
            parts = line.rstrip("\n").split("\t")
            sym = norm(parts[i_symbol])
            if sym not in WANTED:
                continue
            if len(parts) != len(header):
                failures.append(f"{tag} row {row_index} ({sym}): {len(parts)} fields, "
                                f"expected {len(header)}")
                continue
            values = [norm(parts[i]) for i in sample_pos]
            hits.setdefault(sym, []).append({
                "row_index": row_index,
                "ensembl": norm(parts[i_ensembl]).split(".")[0],
                "biotype": norm(parts[i_biotype]),
                "pass_filter": norm(parts[i_filter]) if i_filter is not None else None,
                "values": values,
            })

        check(n_rows == EXPECT_FEATURE_ROWS,
              f"{n_rows} feature rows", f"{n_rows} feature rows, expected {EXPECT_FEATURE_ROWS}")

        log()
        log("--- P4  one row per gene, correct id / biotype / filter ---")
        for g in WANTED:
            rows_g = hits.get(g, [])
            if len(rows_g) != 1:
                check(False, "", f"{g}: {len(rows_g)} rows in {tag}, expected exactly 1")
                continue
            r = rows_g[0]
            ok = (r["ensembl"] == GENE_ENSEMBL[g]
                  and r["biotype"] == "protein_coding"
                  and (i_filter is None or r["pass_filter"] == "TRUE"))
            log(f"  [{'PASS' if ok else 'FAIL'}] {g:<9} row_index={r['row_index']:<6} "
                f"ensembl={r['ensembl']} biotype={r['biotype']} "
                f"TMM_pass_filter={r['pass_filter']}")
            if not ok:
                failures.append(f"{g} in {tag}: id/biotype/filter check failed")
            # P5  row index matches the frozen value, independently re-derived
            if r["row_index"] != FROZEN_ROW_INDEX[g]:
                check(False, "", f"{g}: row_index {r['row_index']} != frozen {FROZEN_ROW_INDEX[g]}")
            else:
                log(f"        row_index matches the frozen value {FROZEN_ROW_INDEX[g]} "
                    f"(independently re-derived, not copied)")

        # P7  value scale
        log()
        log("--- P7  value scale (proves the downstream log2 transform is appropriate) ---")
        scale_rows = []
        for g in WANTED:
            v = [float(x) for x in hits[g][0]["values"]]
            v_sorted = sorted(v)
            n = len(v)
            scale_rows.append({
                "gene": g,
                "min": v_sorted[0],
                "q25": v_sorted[n // 4],
                "median": v_sorted[n // 2],
                "q75": v_sorted[3 * n // 4],
                "max": v_sorted[-1],
                "mean": sum(v) / n,
                "n_zero": sum(1 for x in v if x == 0.0),
                "n_negative": sum(1 for x in v if x < 0.0),
                "all_integer": all(float(x).is_integer() for x in v),
            })
        for s in scale_rows:
            log(f"  {s['gene']:<9} min={s['min']:<12.4g} q25={s['q25']:<12.4g} "
                f"med={s['median']:<12.4g} q75={s['q75']:<12.4g} max={s['max']:<12.4g} "
                f"mean={s['mean']:<12.4g} zeros={s['n_zero']:<4} all_integer={s['all_integer']}")
        mn = min(s["min"] for s in scale_rows)
        mx = max(s["max"] for s in scale_rows)
        check(mn >= 0.0, "no negative values", f"negative values present (min={mn})")
        if tag in ("TMM", "TPM"):
            check(mx > 100,
                  f"max value {mx:.4g} > 100 -> the matrix is NOT log scaled, "
                  f"so log2(x+1) is the manuscript's log-scale analogue",
                  f"max value {mx:.4g} is too small for an unlogged TMM/TPM matrix")
            check(not scale_rows[0]["all_integer"],
                  "values are non-integer -> normalised (not raw counts)",
                  "values look like raw integers; TMM/TPM should be normalised")
        else:
            check(scale_rows[0]["all_integer"] or mx > 1000,
                  f"counts scale (max={mx:.4g})",
                  f"unexpected counts scale (max={mx:.4g})")

        extracted[tag] = {
            "header": header,
            "sample_cols": sample_cols,
            "hits": hits,
            "scale": scale_rows,
        }

# ---------------------------------------------------------------- cross-matrix
log()
log("=" * 78)
log("--- P5  gene row indices identical across the three matrices ---")
log("=" * 78)
for g in WANTED:
    idx = {tag: extracted[tag]["hits"][g][0]["row_index"] for tag in MATRICES}
    check(len(set(idx.values())) == 1,
          f"{g:<9} row_index {list(idx.values())[0]} identical across TMM / TPM / COUNTS  "
          f"(frozen {FROZEN_ROW_INDEX[g]})",
          f"{g}: row indices differ across matrices: {idx}")

log()
log("--- sample order identical across all three matrices ---")
ref = extracted["TMM"]["sample_cols"]
for tag in ("TPM", "COUNTS"):
    check(extracted[tag]["sample_cols"] == ref,
          f"{tag} column order identical to TMM",
          f"{tag} column order differs from TMM")

# ---------------------------------------------------------------- write
log()
log("=" * 78)
log("WRITING expression_3genes.csv")
log("=" * 78)
sample_cols = ref
out_path = OUT / "expression_3genes.csv"
with out_path.open("w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["representation", "gene", "ensembl_gene_id", "gene_biotype",
                "TMM_pass_filter", "feature_row_index_0based"] + sample_cols)
    for tag in ("TMM", "TPM", "COUNTS"):
        for g in WANTED:
            r = extracted[tag]["hits"][g][0]
            w.writerow([tag, g, r["ensembl"], r["biotype"], r["pass_filter"] or "",
                        r["row_index"]] + r["values"])
log(f"wrote {out_path}  ({out_path.stat().st_size} bytes)")
log(f"  rows  : {1 + 3 * 3} (header + 3 representations x 3 genes)")
log(f"  cols  : 6 annotation + {len(sample_cols)} samples")

# a machine-readable provenance block for script 03
provenance = {
    "generated": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    "matrices": {tag: {"file": spec["path"].name, "sha256": spec["sha256"]}
                 for tag, spec in MATRICES.items()},
    "genes": {g: {"ensembl": GENE_ENSEMBL[g], "frozen_row_index": FROZEN_ROW_INDEX[g]}
              for g in WANTED},
    "n_feature_rows": EXPECT_FEATURE_ROWS,
    "n_samples": len(sample_cols),
    "orientation_as_deposited": "genes x samples (289 sample columns)",
    "orientation_consumed_by_manuscript": "samples x genes (frame = genes_by_samples.T)",
    "scale": {tag: extracted[tag]["scale"] for tag in MATRICES},
}
(OUT / "expression_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
log(f"wrote {OUT / 'expression_provenance.json'}")

log()
log("=" * 78)
if failures:
    log(f"RESULT: {len(failures)} CHECK(S) FAILED")
    for f in failures:
        log(f"  - {f}")
else:
    log("RESULT: ALL PROVENANCE CHECKS PASSED")
    log("  P1  3/3 matrix SHA-256 match the frozen record")
    log("  P2  60,554 feature rows in each matrix")
    log("  P3  289 sample columns, identical order across TMM / TPM / COUNTS")
    log("  P4  1 protein_coding row per gene per matrix; TMM_pass_filter == TRUE")
    log("  P5  gene row indices identical across matrices and equal to the frozen values")
    log("  P6  matrix sample columns == the 289 manifest sample titles")
    log("  P7  TMM and TPM are NOT log scaled -> log2(x+1) is the log-scale transform")
log("=" * 78)
(OUT / "expression_extraction_log.txt").write_text("\n".join(log_lines) + "\n", encoding="utf-8")
log(f"wrote {OUT / 'expression_extraction_log.txt'}")

sys.exit(1 if failures else 0)
