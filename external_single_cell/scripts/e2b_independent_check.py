"""E2-B: independent single-cell association check on GSE217845.

Frozen rules are in QA\\gate4\\E2_FROZEN_RULES.md and were written before this ran.

Design, in the order the frozen rules set it:

  * primary pool  = all 10 primary-tumour samples (the `tissue` field is `Primary tumor`
    for every one of them, and the 10 patient IDs are distinct)
  * strata        = PDAC_* unsorted (6), LPDAC_* unsorted (3), LPDAC_15 sorted (1)
  * blood samples = excluded from the tumour myeloid contrast
  * myeloid cells = marker rule, declared in advance, because the deposit has no cell
    annotation (LYZ, AIF1, TYROBP, CD68, CD163, CSF1R, ITGAM, FCER1G, C1QA, C1QB)
  * score         = ACP5 and MMP9 z-standardised within biological sample, then averaged
  * contrast Ki   = mean score in SIGLEC15-detectable minus non-detectable myeloid cells
  * summary       = equal weight per sample, sample-level bootstrap interval

It downloads the public filtered matrices, which is the same data level as the three
series in the manuscript. It therefore cannot provide ambient-RNA correction, and the
report must say so.

Everything that is not estimable is recorded, not dropped.
"""

from __future__ import annotations
import os

import csv
import gzip
import io
import json
import math
import random
import re
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np

V47 = Path(os.path.join(os.environ.get("SIGLEC15_SUBMISSION_ROOT", r"E:\Tex\HMlong\npjPrecOncol\FINAL_SUBMISSION_v6"), "BMC_MEDICINE_SUBMISSION_v47"))
WORK = V47 / "QA" / "gate4" / "_e2"
OUT = V47 / "QA" / "gate4" / "E2B_RESULT.json"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/samples"

SEED = 20260918
N_BOOT = 10000
MYELOID_MARKERS = ["LYZ", "AIF1", "TYROBP", "CD68", "CD163", "CSF1R", "ITGAM",
                   "FCER1G", "C1QA", "C1QB"]
PROGRAM = ["ACP5", "MMP9"]
TARGET = "SIGLEC15"

# GSM -> (patient, label, sorted)
SAMPLES = {
    "GSM6727542": ("LPDAC_15", "tumour", True),
    "GSM6727543": ("LPDAC_25", "tumour", False),
    "GSM6727544": ("LPDAC_26", "tumour", False),
    "GSM6727545": ("LPDAC_30", "tumour", False),
    "GSM6727546": ("PDAC_47", "tumour", False),
    "GSM6727547": ("PDAC_48", "tumour", False),
    "GSM6727548": ("PDAC_50", "tumour", False),
    "GSM6727549": ("PDAC_51", "tumour", False),
    "GSM6727550": ("PDAC_55", "tumour", False),
    "GSM6727551": ("PDAC_60", "tumour", False),
}


def fetch(url: str, tries: int = 3, timeout: int = 600) -> bytes:
    last: Exception | None = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "research-script"})
            return urllib.request.urlopen(req, timeout=timeout).read()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(3 + 3 * k)
    raise RuntimeError(f"fetch failed: {url}") from last


def gsm_suppl(gsm: str) -> list[str]:
    stem = gsm[:-3] + "nnn"
    html = fetch(f"{FTP}/{stem}/{gsm}/suppl/", timeout=120).decode("utf-8", "replace")
    return [f for f in re.findall(r'href="([^"]+)"', html)
            if not f.startswith(("?", "/", "http")) and not f.endswith("/")]


def download(gsm: str) -> dict[str, Path]:
    WORK.mkdir(parents=True, exist_ok=True)
    d = WORK / gsm
    d.mkdir(exist_ok=True)
    paths: dict[str, Path] = {}
    for name in gsm_suppl(gsm):
        low = name.lower()
        key = ("barcodes" if "barcode" in low else
               "features" if "feature" in low or "gene" in low else
               "matrix" if "matrix" in low else None)
        if not key:
            continue
        out = d / name
        if not out.exists():
            out.write_bytes(fetch(f"{FTP}/{gsm[:-3]}nnn/{gsm}/suppl/{name}"))
        paths[key] = out
    return paths


def read_mtx(paths: dict[str, Path]) -> tuple[list[str], list[str], np.ndarray]:
    """Read a 10x matrix market triplet into a dense cell x gene array (samples are small)."""
    with gzip.open(paths["features"], "rt", errors="replace") as f:
        feats = [l.split("\t") for l in f if l.strip()]
    genes = [(c[1].strip() if len(c) > 1 else c[0].strip()) for c in feats]
    with gzip.open(paths["barcodes"], "rt", errors="replace") as f:
        bcs = [l.strip() for l in f if l.strip()]
    rows, cols, vals = [], [], []
    with gzip.open(paths["matrix"], "rt", errors="replace") as f:
        for line in f:
            if line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            try:
                r, c, v = int(parts[0]), int(parts[1]), float(parts[2])
            except ValueError:
                continue  # header
            rows.append(r - 1)
            cols.append(c - 1)
            vals.append(v)
    M = np.zeros((len(bcs), len(genes)), dtype=np.float32)
    M[np.array(cols), np.array(rows)] = np.array(vals, dtype=np.float32)
    return bcs, genes, M


def analyse(gsm: str, patient: str, sorted_sample: bool) -> dict:
    rec: dict = {"gsm": gsm, "patient": patient, "sorted": sorted_sample}
    try:
        paths = download(gsm)
    except Exception as e:  # noqa: BLE001
        rec["status"] = "download_failed"
        rec["reason"] = f"{type(e).__name__}: {e}"
        return rec
    if len(paths) < 3:
        rec["status"] = "incomplete_files"
        rec["reason"] = f"only {sorted(paths)}"
        return rec

    bcs, genes, M = read_mtx(paths)
    gi = {g.upper(): i for i, g in enumerate(genes)}
    need = [TARGET] + PROGRAM + MYELOID_MARKERS
    missing = [g for g in need if g not in gi]
    rec["n_cells_raw"] = len(bcs)
    rec["n_genes"] = len(genes)
    rec["genes_missing"] = missing
    if any(g not in gi for g in [TARGET] + PROGRAM):
        rec["status"] = "gene_missing"
        return rec

    total = M.sum(axis=1)
    n_genes = (M > 0).sum(axis=1)
    keep = (n_genes >= 200) & (total >= 500)
    rec["n_cells_qc"] = int(keep.sum())
    if keep.sum() < 50:
        rec["status"] = "too_few_cells_after_qc"
        return rec
    Mq = M[keep]
    bq = [b for b, k in zip(bcs, keep) if k]

    # myeloid call, declared in advance: number of positive markers >= 3
    mk = [gi[g] for g in MYELOID_MARKERS if g in gi]
    pos = (Mq[:, mk] > 0).sum(axis=1)
    my = pos >= 3
    rec["n_myeloid"] = int(my.sum())
    if my.sum() < 20:
        rec["status"] = "too_few_myeloid"
        return rec
    Mm = Mq[my]

    # library-size normalisation on the myeloid subset, then log1p
    lib = Mm.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1
    Mn = np.log1p(Mm / lib * 1e4)

    # programme score: z within sample, then mean of the two genes
    z = np.zeros(len(PROGRAM), dtype=object)
    parts = []
    for g in PROGRAM:
        v = Mn[:, gi[g]]
        sd = v.std()
        parts.append((v - v.mean()) / sd if sd > 0 else v * 0.0)
    score = np.mean(np.vstack(parts), axis=0)

    det = Mm[:, gi[TARGET]] > 0
    rec["n_detectable"] = int(det.sum())
    rec["n_nondetectable"] = int((~det).sum())
    if det.sum() == 0:
        rec["status"] = "zero_detectable"
        return rec
    if (~det).sum() == 0:
        rec["status"] = "zero_nondetectable"
        return rec

    rec["status"] = "estimable"
    rec["contrast"] = float(score[det].mean() - score[~det].mean())
    rec["prevalence"] = float(det.mean())
    rec["mean_score_detectable"] = float(score[det].mean())
    rec["mean_score_nondetectable"] = float(score[~det].mean())
    return rec


def main() -> None:
    rng = random.Random(SEED)
    out: dict = {"seed": SEED, "n_boot": N_BOOT, "samples": {}}
    for gsm, (patient, kind, srt) in SAMPLES.items():
        print(f"  {gsm}  {patient:<10} sorted={srt}", end="", flush=True)
        rec = analyse(gsm, patient, srt)
        out["samples"][gsm] = rec
        if rec.get("status") == "estimable":
            print(f"  -> Ki={rec['contrast']:+.3f}  myeloid={rec['n_myeloid']:,}  "
                  f"detectable={rec['n_detectable']}")
        else:
            print(f"  -> {rec.get('status')}  {rec.get('reason','')}")

    est = [r for r in out["samples"].values() if r.get("status") == "estimable"]
    out["n_estimable"] = len(est)
    out["n_tumour_samples"] = len(SAMPLES)
    if est:
        k = [r["contrast"] for r in est]
        out["mean_contrast"] = float(np.mean(k))
        out["n_positive"] = int(sum(1 for x in k if x > 0))
        boots = []
        for _ in range(N_BOOT):
            s = [k[rng.randrange(len(k))] for _ in range(len(k))]
            boots.append(float(np.mean(s)))
        boots.sort()
        out["ci95"] = [boots[int(0.025 * len(boots))], boots[int(0.975 * len(boots))]]
        for name, sel in (("PDAC_unsorted", lambda r: r["patient"].startswith("PDAC")),
                          ("LPDAC_unsorted", lambda r: r["patient"].startswith("LPDAC") and not r["sorted"]),
                          ("without_sorted", lambda r: not r["sorted"]),
                          ("with_sorted", lambda r: True)):
            sub = [r["contrast"] for r in est if sel(r)]
            out[f"stratum_{name}"] = {"n": len(sub),
                                      "mean": float(np.mean(sub)) if sub else None}

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")

    print()
    print("=" * 68)
    print(f"可估计 {out['n_estimable']} / {out['n_tumour_samples']} 个肿瘤样本")
    if est:
        print(f"样本等权平均差值: {out['mean_contrast']:+.3f}  "
              f"95% bootstrap CI {out['ci95'][0]:+.3f} 至 {out['ci95'][1]:+.3f}")
        print(f"正向样本: {out['n_positive']} / {out['n_estimable']}")
        for name in ("PDAC_unsorted", "LPDAC_unsorted", "without_sorted", "with_sorted"):
            s = out.get(f"stratum_{name}", {})
            if s.get("n"):
                print(f"  {name:<16} n={s['n']}  均值 {s['mean']:+.3f}")
    print(f"\nwritten: {OUT.relative_to(V47)}")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
