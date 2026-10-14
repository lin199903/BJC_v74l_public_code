"""E1: eligibility check of the three external candidate datasets.

The author's brief specifies six admission criteria and says a dataset that fails them is
not used:

  1. three-gene coverage   SIGLEC15, ACP5 and MMP9 all measurable
  2. demonstrable independence  patient and sample identifiers resolvable, and no overlap
     with the series already analysed (GSE155698, GSE156405, GSE212966)
  3. cell-level measurement  expression data plus metadata sufficient for cell identity
  4. technical quality       background, depth and low-expression signal inspectable
  5. treatment context       untreated / treated and primary / metastatic stated
  6. spatial coordinates     for spatial data, cell boundaries or molecule coordinates

This queries GEO deterministically and inspects file listings. It downloads gene panels or
barcode files where they are small enough to settle a criterion, and it does not download
expression matrices.

Candidates: GSE217845 (scRNA-seq), GSE310352 (CosMx), GSE277782 (CosMx).
"""

from __future__ import annotations
import os

import gzip
import io
import json
import re
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

V47 = Path(os.path.join(os.environ.get("SIGLEC15_SUBMISSION_ROOT", r"E:\Tex\HMlong\npjPrecOncol\FINAL_SUBMISSION_v6"), "BMC_MEDICINE_SUBMISSION_v47"))
OUT = V47 / "QA" / "gate4" / "E1_ELIGIBILITY.json"
CACHE = V47 / "QA" / "gate4" / "_e1_cache"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

CANDIDATES = {
    "GSE217845": {"kind": "scRNA-seq", "need_spatial": False},
    "GSE310352": {"kind": "CosMx", "need_spatial": True},
    "GSE277782": {"kind": "CosMx", "need_spatial": True},
}
ORIGINAL = {"GSE155698", "GSE156405", "GSE212966"}
GENES = ["SIGLEC15", "ACP5", "MMP9"]


def fetch(url: str, tries: int = 3, timeout: int = 120) -> bytes:
    last: Exception | None = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "research-script"})
            return urllib.request.urlopen(req, timeout=timeout).read()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 + 2 * k)
    raise RuntimeError(f"fetch failed: {url}") from last


def get_json(url: str) -> dict:
    return json.loads(fetch(url).decode("utf-8", "replace"))


def listing(gse: str) -> list[dict]:
    """FTP suppl listing plus filelist.txt sizes."""
    short = gse[:-3] + "nnn"
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{short}/{gse}/suppl/"
    sizes: dict[str, int] = {}
    try:
        fl = fetch(base + "filelist.txt").decode("utf-8", "replace")
        for line in fl.splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) >= 4 and parts[3].isdigit():
                sizes[parts[1].strip()] = int(parts[3])
    except Exception:  # noqa: BLE001
        pass
    html = fetch(base).decode("utf-8", "replace")
    names = [f for f in re.findall(r'href="([^"]+)"', html)
             if not f.startswith(("?", "/", "http")) and not f.endswith("/")]
    return [{"name": n, "bytes": sizes.get(n, 0)} for n in names]


def first_tar_members(gse: str, name: str, limit: int = 40) -> tuple[list[str], bytes | None]:
    """Read the member list of the first inner archive; return one small file's bytes."""
    short = gse[:-3] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{short}/{gse}/suppl/{name}"
    CACHE.mkdir(parents=True, exist_ok=True)
    local = CACHE / name
    if not local.exists():
        local.write_bytes(fetch(url, timeout=600))
    members: list[str] = []
    payload: bytes | None = None
    try:
        with tarfile.open(local, "r:*") as t:
            for m in t.getmembers():
                members.append(f"{m.name} ({m.size:,} B)")
            # pull one barcode/feature-like file if small
            for m in t.getmembers():
                if m.size and m.size < 3_000_000 and re.search(
                        r"barcode|feature|gene|panel|meta", m.name, re.I):
                    payload = t.extractfile(m).read()
                    break
    except Exception as e:  # noqa: BLE001
        members = [f"<read error: {type(e).__name__}>"]
    return members, payload


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    report: dict = {"generated": time.strftime("%Y-%m-%d %H:%M"),
                    "criteria": ["three-gene coverage", "demonstrable independence",
                                 "cell-level measurement", "technical quality",
                                 "treatment context", "spatial coordinates"],
                    "candidates": {}}

    for gse, spec in CANDIDATES.items():
        print(f"\n{'=' * 72}\n=== {gse}  ({spec['kind']}) ===")
        e: dict = {"kind": spec["kind"]}

        search = get_json(
            f"{EUTILS}/esearch.fcgi?db=gds&term={gse}[ACCN]+AND+gse[ETYP]&retmode=json")
        uids = search["esearchresult"]["idlist"]
        if not uids:
            print("  !! no UID")
            continue
        rec = get_json(
            f"{EUTILS}/esummary.fcgi?db=gds&id={uids[0]}&retmode=json")["result"][uids[0]]

        e["title"] = rec.get("title")
        e["n_samples"] = rec.get("n_samples")
        e["pdat"] = rec.get("pdat")
        e["gdstype"] = rec.get("gdstype")
        print(f"  {rec.get('accession')}  n={e['n_samples']}  {e.get('pdat')}")
        print(f"  title: {(e['title'] or '')[:92]}")

        rel = rec.get("relation", [])
        if isinstance(rel, str):
            rel = [rel]
        e["has_sra"] = any(re.search(r"SRA|SRP|PRJ|BioProject", str(r)) for r in rel)
        print(f"  SRA/BioProject: {e['has_sra']}")

        samples = [s for s in rec.get("samples", []) if isinstance(s, dict)]
        e["sample_titles"] = [s.get("title", "") for s in samples]
        e["sample_accessions"] = [s.get("accession", "") for s in samples]
        print(f"  samples: {len(samples)}")
        for s in samples[:12]:
            print(f"     {s.get('accession',''):<12} {(s.get('title') or '')[:60]}")

        e["files"] = listing(gse)
        print(f"  suppl files: {len(e['files'])}")
        for f in e["files"][:14]:
            print(f"     {f['bytes']:>14,} B  {f['name'][:74]}")

        # probe an inner archive for a gene panel or barcode list
        e["panel_probe"] = None
        tars = [f["name"] for f in e["files"] if re.search(r"\.tar(\.gz)?$", f["name"])
                and f["bytes"] < 400_000_000]
        if tars:
            print(f"  probing inner archive: {tars[0]}")
            members, payload = first_tar_members(gse, tars[0])
            e["panel_probe"] = {"archive": tars[0], "members": members[:30]}
            for m in members[:14]:
                print(f"     {m}")
            if payload:
                txt = payload.decode("utf-8", "replace")
                head = txt[:400].replace("\n", " | ")
                print(f"  small file head: {head[:300]}")
                hits = {g: (g.upper() in txt.upper()) for g in GENES}
                e["panel_probe"]["gene_hits"] = hits
                print(f"  three-gene presence in that file: {hits}")
                n_lines = txt.count("\n")
                e["panel_probe"]["n_lines"] = n_lines
                print(f"  lines in that file: {n_lines:,}")

        e["overlap_with_original"] = sorted(ORIGINAL & {gse})
        report["candidates"][gse] = e
        time.sleep(1)

    OUT.write_text(json.dumps(report, ensure_ascii=False, indent=1), encoding="utf-8")

    print("\n" + "=" * 72)
    print("E1 概要")
    for gse, e in report["candidates"].items():
        print(f"  {gse:<12} {e.get('kind'):<9} n={e.get('n_samples')}  "
              f"SRA={e.get('has_sra')}  files={len(e.get('files', []))}")
    print(f"\nwritten: {OUT.relative_to(V47)}")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
