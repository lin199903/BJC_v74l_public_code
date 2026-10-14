#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
M_ICGC_PACACA — ICGC 25K PACA-CA public RNA-seq batch download + SIGLEC15 expression extraction
=====================================================================
Data source: s3://icgc25k-open (public bucket, no authentication required)
  - Endpoint: https://object.genomeinformatics.org
  - Path:    release_28/data/PACA-CA/<DOxxx>/exp_seq/part-00000.gz
  - Format:    TSV.gz, 22 columns; col 7=Ensembl Gene ID, col 8=FPKM/RPKM expression
  - Donors:  317
  - SIGLEC15 Ensembl ID = ENSG00000197046

Strategy:
  1) List all 317 donor
  2) Concurrently download all exp_seq/part-00000.gz
  3) Streaming gzip decompression → extract the SIGLEC15 row
  4) Write wide table: rows=donors, cols=[donor_id, sig15_fpkm, n_genes_detected]
  5) Write to disk -> data/ICGC/PACA-CA/
"""
import os, sys, json, time, gzip, io, urllib.request, re, shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

ENDPOINT = 'https://object.genomeinformatics.org'
BUCKET = 'icgc25k-open'
SIG15_ENSG = 'ENSG00000197046'  # SIGLEC15
SIG15_SYMBOL = 'SIGLEC15'

ROOT = Path(r'.\data')
OUT_DIR = ROOT / 'ICGC' / 'PACA-CA'
CACHE_DIR = OUT_DIR / 'gz_cache'
CACHE_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', flush=True)

# ============== 1) List all PACA-CA donors ==============
def list_donors():
    donors = set()
    marker = None
    n_page = 0
    while True:
        n_page += 1
        url = f'{ENDPOINT}/{BUCKET}?list-type=2&delimiter=/&prefix=release_28/data/PACA-CA/&max-keys=1000'
        if marker:
            url += f'&continuation-token={marker}'
        r = urllib.request.urlopen(url, timeout=60)
        body = r.read().decode('utf-8', errors='replace')
        pfx = re.findall(r'<Prefix>(release_28/data/PACA-CA/DO(\d+))/', body)
        for _, do_id in pfx:
            donors.add(do_id)
        m = re.search(r'<NextContinuationToken>([^<]+)</NextContinuationToken>', body)
        if not m:
            break
        marker = m.group(1)
    log(f'  PACA-CA donors listed: {len(donors)} (across {n_page} pages)')
    return sorted(donors)

# ============== 2) Download the gz file of a single donor ==============
def download_one(do_id, force=False):
    gz_path = CACHE_DIR / f'PACA-CA_{do_id}_part-00000.gz'
    if gz_path.exists() and not force:
        return gz_path
    url = f'{ENDPOINT}/{BUCKET}/release_28/data/PACA-CA/DO{do_id}/exp_seq/part-00000.gz'
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=120) as r:
            data = r.read()
        gz_path.write_bytes(data)
        return gz_path
    except Exception as e:
        log(f'  !! {do_id} download failed: {e}')
        return None

# ============== 3) Streaming parse of SIGLEC15 + n_genes ==============
def parse_gz_for_sig15(gz_path, do_id):
    if gz_path is None:
        return None
    try:
        with gzip.open(gz_path, 'rt', encoding='utf-8', errors='replace') as f:
            sig15_val = None
            n_genes = 0
            n_expressed = 0
            for line in f:
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 9:
                    continue
                gene_id = cols[7]
                try:
                    val = float(cols[8])
                except (ValueError, IndexError):
                    continue
                n_genes += 1
                if val > 0:
                    n_expressed += 1
                if gene_id == SIG15_ENSG:
                    sig15_val = val
        return {
            'donor_id': f'DO{do_id}',
            'sig15_fpkm': sig15_val,
            'n_genes_total': n_genes,
            'n_genes_expressed': n_expressed,
        }
    except Exception as e:
        log(f'  !! {do_id} parse failed: {e}')
        return None

# ============== main ==============
if __name__ == '__main__':
    log('=' * 60)
    log('M_ICGC_PACACA — PACA-CA public RNA-seq SIGLEC15 extraction')
    log('=' * 60)

    # 1) Listing
    log('Step 1: Listing PACA-CA donors')
    donors = list_donors()
    log(f'  Total: {len(donors)}')

    # 2) Concurrent download (8 threads)
    log('Step 2: Concurrent download of exp_seq gz')
    t0 = time.time()
    paths = {}
    with ThreadPoolExecutor(max_workers=8) as exe:
        futs = {exe.submit(download_one, d): d for d in donors}
        done = 0
        for fut in as_completed(futs):
            d = futs[fut]
            try:
                p = fut.result()
                if p is not None:
                    paths[d] = p
            except Exception as e:
                log(f'  !! {d} exception: {e}')
            done += 1
            if done % 20 == 0:
                log(f'  Download progress {done}/{len(donors)}  ({time.time()-t0:.1f}s elapsed)')
    log(f'  Download complete: {len(paths)}/{len(donors)}  ({time.time()-t0:.1f}s elapsed)')

    # 3) Concurrent parsing (8 threads)
    log('Step 3: Concurrent parsing of SIGLEC15')
    t1 = time.time()
    rows = []
    with ThreadPoolExecutor(max_workers=8) as exe:
        futs = {exe.submit(parse_gz_for_sig15, paths.get(d), d): d for d in donors if d in paths}
        done = 0
        for fut in as_completed(futs):
            r = fut.result()
            if r:
                rows.append(r)
            done += 1
            if done % 30 == 0:
                log(f'  Parsing progress {done}/{len(paths)}  ({time.time()-t1:.1f}s elapsed)')
    log(f'  Parsing complete: {len(rows)}  ({time.time()-t1:.1f}s elapsed)')

    # 4) Write to disk
    df = pd.DataFrame(rows).sort_values('donor_id').reset_index(drop=True)
    log(f'  SIGLEC15 hits: {df["sig15_fpkm"].notna().sum()}/{len(df)}')
    if df['sig15_fpkm'].notna().any():
        log(f'  SIG15 FPKM: mean={df["sig15_fpkm"].mean():.2f}, median={df["sig15_fpkm"].median():.2f}, max={df["sig15_fpkm"].max():.2f}')

    out_csv = OUT_DIR / 'PACA-CA_SIGLEC15_expression.csv'
    df.to_csv(out_csv, index=False)
    log(f'  → {out_csv}')

    # 5) Clean up gz cache (to save space)
    log('Step 5: Cleaning up gz cache')
    cleared = 0
    for p in CACHE_DIR.glob('*.gz'):
        p.unlink()
        cleared += 1
    log(f'  cleared {cleared} files')

    # 6) Summary
    summary = {
        'n_donors': int(len(donors)),
        'n_with_sig15': int(df['sig15_fpkm'].notna().sum()),
        'n_miss_sig15': int(df['sig15_fpkm'].isna().sum()),
        'sig15_fpkm_mean': float(df['sig15_fpkm'].mean()) if df['sig15_fpkm'].notna().any() else None,
        'sig15_fpkm_median': float(df['sig15_fpkm'].median()) if df['sig15_fpkm'].notna().any() else None,
        'sig15_fpkm_max': float(df['sig15_fpkm'].max()) if df['sig15_fpkm'].notna().any() else None,
        'sig15_detected_frac': float((df['sig15_fpkm'] > 0).mean()) if df['sig15_fpkm'].notna().any() else 0.0,
        'n_genes_per_sample_median': float(df['n_genes_total'].median()),
        'n_genes_expressed_median': float(df['n_genes_expressed'].median()),
    }
    (OUT_DIR / 'PACA-CA_SIGLEC15_summary.json').write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding='utf-8'
    )
    log('  → PACA-CA_SIGLEC15_summary.json')
    log(json.dumps(summary, ensure_ascii=False, indent=2))

    log('===== M_ICGC_PACACA done =====')
