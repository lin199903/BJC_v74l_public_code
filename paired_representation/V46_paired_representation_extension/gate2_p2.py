"""Gate 2 P2 — CPTAC-PDAC paired RNA vs downstream osteoclast-state protein.

Gate P1 (target protein) is already settled: SIGLEC15 protein is undetected in
CPTAC-PDAC. Gate P2 asks whether the downstream state proteins are quantifiable
in the same patients, which would allow a paired tumour-RNA(SIGLEC15) versus
protein(state) test.

Inputs: LinkedOmics .cct matrices in data/analysis/proteogenomic.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PG = Path(r".\data\analysis\proteogenomic")
OUT = Path(r".\data\analysis\V46_validation_extension\Gate2")
OUT.mkdir(parents=True, exist_ok=True)

TARGET = "SIGLEC15"
# frozen osteoclast panel (analysis_contract.json) + receptor + controls
DOWNSTREAM = ["ACP5", "MMP9", "CTSK", "ATP6V0D2", "DCSTAMP", "OCSTAMP",
              "ITGB3", "CALCR", "NFATC1", "TREM2", "CTSB", "MMP14"]
CONTROLS = ["EPCAM", "KRT8", "MUC1", "CD68", "AIF1", "CSF1R", "TYROBP"]


def read_cct(path: Path) -> pd.DataFrame:
    """LinkedOmics .cct: comment lines start with '#'; one header row follows."""
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        lines = []
        for line in fh:
            if line.startswith("#"):
                continue
            lines.append(line.rstrip("\n"))
    if not lines:
        return pd.DataFrame()
    header = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[1:] if ln.strip()]
    # first column carries the feature id
    idx = [r[0] for r in rows]
    data = [r[1:] for r in rows]
    df = pd.DataFrame(data, index=idx, columns=header[1:])
    return df


def detection(df: pd.DataFrame, gene: str) -> dict | None:
    if gene not in df.index:
        return None
    s = pd.to_numeric(df.loc[gene], errors="coerce")
    return {
        "n_samples": int(s.notna().sum()),
        "detection_rate": float(s.notna().mean()),
        "median": float(s.median()) if s.notna().any() else None,
    }


def main() -> None:
    rep: dict = {"inputs": {}, "genes": {}}

    prot_path = PG / "proteome_tumor.cct"
    rna_path = PG / "rna_tumor.cct"
    for p in (prot_path, rna_path):
        rep["inputs"][p.name] = {"bytes": p.stat().st_size if p.exists() else None}

    prot = read_cct(prot_path)
    rna = read_cct(rna_path)
    rep["inputs"]["proteome_shape"] = list(prot.shape)
    rep["inputs"]["rna_shape"] = list(rna.shape)
    rep["inputs"]["proteome_index_head"] = list(prot.index[:5])
    rep["inputs"]["rna_index_head"] = list(rna.index[:5])
    rep["inputs"]["shared_samples"] = int(len(set(prot.columns) & set(rna.columns)))

    for g in [TARGET] + DOWNSTREAM + CONTROLS:
        rep["genes"][g] = {
            "in_proteome": g in prot.index,
            "in_rna": g in rna.index,
            "protein": detection(prot, g),
            "rna": detection(rna, g),
        }

    # Gate P2 decision
    ds_avail = {g: rep["genes"][g]["protein"] for g in DOWNSTREAM
                if rep["genes"][g]["protein"]}
    usable = [g for g, v in ds_avail.items() if v and v["detection_rate"] >= 0.50]
    rep["gate_P2"] = {
        "downstream_quantifiable_ge50pct": usable,
        "n_quantifiable": len(usable),
        "acp5_and_mmp9_both_usable": ("ACP5" in usable and "MMP9" in usable),
        "verdict": ("PROCEED" if ("ACP5" in usable and "MMP9" in usable) else "STOP"),
    }

    # If ACP5 and MMP9 are usable, run the paired test:
    # tumour RNA(SIGLEC15) vs protein(state), same patients.
    if rep["gate_P2"]["verdict"] == "PROCEED":
        common = [c for c in rna.columns if c in set(prot.columns)]
        pro = prot[common].apply(pd.to_numeric, errors="coerce")
        rn = rna[common].apply(pd.to_numeric, errors="coerce")

        def zscore(s: pd.Series) -> pd.Series:
            v = s.astype(float)
            return (v - v.mean()) / v.std(ddof=1)

        state_cols = [g for g in (["ACP5", "MMP9"] + [x for x in usable
                                                      if x not in ("ACP5", "MMP9")])]
        prot_state = pd.concat([zscore(pro.loc[g]) for g in state_cols
                                if g in pro.index], axis=1).mean(axis=1)
        rna_sig = rn.loc[TARGET] if TARGET in rn.index else None

        if rna_sig is not None:
            pair = pd.DataFrame({"rna_SIGLEC15": rna_sig,
                                 "prot_state": prot_state}).dropna()
            rho, p = stats.spearmanr(pair["rna_SIGLEC15"], pair["prot_state"])
            rep["paired_test"] = {
                "n_patients": int(len(pair)),
                "state_proteins_used": state_cols,
                "rho_RNA_SIGLEC15_vs_protein_state": float(rho),
                "p": float(p),
                "note": ("Tumour-level RNA versus downstream protein state in the "
                         "same patients; the target protein itself is undetected."),
            }
            pair.to_csv(OUT / "GATE2_paired_rna_protein_state.csv", index=False)

    (OUT / "_gate2_summary.json").write_text(
        json.dumps(rep, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    print(json.dumps(rep, indent=2, ensure_ascii=False, default=str))


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
