"""S20b: CSF1R second full two-level case (revision-stage decisive analysis).

Runs the CSF1R-centred pan-myeloid programme through the same cell,
aggregate and bulk layers as the locked SIGLEC15 analysis, reusing the
locked functions and the manuscript's own pan-myeloid gene set (the
post hoc myeloid adjustment set, Table S14). Expected verdict: qualified,
demonstrating the framework returns a positive full case in the same
disease and compartment where SIGLEC15 fails.

Outputs go to a writable TEMP directory; summary JSON printed to stdout.
"""
import json
import sys
import tempfile
from pathlib import Path

ANALYSIS = Path(r".\biological_validation\analysis_F_cross_scale")
sys.path.insert(0, str(ANALYSIS))
import run_S20_state_controls as S  # noqa: E402

CSF1R_MYELOID = {
    "anchor": "CSF1R",
    # Manuscript-defined pan-myeloid composition (post hoc myeloid adjustment set);
    # the anchor gene is excluded from every programme score.
    "genes": ["LST1", "AIF1", "FCER1G", "TYROBP", "CTSS", "LILRB1", "CSF1R", "SPI1"],
    "cell": True,
    "aggregate": True,
}


def main():
    out_dir = Path(tempfile.gettempdir()) / "S20b_outputs"
    out_dir.mkdir(exist_ok=True)
    report = {}
    print("[S20b] cell layer: CSF1R_MYELOID", flush=True)
    report["cell"], _, _ = S.cell_layer("CSF1R_MYELOID", CSF1R_MYELOID, out_dir)
    print("[S20b] aggregate layer: CSF1R_MYELOID", flush=True)
    report["aggregate"] = S.aggregate_layer("CSF1R_MYELOID", CSF1R_MYELOID, out_dir)
    print("[S20b] bulk layer", flush=True)
    report["bulk"] = S.bulk_layer({"CSF1R_MYELOID": CSF1R_MYELOID}, out_dir)
    with open(out_dir / "S20b_summary.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)
    print("S20B_SUMMARY_JSON_BEGIN")
    print(json.dumps(report, indent=2))
    print("S20B_SUMMARY_JSON_END")
    print("out_dir:", out_dir)


if __name__ == "__main__":
    main()
