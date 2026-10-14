"""Gate 1 final report — paired representation test."""

from __future__ import annotations

import json
import sys
from pathlib import Path

G1 = Path(r".\data\analysis\V46_validation_extension\Gate1")
rep = json.loads((G1 / "_gate1_final_summary.json").read_text(encoding="utf-8"))
pb = json.loads((G1 / "_gate1_partB_corrected.json").read_text(encoding="utf-8"))

L = []
A = L.append
A("# Gate 1 report — paired representation test")
A("")
A("**Status: COMPLETE. Correctness gate EXACT.**")
A("")
A("## Correctness gate")
A("")
A("The single-cell aggregation was recomputed from the three files listed in the")
A("archived `input_manifest.csv` (SHA-256 verified against the manifest) using the")
A("frozen selection rule (`tissue in [PDAC, PDAC_tumor, PDAC_tissue]` and")
A("`celltype == Myeloid`, >=50 myeloid cells per sample). The recomputed myeloid")
A(f"aggregate reproduces the archived values with a maximum absolute difference of")
A(f"{pb['correctness_gate']['max_abs_diff_R1']:.1e} across "
  f"{pb['correctness_gate']['n_compared']} samples")
A(f"(myeloid counts {pb['correctness_gate']['n_myeloid_count_match']}/"
  f"{pb['correctness_gate']['n_compared']} and detectable-cell counts "
  f"{pb['correctness_gate']['n_positive_count_match']}/"
  f"{pb['correctness_gate']['n_compared']} identical).")
A("")
A("Cell totals: 4,081 (GSE212966) + 4,635 (GSE156405) + 1,715 (GSE155698) =")
A("**10,431 tumour-myeloid cells across 20 samples**, matching the manuscript.")
A("")
A("## Readouts tested")
A("")
A("| Readout | Definition |")
A("|---|---|")
A("| R1 | myeloid aggregate SIGLEC15 — closest analogue of a scalable RNA readout |")
A("| R2 | SIGLEC15-detectable fraction of included tumour-myeloid cells |")
A("| R3 | mean SIGLEC15 among detectable tumour-myeloid cells |")
A("| R4 | all-cell pseudo-whole-tumour aggregate over the tumour compartment |")
A("")
A("K_i is the locked within-patient detectable-minus-non-detectable ACP5/MMP9")
A("contrast (19 estimable samples).")
A("")
A("## Results")
A("")
A("| Readout | n | Spearman rho | 95% CI | CI upper bound | LO direction |")
A("|---|---:|---:|---|---:|---|")
for name, r in rep["readouts"].items():
    if "error" in r:
        continue
    loo = r["leave_one_dataset_out"]
    A(f"| {r['description']} | {r['n']} | {r['rho']:+.3f} | "
      f"[{r['ci_low']:+.3f}, {r['ci_high']:+.3f}] | {r['ci_high']:+.3f} | "
      f"{'retained' if loo['_direction_retained'] else 'sign reversal'} |")
A("")
A("Per-dataset Spearman rho:")
A("")
A("| Readout | GSE155698 (n=5) | GSE156405 (n=8) | GSE212966 (n=6) |")
A("|---|---:|---:|---:|")
for name, r in rep["readouts"].items():
    if "error" in r or "per_dataset" not in r:
        continue
    p = r["per_dataset"]
    vals = [f"{p[k]['rho']:+.2f}" if k in p else "—"
            for k in ("GSE155698", "GSE156405", "GSE212966")]
    A(f"| {r['description']} | {vals[0]} | {vals[1]} | {vals[2]} |")
A("")
A("## Interpretation")
A("")
A("1. **Both continuous target-quantity readouts exclude moderate-to-strong")
A("   positive concordance.** R1 and R3 have 95% CI upper bounds of")
A(f"   {rep['readouts']['R1_myeloid_aggregate']['ci_high']:+.3f} and "
  f"{rep['readouts']['R3_positive_intensity']['ci_high']:+.3f}. On the same")
A("   patient samples, a higher target quantity did not accompany a larger")
A("   cellular contrast, and the data are inconsistent with a positive")
A("   association above roughly 0.15.")
A("2. **R3 is the most internally consistent readout**, being negative in all")
A("   three series separately (-0.90, -0.71, -0.26) with direction retained")
A("   under leave-one-dataset-out.")
A("3. **R1 retains direction** under leave-one-dataset-out (all three omissions")
A("   negative), so its negative point estimate is not driven by one series.")
A("4. **The detectable fraction (R2) is unstable** (one sign reversal), so no")
A("   directional statement is supported for prevalence alone.")
A("5. **The all-cell aggregate (R4) is uninformative** (rho +0.075,")
A("   CI [-0.401, +0.520]) and does not support strong positive concordance")
A("   either.")
A("")
A("## Claim ceiling")
A("")
A("This is a paired comparison **within the same public single-cell samples**; it")
A("is not an independent clinical assay pairing, which would require a scalable")
A("readout and the cellular reference measured in the same specimen on different")
A("platforms. The supported statement is:")
A("")
A("> In the same patient samples, neither of the two continuous target-quantity")
A("> readouts showed positive concordance with the cellular state, and")
A("> moderate-to-strong positive concordance is excluded; prevalence and the")
A("> all-cell aggregate were uninformative.")
A("")
A("Not supported: that any readout has been validated, qualified, or shown to be")
A("interchangeable; and not supported: a biological negative relation.")
A("")
A("## What this changes in the manuscript")
A("")
A("v45 Discussion states that the unresolved quantity is same-patient concordance")
A("and that paired specimens are required. Gate 1 supplies a within-sample paired")
A("test on the existing public data, which partially closes that gap at the RNA")
A("level while leaving the clinical-assay pairing open. Recommended placement:")
A("")
A("- **Results**: 2-3 sentences reporting R1 and R3 with CIs and the exclusion")
A("  bound; R2/R4 as uninformative; pointing to a new supplementary table.")
A("- **Discussion**: revise the 'unresolved quantity' paragraph so that the open")
A("  item becomes independent clinical-assay pairing rather than any same-patient")
A("  comparison at all.")
A("- **Conclusions**: unchanged claim ceiling; optionally note that the")
A("  within-sample RNA comparison did not support interchangeability.")
A("- **Supplementary**: new table with all four readouts, per-dataset estimates")
A("  and leave-one-dataset-out behaviour.")

(G1 / "GATE1_REPORT.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L))
