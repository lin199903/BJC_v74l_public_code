"""FigS5 state-level positive controls (S20) plotting script.

Reads S20_summary.json (written by run_S20_state_controls.py) and renders
FigS5.png / FigS5.svg / FigS5.pdf with three panels:
  a) cell-level sample-effect forest for TREM2 and SPP1 anchors
  b) bulk-layer forest for TREM2, SPP1, LRRC15 with decision lines |r|=0.20, 0.36
  c) aggregate-layer pooled estimates for TREM2 and SPP1

JSON contract (written by the analysis script):
{
  "cell":   {"TREM2": {"n": int, "est": float, "lo": float, "hi": float, "p": float}, "SPP1": {...}},
  "aggregate": {"TREM2": {"est": float, "lo": float, "hi": float, "p": float}, "SPP1": {...}},
  "bulk":   {"TREM2": {"cohort": {"TCGA-PAAD": r, "GSE57495": r, "GSE85916": r},
                       "pooled": {"est": float, "lo": float, "hi": float, "p": float, "I2": float}},
             "SPP1": {...},
             "LRRC15": {...}}
}
"""
import json
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load(path):
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)

def forest(ax, rows, xlabel, title, decision_lines=()):
    labels = [r[0] for r in rows]
    ests = [r[1]["est"] for r in rows]
    los = [r[1]["lo"] for r in rows]
    his = [r[1]["hi"] for r in rows]
    y = list(range(len(rows)))[::-1]
    ax.errorbar(ests, y, xerr=[[e - l for e, l in zip(ests, los)], [h - e for e, h in zip(ests, his)]],
                fmt="o", color="#1a5f7a", capsize=4, markersize=6)
    ax.axvline(0, color="0.4", lw=0.8, ls="--")
    for line, col, lab in decision_lines:
        ax.axvline(line, color=col, lw=1.0, ls=":", alpha=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    ax.grid(axis="x", alpha=0.25)

def main(json_path, out_prefix):
    d = load(json_path)

    def cell_row(anchor):
        r = d["cell"][anchor]
        return {"est": r["mean_sample_difference"], "lo": r["ci_low"], "hi": r["ci_high"]}

    def agg_row(anchor):
        r = d["aggregate"][anchor]
        return {"est": r["pooled_rho"], "lo": r["ci_low"], "hi": r["ci_high"]}

    def bulk_row(anchor):
        r = d["bulk"][anchor]
        return {"est": r["pooled_r"], "lo": r["ci_low"], "hi": r["ci_high"]}

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.6), constrained_layout=True)

    def pretty(anchor):
        return "CSF1R (2nd full case)" if anchor == "CSF1R_MYELOID" else anchor.split("_")[0]

    rows = [(a, cell_row(a)) for a in ("TREM2_TAM", "SPP1_TAM", "CSF1R_MYELOID") if a in d.get("cell", {})]
    forest(axes[0], [(pretty(a), r) for a, r in rows],
           "Sample-level score difference", "a  Cell contrast (positive vs negative)")

    rows = [(a, bulk_row(a)) for a in ("TREM2_TAM", "SPP1_TAM", "CSF1R_MYELOID", "LRRC15_CAF") if a in d.get("bulk", {})]
    forest(axes[1], [(pretty(a), r) for a, r in rows],
           "Pooled Spearman r (3 cohorts)",
           "b  Bulk layer", decision_lines=[(0.20, "#b32f2f", "|r|=0.20"), (0.36, "#7a4f1f", "|r|=0.36")])

    rows = [(a, agg_row(a)) for a in ("TREM2_TAM", "SPP1_TAM", "CSF1R_MYELOID") if a in d.get("aggregate", {})]
    forest(axes[2], [(pretty(a), r) for a, r in rows],
           "Pooled within-series rank r",
           "c  Sample aggregate (n = 20 samples)")

    fig.savefig(out_prefix + ".png", dpi=300)
    fig.savefig(out_prefix + ".svg")
    fig.savefig(out_prefix + ".pdf")
    print("wrote", out_prefix + ".png/.svg/.pdf")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
