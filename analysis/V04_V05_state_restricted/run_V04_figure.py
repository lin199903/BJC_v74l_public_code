"""V0.4 four-panel mechanism discrimination figure.

The four panels correspond to the four analysis blocks:
  a  Gate 0 annotation mapping rate
  b  P1-A1 patient-level state localisation  L_i
  c  P1-A2 within-TREM2_TAM coupling  K_i_T(also compared against the locked K_i from the primary analysis)
  d  P1-B A-O-K scatter plots of the three rank relationships

Reads only the output of run_V04_core.py, does not recompute any biological quantity.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(r".\data\analysis\V04_discriminating_analysis")
COL = {"GSE155698": "#1b6ca8", "GSE156405": "#d1495b", "GSE212966": "#2a9d8f"}


def load(name: str):
    p = OUT / name
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def panel_gate(ax, gate: dict) -> None:
    rows = gate.get("mapping", [])
    if not rows:
        ax.text(0.5, 0.5, "Gate 0 data unavailable", ha="center", va="center")
        ax.axis("off")
        return
    names = [r["dataset"] for r in rows]
    frac = [100 * r["mapping_fraction"] for r in rows]
    ax.bar(names, frac, color=[COL.get(n, "#888") for n in names], width=0.6)
    ax.axhline(95, ls="--", lw=1.2, color="k")
    ax.text(len(names) - 0.5, 95.6, "95% PASS threshold", ha="right", va="bottom", fontsize=8)
    for i, (f, r) in enumerate(zip(frac, rows)):
        ax.text(i, f + 0.15, f"{f:.2f}%\n{r['n_matched']}/{r['n_parent_myeloid']}",
                ha="center", va="bottom", fontsize=8)
    ax.set_ylim(90, 101.5)
    ax.set_ylabel("tumour-myeloid cells with\ndeterministic fine-state mapping (%)")
    ax.set_title(f"a  Gate 0 provenance — {gate.get('verdict', '?')}", loc="left", fontsize=10)


def panel_localisation(ax, loc: pd.DataFrame, summary: dict) -> None:
    if loc.empty:
        ax.text(0.5, 0.5, "No eligible patient", ha="center", va="center")
        ax.axis("off")
        return
    order = loc.sort_values(["dataset", "L_i"]).reset_index(drop=True)
    x = np.arange(len(order))
    ax.bar(x, order["L_i"], color=[COL.get(d, "#888") for d in order["dataset"]], width=0.7)
    ax.axhline(0, color="k", lw=1)
    mean_l = summary.get("mean_L")
    lo, hi = summary.get("ci_low"), summary.get("ci_high")
    if mean_l is not None:
        ax.axhline(mean_l, color="#444", ls="-", lw=1.2)
        ax.axhspan(lo, hi, color="#444", alpha=0.12)
        ax.text(0.02, mean_l, f"mean {mean_l:.3f} [{lo:.3f}, {hi:.3f}]",
                transform=ax.get_yaxis_transform(), fontsize=8, va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels(order["sample_id"], rotation=90, fontsize=7)
    ax.set_ylabel(r"$L_i = p_{i,T} - p_{i,O}$")
    ax.set_title(f"b  P1-A1 localisation (n={summary.get('n_patients', 0)}) — "
                 f"{summary.get('verdict', '?')}", loc="left", fontsize=10)


def panel_within_state(ax, ws: pd.DataFrame, summary: dict) -> None:
    prim = ws.loc[ws["primary_eligible"]].copy()
    if prim.empty:
        ax.text(0.5, 0.5, "No eligible patient", ha="center", va="center")
        ax.axis("off")
        return
    prim = prim.sort_values(["dataset", "K_i_T_within_state"]).reset_index(drop=True)
    x = np.arange(len(prim))
    ax.bar(x - 0.19, prim["K_i_T_within_state"], width=0.36,
           color=[COL.get(d, "#888") for d in prim["dataset"]], label=r"$K_{i,T}$ within TREM2_TAM")
    ax.bar(x + 0.19, prim["K_i_T_depth_adjusted"], width=0.36, color="#bbbbbb",
           edgecolor="#555", label=r"$K_{i,T}$ depth-adjusted (S1)")
    ax.axhline(0, color="k", lw=1)
    m = summary.get("mean")
    lo, hi = summary.get("ci_low"), summary.get("ci_high")
    if m is not None:
        ax.axhline(m, color="#a33", lw=1.4)
        ax.axhspan(lo, hi, color="#a33", alpha=0.12)
        ax.text(0.02, m, f"mean {m:.3f} [{lo:.3f}, {hi:.3f}]",
                transform=ax.get_yaxis_transform(), fontsize=8, va="bottom", color="#a33")
    ax.set_xticks(x)
    ax.set_xticklabels(prim["sample_id"], rotation=90, fontsize=7)
    ax.set_ylabel("ACP5/MMP9 contrast")
    ax.legend(fontsize=7, loc="upper left", frameon=False)
    ax.set_title(f"c  P1-A2 within-state (n={summary.get('n_patients', 0)}) — "
                 f"{summary.get('verdict', '?')}", loc="left", fontsize=10)


def panel_aok(ax, aok: pd.DataFrame, s: dict) -> None:
    if aok.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.axis("off")
        return
    pairs = [
        ("A_i_siglec15_pos_fraction", "O_i_trem2_fraction", s.get("rho_AO", {}), r"$\rho_{AO}$"),
        ("A_i_siglec15_pos_fraction", "K_i_locked_unadjusted", s.get("rho_AK", {}), r"$\rho_{AK}$"),
        ("O_i_trem2_fraction", "K_i_locked_unadjusted", s.get("rho_OK", {}), r"$\rho_{OK}$"),
    ]
    for i, (xc, yc, stat, lbl) in enumerate(pairs):
        sub = aok.dropna(subset=[xc, yc])
        ax.scatter(sub[xc], sub[yc], s=26,
                   c=[COL.get(d, "#888") for d in sub["dataset"]],
                   edgecolor="white", linewidth=0.6, zorder=3)
        r = stat.get("point")
        lo, hi = stat.get("ci_low"), stat.get("ci_high")
        txt = f"{lbl} = {r:.2f}" if r is not None and np.isfinite(r) else f"{lbl} = NA"
        if lo is not None and np.isfinite(lo):
            txt += f"\n[{lo:.2f}, {hi:.2f}]"
        ax.text(0.04, 0.96, txt, transform=ax.transAxes, fontsize=8.5, va="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#999", lw=0.6))
        ax.set_xlabel(xc.replace("_", " "), fontsize=8)
        ax.set_ylabel(yc.replace("_", " "), fontsize=8)
        ax.tick_params(labelsize=7)
        if i < 2:
            ax.set_title("", loc="left")
    ax.set_title(f"d  P1-B A–O–K structure (n={s.get('n_patients', 0)})", loc="left", fontsize=10)


def main() -> None:
    gate = load("V04_GATE_EVALUATION.json")
    loc_sum = load("_v04_p1a_localisation_summary.json")
    within = load("_v04_p1a_within_state_summary.json")
    aok_sum = load("_v04_p1b_aok_summary.json")

    loc = pd.read_csv(OUT / "V04_P1A_LOCALISATION.csv") if (OUT / "V04_P1A_LOCALISATION.csv").exists() else pd.DataFrame()
    ws = pd.read_csv(OUT / "V04_P1A_WITHIN_STATE.csv") if (OUT / "V04_P1A_WITHIN_STATE.csv").exists() else pd.DataFrame()
    aok = pd.read_csv(OUT / "V04_P1B_AOK.csv") if (OUT / "V04_P1B_AOK.csv").exists() else pd.DataFrame()

    fig = plt.figure(figsize=(13.5, 9.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.15], hspace=0.55, wspace=0.42)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1:])
    ax_c = fig.add_subplot(gs[1, 0:2])
    ax_d1 = fig.add_subplot(gs[1, 2])

    panel_gate(ax_a, gate)
    panel_localisation(ax_b, loc, loc_sum)
    panel_within_state(ax_c, ws, within.get("primary_K_i_T", {}))

    # panel d needs three sub-scatter plots: use an inset grid instead
    ax_d1.axis("off")
    inner = gs[1, 2].subgridspec(1, 3, wspace=0.75)
    for i, (xc, yc, stat, lbl) in enumerate([
        ("A_i_siglec15_pos_fraction", "O_i_trem2_fraction", aok_sum.get("rho_AO", {}), r"$\rho_{AO}$"),
        ("A_i_siglec15_pos_fraction", "K_i_locked_unadjusted", aok_sum.get("rho_AK", {}), r"$\rho_{AK}$"),
        ("O_i_trem2_fraction", "K_i_locked_unadjusted", aok_sum.get("rho_OK", {}), r"$\rho_{OK}$"),
    ]):
        axx = fig.add_subplot(inner[0, i])
        sub = aok.dropna(subset=[xc, yc]) if not aok.empty else pd.DataFrame()
        if not sub.empty:
            axx.scatter(sub[xc], sub[yc], s=20,
                        c=[COL.get(d, "#888") for d in sub["dataset"]],
                        edgecolor="white", linewidth=0.5, zorder=3)
        r = stat.get("point")
        lo, hi = stat.get("ci_low"), stat.get("ci_high")
        txt = f"{lbl}={r:.2f}" if r is not None and np.isfinite(r) else f"{lbl}=NA"
        if lo is not None and np.isfinite(lo):
            txt += f"\n[{lo:.2f},{hi:.2f}]"
        axx.set_title(txt, fontsize=8)
        axx.set_xlabel(xc.split("_")[0], fontsize=7)
        if i == 0:
            axx.set_ylabel(yc.split("_")[0], fontsize=7)
        axx.tick_params(labelsize=6)

    fig.suptitle(
        "V0.4 discriminating analysis — state localisation, within-state coupling and A–O–K structure\n"
        "Post hoc mechanistic discrimination. Patient is the inferential unit. "
        f"seed 15072026, 10,000 dataset-stratified patient bootstrap.",
        fontsize=10.5, y=0.985,
    )
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=d) for d, c in COL.items()]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=9)
    fig.subplots_adjust(top=0.88, bottom=0.09, left=0.07, right=0.985)

    target = OUT / "V04_FIG_MECHANISM_DISCRIMINATION.pdf"
    fig.savefig(target)
    fig.savefig(target.with_suffix(".png"), dpi=200)
    plt.close(fig)
    print("wrote", target)


if __name__ == "__main__":
    main()
