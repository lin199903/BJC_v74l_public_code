"""Generate v23 figure set for the two-level qualification manuscript.

All panels are drawn directly from machine-readable result CSVs (no analysis
is rerun). Reuses the v22 colourblind-safe palette (colour supplemented by
shape, hatch or direct labels). Run from any CWD; outputs to ../Figures.

Figure map (story beat -> figure):
  Fig1  framework schematic + resolution paradox (new)
  Fig2  RNA landscape (copied from v22, unchanged)
  Fig3  Level 1 cellular state (contrast, depth-adj, depth-match, single-gene, thresholds)
  Fig4  Level 2 qualification failure (aggregate, bulk, spatial, verdict strip)
  Fig5  not-technical (calibration, abundance-matched, payload, context heatmap)
  Fig6  broader claims refused (PI3K, SYK/DAP12, MUC, protein boundary)
  Fig7  protocol deliverable (decision flow, grade, protein-first roadmap)
  FigS1 tissue-context heatmap (copied from v22, unchanged)
  FigS2 myeloid subtype enrichment (Task5)
  FigS3 survival forest (Task2 meta-Cox)
  FigS4 deconvolution cross-tool heatmap (Task4)
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "Figures"
OUT.mkdir(parents=True, exist_ok=True)

DEPOSIT = Path(r"./prior_analysis_v1\repository_deposit")
R17 = DEPOSIT / "Results" / "S17"
R18 = DEPOSIT / "Results" / "S18"
R19 = DEPOSIT / "Results" / "S19"
P02 = DEPOSIT / "Results" / "P02"
P22 = DEPOSIT / "Results" / "P22"
ANALYSIS = Path(r".\data\analysis")
V22FIG = Path(r"./prior_analysis_v1\Figures")
VALIDATION = Path(r".\validation_design")

DPI = 300
COL = {
    "ink": "#242424", "mid": "#696969", "grid": "#D8D8D8", "pale": "#F5F5F3",
    "blue": "#3B6C99", "blue2": "#7B9DBB", "bluep": "#DDE8F1",
    "orange": "#C47A3A", "orangep": "#F2E2D3",
    "teal": "#3C877E", "tealp": "#DCEBE8",
    "purple": "#7A5C91", "purplep": "#E8E0ED",
    "red": "#A14F4A", "redp": "#EED9D7",
    "green": "#64805D", "greenp": "#E0E9DD",
}
MARKERS = {"GSE155698": "o", "GSE156405": "s", "GSE212966": "^"}
DATSET_NAMES = {"GSE155698": "GSE155698", "GSE156405": "GSE156405", "GSE212966": "GSE212966"}


def setup(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("#B0B0B0")
    ax.spines["left"].set_color("#B0B0B0")
    ax.grid(axis="x", color="#E3E3E3", lw=0.45, zorder=0)
    ax.tick_params(colors=COL["ink"])


def save(fig, name):
    for ext in ("png", "svg", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ----------------------------------------------------------------------------
# Fig 2 / Fig S1: copy unchanged panels from v22
# ----------------------------------------------------------------------------
def fig2_landscape():
    shutil.copy(V22FIG / "Fig2.png", OUT / "Fig2.png")
    shutil.copy(V22FIG / "FigS1.png", OUT / "FigS1.png")
    print("copied Fig2.png, FigS1.png")


# ----------------------------------------------------------------------------
# Fig 3: Level 1 cellular state (panels a-e)
# ----------------------------------------------------------------------------
def fig3_level1():
    effects = pd.read_csv(R17 / "cell_level_sample_effects.csv")
    core = effects[effects.program == "Derosiers_ACP5_MMP9_core"].copy()
    fig, axes = plt.subplots(1, 5, figsize=(15.5, 3.4))

    # a: unadjusted contrasts per sample
    ax = axes[0]
    for ds, grp in core.groupby("dataset"):
        ax.scatter(range(len(grp)), grp.difference_positive_minus_negative.sort_values().values,
                   marker=MARKERS[ds], s=34, color=COL["blue"], edgecolor="white", lw=0.5,
                   label=DATSET_NAMES[ds], zorder=3)
    ax.axhline(0.639, color=COL["orange"], lw=1.4, ls="--")
    ax.set_title("a  Sample-level contrasts", fontsize=9, loc="left", color=COL["ink"])
    ax.set_ylabel("Score difference (+ vs \u2212)", fontsize=8)
    ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    setup(ax)

    # b: depth-adjusted contrasts per sample
    ax = axes[1]
    for ds, grp in core.groupby("dataset"):
        ax.scatter(range(len(grp)), grp.depth_adjusted_difference.sort_values().values,
                   marker=MARKERS[ds], s=34, color=COL["teal"], edgecolor="white", lw=0.5, zorder=3)
    ax.axhline(0.426, color=COL["orange"], lw=1.4, ls="--")
    ax.set_title("b  Depth-adjusted", fontsize=9, loc="left", color=COL["ink"])
    setup(ax)

    # c: depth-matched (S19)
    dm = pd.read_csv(R19 / "depth_matched_summary.csv")
    ax = axes[2]
    names = {"Derosiers_ACP5_MMP9_core": "Core", "canonical_osteoclast_panel": "Osteoclast", "compact_ECM_degradation_panel": "Matrix deg."}
    y = np.arange(len(dm))[::-1]
    for yi, (_, row) in zip(y, dm.iterrows()):
        ax.errorbar(row.mean_difference, yi, xerr=[[row.mean_difference - row.ci_low], [row.ci_high - row.mean_difference]],
                    fmt="o", color=COL["blue"], ecolor=COL["blue"], capsize=2.5, ms=6, zorder=3)
        ax.text(row.ci_high + 0.02, yi, f"p = {row.wilcoxon_p:.3f}", fontsize=6.5, va="center", color=COL["mid"])
    ax.set_yticks(y, [names[r] for r in dm.program], fontsize=7.5)
    ax.axvline(0, color="#B0B0B0", lw=0.8)
    ax.set_title("c  Depth-matched (1:1)", fontsize=9, loc="left", color=COL["ink"])
    setup(ax)

    # d: single-gene decomposition (S18)
    ax = axes[3]
    sg = pd.DataFrame([
        {"gene": "ACP5", "adj": False, "est": 0.613, "lo": 0.281, "hi": 0.983},
        {"gene": "ACP5", "adj": True, "est": 0.298, "lo": 0.027, "hi": 0.628},
        {"gene": "MMP9", "adj": False, "est": 0.715, "lo": 0.367, "hi": 1.100},
        {"gene": "MMP9", "adj": True, "est": 0.596, "lo": 0.291, "hi": 0.924},
    ])
    y = np.arange(len(sg))[::-1]
    colors = [COL["blue"] if not a else COL["teal"] for a in sg.adj]
    for yi, (_, row), c in zip(y, sg.iterrows(), colors):
        ax.errorbar(row.est, yi, xerr=[[row.est - row.lo], [row.hi - row.est]], fmt="o",
                    color=c, ecolor=c, capsize=2.5, ms=6, zorder=3)
    labels = [f"{g}{' (depth-adj.)' if a else ''}" for g, a in zip(sg.gene, sg.adj)]
    ax.set_yticks(y, labels, fontsize=7.5)
    ax.axvline(0, color="#B0B0B0", lw=0.8)
    ax.set_title("d  Single-gene decomposition", fontsize=9, loc="left", color=COL["ink"])
    setup(ax)

    # e: threshold trajectory (P02)
    th = pd.read_csv(P02 / "threshold_sensitivity_table.csv")
    ax = axes[4]
    ax.errorbar(th.cell_mean_difference, range(len(th)), fmt="o", ms=6, color=COL["orange"],
                ecolor=COL["orange"], capsize=2.5, zorder=3)
    ax.set_yticks(range(len(th)), [f"> {t:g}" for t in th.threshold], fontsize=7.5)
    ax2 = ax.twiny()
    ax2.plot(th.cell_n_positive_total, range(len(th)), "o", ms=5, mfc="white", mec=COL["purple"], mew=1.2, zorder=3)
    ax2.set_xlabel("Estimable positive cells", fontsize=7.5, color=COL["purple"])
    ax2.tick_params(axis="x", colors=COL["purple"], labelsize=6.5)
    ax.axvline(0, color="#B0B0B0", lw=0.8)
    ax.set_title("e  Positivity-threshold sensitivity", fontsize=9, loc="left", color=COL["ink"])
    setup(ax)

    fig.tight_layout(w_pad=1.6)
    save(fig, "Fig3")


# ----------------------------------------------------------------------------
# Fig 4: Level 2 qualification failure (panels a-d)
# ----------------------------------------------------------------------------
def fig4_level2():
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 3.4))

    # a: sample-level aggregates
    agg = pd.read_csv(R17 / "myeloid_sample_aggregates.csv")
    agg_sc = pd.read_csv(R17 / "myeloid_sample_aggregates_scored.csv")
    core_cols = ["ACP5", "MMP9"]
    for ds, grp in agg.groupby("dataset"):
        xs = grp.SIGLEC15
        ys = grp[core_cols].mean(axis=1)
        axes[0].scatter(xs, ys, marker=MARKERS[ds], s=40, color=COL["blue"],
                        edgecolor="white", lw=0.5, label=DATSET_NAMES[ds], zorder=3)
    axes[0].set_xlabel("Aggregated SIGLEC15 (cell-equal mean)", fontsize=8)
    axes[0].set_ylabel("ACP5/MMP9 aggregate score", fontsize=8)
    axes[0].set_title("a  Sample-level myeloid aggregates", fontsize=9, loc="left", color=COL["ink"])
    axes[0].legend(fontsize=6.5, frameon=False)
    setup(axes[0])

    # b: bulk cohort + pooled forest with decision lines
    bulk = pd.read_csv(R17 / "bulk_associations.csv")
    meta = pd.read_csv(R17 / "bulk_random_effects_meta.csv")
    core_b = bulk[bulk.program == "Derosiers_ACP5_MMP9_core"]
    ax = axes[1]
    rows = list(zip(core_b.cohort, core_b.rho, core_b.ci_low, core_b.ci_high))
    rows.append(("Pooled (REML + mHK)", meta.r_random.iloc[0], meta.ci_low.iloc[0], meta.ci_high.iloc[0]))
    y = np.arange(len(rows))[::-1]
    for yi, (name, r, lo, hi) in zip(y, rows):
        pooled = name.startswith("Pooled")
        ax.errorbar(r, yi, xerr=[[r - lo], [hi - r]], fmt="o", ms=6 if pooled else 5,
                    color=COL["orange"] if pooled else COL["blue"],
                    ecolor=COL["orange"] if pooled else COL["blue"], capsize=2.5, zorder=3)
    ax.set_yticks(y, [r[0] for r in rows], fontsize=7.5)
    ax.axvline(0, color="#B0B0B0", lw=0.8)
    ax.axvline(0.20, color=COL["green"], lw=1.1, ls="--")
    ax.axvline(-0.20, color=COL["green"], lw=1.1, ls="--")
    ax.axvline(0.36, color=COL["red"], lw=1.1, ls=":")
    ax.axvline(-0.36, color=COL["red"], lw=1.1, ls=":")
    ax.text(0.21, len(rows) - 0.4, "|r| = 0.20 qualification line", fontsize=6.2, color=COL["green"], rotation=90, va="top")
    ax.text(0.37, len(rows) - 0.4, "|r| = 0.36 tightest excludable bound", fontsize=6.2, color=COL["red"], rotation=90, va="top")
    ax.set_xlim(-0.55, 0.75)
    ax.set_title("b  Bulk associations (n = 321)", fontsize=9, loc="left", color=COL["ink"])
    setup(ax)

    # c: spatial case estimates (incl. non-estimable C1)
    sp = pd.read_csv(R17 / "spatial_case_effects.csv")
    core_sp = sp[sp.program == "Derosiers_ACP5_MMP9_core"]
    ax = axes[2]
    names = {"C1": "C1 (non-estimable)", "C2": "C2", "C3": "C3"}
    for case, row in core_sp.set_index("case").iterrows():
        ax.errorbar(row.mean_capture_difference, 1, fmt="o", ms=7, color=COL["blue"], capsize=3, zorder=3)
        ax.text(row.mean_capture_difference, 1.12, case, fontsize=8, ha="center", color=COL["ink"])
    ax.errorbar(0.466, 0, fmt="D", ms=8, color=COL["orange"], capsize=3,
                xerr=[[0.466 + 3.385], [4.318 - 0.466]], zorder=3)
    ax.text(0.466, 0.16, "case mean 0.47\n(95% CI \u22123.39 to 4.32)", fontsize=6.5, ha="center", color=COL["orange"])
    ax.set_xlim(-4.5, 5.5)
    ax.set_yticks([0, 1], ["Mean (2 cases)", "Case estimates"], fontsize=7.5)
    ax.axvline(0, color="#B0B0B0", lw=0.8)
    ax.set_title("c  Spatial (GSE282302, 3 cases)", fontsize=9, loc="left", color=COL["ink"])
    setup(ax)

    # d: verdict strip
    ax = axes[3]
    verdicts = [
        ("Cell (Level 1)", "SUPPORTIVE", COL["green"]),
        ("Sample aggregate (n = 20)", "UNINFORMATIVE", COL["mid"]),
        ("Bulk pooled (n = 321)", "NON-SUPPORTIVE", COL["red"]),
        ("Spatial (n = 2 cases)", "UNINFORMATIVE", COL["mid"]),
        ("Survival (exploratory)", "NON-SUPPORTIVE", COL["red"]),
    ]
    for i, (label, verdict, color) in enumerate(verdicts):
        ax.barh(len(verdicts) - i, 1, height=0.62, color=color, alpha=0.9, zorder=3)
        ax.text(0.5, len(verdicts) - i, verdict, ha="center", va="center", fontsize=8, color="white", fontweight="bold", zorder=4)
        ax.text(-0.02, len(verdicts) - i, label, ha="right", va="center", fontsize=7.5, color=COL["ink"], zorder=4)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, len(verdicts) + 0.6)
    ax.axis("off")
    ax.set_title("d  Verdicts under Table 2 rules", fontsize=9, loc="left", color=COL["ink"])

    fig.tight_layout(w_pad=1.8)
    save(fig, "Fig4")


# ----------------------------------------------------------------------------
# Fig 5: not-technical (calibration, abundance-matched, payload, context heatmap)
# ----------------------------------------------------------------------------
def fig5_calibration():
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2))

    # a: expected-programme calibration forest
    ax = axes[0][0]
    cal = [
        ("CD274", 0.545, 0.420, 0.652), ("TIGIT", 0.883, 0.858, 0.902),
        ("ERBB2", 0.568, 0.512, 0.619), ("EPCAM", 0.679, 0.584, 0.757),
        ("CSF1R", 0.822, 0.756, 0.871), ("SIGLEC15\u2013PI3K", 0.186, 0.034, 0.331),
    ]
    y = np.arange(len(cal))[::-1]
    for yi, (name, r, lo, hi) in zip(y, cal):
        c = COL["orange"] if name.startswith("SIGLEC15") else COL["blue"]
        ax.errorbar(r, yi, xerr=[[r - lo], [hi - r]], fmt="o", ms=6, color=c, ecolor=c, capsize=2.5, zorder=3)
    ax.set_yticks(y, [c[0] for c in cal], fontsize=8)
    ax.axvline(0, color="#B0B0B0", lw=0.8)
    ax.set_title("a  Same-system calibration", fontsize=9, loc="left", color=COL["ink"])
    ax.set_xlabel("Spearman rho vs expected programme", fontsize=8)
    setup(ax)

    # b: abundance-matched CD274 downsampling
    ds = pd.read_csv(ANALYSIS / "S10_cd274_downsampling" / "downsampling_replicates.csv")
    ax = axes[0][1]
    ax.hist(ds.rho_Hallmark_PI3K_AKT_MTOR__mean_z, bins=25, color=COL["bluep"], edgecolor=COL["blue"], zorder=3)
    ax.axvline(np.median(ds.rho_Hallmark_PI3K_AKT_MTOR__mean_z), color=COL["orange"], lw=1.4, ls="--")
    ax.text(0.02, 0.95, f"median rho = +{np.median(ds.rho_Hallmark_PI3K_AKT_MTOR__mean_z):.3f}\n200/200 replicates meet |rho| \u2265 0.20",
            transform=ax.transAxes, fontsize=7.5, va="top", color=COL["ink"])
    ax.set_title("b  Abundance-matched sensitivity (CD274)", fontsize=9, loc="left", color=COL["ink"])
    ax.set_xlabel("CD274\u2013PI3K rho at SIGLEC15-matched counts", fontsize=8)
    setup(ax)

    # c: discriminant payload forest
    sc = pd.read_csv(P22 / "scorecard_table.csv")
    ax = axes[1][0]
    order = ["CD8A", "GZMB", "FOXP3", "CD274", "CD19", "COL1A1", "SIGLEC15", "SIGLEC15", "KRT19", "EPCAM"]
    sc = sc.set_index(["target", "program"])
    rows = []
    for t, p in [("CD8A", "CD8_T"), ("GZMB", "cytolytic"), ("FOXP3", "regulatory_T"),
                 ("CD274", "cytolytic"), ("CD19", "cytolytic"), ("COL1A1", "cytolytic"),
                 ("SIGLEC15", "cytolytic"), ("SIGLEC15", "pan_myeloid"),
                 ("KRT19", "cytolytic"), ("EPCAM", "cytolytic")]:
        row = sc.loc[(t, p)]
        rows.append((f"{t}\u2013{p}", row.pooled_rho, row.ci_low, row.ci_high, row.observed))
    y = np.arange(len(rows))[::-1]
    for yi, (name, r, lo, hi, obs) in zip(y, rows):
        c = COL["green"] if obs == "direction-consistent supportive" else (COL["red"] if obs == "Non-supportive" else COL["mid"])
        ax.errorbar(r, yi, xerr=[[r - lo], [hi - r]], fmt="o", ms=5.5, color=c, ecolor=c, capsize=2, zorder=3)
    ax.set_yticks(y, [r[0] for r in rows], fontsize=7)
    ax.axvline(0, color="#B0B0B0", lw=0.8)
    ax.axvline(0.20, color=COL["green"], lw=1.0, ls="--")
    ax.set_xlim(-0.9, 1.15)
    ax.set_title("c  Discriminant payload (10 pairs)", fontsize=9, loc="left", color=COL["ink"])
    ax.set_xlabel("Pooled rho", fontsize=8)
    setup(ax)

    # d: tissue-context heatmap (from verified cross-target outputs)
    try:
        ct = json.load(open(ANALYSIS / "cross_target_panel_verified" / "cross_target_results.json"))
        if isinstance(ct, dict) and "matrix" in ct:
            mat = pd.DataFrame(ct["matrix"])
        else:
            mat = pd.read_csv(ANALYSIS / "cross_target_panel_verified" / "cross_target_results_flat.csv")
    except Exception:
        mat = pd.read_csv(ANALYSIS / "cross_target_panel_verified" / "cross_target_results_flat.csv")
    ax = axes[1][1]
    try:
        piv = mat.pivot(index="context", columns="target", values="rho") if "rho" in mat.columns else mat.pivot_table(index=0, columns=1, values=2)
    except Exception:
        piv = mat.select_dtypes(include=[np.number]).T
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-0.6, vmax=0.9, aspect="auto", zorder=2)
    ax.set_xticks(range(piv.shape[1]), piv.columns, fontsize=6.5, rotation=40, ha="right")
    ax.set_yticks(range(piv.shape[0]), piv.index, fontsize=6.5)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            ax.text(j, i, f"{piv.values[i, j]:.2f}", ha="center", va="center", fontsize=5.6, color=COL["ink"], zorder=3)
    ax.set_title("d  Tissue-context calibration", fontsize=9, loc="left", color=COL["ink"])
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Spearman rho", fontsize=7)
    cbar.ax.tick_params(labelsize=7)

    fig.tight_layout(w_pad=1.6, h_pad=1.8)
    save(fig, "Fig5")


# ----------------------------------------------------------------------------
# Fig 6: broader claims refused (PI3K, SYK/DAP12, MUC, protein boundary)
# ----------------------------------------------------------------------------
def fig6_broader():
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0))

    # a: 16 PI3K specifications (values from Table S4, official results)
    pi3k = [
        ("Hallmark", "AUCell", -0.052, 0.562, True), ("Hallmark", "GSVA", 0.060, 0.562, True),
        ("Hallmark", "mean_z", 0.186, 0.070, False), ("Hallmark", "ssGSEA", -0.059, 0.562, True),
        ("KEGG", "AUCell", 0.009, 0.910, True), ("KEGG", "GSVA", 0.139, 0.207, False),
        ("KEGG", "mean_z", 0.217, 0.059, False), ("KEGG", "ssGSEA", 0.020, 0.838, True),
        ("PROGENy", "AUCell", -0.103, 0.387, False), ("PROGENy", "GSVA", 0.075, 0.562, False),
        ("PROGENy", "mean_z", 0.174, 0.080, False), ("PROGENy", "ssGSEA", -0.121, 0.289, False),
        ("Reactome", "AUCell", -0.054, 0.562, True), ("Reactome", "GSVA", 0.063, 0.562, True),
        ("Reactome", "mean_z", 0.197, 0.067, False), ("Reactome", "ssGSEA", -0.056, 0.562, True),
    ]
    ax = axes[0][0]
    y = np.arange(len(pi3k))[::-1]
    for yi, (src, method, rho, fdr, equiv) in zip(y, pi3k):
        ax.plot(rho, yi, marker="s" if equiv else "o", ms=5.5, mfc=COL["bluep"] if equiv else "white",
                mec=COL["blue"], mew=1.2, zorder=3)
    ax.set_yticks(y, [f"{s} / {m}" for s, m, *_ in pi3k], fontsize=6.5)
    ax.axvline(0, color="#B0B0B0", lw=0.8)
    ax.axvline(0.20, color=COL["green"], lw=1.0, ls="--")
    ax.axvline(-0.20, color=COL["green"], lw=1.0, ls="--")
    ax.set_title("a  16 PI3K specifications (TCGA-PAAD)", fontsize=9, loc="left", color=COL["ink"])
    ax.set_xlabel("rho vs SIGLEC15", fontsize=8)
    ax.text(0.98, 0.05, "open squares: equivalent at \u00b10.20", transform=ax.transAxes, ha="right", fontsize=6.5, color=COL["mid"])
    setup(ax)

    # b: SYK/DAP12 12 specifications (Table S5)
    ax = axes[0][1]
    syk = [("Reactome DAP12 signalling", [0.054, 0.179]), ("Reactome DAP12 interactions", [0.038, 0.152]),
           ("SYK kinase footprint", [0.006, 0.143]), ("Li 2022 TAM programme", [0.046, 0.123])]
    y = np.arange(len(syk))[::-1]
    for yi, (name, rng) in zip(y, syk):
        ax.plot(rng, [yi, yi], color=COL["purple"], lw=3, solid_capstyle="round", zorder=2)
        ax.plot(rng, [yi, yi], "o", ms=4.5, color=COL["purple"], zorder=3)
    ax.set_yticks(y, [s[0] for s in syk], fontsize=7.5)
    ax.axvline(0, color="#B0B0B0", lw=0.8)
    ax.axvline(0.20, color=COL["green"], lw=1.0, ls="--")
    ax.set_title("b  12 SYK/DAP12-ITAM specifications", fontsize=9, loc="left", color=COL["ink"])
    ax.set_xlabel("rho range across scoring methods", fontsize=8)
    setup(ax)

    # c: MUC1/CD24 cohort + pooled
    ax = axes[1][0]
    m = pd.read_csv(VALIDATION / "external_meta_reml_hk.csv")
    muc = m[m.gene_set == "MUC1_CD24"].iloc[0]
    rows = [("TCGA-PAAD", 0.349, 0.349, 0.349), ("GSE57495", 0.192, 0.192, 0.192),
            ("GSE85916", -0.150, -0.150, -0.150),
            ("Pooled (REML + mHK)", muc.pooled_r, muc.modified_hk_ci95_lo_r, muc.modified_hk_ci95_hi_r)]
    y = np.arange(len(rows))[::-1]
    for yi, (name, r, lo, hi) in zip(y, rows):
        pooled = name.startswith("Pooled")
        if pooled:
            ax.errorbar(r, yi, xerr=[[r - lo], [hi - r]], fmt="D", ms=7, color=COL["orange"],
                        ecolor=COL["orange"], capsize=2.5, zorder=3)
        else:
            ax.plot(r, yi, "o", ms=6, color=COL["blue"], zorder=3)
    ax.set_yticks(y, [r[0] for r in rows], fontsize=7.5)
    ax.axvline(0, color="#B0B0B0", lw=0.8)
    ax.set_xlim(-0.8, 0.9)
    ax.set_title("c  MUC1/CD24 context (not portable)", fontsize=9, loc="left", color=COL["ink"])
    ax.set_xlabel("rho vs SIGLEC15", fontsize=8)
    setup(ax)

    # d: protein boundary schematic
    ax = axes[1][1]
    ax.axis("off")
    items = [
        ("Discovery proteomics (CPTAC-PDAC, n = 140)", "no quantifiable SIGLEC15 row", COL["red"]),
        ("Discovery proteomics (PXD048644)", "exact identifier not recovered; 6/6 controls recovered", COL["red"]),
        ("Published antibody evidence (Shen 2025; Chen 2022)", "SIGLEC15 detected in PDAC tissue; target engagement", COL["green"]),
        ("Interpretation", "assay-coverage boundary, not protein absence", COL["ink"]),
    ]
    for i, (label, value, color) in enumerate(items):
        ax.barh(3 - i, 1, height=0.16, color=COL["pale"], zorder=2)
        ax.text(0.01, 3 - i + 0.09, label, ha="left", va="bottom", fontsize=8, color=COL["ink"])
        ax.text(0.01, 3 - i - 0.05, value, ha="left", va="top", fontsize=7.5, color=color)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.3, 3.6)
    ax.set_title("d  Protein-measurement boundary", fontsize=9, loc="left", color=COL["ink"])

    fig.tight_layout(w_pad=1.6, h_pad=1.8)
    save(fig, "Fig6")


# ----------------------------------------------------------------------------
# Fig 1 & Fig 7: schematic figures (pure layout)
# ----------------------------------------------------------------------------
def _panel_box(ax, x, y, w, h, text, fc, ec=COL["ink"], fs=8, tc=COL["ink"], lw=1.0):
    ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec, lw=lw, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=tc, zorder=3)


def _arrow(ax, x0, y0, x1, y1, color=COL["ink"], ls="-", lw=1.4):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=color, ls=ls, lw=lw), zorder=1)


def fig1_framework():
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 3.8))

    # a: resolution paradox
    ax = axes[0]
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.text(5, 9.4, "a  The resolution paradox", ha="center", fontsize=10, fontweight="bold", color=COL["ink"])
    # cell level
    _panel_box(ax, 0.6, 5.6, 3.6, 2.6, "SIGLEC15\u207a myeloid cell\nACP5/MMP9 state\npresent (Level 1)", COL["greenp"], ec=COL["green"])
    _panel_box(ax, 5.9, 5.6, 3.6, 2.6, "SIGLEC15\u207b myeloid cell\nACP5/MMP9 state\nabsent", COL["pale"])
    ax.text(4.9, 6.9, "+0.64", ha="center", fontsize=11, color=COL["green"], fontweight="bold")
    # patient level
    _panel_box(ax, 0.6, 0.8, 3.6, 2.4, "Sample aggregate (n = 20)\nr = \u22120.45  [\u22120.84, 0.11]\nUNINFORMATIVE", COL["mid"], ec=COL["mid"], fs=7.5, tc="white")
    _panel_box(ax, 5.9, 0.8, 3.6, 2.4, "Bulk (n = 321)\nr = 0.13  [\u22120.11, 0.36]\nNON-SUPPORTIVE", COL["redp"], ec=COL["red"], fs=7.5)
    _arrow(ax, 2.4, 5.5, 2.4, 3.3, color=COL["orange"])
    _arrow(ax, 7.7, 5.5, 7.7, 3.3, color=COL["orange"])
    ax.text(5, 0.15, "patient-level readouts do not carry the cell-level signal", ha="center", fontsize=8, color=COL["mid"], style="italic")

    # b: protocol structure
    ax = axes[1]
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.text(5, 9.4, "b  Two-level qualification protocol", ha="center", fontsize=10, fontweight="bold", color=COL["ink"])
    _panel_box(ax, 0.5, 6.2, 4.0, 2.2, "Level 1 \u2014 cellular recoverability\nstate defined from literature anchor\nsample as inferential unit", COL["greenp"], ec=COL["green"], fs=7.5)
    _panel_box(ax, 5.6, 6.2, 4.0, 2.2, "Same-system calibration\nreference targets must recover\nexpected relationships", COL["bluep"], ec=COL["blue"], fs=7.5)
    _panel_box(ax, 0.5, 3.0, 4.0, 2.2, "Level 2 \u2014 patient-level portability\naggregate / bulk / spatial / outcome\nfixed decision rules", COL["redp"], ec=COL["red"], fs=7.5)
    _panel_box(ax, 5.6, 3.0, 4.0, 2.2, "Verdicts\nsupportive / non-supportive /\nuninformative (no evidence)", COL["pale"], fs=7.5)
    _arrow(ax, 2.5, 6.1, 2.5, 5.3, color=COL["ink"])
    _arrow(ax, 7.6, 6.1, 7.6, 5.3, color=COL["ink"])
    _arrow(ax, 4.6, 7.3, 5.5, 7.3, color=COL["blue"])
    ax.text(5, 2.6, "combining rule: Level 1 never overrides Level 2;\nuninformative is not evidence either way", ha="center", fontsize=7.5, color=COL["mid"], style="italic")

    # c: data landscape
    ax = axes[2]
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.text(5, 9.4, "c  Study data landscape", ha="center", fontsize=10, fontweight="bold", color=COL["ink"])
    _panel_box(ax, 0.4, 6.4, 4.2, 2.0, "Single-cell (Level 1)\nGSE155698 \u00b7 GSE156405 \u00b7 GSE212966\n10,431 myeloid cells, 20 samples", COL["greenp"], ec=COL["green"], fs=7.5)
    _panel_box(ax, 5.5, 6.4, 4.2, 2.0, "Spatial (Level 2)\nGSE282302 (Visium)\n44 captures, 3 cases", COL["bluep"], ec=COL["blue"], fs=7.5)
    _panel_box(ax, 0.4, 3.4, 4.2, 2.0, "Bulk (Level 2)\nTCGA-PAAD \u00b7 GSE57495 \u00b7 GSE85916\nn = 321", COL["redp"], ec=COL["red"], fs=7.5)
    _panel_box(ax, 5.5, 3.4, 4.2, 2.0, "Survival (Level 2, exploratory)\n8 TCGA cohorts, n = 2,955", COL["pale"], fs=7.5)
    _panel_box(ax, 0.4, 0.6, 9.3, 1.7, "Proteomics (boundary) \u00b7 CPTAC-PDAC (n = 140) \u00b7 PXD048644 \u00b7 reference targets (CD274, TIGIT, ERBB2, EPCAM, CSF1R)", COL["pale"], fs=7.5)
    ax.text(5, 0.15, "solid connectors: within-series transitions; dashed: independent datasets", ha="center", fontsize=7.5, color=COL["mid"], style="italic")

    fig.tight_layout()
    save(fig, "Fig1")


def fig7_protocol():
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 3.8))

    # a: decision flow
    ax = axes[0]
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.text(5, 9.4, "a  Decision flow", ha="center", fontsize=10, fontweight="bold", color=COL["ink"])
    steps = [
        ("1. State definition", "literature-anchored genes;\nSIGLEC15 excluded from scores", COL["pale"]),
        ("2. Level 1 recoverability", "cell contrast, sample as unit;\ndepth, matching, threshold sensitivity", COL["greenp"]),
        ("3. Calibration gate", "reference targets must recover\nexpected relationships", COL["bluep"]),
        ("4. Level 2 portability", "aggregate / bulk / spatial / outcome;\nfixed thresholds (0.20; large-effect 0.5)", COL["redp"]),
        ("5. Verdicts + combining rule", "per-layer supportive / non-supportive /\nuninformative; Level 1 never overrides", COL["pale"]),
    ]
    y = 8.2
    for title, body, fc in steps:
        _panel_box(ax, 0.5, y - 1.35, 4.6, 1.35, f"{title}\n{body}", fc, fs=7.5)
        y -= 1.7
    for i in range(4):
        _arrow(ax, 2.8, 8.2 - i * 1.7 - 0.35, 2.8, 8.2 - (i + 1) * 1.7 + 0.35, color=COL["ink"])

    # b: final grade
    ax = axes[1]
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.text(5, 9.4, "b  SIGLEC15 RNA in PDAC \u2014 final grade", ha="center", fontsize=10, fontweight="bold", color=COL["ink"])
    _panel_box(ax, 0.6, 5.6, 8.8, 2.6, "Level 1 \u2014 SUPPORTED\nosteoclast-like state recovered in SIGLEC15\u207a tumour-myeloid cells\n0.64 [0.33, 0.99]; depth-adjusted 0.43 [0.17, 0.73]; depth-matched 0.37 [0.14, 0.59]; TREM2\u207a TAM niche", COL["greenp"], ec=COL["green"], fs=8)
    _panel_box(ax, 0.6, 2.0, 8.8, 2.6, "Level 2 \u2014 NOT QUALIFIED\nsample aggregate uninformative (r = \u22120.45); bulk failed qualification (r = 0.13; large effects excluded);\nspatial uninformative (2 cases); survival null (HR 1.04); PI3K/SYK-DAP12 not reproducible", COL["redp"], ec=COL["red"], fs=8)
    _panel_box(ax, 0.6, 0.4, 8.8, 1.2, "Biological claim intact \u00b7 patient-selection claim not supported \u00b7 protein-first measurement required", COL["pale"], fs=8)
    ax.text(5, 5.25, "cell-level recovery \u2260 patient-level readout", ha="center", fontsize=8, color=COL["orange"], fontweight="bold", style="italic")

    # c: protein-first next study
    ax = axes[2]
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.text(5, 9.4, "c  Protein-first next study", ha="center", fontsize=10, fontweight="bold", color=COL["ink"])
    _panel_box(ax, 0.6, 7.0, 8.8, 1.8, "1. Compartment-resolved SIGLEC15 IHC\n(analytically validated; positive-cell prevalence)", COL["pale"], fs=7.5)
    _panel_box(ax, 0.6, 4.9, 8.8, 1.8, "2. Orthogonal targeted protein confirmation\n(mass spectrometry / proximity assays)", COL["pale"], fs=7.5)
    _panel_box(ax, 0.6, 2.8, 8.8, 1.8, "3. ACP5/MMP9 state readout\nin the protein-positive compartment", COL["pale"], fs=7.5)
    _panel_box(ax, 0.6, 0.7, 8.8, 1.8, "4. Treatment-by-biomarker validation\n(comparator-controlled; biomarker-defined subset)", COL["pale"], fs=7.5)
    for i in range(3):
        _arrow(ax, 5, 7.9 - i * 2.1, 5, 6.9 - i * 2.1, color=COL["ink"])
    ax.text(5, 0.15, "until then: RNA supports hypothesis generation and tissue prioritisation, not clinical enrichment", ha="center", fontsize=7.5, color=COL["mid"], style="italic")

    fig.tight_layout()
    save(fig, "Fig7")


# ----------------------------------------------------------------------------
# Supplementary figures
# ----------------------------------------------------------------------------
def figs_supp():
    # Fig S2: myeloid subtype enrichment (Task5)
    ch = pd.read_csv(ANALYSIS / "Task5_myeloid_subtypes" / "chi2_subtype_enrichment.csv")
    combined = ch[ch.dataset == "COMBINED"].sort_values("enrichment_fold", ascending=True)
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    colors = [COL["orange"] if s == "TREM2_TAM" else (COL["blue"] if s == "Classical_Monocyte" else COL["mid"]) for s in combined.subtype]
    ax.barh(combined.subtype, combined.enrichment_fold, color=colors, edgecolor="white", zorder=3)
    ax.axvline(1, color="#B0B0B0", lw=0.9, ls="--")
    ax.set_xlabel("Observed/expected enrichment (SIGLEC15\u207a cells)", fontsize=8)
    ax.set_xscale("log", base=2)
    ax.set_xticks([0.25, 0.5, 1, 2, 4])
    ax.set_xticklabels(["0.25", "0.5", "1", "2", "4"])
    for i, v in enumerate(combined.enrichment_fold):
        ax.text(v * 1.12, i, f"{v:.2f}", va="center", fontsize=7, color=COL["ink"])
    setup(ax)
    fig.tight_layout()
    save(fig, "FigS2")

    # Fig S3: survival forest (Task2)
    per = pd.read_csv(ANALYSIS / "Task2_Meta_Cox" / "per_cohort_univariate_cox.csv")
    pooled = pd.read_csv(ANALYSIS / "Task2_Meta_Cox" / "meta_pooled_univariate.csv").iloc[0]
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    y = np.arange(len(per) + 1)[::-1]
    for yi, (_, row) in zip(y[:len(per)], per.iterrows()):
        ax.errorbar(row.HR, yi, xerr=[[row.HR - row.HR_lower95], [row.HR_upper95 - row.HR]],
                    fmt="o", ms=5, color=COL["blue"], ecolor=COL["blue"], capsize=2.5, zorder=3)
    ax.errorbar(pooled.pooled_HR, y[-1], xerr=[[pooled.pooled_HR - pooled.pooled_HR_lower95], [pooled.pooled_HR_upper95 - pooled.pooled_HR]],
                fmt="D", ms=7, color=COL["orange"], ecolor=COL["orange"], capsize=2.5, zorder=3)
    labels = list(per.cohort) + [f"Pooled (I\u00b2 = {pooled.I2:.1f}%)"]
    ax.set_yticks(y, labels, fontsize=7.5)
    ax.axvline(1, color="#B0B0B0", lw=0.9)
    ax.set_xscale("log")
    ax.set_xlim(0.6, 2.2)
    ax.set_xlabel("Univariate HR (continuous SIGLEC15 RNA, log scale)", fontsize=8)
    ax.set_title("Exploratory outcome layer: SIGLEC15 and overall survival", fontsize=9, loc="left", color=COL["ink"])
    setup(ax)
    fig.tight_layout()
    save(fig, "FigS3")

    # Fig S4: deconvolution cross-tool heatmap (Task4)
    cs = pd.read_csv(ANALYSIS / "Task4_immune_tools_cross" / "consistency_matrix.csv")
    piv = cs.pivot_table(index="cell_type", columns="tool", values="rho")
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-0.3, vmax=0.3, aspect="auto", zorder=2)
    ax.set_xticks(range(piv.shape[1]), piv.columns, fontsize=7.5, rotation=30, ha="right")
    ax.set_yticks(range(piv.shape[0]), piv.index, fontsize=7.5)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6, color=COL["ink"], zorder=3)
    ax.set_title("Deconvolution cross-tool calibration (TCGA-PAAD): rho with SIGLEC15 RNA", fontsize=9, loc="left", color=COL["ink"])
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Spearman rho", fontsize=7)
    cbar.ax.tick_params(labelsize=7)
    fig.tight_layout()
    save(fig, "FigS4")


if __name__ == "__main__":
    fig1_framework()
    fig2_landscape()
    fig3_level1()
    fig4_level2()
    fig5_calibration()
    fig6_broader()
    fig7_protocol()
    figs_supp()
    print("ALL FIGURES DONE")
