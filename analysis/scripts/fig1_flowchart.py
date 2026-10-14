#!/usr/bin/env python3
"""
Fig. 1 (redrawn) - Two-level qualification protocol applied to SIGLEC15 in PDAC.
Vector PDF, Springer Nature full-page-width compatible (170 mm).
"""
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42      # embed TrueType (editable in AI/Inkscape)
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

# ---------------------------------------------------------------- palette
C = dict(
    ink="#1B2227",
    sub="#4A555D",
    neutral_f="#F4F6F7", neutral_e="#A9B4BB",
    anchor_f="#E9EFF7", anchor_e="#4E70A0",
    l1_f="#E6F1EA",     l1_e="#3C7F5E",
    gate_f="#FBF2E0",   gate_e="#B98A2C",
    l2_f="#F5F1F5",     l2_e="#7C6A87",
    out_f="#EFF1F3",    out_e="#5A666E",
    sup="#2E7D5B", sup_f="#DCEDE4",
    non="#A63A38", non_f="#F5E1DF",
    uni="#5E6A72", uni_f="#E4E8EA",
    solid="#C4632A", dashed="#3D6EA8",
)
FS_H   = 8.2   # box header
FS_B   = 6.7   # body
FS_S   = 6.0   # small / caption-ish
FS_V   = 6.4   # verdict pill

# ---------------------------------------------------------------- canvas
MM = 1 / 25.4
W_IN, H_IN = 170 * MM, 204 * MM    # JTM full width; near-isotropic layout and page-safe height
fig = plt.figure(figsize=(W_IN, H_IN))
fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100)
ax.set_ylim(0, 121)
ax.axis("off")

# ---------------------------------------------------------------- helpers
def box(x0, y0, x1, y1, fc, ec, lw=1.0, r=1.4, z=2, ls="-"):
    p = FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                       boxstyle=f"round,pad=0,rounding_size={r}",
                       facecolor=fc, edgecolor=ec, linewidth=lw,
                       linestyle=ls, zorder=z, mutation_aspect=1)
    ax.add_patch(p)
    return p

def txt(x, y, s, size=FS_B, weight="normal", color=None, ha="center",
        va="center", style="normal", z=5):
    ax.text(x, y, s, fontsize=size, fontweight=weight,
            color=color or C["ink"], ha=ha, va=va, style=style, zorder=z)

def arrow(x0, y0, x1, y1, color=None, dashed=False, lw=1.15, z=1, rad=0.0):
    a = FancyArrowPatch((x0, y0), (x1, y1),
                        arrowstyle="-|>", mutation_scale=8.5,
                        linewidth=lw, color=color or C["sub"],
                        linestyle=(0, (3.5, 2.2)) if dashed else "-",
                        connectionstyle=f"arc3,rad={rad}",
                        shrinkA=0, shrinkB=0, zorder=z,
                        joinstyle="miter")
    ax.add_patch(a)

def pill(xc, y, label, kind, w=15.5, h=3.6):
    fc = {"sup": C["sup"], "non": C["non"], "uni": C["uni"]}[kind]
    box(xc - w / 2, y - h / 2, xc + w / 2, y + h / 2, fc, fc, lw=0.6, r=1.8, z=4)
    txt(xc, y - 0.05, label, size=FS_V, weight="bold", color="white", z=5)

def stage_tag(x, y, s):
    txt(x, y, s, size=FS_S, weight="bold", color=C["sub"], ha="left")

# ================================================================ 1. ANCHOR
box(24, 111.5, 76, 119.5, C["anchor_f"], C["anchor_e"], lw=1.1)
txt(50, 117.4, "STATE DEFINITION", size=FS_H, weight="bold", color=C["anchor_e"])
txt(50, 114.7, "Literature anchor (Derosiers et al.): ACP5 + MMP9", size=FS_B)
txt(50, 112.6, "osteoclast-like myeloid state  ·  SIGLEC15 excluded from every score",
    size=FS_S, color=C["sub"])

arrow(50, 111.5, 50, 107.5)

# ================================================================ 2. LEVEL 1
box(6, 87.0, 94, 107.5, C["l1_f"], C["l1_e"], lw=1.1)
txt(9.5, 104.8, "LEVEL 1", size=FS_H, weight="bold", color=C["l1_e"], ha="left")
ax.plot([20.0, 20.0], [103.6, 106.0], color=C["l1_e"], lw=0.8, zorder=4)
txt(22.5, 104.8, "Cellular recoverability — is the state present in target-bearing cells?",
    size=FS_H, color=C["ink"], ha="left")
txt(9.5, 101.9,
    "3 independent scRNA-seq series  ·  10,431 tumour-myeloid cells  ·  20 samples (19 estimable)  ·  sample = inferential unit",
    size=FS_S, color=C["sub"], ha="left")

# result chips
chip_y0, chip_y1 = 91.35, 99.8
chips = [
    (8.0, 26.5, "Unadjusted contrast", "+0.64", "[0.33, 0.99]", r"$p = 7.9 \times 10^{-4}$"),
    (28.0, 46.5, "Depth-adjusted", "+0.43", "[0.17, 0.73]", "p = 0.007"),
    (48.0, 66.5, "Depth-matched (1:1)", "+0.37", "[0.14, 0.59]", "p = 0.011"),
    (68.0, 92.0, "Niche placement", "TREM2+ TAM", "enrichment 2.45", "monocytes 0.70"),
]
for x0, x1, head, val, ci, p in chips:
    box(x0, chip_y0, x1, chip_y1, "#FFFFFF", C["l1_e"], lw=0.7, r=1.0, z=3)
    xc = (x0 + x1) / 2
    txt(xc, 98.2, head, size=FS_S, weight="bold", color=C["l1_e"])
    txt(xc, 95.6, val, size=8.6, weight="bold")
    txt(xc, 93.6, ci, size=FS_S, color=C["sub"])
    txt(xc, 92.25, p, size=5.5, color=C["sub"])

pill(50, 89.25, "LEVEL 1  VERDICT:  SUPPORTIVE", "sup", w=52, h=3.0)

arrow(50, 87.0, 50, 83.0)

# ================================================================ 3. GATE
box(6, 68.0, 94, 83.0, C["gate_f"], C["gate_e"], lw=1.1)
txt(9.5, 80.4, "CALIBRATION GATE", size=FS_H, weight="bold", color=C["gate_e"], ha="left")
ax.plot([31.5, 31.5], [79.2, 81.6], color=C["gate_e"], lw=0.8, zorder=4)
txt(34.0, 80.4, "No level is scored until the analytical instrument is verified",
    size=FS_H, color=C["ink"], ha="left")

gate_items = [
    (8.0, 28.0, "Expected-programme\ncalibration",
     "5 reference targets\nrecovered  (r 0.55\u20130.88)"),
    (29.3, 49.3, "Abundance-matched\nsensitivity",
     "CD274\u2013PI3K retained in\n200/200 replicates"),
    (50.7, 70.7, "Discriminant\npayload",
     "8/8 control pairs classified\nas expected; SIGLEC15 silent"),
    (72.0, 92.0, "State-level positive\ncontrols",
     "TREM2 0.76 · SPP1 0.41\nLRRC15 0.71 · CSF1R 0.76"),
]
for x0, x1, head, body in gate_items:
    box(x0, 70.6, x1, 78.4, "#FFFFFF", C["gate_e"], lw=0.7, r=1.0, z=3)
    xc = (x0 + x1) / 2
    txt(xc, 76.9, head, size=FS_S, weight="bold", color=C["gate_e"])
    txt(xc, 73.2, body, size=5.7, color=C["ink"])

txt(50, 69.2, "GATE PASSED  \u2014  the analytical system returns positive verdicts when the biology is present",
    size=FS_S, weight="bold", color=C["gate_e"])

# ================================================================ 4. LEVEL 2 header
arrow(50, 68.0, 50, 64.5)
box(6, 58.0, 94, 64.5, C["l2_f"], C["l2_e"], lw=1.1)
txt(9.5, 62.6, "LEVEL 2", size=FS_H, weight="bold", color=C["l2_e"], ha="left")
ax.plot([20.0, 20.0], [61.4, 63.8], color=C["l2_e"], lw=0.8, zorder=4)
txt(22.5, 62.6, "Patient-level portability — does the same state qualify at patient level?",
    size=FS_H, color=C["ink"], ha="left")
txt(9.5, 59.7,
    "Sample is the unit; cells and spots contribute only to scoring.  Qualification |r| \u2265 0.20; large effect |r| \u2265 0.50.",
    size=FS_S, color=C["sub"], ha="left")

# ================================================================ 5. LEVEL 2 layers
lay_y0, lay_y1 = 33.0, 54.0
cols = [(6.0, 26.5), (28.5, 49.0), (51.0, 71.5), (73.5, 94.0)]
layers = [
    dict(head="Sample aggregate", src="20 samples\n(same three series)",
         val="r = \u22120.45", ci="[\u22120.84, 0.11]", note="permutation p = 0.063",
         verdict=("UNINFORMATIVE", "uni"), foot="n too small to adjudicate", dashed=False),
    dict(head="Independent bulk RNA", src="3 cohorts, n = 321\nTCGA-PAAD/57495/85916",
         val="pooled r = 0.13", ci="[\u22120.11, 0.36]", note="|r| \u2265 0.5 excluded",
         verdict=("NOT QUALIFIED", "non"), foot="0.20\u20130.36 remains compatible", dashed=True),
    dict(head="Independent spatial", src="2 of 3 cases estimable\n(post-neoadjuvant)",
         val="\u0394 = 0.47", ci="[\u22123.39, 4.32]", note="both cases positive",
         verdict=("UNINFORMATIVE", "uni"), foot="exploratory", dashed=True),
    dict(head="Survival outcome", src="8 TCGA cohorts\nn = 2,955",
         val="HR = 1.04", ci="[0.96, 1.12]", note="c-index 0.504 (PAAD)",
         verdict=("NON-SUPPORTIVE", "non"), foot="exploratory", dashed=True),
]

for (x0, x1), L in zip(cols, layers):
    xc = (x0 + x1) / 2
    box(x0, lay_y0, x1, lay_y1, C["neutral_f"], C["neutral_e"], lw=0.9, r=1.2)
    txt(xc, 52.0, L["head"], size=7.2, weight="bold")
    ax.plot([x0 + 2.0, x1 - 2.0], [50.5, 50.5], color=C["neutral_e"], lw=0.6, zorder=4)
    txt(xc, 48.2, L["src"], size=5.7, color=C["sub"])
    txt(xc, 44.4, L["val"], size=9.0, weight="bold")
    txt(xc, 42.2, L["ci"], size=6.2, color=C["sub"])
    txt(xc, 40.4, L["note"], size=5.5, color=C["sub"], style="italic")
    pill(xc, 37.2, L["verdict"][0], L["verdict"][1], w=x1 - x0 - 4.0, h=3.5)
    txt(xc, 34.5, L["foot"], size=5.4, color=C["sub"])
    # connector from the Level 2 header bar
    arrow(xc, 58.0, xc, lay_y1,
          color=C["dashed"] if L["dashed"] else C["solid"],
          dashed=L["dashed"], lw=1.05)

# ================================================================ 6. combining rule
for (x0, x1) in cols:
    xc = (x0 + x1) / 2
    arrow(xc, lay_y0, xc, 28.0, lw=0.95)

box(6, 21.5, 94, 28.0, C["out_f"], C["out_e"], lw=1.0)
txt(50, 25.9, "COMBINING RULE", size=FS_H, weight="bold", color=C["out_e"])
txt(50, 23.2,
    "Level 1 never overrides Level 2   ·   uninformative = evidence in neither direction   ·   "
    "excluded effects do not prove absence",
    size=5.6, color=C["ink"])

arrow(50, 21.5, 50, 18.0)

# ================================================================ 7. outcome
box(6, 8.2, 94, 18.0, "#FFFFFF", C["ink"], lw=1.3)
txt(50, 15.9, "Cell-level recovery  \u2260  patient-level portability",
    size=9.4, weight="bold")
txt(50, 13.1,
    "SIGLEC15 in PDAC: the state is recovered in cells; no patient-level RNA readout qualifies for patient selection.",
    size=FS_B)
txt(50, 10.8,
    "Next step \u2014 protein-first:  compartment-resolved IHC  \u2192  orthogonal protein assay",
    size=5.6, color=C["sub"])
txt(50, 9.35,
    "ACP5/MMP9 in protein-positive cells  \u2192  treatment \u00d7 biomarker",
    size=5.6, color=C["sub"])

# ================================================================ 8. legend
ax.plot([6, 94], [6.4, 6.4], color="#C9D0D5", lw=0.7)

for lx, lab, kind in [(6.0, "supportive", "sup"),
                      (26.0, "non-supportive / not qualified", "non"),
                      (60.0, "uninformative", "uni")]:
    box(lx, 3.4, lx + 3.2, 5.2, C[kind], C[kind], lw=0.5, r=0.8, z=4)
    txt(lx + 4.2, 4.3, lab, size=5.7, color=C["sub"], ha="left")

ax.plot([6.0, 11.0], [1.4, 1.4], color=C["solid"], lw=1.15, solid_capstyle="butt")
txt(11.8, 1.4, "solid connector: same scRNA-seq series", size=5.7, color=C["sub"], ha="left")
ax.plot([48.0, 53.0], [1.4, 1.4], color=C["dashed"], lw=1.15,
        linestyle=(0, (3.5, 2.2)))
txt(53.8, 1.4, "dashed connector: independent patient datasets",
    size=5.7, color=C["sub"], ha="left")

# ---------------------------------------------------------------- output
# Files are written next to this script. Change OUTDIR if you prefer elsewhere.
import os
OUTDIR = os.path.dirname(os.path.abspath(__file__))
pdf_path = os.path.join(OUTDIR, "Fig1_flowchart.pdf")
png_path = os.path.join(OUTDIR, "Fig1_flowchart_preview.png")
svg_path = os.path.join(OUTDIR, "Fig1_flowchart.svg")

# Guard against silent clipping after future wording or font changes.
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
canvas = fig.get_window_extent(renderer)
outside = []
for artist in ax.texts:
    extent = artist.get_window_extent(renderer)
    if (extent.x0 < canvas.x0 - 1 or extent.y0 < canvas.y0 - 1 or
            extent.x1 > canvas.x1 + 1 or extent.y1 > canvas.y1 + 1):
        outside.append(artist.get_text())
if outside:
    raise RuntimeError("Text outside figure canvas: " + " | ".join(outside))

fig.savefig(pdf_path, format="pdf", bbox_inches=None, pad_inches=0)  # vector submission file
fig.savefig(svg_path, format="svg", bbox_inches=None, pad_inches=0)  # editable vector source
fig.savefig(png_path, dpi=600, bbox_inches=None, pad_inches=0)       # visual QA
# fig.savefig(os.path.join(OUTDIR, "Fig1_flowchart.tif"), dpi=600)  # if TIFF required
print("written:", pdf_path)
print("written:", svg_path)
print("written:", png_path)
