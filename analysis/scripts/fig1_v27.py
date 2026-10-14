"""Fig 1 (v27): cross-scale construct-to-readout chain with coverage markers."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({"font.family": "sans-serif", "font.size": 8})

BLUE, BLUEP = "#2B5B8C", "#DCE9F5"
TEAL, TEALP = "#1F6F6B", "#DCEDEC"
ORANGE, ORANGEP = "#B5651D", "#F6E4D0"
GREEN, GREENP = "#2E6B34", "#DFEEDF"
RED, REDP = "#A63A2E", "#F5DEDB"
GREY, GREYP = "#4A4A4A", "#E8E8E8"

fig, ax = plt.subplots(figsize=(6.7, 5.3))
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis("off")


def box(x, y, w, h, title, body, fc, ec, fs_t=8.2, fs_b=6.6, yt=0.72, yb=0.35):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04",
                       fc=fc, ec=ec, lw=1.0, mutation_aspect=1)
    ax.add_patch(p)
    ax.text(x + w / 2, y + h * yt, title, ha="center", va="center",
            fontsize=fs_t, fontweight="bold", color="#111111")
    ax.text(x + w / 2, y + h * yb, body, ha="center", va="center",
            fontsize=fs_b, color="#222222")


def arrow(x, y0, y1):
    ax.add_patch(FancyArrowPatch((x, y0), (x, y1), arrowstyle="-|>",
                                 mutation_scale=13, color="#555555", lw=1.1))


BX, BW, BH = 0.55, 5.1, 0.86
ys = [9.0, 7.85, 6.7, 5.55, 4.4, 3.25, 2.1]

box(BX, ys[0], BW, BH, "CELLULAR CONSTRUCT",
    "SIGLEC15-associated ACP5/MMP9-centred relation\nwithin tumour-associated myeloid cells",
    BLUEP, BLUE)
box(BX, ys[1], BW, BH, "PATIENT-LEVEL MEASURAND",
    "abstract within-patient target\u2013state coupling\nin the tumour-myeloid compartment",
    BLUEP, BLUE)
box(BX, ys[2], BW, BH, "MEASUREMENT PROCEDURE",
    "scRNA operationalization: detectable-vs-non-detectable\nwithin-sample ACP5/MMP9 contrast",
    TEALP, TEAL)
box(BX, ys[3], BW, BH, "OBSERVED ESTIMATE",
    "K\u0302\u1d62, sample-specific scRNA estimate\n(locked primary: mean 0.639, CI 0.322\u20130.976)",
    TEALP, TEAL)
box(BX, ys[4], BW, BH, "CANDIDATE SCALABLE PROXY",
    "myeloid aggregate \u00b7 bulk RNA\n(indirect proxy-architecture audit, unpaired data)",
    ORANGEP, ORANGE)
box(BX, ys[5], BW, BH, "REPRESENTATION QUALIFICATION",
    "bulk uninformative for the benchmark;\naggregate/spatial insufficiently precise",
    GREYP, GREY)
box(BX, ys[6], BW, BH, "CLINICAL QUALIFICATION",
    "treatment-by-biomarker \u00b7 patient selection\n(outside this study)",
    REDP, RED)

for k in range(len(ys) - 1):
    arrow(BX + BW / 2, ys[k], ys[k + 1] + BH)

# Right-side coverage markers
cov = [
    (ys[0] + BH / 2, "Directly evaluated"),
    (ys[1] + BH / 2, "Explicitly defined"),
    (ys[2] + BH / 2, "Implemented"),
    (ys[3] + BH / 2, "Generated in locked analysis;\ncharacterised post hoc"),
    (ys[4] + BH / 2, "Indirect architecture audit only"),
    (ys[5] + BH / 2, "Not established;\npaired validation required"),
    (ys[6] + BH / 2, "Not tested"),
]
for y, label in cov:
    ax.text(BX + BW + 0.35, y, label, ha="left", va="center", fontsize=6.6,
            color="#333333")

ax.text(0.55, 9.92, "Cross-scale construct-to-readout chain", fontsize=9.5,
        fontweight="bold", color="#111111")
ax.text(0.55, 0.42, "Direct representation validation would require paired measurement of the coupling "
        "measurand and a scalable readout in the same patients.",
        fontsize=6.6, color="#555555")

for fmt in ("png", "svg", "pdf"):
    fig.savefig(r"./prior_analysis_v2\Figures\Fig1." + fmt,
                dpi=400, bbox_inches="tight")
print("Fig1 schematic written")
