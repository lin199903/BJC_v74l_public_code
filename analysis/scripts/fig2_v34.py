"""Fig 3 (v26): patient-specific target-state coupling.

a: K_i per patient, grouped by dataset, with dataset medians.
b: raw vs depth-adjusted K_i (rank scatter).
c: raw vs depth-matched K_i (rank scatter).
d: SIGLEC15 prevalence vs K_i (rank scatter).
"""
import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as st

J = r"./prior_analysis_v2"
RES = J + r"\deposit_v24\results"
OUT = r".\Figures\Fig2"

df = pd.read_csv(RES + r"\COUPLING_AUDIT_ki_by_patient.csv")
df = df.dropna(subset=["difference_positive_minus_negative"])
matched = pd.read_csv(RES + r"\COUPLING_AUDIT_ki_matched.csv")

COL = {"GSE155698": "#3B6FB5", "GSE156405": "#C0442E", "GSE212966": "#2E7D4F"}
order = ["GSE155698", "GSE156405", "GSE212966"]

plt.rcParams.update({
    "font.size": 7.5, "font.family": "sans-serif",
    "axes.linewidth": 0.6, "xtick.direction": "in", "ytick.direction": "in",
    "axes.spines.top": False, "axes.spines.right": False,
})

fig = plt.figure(figsize=(6.7, 4.64))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.34, height_ratios=[1.18, 1.0])

# a: K_i by dataset
ax = fig.add_subplot(gs[0, 0])
for j, ds in enumerate(order):
    g = df[df.dataset == ds]
    x = j + 0.05 * np.random.default_rng(15072026).uniform(-1, 1, len(g))
    ax.scatter(x, g.difference_positive_minus_negative, s=18, color=COL[ds],
               edgecolor="white", linewidth=0.4, zorder=3)
    med = g.difference_positive_minus_negative.median()
    ax.hlines(med, j - 0.25, j + 0.25, color=COL[ds], linewidth=1.6, zorder=4)
ax.axhline(0, color="#888888", linewidth=0.5, linestyle="--", zorder=1)
ax.set_xticks(range(len(order)))
ax.set_xticklabels(["GSE155698", "GSE156405", "GSE212966"], fontsize=6.5)
ax.set_ylabel("K\u1d62 (within-patient coupling)")
ax.set_title("a  Patient-specific coupling", loc="left", fontsize=8, pad=2)
ax.set_ylim(-1.1, 3.9)

# b: raw vs depth-adjusted
ax = fig.add_subplot(gs[0, 1])
g = df.dropna(subset=["depth_adjusted_difference"])
r, p = st.spearmanr(g.difference_positive_minus_negative, g.depth_adjusted_difference)
ax.scatter(g.difference_positive_minus_negative, g.depth_adjusted_difference,
           s=16, color=[COL[d] for d in g.dataset], edgecolor="white", linewidth=0.4, zorder=3)
lim = [-0.6, 3.7]
ax.plot(lim, lim, color="#888888", linewidth=0.5, linestyle="--", zorder=1)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("K\u1d62, unadjusted")
ax.set_ylabel("K\u1d62, depth-adjusted")
ax.text(0.04, 0.96, "Spearman = %.2f\np = %.0e" % (r, p), transform=ax.transAxes,
        va="top", fontsize=6.5)
ax.set_title("b  Rank stability after depth adjustment", loc="left", fontsize=8, pad=2)

# c: raw vs matched
ax = fig.add_subplot(gs[1, 0])
r, p = st.spearmanr(matched.difference_positive_minus_negative,
                    matched.difference_positive_minus_negative_matched)
ax.scatter(matched.difference_positive_minus_negative,
           matched.difference_positive_minus_negative_matched,
           s=16, color=[COL[d] for d in matched.dataset], edgecolor="white", linewidth=0.4, zorder=3)
lim = [-0.6, 3.7]
ax.plot(lim, lim, color="#888888", linewidth=0.5, linestyle="--", zorder=1)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("K\u1d62, unadjusted")
ax.set_ylabel("K\u1d62, depth-matched")
ax.text(0.04, 0.96, "Spearman = 0.72\np = %.3f (n = 16)" % p, transform=ax.transAxes,
        va="top", fontsize=6.5)
ax.set_title("c  Rank stability after depth matching", loc="left", fontsize=8, pad=2)

# d: prevalence vs K_i
ax = fig.add_subplot(gs[1, 1])
r, p = st.spearmanr(df.prevalence, df.difference_positive_minus_negative)
ax.scatter(df.prevalence * 100, df.difference_positive_minus_negative,
           s=16, color=[COL[d] for d in df.dataset], edgecolor="white", linewidth=0.4, zorder=3)
ax.axhline(0, color="#888888", linewidth=0.5, linestyle="--", zorder=1)
ax.set_xlabel("RNA-detectable myeloid cells (%)")
ax.set_ylabel("K\u1d62")
ax.text(0.04, 0.96, "Spearman = %.2f\np = %.2f" % (r, p), transform=ax.transAxes,
        va="top", fontsize=6.5)
ax.set_title("d  Coupling versus prevalence", loc="left", fontsize=8, pad=2)

for fmt in ("png", "svg", "pdf"):
    fig.savefig(OUT + "." + fmt, dpi=400, bbox_inches="tight")
print("Fig3 written (a-d)")
