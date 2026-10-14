"""Coupling Audit: final assembly (K_i vs intensity; complete report)."""
import io
import json

import pandas as pd
from scipy import stats as st

J = r"./prior_analysis_v2"
RES = J + r"\deposit_v24\results"

h5 = json.load(io.open(r"./temp\COUPLING_AUDIT_h5ad.json", encoding="utf-8"))
rows = []
for ds, recs in h5.items():
    if ds == "obs_columns":
        continue
    rows.extend(recs)
inten = pd.DataFrame(rows)
inten = inten.dropna(subset=["mean_positive_intensity"])

core = pd.read_csv(RES + r"\cell_level_sample_effects.csv")
core = core[core.program == "Derosiers_ACP5_MMP9_core"].copy()
core["prevalence"] = core.n_positive / core.n_myeloid

m = core.merge(inten[["dataset", "sample_id", "mean_positive_intensity"]],
               on=["dataset", "sample_id"], how="left").dropna(subset=["mean_positive_intensity"])
r_i, p_i = st.spearmanr(m.difference_positive_minus_negative, m.mean_positive_intensity)
r_ip, p_ip = st.spearmanr(m.prevalence, m.mean_positive_intensity)

res = {
    "ki_vs_intensity": {"n": int(len(m)), "spearman": float(r_i), "p": float(p_i)},
    "prevalence_vs_intensity": {"n": int(len(m)), "spearman": float(r_ip), "p": float(p_ip)},
}
with open(RES + r"\COUPLING_AUDIT_summary.json", encoding="utf-8") as fh:
    rep = json.load(fh)
rep["ki_vs_intensity"] = res["ki_vs_intensity"]
rep["prevalence_vs_intensity"] = res["prevalence_vs_intensity"]
rep["trem2_fraction_note"] = ("Per-sample TREM2_TAM fractions were not computable: the deposited objects "
                              "carry top-level celltype only; the fine-stratification annotation used for Table S12 "
                              "is not present in the deposited obs.")
with open(RES + r"\COUPLING_AUDIT_summary.json", "w", encoding="utf-8") as fh:
    json.dump(rep, fh, indent=2)
print("ki vs intensity:", round(r_i, 3), "p:", round(p_i, 4))
print("prevalence vs intensity:", round(r_ip, 3), "p:", round(p_ip, 4))
print("n:", int(len(m)))
