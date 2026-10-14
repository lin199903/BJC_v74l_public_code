import io
import json

import pandas as pd
from scipy import stats as st

J = r"./prior_analysis_v2"
RES = J + r"\deposit_v24\results"

core = pd.read_csv(RES + r"\cell_level_sample_effects.csv")
core = core[core.program == "Derosiers_ACP5_MMP9_core"].copy()
matched = pd.read_csv(RES + r"\depth_matched_cell_effects.csv")
matched = matched[matched.program == "Derosiers_ACP5_MMP9_core"].copy()
matched["sample_id"] = matched["sample_id"].astype(str).str.replace("b'", "").str.replace("'", "")

mm = core.merge(matched[["dataset", "sample_id", "difference_positive_minus_negative"]],
                on=["dataset", "sample_id"], suffixes=("", "_matched"))
mm = mm.dropna(subset=["difference_positive_minus_negative", "difference_positive_minus_negative_matched"])
r_m, p_m = st.spearmanr(mm.difference_positive_minus_negative, mm.difference_positive_minus_negative_matched)
print("matched n:", len(mm), "spearman:", round(r_m, 3), "p:", round(p_m, 4))
mm[["dataset", "sample_id", "difference_positive_minus_negative",
    "difference_positive_minus_negative_matched"]].to_csv(RES + r"\COUPLING_AUDIT_ki_matched.csv", index=False)

# update summary json with the matched rank stability
with open(RES + r"\COUPLING_AUDIT_summary.json", encoding="utf-8") as fh:
    rep = json.load(fh)
rep["rank_stability_unadj_vs_matched"] = {"n": int(len(mm)), "spearman": float(r_m), "p": float(p_m)}
with open(RES + r"\COUPLING_AUDIT_summary.json", "w", encoding="utf-8") as fh:
    json.dump(rep, fh, indent=2)
print("summary updated")
