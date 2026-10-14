"""V0.5 Points 3 and 4.

Point 3  estimability selection: systematic comparison of the 13 within-state estimable and 7 non-estimable samples,
       and reported separately for the two missingness mechanisms (too few TREM2 cells / zero positive cells).
Point 4  single-gene decomposition of ACP5 and MMP9 within TREM2, for claim audit only;
       scores follow the within-sample all-myeloid z-score standard of the primary analysis; the primary outcome is not redefined.
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

from run_V04_core import (
    OUT, ANNOTATION, H5ADS, TUMOUR_TISSUES, MIN_MYELOID_PER_SAMPLE,
    TREM2_LABEL, zscore, derive_rng, N_BOOT,
)

P1A2_MIN_TREM2, P1A2_MIN_POS, P1A2_MIN_NEG = 10, 2, 5


def boot_mean(vals: np.ndarray, datasets: np.ndarray, stream: int) -> dict:
    order = np.argsort(datasets)
    vals, datasets = vals[order], datasets[order]
    strata = [vals[datasets == d] for d in np.unique(datasets)]
    strata = [s for s in strata if len(s)]
    if not strata:
        return {"point": float("nan"), "ci_low": float("nan"), "ci_high": float("nan")}
    point = float(np.mean(np.concatenate(strata)))
    rng = derive_rng(stream)
    boot = np.empty(N_BOOT)
    for i in range(N_BOOT):
        boot[i] = np.mean(np.concatenate([rng.choice(s, len(s), replace=True) for s in strata]))
    return {"point": point, "ci_low": float(np.quantile(boot, 0.025)),
            "ci_high": float(np.quantile(boot, 0.975)), "n_patients": int(len(vals)),
            "n_datasets": int(len(strata))}


def main() -> None:
    labels = pd.read_csv(ANNOTATION, index_col=0, low_memory=False)
    rows = []
    for path in H5ADS:
        obj = ad.read_h5ad(path, backed="r")
        dataset = str(obj.obs["dataset"].astype(str).mode().iloc[0])
        mask = (obj.obs["tissue"].astype(str).isin(TUMOUR_TISSUES)
                & obj.obs["celltype"].astype(str).eq("Myeloid"))
        obs = obj.obs.loc[mask].copy()
        obs.index = obs.index.astype(str)
        lab = labels.loc[labels["dataset"].astype(str).eq(dataset)
                         & labels["tissue"].astype(str).isin(TUMOUR_TISSUES),
                         ["myeloid_subtype"]].copy()
        lab.index = lab.index.astype(str)
        overlap = obs.index.intersection(lab.index)
        gene_idx = [obj.var_names.get_loc(g) for g in ["SIGLEC15", "ACP5", "MMP9"]]
        row_idx = obj.obs_names.get_indexer(overlap)
        x = obj.X[row_idx, :][:, gene_idx]
        if sparse.issparse(x):
            x = x.toarray()
        expr = pd.DataFrame(np.asarray(x, dtype=float), index=overlap,
                            columns=["SIGLEC15", "ACP5", "MMP9"])
        obj.file.close()
        obs = obs.loc[overlap].join(lab.loc[overlap])
        expr = expr.loc[overlap]

        for sample_id, ids in obs.groupby("sample_id", observed=True).groups.items():
            ids = pd.Index(ids)
            if len(ids) < MIN_MYELOID_PER_SAMPLE:
                continue
            o, e = obs.loc[ids], expr.loc[ids]
            positive = e["SIGLEC15"].to_numpy(float) > 0
            sub = o["myeloid_subtype"].astype(str).to_numpy()
            trem2 = sub == TREM2_LABEL
            a_score = zscore(e["ACP5"].to_numpy(float))
            m_score = zscore(e["MMP9"].to_numpy(float))
            counts = np.asarray(o["total_counts"].to_numpy(float), dtype=float)

            n_t = int(trem2.sum()); n_tp = int((trem2 & positive).sum())
            elig = (n_t >= P1A2_MIN_TREM2 and n_tp >= P1A2_MIN_POS
                    and int((trem2 & ~positive).sum()) >= P1A2_MIN_NEG)
            k_a = k_m = k_core = np.nan
            if n_tp > 0 and int((trem2 & ~positive).sum()) > 0:
                core = np.nanmean(np.column_stack([a_score, m_score]), axis=1)
                t_pos, t_neg = positive[trem2], ~positive[trem2]
                k_a = float(np.nanmean(a_score[trem2][t_pos]) - np.nanmean(a_score[trem2][t_neg]))
                k_m = float(np.nanmean(m_score[trem2][t_pos]) - np.nanmean(m_score[trem2][t_neg]))
                k_core = float(np.nanmean(core[trem2][t_pos]) - np.nanmean(core[trem2][t_neg]))

            rows.append({
                "dataset": dataset, "sample_id": str(sample_id),
                "n_myeloid": int(len(ids)),
                "median_total_counts": float(np.median(counts)),
                "median_n_genes": float(np.median(np.asarray(o["n_genes_by_counts"].to_numpy(float)))),
                "A_i": float(positive.mean()),
                "O_i": float(trem2.mean()),
                "n_trem2": n_t, "n_trem2_positive": n_tp,
                "n_trem2_negative": int((trem2 & ~positive).sum()),
                "primary_eligible": bool(elig),
                "K_i_T_ACP5": k_a,
                "K_i_T_MMP9": k_m,
                "K_i_T_core": k_core,
            })

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "V05_WITHIN_STATE_GENE_DECOMPOSITION.csv", index=False, encoding="utf-8")

    # ---------- Point 4: single-gene decomposition summary ----------
    prim = df.loc[df["primary_eligible"]]
    decomp = {}
    for col in ("K_i_T_ACP5", "K_i_T_MMP9", "K_i_T_core"):
        b = boot_mean(prim[col].dropna().to_numpy(float),
                      prim.loc[prim[col].notna(), "dataset"].to_numpy(), stream=41)
        vals = prim[col].dropna().to_numpy(float)
        b["n_positive"] = int((vals > 0).sum())
        b["median"] = float(np.median(vals))
        b["verdict_gt0"] = bool(np.isfinite(b["ci_low"]) and b["ci_low"] > 0)
        decomp[col] = b
    # Whether the two genes are comparable in the same set of patients: paired difference
    paired = prim.dropna(subset=["K_i_T_ACP5", "K_i_T_MMP9"])
    diff = paired["K_i_T_MMP9"] - paired["K_i_T_ACP5"]
    decomp["paired_MMP9_minus_ACP5"] = {
        "n_patients": int(len(paired)),
        "mean": float(diff.mean()), "median": float(diff.median()),
        "wilcoxon_p": float(stats.wilcoxon(diff, alternative="two-sided").pvalue)
        if len(diff) > 5 else None,
        "n_mmp9_larger": int((diff > 0).sum()),
    }
    (OUT / "_v05_gene_decomposition.json").write_text(
        json.dumps(decomp, indent=2, ensure_ascii=False), encoding="utf-8")

    # ---------- Point 3: estimability selection ----------
    df["missing_reason"] = np.where(
        df["primary_eligible"], "estimable",
        np.where(df["n_trem2"] < P1A2_MIN_TREM2, "too_few_TREM2_cells",
                 np.where(df["n_trem2_positive"] < P1A2_MIN_POS,
                          "too_few_SIGLEC15pos_TREM2_cells", "too_few_TREM2neg_cells")))
    sel = df.groupby("missing_reason").agg(
        n=("sample_id", "size"),
        median_n_myeloid=("n_myeloid", "median"),
        median_median_counts=("median_total_counts", "median"),
        median_A_i=("A_i", "median"),
        median_O_i=("O_i", "median"),
        median_n_trem2=("n_trem2", "median"),
        median_n_trem2_pos=("n_trem2_positive", "median"),
        datasets=("dataset", lambda s: sorted(set(s))),
    ).reset_index()
    sel.to_csv(OUT / "V05_ESTIMABILITY_SELECTION.csv", index=False, encoding="utf-8")

    est = df.loc[df["primary_eligible"]]
    non = df.loc[~df["primary_eligible"]]
    sel_test = {}
    for col in ("n_myeloid", "median_total_counts", "A_i", "O_i", "n_trem2"):
        a, b = est[col].to_numpy(float), non[col].to_numpy(float)
        try:
            p = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
        except Exception:  # noqa: BLE001
            p = None
        sel_test[col] = {"estimable_median": float(np.median(a)),
                         "non_estimable_median": float(np.median(b)),
                         "mannwhitney_p": p}
    (OUT / "_v05_estimability_selection.json").write_text(
        json.dumps({"by_reason": sel.to_dict(orient="records"),
                    "estimable_vs_non": sel_test,
                    "n_estimable": int(len(est)), "n_non_estimable": int(len(non))},
                   indent=2, ensure_ascii=False), encoding="utf-8")

    print(json.dumps({
        "n_estimable": int(len(est)), "n_non_estimable": int(len(non)),
        "missing_reason_counts": df["missing_reason"].value_counts().to_dict(),
        "gene_decomposition": decomp,
        "selection_tests": sel_test,
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
