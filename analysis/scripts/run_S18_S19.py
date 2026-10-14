"""S18 (single-gene decomposition) and S19 (depth-matched sensitivity) companion.

Reuses the locked S17 pipeline statistical infrastructure but provides its own
HDF5 reader to bypass the anndata<->h5ad version-incompatibility in /uns/log1p.

All random seeds, bootstrap draws, and inferential frameworks match S17.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy
from scipy import sparse, stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import analysis_F_cross_scale as s17  # noqa: E402

from analysis_F_cross_scale import (  # noqa: E402
    ALL_GENES, BULK_PATHS, DESCRIPTIVE, N_RESAMPLES, OUT, POSTHOC_MYELOID_COMPOSITION, PRIMARY,
    PROGRAMS, RNG, SCRNA_PATHS, SEED, SOFT_PATHS, SPATIAL_PATH, TARGET, TUMOUR_TISSUES,
    holm_adjust, json_default, load_bulk, parse_soft_platform,
    permutation_p, pooled_within_dataset_rank, random_effects_meta, required_gene_count,
    residualize_depth, score_program, sha256_file, spearman_bootstrap, stratified_bootstrap_mean,
    stratified_rank_inference, zscore_columns,
)

S18_OUT = HERE / "outputs_S18"
S19_OUT = HERE / "outputs_S19"


# ---- Custom HDF5 h5ad reader (bypasses anndata version incompatibility) ----


def _read_h5ad_custom(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, list[str]]:
    """Read .h5ad via h5py: return (obs_df, obs_index, X_dense, var_names)."""
    with h5py.File(path, "r") as f:
        # read X as dense numpy array from CSR
        if "X" in f:
            data = f["X/data"][:]
            indices = f["X/indices"][:]
            indptr = f["X/indptr"][:]
            shape = (len(indptr) - 1, len(f["var/_index"][:]))
            x_sparse = sparse.csr_matrix((data, indices, indptr), shape=shape)
            x_dense = x_sparse.toarray()
        else:
            raise KeyError("No X matrix in h5ad")

        # read var names
        var_names = [name.decode() if isinstance(name, bytes) else name
                     for name in f["var/_index"][:]]

        # read obs
        obs = f["obs"]
        obs_data = {}
        obs_index = [idx.decode() if isinstance(idx, bytes) else idx
                     for idx in obs["_index"][:]]

        for key in obs.keys():
            if key == "_index":
                continue
            item = obs[key]
            if isinstance(item, h5py.Dataset):
                values = item[:]
                if values.dtype.kind == "S":
                    values = [v.decode() if isinstance(v, bytes) else v for v in values]
                elif values.dtype.kind == "O":
                    try:
                        values = [v.decode() if isinstance(v, bytes) else v for v in values]
                    except Exception:
                        pass
                obs_data[key] = values
            elif isinstance(item, h5py.Group):
                # categorical: categories + codes
                if "categories" in item and "codes" in item:
                    cats = item["categories"][:]
                    codes = item["codes"][:]
                    if cats.dtype.kind == "S":
                        cats = [c.decode() if isinstance(c, bytes) else c for c in cats]
                    obs_data[key] = [cats[code] if 0 <= code < len(cats) else np.nan
                                     for code in codes]
                else:
                    obs_data[key] = np.full(len(obs_index), np.nan)
            else:
                obs_data[key] = np.full(len(obs_index), np.nan)

        obs_df = pd.DataFrame(obs_data, index=obs_index)

    return obs_df, x_dense, var_names


def extract_h5ad_custom(path: Path, spatial: bool = False):
    """Extract myeloid cells from scRNA or spatial h5ad.

    Returns (obs_df, expr_df, dataset_name) replicating the S17 extract_h5ad contract.
    """
    obs_df, x, var_names = _read_h5ad_custom(path)

    if spatial:
        mask = obs_df["spot_type"].astype(str).eq("Myeloid")
        dataset = "GSE282302"
    else:
        mask = (obs_df["tissue"].astype(str).isin(TUMOUR_TISSUES) &
                obs_df["celltype"].astype(str).eq("Myeloid"))
        dataset = str(obs_df.loc[mask, "dataset"].astype(str).iloc[0])

    mask_array = mask.to_numpy()
    present_genes = [g for g in ALL_GENES if g in var_names]
    gene_indices = [var_names.index(g) for g in present_genes]

    expr = pd.DataFrame(
        x[mask_array][:, gene_indices],
        columns=present_genes,
    )
    obs = obs_df.loc[mask].reset_index(drop=True)
    expr = expr.reset_index(drop=True)

    return obs, expr, dataset


# ---- S18: Single-gene decomposition ----


def _score_single_gene(frame: pd.DataFrame, gene: str):
    if gene not in frame.columns:
        return np.full(len(frame), np.nan), []
    values = frame[gene].to_numpy(float)
    mean = np.nanmean(values)
    sd = np.nanstd(values, ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return np.full(len(frame), np.nan), [gene]
    z = (values - mean) / sd
    return np.nan_to_num(z, nan=0.0), [gene]


def run_s18_single_gene_decomposition() -> dict:
    """Run the full S17 pipeline with ACP5 alone and MMP9 alone."""

    S18_OUT.mkdir(parents=True, exist_ok=True)
    results = {}

    single_genes = {"ACP5": ["ACP5"], "MMP9": ["MMP9"]}

    for gene_name, gene_list in single_genes.items():
        print(f"\n{'='*60}")
        print(f"  S18: single-gene run for {gene_name}")
        print(f"{'='*60}")

        # ---- scRNA cell-level ----
        sample_effect_rows = []
        aggregate_rows = []
        for path in SCRNA_PATHS:
            obs, expr, dataset = extract_h5ad_custom(path)
            for sample_id, idx in obs.groupby("sample_id", observed=True).groups.items():
                idx_arr = np.asarray(list(idx), dtype=int)
                if len(idx_arr) < 50:
                    continue
                sample_obs = obs.iloc[idx_arr].reset_index(drop=True)
                sample_expr = expr.iloc[idx_arr].reset_index(drop=True)
                positive = sample_expr[TARGET].to_numpy(float) > 0

                aggregate = np.log1p(np.mean(
                    np.expm1(sample_expr.to_numpy(float)), axis=0
                ))
                agg_row = {
                    "dataset": dataset, "sample_id": str(sample_id),
                    "n_myeloid": len(idx_arr),
                    "n_siglec15_positive": int(positive.sum()),
                }
                agg_row.update(dict(zip(sample_expr.columns, aggregate)))
                aggregate_rows.append(agg_row)

                if positive.sum() == 0 or (~positive).sum() == 0:
                    continue

                score, present = _score_single_gene(sample_expr, gene_name)
                if not np.isfinite(score).any():
                    continue
                residual = s17.residualize_depth(score, sample_obs)
                sample_effect_rows.append({
                    "dataset": dataset, "sample_id": str(sample_id),
                    "program": gene_name, "gene": gene_name,
                    "n_myeloid": len(idx_arr),
                    "n_positive": int(positive.sum()),
                    "n_negative": int((~positive).sum()),
                    "n_genes_present": len(present),
                    "difference_positive_minus_negative": float(
                        np.nanmean(score[positive]) - np.nanmean(score[~positive])
                    ),
                    "depth_adjusted_difference": float(
                        np.nanmean(residual[positive]) - np.nanmean(residual[~positive])
                    ),
                })

        effects = pd.DataFrame(sample_effect_rows)
        aggregates = pd.DataFrame(aggregate_rows)

        # Cell-level summary
        values = effects["difference_positive_minus_negative"].to_numpy(float)
        ci_lo, ci_hi = stratified_bootstrap_mean(
            effects, "difference_positive_minus_negative"
        )
        w_test = stats.wilcoxon(values, alternative="two-sided", zero_method="wilcox")
        adj_vals = effects["depth_adjusted_difference"].to_numpy(float)
        adj_lo, adj_hi = stratified_bootstrap_mean(
            effects, "depth_adjusted_difference"
        )
        adj_w_test = stats.wilcoxon(adj_vals, alternative="two-sided", zero_method="wilcox")

        cell_summary = {"n_samples": len(effects),
                        "mean_sample_difference": float(values.mean()),
                        "ci_low": ci_lo, "ci_high": ci_hi,
                        "wilcoxon_p": float(w_test.pvalue),
                        "depth_adjusted_mean_difference": float(adj_vals.mean()),
                        "depth_adjusted_ci_low": adj_lo,
                        "depth_adjusted_ci_high": adj_hi,
                        "depth_adjusted_wilcoxon_p": float(adj_w_test.pvalue)}

        # ---- Myeloid aggregate ----
        association_rows = []
        for dataset_tag, group in aggregates.groupby("dataset"):
            score, present = _score_single_gene(group, gene_name)
            if not np.isfinite(score).any():
                continue
            rho, p, lo, hi = spearman_bootstrap(
                group[TARGET].to_numpy(float), score
            )
            pp = permutation_p(group[TARGET].to_numpy(float), score, rho)
            association_rows.append({
                "scope": "dataset", "dataset": dataset_tag,
                "program": gene_name, "n_samples": len(group),
                "rho": rho, "ci_low": lo, "ci_high": hi,
                "permutation_p": pp,
            })

        aggregate_scored = []
        for dataset_tag, group in aggregates.groupby("dataset"):
            group = group.copy()
            group[gene_name], _ = _score_single_gene(group, gene_name)
            aggregate_scored.append(group)
        aggregate_scored = pd.concat(aggregate_scored, ignore_index=True)
        rho_p, lo_p, hi_p, p_p = stratified_rank_inference(
            aggregate_scored, TARGET, gene_name
        )
        association_rows.append({
            "scope": "pooled_within_dataset_ranks", "dataset": "ALL",
            "program": gene_name, "n_samples": len(aggregate_scored),
            "rho": rho_p, "ci_low": lo_p, "ci_high": hi_p,
            "permutation_p": p_p,
        })
        agg_assoc = pd.DataFrame(association_rows)
        pooled = agg_assoc[
            agg_assoc["scope"] == "pooled_within_dataset_ranks"
        ].iloc[0]

        # ---- Bulk ----
        cohorts = load_bulk()
        bulk_rows = []
        for cohort_tag, genes_by_samples in cohorts.items():
            frame = genes_by_samples.T
            score, present = _score_single_gene(frame, gene_name)
            if not np.isfinite(score).any() or TARGET not in frame:
                continue
            x = frame[TARGET].to_numpy(float)
            rho, p, lo, hi = spearman_bootstrap(x, score)
            pp = permutation_p(x, score, rho)
            bulk_rows.append({
                "cohort": cohort_tag, "program": gene_name,
                "n_samples": len(frame), "rho": rho,
                "ci_low": lo, "ci_high": hi, "permutation_p": pp,
            })
        bulk_results = pd.DataFrame(bulk_rows)
        meta = pd.DataFrame([
            random_effects_meta(g) for _, g in bulk_results.groupby("program")
        ])

        # ---- Spatial ----
        obs_sp, expr_sp, ds_sp = extract_h5ad_custom(SPATIAL_PATH, spatial=True)
        capture_rows = []
        for sp_id, idx_sp in obs_sp.groupby("sample_id", observed=True).groups.items():
            idx_arr = np.asarray(list(idx_sp), dtype=int)
            capture_obs = obs_sp.iloc[idx_arr].reset_index(drop=True)
            capture_expr = expr_sp.iloc[idx_arr].reset_index(drop=True)
            pos = capture_expr[TARGET].to_numpy(float) > 0
            if len(idx_arr) < 10 or pos.sum() < 2 or (~pos).sum() < 2:
                continue
            case = str(capture_obs["case"].iloc[0])
            score_sp, present_sp = _score_single_gene(capture_expr, gene_name)
            if not np.isfinite(score_sp).any():
                continue
            capture_rows.append({
                "case": case, "sample_id": str(sp_id),
                "program": gene_name,
                "n_myeloid_spots": len(idx_arr),
                "n_positive": int(pos.sum()), "n_negative": int((~pos).sum()),
                "difference_positive_minus_negative": float(
                    np.nanmean(score_sp[pos]) - np.nanmean(score_sp[~pos])
                ),
            })

        sp_captures = pd.DataFrame(capture_rows)
        sp_cases = sp_captures.groupby(
            ["case", "program"], as_index=False
        ).agg(
            n_captures=("sample_id", "nunique"),
            mean_capture_difference=("difference_positive_minus_negative", "mean"),
        )

        sp_row = {}
        if len(sp_cases) > 0:
            val = sp_cases["mean_capture_difference"].to_numpy(float)
            crit = stats.t.ppf(0.975, df=len(val) - 1)
            mean_v = val.mean()
            se_v = val.std(ddof=1) / np.sqrt(len(val))
            sp_row = {
                "n_cases": len(val),
                "mean_case_difference": float(mean_v),
                "ci_low": float(mean_v - crit * se_v),
                "ci_high": float(mean_v + crit * se_v),
                "t_test_p": float(stats.ttest_1samp(val, 0).pvalue),
            }

        # ---- Collate ----
        meta_row = meta.iloc[0]
        results[gene_name] = {
            "cell_myeloid": {
                "mean_difference": cell_summary["mean_sample_difference"],
                "ci_low": cell_summary["ci_low"],
                "ci_high": cell_summary["ci_high"],
                "wilcoxon_p": cell_summary["wilcoxon_p"],
            },
            "cell_depth_adjusted": {
                "mean_difference": cell_summary["depth_adjusted_mean_difference"],
                "ci_low": cell_summary["depth_adjusted_ci_low"],
                "ci_high": cell_summary["depth_adjusted_ci_high"],
                "wilcoxon_p": cell_summary["depth_adjusted_wilcoxon_p"],
            },
            "aggregate_pooled": {
                "rho": float(pooled["rho"]),
                "ci_low": float(pooled["ci_low"]),
                "ci_high": float(pooled["ci_high"]),
                "permutation_p": float(pooled["permutation_p"]),
            },
            "bulk_meta": {
                "r_random": float(meta_row["r_random"]),
                "ci_low": float(meta_row["ci_low"]),
                "ci_high": float(meta_row["ci_high"]),
                "p_mhk": float(meta_row["p_modified_hartung_knapp"]),
                "n_cohorts": int(meta_row["n_cohorts"]),
            },
            "spatial": {
                "n_cases": sp_row.get("n_cases", 0),
                "mean_case_difference": float(sp_row.get("mean_case_difference", np.nan)),
                "ci_low": float(sp_row.get("ci_low", np.nan)),
                "ci_high": float(sp_row.get("ci_high", np.nan)),
                "t_test_p": float(sp_row.get("t_test_p", np.nan)),
            },
        }

        effects.to_csv(S18_OUT / f"{gene_name}_cell_effects.csv", index=False)
        pd.DataFrame([cell_summary]).to_csv(
            S18_OUT / f"{gene_name}_cell_summary.csv", index=False
        )
        agg_assoc.to_csv(S18_OUT / f"{gene_name}_aggregate_associations.csv", index=False)
        bulk_results.to_csv(S18_OUT / f"{gene_name}_bulk_associations.csv", index=False)
        meta.to_csv(S18_OUT / f"{gene_name}_bulk_meta.csv", index=False)
        if len(sp_captures) > 0:
            sp_captures.to_csv(S18_OUT / f"{gene_name}_spatial_captures.csv", index=False)

    (S18_OUT / "S18_summary.json").write_text(
        json.dumps(results, indent=2, ensure_ascii=False, default=json_default),
        encoding="utf-8",
    )
    return results


# ---- S19: Depth-matched sensitivity ----


def run_s19_depth_matched() -> dict:
    """Greedy nearest-neighbour depth-matching within each sample."""

    S19_OUT.mkdir(parents=True, exist_ok=True)
    print(f"\n{'='*60}")
    print(f"  S19: depth-matched sensitivity analysis")
    print(f"{'='*60}")

    sample_effect_rows = []
    all_match_counts = []

    for path in SCRNA_PATHS:
        obs, expr, dataset = extract_h5ad_custom(path)
        for sample_id, idx in obs.groupby("sample_id", observed=True).groups.items():
            idx_arr = np.asarray(list(idx), dtype=int)
            if len(idx_arr) < 50:
                continue
            sample_obs = obs.iloc[idx_arr].reset_index(drop=True)
            sample_expr = expr.iloc[idx_arr].reset_index(drop=True)

            positive_mask = sample_expr[TARGET].to_numpy(float) > 0
            negative_mask = ~positive_mask

            n_pos = positive_mask.sum()
            n_neg = negative_mask.sum()
            if n_pos == 0 or n_neg == 0:
                continue

            if "total_counts" in sample_obs:
                log_depth_pos = np.log1p(
                    pd.to_numeric(
                        sample_obs.loc[positive_mask, "total_counts"],
                        errors="coerce",
                    ).to_numpy(float)
                )
                log_depth_neg = np.log1p(
                    pd.to_numeric(
                        sample_obs.loc[negative_mask, "total_counts"],
                        errors="coerce",
                    ).to_numpy(float)
                )
            else:
                log_depth_pos = np.ones(n_pos)
                log_depth_neg = np.ones(n_neg)

            pos_indices = np.where(positive_mask)[0]
            neg_indices = np.where(negative_mask)[0]

            matched_pairs = []
            available_neg = set(range(len(neg_indices)))

            for i in range(len(pos_indices)):
                if not available_neg:
                    break
                d_pos = log_depth_pos[i]
                best_j, best_dist = None, np.inf
                for j in list(available_neg):
                    dist = abs(d_pos - log_depth_neg[j])
                    if dist < best_dist:
                        best_dist = dist
                        best_j = j
                if best_j is not None:
                    matched_pairs.append((pos_indices[i], neg_indices[best_j]))
                    available_neg.remove(best_j)

            if not matched_pairs:
                continue

            n_pairs = len(matched_pairs)
            all_match_counts.append({
                "dataset": dataset, "sample_id": str(sample_id),
                "n_positive": n_pos, "n_negative_before_match": n_neg,
                "n_matched_pairs": n_pairs,
            })

            pos_matched_idx_list = [p[0] for p in matched_pairs]
            neg_matched_idx_list = [p[1] for p in matched_pairs]

            all_matched_idx = pos_matched_idx_list + neg_matched_idx_list
            all_expr = sample_expr.iloc[all_matched_idx]
            n_pos = len(pos_matched_idx_list)

            for prog_name, prog_genes in PROGRAMS.items():
                all_score, all_present = s17.score_program(all_expr, prog_name)
                if not np.isfinite(all_score).any():
                    continue
                pos_score = all_score[:n_pos]
                neg_score = all_score[n_pos:]
                sample_effect_rows.append({
                    "dataset": dataset, "sample_id": str(sample_id),
                    "program": prog_name,
                    "n_matched_pairs": n_pairs,
                    "n_genes_present": len(all_present),
                    "difference_positive_minus_negative": float(
                        np.nanmean(pos_score) - np.nanmean(neg_score)
                    ),
                })

    effects = pd.DataFrame(sample_effect_rows)
    match_info = pd.DataFrame(all_match_counts)

    summaries = []
    for prog_name, group in effects.groupby("program"):
        values = group["difference_positive_minus_negative"].to_numpy(float)
        ci_lo, ci_hi = stratified_bootstrap_mean(
            group, "difference_positive_minus_negative"
        )
        w_test = stats.wilcoxon(values, alternative="two-sided", zero_method="wilcox")
        total_pairs = int(group["n_matched_pairs"].sum())
        summaries.append({
            "program": prog_name,
            "n_matched_pairs_total": total_pairs,
            "n_samples": len(group),
            "n_datasets": group["dataset"].nunique(),
            "mean_difference": float(values.mean()),
            "ci_low": ci_lo, "ci_high": ci_hi,
            "wilcoxon_p": float(w_test.pvalue),
        })

    depth_summary = pd.DataFrame(summaries)
    secondary = depth_summary["program"].ne(PRIMARY)
    depth_summary.loc[secondary, "holm_p_secondary"] = holm_adjust(
        depth_summary.loc[secondary, "wilcoxon_p"].tolist()
    )

    effects.to_csv(S19_OUT / "depth_matched_cell_effects.csv", index=False)
    depth_summary.to_csv(S19_OUT / "depth_matched_summary.csv", index=False)
    match_info.to_csv(S19_OUT / "match_counts.csv", index=False)

    results = {}
    for _, row in depth_summary.iterrows():
        results[row["program"]] = {
            "n_matched_pairs": int(row["n_matched_pairs_total"]),
            "n_samples": int(row["n_samples"]),
            "mean_difference": float(row["mean_difference"]),
            "ci_low": float(row["ci_low"]),
            "ci_high": float(row["ci_high"]),
            "wilcoxon_p": float(row["wilcoxon_p"]),
        }

    (S19_OUT / "S19_summary.json").write_text(
        json.dumps(results, indent=2, ensure_ascii=False, default=json_default),
        encoding="utf-8",
    )
    return results


def main():
    print(f"SEED = {SEED}  N_RESAMPLES = {N_RESAMPLES}")
    print(f"Running S18 (single-gene decomposition) and S19 (depth-matched sensitivity)")
    print(f"S18 outputs -> {S18_OUT}")
    print(f"S19 outputs -> {S19_OUT}")

    s18 = run_s18_single_gene_decomposition()
    s19 = run_s19_depth_matched()

    print("\n" + "=" * 60)
    print("  S18 SUMMARY (single-gene decomposition)")
    print("=" * 60)
    print(json.dumps(s18, indent=2, default=json_default))

    print("\n" + "=" * 60)
    print("  S19 SUMMARY (depth-matched sensitivity)")
    print("=" * 60)
    print(json.dumps(s19, indent=2, default=json_default))

    print(
        "\nDone. S18 and S19 results ready for population into "
        "Supplementary Tables S28 and S29."
    )


if __name__ == "__main__":
    main()
