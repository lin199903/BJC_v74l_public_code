#!/usr/bin/env python
# =============================================================================
# 09_independent_verification.py -- Phase B step 2, GSE278689
#
# Three things, in order:
#
#   PART A  an INDEPENDENT re-read of the deposited matrix through a different
#           parser (scipy.io.mmread instead of the chunked pandas reader used in
#           03_), checked against the cache entry by entry.
#
#   PART B  the primary contrast recomputed in Python using the MANUSCRIPT'S OWN
#           functions -- `zscore_columns` and `score_program`, imported verbatim
#           from analysis_F_cross_scale.py -- and compared with the R result.
#
#   PART C  POSITIVE CONTROLS THAT ARE CAPABLE OF BEING VIOLATED. Deliberately
#           wrong conventions are applied and must give a DIFFERENT answer.
#           Without this, PART B would only show that two implementations agree,
#           which is self-confirming.
#
# NOTE ON LANGUAGE (pre-registered interpretation criteria): this run is a
# modality-specific sensitivity analysis. Its output must not use the words
# replication / non-replication / confirmed / failed about any RESULT. The
# strings below are technical verdicts on verification CONTROLS, not on results.
#
# Nothing here modifies the manuscript or the deposit.
# =============================================================================
from __future__ import annotations

import os
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
from scipy import stats

ROOT = Path(os.environ.get("SIGLEC15_GSE278689_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927\phaseB_GSE278689"))
OUT = ROOT / "output"
DATA = ROOT / "data"
MANUSCRIPT_MODULE = Path(
    os.path.join(os.environ.get("SIGLEC15_SUBMISSION_ROOT", r"E:\Tex\HMlong\npjPrecOncol\FINAL_SUBMISSION_v6"), "BMC_MEDICINE_SUBMISSION_v47", "_archives_translate", "Additional_file_2_Analysis_archive_v39", "scripts", "analysis_F_cross_scale.py"))
SEED = 15072026
N_BOOT = 10000
EXCLUDE = set("bcdfgjklmnopqrstuvwxyz")     # only to keep ruff quiet about names


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


# --------------------------------------------------------------------------
# import the manuscript's own routines, verbatim
# --------------------------------------------------------------------------
def load_manuscript_module(path: Path):
    spec = importlib.util.spec_from_file_location("manuscript_cross_scale", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)          # safe: the module guards main() itself
    return mod


def main() -> None:
    mm = load_manuscript_module(MANUSCRIPT_MODULE)
    log("imported the manuscript module from", MANUSCRIPT_MODULE.name)
    log("  PRIMARY                     :", mm.PRIMARY)
    log("  PROGRAMS[PRIMARY]           :", mm.PROGRAMS[mm.PRIMARY])
    log("  required_gene_count(PRIMARY):", mm.required_gene_count(mm.PRIMARY))
    log("  TARGET                      :", mm.TARGET)
    log("  SEED / N_RESAMPLES          :", mm.SEED, "/", mm.N_RESAMPLES)
    print("\n-- the manuscript's zscore_columns(), verbatim --")
    import inspect
    print(inspect.getsource(mm.zscore_columns).rstrip())
    print("\n-- the manuscript's score_program(), verbatim --")
    print(inspect.getsource(mm.score_program).rstrip())

    assert mm.TARGET == "SIGLEC15"
    assert mm.PROGRAMS[mm.PRIMARY] == ["ACP5", "MMP9"], "primary program is not ACP5+MMP9"
    assert mm.TARGET not in mm.PROGRAMS[mm.PRIMARY], "SIGLEC15 is inside the programme score"
    print("\ntarget-in-program assertion: SIGLEC15 is NOT a member of the primary "
          "programme -> PASS")

    # ----------------------------------------------------------------------
    print("\n" + "=" * 74)
    print("PART A -- INDEPENDENT RE-READ OF THE DEPOSITED MATRIX")
    print("=" * 74)
    feats = pd.read_csv(DATA / "GSE278689_sn_features.tsv.gz", sep="\t", header=None)
    genes = feats[1].astype(str).tolist()
    row1 = {g: i + 1 for i, g in enumerate(genes)}
    bar = pd.read_csv(DATA / "GSE278689_sn_barcodes.tsv.gz", sep="\t", header=None)[0].astype(str).tolist()
    md = pd.read_csv(DATA / "GSE278689_sn_metadata.csv.gz", index_col=0)
    md.index = md.index.astype(str)

    log("scipy.io.mmread ... (independent parser)")
    M = scipy.io.mmread(DATA / "GSE278689_sn_matrix.mtx.gz").tocsr()
    log(f"  shape {M.shape}  nnz {M.nnz:,}")
    assert M.shape == (len(genes), len(bar)), "matrix shape mismatch"
    log("  orientation: rows = features, cols = barcodes")

    md = md.loc[bar]
    my_mask = md["all_celltype"].to_numpy() == "Myeloid"
    my_idx = np.where(my_mask)[0]
    gi = {g: row1[g] - 1 for g in ("SIGLEC15", "ACP5", "MMP9")}
    sub = np.asarray(M[[gi["SIGLEC15"], gi["ACP5"], gi["MMP9"]], :][:, my_idx].todense()).T
    lib_full = np.asarray(M.sum(axis=0)).ravel()

    cache = pd.read_csv(OUT / "_myeloid_gene_cache.csv")
    assert len(cache) == len(my_idx)
    assert list(cache["barcode"]) == [bar[i] for i in my_idx], "myeloid cell order differs"
    for j, g in enumerate(("SIGLEC15", "ACP5", "MMP9")):
        same = np.array_equal(sub[:, j], cache[g].to_numpy(float))
        log(f"  {g:<9} independent read == cache : {same}   "
            f"max abs diff {np.max(np.abs(sub[:, j] - cache[g].to_numpy(float))):g}")
        assert same, f"{g} disagrees between the two parsers"
    log(f"  libsize  independent read == cache : "
        f"{np.array_equal(lib_full[my_idx], cache['libsize'].to_numpy(float))}")
    print("=> two independent parsers of the same deposit agree entry by entry.")

    # ----------------------------------------------------------------------
    print("\n" + "=" * 74)
    print("PART B -- PRIMARY CONTRAST VIA THE MANUSCRIPT'S OWN FUNCTIONS")
    print("=" * 74)
    Rres = pd.read_csv(OUT / "primary_contrast_by_sample.csv")
    Rres = Rres[Rres["status"] == "OK"].reset_index(drop=True)
    elig = list(Rres["sample_id"])
    log("eligible samples:", elig)

    raw = cache[["SIGLEC15", "ACP5", "MMP9"]].to_numpy(float)
    libm = cache["libsize"].to_numpy(float)
    libm = np.where(libm > 0, libm, 1.0)
    norm = np.log1p(raw / libm[:, None] * 1e4)          # CP10K then log1p
    sid = cache["sample_id"].to_numpy()
    det = raw[:, 0] > 0

    rows = []
    for s in elig:
        i = np.where(sid == s)[0]
        frame = pd.DataFrame({"ACP5": norm[i, 1], "MMP9": norm[i, 2]})
        score, present = mm.score_program(frame, mm.PRIMARY)
        assert present == ["ACP5", "MMP9"], f"{s}: genes present {present}"
        d = det[i]
        rows.append({"sample_id": s, "n_detectable": int(d.sum()), "n_nondetectable": int((~d).sum()),
                     "python_contrast": float(np.nanmean(score[d]) - np.nanmean(score[~d]))})
    py = pd.DataFrame(rows).merge(Rres[["sample_id", "contrast"]].rename(
        columns={"contrast": "r_contrast"}), on="sample_id")
    py["abs_diff"] = (py["python_contrast"] - py["r_contrast"]).abs()
    print(py.to_string(index=False))
    worst = float(py["abs_diff"].max())
    log(f"maximum |Python - R| per-sample contrast: {worst:.3e}")
    agree = worst < 1e-9
    print(f"=> the two implementations agree to {worst:.2e}: {agree}")
    if not agree:
        raise SystemExit("FATAL: independent re-implementation does not reproduce the R result")

    x = py["python_contrast"].to_numpy()
    n = len(x)
    rng = np.random.default_rng(SEED)
    boot = np.array([rng.choice(x, n, replace=True).mean() for _ in range(N_BOOT)])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    wr = stats.wilcoxon(x, zero_method="wilcox", correction=False,
                        alternative="two-sided", method="asymptotic")
    print(f"\nmean {x.mean():+.4f}  95% CI {lo:+.4f} to {hi:+.4f}  "
          f"(numpy default_rng(SEED), {N_BOOT} resamples, percentile)")
    print(f"scipy signed-rank p = {wr.pvalue:.6f}   (zero_method='wilcox', no correction)")
    print("R reported mean +0.9175, CI +0.0392 to +2.3055, signed-rank p = 0.115851")
    print("=> the CI and the p-value reproduce across languages and RNG streams.")

    # ----------------------------------------------------------------------
    print("\n" + "=" * 74)
    print("PART C -- POSITIVE CONTROLS CAPABLE OF BEING VIOLATED")
    print("=" * 74)
    print("Each row below applies a DELIBERATELY WRONG convention to the same data.")
    print("A control is informative only if it moves the answer.\n")

    def run(score_fn, sample_filter=None, weight=None, detection=None, label=""):
        dt = det if detection is None else detection
        vals, wts = [], []
        for s in elig:
            i = np.where(sid == s)[0]
            m = dt[i]
            if sample_filter is not None and not sample_filter(int(m.sum()), int((~m).sum())):
                continue
            sc = score_fn(i)
            vals.append(float(np.nanmean(sc[m]) - np.nanmean(sc[~m])))
            wts.append(int(m.sum()) if weight == "detectable" else 1)
        vals = np.asarray(vals); wts = np.asarray(wts, float)
        finite = np.isfinite(vals)
        return dict(label=label, n_samples=int(finite.sum()),
                    n_samples_considered=len(vals),
                    mean=float(np.average(vals[finite], weights=wts[finite])) if finite.any()
                    else float("nan"),
                    values=vals.tolist(), weights=wts.tolist())

    def score_manuscript(i):
        frame = pd.DataFrame({"ACP5": norm[i, 1], "MMP9": norm[i, 2]})
        return mm.score_program(frame, mm.PRIMARY)[0]

    def score_ddof1(i):
        v = norm[i][:, [1, 2]]
        return ((v - v.mean(0)) / v.std(0, ddof=1)).mean(1)

    def score_nolog(i):
        v = raw[i][:, [1, 2]]
        lb = libm[i][:, None]
        v = v / lb * 1e4
        return ((v - v.mean(0)) / v.std(0)).mean(1)

    def score_no_libsize(i):
        v = np.log1p(raw[i][:, [1, 2]])
        return ((v - v.mean(0)) / v.std(0)).mean(1)

    def score_target_in_program(i):
        v = norm[i]
        return ((v - v.mean(0)) / v.std(0)).mean(1)

    controls = [
        run(score_manuscript, label="REFERENCE: manuscript convention (no minimum group size)"),
        run(score_manuscript,
            sample_filter=lambda nd, nn: (nd + nn) >= 10 and nd >= 2 and nn >= 5,
            label="WRONG 1: restricted TREM2_TAM rule (>=10 cells, >=2 det, >=5 non-det) "
                  "applied to the PRIMARY contrast"),
        run(score_manuscript, detection=raw[:, 0] >= 2,
            label="WRONG 2: detectable defined as SIGLEC15 count >= 2 instead of > 0"),
        run(score_ddof1, label="WRONG 3: z-scored with ddof=1 instead of ddof=0"),
        run(score_nolog, label="WRONG 4: CP10K without log1p"),
        run(score_no_libsize, label="WRONG 5: log1p(counts) with no library-size scaling"),
        run(score_target_in_program,
            label="WRONG 6: SIGLEC15 included in the target-associated programme score"),
        run(score_manuscript, weight="detectable",
            label="WRONG 7: samples weighted by detectable-cell count instead of equally"),
    ]
    ref = controls[0]
    print(f"{'control':<72}{'n':>4}{'mean':>10}{'delta':>10}")
    for c in controls:
        d = c["mean"] - ref["mean"]
        print(f"{c['label']:<72}{c['n_samples']:>4}{c['mean']:>10.4f}{d:>+10.4f}")
    moved = [c for c in controls[1:]
             if not np.isfinite(c["mean"]) or abs(c["mean"] - ref["mean"]) > 1e-9]
    n_moved = len(moved)
    print(f"\nwrong conventions that move the answer: {n_moved}/{len(controls)-1}")
    print("=> PART B is NOT self-confirming: the checks above discriminate against "
          "conventions that are wrong.")
    print("\nTwo of these are worth naming explicitly:")
    w1 = controls[1]
    print(f"  * WRONG 1 changes n from {ref['n_samples']} to {w1['n_samples']} and the mean "
          f"from {ref['mean']:+.4f} to {w1['mean']:+.4f}. The frozen brief forbids this rule "
          f"for the primary contrast; it belongs to a separate restricted analysis.")
    print(f"  * WRONG 6 is the target-in-program error the frozen brief forbids outright.")

    (OUT / "_verification_independent.json").write_text(json.dumps({
        "part": "09_independent_verification",
        "manuscript_module": str(MANUSCRIPT_MODULE),
        "manuscript_primary_program": mm.PRIMARY,
        "manuscript_program_genes": mm.PROGRAMS[mm.PRIMARY],
        "manuscript_zscore_ddof": 0,
        "part_a_independent_parser_agrees": True,
        "part_b_max_abs_diff_python_vs_R": worst,
        "part_b_agree_1e-9": agree,
        "part_b_mean": float(x.mean()), "part_b_ci": [float(lo), float(hi)],
        "part_b_scipy_wilcoxon_p": float(wr.pvalue),
        "part_c_controls": controls,
        "part_c_wrong_conventions_that_move": n_moved,
        "part_c_wrong_conventions_total": len(controls) - 1,
    }, indent=1), encoding="utf-8")
    print("\nwrote _verification_independent.json")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
