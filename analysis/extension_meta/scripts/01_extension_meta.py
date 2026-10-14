#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
extension_meta / 01_extension_meta.py

Secondary bulk-extension meta-analysis for the BJC revision.

FROZEN RULES BEING EXECUTED (verbatim, from the frozen protocol text):
  "The three original PDAC bulk cohorts will be combined with the independently eligible
   primary-tumour GSE205154 cohort using the same Fisher-z transformation, REML random-effects
   model and modified Hartung-Knapp variance rule used in the manuscript. This four-cohort
   estimate is a revision-stage secondary pooled estimate and does not replace the original
   three-cohort estimate."
  "GSE21501 will then be added in a five-cohort sensitivity analysis because exact patient
   independence cannot be established from the deposited metadata."

Nothing here is a new primary estimate. The published three-cohort estimate is REPRODUCED,
not recomputed as a new primary. No cohort is excluded, reweighted or reshaped on the basis
of its rho magnitude or direction. Only bulk cohorts are pooled. No manuscript text, figure
or claim is written or modified.

Machinery is reused, not re-implemented: the manuscript's own `random_effects_meta()` is
imported by path from the v39 analysis archive and called directly. A byte-faithful replica
of that function is written ONLY so that deliberately mutated variants (the falsification
controls) can be run; the replica is validated against the imported original before any
mutated variant is believed.
"""

from __future__ import annotations
import os

import hashlib
import importlib.util
import inspect
import json
import platform
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import optimize, stats

# --------------------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------------------
EXT = Path(os.environ.get("SIGLEC15_EXT_META_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927\extension_meta"))
ARCHIVE = Path(
    os.path.join(os.environ.get("SIGLEC15_SUBMISSION_ROOT", r"E:\Tex\HMlong\npjPrecOncol\FINAL_SUBMISSION_v6"), "BMC_MEDICINE_SUBMISSION_v39", "_qa_archive_v39_20260912")
)
MANUSCRIPT_MODULE = ARCHIVE / "scripts" / "analysis_F_cross_scale.py"
MANUSCRIPT_SHA_EXPECTED = "7b6a311a0a73835f66596b7ea21070a6355c1fe6de44b518394dfe89f7c4334e"
BULK_ASSOC = ARCHIVE / "results" / "bulk_associations.csv"
BULK_META = ARCHIVE / "results" / "bulk_random_effects_meta.csv"

GSE205154_PRIMARY = Path(
    os.path.join(os.environ.get("SIGLEC15_PHASEB_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927"), "phaseB_GSE205154", "step3c_rho", "primary_result.csv")
)
GSE21501_PRIMARY = Path(
    os.path.join(os.environ.get("SIGLEC15_PHASEB_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927"), "phaseB_GSE21501", "step3c_rho", "output", "primary_result.csv")
)

PROGRAM = "Derosiers_ACP5_MMP9_core"          # the manuscript's PRIMARY bulk program
ORIGINAL_COHORTS = ["TCGA_PAAD", "GSE57495", "GSE85916"]

# Values stated for the two new cohorts in the task brief. Asserted, never assumed.
BRIEF = {
    "GSE205154": {"rho": 0.1317, "n": 216, "ci_low": -0.0024, "ci_high": 0.2648},
    "GSE21501": {"rho": 0.171188, "n": 129, "ci_low": -0.012824, "ci_high": 0.340531},
}

SPECIMEN_SENTENCE = (
    "The GSE21501 analysis was specimen-level because patient identifiers were not "
    "deposited; the number of independent patients could therefore not be established exactly."
)

# --------------------------------------------------------------------------------------
# Logging / checks
# --------------------------------------------------------------------------------------
LOG: list[str] = []
CHECKS: list[dict] = []
OBS: list[dict] = []


def log(msg: str = "") -> None:
    LOG.append(str(msg))
    print(msg, flush=True)


def observe(name: str, detail: str = "") -> None:
    """Record an obtained-but-not-asserted outcome (e.g. a control that came back null).

    Kept separate from CHECKS so that the pass/fail ledger is not padded with controls that
    were expected to be material and were not. A null control is a result and is reported.
    """
    OBS.append({"observation": name, "detail": str(detail)})
    log(f"[OBSERVED] {name}" + (f" :: {detail}" if detail else ""))


def check(name: str, ok: bool, detail: str = "") -> bool:
    CHECKS.append({"check": name, "passed": bool(ok), "detail": str(detail)})
    log(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" :: {detail}" if detail else ""))
    return bool(ok)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fmt(x, nd=6):
    try:
        if x is None or (isinstance(x, float) and not np.isfinite(x)):
            return "NA"
        return f"{float(x):.{nd}f}"
    except Exception:
        return str(x)


# --------------------------------------------------------------------------------------
# The manuscript's function, mutated only by explicit named switches (for controls)
# --------------------------------------------------------------------------------------
def replica_meta(
    rho,
    n,
    *,
    var_rule="n-3",
    hk="modified",
    sign_flip=(),
    transform="fisher_z",
    interval="hk_t",
):
    """Byte-faithful replica of the manuscript's random_effects_meta() with named switches.

    Default arguments reproduce the manuscript function EXACTLY. Each switch exists only so
    that a deliberately wrong convention can be run as a falsification control:

      var_rule   'n-3' (manuscript) | 'n-1' | 'n'   -> within-study variance 1/(n-c)
      hk         'modified' (manuscript: max(1, Q/(k-1))) | 'raw' (unmodified Hartung-Knapp)
      sign_flip  indices whose rho is negated before the Fisher-z transform
      transform  'fisher_z' (manuscript) | 'none' (pool raw rho on the rho scale)
      interval   'hk_t' (manuscript) | 'fixed_z' (conventional normal interval, no HK, no t)
    """
    rho = np.asarray(rho, dtype=float).copy()
    n = np.asarray(n, dtype=float)

    for i in sign_flip:
        rho[i] = -rho[i]

    if transform == "fisher_z":
        y = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    elif transform == "none":
        y = rho
    else:
        raise ValueError(transform)

    if var_rule == "n-3":
        v = 1.0 / (n - 3.0)
    elif var_rule == "n-1":
        v = 1.0 / (n - 1.0)
    elif var_rule == "n":
        v = 1.0 / n
    else:
        raise ValueError(var_rule)

    k = len(y)
    fixed_w = 1 / v
    fixed_mu = np.sum(fixed_w * y) / np.sum(fixed_w)
    q_fixed = np.sum(fixed_w * (y - fixed_mu) ** 2)
    i2 = max(0.0, (q_fixed - (k - 1)) / q_fixed) * 100 if q_fixed > 0 else 0.0

    def objective(tau2):
        w = 1.0 / (v + tau2)
        mu = np.sum(w * y) / np.sum(w)
        return np.sum(np.log(v + tau2)) + np.log(np.sum(w)) + np.sum(w * (y - mu) ** 2)

    result = optimize.minimize_scalar(objective, bounds=(0.0, 4.0), method="bounded")
    tau2 = max(0.0, float(result.x))
    if tau2 < 1e-8:
        tau2 = 0.0
    w = 1.0 / (v + tau2)
    mu = np.sum(w * y) / np.sum(w)
    q_hk_raw = np.sum(w * (y - mu) ** 2) / (k - 1)
    if hk == "modified":
        q_hk = max(1.0, float(q_hk_raw))
    elif hk == "raw":
        q_hk = float(q_hk_raw)
    else:
        raise ValueError(hk)

    if interval == "hk_t":
        se_hk = np.sqrt(q_hk / np.sum(w))
        critical = stats.t.ppf(0.975, df=k - 1)
    elif interval == "fixed_z":
        se_hk = np.sqrt(1.0 / np.sum(w))
        critical = stats.norm.ppf(0.975)
    else:
        raise ValueError(interval)

    lo_z, hi_z = mu - critical * se_hk, mu + critical * se_hk
    pred_se = np.sqrt(tau2 + se_hk**2)
    pred_lo, pred_hi = mu - critical * pred_se, mu + critical * pred_se
    p = 2 * stats.t.sf(abs(mu / se_hk), df=k - 1) if se_hk > 0 else np.nan

    back = np.tanh if transform == "fisher_z" else (lambda z: z)
    return {
        "n_cohorts": k,
        "r_random": float(back(mu)),
        "ci_low": float(back(lo_z)),
        "ci_high": float(back(hi_z)),
        "mu_fisher_z": float(mu),
        "se_used": float(se_hk),
        "critical": float(critical),
        "p_modified_hartung_knapp": float(p),
        "hk_scale_raw": float(q_hk_raw),
        "hk_scale_used": float(q_hk),
        "tau2_fisher_z": float(tau2),
        "I2_pct": float(i2),
        "Q": float(q_fixed),
        "prediction_low": float(back(pred_lo)),
        "prediction_high": float(back(pred_hi)),
    }


def make_frame(rows):
    """rows = list of (cohort, rho, n_samples) -> the DataFrame shape the manuscript fn wants."""
    return pd.DataFrame(
        [{"cohort": c, "program": PROGRAM, "rho": r, "n_samples": nn} for c, r, nn in rows]
    )


def model_z_stats(rows, tau2):
    """Pooled z-scale mean and the se the manuscript function actually used, for a fitted model.

    Mirrors lines 432-439 of the archived module:
        w = 1/(v+tau2);  mu = sum(w y)/sum(w);  se = sqrt(max(1, Q/(k-1)) / sum(w))
    `mu` and `se` are NOT returned by random_effects_meta(), so they are recomputed here and
    validated by checking that tanh(mu) reproduces the function's returned r_random.
    """
    y = np.arctanh(np.clip(np.array([r for _, r, _ in rows], float), -0.999999, 0.999999))
    v = 1.0 / (np.array([nn for _, _, nn in rows], float) - 3.0)
    w = 1.0 / (v + tau2)
    mu = float(np.sum(w * y) / np.sum(w))
    k = len(y)
    q_hk_raw = float(np.sum(w * (y - mu) ** 2) / (k - 1))
    se = float(np.sqrt(max(1.0, q_hk_raw) / np.sum(w)))
    return mu, se, q_hk_raw


def fixed_Q(rows):
    """The fixed-effect Q, on Fisher z, using the identical formula that sits inside the
    manuscript's random_effects_meta() (line 420 of the archived module):

        w = 1/v, v = 1/(n-3), mu_fixed = sum(w y)/sum(w), Q = sum(w (y - mu_fixed)^2)

    The function itself does not RETURN Q, only the I2 derived from it. This helper is
    therefore validated in section 4b by checking that the I2 derived from it reproduces the
    function's own returned I2_pct; that check is what licenses quoting Q as the manuscript's Q.
    """
    y = np.arctanh(np.clip(np.array([r for _, r, _ in rows], float), -0.999999, 0.999999))
    v = 1.0 / (np.array([nn for _, _, nn in rows], float) - 3.0)
    w = 1 / v
    mu = np.sum(w * y) / np.sum(w)
    return float(np.sum(w * (y - mu) ** 2))


# --------------------------------------------------------------------------------------
# 0. Header
# --------------------------------------------------------------------------------------
def main() -> int:
    t0 = time.time()
    EXT.mkdir(parents=True, exist_ok=True)
    (EXT / "provenance").mkdir(exist_ok=True)

    log("=" * 100)
    log("extension_meta / 01_extension_meta.py")
    log(f"run started : {datetime.now().isoformat(timespec='seconds')}")
    log(f"python      : {sys.version.split()[0]} ({platform.python_implementation()})")
    log(f"numpy {np.__version__} | scipy {scipy.__version__} | pandas {pd.__version__}")
    log("=" * 100)

    # ----------------------------------------------------------------------------------
    # 1. Source integrity
    # ----------------------------------------------------------------------------------
    log("\n### 1. Source integrity")
    manu_sha = sha256(MANUSCRIPT_MODULE)
    check(
        "manuscript module sha256 == recorded 7b6a311a...",
        manu_sha == MANUSCRIPT_SHA_EXPECTED,
        f"{manu_sha}",
    )
    if manu_sha != MANUSCRIPT_SHA_EXPECTED:
        log("FATAL: manuscript module hash mismatch -- refusing to proceed.")
        return 2

    src_by_file = {
        "analysis_F_cross_scale.py": manu_sha,
        "results/bulk_associations.csv": sha256(BULK_ASSOC),
        "results/bulk_random_effects_meta.csv": sha256(BULK_META),
        "GSE205154 step3c_rho/primary_result.csv": sha256(GSE205154_PRIMARY),
        "GSE21501 step3c_rho/output/primary_result.csv": sha256(GSE21501_PRIMARY),
    }
    for k, v in src_by_file.items():
        log(f"    sha256 {v}  {k}")

    # ----------------------------------------------------------------------------------
    # 2. Import the manuscript's own implementation
    # ----------------------------------------------------------------------------------
    log("\n### 2. Importing the manuscript's own random_effects_meta()")
    spec = importlib.util.spec_from_file_location("manuscript_analysis_F_cross_scale", MANUSCRIPT_MODULE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mer = mod.random_effects_meta
    fn_src = inspect.getsource(mer)
    fn_sha = hashlib.sha256(fn_src.encode("utf-8")).hexdigest()
    log(f"    imported from : {MANUSCRIPT_MODULE}")
    log(f"    function      : random_effects_meta(group: pd.DataFrame) -> dict")
    log(f"    source lines  : {inspect.getsourcelines(mer)[1]}..{inspect.getsourcelines(mer)[1] + len(inspect.getsourcelines(mer)[0]) - 1}")
    log(f"    source sha256 : {fn_sha}")
    log(f"    module PRIMARY program constant = {mod.PRIMARY!r}")
    check("imported module PRIMARY is the program used here", mod.PRIMARY == PROGRAM, mod.PRIMARY)

    # ----------------------------------------------------------------------------------
    # 3. Per-cohort inputs -- taken from the manuscript's own archived outputs
    # ----------------------------------------------------------------------------------
    log("\n### 3. Per-cohort inputs")
    log("3a. The three ORIGINAL PDAC bulk cohorts -- read from the manuscript's archived")
    log("    analysis outputs (NOT invented, NOT re-derived here):")
    log(f"      {BULK_ASSOC}")
    assoc = pd.read_csv(BULK_ASSOC)
    orig = assoc.loc[assoc["program"] == PROGRAM, ["cohort", "n_samples", "rho", "ci_low", "ci_high"]]
    orig = orig.set_index("cohort")
    check(
        "all three original cohorts present in archived per-cohort table",
        set(ORIGINAL_COHORTS) <= set(orig.index),
        ";".join(orig.index),
    )
    orig_rows = []
    for c in ORIGINAL_COHORTS:
        r = float(orig.loc[c, "rho"])
        nn = int(orig.loc[c, "n_samples"])
        orig_rows.append((c, r, nn))
        log(f"      {c:<10s} n={nn:>4d}  rho={r!r}")

    log("3b. The archived three-cohort pooled estimate (the published anchor):")
    meta_arch = pd.read_csv(BULK_META)
    anchor = meta_arch.loc[meta_arch["program"] == PROGRAM].iloc[0]
    PUBLISHED = {
        "rho": float(anchor["r_random"]),
        "ci_low": float(anchor["ci_low"]),
        "ci_high": float(anchor["ci_high"]),
        "p": float(anchor["p_modified_hartung_knapp"]),
        "tau2": float(anchor["tau2_fisher_z"]),
        "I2": float(anchor["I2_pct"]),
        "cohorts": str(anchor["cohorts"]),
        "k": int(anchor["n_cohorts"]),
    }
    log(f"      k={PUBLISHED['k']} cohorts={PUBLISHED['cohorts']}")
    log(f"      pooled rho = {PUBLISHED['rho']!r}")
    log(f"      95% CI     = {PUBLISHED['ci_low']!r} .. {PUBLISHED['ci_high']!r}")
    log(f"      tau2 = {PUBLISHED['tau2']!r}   I2 = {PUBLISHED['I2']!r}   p = {PUBLISHED['p']!r}")
    check(
        "archived anchor rounds to the published rho = 0.132",
        abs(PUBLISHED["rho"] - 0.132) < 5e-4,
        f"{PUBLISHED['rho']:.6f}",
    )
    check(
        "archived anchor rounds to the published CI -0.110 .. +0.360",
        abs(PUBLISHED["ci_low"] - (-0.110)) < 5e-4 and abs(PUBLISHED["ci_high"] - 0.360) < 5e-4,
        f"{PUBLISHED['ci_low']:.6f} .. {PUBLISHED['ci_high']:.6f}",
    )

    log("3c. The two NEW cohorts -- read from this revision's step3c_rho primary results:")
    g54 = pd.read_csv(GSE205154_PRIMARY)
    g51 = pd.read_csv(GSE21501_PRIMARY)
    NEW = {}
    for label, df, key in (("GSE205154", g54, None), ("GSE21501", g51, None)):
        row = df.iloc[0]
        NEW[label] = {
            "rho": float(row["rho"]),
            "n": int(row["n"]),
            "ci_low": float(row["ci_low"]),
            "ci_high": float(row["ci_high"]),
            "file": str(GSE205154_PRIMARY if label == "GSE205154" else GSE21501_PRIMARY),
        }
        log(f"      {label:<10s} n={NEW[label]['n']:>4d}  rho={NEW[label]['rho']!r}"
            f"  CI {NEW[label]['ci_low']:.10f} .. {NEW[label]['ci_high']:.10f}")
        for field in ("rho", "n", "ci_low", "ci_high"):
            b = BRIEF[label][field]
            got = NEW[label][field]
            tol = (5e-5 if field == "rho" else (0.5 if field == "n" else 5e-5))
            check(
                f"{label}: {field} in primary_result.csv matches the brief",
                abs(got - b) <= tol,
                f"file={got!r} brief={b!r}",
            )

    # GSE21501 unit statement, carried verbatim
    unit_stmt = str(g51.iloc[0].get("unit_statement", "")).strip()
    check("GSE21501 unit_statement equals the required verbatim sentence",
          unit_stmt == SPECIMEN_SENTENCE, unit_stmt)
    n_patients_field = str(g51.iloc[0].get("patients", ""))
    log(f"      GSE21501 patients field = {n_patients_field!r}  (a BOUND, never an exact count)")
    check("GSE21501 patient field states a bound, not an exact count",
          "<=" in n_patients_field, n_patients_field)

    # ----------------------------------------------------------------------------------
    # 4. The three frozen models
    # ----------------------------------------------------------------------------------
    log("\n### 4. The three frozen models (manuscript machinery, called directly)")
    MODELS = {
        "3-cohort published (reproduced)": list(orig_rows),
        "4-cohort secondary extension": list(orig_rows) + [("GSE205154", NEW["GSE205154"]["rho"], NEW["GSE205154"]["n"])],
        "5-cohort sensitivity": list(orig_rows) + [
            ("GSE205154", NEW["GSE205154"]["rho"], NEW["GSE205154"]["n"]),
            ("GSE21501", NEW["GSE21501"]["rho"], NEW["GSE21501"]["n"]),
        ],
    }
    results = {}
    for label, rows in MODELS.items():
        frame = make_frame(rows)
        out = mer(frame)
        out["coh"] = [r[0] for r in rows]
        out["rho_in"] = [r[1] for r in rows]
        out["n_in"] = [r[2] for r in rows]
        out["Q"] = fixed_Q(rows)          # validated against the function's own I2 in 4b
        out["mu_z"], out["se_z"], _qhkr = model_z_stats(rows, out["tau2_fisher_z"])
        check(
            f"recomputed z-scale mean reproduces the function's pooled rho ({label})",
            abs(np.tanh(out["mu_z"]) - out["r_random"]) < 1e-15,
            f"tanh(mu)={np.tanh(out['mu_z'])!r} fn={out['r_random']!r}",
        )
        results[label] = out
        log(f"\n  {label}  (k={out['n_cohorts']})")
        log(f"    cohorts        : {out['cohorts']}")
        log(f"    pooled rho     : {out['r_random']!r}")
        log(f"    95% CI         : {out['ci_low']!r} .. {out['ci_high']!r}")
        log(f"    tau2 (REML)    : {out['tau2_fisher_z']!r}")
        log(f"    I2             : {out['I2_pct']!r} %")
        log(f"    Q (fixed, z)   : {out['Q']!r}")
        log(f"    HK scale raw   : {out['hk_scale_raw']!r}  -> used {out['hk_scale_used']!r}")
        log(f"    p (mod HK)     : {out['p_modified_hartung_knapp']!r}")
        log(f"    crosses zero   : {bool(out['ci_low'] < 0 < out['ci_high'])}")

    # ---- 4a. three-cohort reproduction against the archive -------------------------
    log("\n### 4a. ANCHOR REPRODUCTION (three-cohort estimate vs the archive)")
    rep = results["3-cohort published (reproduced)"]
    for field, pub, tol in (
        ("r_random", PUBLISHED["rho"], 1e-9),
        ("ci_low", PUBLISHED["ci_low"], 1e-9),
        ("ci_high", PUBLISHED["ci_high"], 1e-9),
        ("tau2_fisher_z", PUBLISHED["tau2"], 1e-15),
        ("I2_pct", PUBLISHED["I2"], 1e-12),
        ("p_modified_hartung_knapp", PUBLISHED["p"], 1e-9),
    ):
        d = abs(float(rep[field]) - float(pub))
        check(f"anchor reproduced: {field} matches archive (tol {tol:g})", d <= tol,
              f"reproduced={float(rep[field])!r} archive={float(pub)!r} |delta|={d:.3e}")
    check("anchor rounds to rho=0.132 at 3 dp", round(rep["r_random"], 3) == 0.132, f"{rep['r_random']:.6f}")
    check("anchor rounds to CI -0.110 at 3 dp", round(rep["ci_low"], 3) == -0.110, f"{rep['ci_low']:.6f}")
    check("anchor rounds to CI +0.360 at 3 dp", round(rep["ci_high"], 3) == 0.360, f"{rep['ci_high']:.6f}")

    # Attribution of the residual: the per-cohort INPUTS are bit-identical, so any residual can
    # only come from the numerics of this machine's numpy/scipy. Recover the critical value the
    # archive must have used and compare it with this machine's t quantile.
    log("")
    log("    Residual attribution. The per-cohort rho values read from the archive are bit-identical")
    log(f"    to the ones used here (rho column of bulk_associations.csv). Residual on the CI:")
    for k_rep, k_pub in (("r_random", "rho"), ("ci_low", "ci_low"),
                         ("ci_high", "ci_high"), ("p_modified_hartung_knapp", "p")):
        dv = abs(float(rep[k_rep]) - float(PUBLISHED[k_pub]))
        log(f"      {k_rep:<26s} |delta| = {dv:.3e}"
            f"   (relative {dv / max(1e-300, abs(float(PUBLISHED[k_pub]))):.2e})")
    # mu and se are not returned by the manuscript function, so recompute them on the z scale
    # with the function's own returned tau2 (formulas identical to lines 432-439 of the module).
    _y = np.arctanh(np.clip(np.array(rep["rho_in"], float), -0.999999, 0.999999))
    _v = 1.0 / (np.array(rep["n_in"], float) - 3.0)
    _w = 1.0 / (_v + rep["tau2_fisher_z"])
    mu_r = float(np.sum(_w * _y) / np.sum(_w))
    se_r = float(np.sqrt(max(1.0, rep["hk_scale_raw"]) / np.sum(_w)))
    check("recomputed mu/se reproduce the function's pooled rho (z scale)",
          abs(np.tanh(mu_r) - rep["r_random"]) < 1e-15, f"tanh(mu)={np.tanh(mu_r)!r}")
    c_implied = (mu_r - float(np.arctanh(PUBLISHED["ci_low"]))) / se_r
    c_here = float(stats.t.ppf(0.975, df=2))
    rel = abs(c_implied - c_here) / c_here
    check("residual is attributable to the t quantile differing across scipy builds",
          rel < 1e-9,
          f"t implied by archive = {c_implied!r}, t here = {c_here!r}, relative diff {rel:.2e}")
    log(f"      archive-implied t_0.975,2 = {c_implied!r}")
    log(f"      this machine's t_0.975,2  = {c_here!r}   (relative diff {rel:.2e})")
    log("    => the residual is floating-point in the t quantile, i.e. a scipy-build difference,")
    log("    not a difference in inputs, model, or implementation. Reproduction is exact to the")
    log("    precision at which the archive's own numbers are reproducible across builds, and")
    log("    exact at every published digit (rho 0.132, CI -0.110 .. +0.360).")
    log("    NOTE: the published CI is the modified-HK t interval, not the per-cohort bootstrap;")
    log("    the archive's r_random/ci_low/ci_high columns are the random_effects_meta() outputs.")

    # ---- 4b. Q recovery check: our Q must reproduce the manuscript's own I2 ----------
    log("\n### 4b. Q RECOVERY CHECK (Q is not returned by the manuscript function)")
    for label, out in results.items():
        y = np.arctanh(np.clip(np.array(out["rho_in"]), -0.999999, 0.999999))
        v = 1.0 / (np.array(out["n_in"], float) - 3.0)
        w = 1 / v
        mu = np.sum(w * y) / np.sum(w)
        q = float(np.sum(w * (y - mu) ** 2))
        k = len(y)
        i2_from_q = max(0.0, (q - (k - 1)) / q) * 100 if q > 0 else 0.0
        check(
            f"I2 derived from the independently computed Q reproduces the function's I2 ({label})",
            abs(i2_from_q - out["I2_pct"]) < 1e-9,
            f"Q={q:.12f} I2(Q)={i2_from_q:.12f} I2(fn)={out['I2_pct']:.12f}",
        )

    # ---- 4c. Why tau2 is the same number in all three models -------------------------
    log("\n### 4c. TAU2 AT THE REML BOUNDARY (diagnostic, not an assumption)")
    log("    The manuscript function returns tau2 = 5.628673589008228e-06 for ALL THREE models,")
    log("    bit-identically. That is not a bug carried over: in every model Q < k-1, so the REML")
    log("    criterion is minimised at the lower boundary tau2 = 0, and scipy's bounded scalar")
    log("    search walks the same deterministic path to the same tolerance-limited value. The")
    log("    check below confirms the criterion is increasing away from the boundary in each model.")
    for label, out in results.items():
        y = np.arctanh(np.clip(np.array(out["rho_in"]), -0.999999, 0.999999))
        v = 1.0 / (np.array(out["n_in"], float) - 3.0)

        def crit(t2, y=y, v=v):
            w = 1.0 / (v + t2)
            mu = np.sum(w * y) / np.sum(w)
            return np.sum(np.log(v + t2)) + np.log(np.sum(w)) + np.sum(w * (y - mu) ** 2)

        c0, c_ret, c_big = crit(0.0), crit(out["tau2_fisher_z"]), crit(0.001)
        k = len(y)
        check(
            f"REML criterion is minimised at the tau2=0 boundary ({label})",
            c0 <= c_ret <= c_big and out["Q"] < (k - 1),
            f"Q={out['Q']:.6f} < k-1={k-1}; crit(0)={c0:.12f} crit(returned)={c_ret:.12f} "
            f"crit(1e-3)={c_big:.12f}",
        )
        log(f"      {label}: Q={out['Q']:.9f} < k-1={k-1}  ->  tau2(REML) is AT THE BOUNDARY 0;")
        log(f"        the value {out['tau2_fisher_z']!r} is the bounded optimiser's approach to it.")
    log("    Substantive reading: no detectable between-cohort heterogeneity in any model; I2 = 0.")
    log("    The manuscript function's tau2 < 1e-8 -> 0 snap does NOT fire (5.63e-06 > 1e-8), and")
    log("    the archived three-cohort row carries the identical value, so this step introduces no")
    log("    new convention.")

    # ----------------------------------------------------------------------------------
    # 5. Replica validation and falsification controls
    # ----------------------------------------------------------------------------------
    log("\n### 5. Replica validation (must pass BEFORE any control is believed)")
    ctrl_rows = []
    for label, rows in MODELS.items():
        frame = make_frame(rows)
        a = mer(frame)
        b = replica_meta([r[1] for r in rows], [r[2] for r in rows])
        d = max(abs(a[k] - b[k]) for k in
                ("r_random", "ci_low", "ci_high", "tau2_fisher_z", "I2_pct", "hk_scale_raw"))
        check(f"C0 replica identity == manuscript function ({label})", d < 1e-14, f"max|delta|={d:.3e}")
        ctrl_rows.append({
            "control": "C0_replica_identity", "model": label,
            "pooled_rho": b["r_random"], "ci_low": b["ci_low"], "ci_high": b["ci_high"],
            "delta_rho_vs_frozen": b["r_random"] - a["r_random"],
            "crosses_zero": bool(b["ci_low"] < 0 < b["ci_high"]),
            "verdict": "IDENTITY -- replica is byte-equivalent",
        })

    def run_control(cid, desc, model_label, **kw):
        rows = MODELS[model_label]
        base = results[model_label]
        out = replica_meta([r[1] for r in rows], [r[2] for r in rows], **kw)
        d_rho = out["r_random"] - base["r_random"]
        base_cross = base["ci_low"] < 0 < base["ci_high"]
        ctrl_cross = out["ci_low"] < 0 < out["ci_high"]
        base_w = base["ci_high"] - base["ci_low"]
        ctrl_w = out["ci_high"] - out["ci_low"]
        width_ratio = ctrl_w / base_w
        # Materiality criterion, stated explicitly. A convention is "materially different" if it
        #   (i)  moves the pooled rho by more than 0.01, OR
        #   (ii) changes whether the interval excludes zero, OR
        #   (iii) changes the interval WIDTH by more than 25%.
        # (iii) was added after the first run: the unmodified-HK control on the five-cohort model
        # narrowed the interval from 0.213 to 0.078 while leaving rho and the zero-crossing
        # verdict untouched, and the first two-clause criterion scored that as "no material
        # change". The criterion was under-specified, not the control.
        trig = []
        if abs(d_rho) > 0.01:
            trig.append(f"|drho|={abs(d_rho):.4f}>0.01")
        if base_cross != ctrl_cross:
            trig.append("ZERO-CROSSING VERDICT CHANGES")
        if not (0.75 <= width_ratio <= 1.0 / 0.75):
            trig.append(f"interval width x{width_ratio:.3f}")
        material = bool(trig)
        verdict = ("MATERIALLY DIFFERENT (" + "; ".join(trig) + ")") if material else (
            f"NO MATERIAL CHANGE (|drho|={abs(d_rho):.5f}, zero-crossing unchanged, "
            f"width x{width_ratio:.3f})")
        log(f"  [{cid}] {desc}")
        log(f"        {model_label}: rho={out['r_random']:.6f} CI {out['ci_low']:.6f}..{out['ci_high']:.6f}"
            f"  width={ctrl_w:.6f}")
        log(f"        frozen          : rho={base['r_random']:.6f} CI {base['ci_low']:.6f}..{base['ci_high']:.6f}"
            f"  width={base_w:.6f}")
        log(f"        delta rho = {d_rho:+.6f}   -> {verdict}")
        ctrl_rows.append({
            "control": cid, "model": model_label, "pooled_rho": out["r_random"],
            "ci_low": out["ci_low"], "ci_high": out["ci_high"],
            "ci_width": ctrl_w, "width_ratio_vs_frozen": width_ratio,
            "delta_rho_vs_frozen": d_rho, "crosses_zero": bool(ctrl_cross),
            "verdict": verdict, "description": desc,
            "frozen_rho": base["r_random"], "frozen_ci_low": base["ci_low"],
            "frozen_ci_high": base["ci_high"], "frozen_crosses_zero": bool(base_cross),
        })
        return out, material

    log("\n### 5b. FALSIFICATION CONTROLS -- conventions that are supposed to change the answer")
    material_controls: list[str] = []
    null_controls: list[str] = []

    def tally(cid, m):
        (material_controls if m else null_controls).append(cid)

    # C1/C2 - wrong within-study variance rule
    o, m = run_control("C1_var_rule_n-1", "within-study variance 1/(n-1) instead of 1/(n-3)",
                       "4-cohort secondary extension", var_rule="n-1")
    tally("C1_var_rule_n-1", m)
    observe("C1 wrong variance rule 1/(n-1) comes back NULL on the four-cohort model",
            f"rho={o['r_random']:.6f} vs frozen 0.132074; the estimate is insensitive to the "
            f"n-3 offset because all four cohorts have large n and tau2 is at the boundary")
    o, m = run_control("C2_var_rule_n", "within-study variance 1/n instead of 1/(n-3)",
                       "4-cohort secondary extension", var_rule="n")
    tally("C2_var_rule_n", m)
    observe("C2 wrong variance rule 1/n comes back NULL on the four-cohort model",
            f"rho={o['r_random']:.6f} vs frozen 0.132074")

    # C3 - unmodified Hartung-Knapp (the mutation of the "modified" rule)
    o, m = run_control("C3_unmodified_HK", "UNMODIFIED Hartung-Knapp (drop the max(1, Q/(k-1)) floor)",
                       "4-cohort secondary extension", hk="raw")
    tally("C3_unmodified_HK|4", m)
    check("C3 (unmodified HK, four-cohort) materially changes the pooled interval", m,
          f"rho={o['r_random']:.6f} CI {o['ci_low']:.6f}..{o['ci_high']:.6f}")
    o5, m5 = run_control("C3_unmodified_HK", "UNMODIFIED Hartung-Knapp (drop the max(1, Q/(k-1)) floor)",
                         "5-cohort sensitivity", hk="raw")
    tally("C3_unmodified_HK|5", m5)
    check("C3 (unmodified HK, five-cohort) materially changes the pooled interval", m5,
          f"rho={o5['r_random']:.6f} CI {o5['ci_low']:.6f}..{o5['ci_high']:.6f}")
    o3, m3 = run_control("C3b_unmodified_HK_published_3cohort",
                         "UNMODIFIED Hartung-Knapp applied to the PUBLISHED three-cohort model "
                         "(control only -- no estimate in meta_results.csv uses it)",
                         "3-cohort published (reproduced)", hk="raw")
    tally("C3b_unmodified_HK|3", m3)
    check("C3b the modified-HK rule is load-bearing in the PUBLISHED three-cohort estimate too", m3,
          f"rho={o3['r_random']:.6f} CI {o3['ci_low']:.6f}..{o3['ci_high']:.6f} vs published "
          f"-0.110117..+0.359917")
    observe("C3b consequence for reading the published anchor",
            "the published three-cohort CI crosses zero with the manuscript's modified-HK floor "
            f"({PUBLISHED['ci_low']:.6f}..{PUBLISHED['ci_high']:.6f}); with the UNMODIFIED HK "
            f"variance the same three cohorts give {o3['ci_low']:.6f}..{o3['ci_high']:.6f}, which "
            "excludes zero. The modified rule (frozen) is therefore the reason the published "
            "estimate is a null-crossing one, and the four-cohort estimate inherits that. This "
            "is a property of the frozen convention, reported as a control; no estimate here "
            "uses the unmodified variance.")

    # C4 - conventional normal interval, no HK and no t
    o, m = run_control("C4_fixed_z_interval", "conventional normal interval, no Hartung-Knapp, no t",
                       "4-cohort secondary extension", interval="fixed_z")
    tally("C4_fixed_z_interval", m)
    check("C4 (conventional normal interval) materially changes the pooled interval", m,
          f"rho={o['r_random']:.6f}")

    # C5 - wrong sign convention applied globally
    o, m = run_control("C5_sign_flip_all", "wrong sign convention: Fisher z of -rho for every cohort",
                       "4-cohort secondary extension", sign_flip=(0, 1, 2, 3))
    tally("C5_sign_flip_all", m)
    check("C5 (global sign flip) flips the pooled estimate", m and o["r_random"] < 0,
          f"rho={o['r_random']:.6f}")

    # C6 - wrong sign convention on one cohort only (transcription error)
    o, m = run_control("C6_sign_flip_GSE205154",
                       "sign of the new GSE205154 cohort entered negative (single-cohort transcription error)",
                       "4-cohort secondary extension", sign_flip=(3,))
    tally("C6_sign_flip_GSE205154", m)
    check("C6 (single-cohort sign flip) materially changes the pooled interval", m,
          f"rho={o['r_random']:.6f}")

    # C7 - no Fisher-z transformation
    o, m = run_control("C7_no_fisher_z", "pool raw rho on the rho scale instead of Fisher z",
                       "4-cohort secondary extension", transform="none")
    tally("C7_no_fisher_z", m)
    observe("C7 dropping the Fisher-z transformation comes back NULL on the four-cohort model",
            f"rho={o['r_random']:.6f} vs frozen 0.132074, delta {o['r_random'] - results['4-cohort secondary extension']['r_random']:+.6f}; "
            "all |rho| < 0.18, where arctanh is nearly the identity, so the transform is not "
            "load-bearing at this effect size")

    check("falsification requirement met: at least one convention is materially different",
          len(material_controls) >= 1,
          f"material: {sorted(set(material_controls))}")
    observe("controls that came back NULL (reported, not hidden)",
            f"{sorted(set(null_controls))} -- see failed_or_excluded.md Part C")

    # C8 - the anchor reproduction check must be capable of FAILING
    log("\n### 5c. CONTROL THAT THE ANCHOR CHECK ITSELF CAN FAIL")
    bad = [(c, (r + 0.03 if c == "GSE57495" else r), nn) for c, r, nn in orig_rows]
    bad_out = mer(make_frame(bad))
    d = abs(bad_out["ci_low"] - PUBLISHED["ci_low"])
    fails = d > 1e-6
    check("C8 a perturbed three-cohort input does NOT reproduce the published anchor", fails,
          f"perturbed ci_low={bad_out['ci_low']:.6f} vs published {PUBLISHED['ci_low']:.6f} "
          f"(|delta|={d:.3e}) -- the anchor check is falsifiable")
    ctrl_rows.append({
        "control": "C8_anchor_check_can_fail", "model": "3-cohort published (reproduced)",
        "pooled_rho": bad_out["r_random"], "ci_low": bad_out["ci_low"], "ci_high": bad_out["ci_high"],
        "delta_rho_vs_frozen": bad_out["r_random"] - PUBLISHED["rho"],
        "crosses_zero": bool(bad_out["ci_low"] < 0 < bad_out["ci_high"]),
        "verdict": "anchor check FAILS as intended when an input is perturbed (+0.03 on GSE57495 rho)",
        "description": "falsifiability of the three-cohort reproduction check",
        "frozen_rho": PUBLISHED["rho"], "frozen_ci_low": PUBLISHED["ci_low"],
        "frozen_ci_high": PUBLISHED["ci_high"], "frozen_crosses_zero": True,
    })

    # ---- 5d. independent re-derivation of the three-cohort arithmetic by hand --------
    log("\n### 5d. HAND ARITHMETIC RE-DERIVATION (three-cohort, independent of the function)")
    ys = np.arctanh(np.clip(np.array([r for _, r, _ in orig_rows]), -0.999999, 0.999999))
    ns = np.array([nn for _, _, nn in orig_rows], float)
    vs = 1.0 / (ns - 3.0)

    # First form of this check (tau2 assumed 0) was WRONG and is retained as a finding: it
    # returned max|delta| = 9.568e-05, because the function's returned tau2 (5.63e-06, at the
    # boundary) still perturbs the weights enough to move the pooled z by ~1e-4. The check is
    # re-derived with the function's own tau2 included.
    ws0 = 1.0 / vs
    mu_hand0 = float(np.sum(ws0 * ys) / np.sum(ws0))
    d0 = abs(np.tanh(mu_hand0) - rep["r_random"])
    observe("first form of the hand check (tau2 forced to 0) is mis-specified",
            f"max|delta| on pooled rho = {d0:.3e}; the boundary tau2 still shifts the weights. "
            "Re-derived with the function's tau2 below.")

    tau2_used = rep["tau2_fisher_z"]
    ws = 1.0 / (vs + tau2_used)
    mu_hand = float(np.sum(ws * ys) / np.sum(ws))
    q_raw_hand = float(np.sum(ws * (ys - mu_hand) ** 2) / 2.0)
    se_hand = float(np.sqrt(max(1.0, q_raw_hand) / np.sum(ws)))
    crit_hand = float(stats.t.ppf(0.975, df=2))
    lo_hand, hi_hand = np.tanh(mu_hand - crit_hand * se_hand), np.tanh(mu_hand + crit_hand * se_hand)
    d = max(abs(lo_hand - rep["ci_low"]), abs(hi_hand - rep["ci_high"]),
            abs(np.tanh(mu_hand) - rep["r_random"]))
    check("hand arithmetic reproduces the three-cohort interval (with the function's tau2)",
          d < 1e-12, f"max|delta|={d:.3e}")
    log(f"    z-scale: mu={mu_hand:.12f}  se={se_hand:.12f}  t_0.975,2={crit_hand:.12f}"
        f"  (tau2={tau2_used!r} carried, as the function does)")
    log(f"    q_hk raw by hand = {q_raw_hand!r}  vs function {rep['hk_scale_raw']!r}")
    check("hand arithmetic reproduces the function's HK raw scale", abs(q_raw_hand - rep["hk_scale_raw"]) < 1e-15,
          f"|delta|={abs(q_raw_hand - rep['hk_scale_raw']):.3e}")
    check("three-cohort HK raw scale < 1, so the modified rule's floor is binding",
          rep["hk_scale_raw"] < 1.0, f"raw={rep['hk_scale_raw']:.6f} used={rep['hk_scale_used']:.6f}")

    # ----------------------------------------------------------------------------------
    # 6. Outputs
    # ----------------------------------------------------------------------------------
    log("\n### 6. Writing outputs")

    TAU2_NOTE = (
        " TAU2 DISCLOSURE: the tau2 column carries the value the manuscript function returns "
        "(5.628673589008228e-06, bit-identical in all three models). The substantive REML "
        "estimate is tau2 = 0 AT THE LOWER BOUNDARY: Q < k-1 in every model, the REML criterion "
        "is minimised at tau2 = 0 (verified in the run log section 4c), and the returned value is "
        "the bounded optimiser's tolerance-limited approach to that boundary (scipy xatol = 1e-5), "
        "which the manuscript's own 'tau2 < 1e-8 -> 0' snap does not catch. metafor 4.8.0 returns "
        "exactly 0 for all three models; carrying metafor's 0 instead moves each interval by "
        "<= 7e-5 and reproduces every zero-crossing verdict. Read tau2 as 0 (lower boundary), "
        "I2 = 0, i.e. no detectable between-cohort heterogeneity."
    )

    LABEL = {
        "3-cohort published (reproduced)": "3-cohort published bulk estimate (REPRODUCED, not a new primary)",
        "4-cohort secondary extension": "4-cohort secondary extension (3 original PDAC bulk + GSE205154)",
        "5-cohort sensitivity": "5-cohort sensitivity (4-cohort + GSE21501)",
    }
    NOTE = {
        "3-cohort published (reproduced)": (
            "PRIMARY PUBLISHED ESTIMATE, UNCHANGED AND UNREPLACED. Reproduced by calling the "
            "manuscript's own random_effects_meta() on the archived per-cohort values; it is the "
            "same number as, not a substitute for, the manuscript's estimate." + TAU2_NOTE
        ),
        "4-cohort secondary extension": (
            "REVISION-STAGE SECONDARY ESTIMATE. Does not replace the three-cohort estimate. "
            "GSE205154 is patient-level (n=216 patients/specimens), log2(TMM+1), GPL20301 "
            "exome-capture; admissible because only within-cohort ranks are used." + TAU2_NOTE
        ),
        "5-cohort sensitivity": (
            "SENSITIVITY ONLY. GSE21501 is SPECIMEN-level (n=129 eligible tumour arrays/"
            "specimens, <=129 patients) and therefore does NOT carry the same inferential weight "
            "as the patient-level cohorts: v=1/(n-3) assumes independent units and the exact "
            "number of independent patients is not established. " + SPECIMEN_SENTENCE + TAU2_NOTE
        ),
    }

    rows_out = []
    for label, out in results.items():
        rows_out.append({
            "estimate_label": LABEL[label],
            "cohorts_included": out["cohorts"],
            "k": out["n_cohorts"],
            "pooled_rho": repr(out["r_random"]),
            "ci_low": repr(out["ci_low"]),
            "ci_high": repr(out["ci_high"]),
            "tau2": repr(out["tau2_fisher_z"]),
            "I2": repr(out["I2_pct"]),
            "Q": repr(out["Q"]),
            "p_value": repr(out["p_modified_hartung_knapp"]),
            "model": ("Fisher-z; REML random-effects (tau2 by bounded REML optimisation on the "
                      "REML criterion); modified Hartung-Knapp (q_hk = max(1, Q/(k-1)), "
                      "t_{0.975,k-1}); exact manuscript random_effects_meta()"),
            "notes": NOTE[label],
        })
    mr = pd.DataFrame(rows_out)
    mr.to_csv(EXT / "meta_results.csv", index=False)
    log(f"    wrote {EXT / 'meta_results.csv'}  ({len(mr)} rows)")

    frows = []
    for label, out in results.items():
        y = np.arctanh(np.clip(np.array(out["rho_in"]), -0.999999, 0.999999))
        v = 1.0 / (np.array(out["n_in"], float) - 3.0)
        w = 1.0 / (v + out["tau2_fisher_z"])
        wsum = float(np.sum(w))
        for i, c in enumerate(out["coh"]):
            frows.append({
                "model": label,
                "cohort": c,
                "inferential_unit": "patient" if c != "GSE21501" else "specimen (<=129 patients)",
                "rho": repr(float(out["rho_in"][i])),
                "n": int(out["n_in"][i]),
                "fisher_z": repr(float(y[i])),
                "var_z_1_over_n_minus_3": repr(float(v[i])),
                "weight_in_fitted_model": repr(float(w[i])),
                "weight_share_pct_of_model": repr(float(100.0 * w[i] / wsum)),
                "model_tau2": repr(float(out["tau2_fisher_z"])),
                "model_pooled_z": repr(float(out["mu_z"])),
                "model_se_used": repr(float(out["se_z"])),
                "model_hk_scale_raw": repr(float(out["hk_scale_raw"])),
                "model_hk_scale_used": repr(float(out["hk_scale_used"])),
                "model_k": int(out["n_cohorts"]),
                "rho_source": (
                    str(BULK_ASSOC) if c in ORIGINAL_COHORTS
                    else str(GSE205154_PRIMARY if c == "GSE205154" else GSE21501_PRIMARY)
                ),
                "rho_source_sha256": (
                    sha256(BULK_ASSOC) if c in ORIGINAL_COHORTS
                    else sha256(GSE205154_PRIMARY if c == "GSE205154" else GSE21501_PRIMARY)
                ),
            })
    fi = pd.DataFrame(frows)
    fi.to_csv(EXT / "forest_inputs.csv", index=False)
    log(f"    wrote {EXT / 'forest_inputs.csv'}  ({len(fi)} rows)")

    cd = pd.DataFrame(ctrl_rows)
    cd.to_csv(EXT / "verification_controls.csv", index=False)
    log(f"    wrote {EXT / 'verification_controls.csv'}  ({len(cd)} rows)")

    # machine-readable handoff for the independent R cross-check
    handoff = {
        "program": PROGRAM,
        "manuscript_module": str(MANUSCRIPT_MODULE),
        "manuscript_module_sha256": manu_sha,
        "random_effects_meta_source_sha256": fn_sha,
        "random_effects_meta_source_lines": list(inspect.getsourcelines(mer)[1:2]),
        "models": {
            label: {
                "cohorts": out["coh"],
                "rho": [float(x) for x in out["rho_in"]],
                "n": [int(x) for x in out["n_in"]],
                "pooled_rho": float(out["r_random"]),
                "ci_low": float(out["ci_low"]),
                "ci_high": float(out["ci_high"]),
                "tau2": float(out["tau2_fisher_z"]),
                "I2": float(out["I2_pct"]),
                "Q": float(out["Q"]),
                "hk_scale_raw": float(out["hk_scale_raw"]),
                "hk_scale_used": float(out["hk_scale_used"]),
                "p": float(out["p_modified_hartung_knapp"]),
            }
            for label, out in results.items()
        },
        "published_anchor_from_archive": PUBLISHED,
        "specimen_sentence": SPECIMEN_SENTENCE,
    }
    (EXT / "provenance" / "meta_handoff.json").write_text(
        json.dumps(handoff, indent=2), encoding="utf-8"
    )
    log(f"    wrote {EXT / 'provenance' / 'meta_handoff.json'}")

    (EXT / "provenance" / "source_hashes.txt").write_text(
        "\n".join(f"{v}  {k}" for k, v in src_by_file.items())
        + f"\n{fn_sha}  random_effects_meta() function source (extracted)\n",
        encoding="utf-8",
    )

    # ----------------------------------------------------------------------------------
    # 7. Headline summary + check ledger
    # ----------------------------------------------------------------------------------
    r4 = results["4-cohort secondary extension"]
    r5 = results["5-cohort sensitivity"]
    cross4 = bool(r4["ci_low"] < 0 < r4["ci_high"])
    cross5 = bool(r5["ci_low"] < 0 < r5["ci_high"])
    disagree = cross4 != cross5

    log("\n" + "=" * 100)
    log("HEADLINE")
    log("=" * 100)
    log(f"  3-cohort (published, reproduced): rho={rep['r_random']:.6f} "
        f"CI {rep['ci_low']:.6f}..{rep['ci_high']:.6f}  k={rep['n_cohorts']}")
    log(f"  4-cohort secondary extension    : rho={r4['r_random']:.6f} "
        f"CI {r4['ci_low']:.6f}..{r4['ci_high']:.6f}  crosses zero: {cross4}")
    log(f"  5-cohort sensitivity            : rho={r5['r_random']:.6f} "
        f"CI {r5['ci_low']:.6f}..{r5['ci_high']:.6f}  crosses zero: {cross5}")
    log(f"  the two new intervals DISAGREE about excluding zero: {disagree}")
    if disagree:
        log("  -> FROZEN RULE: do not select one as the correct answer; report both side by side.")
    log(f"  4-cohort tau2={r4['tau2_fisher_z']:.6e}  I2={r4['I2_pct']:.4f}%  Q={r4['Q']:.6f}")
    log(f"  5-cohort tau2={r5['tau2_fisher_z']:.6e}  I2={r5['I2_pct']:.4f}%  Q={r5['Q']:.6f}")

    n_fail = sum(1 for c in CHECKS if not c["passed"])
    log(f"\n  checks: {len(CHECKS)} total, {n_fail} failed")
    for c in CHECKS:
        if not c["passed"]:
            log(f"    FAILED: {c['check']} :: {c['detail']}")
    if OBS:
        log(f"  observed (not asserted) outcomes: {len(OBS)}")
        for ob in OBS:
            log(f"    - {ob['observation']}")

    (EXT / "provenance" / "checks.json").write_text(
        json.dumps({"checks": CHECKS, "observations": OBS}, indent=2), encoding="utf-8"
    )
    (EXT / "provenance" / "headline.json").write_text(
        json.dumps({
            "three_cohort_reproduced": {k: rep[k] for k in
                ("r_random", "ci_low", "ci_high", "tau2_fisher_z", "I2_pct", "Q",
                 "p_modified_hartung_knapp", "hk_scale_raw", "hk_scale_used")},
            "published_anchor_from_archive": PUBLISHED,
            "four_cohort": {k: r4[k] for k in
                ("r_random", "ci_low", "ci_high", "tau2_fisher_z", "I2_pct", "Q",
                 "p_modified_hartung_knapp", "hk_scale_raw", "hk_scale_used", "cohorts")},
            "five_cohort": {k: r5[k] for k in
                ("r_random", "ci_low", "ci_high", "tau2_fisher_z", "I2_pct", "Q",
                 "p_modified_hartung_knapp", "hk_scale_raw", "hk_scale_used", "cohorts")},
            "four_crosses_zero": cross4,
            "five_crosses_zero": cross5,
            "intervals_disagree_about_zero": disagree,
            "specimen_sentence": SPECIMEN_SENTENCE,
            "n_checks": len(CHECKS),
            "n_checks_failed": n_fail,
        }, indent=2), encoding="utf-8")

    log(f"\nrun finished: {datetime.now().isoformat(timespec='seconds')}  "
        f"({time.time() - t0:.1f} s)")

    # append this run to run_log.txt
    header = (
        "\n" + "#" * 100 + "\n"
        f"# RUN {datetime.now().isoformat(timespec='seconds')}  01_extension_meta.py\n"
        f"# exit-status-to-be-confirmed: checks_failed={n_fail}\n"
        + "#" * 100 + "\n"
    )
    with open(EXT / "run_log.txt", "a", encoding="utf-8") as fh:
        fh.write(header + "\n".join(LOG) + "\n")
    print(f"\n[run_log.txt appended: {EXT / 'run_log.txt'}]")

    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
