"""
Cross-target reference panel analysis
=====================================
Implements PRE_SPECIFICATION.md (frozen 2026-07-14, IMMUTABLE).

Three-layer analysis:
  Layer 1: Target existence (expression stats)
  Layer 2: Expected program (one per target, leave-target-out)
  Layer 3: Competitive ecology (5 modules, leave-target-out per target)

Test families:
  Family A: 6 expected-program tests → BH-FDR within family
  Family B: 30 ecology tests (6×5) → BH-FDR within family
  Family C: Orthogonal anchoring (descriptive)

All 6 target-context pairs run in one pass. No result-based filtering.
"""
import warnings
warnings.filterwarnings('ignore')
import json, os, sys
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

np.random.seed(42)

# ── Paths ──
# Script: e:\Tex\HMlong\data\analysis\cross_target_panel\run_cross_target.py
# Need 4x dirname to reach project root: cross_target_panel -> analysis -> data -> HMlong
#
# In the original tree the inputs sat at <project root>/data/analysis and
# <project root>/data/TCGA, and walking four directories up from __file__ reached the
# project root. That layout does not exist in the public code deposit, so the root is
# now overridable: set SIGLEC15_PROJECT_ROOT to a tree containing
# data/analysis/*.csv and data/TCGA/*.csv to run this script elsewhere. When the
# variable is unset the original four-directory walk is used unchanged.
BASE = os.environ.get("SIGLEC15_PROJECT_ROOT") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
EXPR_DIR = os.path.join(BASE, 'data', 'analysis')
CLINICAL_DIR = os.path.join(BASE, 'data', 'TCGA')
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# ── 1. Gene sets ──

# Hallmark PI3K/AKT/MTOR (from analysis1_pi3k_multisource.py, MSigDB H: HALLMARK_PI3K_AKT_MTOR_SIGNALING)
HALLMARK_PI3K = [
    "AKT1","AKT2","AKT3","AKT1S1","BCL2","BCL2L1","BID","BIRC2","BIRC3",
    "CASP9","CDC37","CREB1","DEPTOR","EIF4B","EIF4E","EIF4EBP1","EIF4G1",
    "F2R","FKBP1A","FOXO1","FOXO3","GSK3A","GSK3B","IL6","ITGA1","ITGA2",
    "ITGA2B","ITGA3","ITGA4","ITGA5","ITGA6","ITGAV","ITGB1","ITGB3","ITGB5",
    "JUN","MAPK13","MAPK1","MAPK3","MAPKAP1","MLST8","MMP9","MTOR","MYC",
    "NFKB1","PDPK1","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2",
    "PIK3R3","PPARA","PPARD","PPARG","PRKAA1","PRKAA2","PRKCA","PRKCB",
    "PRKCG","PTEN","RHEB","RHOA","RPS6","RPS6KB1","RPTOR","SGK1","SGK3",
    "SREBF1","SYK","TSC1","TSC2","ULK1","VEGFA","VEGFB","VEGFC","YAP1",
    "YY1","ZFP36"
]

# Cytolytic immune (CD274 expected program)
CYTOLYTIC_IMMUNE = ['GZMA','GZMB','PRF1','IFNG','NKG7','GNLY','KLRD1','KLRB1']

# T-cell exhaustion / checkpoint (TIGIT expected program)
TCELL_CHECKPOINT = ['PDCD1','CTLA4','HAVCR2','LAG3','TIGIT','CD28','ICOS']

# HER2/ERBB signaling (ERBB2 expected program; ERBB2 itself excluded from definition)
HER2_ERBB = ['ERBB3','ERBB4','NRG1','NDRG1','GRB7','MED1']

# Epithelial differentiation (EPCAM expected program)
EPITHELIAL_DIFF = ['KRT18','KRT19','KRT8','EPCAM','CDH1','CLDN4','CLDN7']

# Macrophage / myeloid (CSF1R expected program)
MYELOID_MACRO = ['CD68','CSF1R','CD14','ITGAM','AIF1','LYZ','FCGR3A']

# Stromal / fibroblast
STROMAL = ['DCN','LUM','COL1A1','COL3A1','PDGFRB','ACTA2']

# MUC1/CD24 mucin-state module (from main manuscript modules_definition.json, 12 genes)
# NOTE: PRE_SPECIFICATION.md listed "MUC1, CD24, MUC3A, MUC4" as a clerical error.
# The spec states "from main manuscript"; the main manuscript module has 12 genes.
# Using the 12-gene module for cross-target comparability with main SIGLEC15 findings.
MUC_MODULE = ['CD24','CD44','EPCAM','KRT8','MUC1','MUC13','MUC16','MUC20',
              'MUC4','MUC5AC','MUC5B','MUC6']

# ── 2. Target-context panel (frozen) ──

TARGETS = [
    {
        'target': 'CD274',
        'cohorts': [('PAAD', 'TCGA_PAAD_178_log2cpm.csv')],
        'expected_program': 'cytolytic_immune',
        'expected_genes': CYTOLYTIC_IMMUNE,
        'reference_class': 'Checkpoint-centred (same-cancer positive calibration)',
    },
    {
        'target': 'TIGIT',
        'cohorts': [('LUAD', 'TCGA_LUAD_full_log2cpm.csv'),
                    ('LUSC', 'TCGA_LUSC_full_log2cpm.csv')],
        'expected_program': 'tcell_checkpoint',
        'expected_genes': TCELL_CHECKPOINT,
        'reference_class': 'Checkpoint-centred (cross-cancer calibration)',
    },
    {
        'target': 'ERBB2',
        'cohorts': [('BRCA', 'TCGA_BRCA_full_log2cpm.csv')],
        'expected_program': 'her2_erbb',
        'expected_genes': HER2_ERBB,
        'reference_class': 'Pathway-coupled',
    },
    {
        'target': 'EPCAM',
        'cohorts': [('PAAD', 'TCGA_PAAD_178_log2cpm.csv')],
        'expected_program': 'epithelial_diff',
        'expected_genes': EPITHELIAL_DIFF,
        'reference_class': 'Epithelial-state',
    },
    {
        'target': 'CSF1R',
        'cohorts': [('PAAD', 'TCGA_PAAD_178_log2cpm.csv')],
        'expected_program': 'myeloid_macro',
        'expected_genes': MYELOID_MACRO,
        'reference_class': 'Myeloid-state',
    },
    {
        'target': 'SIGLEC15',
        'cohorts': [('PAAD', 'TCGA_PAAD_178_log2cpm.csv')],
        'expected_program': 'pi3k_hallmark',
        'expected_genes': HALLMARK_PI3K,
        'reference_class': 'Identity-mismatch (core case)',
    },
]

ECOLOGY_MODULES = {
    'T_cell_cytolytic': CYTOLYTIC_IMMUNE,
    'Myeloid': MYELOID_MACRO,
    'Epithelial': EPITHELIAL_DIFF,
    'Stromal': STROMAL,
    'MUC_state': MUC_MODULE,
}

B_BOOTSTRAP = 5000


# ── 3. Scoring ──

def score_mean_z(expr_sub):
    """Mean of z-scored expression across genes (genes × samples).
    Input: pandas DataFrame (genes × samples). Returns pandas Series indexed by sample."""
    # Convert to numpy array to support keepdims (pandas mean/std don't accept keepdims)
    arr = expr_sub.values if hasattr(expr_sub, 'values') else np.asarray(expr_sub, dtype=float)
    mu = arr.mean(axis=1, keepdims=True)
    sd = arr.std(axis=1, ddof=1, keepdims=True)
    sd[sd < 1e-10] = 1e-10
    z = (arr - mu) / sd
    scores = z.mean(axis=0)
    # Return as pandas Series to preserve sample labels
    if hasattr(expr_sub, 'columns'):
        return pd.Series(scores, index=expr_sub.columns)
    return pd.Series(scores)


def leave_target_out(genes, target):
    """Remove target gene from gene set if present."""
    return [g for g in genes if g != target]


def filter_to_matrix(genes, expr_index):
    """Keep only genes present in expression matrix."""
    return [g for g in genes if g in expr_index]


# ── 4. Data loading ──

_expr_cache = {}

def load_expr(filename):
    if filename not in _expr_cache:
        path = os.path.join(EXPR_DIR, filename)
        df = pd.read_csv(path, index_col=0)
        # Aggregate duplicate gene names (take mean) — required for .loc indexing
        # and for cross-cohort pooling via pd.concat
        if df.index.duplicated().any():
            n_dup = int(df.index.duplicated().sum())
            print(f"  Aggregating {n_dup} duplicate gene names in {filename}")
            df = df.groupby(level=0, sort=False).mean()
        _expr_cache[filename] = df
        print(f"  Loaded {filename}: {df.shape[0]} genes × {df.shape[1]} samples")
    return _expr_cache[filename]


def get_pooled_expr(cohorts):
    """Load and pool expression from multiple cohorts (column-wise concat)."""
    frames = []
    for cohort_name, filename in cohorts:
        df = load_expr(filename)
        frames.append(df)
    if len(frames) == 1:
        return frames[0]
    # Pool: align on gene rows (intersection)
    common_genes = frames[0].index
    for f in frames[1:]:
        common_genes = common_genes.intersection(f.index)
    pooled = pd.concat([f.loc[common_genes] for f in frames], axis=1)
    return pooled


# ── 5. Correlation + bootstrap ──

def spearman_with_ci(x, y, n_boot=B_BOOTSTRAP):
    """Spearman ρ with patient-level bootstrap 95% CI."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    rho, p = stats.spearmanr(x, y)
    # Patient-level bootstrap
    rng = np.random.default_rng(42)
    boot_rhos = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        xb, yb = x[idx], y[idx]
        if np.std(xb) > 0 and np.std(yb) > 0:
            boot_rhos[i] = stats.spearmanr(xb, yb)[0]
        else:
            boot_rhos[i] = np.nan
    boot_rhos = boot_rhos[~np.isnan(boot_rhos)]
    ci_lo = np.percentile(boot_rhos, 2.5)
    ci_hi = np.percentile(boot_rhos, 97.5)
    return {
        'rho': float(rho),
        'p': float(p),
        'ci_lo': float(ci_lo),
        'ci_hi': float(ci_hi),
        'n': int(n),
        'n_boot_valid': int(len(boot_rhos)),
    }


# ── 6. Main analysis ──

def run_analysis():
    results = {
        'spec_version': 'PRE_SPECIFICATION.md frozen 2026-07-14',
        'method_notes': {
            'scoring': 'mean-z (consistent with main manuscript primary method)',
            'leave_target_out': 'Applied to BOTH expected program AND ecology modules '
                                '(target gene removed from any gene set containing it)',
            'muc_module': '12-gene MUC1_CD24_Module from main manuscript '
                          '(PRE_SPECIFICATION listed 4-gene subset as clerical error; '
                          'spec states "from main manuscript")',
            'tigit_pooling': 'LUAD + LUSC pooled for single NSCLC test (per spec: 1 test)',
            'bootstrap': f'Patient-level non-parametric B={B_BOOTSTRAP}',
            'fdr': 'BH within each pre-specified family',
        },
        'layer1_existence': {},
        'layer2_expected_program': {},
        'layer3_ecology': {},
        'family_a_fdr': {},
        'family_b_fdr': {},
        'erbb2_sensitivity': {},
    }

    family_a_records = []  # (target, rho, p)
    family_b_records = []  # (target, ecology, rho, p)

    for tinfo in TARGETS:
        target = tinfo['target']
        print(f"\n{'='*60}")
        print(f"Processing: {target} ({tinfo['reference_class']})")

        expr = get_pooled_expr(tinfo['cohorts'])
        cohort_label = '+'.join([c[0] for c in tinfo['cohorts']])

        # --- Layer 1: Existence ---
        if target not in expr.index:
            print(f"  WARNING: {target} not found in expression matrix!")
            results['layer1_existence'][target] = {
                'error': f'{target} not in expression matrix',
                'cohort': cohort_label,
            }
            continue

        target_expr = expr.loc[target]
        det_rate = float((target_expr > 0).mean())
        results['layer1_existence'][target] = {
            'cohort': cohort_label,
            'n_samples': int(len(target_expr)),
            'detection_rate': det_rate,
            'median_log2cpm': float(target_expr.median()),
            'mean_log2cpm': float(target_expr.mean()),
            'q25': float(target_expr.quantile(0.25)),
            'q75': float(target_expr.quantile(0.75)),
        }
        print(f"  n={len(target_expr)}, det_rate={det_rate:.3f}, "
              f"median={target_expr.median():.2f}")

        # --- Layer 2: Expected program (leave-target-out) ---
        prog_genes = leave_target_out(tinfo['expected_genes'], target)
        prog_genes = filter_to_matrix(prog_genes, expr.index)
        if len(prog_genes) < 3:
            print(f"  WARNING: Only {len(prog_genes)} genes for expected program!")
        prog_score = score_mean_z(expr.loc[prog_genes])
        prog_result = spearman_with_ci(target_expr.values, prog_score.values)
        prog_result['program'] = tinfo['expected_program']
        prog_result['n_genes_used'] = len(prog_genes)
        prog_result['target_removed'] = target in tinfo['expected_genes']
        results['layer2_expected_program'][target] = prog_result
        family_a_records.append((target, prog_result['rho'], prog_result['p']))
        print(f"  Expected program ({tinfo['expected_program']}): "
              f"ρ={prog_result['rho']:.3f}, p={prog_result['p']:.2e}, "
              f"CI=[{prog_result['ci_lo']:.3f}, {prog_result['ci_hi']:.3f}]")

        # --- Layer 3: Competitive ecology (leave-target-out per target) ---
        results['layer3_ecology'][target] = {}
        for eco_name, eco_genes_full in ECOLOGY_MODULES.items():
            eco_genes = leave_target_out(eco_genes_full, target)
            eco_genes = filter_to_matrix(eco_genes, expr.index)
            if len(eco_genes) < 3:
                print(f"    WARNING: {eco_name} has only {len(eco_genes)} genes!")
            eco_score = score_mean_z(expr.loc[eco_genes])
            eco_result = spearman_with_ci(target_expr.values, eco_score.values)
            eco_result['n_genes_used'] = len(eco_genes)
            eco_result['target_removed'] = target in eco_genes_full
            results['layer3_ecology'][target][eco_name] = eco_result
            family_b_records.append((target, eco_name, eco_result['rho'],
                                     eco_result['p']))
            print(f"    {eco_name}: ρ={eco_result['rho']:.3f}, "
                  f"p={eco_result['p']:.2e}")

    # ── 7. BH-FDR within families ──

    # Family A: 6 expected program tests
    if family_a_records:
        p_vals_a = np.array([r[2] for r in family_a_records])
        _, fdr_a, _, _ = multipletests(p_vals_a, method='fdr_bh')
        for i, (target, rho, p) in enumerate(family_a_records):
            results['family_a_fdr'][target] = float(fdr_a[i])
            results['layer2_expected_program'][target]['fdr_bh'] = float(fdr_a[i])
        print(f"\n{'='*60}")
        print("Family A (expected program) BH-FDR:")
        for i, (target, rho, p) in enumerate(family_a_records):
            print(f"  {target}: ρ={rho:.3f}, p={p:.2e}, FDR={fdr_a[i]:.3f}")

    # Family B: 30 ecology tests
    if family_b_records:
        p_vals_b = np.array([r[3] for r in family_b_records])
        _, fdr_b, _, _ = multipletests(p_vals_b, method='fdr_bh')
        for i, (target, eco, rho, p) in enumerate(family_b_records):
            key = f"{target}__{eco}"
            results['family_b_fdr'][key] = float(fdr_b[i])
            results['layer3_ecology'][target][eco]['fdr_bh'] = float(fdr_b[i])
        print(f"\nFamily B (ecology) BH-FDR:")
        for i, (target, eco, rho, p) in enumerate(family_b_records):
            print(f"  {target}×{eco}: ρ={rho:.3f}, p={p:.2e}, FDR={fdr_b[i]:.3f}")

    # ── 8. ERBB2-BRCA sensitivity analysis ──
    print(f"\n{'='*60}")
    print("ERBB2-BRCA sensitivity analysis...")
    erbb2_sensitivity = run_erbb2_sensitivity()
    results['erbb2_sensitivity'] = erbb2_sensitivity

    # ── 9. Save ──
    out_json = os.path.join(OUT_DIR, 'cross_target_results.json')
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {out_json}")

    # Also save flat CSV for Table S22
    save_flat_csv(results, family_b_records)
    return results


def run_erbb2_sensitivity():
    """ERBB2-BRCA: full cohort + HER2-enriched (if PAM50) + ERBB2-high (expression proxy)."""
    sens = {}
    expr_path = os.path.join(EXPR_DIR, 'TCGA_BRCA_full_log2cpm.csv')
    expr = pd.read_csv(expr_path, index_col=0)

    if 'ERBB2' not in expr.index:
        sens['error'] = 'ERBB2 not in BRCA expression matrix'
        return sens

    erbb2_expr = expr.loc['ERBB2']
    n_total = len(erbb2_expr)
    sens['full_cohort'] = {'n': n_total}

    # ERBB2-high subgroup (top quartile as expression proxy for amplification)
    q75 = erbb2_expr.quantile(0.75)
    high_mask = erbb2_expr >= q75
    sens['erbb2_high_quartile'] = {
        'n': int(high_mask.sum()),
        'threshold_log2cpm': float(q75),
        'note': 'Top quartile of ERBB2 expression as proxy for amplification '
                 '(CNV data not loaded)',
    }

    # Try PAM50 subgroup
    pam50_path = os.path.join(BASE, 'data', 'analysis',
                              'Task3_BRCA_ER_subgroup', 'brca_clinical_wide.csv')
    if os.path.exists(pam50_path):
        try:
            clin = pd.read_csv(pam50_path, index_col=0)
            # Look for PAM50 column
            pam50_col = None
            for col in clin.columns:
                if 'pam' in col.lower() or 'subtype' in col.lower():
                    pam50_col = col
                    break
            if pam50_col and 'HER2' in str(clin[pam50_col].unique()):
                her2_samples = clin[clin[pam50_col].astype(str).str.contains(
                    'HER2', case=False, na=False)].index
                # Match to expression columns
                her2_in_expr = [s for s in her2_samples if s in expr.columns]
                sens['her2_enriched_pam50'] = {
                    'n': len(her2_in_expr),
                    'pam50_column': pam50_col,
                    'note': 'PAM50 HER2-enriched subgroup',
                }
                # Compute correlations in HER2-enriched subgroup
                if len(her2_in_expr) >= 10:
                    expr_her2 = expr[her2_in_expr]
                    erbb2_her2 = expr_her2.loc['ERBB2']
                    # Expected program (HER2/ERBB, leave-ERBB2-out)
                    prog_genes = [g for g in HER2_ERBB if g != 'ERBB2'
                                  and g in expr_her2.index]
                    if len(prog_genes) >= 3:
                        prog_score = score_mean_z(expr_her2.loc[prog_genes])
                        corr = spearman_with_ci(erbb2_her2.values,
                                                prog_score.values,
                                                n_boot=min(B_BOOTSTRAP, 2000))
                        corr['n_genes_used'] = len(prog_genes)
                        sens['her2_enriched_pam50']['expected_program_corr'] = corr
            else:
                sens['her2_enriched_pam50'] = {
                    'note': 'PAM50 HER2-enriched subgroup not identifiable '
                            'in clinical data',
                }
        except Exception as e:
            sens['her2_enriched_pam50'] = {'error': str(e)}
    else:
        sens['her2_enriched_pam50'] = {
            'note': 'PAM50 clinical data file not found',
        }

    # Compute expected program correlation in ERBB2-high subgroup
    if high_mask.sum() >= 10:
        expr_high = expr.loc[:, high_mask]
        erbb2_high = expr_high.loc['ERBB2']
        prog_genes = [g for g in HER2_ERBB if g != 'ERBB2'
                      and g in expr_high.index]
        if len(prog_genes) >= 3:
            prog_score = score_mean_z(expr_high.loc[prog_genes])
            corr = spearman_with_ci(erbb2_high.values, prog_score.values,
                                    n_boot=min(B_BOOTSTRAP, 2000))
            corr['n_genes_used'] = len(prog_genes)
            sens['erbb2_high_quartile']['expected_program_corr'] = corr

    return sens


def save_flat_csv(results, family_b_records):
    """Save flat CSV for Supplementary Table S22."""
    rows = []
    # Layer 2
    for target, vals in results['layer2_expected_program'].items():
        if 'error' in vals:
            continue
        rows.append({
            'Target': target,
            'Layer': 'Expected program',
            'Program/Ecology': vals['program'],
            'Spearman_rho': vals['rho'],
            'p_value': vals['p'],
            'FDR_BH': vals.get('fdr_bh', np.nan),
            'CI_lo': vals['ci_lo'],
            'CI_hi': vals['ci_hi'],
            'n': vals['n'],
            'n_genes': vals.get('n_genes_used', np.nan),
            'target_removed': vals.get('target_removed', False),
        })
    # Layer 3
    for target, ecos in results['layer3_ecology'].items():
        for eco_name, vals in ecos.items():
            if 'error' in vals:
                continue
            rows.append({
                'Target': target,
                'Layer': 'Competitive ecology',
                'Program/Ecology': eco_name,
                'Spearman_rho': vals['rho'],
                'p_value': vals['p'],
                'FDR_BH': vals.get('fdr_bh', np.nan),
                'CI_lo': vals['ci_lo'],
                'CI_hi': vals['ci_hi'],
                'n': vals['n'],
                'n_genes': vals.get('n_genes_used', np.nan),
                'target_removed': vals.get('target_removed', False),
            })
    df = pd.DataFrame(rows)
    out_csv = os.path.join(OUT_DIR, 'cross_target_results_flat.csv')
    df.to_csv(out_csv, index=False)
    print(f"Flat CSV saved to {out_csv}")


if __name__ == '__main__':
    print("=" * 60)
    print("CROSS-TARGET REFERENCE PANEL ANALYSIS")
    print("Per PRE_SPECIFICATION.md (frozen 2026-07-14, IMMUTABLE)")
    print("=" * 60)
    run_analysis()
