# =============================================================================
# 08_final_analysis.R -- Phase B step 2, GSE278689 (FROZEN RUN)
#
# MODALITY: snRNA (single NUCLEI). This is a CROSS-MODALITY SENSITIVITY ANALYSIS.
# It is NOT pooled with the manuscript's 19-sample scRNA estimate and NOT pooled
# with GSE278688. It is a single stratum and is reported as one.
#
# LANGUAGE CONVENTION (author's pre-registered interpretation criteria): the four
# words named in those criteria are never applied to a sample, estimate, interval
# or contrast. This run answers only whether the target-state relation remains
# estimable, and in what direction, under a different transcript-capture process.
#
# Every convention below is either a FROZEN constant or read from the manuscript's
# own pipeline -- none was chosen by this script:
#
#   eligible tumour label : tissue == "Tumor"  (mapping fixed before any
#                           expression-based eligibility was evaluated; 02_)
#   cell type             : deposited all_celltype == "Myeloid" (no re-annotation)
#   >= 50 included tumour-myeloid cells per sample
#   both groups present   : NO minimum group size. The `>=10 cells with >=2
#                           detectable and >=5 non-detectable` rule belongs to the
#                           SEPARATE restricted TREM2_TAM sensitivity analysis and
#                           is NOT applied here.
#   programme score       : PROGRAM = ACP5 + MMP9. Each gene z-standardised WITHIN
#                           the biological sample, then averaged.
#                           zscore_columns() uses ddof=0 (population SD).
#                           An invariant gene would be OMITTED for that sample;
#                           the primary score REQUIRES BOTH genes, so such a
#                           sample would leave the primary contrast. (None did.)
#                           SIGLEC15 NEVER enters the programme score.
#   sample contrasts      : equally weighted
#   detection             : non-zero RAW deposited count (invariant to normalisation)
#   normalisation         : deposit is RAW INTEGER COUNTS -> CP10K + log1p,
#                           matching the manuscript pipeline (see 05_)
#   aggregation of RNA    : inverse-log = log1p(mean(expm1(x))), the manuscript's
#                           own rule (analysis_F_cross_scale.py line 324)
#   bootstrap             : 10000 resamples, seed 15072026, resampling unit =
#                           biological sample WITHIN the single dataset stratum
#
# SEED 15072026   N_RESAMPLES 10000
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

SEED <- 15072026L
N_RESAMPLES <- 10000L
set.seed(SEED)
ROOT <- "E:/Tex/HMlong/BJC_EXT_v1_20260927/phaseB_GSE278689"
OUT  <- file.path(ROOT, "output")
P <- function(...) cat(..., "\n", sep = "")
P("SEED=", SEED, "  N_RESAMPLES=", N_RESAMPLES)
P("time=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))

d    <- fread(file.path(OUT, "_myeloid_gene_cache.csv"))
allc <- fread(file.path(OUT, "_allcells_siglec15_cache.csv"))
man  <- fread(file.path(OUT, "sample_manifest.csv"))
P("myeloid cache : ", nrow(d), " x ", ncol(d), "   all-cell cache: ", nrow(allc))
elig <- man[final_eligible == TRUE, sample_id]
P("final estimable samples: ", length(elig), " -> ", paste(elig, collapse = ", "))

## ---- CP10K + log1p reconstruction (deposit is raw integer counts) ----------
d[, c("SIGLEC15_n", "ACP5_n", "MMP9_n") :=
    lapply(.SD, function(v) log1p(v * 1e4 / libsize)),
  .SDcols = c("SIGLEC15", "ACP5", "MMP9")]
allc[, SIGLEC15_n := log1p(SIGLEC15 * 1e4 / libsize)]
allc[, SIGLEC15_cpm := SIGLEC15 * 1e4 / libsize]

## zscore_columns() from analysis_F_cross_scale.py: ddof = 0 (population SD).
## stats::sd uses ddof = 1, hence the correction. Zero/undefined SD -> NA, which
## np.nanmean then skips -- reproduced by returning NA here.
zpop <- function(x) {
  s <- stats::sd(x) * sqrt((length(x) - 1) / length(x))
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x)) / s
}

## ===================== 1. PER-SAMPLE PRIMARY CONTRAST ======================
P("\n===== 1. PER-SAMPLE PRIMARY CONTRAST =====")
P("score   : mean of z(ACP5) and z(MMP9), each z-standardised within the sample (ddof=0)")
P("contrast: mean score in SIGLEC15 RNA-detectable MINUS non-detectable myeloid cells")
recs <- list()
for (s in elig) {
  dd <- d[sample_id == s]
  za <- zpop(dd$ACP5_n); zm <- zpop(dd$MMP9_n)
  acp5_inv <- all(is.na(za)); mmp9_inv <- all(is.na(zm))
  if (acp5_inv || mmp9_inv) {
    recs[[s]] <- data.table(sample_id = s, patient = dd$patient[1],
      detectable_cells = sum(dd$SIGLEC15 > 0), nondetectable_cells = sum(dd$SIGLEC15 == 0),
      mean_score_detectable = NA_real_, mean_score_nondetectable = NA_real_,
      contrast = NA_real_, acp5_invariant = acp5_inv, mmp9_invariant = mmp9_inv,
      genes_in_score = NA_character_,
      status = "EXCLUDED: programme gene invariant within sample")
    next
  }
  score <- (za + zm) / 2                      # both genes present
  det <- dd$SIGLEC15 > 0
  recs[[s]] <- data.table(sample_id = s, patient = dd$patient[1],
    detectable_cells = sum(det), nondetectable_cells = sum(!det),
    mean_score_detectable = mean(score[det]), mean_score_nondetectable = mean(score[!det]),
    contrast = mean(score[det]) - mean(score[!det]),
    acp5_invariant = FALSE, mmp9_invariant = FALSE, genes_in_score = "ACP5;MMP9",
    status = "OK")
}
res <- rbindlist(recs)
res[, direction_positive := contrast > 0]
setorder(res, -contrast, na.last = TRUE)
fwrite(res, file.path(OUT, "primary_contrast_by_sample.csv"))
print(res, nrows = 20)

x <- res$contrast[is.finite(res$contrast)]
n <- length(x)
P("\nmodality                       : snRNA (single nuclei)")
P("samples with an estimable contrast: ", n, " of ", length(elig), " eligible")
P("score genes used in every sample  : ACP5;MMP9 (both required)")
P("SIGLEC15 in the programme score   : NO (asserted below)")
stopifnot(!any(c("SIGLEC15_n") %in% c("ACP5_n", "MMP9_n")))     # target never a program gene
P("target-in-program assertion       : PASS (SIGLEC15 is not an argument of the score)")

## ===================== 2. SUMMARY ==========================================
P("\n===== 2. SUMMARY OF THE PRIMARY CONTRAST =====")
P("STRATUM: this is a SINGLE stratum (GSE278689 snRNA tumour samples). It is NOT")
P("         pooled with the manuscript's 19-sample scRNA estimate and NOT pooled")
P("         with GSE278688. No cross-dataset combination is computed anywhere.")
P("n                  : ", n)
P("mean contrast      : ", sprintf("%+.4f", mean(x)))
P("sd                 : ", sprintf("%.4f", stats::sd(x)))
P("median             : ", sprintf("%+.4f", stats::median(x)))
P("range              : ", sprintf("%+.4f to %+.4f", min(x), max(x)))
P("positive direction : ", sum(x > 0), " / ", n)
P("negative direction : ", sum(x < 0), " / ", n)
P("zero               : ", sum(x == 0), " / ", n)

set.seed(SEED)
boot <- replicate(N_RESAMPLES, mean(sample(x, n, replace = TRUE)))
ci   <- unname(stats::quantile(boot, c(0.025, 0.975), type = 7))
P("\nbootstrap          : ", N_RESAMPLES, " resamples, seed ", SEED,
  ", type = 7 quantiles")
P("resampling unit    : biological sample, within the single dataset stratum")
P("95% CI             : ", sprintf("%+.4f to %+.4f", ci[1], ci[2]))
P("CI excludes zero   : ", ci[1] > 0 || ci[2] < 0)

wr  <- suppressWarnings(stats::wilcox.test(x, mu = 0, alternative = "two.sided",
                                           correct = FALSE, exact = FALSE))
wr0 <- suppressWarnings(stats::wilcox.test(x, mu = 0, alternative = "two.sided",
                                           correct = TRUE,  exact = FALSE))
st  <- stats::binom.test(sum(x > 0), n, 0.5)
P("signed-rank p (no continuity correction, = scipy default) : ",
  format.pval(wr$p.value, digits = 4), "   V = ", wr$statistic)
P("signed-rank p (continuity-corrected)                      : ",
  format.pval(wr0$p.value, digits = 4))
P("sign test p                                               : ",
  format.pval(st$p.value, digits = 4))

sm <- data.table(
  dataset = "GSE278689", modality = "snRNA (single nuclei)",
  analysis_role = "cross-modality sensitivity analysis (single stratum; not pooled)",
  stratum = "GSE278689 snRNA tumour samples (single stratum)",
  n_samples = n, n_eligible = length(elig),
  mean_contrast = mean(x), sd_contrast = stats::sd(x),
  median_contrast = stats::median(x), min_contrast = min(x), max_contrast = max(x),
  n_positive = sum(x > 0), n_negative = sum(x < 0),
  score = "ACP5;MMP9 (both required; gene z-standardised within sample, ddof=0)",
  normalisation = "CP10K then log1p (deposit verified as raw integer counts)",
  detection = "non-zero raw deposited SIGLEC15 count",
  boot_resamples = N_RESAMPLES, boot_seed = SEED,
  ci_lower = ci[1], ci_upper = ci[2],
  wilcoxon_V = as.numeric(wr$statistic), wilcoxon_p = wr$p.value,
  sign_test_p = st$p.value)
fwrite(sm, file.path(OUT, "primary_summary.csv"))
P(""); print(sm)

## ===================== 3. FOUR SPECIMEN-LEVEL SUMMARIES ====================
P("\n===== 3. FOUR SAME-SAMPLE SPECIMEN-LEVEL SIGLEC15 SUMMARIES =====")
invlog <- function(v) log1p(mean(expm1(v)))     # manuscript aggregation rule
spec <- rbindlist(lapply(elig, function(s) {
  dd  <- d[sample_id == s]
  ac  <- allc[sample_id == s]                   # EVERY cell in the retained sample
  det <- dd$SIGLEC15 > 0
  data.table(sample_id = s, patient = dd$patient[1],
    agg_myeloid_siglec15        = invlog(dd$SIGLEC15_n),        # manuscript rule
    agg_myeloid_siglec15_cp10k  = mean(dd$SIGLEC15_n),          # convention control
    detectable_fraction         = mean(det),
    mean_expr_detectable        = if (any(det)) mean(dd$SIGLEC15_n[det]) else NA_real_,
    allcell_tumour_agg          = invlog(ac$SIGLEC15_n),        # manuscript rule
    allcell_tumour_agg_cp10k    = mean(ac$SIGLEC15_cpm),        # convention control
    n_all_cells                 = nrow(ac),
    n_detectable_allcell        = sum(ac$SIGLEC15 > 0))
}))
spec <- merge(spec, res[, .(sample_id, contrast)], by = "sample_id")
setorder(spec, -contrast)
fwrite(spec, file.path(OUT, "specimen_level_summaries.csv"))
print(spec, nrows = 20)

spearman_boot <- function(xv, yv) {
  ok <- is.finite(xv) & is.finite(yv); xv <- xv[ok]; yv <- yv[ok]; m <- length(xv)
  if (m < 4) return(list(rho = NA_real_, p = NA_real_, lo = NA_real_, hi = NA_real_, n = m))
  rho <- suppressWarnings(stats::cor(xv, yv, method = "spearman"))
  pt  <- suppressWarnings(stats::cor.test(xv, yv, method = "spearman", exact = FALSE))
  set.seed(SEED)
  bs <- vapply(seq_len(N_RESAMPLES), function(b) {
    i <- sample.int(m, m, replace = TRUE)
    if (length(unique(xv[i])) < 2 || length(unique(yv[i])) < 2) return(NA_real_)
    suppressWarnings(stats::cor(xv[i], yv[i], method = "spearman"))
  }, numeric(1))
  q <- unname(stats::quantile(bs, c(0.025, 0.975), na.rm = TRUE, type = 7))
  list(rho = unname(rho), p = pt$p.value, lo = q[1], hi = q[2], n = m)
}
summ_vars <- c("agg_myeloid_siglec15", "detectable_fraction",
               "mean_expr_detectable", "allcell_tumour_agg")
cors <- rbindlist(lapply(summ_vars, function(v) {
  r <- spearman_boot(spec[[v]], spec$contrast)
  data.table(summary_statistic = v, n = r$n, spearman_rho = r$rho,
             ci_lower = r$lo, ci_upper = r$hi, p_value = r$p,
             ci_method = sprintf("sample bootstrap, %d resamples, seed %d, type=7",
                                 N_RESAMPLES, SEED))
}))
## convention controls: the other aggregation rule for the two aggregate measures
conv <- rbindlist(lapply(c("agg_myeloid_siglec15_cp10k", "allcell_tumour_agg_cp10k"),
  function(v) {
    r <- spearman_boot(spec[[v]], spec$contrast)
    data.table(summary_statistic = v, n = r$n, spearman_rho = r$rho,
               ci_lower = r$lo, ci_upper = r$hi, p_value = r$p,
               ci_method = "convention control: plain CP10K mean instead of inverse-log")
  }))
cors_all <- rbind(cors, conv, fill = TRUE)
fwrite(cors_all, file.path(OUT, "specimen_level_correlations.csv"))
P("\n-- Spearman correlations with the sample-specific contrast --")
print(cors_all, nrows = 20)
P("\nattainable Spearman grid at n=", n, " is k/", n * (n - 1) / 2,
  " -- these estimates carry very little information at this n.")
boot_neg1 <- sum(vapply(seq_len(N_RESAMPLES), function(b) {
  i <- sample.int(n, n, replace = TRUE)
  suppressWarnings(stats::cor(spec$agg_myeloid_siglec15[i], spec$contrast[i],
                              method = "spearman")) == -1
}, logical(1)), na.rm = TRUE)
P("of ", N_RESAMPLES, " bootstrap resamples, ", boot_neg1,
  " land on exactly rho = -1 (a degeneracy of small n, not an estimate).")

## ===================== 4. DETECTABILITY / COMPARABILITY ====================
P("\n===== 4. DETECTABILITY AND MODALITY COMPARABILITY =====")
P("MODALITY: snRNA (nuclei). The published primary series and the GSE278688")
P("extension are scRNA (whole cells). These are DIFFERENT MEASUREMENT PROCESSES.")
P("The comparison below is a modality/dataset difference statement. It is NOT a")
P("quality ranking. Only the phrasings permitted by the author's pre-registered")
P("interpretation criteria are used; see results_summary.md.")
det_my <- d[, .(myeloid_cells = .N, detectable = sum(SIGLEC15 > 0)), by = sample_id]
det_my <- merge(man[, .(sample_id, patient, final_eligible)], det_my, by = "sample_id")
setorder(det_my, -final_eligible, patient)
print(det_my)
elig_det <- det_my[final_eligible == TRUE]
P("\n-- eligible tumour-myeloid cells --")
P("myeloid cells                       : ", sum(elig_det$myeloid_cells))
P("SIGLEC15 RNA-detectable             : ", sum(elig_det$detectable))
P("detectable fraction (eligible)      : ",
  sprintf("%.4f%%", 100 * sum(elig_det$detectable) / sum(elig_det$myeloid_cells)))
P("detectable fraction (all 8 samples) : ",
  sprintf("%.4f%%", 100 * d[, sum(SIGLEC15 > 0)] / nrow(d)))
P("\n-- SIGLEC15 across ALL cells of the series (any cell type) --")
P("all cells                           : ", nrow(allc))
P("SIGLEC15-detectable cells           : ", allc[, sum(SIGLEC15 > 0)], " (",
  sprintf("%.4f%%", 100 * allc[, sum(SIGLEC15 > 0)] / nrow(allc)), ")")
P("myeloid share of all cells          : ",
  sprintf("%.2f%%", 100 * nrow(d) / nrow(allc)))
P("myeloid share of SIGLEC15+ cells    : ",
  sprintf("%.2f%%", 100 * d[, sum(SIGLEC15 > 0)] / allc[, sum(SIGLEC15 > 0)]))

cmp <- data.table(
  dataset = c("Published scRNA series (manuscript primary)",
              "GSE278688 (scRNA, whole cells)",
              "GSE278689 (snRNA, nuclei) -- THIS RUN"),
  modality = c("scRNA", "scRNA", "snRNA"),
  role = c("primary 19-sample estimate", "extension", "cross-modality sensitivity analysis"),
  detectable_myeloid_fraction_pct = c(2.0, 0.57,
    100 * sum(elig_det$detectable) / sum(elig_det$myeloid_cells)),
  n_myeloid = c(9918L, 9102L, as.integer(sum(elig_det$myeloid_cells))),
  n_detectable = c(197L, 52L, as.integer(sum(elig_det$detectable))))
cmp[, ratio_vs_published_scrna := detectable_myeloid_fraction_pct /
      detectable_myeloid_fraction_pct[1]]
fwrite(cmp, file.path(OUT, "detectability_comparability.csv"))
P("\n-- detectability comparison (modality/dataset difference, not a ranking) --")
print(cmp)
P("\nThis run is NOT pooled with either scRNA figure and no combined estimate is formed.")
P("The denominators are not equated: the published fraction is over that series'")
P("myeloid cells in its own eligible samples; this one is over this series' 6")
P("estimable tumour samples' myeloid cells.")

P("\nDONE ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
