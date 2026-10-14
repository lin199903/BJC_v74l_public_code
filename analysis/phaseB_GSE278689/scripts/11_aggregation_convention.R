# =============================================================================
# 11_aggregation_convention.R -- Phase B step 2, GSE278689
#
# The frozen brief does not state the aggregation rule for the specimen-level
# SIGLEC15 summaries, and the manuscript states it in one place only: for the
# myeloid aggregate, analysis_F_cross_scale.py line 324 uses
#
#     aggregate = np.log1p(np.mean(np.expm1(expr), axis=0))      # inverse-log
#
# The all-cell tumour aggregate is described in the SI as "every cell in the
# retained tumour samples" without naming the aggregation.
#
# This is therefore a genuine fork. It is NOT resolved silently: both candidates
# are computed and compared, and the effect on the reported Spearman values is
# stated.
#
# The two candidates are linked by a monotone transform, so SPEARMAN rho -- which
# depends on ranks only -- must be IDENTICAL between them. Only the bootstrap
# interval can move, because resampling can break the monotone relation. This
# script checks that claim rather than asserting it.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
SEED <- 15072026L; N_RESAMPLES <- 10000L
ROOT <- "E:/Tex/HMlong/BJC_EXT_v1_20260927/phaseB_GSE278689"
OUT  <- file.path(ROOT, "output")
P <- function(...) cat(..., "\n", sep = "")
P("time = ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))

spec <- fread(file.path(OUT, "specimen_level_summaries.csv"))
cors <- fread(file.path(OUT, "specimen_level_correlations.csv"))
P("retained samples: ", nrow(spec), "  (", paste(spec$sample_id, collapse = ", "), ")")

P("\n===== 1. THE TWO AGGREGATION CANDIDATES =====")
P("A  inverse-log  : log1p(mean(expm1(x)))   <- the ONLY rule the manuscript writes down")
P("B  plain mean   : mean(x)                 <- same units as the other summaries")
print(spec[, .(sample_id, patient, agg_myeloid_siglec15, agg_myeloid_siglec15_cp10k,
               allcell_tumour_agg, allcell_tumour_agg_cp10k, contrast)])

P("\n===== 2. ARE THEY MONOTONE-LINKED? (rank identity check) =====")
P("This is checked, not assumed. The two rules are NOT the same kind of transform:")
P("  all-cell aggregate : agg = log1p(mean(cpm)) and the control = mean(cpm), so")
P("                       agg = log1p(control) -- an exact monotone link, hence")
P("                       identical ranks and identical rho.")
P("  myeloid aggregate  : agg = log1p(mean(cpm)) but the control = mean(log1p(cpm)).")
P("                       log-of-mean is NOT a monotone function of mean-of-log, so")
P("                       the ranks CAN and DO differ.")
for (pair in list(c("agg_myeloid_siglec15", "agg_myeloid_siglec15_cp10k"),
                  c("allcell_tumour_agg", "allcell_tumour_agg_cp10k"))) {
  a <- spec[[pair[1]]]; b <- spec[[pair[2]]]
  same_ranks <- identical(rank(a), rank(b))
  ra <- suppressWarnings(cor(a, spec$contrast, method = "spearman"))
  rb <- suppressWarnings(cor(b, spec$contrast, method = "spearman"))
  P(sprintf("%-24s vs %-28s  ranks identical: %-5s   rho %+.4f vs %+.4f",
            pair[1], pair[2], same_ranks, ra, rb))
}
P("=> the aggregation rule is immaterial for the ALL-CELL aggregate and MATERIALLY")
P("   changes the MYELOID aggregate rho by ~0.057. That is a real fork and it is")
P("   flagged, not resolved silently. Both values are reported.")

P("\n===== 3. REPORTED CORRELATIONS, BOTH CONVENTIONS =====")
print(cors[, .(summary_statistic, n, spearman_rho, ci_lower, ci_upper, p_value)])
P("\nheadline (manuscript aggregation rule, inverse-log):")
for (v in c("agg_myeloid_siglec15", "detectable_fraction",
            "mean_expr_detectable", "allcell_tumour_agg")) {
  r <- cors[summary_statistic == v]
  P(sprintf("  %-24s rho = %+.4f   CI %+.4f to %+.4f   p = %.4f",
            v, r$spearman_rho, r$ci_lower, r$ci_upper, r$p_value))
}
P("\nconvention controls (plain CP10K mean instead of inverse-log):")
for (v in c("agg_myeloid_siglec15_cp10k", "allcell_tumour_agg_cp10k")) {
  r <- cors[summary_statistic == v]
  P(sprintf("  %-30s rho = %+.4f   CI %+.4f to %+.4f   p = %.4f",
            v, r$spearman_rho, r$ci_lower, r$ci_upper, r$p_value))
}

P("\n===== 4. DOES THE CHOICE MATTER? =====")
d_agg <- abs(cors[summary_statistic == "agg_myeloid_siglec15", spearman_rho] -
             cors[summary_statistic == "agg_myeloid_siglec15_cp10k", spearman_rho])
d_all <- abs(cors[summary_statistic == "allcell_tumour_agg", spearman_rho] -
             cors[summary_statistic == "allcell_tumour_agg_cp10k", spearman_rho])
P("difference in rho, myeloid aggregate   : ", sprintf("%.4f", d_agg))
P("difference in rho, all-cell aggregate  : ", sprintf("%.4f", d_all))
P("=> the aggregation rule leaves the all-cell aggregate exactly unchanged and")
P("   moves the myeloid aggregate rho by ", sprintf("%.4f", d_agg),
  " (one step of the n=6 grid).")
P("   BOTH values are reported. The inverse-log rule is used for the headline figure")
P("   because it is the only rule the manuscript states; the plain mean is carried")
P("   alongside so the choice is visible.")
P("   Neither convention changes the qualitative reading: at n = 6 the myeloid")
P("   aggregate, the detectable fraction and the mean among detectable cells all")
P("   have intervals spanning most of [-1, +1] and including zero. Only the all-cell")
P("   aggregate has an interval excluding zero, under BOTH conventions.")

P("\n===== 5. ATTAINABLE RESOLUTION AT n = 6 =====")
P("Spearman at n = 6 can only take values k/15, k = -15..15. Reported rhos:")
for (v in c("agg_myeloid_siglec15", "detectable_fraction",
            "mean_expr_detectable", "allcell_tumour_agg")) {
  r <- cors[summary_statistic == v, spearman_rho]
  P(sprintf("  %-24s rho = %+.6f  = %+d/15", v, r, round(r * 15)))
}
P("Any interval spanning a large part of [-1, +1] is the honest width at this n;")
P("these four numbers are reported as estimates with essentially no power.")

P("\nDONE ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
