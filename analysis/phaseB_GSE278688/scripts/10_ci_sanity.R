# 10_ci_sanity.R -- check the bootstrap CI boundary of exactly -1 reported for
# the aggregated-myeloid correlation, so that an artefact is not mistaken for a result.
suppressPackageStartupMessages({ library(data.table) })
ROOT <- "E:/Tex/HMlong/BJC_EXT_v1_20260927/phaseB_GSE278688"
sp <- fread(file.path(ROOT,"output","specimen_level_summaries.csv"))
for (v in c("agg_myeloid_siglec15","detectable_fraction","mean_expr_detectable","allcell_tumour_agg")) {
  x <- sp[[v]]; y <- sp$contrast; n <- length(x)
  set.seed(15072026)
  bs <- replicate(10000, {
    i <- sample.int(n, n, replace=TRUE)
    if (length(unique(x[i])) < 2 || length(unique(y[i])) < 2) return(NA_real_)
    suppressWarnings(cor(x[i], y[i], method="spearman"))
  })
  obs <- suppressWarnings(cor(x, y, method="spearman"))
  q <- unname(quantile(bs, c(0.025, 0.975), na.rm=TRUE, type=7))
  cat(sprintf("%-22s obs=%+.4f  boot_min=%+.4f boot_max=%+.4f  q2.5=%+.4f q97.5=%+.4f  n_at_-1=%d  n_NA=%d\n",
              v, obs, min(bs,na.rm=TRUE), max(bs,na.rm=TRUE), q[1], q[2],
              sum(bs == -1, na.rm=TRUE), sum(is.na(bs))))
}
cat("\nNOTE: with n=8 the attainable Spearman rho grid is coarse (k/28).")
cat("\n      A resample can therefore land on exactly -1.0, which is a property")
cat("\n      of the small sample and the bootstrap, not an anomaly.\n")
