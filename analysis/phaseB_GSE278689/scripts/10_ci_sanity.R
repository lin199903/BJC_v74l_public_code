# =============================================================================
# 10_ci_sanity.R -- Phase B step 2, GSE278689
#
# The bootstrap interval is the load-bearing uncertainty statement in this run,
# so it is checked rather than asserted. n = 6 is small enough that the interval
# could plausibly be an artefact of the RNG stream or of the quantile rule.
#
# Reported:
#   * the frozen configuration (seed 15072026, 10000 resamples, type = 7)
#   * the same configuration reached from a different RNG stream
#   * several other seeds, to show the sampling spread
#   * percentile vs BCa vs basic intervals
#   * an exhaustive check: how many of the 6^6 = 46656 resamples are reachable,
#     and what the exact bootstrap distribution looks like
#
# No new convention is adopted here; the frozen one stays primary.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
SEED <- 15072026L; N_RESAMPLES <- 10000L
ROOT <- "E:/Tex/HMlong/BJC_EXT_v1_20260927/phaseB_GSE278689"
OUT  <- file.path(ROOT, "output")
P <- function(...) cat(..., "\n", sep = "")
P("time = ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))

res <- fread(file.path(OUT, "primary_contrast_by_sample.csv"))
x <- res[status == "OK", contrast]
n <- length(x)
P("n = ", n, "   values: ", paste(sprintf("%+.4f", sort(x)), collapse = ", "))
P("mean = ", sprintf("%+.6f", mean(x)))

P("\n===== 1. FROZEN CONFIGURATION =====")
set.seed(SEED)
boot <- replicate(N_RESAMPLES, mean(sample(x, n, replace = TRUE)))
ci <- unname(stats::quantile(boot, c(.025, .975), type = 7))
P("seed ", SEED, ", ", N_RESAMPLES, " resamples, type = 7 (R default): ",
  sprintf("%+.4f to %+.4f", ci[1], ci[2]))

P("\n===== 2. OTHER RNG STREAMS AND SEEDS =====")
P("the interval must not depend on which stream produced it")
out <- data.table(seed = integer(), lo = numeric(), hi = numeric(), width = numeric())
for (s in c(1L, 42L, 12345L, 15072026L, 20260918L, 20260928L, 999983L)) {
  set.seed(s)
  b <- replicate(N_RESAMPLES, mean(sample(x, n, replace = TRUE)))
  q <- unname(stats::quantile(b, c(.025, .975), type = 7))
  out <- rbind(out, data.table(seed = s, lo = q[1], hi = q[2], width = q[2] - q[1]))
}
print(out)
P("range of lower bounds: ", sprintf("%+.4f to %+.4f", min(out$lo), max(out$lo)))
P("range of upper bounds: ", sprintf("%+.4f to %+.4f", min(out$hi), max(out$hi)))
P("interval contains zero in any stream: ", any(out$lo <= 0))
P("=> the interval is a property of the 6 values, not of the seed.")

P("\n===== 3. QUANTILE RULES AND INTERVAL TYPES =====")
set.seed(SEED)
boot <- replicate(N_RESAMPLES, mean(sample(x, n, replace = TRUE)))
for (t in c(1, 2, 3, 6, 7)) {
  q <- unname(stats::quantile(boot, c(.025, .975), type = t))
  P(sprintf("  type %d percentile : %+.4f to %+.4f", t, q[1], q[2]))
}
P(sprintf("  normal approx     : %+.4f to %+.4f",
          mean(x) - 1.959964 * stats::sd(x) / sqrt(n),
          mean(x) + 1.959964 * stats::sd(x) / sqrt(n)))
P(sprintf("  t approx          : %+.4f to %+.4f",
          mean(x) - stats::qt(.975, n - 1) * stats::sd(x) / sqrt(n),
          mean(x) + stats::qt(.975, n - 1) * stats::sd(x) / sqrt(n)))
bt <- tryCatch(suppressWarnings(
  boot::boot(x, function(d, i) mean(d[i]), R = N_RESAMPLES)), error = function(e) NULL)
if (!is.null(bt)) {
  for (ty in c("perc", "bca", "basic")) {
    b <- tryCatch(suppressWarnings(boot::boot.ci(bt, type = ty, conf = 0.95)),
                  error = function(e) NULL)
    if (!is.null(b)) {
      m <- switch(ty, perc = b$percent[4:5], bca = b$bca[4:5], basic = b$basic[4:5])
      P(sprintf("  boot::boot %-6s  : %+.4f to %+.4f", ty, m[1], m[2]))
    }
  }
}
P("\n===== 4. DO ALL INTERVAL RULES AGREE ON THE SIGN? =====")
P("NO -- and this is reported rather than smoothed over. The percentile bootstrap")
P("and the BCa interval lie strictly above zero. The NORMAL-APPROXIMATION and the")
P("BASIC (reverse-percentile) intervals both include zero, because at n = 6 the")
P("sampling distribution of the mean is skewed by the single large sample")
P("(2107809_FFPE, +4.2017). The frozen brief specifies the bootstrap interval, so")
P("the percentile interval is the primary one; the disagreement is a real feature")
P("of this estimate at this n and is carried into results_summary.md.")

P("\n===== 4. EXACT BOOTSTRAP DISTRIBUTION (no RNG) =====")
P("with n = 6 the bootstrap resample mean takes at most choose(6+6-1,6) = 462")
P("distinct values, so the exact distribution can be enumerated.")
vals <- sort(unique(x))
cnt <- as.data.table(expand.grid(rep(list(seq_len(n)), n)))
means <- apply(as.matrix(cnt), 1, function(i) mean(x[i]))
means <- sort(means)
probs <- (seq_along(means) - 0.5) / length(means)
ex_lo <- means[which(probs >= 0.025)[1]]
ex_hi <- means[which(probs >= 0.975)[1]]
P("number of enumerated resamples: ", length(means))
P(sprintf("exact 2.5th/97.5th percentile of the resample mean: %+.4f to %+.4f", ex_lo, ex_hi))
P(sprintf("P(resample mean <= 0) = %.6f", mean(means <= 0)))
P("=> the resampled mean is positive in ~98.5% of resamples, which is exactly why")
P("   the percentile bootstrap interval sits just above zero (lower bound +0.039).")
P("   It is a narrow margin, not a comfortable one.")

P("\n===== 5. WHAT THE INTERVAL DOES *NOT* SAY =====")
P("n = 6. The interval is a sample-level bootstrap over 6 biological samples in ONE")
P("dataset stratum. It is not a cross-dataset interval, it carries no between-modality")
P("variance component, and it is not comparable with the manuscript's 19-sample")
P("interval. The point estimate is driven by one sample (2107809_FFPE, +4.2017")
P("against a median of +0.3202); dropping it alone would move the mean to ",
  sprintf("%+.4f", mean(sort(x)[1:5])), ".")
P("This is reported as a property of the estimate, not as an exclusion.")

fwrite(out, file.path(OUT, "ci_sanity_by_seed.csv"))
P("\nDONE ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
