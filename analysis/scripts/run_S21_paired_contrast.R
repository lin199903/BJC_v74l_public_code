## Post hoc paired contrast for S21 depth strata.
## Uses only the archived sample-level S21 output; the patient sample remains
## the inferential unit. The bootstrap resamples paired samples within dataset.

set.seed(15072026)

input_path <- file.path("deposit_v24", "results", "S21_depth_strata.csv")
output_pairs <- file.path("deposit_v24", "results", "S21_depth_strata_paired.csv")
output_summary <- file.path("deposit_v24", "results", "S21_depth_strata_paired_summary.csv")

x <- read.csv(input_path, check.names = FALSE)
wide <- reshape(
  x[, c("dataset", "sample_id", "stratum", "difference_positive_minus_negative")],
  idvar = c("dataset", "sample_id"), timevar = "stratum", direction = "wide"
)
names(wide) <- sub("difference_positive_minus_negative\\.", "", names(wide))
wide <- wide[complete.cases(wide[, c("deep", "shallow")]), ]
wide$deep_minus_shallow <- wide$deep - wide$shallow

datasets <- unique(wide$dataset)
boot_mean <- replicate(10000, {
  sampled <- do.call(rbind, lapply(datasets, function(d) {
    z <- wide[wide$dataset == d, , drop = FALSE]
    z[sample(seq_len(nrow(z)), nrow(z), replace = TRUE), , drop = FALSE]
  }))
  mean(sampled$deep_minus_shallow)
})

ci <- unname(quantile(boot_mean, c(0.025, 0.975), type = 7))
p_value <- wilcox.test(wide$deep_minus_shallow, mu = 0, exact = FALSE)$p.value

summary <- data.frame(
  contrast = "Deep minus shallow cell-level effect",
  n_paired_samples = nrow(wide),
  mean_difference = mean(wide$deep_minus_shallow),
  ci_low = ci[1],
  ci_high = ci[2],
  wilcoxon_p = p_value,
  bootstrap_resamples = 10000,
  seed = 15072026
)

write.csv(wide, output_pairs, row.names = FALSE)
write.csv(summary, output_summary, row.names = FALSE)
print(summary)
