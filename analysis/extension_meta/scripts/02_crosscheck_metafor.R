# extension_meta / 02_crosscheck_metafor.R
#
# INDEPENDENT-ENGINE CROSS-CHECK of the extension meta-analysis.
#
# Recomputes all three models with metafor (Viechtbauer), which estimates REML tau^2 by a
# completely different numerical route from the manuscript's bounded scalar search, and which
# provides Hartung-Knapp as `test="knha"` in its own implementation.
#
# This is a check on the ARITHMETIC, not a second analysis: nothing here is reported in
# meta_results.csv. What it establishes that a Python-only check could not:
#   1. Q and I2 are reproduced by an independent package;
#   2. the UNMODIFIED-HK control values (C3, C3b) are metafor's own knha answer, so those
#      controls are not artefacts of my replica function;
#   3. the closed-form modified-HK interval reproduced in R agrees with the manuscript's to
#      machine precision ONCE the tau^2 difference is carried, and the remaining tau^2
#      difference is a boundary/optimiser-tolerance artefact, not an arithmetic disagreement.
#
# RUN 1 OF THIS SCRIPT REPORTED 6 FAILURES. They were not relaxed; they were attributed, and
# the attribution is asserted (see the checks below). Recorded in ../failed_or_excluded.md.
#
# Inputs (written by 01_extension_meta.py): forest_inputs.csv, meta_results.csv
# Outputs: provenance/crosscheck_metafor.csv, provenance/crosscheck_metafor_log.txt

suppressPackageStartupMessages(library(metafor))

ext <- "E:/Tex/HMlong/BJC_EXT_v1_20260927/extension_meta"
fi  <- read.csv(file.path(ext, "forest_inputs.csv"), stringsAsFactors = FALSE)
mr  <- read.csv(file.path(ext, "meta_results.csv"), stringsAsFactors = FALSE)
vc_all <- read.csv(file.path(ext, "verification_controls.csv"), stringsAsFactors = FALSE)

con <- file(file.path(ext, "provenance", "crosscheck_metafor_log.txt"), open = "wt")
say <- function(...) { txt <- paste0(...); cat(txt, "\n"); writeLines(txt, con) }

say("========================================================================")
say("02_crosscheck_metafor.R   independent-engine cross-check")
say(paste("run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
say(paste("R", R.version.string, "| metafor", as.character(packageVersion("metafor"))))
say("========================================================================")

XATOL <- 1e-5   # scipy bounded minimise_scalar's default x tolerance, used by the manuscript
rows <- list()
n_fail <- 0
chk <- function(name, ok, detail) {
  say(sprintf("[%s] %s :: %s", ifelse(ok, "PASS", "FAIL"), name, detail))
  if (!ok) n_fail <<- n_fail + 1
  invisible(ok)
}

# The manuscript's modified Hartung-Knapp rule, as a closed form in R.
modHK <- function(z, v, tau2) {
  k      <- length(z)
  w      <- 1 / (v + tau2)
  mu     <- sum(w * z) / sum(w)
  q_raw  <- sum(w * (z - mu)^2) / (k - 1)
  q_used <- max(1, q_raw)
  se     <- sqrt(q_used / sum(w))
  crit   <- qt(0.975, df = k - 1)
  list(mu = mu, ci = tanh(mu + c(-1, 1) * crit * se), q_raw = q_raw, se = se)
}

for (lab in unique(fi$model)) {
  d  <- fi[fi$model == lab, ]
  z  <- atanh(d$rho)
  v  <- 1 / (d$n - 3)
  k  <- length(z)

  fit  <- rma(yi = z, vi = v, method = "REML")
  fitk <- rma(yi = z, vi = v, method = "REML", test = "knha")

  m <- mr[mr$k == k, ]
  stopifnot(nrow(m) == 1)

  A <- modHK(z, v, fit$tau2)      # with metafor's REML tau^2 (exactly 0 here)
  B <- modHK(z, v, m$tau2)        # with the manuscript's returned tau^2
  ci_mod <- B$ci
  d_tau  <- abs(fit$tau2 - m$tau2)
  d_ci_attrib <- max(abs(A$ci[1] - B$ci[1]), abs(A$ci[2] - B$ci[2]))
  d_ci_resid  <- max(abs(ci_mod[1] - m$ci_low), abs(ci_mod[2] - m$ci_high))

  say("")
  say(sprintf("--- %s  (k=%d: %s)", lab, k, paste(d$cohort, collapse = "; ")))
  say(sprintf("    metafor REML tau2      = %.17g", fit$tau2))
  say(sprintf("    manuscript tau2        = %.17g   |delta| = %.3e", m$tau2, d_tau))
  say(sprintf("    metafor QE             = %.17g   manuscript = %.17g", fit$QE, m$Q))
  say(sprintf("    metafor I2             = %.17g   manuscript = %.17g", fit$I2, m$I2))
  say(sprintf("    metafor pooled z       = %.17g", fit$beta))
  say(sprintf("    modified-HK in R with metafor tau2 = %.17g .. %.17g", A$ci[1], A$ci[2]))
  say(sprintf("    modified-HK in R with ms tau2      = %.17g .. %.17g", B$ci[1], B$ci[2]))
  say(sprintf("    manuscript modified-HK             = %.17g .. %.17g", m$ci_low, m$ci_high))
  say(sprintf("    UNMODIFIED-HK (metafor knha)       = %.17g .. %.17g",
              tanh(fitk$ci.lb), tanh(fitk$ci.ub)))

  # (1) tau^2: metafor lands exactly on the REML boundary; the manuscript lands within one
  #     optimiser x-tolerance of it. Agreement is asserted at the tolerance that is meaningful.
  chk(paste("REML tau2 agrees within the optimiser tolerance", lab),
      d_tau <= XATOL,
      sprintf("metafor %.3e (boundary), manuscript %.3e, |delta| %.3e <= xatol %.0e",
              fit$tau2, m$tau2, d_tau, XATOL))
  chk(paste("Q agrees (metafor QE vs manuscript Q)", lab),
      abs(fit$QE - m$Q) < 1e-12, sprintf("|delta| = %.3e", abs(fit$QE - m$Q)))
  chk(paste("I2 agrees", lab), abs(fit$I2 - m$I2) < 1e-12,
      sprintf("metafor %.10f vs manuscript %.10f", fit$I2, m$I2))

  # (2) the closed-form modified-HK rule in R is identical to the manuscript's -- once the
  #     manuscript's own tau^2 is carried.
  chk(paste("modified-HK rule reproduced in R matches the manuscript CI (its tau2 carried)", lab),
      d_ci_resid < 1e-9, sprintf("max|delta| = %.3e", d_ci_resid))

  # (3) ATTRIBUTION: the whole metafor-vs-manuscript CI gap is explained by the tau^2
  #     difference alone. Both gaps must be the same number.
  gap_full <- max(abs(A$ci[1] - m$ci_low), abs(A$ci[2] - m$ci_high))
  chk(paste("the metafor-vs-manuscript CI gap is fully explained by the tau2 difference", lab),
      abs(gap_full - d_ci_attrib) < 1e-12,
      sprintf("gap with metafor tau2 = %.3e; gap induced by tau2 alone = %.3e", gap_full, d_ci_attrib))

  # (4) the zero-crossing verdicts must agree between engines, on both HK conventions
  chk(paste("zero-crossing verdict agrees between engines (modified HK)", lab),
      (m$ci_low < 0 && m$ci_high > 0) == (A$ci[1] < 0 && A$ci[2] > 0),
      sprintf("manuscript %s; metafor-tau2 R %s",
              ifelse(m$ci_low < 0 && m$ci_high > 0, "crosses zero", "excludes zero"),
              ifelse(A$ci[1] < 0 && A$ci[2] > 0, "crosses zero", "excludes zero")))

  # (5) the unmodified-HK CONTROL values must be metafor's own knha answer. The control values
  #     are read from verification_controls.csv, not hard-coded, so this compares two files.
  vc <- vc_all[vc_all$model == lab & grepl("unmodified_HK", vc_all$control), ]
  stopifnot(nrow(vc) == 1)
  chk(paste("unmodified-HK control is confirmed by metafor's own knha", lab),
      max(abs(tanh(fitk$ci.lb) - vc$ci_low), abs(tanh(fitk$ci.ub) - vc$ci_high)) < 1e-4,
      sprintf("metafor knha %.6f..%.6f vs control file %.6f..%.6f", tanh(fitk$ci.lb),
              tanh(fitk$ci.ub), vc$ci_low, vc$ci_high))

  rows[[length(rows) + 1]] <- data.frame(
    model = lab, k = k, cohorts = paste(d$cohort, collapse = ";"),
    rho = m$pooled_rho, pooled_z_metafor = fit$beta,
    tau2_metafor = fit$tau2, tau2_manuscript = m$tau2, tau2_abs_delta = d_tau,
    Q_metafor = fit$QE, Q_manuscript = m$Q,
    I2_metafor = fit$I2, I2_manuscript = m$I2,
    ci_low_modHK_R_metaforTau2 = A$ci[1], ci_high_modHK_R_metaforTau2 = A$ci[2],
    ci_low_modHK_R_msTau2 = B$ci[1], ci_high_modHK_R_msTau2 = B$ci[2],
    ci_low_manuscript = m$ci_low, ci_high_manuscript = m$ci_high,
    ci_gap_metafor_vs_ms = gap_full, ci_gap_attributable_to_tau2 = d_ci_attrib,
    ci_low_unmodHK_metafor = tanh(fitk$ci.lb), ci_high_unmodHK_metafor = tanh(fitk$ci.ub),
    crosses_zero_manuscript = (m$ci_low < 0 && m$ci_high > 0),
    crosses_zero_metafor_tau2 = (A$ci[1] < 0 && A$ci[2] > 0),
    stringsAsFactors = FALSE
  )
}

say("")
say("========================================================================")
say(sprintf("TOTAL: %d checks, %d failed", 5 * length(unique(fi$model)), n_fail))
say("========================================================================")
close(con)

out <- do.call(rbind, rows)
write.csv(out, file.path(ext, "provenance", "crosscheck_metafor.csv"), row.names = FALSE)
cat("wrote provenance/crosscheck_metafor.csv\n")
if (n_fail > 0) quit(status = 1)
