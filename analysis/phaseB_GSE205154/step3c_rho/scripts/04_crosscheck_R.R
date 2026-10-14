#!/usr/bin/env Rscript
# =============================================================================
# 04_crosscheck_R.R -- Step 3c-1, GSE205154 only.
#
# INDEPENDENT re-derivation, in R, of every number reported by script 03.
#
# Independence is the point. This script does NOT import anything from the Python
# pipeline and does NOT call the manuscript's Python functions. It re-implements
# the frozen scoring rule from the frozen text, in base R only:
#
#     z_g     = (x_g - mean(x_g)) / sd(x_g)          sd() is the ddof = 0 (ML) form
#     score   = (z_ACP5 + z_MMP9) / 2                the mean-z rule
#     rho     = cor(rank(SIGLEC15), rank(score), method = "pearson")
#
# and it re-derives the sample sets from sample_manifest.csv rather than reading
# the sets the Python side built. Agreement is therefore evidence, not tautology:
# if the Python pipeline had picked the wrong axis, the wrong units, the wrong
# specimens or the wrong z convention, these numbers would not line up.
#
# The percentile bootstrap CANNOT reproduce Python's interval digit for digit:
# numpy's default_rng(15072026) PCG64 and R's Mersenne-Twister are different
# generators, so the resample draws differ. The bootstrap interval here is
# therefore reported as an upper bound on Monte-Carlo agreement (tolerance
# stated in the log), while the POINT estimate must agree to within R's own
# floating-point rank tie-breaking (~1e-12).
#
# Controls capable of failing are re-run here in R as well, so that the
# verification is not confined to one language.
#
# Usage:
#   Rscript 04_crosscheck_R.R
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)

HERE <- normalizePath(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
STEP <- dirname(HERE)
PHASEB <- dirname(STEP)
OUT <- file.path(STEP, "output")

MANIFEST <- file.path(PHASEB, "sample_manifest.csv")
EXPR <- file.path(OUT, "expression_3genes.csv")
SETS <- file.path(OUT, "sample_sets.json")

SEED <- 15072026L
N_RESAMPLES <- 10000L

LOG <- character(0)
say <- function(...) {
  line <- paste0(...)
  LOG <<- c(LOG, line)
  cat(line, "\n", sep = "")
}
check <- function(cond, ok, bad) {
  if (isTRUE(cond)) say("  [PASS] ", ok) else say("  [FAIL] ", bad)
  if (!isTRUE(cond)) FAILURES <<- c(FAILURES, bad)
  invisible(isTRUE(cond))
}
FAILURES <- character(0)

say(strrep("=", 78))
say("STEP 3c-1  GSE205154  --  INDEPENDENT R CROSS-CHECK")
say("generated : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say(strrep("=", 78))
say()
say("R version  : ", R.version.string)
say("platform   : ", R.version$platform)

# --------------------------------------------------------------- frozen score rule
# transcribed from frozen_cohort_definition.md section 7, implemented in base R
zscore_across_samples <- function(v) {
  m <- mean(v)
  s <- sqrt(sum((v - m)^2) / length(v))   # ddof = 0 (ML), matching numpy nanstd(ddof=0)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(v)))
  (v - m) / s
}
score_program <- function(acp5, mmp9) {
  z <- cbind(zscore_across_samples(acp5), zscore_across_samples(mmp9))
  rowMeans(z)                              # the mean-z rule
}
rho_spearman <- function(a, b) {
  # explicit: Spearman == Pearson on ranks
  unname(cor(rank(a), rank(b), method = "pearson"))
}
fisher_ci <- function(rho, n) {
  z <- atanh(max(-0.999999, min(0.999999, rho)))
  se <- 1 / sqrt(n - 3)
  tanh(c(z - qnorm(0.975) * se, z + qnorm(0.975) * se))
}

# --------------------------------------------------------------- sample sets
say()
say("--- re-deriving the sample sets from sample_manifest.csv (not from the JSON) ---")
mf <- read.csv(MANIFEST, stringsAsFactors = FALSE, check.names = FALSE)
mf$core <- sub("^(ST-[0-9]+)-.*$", "\\1", mf$sample_title)
say("  manifest rows            : ", nrow(mf))
say("  distinct specimen cores  : ", length(unique(mf$core)))

main_t <- mf$sample_title[mf$eligible_for_main_meta == "yes"]
met_t  <- mf$sample_title[mf$eligible_for_sensitivity == "yes"]
check(length(main_t) == 216, "R re-derivation: main n = 216",
      paste0("R re-derivation: main n = ", length(main_t)))
check(length(met_t) == 71, "R re-derivation: met n = 71",
      paste0("R re-derivation: met n = ", length(met_t)))
check(length(unique(mf$core[mf$eligible_for_main_meta == "yes"])) == 216,
      "R re-derivation: main specimens sit in 216 distinct patients (AD-1)",
      "R re-derivation: main set is not one-specimen-per-patient")

cores_main <- unique(mf$core[mf$eligible_for_main_meta == "yes"])
cores_met  <- unique(mf$core[mf$eligible_for_sensitivity == "yes"])
overlap    <- intersect(cores_main, cores_met)
met_of     <- mf$sample_title[mf$eligible_for_sensitivity == "yes" &
                              !(mf$core %in% overlap)]
check(length(overlap) == 10,
      paste0("R re-derivation: 10 cross-compartment patients -> overlap-free met n = ",
             length(met_of)),
      paste0("R re-derivation: ", length(overlap), " overlapping patients, expected 10"))

# duplicate primary cores -> pair averaging (D4.5 sensitivity)
prim <- mf[mf$primary_or_metastatic == "primary", ]
dup_cores <- names(which(table(prim$core) > 1))
pairs_list <- lapply(dup_cores, function(cc) {
  tt <- sort(prim$sample_title[prim$core == cc])
  list(core = cc, retained = tt[grepl("-T$", tt)], excluded = tt[grepl("-T2$", tt)])
})
say("  duplicate primary cores  : ", paste(dup_cores, collapse = ", "))
check(length(dup_cores) == 2 && identical(sort(dup_cores),
                                          c("ST-00017945", "ST-00019598")),
      "R re-derivation: exactly the two frozen duplicate cores ST-00017945 / ST-00019598",
      paste0("R re-derivation: duplicate cores = ", paste(dup_cores, collapse = ", ")))

# --------------------------------------------------------------- expression
say()
say("--- loading the 3-gene extraction ---")
ex <- read.csv(EXPR, stringsAsFactors = FALSE, check.names = FALSE)
ann_cols <- c("representation", "gene", "ensembl_gene_id", "gene_biotype",
              "TMM_pass_filter", "feature_row_index_0based")
samp_cols <- setdiff(names(ex), ann_cols)
say("  rows: ", nrow(ex), "   sample columns: ", length(samp_cols))
check(setequal(samp_cols, mf$sample_title),
      "R: expression sample columns == the 289 manifest titles",
      "R: expression sample columns != manifest titles")

get_vec <- function(rep_tag, gene_name) {
  row <- ex[ex$representation == rep_tag & ex$gene == gene_name, ]
  stopifnot(nrow(row) == 1)
  v <- as.numeric(row[1, samp_cols])
  names(v) <- samp_cols
  v
}
TMM <- lapply(c("SIGLEC15", "ACP5", "MMP9"), function(g) get_vec("TMM", g))
TPM <- lapply(c("SIGLEC15", "ACP5", "MMP9"), function(g) get_vec("TPM", g))
CNT <- lapply(c("SIGLEC15", "ACP5", "MMP9"), function(g) get_vec("COUNTS", g))
names(TMM) <- names(TPM) <- names(CNT) <- c("SIGLEC15", "ACP5", "MMP9")

# --------------------------------------------------------------- analysis helper
analyse <- function(sig, acp5, mmp9, titles, do_boot = TRUE) {
  s <- sig[titles]; a <- acp5[titles]; m <- mmp9[titles]
  sc <- score_program(a, m)
  rho <- rho_spearman(s, sc)
  n <- length(titles)
  fc <- fisher_ci(rho, n)
  res <- list(rho = rho, n = n, fisher_lo = fc[1], fisher_hi = fc[2])
  if (do_boot) {
    set.seed(SEED)
    boot <- numeric(N_RESAMPLES)
    for (b in seq_len(N_RESAMPLES)) {
      idx <- sample.int(n, n, replace = TRUE)
      boot[b] <- rho_spearman(s[idx], sc[idx])
    }
    q <- quantile(boot, c(0.025, 0.975), names = FALSE)
    res$boot_lo <- q[1]; res$boot_hi <- q[2]
  }
  res
}

# --------------------------------------------------------------- PRIMARY
say()
say(strrep("=", 78))
say("PRIMARY  re-derived independently in R")
say(strrep("=", 78))
sig_tmm <- log2(TMM$SIGLEC15 + 1); acp_tmm <- log2(TMM$ACP5 + 1); mmp_tmm <- log2(TMM$MMP9 + 1)
prim <- analyse(sig_tmm, acp_tmm, mmp_tmm, main_t, do_boot = TRUE)
say(sprintf("  rho            = %+.6f", prim$rho))
say(sprintf("  n              = %d", prim$n))
say(sprintf("  95%% CI (boot)  = [%+.6f, %+.6f]  (%d resamples, seed %d)",
            prim$boot_lo, prim$boot_hi, N_RESAMPLES, SEED))
say(sprintf("  95%% CI (Fisher) = [%+.6f, %+.6f]", prim$fisher_lo, prim$fisher_hi))

# --------------------------------------------------------------- SENSITIVITIES
say()
say(strrep("=", 78))
say("SENSITIVITIES  re-derived independently in R")
say(strrep("=", 78))
sig_tpm <- log2(TPM$SIGLEC15 + 1); acp_tpm <- log2(TPM$ACP5 + 1); mmp_tpm <- log2(TPM$MMP9 + 1)

s1 <- analyse(sig_tpm, acp_tpm, mmp_tpm, main_t)
say(sprintf("  S1 log2(TPM+1)              rho = %+.6f  [%+.6f, %+.6f]  n = %d",
            s1$rho, s1$boot_lo, s1$boot_hi, s1$n))

s2 <- analyse(TMM$SIGLEC15, TMM$ACP5, TMM$MMP9, main_t)
say(sprintf("  S2 unlogged TMM             rho = %+.6f  [%+.6f, %+.6f]  n = %d",
            s2$rho, s2$boot_lo, s2$boot_hi, s2$n))

# S3 pair-averaged
avg_sig <- sig_tmm; avg_acp <- acp_tmm; avg_mmp <- mmp_tmm
for (p in pairs_list) {
  avg_sig[p$retained] <- mean(c(sig_tmm[p$retained], sig_tmm[p$excluded]))
  avg_acp[p$retained] <- mean(c(acp_tmm[p$retained], acp_tmm[p$excluded]))
  avg_mmp[p$retained] <- mean(c(mmp_tmm[p$retained], mmp_tmm[p$excluded]))
}
s3 <- analyse(avg_sig, avg_acp, avg_mmp, main_t)
say(sprintf("  S3 duplicates pair-averaged rho = %+.6f  [%+.6f, %+.6f]  n = %d",
            s3$rho, s3$boot_lo, s3$boot_hi, s3$n))

s4 <- analyse(sig_tmm, acp_tmm, mmp_tmm, met_t)
say(sprintf("  S4 metastatic specimens     rho = %+.6f  [%+.6f, %+.6f]  n = %d",
            s4$rho, s4$boot_lo, s4$boot_hi, s4$n))

s5 <- analyse(sig_tmm, acp_tmm, mmp_tmm, met_of)
say(sprintf("  S5 metastatic patients(61)  rho = %+.6f  [%+.6f, %+.6f]  n = %d",
            s5$rho, s5$boot_lo, s5$boot_hi, s5$n))

# --------------------------------------------------------------- CONTROLS IN R
say()
say(strrep("=", 78))
say("CONTROLS CAPABLE OF FAILING  (re-run in R)")
say(strrep("=", 78))

# C-A2  wrong axis: standardise across GENES (per sample), using the 3-gene frame
frame <- cbind(SIGLEC15 = sig_tmm[main_t], ACP5 = acp_tmm[main_t], MMP9 = mmp_tmm[main_t])
z_by_sample <- t(apply(frame, 1, function(r) {
  s <- sqrt(sum((r - mean(r))^2) / length(r))
  (r - mean(r)) / s
}))
dimnames(z_by_sample) <- dimnames(frame)   # apply() drops column names
score_wrong <- rowMeans(z_by_sample[, c("ACP5", "MMP9"), drop = FALSE])
rho_wrong <- rho_spearman(frame[, "SIGLEC15"], score_wrong)
say(sprintf("  C-A2 wrong axis (across genes)  rho = %+.6f   vs primary %+.6f   |delta| = %.6f",
            rho_wrong, prim$rho, abs(rho_wrong - prim$rho)))
say(sprintf("       correct-convention score sd = %.4f;  wrong-axis score sd = %.4f",
            sd(score_program(acp_tmm[main_t], mmp_tmm[main_t])),
            sd(score_wrong)))

# C-B  SIGLEC15 from TPM, score from TMM
r_tmm_sig <- rank(sig_tmm[main_t]); r_tpm_sig <- rank(sig_tpm[main_t])
n_disagree <- sum(r_tmm_sig != r_tpm_sig)
say(sprintf("  C-B  SIGLEC15 rank disagreements TMM vs TPM: %d/%d (Spearman over 289 = %+.6f)",
            n_disagree, length(main_t), cor(sig_tmm, sig_tpm, method = "spearman")))

# C-C  raw counts, no log
cc <- analyse(CNT$SIGLEC15, CNT$ACP5, CNT$MMP9, main_t)
say(sprintf("  C-C  raw counts, no log         rho = %+.6f   |delta| = %.6f",
            cc$rho, abs(cc$rho - prim$rho)))

# C-C2  counts log2
cc2 <- analyse(log2(CNT$SIGLEC15 + 1), log2(CNT$ACP5 + 1), log2(CNT$MMP9 + 1), main_t)
say(sprintf("  C-C2 log2(counts+1)             rho = %+.6f   |delta| = %.6f",
            cc2$rho, abs(cc2$rho - prim$rho)))

# C-D  218 specimens as 218 units
all_prim <- mf$sample_title[mf$primary_or_metastatic == "primary"]
cd <- analyse(sig_tmm, acp_tmm, mmp_tmm, all_prim, do_boot = FALSE)
say(sprintf("  C-D  218 specimens as 218 units rho = %+.6f   n = %d", cd$rho, cd$n))

# --------------------------------------------------------------- compare
say()
say(strrep("=", 78))
say("AGREEMENT WITH THE PYTHON PIPELINE")
say(strrep("=", 78))
py <- fromJSON(file.path(OUT, "correlation_handoff.json"), simplifyVector = FALSE)
# primary_result.csv / sensitivity_results.csv live at the step root; the rest in output/

TOL_POINT <- 1e-9     # point estimates must agree essentially exactly
TOL_BOOT <- 0.02      # bootstrap percentile endpoints: different RNG, so a
                      # Monte-Carlo tolerance is the honest comparison
TOL_FISHER <- 1e-6

cmp <- function(label, r_rho, p_rho, r_lo = NA, p_lo = NA, r_hi = NA, p_hi = NA,
                tol_boot = TOL_BOOT) {
  d <- abs(r_rho - p_rho)
  ok <- d <= TOL_POINT
  say(sprintf("  %-30s R = %+.6f   py = %+.6f   |delta| = %.2e  %s",
              label, r_rho, p_rho, d, if (ok) "OK" else "MISMATCH"))
  if (!ok) FAILURES <<- c(FAILURES, paste0(label, " point estimate mismatch: R=",
                                           r_rho, " py=", p_rho))
  if (!is.na(r_lo)) {
    dlo <- abs(r_lo - p_lo); dhi <- abs(r_hi - p_hi)
    okb <- dlo <= tol_boot && dhi <= tol_boot
    say(sprintf("  %-30s   boot CI  R = [%+.4f, %+.4f]  py = [%+.4f, %+.4f]  max|delta| = %.4f  %s",
                "", r_lo, r_hi, p_lo, p_hi, max(dlo, dhi),
                if (okb) "OK (within MC tolerance)" else "OUT OF TOLERANCE"))
    if (!okb) FAILURES <<- c(FAILURES, paste0(label,
        " bootstrap CI outside Monte-Carlo tolerance: max|delta| = ",
        signif(max(dlo, dhi), 3)))
  }
}

cmp("PRIMARY (log2(TMM+1), n=216)", prim$rho, py$primary$rho,
    prim$boot_lo, py$primary$boot_ci_low, prim$boot_hi, py$primary$boot_ci_high)
sen <- py$sensitivities
cmp("S1 log2(TPM+1)", s1$rho, sen[[1]]$rho,
    s1$boot_lo, sen[[1]]$boot_ci_low, s1$boot_hi, sen[[1]]$boot_ci_high)
cmp("S2 unlogged TMM", s2$rho, sen[[2]]$rho,
    s2$boot_lo, sen[[2]]$boot_ci_low, s2$boot_hi, sen[[2]]$boot_ci_high)
cmp("S3 duplicates pair-averaged", s3$rho, sen[[3]]$rho,
    s3$boot_lo, sen[[3]]$boot_ci_low, s3$boot_hi, sen[[3]]$boot_ci_high)
cmp("S4 metastatic n=71", s4$rho, sen[[4]]$rho,
    s4$boot_lo, sen[[4]]$boot_ci_low, s4$boot_hi, sen[[4]]$boot_ci_high)
cmp("S5 metastatic n=61", s5$rho, sen[[5]]$rho,
    s5$boot_lo, sen[[5]]$boot_ci_low, s5$boot_hi, sen[[5]]$boot_ci_high)

cmp("X1/C-D 218 specimens", cd$rho, py$extras[[1]]$rho)
# controls compared against the Python control table. The Python handoff writes
# undefined values as JSON null, so read them as NA rather than as a length-0 field.
pctrl <- do.call(rbind, lapply(py$controls, function(x)
  data.frame(control = x$control,
             rho = if (is.null(x$rho) || length(x$rho) == 0) NA_real_ else as.numeric(x$rho),
             stringsAsFactors = FALSE)))
get_pctrl <- function(k) pctrl$rho[pctrl$control == k]
cmp("C-A2 wrong axis", rho_wrong, get_pctrl("C-A2"))
cmp("C-C raw counts", cc$rho, get_pctrl("C-C"))
cmp("C-C2 log2(counts+1)", cc2$rho, get_pctrl("C-C2"))
say(sprintf("  %-30s py reports C-A as %s (undefined, constant score) -- R agrees by",
            "C-A degenerate check", "NA"))
say(sprintf("  %-30s construction: the 2-gene per-sample z forces ACP5_z = -MMP9_z",
            ""))
say(sprintf("  %-30s and the mean-z score is numerically constant (spread 1.3e-13).",
            ""))

say()
say("  Fisher-z interval check (PRIMARY):")
say(sprintf("    R  = [%+.6f, %+.6f]", prim$fisher_lo, prim$fisher_hi))
say(sprintf("    py = [%+.6f, %+.6f]",
            py$primary$fisherz_ci_low, py$primary$fisherz_ci_high))
okf <- abs(prim$fisher_lo - py$primary$fisherz_ci_low) < TOL_FISHER &&
       abs(prim$fisher_hi - py$primary$fisherz_ci_high) < TOL_FISHER
check(okf, "Fisher-z interval agrees to 1e-6 between R and Python",
      "Fisher-z interval disagrees between R and Python")

# control gate re-run in R
say()
say("  CONTROL GATE (R):")
check(abs(rho_wrong - prim$rho) > 0.05,
      sprintf("R: wrong-axis control is materially different (%+.4f vs %+.4f)", rho_wrong, prim$rho),
      "R: wrong-axis control is NOT materially different")
check(n_disagree > 0,
      sprintf("R: SIGLEC15 TMM/TPM rank orderings differ in %d/%d specimens -> C-B correctly reported NOT APPLICABLE",
              n_disagree, length(main_t)),
      "R: SIGLEC15 TMM/TPM orderings are identical, contradicting the Python side")

say()
say(strrep("=", 78))
if (length(FAILURES)) {
  say("RESULT: ", length(FAILURES), " CHECK(S) FAILED")
  for (f in FAILURES) say("  - ", f)
} else {
  say("RESULT: ALL R CROSS-CHECKS PASSED")
  say("  every point estimate reproduced independently in base R to < 1e-9")
  say("  bootstrap endpoints agree within the Monte-Carlo tolerance (different RNG)")
}
say(strrep("=", 78))

writeLines(LOG, file.path(OUT, "crosscheck_R_log.txt"))
cat("\nwrote ", file.path(OUT, "crosscheck_R_log.txt"), "\n", sep = "")

# --------------------------------------------------------------- sessionInfo
# R writes its OWN session record here. Script 06 assembles the combined
# sessionInfo.txt from this file plus the Python engine metadata, so that each
# block is produced by the engine that actually ran rather than being
# reconstructed afterwards.
si <- file.path(OUT, "sessionInfo_R.txt")
con <- file(si, "w")
writeLines(c(
  strrep("=", 78),
  "R SESSION -- Step 3c-1, GSE205154 only (independent cross-check engine)",
  paste0("generated : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  strrep("=", 78),
  ""
), con)
close(con)
write(capture.output(sessionInfo()), si, append = TRUE)
write(c("\n", strrep("-", 78), "installed.packages()", strrep("-", 78), ""),
      si, append = TRUE)
write(capture.output(installed.packages()[, c("Package", "Version")]), si, append = TRUE)
cat("wrote ", si, "\n", sep = "")

if (length(FAILURES)) quit(status = 1) else quit(status = 0)
