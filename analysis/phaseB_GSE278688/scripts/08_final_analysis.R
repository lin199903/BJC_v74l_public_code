# =============================================================================
# 08_final_analysis.R -- Phase B step 1, GSE278688 (FROZEN RUN)
#
# FROZEN RULES, applied verbatim. Every convention below was READ FROM the
# manuscript's own pipeline (analysis_F_cross_scale.py,
# run_S20_state_controls.py, run_S21_ambient.py) and from the frozen brief --
# none was chosen by this script:
#
#  eligible tumour label : tissue == "Tumor"
#  cell type             : deposited all_celltype == "Myeloid" (no re-annotation)
#  >= 50 included tumour-myeloid cells per sample
#  SIGLEC15 detectability: stored log1p-normalised expression > 0
#  primary contrast      : needs BOTH detectable and non-detectable groups
#                          (NO >=2/>=5 minimum -- that is the separate
#                           restricted TREM2_TAM sensitivity rule)
#  programme score       : PROGRAM = ACP5 + MMP9. Each gene z-standardised
#                          WITHIN the biological sample, then averaged.
#                          zscore_columns() uses ddof=0 (population SD).
#                          A gene invariant within a sample contributes no
#                          information and is OMITTED for that sample; the
#                          primary score REQUIRES BOTH genes.
#                          SIGLEC15 NEVER enters the programme score.
#  sample contrasts      : equally weighted
#  aggregation of RNA    : inverse-log = log1p(mean(expm1(x)))
#  bootstrap             : dataset-stratified, 10000 resamples, seed 15072026
#  normality of the matrix: deposited values are RAW INTEGER UMI COUNTS.
#                          Detectability (>0) is invariant to this. The
#                          manuscript pipeline z-scores stored log1p(CP10K)
#                          values, so CP10K+log1p is reconstructed here.
#  all-cell tumour aggregate: mean of SIGLEC15 over EVERY cell in the
#                          retained tumour samples (manuscript SI section 15).
#
# SEED 15072026 (manuscript RNG seed)   N_RESAMPLES 10000
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

SEED <- 15072026L
N_RESAMPLES <- 10000L
set.seed(SEED)

ROOT <- "E:/Tex/HMlong/BJC_EXT_v1_20260927/phaseB_GSE278688"
OUT  <- file.path(ROOT,"output")
P <- function(...) cat(..., "\n", sep="")
P("SEED=", SEED, "  N_RESAMPLES=", N_RESAMPLES)
P("time=", format(Sys.time(), "%Y-%m-%d %H:%M:%S K"))

d    <- fread(file.path(OUT,"_myeloid_gene_cache.csv"))
allc <- fread(file.path(OUT,"_allcells_siglec15_cache.csv"))
P("myeloid cache : ", nrow(d), " x ", ncol(d))
P("all-cell cache: ", nrow(allc), " x ", ncol(allc))

## ===================== 1. SAMPLE MANIFEST ==================================
P("\n===== 1. SAMPLE MANIFEST =====")
agg <- d[, .(myeloid_cells = .N,
             sig_det    = sum(SIGLEC15 > 0),
             sig_non    = sum(SIGLEC15 == 0)), by = sample_id]
mn <- unique(d[, .(sample_id, tissue, patient)])
mn <- merge(mn, agg, by="sample_id", all.x=TRUE)
tc <- allc[, .(total_cells = .N), by=sample_id]
mn <- merge(mn, tc, by="sample_id", all.x=TRUE)
setnames(mn, c("sig_det","sig_non"),
             c("siglec15_detectable_myeloid","siglec15_nondetectable_myeloid"))
mn[, passes_min50_myeloid := myeloid_cells >= 50]
mn[, passes_both_groups   := passes_min50_myeloid &
                              siglec15_detectable_myeloid   >= 1 &
                              siglec15_nondetectable_myeloid >= 1]
mn[, eligible_tumour_label := tissue == "Tumor"]
mn[, exclusion_reason := ""]
mn[eligible_tumour_label == FALSE,
   exclusion_reason := sprintf("ineligible tissue label '%s' (not primary PDAC tumour)", tissue)]
mn[eligible_tumour_label == TRUE & passes_min50_myeloid == FALSE,
   exclusion_reason := sprintf("only %d tumour-myeloid cells (<50)", myeloid_cells)]
mn[eligible_tumour_label == TRUE & passes_min50_myeloid == TRUE & passes_both_groups == FALSE,
   exclusion_reason := sprintf("0 SIGLEC15 RNA-detectable tumour-myeloid cells (of %d); detectable group empty",
                               myeloid_cells)]
mn[, final_eligible := eligible_tumour_label == TRUE & passes_min50_myeloid == TRUE &
                        passes_both_groups == TRUE & exclusion_reason == ""]
setcolorder(mn, c("sample_id","patient","tissue","total_cells","myeloid_cells",
                  "siglec15_detectable_myeloid","siglec15_nondetectable_myeloid",
                  "passes_min50_myeloid","passes_both_groups","eligible_tumour_label",
                  "final_eligible","exclusion_reason"))
setorder(mn, -final_eligible, patient, tissue)
fwrite(mn, file.path(OUT,"sample_manifest.csv"))
P("manifest rows=", nrow(mn), "  final n=", sum(mn$final_eligible))
print(mn[, .(sample_id,patient,tissue,myeloid_cells,siglec15_detectable_myeloid,
             siglec15_nondetectable_myeloid,passes_min50_myeloid,passes_both_groups,
             final_eligible)], nrows=35)

P("\n-- EXCLUSION ARITHMETIC --")
P(sprintf("%-46s %d", "GSE278688 samples deposited", nrow(mn)))
P(sprintf("%-46s -%d", "  ineligible tissue label (!= Tumor)", sum(mn$eligible_tumour_label == FALSE)))
P(sprintf("%-46s %d", "  = tumour samples", sum(mn$eligible_tumour_label == TRUE)))
P(sprintf("%-46s -%d", "  <50 tumour-myeloid cells",
          sum(mn$eligible_tumour_label == TRUE & mn$passes_min50_myeloid == FALSE)))
P(sprintf("%-46s -%d", "  detectable group empty (sig_det==0)",
          sum(mn$eligible_tumour_label == TRUE & mn$passes_min50_myeloid == TRUE &
              mn$passes_both_groups == FALSE)))
P(sprintf("%-46s %d", "  = FINAL ESTIMABLE n", sum(mn$final_eligible)))

## ===================== 2. PRIMARY CONTRAST =================================
P("\n===== 2. PRIMARY CONTRAST =====")
elig <- mn[final_eligible == TRUE, sample_id]

## CP10K + log1p reconstruction. Detectability is unchanged by this
## (monotone, f(0)=0) -- verified in 05_normalisation_fork_test.py.
d[, c("SIGLEC15_cpm","ACP5_cpm","MMP9_cpm") :=
    lapply(.SD, function(v) log1p(v * 1e4 / libsize)),
  .SDcols = c("SIGLEC15","ACP5","MMP9")]

zpop <- function(x) {           # matches zscore_columns(values): ddof=0
  s <- stats::sd(x) * sqrt((length(x)-1)/length(x))   # population SD
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x)) / s
}

recs <- list()
for (s in elig) {
  i  <- which(d$sample_id == s)
  dd <- d[i]
  za <- zpop(dd$ACP5_cpm); zm <- zpop(dd$MMP9_cpm)
  inv <- c(if (all(is.na(za))) "ACP5", if (all(is.na(zm))) "MMP9")
  if (all(is.na(za)) || all(is.na(zm))) {
    ## primary score REQUIRES both genes
    recs[[s]] <- data.table(sample_id=s, patient=dd$patient[1],
      contrast=NA_real_, detectable_cells=sum(dd$SIGLEC15>0),
      nondetectable_cells=sum(dd$SIGLEC15==0),
      acp5_invariant=all(is.na(za)), mmp9_invariant=all(is.na(zm)),
      genes_in_score=NA_character_, status="EXCLUDED: programme gene invariant")
    next
  }
  base <- (za + zm)/2                                   # both genes present
  det  <- dd$SIGLEC15 > 0
  recs[[s]] <- data.table(sample_id=s, patient=dd$patient[1],
    contrast = mean(base[det]) - mean(base[!det]),
    detectable_cells = sum(det), nondetectable_cells = sum(!det),
    acp5_invariant=FALSE, mmp9_invariant=FALSE,
    genes_in_score="ACP5;MMP9", status="OK")
}
res <- rbindlist(recs)
res[, direction_positive := contrast > 0]
setorder(res, -contrast, na.last=TRUE)
fwrite(res, file.path(OUT,"primary_contrast_by_sample.csv"))
print(res, nrows=30)

## single-gene fallback, reported for transparency (NOT the primary score)
fb <- rbindlist(lapply(elig, function(s){
  dd <- d[which(d$sample_id==s)]
  za <- zpop(dd$ACP5_cpm); zm <- zpop(dd$MMP9_cpm)
  det <- dd$SIGLEC15>0
  if (all(is.na(za)) && all(is.na(zm))) return(NULL)
  v <- if (all(is.na(za))) zm else if (all(is.na(zm))) za else NULL
  if (is.null(v)) return(NULL)
  data.table(sample_id=s, single_gene_contrast = mean(v[det])-mean(v[!det]))
}))
if (nrow(fb)) { P("\n-- single-gene fallback (transparency only, NOT primary) --"); print(fb) }

x  <- res$contrast[is.finite(res$contrast)]
n  <- length(x)
P("\nsamples with an estimable primary contrast: ", n, " of ", nrow(res))

## ===================== 3. SUMMARY ==========================================
P("\n===== 3. SUMMARY OF THE PRIMARY CONTRAST =====")
P("n                  : ", n)
P("mean contrast      : ", sprintf("%.4f", mean(x)))
P("sd                 : ", sprintf("%.4f", stats::sd(x)))
P("median             : ", sprintf("%.4f", stats::median(x)))
P("range              : ", sprintf("%.4f to %.4f", min(x), max(x)))
P("positive direction : ", sum(x > 0), " / ", n)
P("negative direction : ", sum(x < 0), " / ", n)

## dataset-stratified bootstrap: ONE stratum -> resample samples within it
set.seed(SEED)
boot <- replicate(N_RESAMPLES, mean(sample(x, n, replace=TRUE)))
ci   <- unname(stats::quantile(boot, c(0.025, 0.975), type=7))
P("bootstrap          : ", N_RESAMPLES, " resamples, seed ", SEED,
  ", resampling unit = sample within the single dataset stratum")
P("95% CI             : ", sprintf("%.4f to %.4f", ci[1], ci[2]))

## signed-rank, matching the manuscript (scipy zero_method='wilcox')
wr  <- suppressWarnings(stats::wilcox.test(x, mu=0, alternative="two.sided",
                                           correct=FALSE, exact=FALSE))
wr0 <- suppressWarnings(stats::wilcox.test(x, mu=0, alternative="two.sided",
                                           correct=TRUE,  exact=FALSE))
P("signed-rank p (no continuity corr, = scipy default) : ",
  format.pval(wr$p.value, digits=4), "  V=", wr$statistic)
P("signed-rank p (continuity-corrected)                : ",
  format.pval(wr0$p.value, digits=4))
P("sign test  p                                        : ",
  format.pval(stats::binom.test(sum(x>0), n, 0.5)$p.value, digits=4))

sm <- data.table(dataset="GSE278688", stratum="GSE278688 (single stratum)",
  n_samples=n, mean_contrast=mean(x), sd_contrast=stats::sd(x),
  median_contrast=stats::median(x), min_contrast=min(x), max_contrast=max(x),
  score="ACP5;MMP9 (both required, z-standardised within sample, ddof=0)",
  normalisation="CP10K then log1p (reconstructed; deposit is raw UMI counts)",
  n_positive=sum(x>0), n_negative=sum(x<0),
  boot_resamples=N_RESAMPLES, boot_seed=SEED, ci_lower=ci[1], ci_upper=ci[2],
  wilcoxon_V=as.numeric(wr$statistic), wilcoxon_p=wr$p.value,
  sign_test_p=stats::binom.test(sum(x>0), n, 0.5)$p.value)
fwrite(sm, file.path(OUT,"primary_summary.csv"))
print(sm)

## ===================== 4. FOUR SPECIMEN-LEVEL SUMMARIES ====================
P("\n===== 4. FOUR SAME-SAMPLE SPECIMEN-LEVEL SIGLEC15 SUMMARIES =====")
retained <- elig                       # retained tumour samples = eligible set
spec <- rbindlist(lapply(elig, function(s){
  dd  <- d[which(d$sample_id==s)]
  ac  <- allc[sample_id == s]          # EVERY cell in the retained tumour sample
  det <- dd$SIGLEC15 > 0
  data.table(sample_id=s,
    agg_myeloid_siglec15 = mean(dd$SIGLEC15_cpm),                     # CP10K mean
    detectable_fraction  = mean(det),
    mean_expr_detectable = if (any(det)) mean(dd$SIGLEC15_cpm[det]) else NA_real_,
    allcell_tumour_agg   = mean(ac$SIGLEC15 * 1e4 / ac$libsize),
    n_all_cells          = nrow(ac),
    n_detectable_allcell = sum(ac$SIGLEC15 > 0))
}))
spec <- merge(spec, res[, .(sample_id, contrast)], by="sample_id")
setorder(spec, -contrast)
fwrite(spec, file.path(OUT,"specimen_level_summaries.csv"))
print(spec)

## --- Spearman with dataset-stratified sample bootstrap (manuscript method) --
## Single stratum => resample samples within GSE278688.
spearman_boot <- function(xv, yv) {
  ok <- is.finite(xv) & is.finite(yv); xv <- xv[ok]; yv <- yv[ok]; m <- length(xv)
  if (m < 4) return(list(rho=NA_real_, p=NA_real_, lo=NA_real_, hi=NA_real_, n=m))
  rho <- suppressWarnings(stats::cor(xv, yv, method="spearman"))
  pt  <- suppressWarnings(stats::cor.test(xv, yv, method="spearman", exact=FALSE))
  set.seed(SEED)
  bs <- vapply(seq_len(N_RESAMPLES), function(b){
    i <- sample.int(m, m, replace=TRUE)
    if (length(unique(xv[i])) < 2 || length(unique(yv[i])) < 2) return(NA_real_)
    suppressWarnings(stats::cor(xv[i], yv[i], method="spearman"))
  }, numeric(1))
  q <- unname(stats::quantile(bs, c(0.025,0.975), na.rm=TRUE, type=7))
  list(rho=unname(rho), p=pt$p.value, lo=q[1], hi=q[2], n=m)
}
cors <- rbindlist(lapply(c("agg_myeloid_siglec15","detectable_fraction",
                           "mean_expr_detectable","allcell_tumour_agg"), function(v){
  r <- spearman_boot(spec[[v]], spec$contrast)
  data.table(summary_statistic=v, n=r$n, spearman_rho=r$rho,
             ci_lower=r$lo, ci_upper=r$hi, p_value=r$p,
             ci_method=sprintf("dataset-stratified sample bootstrap, %d resamples, seed %d",
                               N_RESAMPLES, SEED))
}))
fwrite(cors, file.path(OUT,"specimen_level_correlations.csv"))
P("\n-- Spearman correlations with the sample-specific contrast --")
print(cors)

P("\nDONE ", format(Sys.time(), "%Y-%m-%d %H:%M:%S K"))
