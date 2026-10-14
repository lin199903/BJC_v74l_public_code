# 11_allcell_convention.R -- the all-cell tumour aggregate is defined in the
# manuscript SI as "every cell in the retained tumour samples", but SI section 15
# states the aggregation only for the myeloid aggregate (inverse-log).
# This reports BOTH candidate conventions so the choice is visible and auditable.
suppressPackageStartupMessages({ library(data.table) })
ROOT <- "E:/Tex/HMlong/BJC_EXT_v1_20260927/phaseB_GSE278688"
OUT  <- file.path(ROOT,"output"); SEED <- 15072026L; N <- 10000L
sp   <- fread(file.path(OUT,"specimen_level_summaries.csv"))
allc <- fread(file.path(OUT,"_allcells_siglec15_cache.csv"))
elig <- sp$sample_id

ac <- rbindlist(lapply(elig, function(s){
  g <- allc[sample_id == s]
  cpm <- g$SIGLEC15 * 1e4 / g$libsize          # back-transformed CP10K
  data.table(sample_id = s,
             plain_cpm_mean    = mean(cpm),                 # plain mean, CP10K space
             inverse_log_cpm   = expm1(mean(log1p(cpm))))   # log1p(mean(expm1(.)))
}))
ac <- merge(ac, sp[, .(sample_id, contrast)], by="sample_id")

sb <- function(x, y){
  n <- length(x); set.seed(SEED)
  rho <- suppressWarnings(cor(x, y, method="spearman"))
  p   <- suppressWarnings(cor.test(x, y, method="spearman", exact=FALSE)$p.value)
  bs  <- replicate(N, { i <- sample.int(n,n,replace=TRUE)
    if (length(unique(x[i]))<2 || length(unique(y[i]))<2) return(NA_real_)
    suppressWarnings(cor(x[i],y[i],method="spearman")) })
  q <- unname(quantile(bs, c(.025,.975), na.rm=TRUE, type=7))
  c(rho=unname(rho), p=p, lo=q[1], hi=q[2])
}
cat("=== all-cell tumour aggregate: the two candidate conventions ===\n")
r1 <- sb(ac$plain_cpm_mean,  ac$contrast)
r2 <- sb(ac$inverse_log_cpm, ac$contrast)
cat(sprintf("%-20s rho=%+.4f  p=%.4f  95%% CI %+.4f to %+.4f\n",
            "plain CP10K mean", r1[1], r1[2], r1[3], r1[4]))
cat(sprintf("%-20s rho=%+.4f  p=%.4f  95%% CI %+.4f to %+.4f\n",
            "inverse-log", r2[1], r2[2], r2[3], r2[4]))
cat("\nSpearman over samples depends only on RANKS. The two conventions are")
cat("\nmonotone transforms of one another, so the point estimate is identical:\n")
cat(sprintf("  identical point estimate: %s\n", isTRUE(all.equal(r1[1], r2[1]))))
cat("Only the bootstrap interval differs slightly, because resampling with")
cat("replacement can break the monotone relation.\n")

## the same logic applies to the myeloid aggregate: check rank equality
myc <- rbindlist(lapply(elig, function(s){
  g <- fread(file.path(OUT,"_myeloid_gene_cache.csv"))[sample_id == s]
  cpm <- g$SIGLEC15 * 1e4 / g$libsize
  data.table(sample_id=s, plain=mean(cpm), invlog=expm1(mean(log1p(cpm))))
}))
cat(sprintf("\nmyeloid aggregate, rank correlation between conventions: %+.6f\n",
            cor(myc$plain, myc$invlog, method="spearman")))
cat(sprintf("myeloid aggregate, Pearson  correlation between conventions: %+.6f\n",
            cor(myc$plain, myc$invlog)))
cat("\nper-sample values:\n"); print(ac)

