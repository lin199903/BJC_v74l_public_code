# 07_eligibility_check.R -- Phase B step 1, GSE278688 (FROZEN RUN)
# Applies the frozen eligibility rules to the TUMOUR samples only and prints the
# 29 -> ... -> n arithmetic. Uses the myeloid cache (raw counts); detectability
# is invariant to normalisation (verified: every candidate normalisation is
# monotone with f(0)=0), so this eligibility set is normalisation-independent.
suppressPackageStartupMessages({ library(data.table) })
ROOT <- "E:/Tex/HMlong/BJC_EXT_v1_20260927/phaseB_GSE278688"
P <- function(...) cat(..., "\n", sep="")
d <- fread(file.path(ROOT,"output","_myeloid_gene_cache.csv"))
meta <- fread(file.path(ROOT,"data","GSE278688_sc_metadata.csv.gz"), header=TRUE, sep=",", quote="\"")
meta$sample_id <- sub("[_-]GEX_out_.*$", "", meta[[1]])

tot <- as.data.table(table(meta$sample_id)); setnames(tot, c("sample_id","total_cells"))
## NOTE: compute all per-sample summaries FIRST, then drop the per-cell rows,
## to avoid data.table `by=` scoping problems.
d[, sig_pos := as.integer(SIGLEC15 > 0)]
agg <- d[, .(myeloid_cells = .N, sig_det = sum(sig_pos),
             sig_nondet = sum(1L - sig_pos)), by = sample_id]
s <- unique(d[, .(sample_id, tissue, patient)])
s <- merge(s, agg, by="sample_id", all.x=TRUE)
s <- merge(s, tot, by="sample_id", all.x=TRUE)
setorder(s, patient, tissue)

P("=== ALL 29 SAMPLES: myeloid + SIGLEC15 detectability ===")
print(s[, .(sample_id, tissue, myeloid_cells, sig_det, sig_nondet)], nrows=30)

P("\n=== FROZEN RULE ARITHMETIC ===")
P("deposited samples                                   : ", nrow(s))
P("  minus tissue != 'Tumor'  (Adjacent_normal + PBMC) : ", sum(s$tissue != "Tumor"))
P("  = tumour samples (eligible label, pre-cell-count) : ", sum(s$tissue == "Tumor"))
tm <- s[tissue == "Tumor"]
P("  minus tumour samples with <50 myeloid cells       : ", sum(tm$myeloid_cells < 50))
tm2 <- tm[myeloid_cells >= 50]
P("  = tumour samples with >=50 myeloid                : ", nrow(tm2))
P("  minus samples with NO SIGLEC15-detectable myeloid : ", sum(tm2$sig_det == 0))
tm3 <- tm2[sig_det >= 1]
P("  minus samples with NO SIGLEC15-non-detectable     : ", sum(tm3$sig_nondet == 0))
P("  = FINAL ESTIMABLE n                               : ", nrow(tm3))

P("\n=== ELIGIBLE SAMPLES (final) ===")
print(tm3[, .(sample_id, patient, tissue, myeloid_cells, siglec15_detectable_myeloid=sig_det,
              siglec15_nondetectable_myeloid=sig_nondet)], nrows=20)

P("\n=== TUMOUR SAMPLES EXCLUDED FOR MISSING DETECTABLE GROUP ===")
ex <- tm2[sig_det == 0]
print(ex[, .(sample_id, patient, myeloid_cells, sig_det, sig_nondet)], nrows=20)

P("\n=== PATIENT OVERLAP AMONG ELIGIBLE SAMPLES ===")
pt <- tm3[, .(n_eligible_samples = .N), by=patient]
print(pt, nrows=20)
P("patients contributing >1 eligible tumour sample: ", sum(pt$n_eligible_samples > 1))
P("duplicate samples (same sample_id twice)       : ", sum(duplicated(s$sample_id)))
P("patients represented in the eligible set       : ", nrow(pt), " of ", length(unique(s$patient)))

P("\n=== SIGLEC15 SPARSITY, WHOLE DATASET ===")
allc <- fread(file.path(ROOT,"output","_allcells_siglec15_cache.csv"))
P("all cells                      : ", nrow(allc))
P("all cells SIGLEC15 > 0         : ", sum(allc$SIGLEC15 > 0))
P("  in Myeloid                   : ", sum(allc$SIGLEC15 > 0 & allc$all_celltype=="Myeloid"))
P("  NOT in Myeloid               : ", sum(allc$SIGLEC15 > 0 & allc$all_celltype!="Myeloid"))
P("myeloid cells                  : ", sum(allc$all_celltype=="Myeloid"))
P("myeloid SIGLEC15+ rate         : ", sprintf("%.3f%%", 100*sum(allc$SIGLEC15>0 & allc$all_celltype=="Myeloid")/sum(allc$all_celltype=="Myeloid")))
P("mean SIGLEC15 over all cells   : ", sprintf("%.6f", mean(allc$SIGLEC15)))
