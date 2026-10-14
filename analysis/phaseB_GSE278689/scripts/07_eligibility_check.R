# =============================================================================
# 07_eligibility_check.R -- Phase B step 2, GSE278689
#
# Builds output/sample_manifest.csv and the exclusion arithmetic, applying the
# FROZEN rules exactly and nothing else.
#
# The rules, in the order they are applied:
#   1. eligible tumour label : tissue == "Tumor"
#   2. cell type             : all_celltype == "Myeloid"      (deposited; no re-annotation)
#   3. >= 50 included tumour-myeloid cells per sample
#   4. BOTH groups present   : >= 1 SIGLEC15 RNA-detectable AND >= 1 non-detectable
#                              *** THERE IS NO >=2 / >=5 / >=10 MINIMUM HERE. ***
#                              The `>=10 cells with >=2 detectable and >=5 non-detectable`
#                              rule belongs to the SEPARATE restricted TREM2_TAM
#                              sensitivity analysis and is NOT applied in this file.
#   5. exclusions            : patient overlap, duplicate samples, integrated objects
#
# Detectability is read from the RAW deposited count (> 0). This is invariant to
# the CP10K + log1p reconstruction used downstream (verified in 05_).
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

ROOT <- "E:/Tex/HMlong/BJC_EXT_v1_20260927/phaseB_GSE278689"
OUT  <- file.path(ROOT, "output")
P <- function(...) cat(..., "\n", sep = "")
P("time = ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))

d    <- fread(file.path(OUT, "_myeloid_gene_cache.csv"))
allc <- fread(file.path(OUT, "_allcells_siglec15_cache.csv"))
P("myeloid cache : ", nrow(d), " x ", ncol(d))
P("all-cell cache: ", nrow(allc), " x ", ncol(allc))

## ===================== 1. SAMPLE MANIFEST ==================================
P("\n===== 1. SAMPLE MANIFEST =====")
agg <- d[, .(myeloid_cells = .N,
             sig_det = sum(SIGLEC15 > 0),
             sig_non = sum(SIGLEC15 == 0)), by = sample_id]
mn  <- unique(d[, .(sample_id, patient, tissue)])
mn  <- merge(mn, agg, by = "sample_id", all.x = TRUE)
tc  <- allc[, .(total_cells = .N), by = sample_id]
mn  <- merge(mn, tc, by = "sample_id", all.x = TRUE)
setnames(mn, c("sig_det", "sig_non"),
             c("siglec15_detectable_myeloid", "siglec15_nondetectable_myeloid"))

mn[, passes_min50_myeloid := myeloid_cells >= 50]
mn[, passes_both_groups   := passes_min50_myeloid &
                              siglec15_detectable_myeloid    >= 1 &
                              siglec15_nondetectable_myeloid >= 1]
mn[, eligible_tumour_label := tissue == "Tumor"]

## ---- exclusions that are properties of the DEPOSIT, not of expression ------
dupe_patient <- mn[, .N, by = patient][N > 1, patient]
dupe_sample  <- mn[duplicated(sample_id), sample_id]
P("patient overlap (same patient, >1 sample) : ", length(dupe_patient),
  if (length(dupe_patient)) paste0(" -> ", paste(dupe_patient, collapse = ",")) else "")
P("duplicate sample_id                       : ", length(dupe_sample))

mn[, exclusion_reason := ""]
mn[eligible_tumour_label == FALSE,
   exclusion_reason := sprintf("ineligible tissue label '%s' (frozen rule requires primary PDAC tumour tissue)", tissue)]
mn[eligible_tumour_label == TRUE & patient %in% dupe_patient,
   exclusion_reason := sprintf("patient overlap: patient %s contributes more than one sample", patient)]
mn[eligible_tumour_label == TRUE & sample_id %in% dupe_sample,
   exclusion_reason := "duplicate sample"]
mn[eligible_tumour_label == TRUE & exclusion_reason == "" & passes_min50_myeloid == FALSE,
   exclusion_reason := sprintf("only %d tumour-myeloid cells (<50)", myeloid_cells)]
mn[eligible_tumour_label == TRUE & exclusion_reason == "" & passes_min50_myeloid == TRUE &
   passes_both_groups == FALSE,
   exclusion_reason := sprintf(
     "SIGLEC15-detectable group empty (0 of %d tumour-myeloid cells carry a non-zero SIGLEC15 count); the primary contrast needs BOTH groups",
     myeloid_cells)]
mn[, final_eligible := eligible_tumour_label == TRUE & passes_min50_myeloid == TRUE &
                        passes_both_groups == TRUE & exclusion_reason == ""]

setcolorder(mn, c("sample_id","patient","tissue","total_cells","myeloid_cells",
                  "siglec15_detectable_myeloid","siglec15_nondetectable_myeloid",
                  "passes_min50_myeloid","passes_both_groups","eligible_tumour_label",
                  "final_eligible","exclusion_reason"))
setorder(mn, -final_eligible, patient)
fwrite(mn, file.path(OUT, "sample_manifest.csv"))
P("manifest rows = ", nrow(mn), "   final eligible = ", sum(mn$final_eligible))
print(mn[, .(sample_id, patient, tissue, total_cells, myeloid_cells,
             siglec15_detectable_myeloid, siglec15_nondetectable_myeloid,
             passes_min50_myeloid, passes_both_groups, final_eligible)], nrows = 40)

## ===================== 2. EXCLUSION ARITHMETIC =============================
P("\n===== 2. EXCLUSION ARITHMETIC (8 -> final estimable n) =====")
steps <- data.table(
  step = c("GSE278689 samples deposited",
           "- ineligible tissue label (tissue != 'Tumor')",
           "= tumour-labelled samples",
           "- patient overlap between samples",
           "- duplicate samples",
           "- integrated objects",
           "- fewer than 50 tumour-myeloid cells",
           "- SIGLEC15-detectable group empty (no non-zero count in any myeloid cell)",
           "- SIGLEC15-non-detectable group empty",
           "= FINAL ESTIMABLE n"),
  n = c(nrow(mn),
        -sum(!mn$eligible_tumour_label),
        sum(mn$eligible_tumour_label),
        -sum(mn$eligible_tumour_label & mn$patient %in% dupe_patient),
        -length(dupe_sample),
        0L,
        -sum(mn$eligible_tumour_label & mn$passes_min50_myeloid == FALSE),
        -sum(mn$eligible_tumour_label & mn$passes_min50_myeloid &
             mn$siglec15_detectable_myeloid == 0),
        -sum(mn$eligible_tumour_label & mn$passes_min50_myeloid &
             mn$siglec15_nondetectable_myeloid == 0),
        sum(mn$final_eligible)))
for (i in seq_len(nrow(steps)))
  P(sprintf("%-64s %s", steps$step[i],
            ifelse(steps$n[i] < 0, sprintf("-%d", -steps$n[i]), sprintf("%d", steps$n[i]))))
fwrite(steps, file.path(OUT, "exclusion_arithmetic.csv"))

P("\n-- the primary rule carries NO minimum group size --")
P("eligible samples and their detectable-cell counts:")
print(mn[final_eligible == TRUE, .(sample_id, patient,
      detectable = siglec15_detectable_myeloid,
      nondetectable = siglec15_nondetectable_myeloid)][order(detectable)])
P("samples with <5 detectable cells that ARE retained by the primary rule: ",
  mn[final_eligible == TRUE & siglec15_detectable_myeloid < 5, .N])
P("  -> these would be dropped by the SEPARATE restricted TREM2_TAM rule")
P("     (>=10 cells with >=2 detectable and >=5 non-detectable). That rule is NOT")
P("     applied here; applying it to the primary contrast would be a convention error.")

## ===================== 3. INVARIANT PROGRAMME GENES ========================
P("\n===== 3. INVARIANT PROGRAMME GENES WITHIN SAMPLE =====")
inv <- d[sample_id %in% mn[final_eligible == TRUE, sample_id],
         .(acp5_sd = sd(ACP5), mmp9_sd = sd(MMP9),
           acp5_nz = sum(ACP5 > 0), mmp9_nz = sum(MMP9 > 0),
           n = .N), by = sample_id]
inv[, acp5_invariant := acp5_sd == 0]
inv[, mmp9_invariant := mmp9_sd == 0]
print(inv)
P("samples with an invariant programme gene : ",
  sum(inv$acp5_invariant | inv$mmp9_invariant))
P("(An invariant gene would be omitted for that sample and, because the primary score")
P(" REQUIRES both ACP5 and MMP9, the sample would drop out of the primary contrast.)")
fwrite(inv, file.path(OUT, "programme_gene_variance_by_sample.csv"))

P("\nDONE ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
