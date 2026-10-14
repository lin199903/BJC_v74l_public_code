# =============================================================================
# 00_audit_inputs.R  --  Phase B step 1, GSE278688 (FROZEN RUN)
# READ-ONLY audit of the deposited files. No analysis. No modification of raw.
# Output: stdout log only.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

ROOT <- "E:/Tex/HMlong/BJC_EXT_v1_20260927/phaseB_GSE278688"
cat("R:", R.version.string, "\n")
cat("time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

## ---- 1. barcodes -----------------------------------------------------------
bar_p <- file.path(ROOT, "data", "GSE278688_sc_barcodes.tsv.gz")
bar <- fread(bar_p, header = FALSE, sep = "\t", quote = "", data.table = FALSE)[[1]]
cat("== barcodes ==\n")
cat("file            :", bar_p, "\n")
cat("n               :", length(bar), "\n")
cat("n unique        :", length(unique(bar)), "\n")
cat("duplicated      :", sum(duplicated(bar)), "\n")
cat("head            :", paste(head(bar, 3), collapse = " | "), "\n")
cat("tail            :", paste(tail(bar, 2), collapse = " | "), "\n")

## Derive sample prefix = everything before the trailing -1 (strip _out suffix text)
pref <- sub("-[ACGT]+-1$", "", bar)
cat("n distinct prefixes:", length(unique(pref)), "\n")
cat("prefix pattern ex  :", paste(head(unique(pref), 5), collapse = " | "), "\n")
cat("\n-- cells per prefix (all) --\n")
print(sort(table(pref), decreasing = TRUE))
cat("\n")

## ---- 2. features -----------------------------------------------------------
feat_p <- file.path(ROOT, "data", "GSE278688_sc_features.tsv.gz")
feat <- fread(feat_p, header = FALSE, sep = "\t", quote = "", data.table = FALSE)
cat("== features ==\n")
cat("nrow x ncol     :", nrow(feat), "x", ncol(feat), "\n")
cat("first col head  :", paste(head(feat[[1]], 3), collapse = " | "), "\n")
if (ncol(feat) >= 2) cat("second col head :", paste(head(feat[[2]], 3), collapse = " | "), "\n")
cat("colnames        :", paste(colnames(feat), collapse = " | "), "\n")
for (g in c("SIGLEC15", "ACP5", "MMP9")) {
  hit <- which(feat[[1]] == g)
  cat(sprintf("  %-9s rows in col1: %s\n", g, if (length(hit)) paste(hit, collapse = ",") else "ABSENT"))
  if (ncol(feat) >= 2) {
    hit2 <- which(feat[[2]] == g)
    cat(sprintf("  %-9s rows in col2: %s\n", g, if (length(hit2)) paste(hit2, collapse = ",") else "ABSENT"))
  }
}
cat("\n-- any duplicated gene symbols in col2? --\n")
if (ncol(feat) >= 2) {
  tb <- table(feat[[2]]); cat("n dup symbols:", sum(tb > 1), "\n")
  cat("SIGLEC15 count in col2:", sum(feat[[2]] == "SIGLEC15"), "\n")
  cat("ACP5 count in col2    :", sum(feat[[2]] == "ACP5"), "\n")
  cat("MMP9 count in col2    :", sum(feat[[2]] == "MMP9"), "\n")
}
cat("\n")

## ---- 3. metadata -----------------------------------------------------------
meta_p <- file.path(ROOT, "data", "GSE278688_sc_metadata.csv.gz")
if (!file.exists(meta_p)) meta_p <- "E:/Tex/HMlong/BJC_EXT_v1_20260927/annotations/GSE278688_sc_metadata.csv.gz"
cat("== metadata ==\n")
cat("file            :", meta_p, "\n")
meta <- fread(meta_p, header = TRUE, sep = ",", quote = "\"", data.table = FALSE)
cat("dim             :", nrow(meta), "x", ncol(meta), "\n")
cat("colnames        :", paste(sprintf("[%s]", colnames(meta)), collapse = " "), "\n")
cat("head rownames   :", paste(head(rownames(meta), 3), collapse = " | "), "\n")
print(utils::head(meta, 3))
cat("\n-- tissue --\n"); print(table(meta$tissue, useNA = "ifany"))
cat("\n-- patients --\n"); print(table(meta$patients, useNA = "ifany"))
cat("\n-- all_celltype --\n"); print(sort(table(meta$all_celltype, useNA = "ifany"), decreasing = TRUE))
cat("\n")
cat("n cells with all_celltype == 'Myeloid':", sum(meta$all_celltype == "Myeloid", na.rm = TRUE), "\n")

## ---- 4. barcode <-> metadata linkage --------------------------------------
cat("\n== linkage barcodes vs metadata ==\n")
cat("meta rows not in barcodes :", sum(!rownames(meta) %in% bar), "\n")
cat("barcodes not in meta rows :", sum(!bar %in% rownames(meta)), "\n")
cat("\n== sessionInfo ==\n"); print(sessionInfo())
