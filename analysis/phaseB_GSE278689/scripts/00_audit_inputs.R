# =============================================================================
# 00_audit_inputs.R -- Phase B step 2, GSE278689
#
# READ-ONLY input audit. Establishes, before anything is computed:
#   * the four deposited files, their dimensions and their internal consistency
#   * the barcode -> sample -> patient mapping
#   * that `all_celltype` contains a literal top-level `Myeloid`
#   * that every metadata row maps to a barcode present in the barcode file
#
# Nothing here chooses a threshold. It only reports what the deposit contains.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

ROOT <- "E:/Tex/HMlong/BJC_EXT_v1_20260927/phaseB_GSE278689"
DATA <- file.path(ROOT, "data")
P <- function(...) cat(..., "\n", sep = "")
P("time = ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
P("R    = ", R.version.string)

## ------------------------------- files --------------------------------------
P("\n===== 1. DEPOSITED FILES =====")
fs <- file.info(list.files(DATA, full.names = TRUE))
print(data.frame(file = basename(rownames(fs)), bytes = fs$size,
                 mtime = format(fs$mtime, "%Y-%m-%d %H:%M:%S"), row.names = NULL))

## ------------------------------- barcodes -----------------------------------
bc <- fread(file.path(DATA, "GSE278689_sn_barcodes.tsv.gz"), header = FALSE)
setnames(bc, "barcode")
P("\n===== 2. BARCODES =====")
P("barcodes            : ", nrow(bc))
P("distinct            : ", uniqueN(bc$barcode))
P("duplicated          : ", sum(duplicated(bc$barcode)))
bc[, sample_id := sub("_[^_]+$", "", barcode)]
P("distinct prefixes   : ", uniqueN(bc$sample_id))
print(bc[, .(n_cells = .N), by = sample_id][order(sample_id)])

## ------------------------------- features -----------------------------------
ft <- fread(file.path(DATA, "GSE278689_sn_features.tsv.gz"), header = FALSE)
P("\n===== 3. FEATURES =====")
P("feature table dims  : ", nrow(ft), " x ", ncol(ft))
P("column classes      : ", paste(sapply(ft, class), collapse = ", "))
print(head(ft, 3))
genes <- ft[[2]]
P("gene symbols        : ", length(genes), " distinct ", uniqueN(genes))
for (g in c("SIGLEC15", "ACP5", "MMP9")) {
  P(sprintf("  %-9s present=%s  rows=%s", g, g %in% genes,
            paste(which(genes == g), collapse = ",")))
}
## duplicate symbols would break a name-based row lookup -- report, do not fix
dup <- genes[duplicated(genes)]
P("duplicated symbols  : ", uniqueN(dup), if (uniqueN(dup)) paste0(" -> ", paste(unique(dup), collapse = ",")) else "")

## ------------------------------- metadata -----------------------------------
md <- fread(file.path(DATA, "GSE278689_sn_metadata.csv.gz"), header = TRUE)
P("\n===== 4. METADATA =====")
P("metadata dims       : ", nrow(md), " x ", ncol(md))
P("column names        : ", paste(sprintf("'%s'", names(md)), collapse = ", "))
print(head(md, 3))
## the deposited header for the barcode column is empty -> fread names it V1
bc_col <- names(md)[1]
setnames(md, bc_col, "barcode")
md[, sample_id := sub("_[^_]+$", "", barcode)]

P("\n-- tissue --"); print(md[, .N, by = tissue][order(-N)])
P("\n-- patients --"); print(md[, .N, by = patients][order(patients)])
P("\n-- all_celltype --"); print(md[, .N, by = all_celltype][order(-N)])

P("\n-- sample_id x patient (must be 1:1) --")
xt <- dcast(md[, .N, by = .(sample_id, patients)], sample_id ~ patients,
            value.var = "N", fill = 0L)
print(xt)
P("prefixes per patient: ", paste(unique(md[, .(n = uniqueN(sample_id)), by = patients]$n), collapse = ","))
P("patients per prefix: ", paste(unique(md[, .(n = uniqueN(patients)), by = sample_id]$n), collapse = ","))

## --------------------- cross-file consistency -------------------------------
P("\n===== 5. CROSS-FILE CONSISTENCY =====")
P("barcodes == metadata rows            : ", nrow(bc) == nrow(md))
P("barcode multisets identical          : ", identical(sort(bc$barcode), sort(md$barcode)))
P("barcodes in metadata absent from file: ", sum(!md$barcode %in% bc$barcode))
P("barcodes in file absent from metadata: ", sum(!bc$barcode %in% md$barcode))
P("metadata row order == barcode order  : ", identical(bc$barcode, md$barcode))
P("sample_id agrees across the two files: ",
  identical(bc[order(barcode), sample_id], md[order(barcode), sample_id]))

## --------------------- what the frozen rules need ---------------------------
P("\n===== 6. FROZEN-RULE PRE-CHECK (counts only, no filtering) =====")
P("samples deposited                    : ", uniqueN(md$sample_id))
P("samples with tissue == 'Tumor'       : ", uniqueN(md[tissue == "Tumor", sample_id]))
P("patients                             : ", uniqueN(md$patients))
my <- md[all_celltype == "Myeloid"]
P("cells with all_celltype == 'Myeloid' : ", nrow(my))
P("\n-- tumour-myeloid cells per sample --")
per <- my[tissue == "Tumor", .(myeloid_cells = .N), by = .(sample_id, patients)]
print(per[order(myeloid_cells)])
P("samples with >=50 tumour-myeloid     : ", per[myeloid_cells >= 50, .N], " of ", nrow(per))
P("minimum                              : ", min(per$myeloid_cells))
P("maximum                              : ", max(per$myeloid_cells))

P("\nDONE ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
