# =============================================================================
# 01b_audit_design_fixed.R -- Phase B step 1, GSE278688 (FROZEN RUN)
# Corrects the sample-ID parser (two naming conventions coexist) and prints
# compactly. READ-ONLY.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
ROOT <- "E:/Tex/HMlong/BJC_EXT_v1_20260927/phaseB_GSE278688"
P <- function(...) cat(..., "\n", sep="")

bar  <- fread(file.path(ROOT,"data","GSE278688_sc_barcodes.tsv.gz"),
              header=FALSE, sep="\t", quote="", data.table=FALSE)[[1]]
meta <- fread(file.path(ROOT,"data","GSE278688_sc_metadata.csv.gz"),
              header=TRUE, sep=",", quote="\"", data.table=FALSE)
mb <- meta[[1]]

## ---- canonical sample_id parser --------------------------------------------
## barcodes come in two conventions:
##   PDAC01-Tumor-GEX_out_<BC>-1     PDAC05_PBMC_GEX_out_<BC>-1
sid_parse <- function(x) sub("[_-]GEX_out_.*$", "", x)
## tissue = token after the last - or _
tissue_parse <- function(x) sub(".*[-_]", "", x)
patient_parse <- function(x) sub("[-_].*$", "", x)

P("===== A. BARCODE FILE vs METADATA (reconciliation) =====")
P("barcode file rows         : ", length(bar))
P("barcode rows unique       : ", length(unique(bar)))
P("barcode rows duplicated   : ", sum(duplicated(bar)))
P("metadata rows             : ", nrow(meta))
P("meta barcode col name     : [", colnames(meta)[1], "]")
P("meta barcodes unique      : ", length(unique(mb)))
P("meta barcodes duplicated  : ", sum(duplicated(mb)))
P("setequal(barcodes, metabc): ", setequal(bar, mb))
P("in meta not in barcode    : ", length(setdiff(mb, bar)))
P("in barcode not in meta    : ", length(setdiff(bar, mb)))

## multiset comparison (order-insensitive, keeps duplicates)
msd <- function(a,b) sum(table(a) != table(b))
P("\nmultiset equal            : ", identical(sort(bar), sort(mb)))
dups <- unique(mb[duplicated(mb)])
P("distinct dup barcode strs : ", length(dups))
if (length(dups)) { P("examples:"); print(head(sort(dups), 5)) }

P("\n===== B. SAMPLE STRUCTURE =====")
meta$sample_id <- sid_parse(mb)
meta$tissue_p  <- tissue_parse(meta$sample_id)
meta$patient_p <- patient_parse(meta$sample_id)
smap <- unique(meta[, c("sample_id","tissue_p","patient_p","tissue","patients")])
smap <- smap[order(smap$sample_id), ]
P("n sample_id               : ", nrow(smap))
P("n patients (barcode)      : ", length(unique(smap$patient_p)))
P("n patients (metadata col) : ", length(unique(meta$patients)))
P("n tissues (from barcode)  : ", paste(sort(unique(smap$tissue_p)), collapse=" | "))
P("n tissues (metadata col)  : ", paste(sort(unique(smap$tissue)), collapse=" | "))
P("barcode-tissue == meta-tissue for all rows : ",
  all(tolower(smap$tissue_p) == tolower(gsub("Adjacent_normal","Normal",smap$tissue))))
P("barcode-patient == meta-patient for all rows: ",
  all(paste0("PA", sub("^PDAC","",smap$patient_p)) == smap$patients))

P("\n-- the 29 samples --")
smap$n_cells <- as.integer(table(meta$sample_id)[smap$sample_id])
smap$n_myeloid <- as.integer(table(meta$sample_id[meta$all_celltype=="Myeloid"])[smap$sample_id])
smap$n_myeloid[is.na(smap$n_myeloid)] <- 0L
print(smap[, c("sample_id","patient_p","patients","tissue","n_cells","n_myeloid")], row.names=FALSE)

P("\n===== C. TISSUE x SAMPLE COUNTS =====")
P("samples by tissue (meta)  :")
print(table(smap$tissue))
P("cells by tissue (meta)    :")
print(table(meta$tissue))
P("\nsamples with >=50 myeloid, by tissue:")
print(table(smap$tissue[smap$n_myeloid >= 50]))
P("TOTAL samples with >=50 myeloid: ", sum(smap$n_myeloid >= 50))

P("\n===== D. PATIENT OVERLAP =====")
pt <- aggregate(cbind(n_samples=rep(1L,nrow(smap)),
                      n_tumor=as.integer(smap$tissue=="Tumor"),
                      n_tumor_my50=as.integer(smap$tissue=="Tumor" & smap$n_myeloid>=50)) ~ patient_p,
                data=smap, FUN=sum)
print(pt, row.names=FALSE)
P("\npatients with >1 sample in dataset      : ", sum(pt$n_samples > 1))
P("patients with >1 ELIGIBLE TUMOUR sample : ", sum(pt$n_tumor_my50 > 1))
P("patients contributing >=1 eligible tumour: ", sum(pt$n_tumor_my50 >= 1))

P("\n===== E. GENE PRESENCE =====")
feat <- fread(file.path(ROOT,"data","GSE278688_sc_features.tsv.gz"),
              header=FALSE, sep="\t", quote="", data.table=FALSE)
P("features nrow x ncol      : ", nrow(feat), " x ", ncol(feat))
for (g in c("SIGLEC15","ACP5","MMP9")) P(sprintf("  %-9s n_rows=%d", g, sum(feat[[2]]==g)))
P("duplicated symbols        : ", sum(table(feat[[2]])>1))
