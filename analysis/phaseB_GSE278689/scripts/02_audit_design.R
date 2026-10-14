# =============================================================================
# 02_audit_design.R -- Phase B step 2, GSE278689
#
# READ-ONLY. States, before any expression-based eligibility is evaluated, how
# each FROZEN rule is realised in THIS deposit, and flags anything the deposit
# cannot support. This is a design audit: it decides nothing by expression.
#
# Frozen rules realised here:
#   eligible tumour label : literal `tissue == "Tumor"` (the deposit supplies the
#                           literal string `Tumor`; the frozen rule admits labels
#                           semantically equivalent to primary PDAC tumour tissue
#                           provided the accompanying metadata identify the
#                           specimen as primary PDAC tumour tissue -- GEO
#                           `!Sample_source_name_ch1 = Tumor` and
#                           `!Sample_characteristics_ch1 = tissue: Tumor` do so)
#   cell type             : deposited `all_celltype == "Myeloid"` -- no re-annotation
#   cell count            : >= 50 included tumour-myeloid cells per sample
#   both groups           : >= 1 SIGLEC15-detectable AND >= 1 non-detectable
#                           (NO >=2/>=5 minimum -- that is the separate restricted
#                            TREM2_TAM sensitivity rule)
#   exclusions            : patient overlap, duplicate samples, integrated objects
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

ROOT <- "E:/Tex/HMlong/BJC_EXT_v1_20260927/phaseB_GSE278689"
DATA <- file.path(ROOT, "data")
P <- function(...) cat(..., "\n", sep = "")
P("time = ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))

md <- fread(file.path(DATA, "GSE278689_sn_metadata.csv.gz"), header = TRUE)
setnames(md, names(md)[1], "barcode")
md[, sample_id := sub("_[^_]+$", "", barcode)]

P("\n===== RULE 1: ELIGIBLE TUMOUR LABEL =====")
P("deposited `tissue` values          : ", paste(sort(unique(md$tissue)), collapse = ", "))
P("rule for this dataset              : tissue == \"Tumor\"")
P("samples satisfying it              : ", uniqueN(md[tissue == "Tumor", sample_id]),
  " of ", uniqueN(md$sample_id))
P("cells carrying the label           : ", md[tissue == "Tumor", .N], " of ", nrow(md))
P("GEO counterpart (provenance probe) : !Sample_source_name_ch1 = 'Tumor' for all 8 GSMs;")
P("                                     !Sample_characteristics_ch1 = 'tissue: Tumor'")
P("NOTE the frozen rule also admits PDAC / PDAC_tumor / PDAC_tissue. The deposit uses")
P("     none of those strings; it uses the literal 'Tumor' branch of the same rule.")
P("     Mapping fixed HERE, before any expression-based eligibility is evaluated.")

P("\n===== RULE 2: CELL TYPE (no re-annotation) =====")
P("deposited `all_celltype` levels    : ", uniqueN(md$all_celltype))
P("literal top-level `Myeloid` present: ", "Myeloid" %in% md$all_celltype)
P("cells labelled Myeloid             : ", md[all_celltype == "Myeloid", .N])
P("no marker-based re-call is performed anywhere in this run.")

P("\n===== RULE 3-4: SAMPLE IDENTIFIERS AND EXCLUSIONS =====")
smp <- unique(md[, .(sample_id, patient = patients, tissue)])
P("samples                            : ", nrow(smp))
P("distinct patients                  : ", uniqueN(smp$patient))
P("samples per patient max            : ", max(smp[, .N, by = patient]$N),
  "   (patient overlap would exclude)")
P("duplicated (patient, tissue) pairs : ",
  sum(duplicated(smp[, .(patient, tissue)])))
P("samples per (patient, tissue) max  : ",
  max(smp[, .N, by = .(patient, tissue)]$N))
P("integrated objects                 : none -- the deposit is one snRNA series;")
P("                                     no batch-corrected/integrated embedding is used,")
P("                                     and no expression matrix other than the deposited")
P("                                     count matrix is read.")

P("\n-- sample table as deposited --")
print(smp[order(patient)])

P("\n===== WHAT THIS DEPOSIT DOES *NOT* CONTAIN =====")
P("* no sample_id column -- sample identity is encoded in the barcode prefix")
P("  (8 distinct prefixes, each mapping 1:1 onto one patient; verified in 00_)")
P("* no raw FASTQ / ambient-RNA correction -- GEO supplement is a Cell Ranger")
P("  filtered count matrix (10x Genomics, Cell Ranger 7.1.0, GRCh38)")
P("* no nuclei-vs-cell flag -- this whole series is snRNA by construction")
P("* no second modality in this series -- comparability with the scRNA series is")
P("  therefore a between-dataset statement, never a within-dataset one.")

P("\n===== THRESHOLDS / CONVENTIONS CHOSEN BY THIS RUN =====")
P("none. Every value used downstream is either a frozen constant (>=50 cells;")
P("both groups present; CP10K+log1p; ddof=0; seed 15072026; 10000 resamples) or")
P("the deposit's own labels. Where the deposit forces a choice not covered by the")
P("brief it is reported in failed_or_excluded.md rather than resolved silently.")

P("\nDONE ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
