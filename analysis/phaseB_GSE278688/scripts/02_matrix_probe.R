# =============================================================================
# 02_matrix_probe.R -- READ-ONLY probe of the deposited matrix.
# Purpose: determine (a) dimensions, (b) whether stored values are raw integer
# counts or already log1p-normalised. This determines the ONLY branch point in
# the frozen pipeline. Streams the .mtx.gz without loading it.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
ROOT <- "E:/Tex/HMlong/BJC_EXT_v1_20260927/phaseB_GSE278688"
P <- function(...) cat(..., "\n", sep="")

mtx_p <- file.path(ROOT,"data","GSE278688_sc_matrix.mtx.gz")
P("file      : ", mtx_p)
P("bytes     : ", file.info(mtx_p)$size)

con <- gzfile(mtx_p, "rt")
on.exit(close(con))

## ---- header + comment lines ----
hdr <- NULL; ln <- 0L
repeat {
  l <- readLines(con, n=1L, warn=FALSE)
  if (length(l)==0) break
  ln <- ln + 1L
  if (substr(l,1,1) == "%") { P("comment   : ", l); next }
  hdr <- l; break
}
P("dims line : ", hdr)
d <- as.integer(strsplit(trimws(hdr), "\\s+")[[1]])
P("nrow(genes)=", d[1], " ncol(cells)=", d[2], " nnz=", d[3])

## ---- read first N data lines to judge integer-ness ----
N <- 200000L
chunk <- readLines(con, n=N, warn=FALSE)
P("scanned data lines: ", length(chunk))
vals <- as.numeric(sub("^\\S+\\s+\\S+\\s+", "", chunk))
P("value summary (first ", length(vals), "):")
print(summary(vals))
P("all integer-valued : ", all(vals == round(vals)))
P("min / max          : ", min(vals), " / ", max(vals))
P("n distinct values  : ", length(unique(vals)))
P("first 5 distinct   : ", paste(head(sort(unique(vals)),5), collapse=", "))
P("frac == 1          : ", round(mean(vals==1),4))
P("frac == 2          : ", round(mean(vals==2),4))
P("max value seen     : ", max(vals))
close(con)

## ---- full-file streaming check: integer-ness across whole file ----
con <- gzfile(mtx_p, "rt"); on.exit(close(con))
l <- readLines(con, n=1L, warn=FALSE)              # header
while (substr(l,1,1) == "%") l <- readLines(con, n=1L, warn=FALSE)
P("\n== full streaming scan ==")
n <- 0L; nonint <- 0L; mx <- 0; mn <- Inf
repeat {
  ch <- readLines(con, n=500000L, warn=FALSE)
  if (length(ch)==0) break
  v <- as.numeric(sub("^\\S+\\s+\\S+\\s+", "", ch))
  n <- n + length(v)
  nonint <- nonint + sum(v != round(v))
  mx <- max(mx, max(v)); mn <- min(mn, min(v))
}
P("total data lines read : ", n)
P("expected nnz          : ", d[3])
P("match                 : ", n == d[3])
P("non-integer values    : ", nonint)
P("global min / max      : ", mn, " / ", mx)
P("\nVERDICT: stored values are ",
  if (nonint > 0) "NON-INTEGER -> already normalised" else
    "ALL INTEGER -> raw counts", sep="")
