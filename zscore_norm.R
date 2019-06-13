# z-score normalizartion
# for the esets

library(data.table)
library(dplyr)

# ---
# functions
# ---

# general zscore
zscore <- function(x) {
  (x - mean(x))/sd(x)
}

# zscore per sample (column)
zscorew <- function(eset, sampattern) {
  cns <- colnames(eset)
  ze <- apply(dplyr::select(eset, matches(sampattern)), 2, zscore)
  ccns <- cns[which(!seq_along(cns) %in% grep(sampattern, cns))]
  cbind(dplyr::select(eset, matches(paste(ccns, collapse = "|"))), ze)
}

filetype <- "raw_eset"

efiles <- list.files(path = "../data/rawdata/", pattern = filetype)
efiles

es <- lapply(efiles, function(...) fread(paste0("../data/rawdata/", ...)))
zse <- lapply(es, function(...) zscorew(..., "GSM"))

# write esets

mapply(function(x, y) { # x = raw file name, y = yugene normalized matrix
  outname <- gsub("_raw_", "_zscore_", x)
  write.table(y, paste0("../data/rawdata/", outname), sep = "\t", quote = FALSE, row.names = FALSE)
}, efiles, zse)
