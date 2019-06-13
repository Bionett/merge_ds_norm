# YuGene normalization

library(YuGene)


efiles <- list.files(path = "../data/rawdata/", pattern = "raw_eset")
efiles

es <- lapply(efiles, function(f) fread(paste0("../data/rawdata/", f)))

# Yugene
do_yugen <- function(eset) {
  e <- dplyr::select(eset, matches("GSM"))
  y <- YuGene(e)
  y <- apply(y, 2, function(x) x)
  cbind(data.frame("probe" = eset$probe), y)
}

yge <- lapply(es, do_yugen)

# write esets

mapply(function(x, y) { # x = raw file name, y = yugene normalized matrix
  outname <- gsub("_raw_", "_yugene_", x)
  write.table(y, paste0("../data/rawdata/", outname), sep = "\t", quote = FALSE, row.names = FALSE)
}, efiles, yge)

