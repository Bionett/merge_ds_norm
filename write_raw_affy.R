# write affymetrix raw eset

library(affy)

# files
affyfiles <- list.files(path = "../data/rawdata/", pattern = "GSM")
affyfiles <- affyfiles[-grep(".txt", affyfiles)]
affyfiles

# reading cel files into a eset object
eset <- ReadAffy(filenames = paste0("../data/rawdata/", affyfiles))

# extracting scan date data from eset object
scandat <- pData(protocolData(eset))
sd <- Reduce(rbind, lapply(scandat$ScanDate, function(d) unlist(strsplit(d, "T| "))))

# original long date format
date <- sapply(sd[, 1], function(d) {
  if(grepl("/", d)) {
    dd <- unlist(strsplit(d, "/"))
    paste(paste0("20", dd[3]), dd[1], dd[2], sep = "-")
  } else {
    d
  }
}) 

sd <- data.frame(sd)
sd[["Sample_geo_accession"]] <- unlist(regmatches(rownames(scandat), regexec('GSM\\d+', rownames(scandat))))
sd[["scandate"]] <- date
# shorter date format (YYYY-MM)
date_yymm <- sapply(as.character(sd[["scandate"]]), function(d) substring(d, 1, 7))
sd[["scandate_yymm"]] <- date_yymm
# writing complementary metadata file
write.table(select(sd, matches("scan|geo")), "../data/metadata_scandate_affy_samples_sjia_pjia_normal.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
# -------------------------------------------------------

eset <- rma(eset, subset = NULL, verbose = TRUE, destructive = TRUE, 
            normalize = FALSE, background = TRUE)

emat <- exprs(eset)
cns <- colnames(emat)
colnames(emat) <- unlist(regmatches(cns, regexec("GSM\\d+", cns)))
emat <- cbind(data.frame("probe" = rownames(emat)), emat)
# writing raw eset
write.table(emat, "../data/rawdata/affymetrix_blood_sjia_pjia_raw_eset.tsv", sep = "\t", 
            quote = FALSE, row.names = FALSE)

