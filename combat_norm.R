# ComBat batch correction

library(sva)


run_combat <- function(filetype, sampattern, genecol, inputfolder, metafile, batchcol) {
  esetfile <- paste0("blood_sjia_pjia_", filetype, ".tsv")
  eset <- fread(paste(inputfolder, esetfile, sep = "/"))
  meta <- read.delim(paste0(inputfolder, metafile))
  rownames(meta) <- meta[["Sample_geo_accession"]]
  batch <- setNames(meta[[batchcol]], meta[["Sample_geo_accession"]])
  
  esetmat <- dplyr::select(eset, matches(sampattern))
  esetmat <- data.matrix(esetmat)
  rownames(esetmat) <- eset[[genecol]] 
  
  combat_eset <- ComBat(dat = esetmat, batch = as.character(batch[colnames(esetmat)]), 
                        mod = NULL, par.prior = TRUE, 
                        prior.plots = FALSE)
  
  cbeset <- cbind(data.frame("Symbol" = rownames(combat_eset)), combat_eset)
  outname <- gsub(filetype, paste0(filetype, "_naive_combat"), esetfile)
  write.table(cbeset, paste0(inputfolder, outname), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # 
  meta <- meta[colnames(esetmat), ]
  mod <- model.matrix(~as.factor(Characteristics..disease.), data=meta)
  combat_eset <- ComBat(dat = esetmat, batch = as.character(batch[colnames(esetmat)]), 
                        mod = mod, par.prior = TRUE, 
                        prior.plots = FALSE)
  
  cbeset <- cbind(data.frame("Symbol" = rownames(combat_eset)), combat_eset)
  outname <- gsub(filetype, paste0(filetype, "_superv_combat"), esetfile)
  write.table(cbeset, paste0(inputfolder, outname), sep = "\t", row.names = FALSE, quote = FALSE)
}  

filetype <- "quantile_eset"
inputfolder <- "../data/"
sampattern <- "GSM"
genecol <- "Symbol"
metafile <- "metadata_blood_sjia_pjia_normal_samples.tsv"
batchcol <- "batch"

run_combat(filetype, sampattern, genecol, inputfolder, metafile, batchcol)
