# UPC normalized affy and Illu raw files
setwd("/Users/cesar/Documents/jia_gxa_project/src/")

library(SCAN.UPC)
library(data.table)
library(dplyr)

# affy.CEL files
upcaffy <- UPCfast("../data/rawdata/*.CEL")
upcaffym <- exprs(upcaffy)
upcaffym <- cbind(data.frame('probe' = rownames(upcaffym)), upcaffym)
write.table(upcaffym, "../data/rawdata/affymetrix_UPC_eset.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# Illumina files
# GPL6106 Sentrix Human-6 v2 Expression BeadChip
fs <- list.files(path = "../data/rawdata/", pattern = ".txt")
fs <- fs[grep("GSM438", fs)] 
names(fs) <- gsub(".txt", "", fs)
fs

exps <- lapply(fs, function(...) fread(paste0("../data/rawdata/", ...)))
sapply(exps, dim)

exps <- lapply(names(exps), function(en) {
          e <- exps[[en]]
          cn <- colnames(e)[2:3]
          cn <- paste(en, cn, sep = "_")
          colnames(e)[2:3] <- cn
          e
      })

exptab <- Reduce(function(x, y) merge(x, y, by = "ID"), exps)

# Background correction
exptab[['nd']] <- apply(select(exptab, matches("Detection")), 1, function(...) length(which(... > 0.05))/84)
exptab <- filter(exptab, nd < 1)
eset <- select(exptab, matches("Signal|ID"))
id <- eset$ID
eset <- select(eset, matches("GSM"))
teset <- apply(eset, 2, function(c) {
          idx <- which(c < 3) 
          c[idx] <- rnorm(length(idx), 5, 1)
          c
        })
teset <- log2(teset)
colnames(teset) <- gsub("_Signal", "", colnames(teset))
meset <- data.matrix(teset)
rownames(meset) <- id
# write raw eset
raweset <- meset
raweset <- cbind(data.frame('probe' = rownames(raweset)), raweset)
write.table(raweset, "../data/rawdata/GPL6106_illumina_raw_eset.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
# normalize
eset <- ExpressionSet(meset)
upceset <- UPC_Generic_ExpressionSet(eset)
# writing UPC illu eset
upcesetm <- exprs(upceset)
upcesetm <- cbind(data.frame('probe' = rownames(upcesetm)), upcesetm)
write.table(upcesetm, "../data/rawdata/GPL6106_illumina_UPC_eset.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# GPL6884 
fs <- "../data/rawdata/GSE55319_non-normalized.txt"
exp <- fread(fs)

exp[['nd']] <- apply(dplyr::select(exp, matches("Pval")), 1, function(...) length(which(... > 0.05))/112)

exp <- filter(exp, nd < 1)
eset <- dplyr::select(exp, matches("Average|ID_REF"))
id <- eset$ID_REF
eset <- dplyr::select(eset, matches("Average"))
# min values summary
summary(apply(eset, 2, summary)[1,])
# --
# log2 scale 
eset <- log2(eset)
# colnames formatting
colnames(eset) <- gsub(" Average", "", colnames(eset))
## ref colnmaes to GSM accessions
meset <- data.matrix(eset)
gsmref <- fread("../data/GSE55319.txt", h = FALSE)
gsmref$V2 == colnames(meset)
colnames(meset) <- gsmref$V1
rownames(meset) <- id
# write raw eset
raweset <- meset
raweset <- cbind(data.frame('probe' = rownames(raweset)), raweset)
write.table(raweset, "../data/rawdata/GPL6884_illumina_raw_eset.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
# normalize
eset <- ExpressionSet(meset)
upceset <- UPC_Generic_ExpressionSet(eset)
# writing UPC illu eset
upcesetm <- exprs(upceset)
upcesetm <- cbind(data.frame('probe' = rownames(upcesetm)), upcesetm)
write.table(upcesetm, "../data/rawdata/GPL6884_illumina_UPC_eset.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

