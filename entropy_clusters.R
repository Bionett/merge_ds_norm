# Calculate the mean Shannon entropy for each clustering result
# based on the disease diversity per cluster
# if all clusters contain samples from a unique disease class, then H = 0,
# otherwise, H > 0. H is maximum when all disease classes are uniformely 
# distributed inside each cluster

library(data.table)
library(dplyr)
library(entropy)



library(philentropy)
#BiocManager::install("philentropy")



entrop <- function(file, inputfolder, metadat) {
  clustr <- fread(paste0(inputfolder, "/", file))
  clustr <- merge(clustr, dplyr::select(metadat, matches("disease|series_|Sample_geo_")), by = "Sample_geo_accession")
  print(head(clustr))
  clustr[["disease"]] <- as.numeric(as.factor(clustr[["Characteristics..disease."]]))
  cls <- split(clustr, clustr$`cluster`)
  cls <- lapply(cls, function(cl) cl[["disease"]])
  mentr <- mean(sapply(cls, entropy))
  ment <-  mean(sapply(cls, function(c) H(table(c)/sum(table(c)))))
  setNames(c(mentr, ment), c("entropy", "ShannonEntropy"))
}

inputfolder <- "../result/sampl_clust"
files <- list.files(path = inputfolder, pattern = ".tsv")
files <- files[grep("_\\d+cluster", files)]
files


metadat <- fread("../data/metadata_blood_sjia_pjia_normal_samples.tsv")


entab <- Reduce(rbind,
                lapply(files, function(file) entrop(file, inputfolder, metadat))
                )
entab <- data.frame(entab)
entab[["clust_ds"]] <- gsub(".tsv", "", files)

head(arrange(entab, ShannonEntropy))
head(arrange(entab, entropy))

write.table(entab, "../result/entropy_clusters.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
