# visualize dataset clustering
# using PCA, tSNE, and UMAP


#BiocManager::install("umap")
#BiocManager::install("Rtsne")

library(data.table)
library(dplyr)
library(ggfortify)
library(umap)
library(Rtsne)


inputfolder <- "../data"
outfolder <- "../result"
esetfiles <- list.files("../data", pattern = "_eset")
esetfiles <- esetfiles[grep(".tsv", esetfiles)]

metadat <- read.delim("../data/metadata_blood_sjia_pjia_normal_samples.tsv")
metadat[["manufacturer"]] <- ifelse(metadat$Sample_platform_id == "GPL570", "Affymetrix", "Illumina")

# variables of interest
info.var <- colnames(metadat)[c(5, 8, 10, 11, 17, 18, 20, 21, 38, 40, 58, 59)]
info.var <- info.var[c(3, 12)]


draw_pca <- function(pca_dat, esetfile, info.var, outfolder) {
  #pdf(paste0(outfolder, "/dimred_pca_", gsub(".tsv", ".pdf", esetfile)))
  #for(v in info.var) {
  #  plot(
  #    autoplot(prcomp(pca_dat[["x"]]), data = pca_dat[["metax"]], colour = v)
  #  )
  #}
  #dev.off()
  lapply(info.var, function(v) {
    autoplot(prcomp(pca_dat[["x"]]), data = pca_dat[["metax"]], colour = v) + theme_classic() +theme(legend.position="bottom", legend.title = element_blank())
  })
}

format_pca <- function(eset, metadat) {
  # setting eset format to the PCA visualization
  em <- dplyr::select(eset, matches("GSM"))
  rownames(em) <- eset$Symbol
  em <- t(em)
  colnames(em) <- gsub("\\.|\\-", "_", colnames(em))
  dat <- cbind(data.frame("Sample_geo_accession" = rownames(em)), em)
  print(dat[1:5, 1:10])
  dat <- data.frame(merge(metadat, dat, by = "Sample_geo_accession"))
  rownames(dat) <- as.character(dat[, 'Sample_geo_accession'])
  print(dat[1:5, 1:10])
  list("x" = em, "metax" = dat)
}


draw_tsne <- function(pca_dat, metadat) {
  tsne <- Rtsne(pca_dat[["x"]], dims = 2, perplexity = 30, verbose=TRUE, max_iter = 500)
  tsnedf <- as.data.frame(tsne$Y)
  colnames(tsnedf) <- c("tSNE1", "tSNE2")
  tsnedf[["Sample_geo_accession"]] <- rownames(pca_dat[["x"]])
  tsnedf <- merge(tsnedf, metadat, by = "Sample_geo_accession")
  
  g1 <- ggplot(tsnedf, aes(x = tSNE1, y = tSNE2)) +
    geom_point(aes(color = Characteristics..disease.)) +
    theme_classic() +
    theme(legend.position="bottom", legend.title = element_blank())
  
  g2 <- ggplot(tsnedf, aes(x = tSNE1, y = tSNE2)) +
    geom_point(aes(color = batch)) +
    theme_classic() + 
    theme(legend.position="bottom", legend.title = element_blank())
  
  list("tsnedisease" = g1, "tsnebatch" = g2)
}

draw_umap <- function(pca_dat, metadat) {
  umapres <- umap(pca_dat[["x"]], dims = 2, perplexity = 30, verbose=TRUE, max_iter = 500)
  tsnedf <- as.data.frame(umapres$layout)
  colnames(tsnedf) <- c("UMAP1", "UMAP2")
  tsnedf[["Sample_geo_accession"]] <- rownames(pca_dat[["x"]])
  tsnedf <- merge(tsnedf, metadat, by = "Sample_geo_accession")
  
  g1 <- ggplot(tsnedf, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = Characteristics..disease.)) +
    theme_classic() +
    theme(legend.position="bottom", legend.title = element_blank())
  
  g2 <- ggplot(tsnedf, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = batch)) +
    theme_classic() + 
    theme(legend.position="bottom", legend.title = element_blank())
  
  list("umapdisease" = g1, "umapbatch" = g2)
}

draw_dimred <- function(esetfile, inputfolder, outfolder, metadat, info.var) {
  eset <- fread(paste0(inputfolder, "/", esetfile))
  # PCA
  pca_dat <- format_pca(eset, metadat)
  pcas <- draw_pca(pca_dat, esetfile, info.var, outfolder)
  # tSNE
  tsnes <- draw_tsne(pca_dat, metadat)
  # UMAP
  umaps <- draw_umap(pca_dat, metadat)
  # writing plots
  dimreds <- append(append(pcas, tsnes), umaps)
  pdf(paste0(outfolder, "/dimred_pca_tsne_umap_", gsub(".tsv", ".pdf", esetfile)), width = 10)
    multiplot(dimreds[[1]], dimreds[[2]], dimreds[[3]], dimreds[[4]], dimreds[[5]], dimreds[[6]], cols = 3)
  dev.off()
}

mclapply(esetfiles, function(esetfile) {
  draw_dimred(esetfile, inputfolder, outfolder, metadat, info.var)
}, mc.cores = 3)

