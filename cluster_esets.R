# clustering eset datasets
# using HC, K-NN and density based methods
# for scRNA-seq

library(data.table)
library(dplyr)
library(parallel)
library(cluster)
library(factoextra)



# ----
# - Functions
# ----

draw_sample_dist <- function(esemat, esetfile, outfolder) {
  # distance calculation (3 availble methods)
  distmat <- get_dist(esemat, stand = FALSE, method = "spearman")
  pdf(paste0(outfolder, "/sample_dist_", gsub(".tsv", ".pdf", esetfile)), width = 15, height = 15)
  plot(
    fviz_dist(distmat, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
  )
  distmat <- get_dist(esemat, stand = FALSE, method = "euclidean")
  plot(
    fviz_dist(distmat, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
  )
  dev.off()
}

# -----
# --
# -----

do_kmeans <- function(esemat, esetfile, outfolder, metadat, samplecol, nclusts) {
  # determinig the best number of cluster by consensus of the silhouette and gap_stat methods
  #sil <- fviz_nbclust(esemat, kmeans, method = "silhouette")
  #gap <- fviz_nbclust(esemat, kmeans, method = "gap_stat")
  #best <- unlist(mapply(mean, sil$data[, 2], gap$data[, 2]))
  
  # visualizing best nclusters
  #pdf(paste0(outfolder, "/best_kmeans_nclust_", gsub(".tsv", ".pdf", esetfile)), width = 6, height = 6)
  #  plot(sil)
  #  plot(gap)
  #dev.off()
  #print(best)
  
  # "best" number of cluster for Kmeans
  kmres <- kmeans(esemat, nclusts[1], nstart = 25)
  # writing clustering result
  res <- data.frame("Sample_geo_accession" = names(kmres$cluster), "cluster" = kmres$cluster)
  write.table(res, paste0(outfolder, "/kmeans_", nclusts[1], "clusters_", esetfile), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Visualizing kmens clusters and disease and technical metadata
  pdf(paste0(outfolder, "/kmeans_", gsub(".tsv", ".pdf", esetfile)), width = 10, height = 10)
  plot(
    fviz_cluster(kmres, data = esemat, ellipse.type = "convex") +
      theme_minimal()
  )
  # best nclust+1 clusters
  kmres <- kmeans(esemat, nclusts[2], nstart = 25)
  # writing clustering result
  res <- data.frame("Sample_geo_accession" = names(kmres$cluster), "cluster" = kmres$cluster)
  write.table(res, paste0(outfolder, "/kmeans_", nclusts[2], "clusters_", esetfile), sep = "\t", row.names = FALSE, quote = FALSE)
  # Visualizing metadata
  plot(
    fviz_cluster(kmres, data = esemat, ellipse.type = "convex") +
      theme_minimal()
  )
  kmres$cluster <- setNames(metadat[["Characteristics..disease."]], metadat[[samplecol]])
  plot(
    fviz_cluster(kmres, data = esemat, ellipse.type = "convex") +
      theme_minimal()
  )
  kmres$cluster <- setNames(metadat[["Series_geo_accession"]], metadat[[samplecol]])
  plot(
    fviz_cluster(kmres, data = esemat, ellipse.type = "convex") +
      theme_minimal()
  )
  dev.off()
}

# -----
# --
# -----

do_clara <- function(esemat, esetfile, outfolder, metadat, samplecol, nclusts) {
  # determinig the best number of cluster by consensus of the silhouette and gap_stat methods
  #sil <- fviz_nbclust(esemat, clara, method = "silhouette")
  #gap <- fviz_nbclust(esemat, clara, method = "gap_stat")
  #best <- unlist(mapply(mean, sil$data[, 2], gap$data[, 2]))
  
  # visualizing best nclusters
  #pdf(paste0(outfolder, "/best_clara_nclust_", gsub(".tsv", ".pdf", esetfile)), width = 6, height = 6)
  #  plot(sil)
  #  plot(gap)
  #dev.off()
  
  #print(best)
  
  # "best" number of cluster for Kmeans
  clarares <- clara(esemat, nclusts[1])
  res <- data.frame("Sample_geo_accession" = names(clarares$cluster), "cluster" = clarares$cluster)
  write.table(res, paste0(outfolder, "/clara_", nclusts[1], "clusters_", esetfile), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Visualizing kmens clusters and disease and technical metadata
  pdf(paste0(outfolder, "/clara_", gsub(".tsv", ".pdf", esetfile)), width = 10, height = 10)
  plot(
    fviz_cluster(clarares, data = esemat, ellipse.type = "convex") +
      theme_minimal()
  )
  # best nclust+1 clusters
  clarares <- clara(esemat, nclusts[2])
  # writing clustering result
  res <- data.frame("Sample_geo_accession" = names(clarares$cluster), "cluster" = clarares$cluster)
  write.table(res, paste0(outfolder, "/clara_", nclusts[2], "clusters_", esetfile), sep = "\t", row.names = FALSE, quote = FALSE)
  # Visualizing metadata
  plot(
    fviz_cluster(clarares, data = esemat, ellipse.type = "convex") +
      theme_minimal()
  )
  clarares$cluster <- setNames(metadat[["Characteristics..disease."]], metadat[[samplecol]])
  plot(
    fviz_cluster(clarares, data = esemat, ellipse.type = "convex") +
      theme_minimal()
  )
  clarares$cluster <- setNames(metadat[["Series_geo_accession"]], metadat[[samplecol]])
  plot(
    fviz_cluster(clarares, data = esemat, ellipse.type = "convex") +
      theme_minimal()
  )
  dev.off()
}

# -----
# --
# -----

do_hc <- function(esemat, esetfile, outfolder, nclusts) {
  res <- hcut(esemat, k = as.numeric(nclusts[1]), stand = FALSE)
  # writing clustering result
  restab <- data.frame("Sample_geo_accession" = names(res$cluster), "cluster" = res$cluster)
  write.table(restab, paste0(outfolder, "/hc_", nclusts[1], "clusters_", esetfile), sep = "\t", row.names = FALSE, quote = FALSE)
  # visualazing tree
  pdf(paste0(outfolder, "/hc_", gsub(".tsv", ".pdf", esetfile)), width = 10, height = 6)
  plot(
    fviz_dend(res, rect = TRUE, cex = 0.5)
  )
  plot(
    fviz_silhouette(res)
  )
  res <- hcut(esemat, k = as.numeric(nclusts[2]), stand = FALSE)
  # writing clustering result
  restab <- data.frame("Sample_geo_accession" = names(res$cluster), "cluster" = res$cluster)
  write.table(restab, paste0(outfolder, "/hc_", nclusts[2], "clusters_", esetfile), sep = "\t", row.names = FALSE, quote = FALSE)
  plot(
    fviz_dend(res, rect = TRUE, cex = 0.5)
  )
  plot(
    fviz_silhouette(res)
  )
  dev.off()
}

# -----
# --
# -----

cluster_sampls <- function(esetfile, gencol, sampattern, metadat, samplecol, inputfolder, outfolder) {
  eset <- fread(paste0(inputfolder, "/", esetfile))
  esemat <- data.matrix(dplyr::select(eset, matches(sampattern)))
  rownames(esemat) <- eset[[gencol]]
  esemat <- t(esemat)
  metadat <- dplyr::filter(metadat, Sample_geo_accession %in% rownames(esemat))
  
  # best number of clusters
  print(dim(esemat))
  #print(which(esemat < 0))
  print(esemat[1:5, 1:10])
  
  nclusts <- tryCatch(
    {
      NbClust(esemat, distance = "euclidean",
              min.nc = 2, max.nc = 10,
              method = "average", index ="all")
      
      pdf(paste0(outfolder, "/best_ncluster_", gsub(".tsv", ".pdf", esetfile)), width = 10, height = 10)
      plot(
        fviz_nbclust(res.nbclust) +
          theme_minimal()
      )
      dev.off()
    }, error = function(e) {
      print(e)
      return(c(3, 4))
    }
  )
  
  
  # sample distance
  draw_sample_dist(esemat, esetfile, outfolder)
  
  # cluster using k-means
  do_kmeans(esemat, esetfile, outfolder, metadat, samplecol, nclusts)
  
  # cluster using clara
  do_clara(esemat, esetfile, outfolder, metadat, samplecol, nclusts)
  
  # cluster using heriarquical clustering
  do_hc(esemat, esetfile, outfolder, nclusts)
}

# -----
# -- Inputs
# -----

inputfolder <- "../data"
esetfiles <- list.files("../result", pattern = "^MNN_corrected_")

metadat <- read.delim("../data/metadata_blood_sjia_pjia_normal_samples.tsv")

lapply(esetfiles, function(esetfile) {
  print(esetfile)
  cluster_sampls(esetfile, "Symbol", "GSM",  metadat, "Sample_geo_accession", inputfolder, "../result/sampl_clust")
})
# }, mc.cores = 3)


# visualize for only the best ds + fs + cls combination

#res <- fread("../result/entropy_clusters.tsv")
#res <- melt.data.table(res, 'clust_ds')

#esetfiles <- unique(c(
#  arrange(filter(res, grepl("3clusters", clust_ds) & variable == "ShannonEntropy"), value)[["clust_ds"]][1:5],
#  arrange(filter(res, grepl("3clusters", clust_ds) & variable == "entropy"), value)[["clust_ds"]][1:5]
#)
#)

#esetfiles <- gsub(".+_\\d+clusters_", "", esetfiles)
#esetfiles <- paste0(esetfiles, ".tsv")

