# unsupervised clustering methods

# semiunsupervised method

## Matching Mutual Nearest Neighbors
#BiocManager::install("scran")
library(scran)
library(data.table)
library(dplyr)
library(scater)


file2esets <- function(esetfile, inputfolder, sampattern, metafile, batch) {
  eset <- fread(paste0(inputfolder, esetfile))
  meta <- fread(paste0(inputfolder, metafile))
  bs <- unique(meta[[batch]])
  setNames(
    lapply(bs, function(b) {
      bsampls <- dplyr::filter(meta, batch == b)[["Sample_geo_accession"]]
      select(eset, matches(paste(c("Symbol", bsampls), collapse = "|")))
    }), bs)
}

#sapply(file2esets(fs[1], inputfolder, sampattern, metafile, batch), dim)

eset2obj <- function(eset, sampattern) {
  emat <- data.matrix(select(eset, matches(sampattern)))
  rownames(emat) <- eset$Symbol
  sce <- SingleCellExperiment(list(counts = emat, logcounts = emat))
  sce
}



do_mnn_corr <- function(esetfiles, inputfolder, sampattern, metafile, batch, ncores) {
  #mclapply(esetfiles, function(esetfile) {
  lapply(esetfiles, function(esetfile) {
    if(esetfile %in% c("blood_sjia_pjia_UPC_eset_sjia_sex_difF_REM_deg.tsv", 
                       "blood_sjia_pjia_yugene_eset_sjia_sex_difF_REM_deg.tsv",
                       "blood_sjia_pjia_zscore_eset_sjia_sex_difF_REM_deg.tsv")) {
      NULL
    } else {
      print(esetfile)
      # get merged eset file and split it in several tables based on the GSE accessions
      print("starting...")
      esets <- file2esets(esetfile, inputfolder, sampattern, metafile, batch)
      print("esetfile splitted...")
      # from tables to SingleCellExpression objects, deprecating the GSE58667 that contains just 4 samples
      print(names(esets))
      esetobjs <- lapply(esets[c("GSE55319", "GSE58667", "GSE80060_2011-03", "GSE80060", "GSE17590")], 
                         function(e) eset2obj(e, sampattern))
      
      # MNN batch correction 
      set.seed(224033911)
      print("correcting...")
      mnn.out <- do.call(fastMNN, c(esetobjs, list(k = 10, d = 50, approximate = TRUE, auto.order=TRUE)))
      print("corrected...")
      
      # Assembling original datset 
      omat <- Reduce(cbind, lapply(esetobjs, function(eo) logcounts(eo)))
      
      # writing MNN corrected matrix
      print("writing MNN corrected eset...")
      cor.exp <- tcrossprod(mnn.out$rotation, mnn.out$corrected)
      colnames(cor.exp) <- colnames(omat)
      cor.exp <- cbind(data.frame("Symbol" = rownames(omat)), cor.exp)
      write.table(cor.exp, paste0("../result/MNN_corrected_", gsub("\\.tsv", "", esetfile), ".tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      
      print("visualizing...")
      # visualize cluster differences among the original and corrected datasets through tSNE plots
      sce <- SingleCellExperiment(list(counts = omat, logcounts = omat))
      reducedDim(sce, "MNN") <- mnn.out$corrected
      sce$Batch <- as.character(mnn.out$batch)
      
      # Including disease metadata 
      msce <- data.frame("Sample_geo_accession" = rownames(sce@colData)) 
      
      # Reference metadata
      meta <- read.delim(paste0(inputfolder, metafile))
      rownames(meta) <- meta$Sample_geo_accession
      msce <- merge(msce, meta, by = "Sample_geo_accession")
      sce$Disease <- as.character(msce$Characteristics..disease.)
      sce$scandate <- as.character(msce$scandate_yymm)
      
      set.seed(100)
      # Using irlba to set up the t-SNE, for speed.
      # visualization of the original data
      osce <- runPCA(sce, ntop=Inf, method="irlba")
      osce <- runTSNE(osce, use_dimred="PCA")
      ot <- plotTSNE(osce, colour_by="Batch") + ggtitle("Original")
      dot <- plotTSNE(osce, colour_by="Disease") + ggtitle("Original")
      sdot <- plotTSNE(osce, colour_by="scandate") + ggtitle("scandate")
      # Visualizartion of the MNN transformed data
      set.seed(100)
      csce <- runTSNE(sce, use_dimred="MNN")
      ct <- plotTSNE(csce, colour_by="Batch") + ggtitle("Corrected")
      dct <- plotTSNE(csce, colour_by="Disease") + ggtitle("Corrected")
      sdct <- plotTSNE(csce, colour_by="scandate") + ggtitle("scandate")
      # Clustering
      # The aim is to use the SNN graph to perform clustering of cells 
      # via community detection algorithms in the igraph package
      snn.gr <- buildSNNGraph(sce, use.dimred = "MNN")
      clusters <- igraph::cluster_walktrap(snn.gr)
      table(clusters$membership, sce$Batch)
      csce$Cluster <- factor(clusters$membership)
      # Ploting cluster on tSNE plot
      cls <- plotTSNE(csce, colour_by="Cluster")
           
      # Ploting batch, disease and ad-hoc clustering on the original and MNN transformed data
      pdf(paste0("../result/tSNE_MNN_corrected_", gsub("\\.tsv", "", esetfile), ".pdf"), width = 15, height = 15)
        multiplot(ot, ct, sdot, cls, dot, dct, sdct, cls, cols=2)
      dev.off()
    #}, mc.cores = ncores)
    }
  })
}

# Inputs

inputfolder <- "../data/"
esetpattern <- "^blood_sjia_pjia_"
sampattern <- "GSM"
metafile <- "metadata_blood_sjia_pjia_normal_samples.tsv"
batch <- "batch" #"GSE_accession_wfound"
esetfiles <- list.files(path = inputfolder, pattern = esetpattern)
ncores <- 3

do_mnn_corr(esetfiles[1:4], inputfolder, sampattern, metafile, batch, ncores)

