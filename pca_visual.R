# Visualize sample distribution
# PCA with metadata labels

library(data.table)
library(dplyr)
library(ggfortify)

# Input data

# --------
# Function
# --------
filetypes <- c("raw_eset", "quantile_eset", "quantile_eset_naive_combat", 
               "yugene_eset", "zscore_eset", "UPC_eset", "quantile_eset_superv_combat")

#filetypes <- "quantile_eset_superv_combat"
# Input eset data
filetype <- filetypes[3]
for(filetype in as.vector(filetypes)[1:43]) {
  esetfile <- paste0("blood_sjia_pjia_", filetype, ".tsv")
  eset <- read.delim(paste0("../data/", esetfile))
  
  # metadata
  md <- read.delim("../data/metadata_blood_sjia_pjia_normal_samples.tsv")
  md[["manufacturer"]] <- ifelse(md$Sample_platform_id == "GPL570", "Affymetrix", "Illumina")
  
  # variables of interest
  info.var <- colnames(md)[c(5, 8, 10, 11, 17, 18, 20, 21, 37, 38, 40, 57, 58)]
  print(head(info.var))
  
  # setting eset format to the PCA visualization
  em <- dplyr::select(eset, matches("GSM"))
  rownames(em) <- eset$Symbol
  em <- t(em)
  colnames(em) <- gsub("\\.|\\-", "_", colnames(em))
  
  dat <- cbind(data.frame("Sample_geo_accession" = rownames(em)), em)
  print(dat[1:5, 1:10])
  dat <- data.frame(merge(md, dat, by = "Sample_geo_accession"))
  rownames(dat) <- as.character(dat[, 'Sample_geo_accession'])
  print(dat[1:5, 1:10])
  
  # writing .PDF visualizations
  pdf(paste0("../result/pca_blood_sjia_pjia_", filetype, ".pdf"))
  for(v in info.var) {
    plot(
    autoplot(prcomp(em), data = dat, colour = v)
    )
  }
  dev.off()
}
