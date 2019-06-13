# Assess gene expression dataset variance
# How does metadata explain gene expression variation?

library(data.table)
library(dplyr)
library(pvca)
library(arrayQualityMetrics)
library(Biobase)
library(ggplot2)

# ------------
# Functions
# ------------

# replace points with _ from string vectors
cleaname <- function(names) {
  gsub("\\.", "_", names)
}

# write .pdf barplot with pvca results
draw_pvca <- function(pvcares, filename) {
  pdf(filename)
  plot(
    ggplot(pvcares, aes(x = label, y = dat)) +
      geom_bar(stat = "identity") +
      scale_x_discrete(limits = dplyr::arrange(pvcares, dat)[["label"]]) +
      theme_classic() +
      coord_flip()
  )
  dev.off()
} 

# running either or both pvca and AQM
do_QC <- function (exprss, sampattern, pheno, info.var, inputfolder, outfolder, pvca, AQM) {
  # exprss: eset file name
  # sampattern: common string among variables names that identifies samples
  # pheno: pheno file name
  # info.var: string vector with pheno variables to consider
  # inputfolder: folder where express and pheno are located
  # outfolder: where results will be saved
  # pvca: logic, if TRUE pvca will be performed
  # AQM: logic, if TRUE AQM will be performed
  phenodata <- read.delim(paste(inputfolder, "/", pheno, sep = ""))
  exp.data <- read.delim(paste(inputfolder, "/", exprss, sep = ""))
  #exp.data <- read.delim(paste(outfolder, "/", exprss, sep = ""))
  
  # cleanames
  colnames(phenodata) <- cleaname(colnames(phenodata))
  colnames(exp.data) <- cleaname(colnames(exp.data))
  
  # defining eset object
  emat <- data.matrix(dplyr::select(exp.data, matches(sampattern)))
  print(dim(emat))
  if(nrow(emat) < 50) {
    NULL
  } else {
    rownames(emat) <- exp.data$Symbol
    eset <- ExpressionSet(emat)
    
    rownames(phenodata) <- phenodata$Sample_geo_accession
    # filtering phenodata based on samples in the eset
    print("this sample is absent in the expression set:")
    print(rownames(pData(eset))[which(!rownames(pData(eset)) %in% colnames(exprs(eset)))])
    phenodata <- dplyr::filter(phenodata, Sample_geo_accession %in% colnames(exprs(eset)))
    pData(eset) <- phenodata
    # printing info on upoad data
    print(eset)
    print(exprs(eset)[1:3, 1:10])
    print(head(pData(eset), 3))
    print(info.var)
    
    if (pvca == TRUE) {
      pvcaObj <- pvcaBatchAssess(eset, info.var, 0.6)
      print("pvca run!")
      data_var_expl = data.frame("dat" = as.vector(pvcaObj$dat), "label" = pvcaObj$label)
      print(head(data_var_expl, 3))
      out.file.name <- paste(outfolder, "/", "var_expl_", exprss, sep = "")
      print(out.file.name)
      write.table(data_var_expl, out.file.name, quote = FALSE, sep = "\t", row.names = FALSE)
      draw_pvca(data_var_expl, gsub("\\.tsv", ".pdf", out.file.name))
    }
    if (AQM == TRUE) {
      arrayQualityMetrics(eset, 
                          outdir = paste(outfolder, "/", "AQM_", gsub("\\.tsv", "", exprss), sep = ""), 
                          force = TRUE, intgroup = info.var)
    }
  }
}


# Input data
#filetype <- "yugene_eset" #"zscore_eset" #"yugene_eset" #"UPC_eset"
#filetypes <- sapply(c("raw_eset", "quantile_eset", "quantile_eset_naive_combat", 
#                     "quantile_eset_superv_combat", "zscore_eset", 
#                      "yugene_eset", "UPC_eset"), function(x) {
#  sapply(c("CEMi", "jia_citedinpaper", "sjia_systrev_deg1", "REM_sjia_deg", 
#           "meta_coexp_sjia_normal", "nonzero", "PCA", # "sjia_sex_difF_REM_deg"           
#           "sjia_systrev_deg1", "sPCA1000", "var25", "var50"), function(y) {
#               paste(x, y, sep = "_")
#             })
#})

filetypes <- c("raw_eset", "quantile_eset", "quantile_eset_naive_combat", 
               "quantile_eset_superv_combat", "zscore_eset", 
               "yugene_eset", "UPC_eset")

class(filetypes)

# meta variables
md <- read.delim("../data/metadata_sjia_pjia_samples.tsv")

info.var <- colnames(md)[c(2, 3, 9)]
info.var <- cleaname(info.var)
info.var

#mclapply(sapply(filetypes, function(f) f)[c(30, 32, 34:35)], function(filetype) {
#  print(filetype)
#  esetfile <- paste0("blood_sjia_pjia_", filetype, ".tsv")
#  print(esetfile)
#  do_QC(esetfile, "GSM", "metadata_sjia_pjia_samples.tsv", info.var, "../data", "../result", TRUE, FALSE)
#}, mc.cores = 4)

files <- list.files(path = "../result/", ".tsv")
files <- files[grep("^MNN_corrected_", files)]

for(esetfile in files) {
  #esetfile <- paste0("blood_sjia_pjia_", filetype, ".tsv")
  print(esetfile)
  do_QC(esetfile, "GSM", "metadata_sjia_pjia_samples.tsv", info.var, "../result", "../result", TRUE, FALSE)
}

for(filetype in as.vector(filetypes)[c(33:36)]) {
  print(filetype)
  esetfile <- paste0("MNN_corrected_blood_sjia_pjia_", filetype, ".tsv")
  print(esetfile)
  do_QC(esetfile, "GSM", "metadata_sjia_pjia_samples.tsv", info.var, "../data", "../result", TRUE, FALSE)
}
