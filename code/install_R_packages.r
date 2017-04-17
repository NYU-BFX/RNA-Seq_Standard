#!/usr/bin/env Rscript

list.of.packages=c("EDASeq","ggthemes","TxDb.Hsapiens.UCSC.hg19.knownGene","TCGAbiolinks","tools","DESeq2","GenomicFeatures","dplyr","BiocParallel","pheatmap","RColorBrewer","ggplot2","ReportingTools","hwriter","reshape2","preprocessCore","cowplot","diagram","GGally","tidyr","animation")
missing.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
for (p in missing.packages){
 if(!file.exists(Sys.getenv("R_LIBS_USER"))){
    system(paste0("mkdir -p ",Sys.getenv("R_LIBS_USER")))
 }
 source("http://bioconductor.org/biocLite.R") 
 biocLite(p,lib=Sys.getenv("R_LIBS_USER"),suppressUpdates=TRUE)
}
