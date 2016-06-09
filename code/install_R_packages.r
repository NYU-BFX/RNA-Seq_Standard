#!/usr/bin/env Rscript

list.of.packages=c("DESeq2","GenomicFeatures","dplyr","BiocParallel","pheatmap","RColorBrewer","ggplot2","ReportingTools","hwriter","reshape2","preprocessCore","cowplot","diagram")
missing.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
for (p in missing.packages){
 if(!file.exists(Sys.getenv("R_LIBS_USER"))){
    system(paste0("mkdir ",Sys.getenv("R_LIBS_USER")))
 }
 source("https://bioconductor.org/biocLite.R") 
 biocLite(p,lib=Sys.getenv("R_LIBS_USER"))
}