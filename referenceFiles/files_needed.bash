#!/bin/bash


# genome.fa
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

twoBitToFa hg19.2bit genome.fa

samtools index genome.fa

# gene.gtf
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

gunzip gencode.v19.annotation.gtf.gz

# chromosome sizes
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
mv hg19.chrom.sizes chromInfo.txt

# STAR Index
cd ../
code/generate_STAR-Index 6 params referenceFiles/STAR_Reference

