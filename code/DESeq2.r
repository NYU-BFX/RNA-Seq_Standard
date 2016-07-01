#!/usr/bin/env Rscript
usage = "\
	Rscript DESeq2.r [GTF File] [Gene ID to Name Match File] [Group Info File] [Group_By Index] [RAW Count File]
"
Sys.time()
args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=5) { write(usage,stderr()); quit(save='no'); }

my.packages=c("DESeq2","GenomicFeatures","dplyr","BiocParallel","tools")
for (p in my.packages){
    suppressPackageStartupMessages(library(p,character.only=TRUE,verbose=FALSE))
}

sessionInfo()
register(MulticoreParam(Sys.getenv(x = "NSLOTS", unset = "6", names = NA)))


(gtffile <- file.path(args[1]))
id2name=as.matrix(read.table(args[2],header=F,row.names = 1))
id2name=as.matrix(id2name[order(rownames(id2name)),])
group=read.table(args[3],header=T,fill=T)
group_by=args[4]
outdir=file_path_sans_ext(basename(args[3]))
grp=factor(group[,group_by])
countsData <-read.table(args[5],header=T,row.names=1)
# raw count must be sorted by rowname before run DESeq2 because when matching with GTF file, it doesn't sort automaticlly and will cause mistake
countsData <- countsData[order(rownames(countsData)),]
countsData <- countsData[,order(colnames(countsData))]
sampleName=colnames(countsData)
sampleData <- data.frame(row.names=sampleName,group[order(group[,group_by]),group_by],sampleName)
colnames(sampleData) <- c("grp","sampleName")
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))
(ebg <- exonsBy(txdb, by="gene"))

# construct DESeq2 matrix from count data
dds <- DESeqDataSetFromMatrix(countsData, colData=data.frame(sampleData), design=~grp)

# add GTF file into DESeq2 Matrix
rowRanges(dds) <- ebg
dds <- DESeq(dds)
# Regularized-logarithm transformation
rld <- rlog(dds, blind=FALSE)
rownames(rld)=cbind(id2name[,1][match(rownames(rld),rownames(id2name))],rld)[,1]

dir.create(file.path(outdir,group_by),recursive=T,showWarnings=F)

# extract norm count
norm_count=counts(dds, normalized=T)
# add gene name
norm_count=data.frame(cbind(Gene=id2name[,1][match(rownames(norm_count),rownames(id2name))],norm_count))
# convert to Numeric for the count field, for ReportingTools
if(!file.exists(file.path(outdir,"norm_count.txt"))){write.table(norm_count,file=file.path(outdir,"norm_count.txt"),quote=F,sep="\t",row.names=F,col.names=T)}
rownames(norm_count)=norm_count$Gene
norm_count=norm_count[,-1]


# extract fpkm 
fpkm_table=fpkm(dds,robust=T)
# add gene name
fpkm_table=data.frame(cbind(Gene=id2name[,1][match(rownames(fpkm_table),rownames(id2name))],fpkm_table))
# convert to Numeric for the fpkm field, for ReportingTools
if(!file.exists(file.path(outdir,"fpkm_table.txt"))){write.table(fpkm_table,file=file.path(outdir,"fpkm_table.txt"),quote=F,sep="\t",row.names=F,col.names=T)}
rownames(fpkm_table)=fpkm_table$Gene
fpkm_table=fpkm_table[,-1]

comparisons=cbind("grp",t(combn(levels(colData(dds)$grp),2)))

ToCharacter=function(x){x}

for (i in 1:dim(comparisons)[1]){
  res=results(dds,contrast=comparisons[i,])
  id2name=as.matrix(id2name[order(rownames(id2name)),])
  res$symbol=apply(id2name,1,ToCharacter)
  mytable <- res %>% data.frame() %>% select(symbol,log2FoldChange,padj,pvalue) %>% mutate(log2FoldChange=round(log2FoldChange,digits=4),padj=round(padj,digits=4),pvalue=round(pvalue,digits=4))
  colnames(mytable)=c("_ID_",paste(comparisons[i,2],"_vs_",comparisons[i,3],"_LFC",sep=""),paste(comparisons[i,2],"_vs_",comparisons[i,3],"_FDR",sep=""),paste(comparisons[i,2],"_vs_",comparisons[i,3],"_PVAL",sep=""))
  write.table(mytable,file=paste(outdir,"/",group_by,"/",comparisons[i,2],"_vs_",comparisons[i,3],"_DE.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
}

save.image(paste0(outdir,"/",group_by,"/",group_by,"_Rscript.RData"))
Sys.time()
