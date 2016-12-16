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

ReadAsMatrix<-function(filename){
  myMatrix <-as.matrix(read.table(filename,row.names=1))
  colnames(myMatrix)=myMatrix[1,]
  myMatrix=myMatrix[-1,]
  temp=apply(myMatrix,2,as.numeric)
  rownames(temp)=rownames(myMatrix)
  colnames(temp)=colnames(myMatrix)
  myMatrix=temp
  return(myMatrix)
}

(gtffile <- file.path(args[1]))
id2name=as.matrix(read.table(args[2],header=F,row.names = 1))
id2name=as.matrix(id2name[order(rownames(id2name)),])
group=read.table(args[3],header=T,fill=T)
group_by=args[4]
outdir=file_path_sans_ext(basename(args[3]))
grp=factor(group[,group_by])
countsData=ReadAsMatrix(args[5])
# raw count must be sorted by rowname before run DESeq2 because when matching with GTF file, it doesn't sort automaticlly and will cause mistake
countsData <- countsData[order(rownames(countsData)),]
countsData <- countsData[,order(colnames(countsData))]
sampleName=colnames(countsData)
sampleData <- group[,c(which(colnames(group)==group_by),1)]
rownames(sampleData)=group[,1]
colnames(sampleData) <- c("grp","sampleName")
sampleData <- sampleData[order(rownames(sampleData)),]
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))
(ebg <- exonsBy(txdb, by="gene"))

# construct DESeq2 matrix from count data
dds <- DESeqDataSetFromMatrix(countsData, colData=data.frame(sampleData), design=~grp)
# add GTF file into DESeq2 Matrix
rowRanges(dds) <- ebg

dir.create(file.path(outdir,group_by),recursive=T,showWarnings=F)

# calculate size factor
dds <- estimateSizeFactors(dds)

# extract norm count
raw_count=counts(dds)
# add gene name
rownames(raw_count)=id2name[,1][match(rownames(raw_count),rownames(id2name))]
# convert to Numeric for the count field, for ReportingTools
if(!file.exists(file.path(outdir,"raw_count.txt"))){
	write.table(matrix(c("Gene",colnames(raw_count)),nrow=1),file=file.path(outdir,"raw_count.txt"),quote=F,sep="\t",col.names=F,row.names=F)
	write.table(raw_count,file=file.path(outdir,"raw_count.txt"),quote=F,sep="\t",row.names=T,col.names=F,append=T)
}

# extract norm count
norm_count=counts(dds, normalized=T)
# add gene name
rownames(norm_count)=id2name[,1][match(rownames(norm_count),rownames(id2name))]
# convert to Numeric for the count field, for ReportingTools
if(!file.exists(file.path(outdir,"norm_count.txt"))){
	write.table(matrix(c("Gene",colnames(norm_count)),nrow=1),file=file.path(outdir,"norm_count.txt"),quote=F,sep="\t",col.names=F,row.names=F)
	write.table(norm_count,file=file.path(outdir,"norm_count.txt"),quote=F,sep="\t",row.names=T,col.names=F,append=T)
}


# extract fpkm 
fpkm_table=fpkm(dds,robust=T)
# add gene name
rownames(fpkm_table)=id2name[,1][match(rownames(fpkm_table),rownames(id2name))]
# convert to Numeric for the fpkm field, for ReportingTools
if(!file.exists(file.path(outdir,"fpkm_table.txt"))){
	write.table(matrix(c("Gene",colnames(fpkm_table)),nrow=1),file=file.path(outdir,"fpkm_table.txt"),quote=F,sep="\t",col.names=F,row.names=F)
        write.table(fpkm_table,file=file.path(outdir,"fpkm_table.txt"),quote=F,sep="\t",row.names=T,col.names=F,append=T)
}

comparisons=cbind("grp",t(combn(levels(colData(dds)$grp),2)))

if(dim(comparisons)[1]<=36){

# clean large objects before running DESeq() in parallel to help avoid memory issue (R's internal memory garbage collector, gc(), may be making copies of the large objects in the environment in all of the workers.)
rm(id2name,countsData,txdb,ebg,raw_count,norm_count,fpkm_table)
dds <- DESeq(dds, parallel = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = Sys.getenv(x = "NSLOTS", unset = "6", names = NA)))
# reload cleaned objects
id2name=as.matrix(read.table(args[2],header=F,row.names = 1))
countsData=ReadAsMatrix(args[5])
countsData <- countsData[order(rownames(countsData)),]
countsData <- countsData[,order(colnames(countsData))]
raw_count=counts(dds)
sampleName=colnames(countsData)
norm_count=counts(dds, normalized=T)
rownames(norm_count)=id2name[,1][match(rownames(norm_count),rownames(id2name))]
fpkm_table=fpkm(dds,robust=T)
rownames(fpkm_table)=id2name[,1][match(rownames(fpkm_table),rownames(id2name))]


# Regularized-logarithm transformation
rld <- rlog(dds, blind=FALSE)
rownames(rld)=id2name[,1][match(rownames(rld),rownames(id2name))]


ToCharacter=function(x){x}


for (i in 1:dim(comparisons)[1]){
  res=results(dds,contrast=comparisons[i,])
  id2name=as.matrix(id2name[order(rownames(id2name)),])
  res$symbol=apply(id2name,1,ToCharacter)
  mytable <- res %>% data.frame() %>% select(symbol,log2FoldChange,padj,pvalue) %>% mutate(log2FoldChange=round(log2FoldChange,digits=4),padj=round(padj,digits=4),pvalue=round(pvalue,digits=4))
  colnames(mytable)=c(".ID_",paste(comparisons[i,2],"_vs_",comparisons[i,3],"_LFC",sep=""),paste(comparisons[i,2],"_vs_",comparisons[i,3],"_FDR",sep=""),paste(comparisons[i,2],"_vs_",comparisons[i,3],"_PVAL",sep=""))
  write.table(mytable,file=paste(outdir,"/",group_by,"/",comparisons[i,2],"_vs_",comparisons[i,3],"_DE.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
}
}

save.image(paste0(outdir,"/",group_by,"/",group_by,"_Rscript.RData"))
Sys.time()
