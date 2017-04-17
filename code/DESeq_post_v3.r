#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
dataset=args[1]
align_summary=args[2]
signature_list=args[3]
target_list=args[4]
load(dataset)
source(paste0(dirname(unlist(strsplit(commandArgs(trailingOnly = FALSE)[4],"="))[2]),"/DESeq_Comparison_Report.r"))

my.packages=c("DESeq2","pheatmap","RColorBrewer","ggplot2","ReportingTools","hwriter","reshape2","dplyr","preprocessCore","cowplot")
for (p in my.packages){
    suppressPackageStartupMessages(library(p,character.only=TRUE,verbose=FALSE))
}

toNumeric<-function(myMatrix){
  temp=apply(myMatrix,2,as.numeric)
  rownames(temp)=rownames(myMatrix)
  colnames(temp)=colnames(myMatrix)
  myMatrix=temp
  return(myMatrix)
}

fpkm_table=toNumeric(fpkm_table)
norm_count=toNumeric(norm_count)

longdir=dirname(dataset)

if(grep("/",longdir)){outdir=basename(longdir)}else{outdir=longdir}
dir.create(longdir,showWarnings=F)
htmlRep <- HTMLReport(shortName = outdir, title="RNA-Seq Experiment Standard Report", reportDirectory = longdir)


# signature as goldern standard to be tested
signature=read.table(signature_list)[,1]
mylist=""
mydata=NULL
for (i in signature){mylist=c(mylist,grep(paste("^",i,"$|:",i,"$",sep=""), rownames(fpkm_table), ignore.case = T))}
signature=rownames(fpkm_table)[as.numeric(mylist[-1])]

target=read.table(target_list)[,1]
mylist=""
mydata=NULL
for (i in target){mylist=c(mylist,grep(paste("^",i,"$|:",i,"$",sep=""), rownames(fpkm_table), ignore.case = T))}
target=rownames(fpkm_table)[as.numeric(mylist[-1])]

# plot the alignment summary bar plot
align_reads <- function(x){
  x = x %>% filter(V2!="Uniquely_Mapped") %>% mutate(V4=ifelse(grepl("NonDuplicated|Duplicated",V2, ignore.case=T),"Uniquely_Mapped",as.character(V2))) %>% arrange(desc(V1)) %>% mutate(V3=V3/1000000)
  x$V2 <- factor(x$V2, levels = c("NonDuplicated","Duplicated","Multi_Mapped","Unmapped"))
  x$V4 <- factor(x$V4, levels = c("Uniquely_Mapped","Multi_Mapped","Unmapped"))
  x=x[order(x$V1,x$V2),]
  myplot <-
    ggplot(x,aes(x=V1,y=V3,fill=V2)) +
    geom_bar(stat="identity",position=position_stack(vjust = 1, reverse = T)) +
    scale_fill_manual(values=c("limegreen","forestgreen","orchid","indianred1")) +
    ggtitle("B. Alignment Category Reads") +
    xlab("") + ylab("Number of Reads in Millions") +
    theme(axis.title = element_text(size=8),axis.text.y = element_blank(),axis.text.x = element_text(size=8),plot.title = element_text(size=8),legend.title=element_blank(),legend.text=element_text(size=8),axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
    coord_flip()
    return(myplot)
}

align_PCT <- function(x){
  x = x %>% filter(V2!="Uniquely_Mapped") %>% group_by(V1,add=F) %>% mutate(V3=V3/sum(V3)) %>% mutate(V4=ifelse(grepl("NonDuplicated|Duplicated",V2, ignore.case=T),"Uniquely_Mapped",as.character(V2))) %>% arrange(desc(V1)) %>% mutate(V3=V3*100)
  x$V2 <- factor(x$V2, levels = c("NonDuplicated","Duplicated","Multi_Mapped","Unmapped"))
  x$V4 <- factor(x$V4, levels = c("Uniquely_Mapped","Multi_Mapped","Unmapped"))
  x=x[order(x$V1,x$V2),]
    myplot <- ggplot(x,aes(x=V1,y=V3,fill=V2)) +
    geom_bar(stat="identity",position=position_stack(vjust = 1, reverse = T)) +
    ggtitle("A. Alignment Category Percentage") +
    scale_fill_manual(values=c("limegreen","forestgreen","orchid","indianred1")) +
    xlab("") + ylab("Percentage of Reads") +
    theme(axis.title = element_text(size=8),axis.text.y = element_text(size=8,angle=0),axis.text.x = element_text(size=8),plot.title = element_text(size=8)) +
    guides(fill=FALSE) +
    coord_flip()
    return(myplot)
}

dir.create(paste(longdir,"/figures",outdir,sep=""))
pdf(paste0(longdir,"/figures",outdir,"/Align_Summary.pdf",sep=""),width = 7,height=6,onefile = F)
plot_grid(align_PCT(read.table(align_summary)),align_reads(read.table(align_summary)))
dev.off()
png(filename=paste0(longdir,"/figures",outdir,"/Align_Summary.png",sep=""),width = 7,height=6,units = "in",res=300)
plot_grid(align_PCT(read.table(align_summary)),align_reads(read.table(align_summary)))
dev.off()
himg <- hwriteImage(paste("figures",outdir,"/Align_Summary.png",sep=""),link=paste("figures",outdir,"/Align_Summary.pdf",sep=""))
publish(hwrite(himg, br=TRUE), htmlRep,name=paste("figures",outdir,"/Align_Summary",sep=""))

# compute number of reads total for each samples
sample_summary<-colData(dds)[,1:2]%>%data.frame()%>%left_join(data.frame(sample=rownames(data.frame(colSums(assay(dds)))),NumofReads=colSums(assay(dds))),by=c("sampleName"="sample"))
publish(paste("Sample Group Information ",sep=""), htmlRep)
publish(sample_summary, htmlRep)

# Another transformation, the variance stabilizing transformation [@Anders2010Differential], is discussed alongside the rlog in the DESeq2 vignette
# 
# calculate the Euclidean distance between samples on the transposed matrix of log-transformed data, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns
# plot the sample distance
splDistPlot <- function(mymatrix){
  sampleDists <- dist( t( mymatrix ) )
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- colnames(mymatrix)
  colnames(sampleDistMatrix) <- colnames(mymatrix)
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  png(filename=paste(longdir,"/figures",outdir,"/Sample_Distance.png",sep=""))
  pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, main="Distance Between Samples")
  dev.off()
  pdf(file=paste(longdir,"/figures",outdir,"/Sample_Distance.pdf",sep=""),onefile=F)
  pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, main="Distance Between Samples")
  dev.off()
}
  splDistPlot(assay(rld))
  himg <- hwriteImage(paste("figures",outdir,"/Sample_Distance.png",sep=""),link=paste("figures",outdir,"/Sample_Distance.pdf",sep=""))
  publish(hwrite(himg, br=TRUE), htmlRep,name=paste("figures",outdir,"/Sample_Distance",sep=""))

# Another option for calculating sample distances is to use the Poisson Distance [@Witten2011Classification], implemented in the PoiClaClu package. This measure of dissimilarity between counts also takes the inherent variance structure of counts into consideration when calculating the distances between samples. The PoissonDistance function takes the original count matrix (not normalized) with samples as rows instead of columns, so we need to transpose the counts in dds
#poisd <- PoissonDistance(t(counts(dds)))
#samplePoisDistMatrix <- as.matrix( poisd$dd )
#rownames(samplePoisDistMatrix) <- rld$sampleName
#colnames(samplePoisDistMatrix) <- rld$sampleName
#pheatmap(samplePoisDistMatrix,clustering_distance_rows=poisd$dd,clustering_distance_cols=poisd$dd,col=colors)

# PCA on rlog transofrmed data
data <- plotPCA(rld, intgroup = c("grp","sampleName"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
myplot=ggplot(data, aes(PC1, PC2, color=grp,label=sampleName)) + geom_point(size=3) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  labs(title = "PCA Plot")+geom_text(size=3.5,hjust=0.5, vjust=1.5)
publish(myplot, htmlRep,name="PCA Plot on Known Genes")


# Plot FPKM for signature gene
mydata=fpkm_table[signature,]
mydata.colnames=c("Gene",colnames(mydata))
mydata=data.frame(cbind(Gene=rownames(mydata),mydata))
colnames(mydata)=mydata.colnames
mydata=melt(mydata,id.vars=1, factorsAsStrings=F)
mydata<-mydata%>%left_join(data.frame(colData(dds)),by=c("variable"="sampleName"))%>%select(Gene,variable,value,grp,sizeFactor)
colnames(mydata)=c("gene","sampleName","count","grp","sizeFactor")
mydata<-mydata%>%arrange(gene)
mydata=mydata[,c(1,3,4,2)]
mydata<-mydata %>% mutate(count=as.numeric(as.character(count)),sampleName=as.factor(sampleName))

geneBarplot <- function(mytitle,data){
  if(dim(data)[1]!=0){
    myplot <- ggplot(data,aes(x=sampleName,y=count,fill=grp)) + 
      geom_bar(stat="identity",position="dodge")  + 
      facet_grid(gene~.)+
      theme(strip.text.y = element_text(size=8, angle=0), legend.title=element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank()) +
      labs(title = paste("FPKM for ",mytitle," Signature Genes",sep="")) +
      guides(color=FALSE) + 
      xlab("Gene") + ylab("FPKM") +
      coord_flip()
    return(myplot)
  }else{return(NULL)}
}


geneBoxplot <- function(mytitle,data){
  if(dim(data)[1]!=0){
    myplot <- ggplot(data, aes(x=grp, y=count)) +
      geom_boxplot(aes(fill=grp,color=grp)) +
      facet_grid(gene~.) +
      theme(strip.text.y = element_text(size=8, angle=0), legend.title=element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
      labs(title = paste("FPKM for ", mytitle," Signature Genes",sep="")) +
      guides(color=FALSE) +
      xlab("Gene") + ylab("FPKM") +
      coord_flip()+ 
      stat_summary(geom = "crossbar", fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) })
    return(myplot)
  }else{return(NULL)}
}

highexp <- mydata %>% group_by(gene) %>% filter(max(count)>100) %>% ungroup
mytitle="High Exp"
myplot <- geneBoxplot(mytitle,highexp)
if(!is.null(myplot)){publish(myplot, htmlRep,name=paste("FPKM for ",mytitle," Signature Genes",sep=""))}

lowexp <- mydata %>% group_by(gene) %>% filter(max(count)<100) %>% ungroup
mytitle="Low Exp"
myplot <- geneBoxplot(mytitle,lowexp)
if(!is.null(myplot)){publish(myplot, htmlRep,name=paste("FPKM for ",mytitle," Signature Genes",sep=""))}

publish("FPKM Table:", htmlRep)
publish(hwrite("Download", link = paste("../fpkm_table.txt",sep="")), htmlRep)
publish("          ", htmlRep)
publish("          ", htmlRep)
publish("          ", htmlRep)
publish("          ", htmlRep)
publish("          ", htmlRep)
publish("          ", htmlRep)
publish("Normalized Count Table:", htmlRep)
publish(hwrite("Download", link = paste("../norm_count.txt",sep="")), htmlRep)
publish("          ", htmlRep)
publish("          ", htmlRep)
publish("          ", htmlRep)
publish("          ", htmlRep)
publish("          ", htmlRep)
publish("          ", htmlRep)
publish("Raw Count Table:", htmlRep)
publish(hwrite("Download", link = paste("../raw_count.txt",sep="")), htmlRep)

# Volcano Plot
#plot(res$log2FoldChange,-log(res$padj,10),main="Volcano Plot of FDR-adjusted P-values",pch=20,cex=0.4,xlab="log2(FC of mean of normalized counts)",ylab="-log10(adj P-value)")

# Dispersion Plot. First, gene-wise MLEs are obtained using only the respective gene’s data (black dots). Then, a curve (red) is fit to the MLEs to capture the overall trend of dispersion-mean dependence. This fit is used as a prior mean for a second estimation round, which results in the final MAP estimates of dispersion (arrow heads). This can be understood as a shrinkage (along the blue arrows) of the noisy gene-wise estimates toward the consensus represented by the red line. The black points circled in blue are detected as dispersion outliers and not shrunk toward the prior (shrinkage would follow the dotted line). 
png(filename=paste(longdir,"/figures",outdir,"/Dispersion_Plot.png",sep=""))
plotDispEsts(dds, main="Dispersion Plot")
dev.off()
pdf(file=paste(longdir,"/figures",outdir,"/Dispersion_Plot.pdf",sep=""))
plotDispEsts(dds, main="Dispersion Plot")
dev.off()
himg <- hwriteImage(paste("figures",outdir,"/Dispersion_Plot.png",sep=""),link=paste("figures",outdir,"/Dispersion_Plot.pdf",sep=""))
publish(hwrite(himg, br=TRUE), htmlRep,name=paste("figures",outdir,"/Dispersion_Plot",sep=""))

if(dim(comparisons)[1]<=36){
publish(paste("Analysis for Different Comparison Groups ",sep=""), htmlRep)
for (i in 1:dim(comparisons)[1]){
  preprocessing="rld"
  print(comparisons[i,])
  myreport=summarize_comparison(longdir,signature,target,"DE",comparisons[i,],preprocessing,fpkm_table)
  publish(hwrite(myreport, link = paste(myreport,"/","report.html",sep="")), htmlRep)
}
}
finish(htmlRep)


