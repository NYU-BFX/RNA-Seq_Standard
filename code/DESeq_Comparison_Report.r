#!/usr/bin/env Rscript

my.packages=c("DESeq2","pheatmap","RColorBrewer","ggplot2","ReportingTools","hwriter","reshape2","dplyr","preprocessCore")
for (p in my.packages){
    suppressPackageStartupMessages(library(p,character.only=TRUE,verbose=FALSE))
}

summarize_comparison=function(outdir,signature,target,DE,comparison,preprocessing,myMatrix){
  
  htmlRep1 <- HTMLReport(shortName = "report", title=paste(comparison[2],"_vs_",comparison[3],sep=""), reportDirectory = paste(outdir,"/",comparison[2],"_vs_",comparison[3],sep=""))
  
  ToCharacter=function(x){x}
  res=results(dds,contrast = comparison)
  res$symbol=apply(id2name,1,ToCharacter)
  resPcutoff <- subset(res,padj<0.1)
  if(dim(resPcutoff)[1]==0){resPcutoff <- subset(res,pvalue<0.05)}
  resOrderedDF <- resPcutoff[order(abs(resPcutoff$log2FoldChange),decreasing = T),]
  if(dim(resOrderedDF)[1]<1){return(paste("No Significant Differnece Between ",comparison[2],"_vs_",comparison[3],sep=""))}
  
  sampleSheet=colData(dds)[(colData(dds)$grp==comparison[2]|colData(dds)$grp==comparison[3]),1:2]%>%data.frame()%>%left_join(data.frame(sample=rownames(data.frame(colSums(assay(dds)))),NumofReads=colSums(assay(dds))),by=c("sampleName"="sample"))
  publish(paste("Sample Group Information ",sep=""), htmlRep1)
  publish(sampleSheet, htmlRep1)
  
  prefix="Signature_Genes"
  
  markMAplot(outdir,signature,comparison,prefix,res)  
  himg <- hwriteImage(paste(comparison[2],"_vs_",comparison[3],"/",prefix,"_MA_plot.png",sep=""),link=paste(comparison[2],"_vs_",comparison[3],"/",prefix,"_MA_plot.pdf",sep=""))
  publish(hwrite(himg, br=TRUE), htmlRep,name=paste("MA Plot ", comparison[2],"_vs_",comparison[3],sep=""))
  himg <- hwriteImage(paste(prefix,"_MA_plot.png",sep=""),link=paste(prefix,"_MA_plot.pdf",sep=""))
  publish(hwrite(himg, br=TRUE), htmlRep1,name=paste("MA Plot ", comparison[2],"_vs_",comparison[3],sep=""))
  

  if(sum(!(target %in% signature))==0){
    prefix="Signature_Genes"
    target=signature;
    deGenes=head(resOrderedDF[resOrderedDF$symbol %in% target,]$symbol,100)
    myGenes=unique(c(deGenes))
  }else{
    prefix="Target_Genes"
    if(DE=="DE"){deGenes=head(resOrderedDF[resOrderedDF$symbol %in% target,]$symbol,50)}
    else{deGenes=as.character(head(target,DE))}
    myGenes=unique(c(deGenes))
    markMAplot(outdir,myGenes,comparison,prefix,res)  
    himg <- hwriteImage(paste(prefix,"_MA_plot.png",sep=""),link=paste(prefix,"_MA_plot.pdf",sep=""))
    publish(hwrite(himg, br=TRUE), htmlRep1,name=paste(prefix," MA Plot ", comparison[2],"_vs_",comparison[3],sep=""))
  }
  
  sampleGroup <- data.frame(colData(dds)[,c("grp")])
  rownames(sampleGroup)=rownames(colData(dds))
  colnames(sampleGroup)=c("group")
  if(!topGeneHeatmap(outdir,myGenes,comparison,prefix,sampleGroup,preprocessing,myMatrix)){
    himg <- hwriteImage(paste(prefix,"png",sep="."),link=paste(prefix,"pdf",sep="."))
    publish(paste(prefix," in Comparison ", comparison[2],"_vs_",comparison[3],sep=""), htmlRep1)
    publish(hwrite(himg, br=TRUE), htmlRep1,name=paste(prefix," in Comparison ", comparison[2],"_vs_",comparison[3],sep=""))
    resSubset <- res %>% data.frame %>% filter(symbol %in% target) %>% select(symbol,log2FoldChange,pvalue,padj) %>% arrange(symbol)
    fpkmSubset <- fpkm_table[resSubset$symbol, colnames(fpkm_table) %in% sampleSheet$sampleName]
    fpkmSubset <- fpkmSubset %>% data.frame %>% mutate(symbol=rownames(fpkmSubset)) 
    fpkmSubset <- fpkmSubset[,c(dim(fpkmSubset)[2],seq(1,dim(fpkmSubset)[2]-1))]
    mytable <- inner_join(resSubset,fpkmSubset,by="symbol")
    publish(paste(prefix," ", comparison[2],"_vs_",comparison[3],sep=""), htmlRep1)
    publish(mytable, htmlRep1)
  }else{
    publish(paste("No Significant Differnetial Expressed ",prefix," Between ",comparison[2],"_vs_",comparison[3],sep=""),htmlRep1)
  }

  
  deGenes=head(resOrderedDF$symbol,50)
  myGenes=unique(c(deGenes))
  sampleGroup <- data.frame(colData(dds)[,c("grp")])
  rownames(sampleGroup)=rownames(colData(dds))
  colnames(sampleGroup)=c("group")
  prefix="Top_50_LFC_DE_Genes"
  if(!topGeneHeatmap(outdir,myGenes,comparison,prefix,sampleGroup,preprocessing,myMatrix)){
    himg <- hwriteImage(paste(prefix,"png",sep="."),link=paste(prefix,"pdf",sep="."))
    publish(paste("Top 50 LFC DE Gene in Comparison ", comparison[2],"_vs_",comparison[3],sep=""), htmlRep1)
    publish(hwrite(himg, br=TRUE), htmlRep1,name=paste(prefix," in Comparison ", comparison[2],"_vs_",comparison[3],sep=""))
    resSubset <- resOrderedDF %>% data.frame %>% select(symbol,log2FoldChange,pvalue,padj) %>% arrange(symbol)
    fpkmSubset <- fpkm_table[resSubset$symbol, colnames(fpkm_table) %in% sampleSheet$sampleName]
    fpkmSubset <- fpkmSubset %>% data.frame %>% mutate(symbol=rownames(fpkmSubset)) 
    fpkmSubset <- fpkmSubset[,c(dim(fpkmSubset)[2],seq(1,dim(fpkmSubset)[2]-1))]
    mytable <- inner_join(resSubset,fpkmSubset,by="symbol")
    publish(paste("Differential Expressed Genes ", comparison[2],"_vs_",comparison[3],sep=""), htmlRep1)
    publish(hwrite("Download", link = paste("../",comparison[2],"_vs_",comparison[3],"_DE.txt",sep="")), htmlRep1)
    publish(mytable, htmlRep1)
  }else{
    publish(paste("No Significant Differnece Between ",comparison[2],"_vs_",comparison[3],sep=""),htmlRep1)
  }
  

  
  finish(htmlRep1)
  
  return(paste(comparison[2],"_vs_",comparison[3],sep=""))
}

# MA-plot. The log2 fold change for a particular comparison is plotted on the y-axis and the average of the counts normalized by size factor is shown on the x-axis (“M” for minus, because a log ratio is equal to log minus log, and “A” for average). Each gene is represented with a dot. Genes with an adjusted p value below a threshold (here 0.1, the default) are shown in red
markMAplot <- function(outdir,signature,comparison,prefix,res){
  if(outdir==""){outdir="."}
  png(filename=paste(outdir,"/",comparison[2],"_vs_",comparison[3],"/",prefix,"_MA_plot.png",sep=""))
  plotMA(res, ylim=c(-5,5), main=paste(prefix," ", comparison[2],"_vs_",comparison[3],sep=""))
  mylist=""
  for (i in signature){mylist=c(mylist,grep(paste("^",i,"$|:",i,"$",sep=""), res$symbol))}
  mylist=rownames(res)[as.numeric(mylist[-1])]
  for (topGene in mylist){
    topGeneName=res[topGene,]$symbol
    with(res[topGene, ], {
      points(baseMean, log2FoldChange, col="dodgerblue", cex=0.5, lwd=2)
      text(baseMean, log2FoldChange, topGeneName, pos=2, col="dodgerblue",cex=0.6)
    })
  }
  dev.off()
  pdf(file=paste(outdir,"/",comparison[2],"_vs_",comparison[3],"/",prefix,"_MA_plot.pdf",sep=""),onefile=FALSE)
  plotMA(res, ylim=c(-5,5), main=paste(prefix," ", comparison[2],"_vs_",comparison[3],sep=""))
  mylist=""
  for (i in signature){mylist=c(mylist,grep(paste("^",i,"$|:",i,"$",sep=""), res$symbol))}
  mylist=rownames(res)[as.numeric(mylist[-1])]
  for (topGene in mylist){
    topGeneName=res[topGene,]$symbol
    with(res[topGene, ], {
      points(baseMean, log2FoldChange, col="dodgerblue", cex=0.5, lwd=2)
      text(baseMean, log2FoldChange, topGeneName, pos=2, col="dodgerblue",cex=0.6)
    })
  }
  dev.off()
}

# Plot the Differential Expression Test Result
topGeneHeatmap <- function(outdir,myGenes,comparison,prefix,sampleGroup,preprocessing,myMatrix){
  if(outdir==""){outdir="."}
  if(preprocessing=="rld"){
    myMatrix=assay(rld)
  }else{
    temp=normalize.quantiles(myMatrix)
    myMatrix=log2(temp)
  }
  if(length(myGenes)<2){return(1)}
  mat <- myMatrix[ myGenes, rownames((sampleGroup))]
#  mat <- mat - rowMeans(mat)
  mywidth=0.5*dim(mat)[2]
  myheight=0.1*dim(mat)[1]
  if(mywidth<5){mywidth=5}
  if(myheight<5){myheight=5}
  png(filename=paste(outdir,"/",comparison[2],"_vs_",comparison[3],"/",prefix,".png",sep=""),units="in",width=mywidth,height=myheight,res=mywidth*myheight/0.25)
  pheatmap(mat, annotation=sampleGroup,col=colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize_row=6,fontsize_col=5,cellwidth = 5, cellheight = 5,cluster_cols = F,scale="row")
  dev.off()
  pdf(file=paste(outdir,"/",comparison[2],"_vs_",comparison[3],"/",prefix,".pdf",sep=""),onefile=FALSE,width=mywidth,height=myheight)
  pheatmap(mat, annotation=sampleGroup,col=colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize_row=6,fontsize_col=5,cellwidth = 5, cellheight = 5,cluster_cols = F,scale="row")
  dev.off()
  return(0)
}

