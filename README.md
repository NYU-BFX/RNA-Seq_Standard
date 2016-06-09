# RNA-Seq Standard

Author: Yixiao Gong
Date: 06/08/2016
Version 1.0

RNA-Seq Standard Pipeline is designed to perform standard RNA-Seq analysis for Illimina TrueSeq sequencing data. It generates a standard RNA-Seq HTML report which include alignment quality assessment, sample distance evaluation, PCA plots, dispersion plots, MA plots, raw count/normalized count/fpkm table, pairwised differential expression analysis, top differential expressed genes heatmap, user-defined targets differential expression analysis/heatmap. 

Manual:

1. You first need to have pre-built STAR reference in your [referenceFiles](referenceFiles/) folder. If you do not have it, please download a "genome.fa" file from public sites. And you also need to have a GTF file in your [referenceFiles](referenceFiles/) folder. This file should be the same GTF file you used to built STAR reference.
   Please use the following command to generate STAR Index:    
   
	```
	code/generate_STAR-Index params referenceFiles/STAR_Reference
	```
	
2. You need to have a "chromInfo.txt" file in your [referenceFiles](referenceFiles/) folder in order to produce RNA-Seq signal tracks. 

3. This pipeline can process GTC generated RNA-Seq data in the "/ifs/data/sequence/results/" folder. Please see the [template](meta_data/20160224.txt). Sample name need to be defined in this file for each of your samples. You will need to generate a .txt file for each of your fastq files directories put them into [meta data](meta_data/) directory.

4. This pipeline can also perform automatic download from SRA and process the samples. Please see the [template](meta_data/sra_info.txt). Sample names need to be defined in this file for corresponding SRX number. 

5. In [group information file](meta_data/group_info.txt) file, you need to categorize your sample into different groups to perform differential expression analysis.

6. If you have a signature gene list (for example you are doing knockdown of SPOP, then the list should include SPOP and some other related genes) which you would like to check the expression vaule immediately you can define it in the [signature gene list](meta_data/signature.txt). If you have a list of defined targets (for example some ChIP-Seq experiment define a list of direct targets of your knockdown gene), you can define them in the [target gene list](meta_data/target.txt). Or you can keep the existing link to the [signature gene list](meta_data/signature.txt). Please notice the name of the gene in these two file must be consistant with the GTC you use to generate the STAR genome. 


How to run it:

1. Once everything has been set up, you can run the pipeline in two stages:

   a. Submit your jobs for each of your sample using the following command:
   
	```
	cut -f1 meta_data/group_info.txt | code/skipn 1 | xargs -n1 -I {} qsub -hard -pe threaded 8 -l tmp_free=190G -l tmp_token=24G -j Y -b Y -cwd -N {} -o {}.out code/By_Sample {}
	```
	
   b. After all the samples are processed, you can run the following command to summarize your results:   
   
	```
        code/Summarize params pipeline/summarize `head -1 meta_data/group_info.txt | cut -f2`
	```
	
   c. The summary html report will be generated in the program directory and you need to unzip and open the html file for the report. 
