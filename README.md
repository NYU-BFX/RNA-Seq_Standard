# RNA-Seq Standard

### Author: Yixiao Gong

RNA-Seq Standard Pipeline is designed to perform standard RNA-Seq analysis for Illumina TruSeq sequencing data. It generates a standard RNA-Seq HTML report which includes alignment quality assessment, sample distance evaluation, PCA plots, dispersion plots, MA plots, raw count/normalized count/fpkm table, pairwise differential expression analysis, top differential expressed genes heatmap, user-defined targets differential expression analysis/heatmap. 

#### Setting up:

1. You first need to have pre-built STAR reference in your [referenceFiles](referenceFiles/) folder. If you do not have it, please download a "genome.fa" file from public sites.
   Please use the following command to generate STAR Index:    
   
	```
	code/generate_STAR-Index params referenceFiles/STAR_Reference
	```
	
2. You need to have a "chromInfo.txt" file in your [referenceFiles](referenceFiles/) folder in order to produce RNA-Seq signal tracks. 

3. This pipeline can process GTC generated RNA-Seq data in the "/ifs/data/sequence/results/" folder. Please see the [template](meta_data/20160224.txt). Sample name need to be defined in this file for each of your samples. You will need to generate a .txt file for each of your fastq files directories put them into [meta data](meta_data/) directory.

4. This pipeline can also perform automatic download from SRA and process the samples. Please see the [template](meta_data/sra_info.txt). Sample names need to be defined in this file for corresponding SRX number. 

5. In [group information file](meta_data/group_info.txt) file, you need to categorize your sample into different groups to perform differential expression analysis.

6. This pipeline use STAR aligner to align the reads. An transcripts annotation GTF format file ("--sjdbGTFfile" option for STAR aligner in the [params](params) file) is required for STAR to extract splice junctions and use them to greatly improve the accuracy of the mapping. And this step also counting number of reads per gene for all the genes in the transcripts annotation GTF file using in the alignment step. The downstream RNA-Seq standard report are all based on this transcript annotation GTF. 

7. If you have a signature gene list (for example you are doing knockdown of SPOP, then the list should include SPOP and some other related genes) which you would like to check the expression vaule immediately you can define it in the [signature gene list](meta_data/signature.txt). If you have a list of defined targets (for example some ChIP-Seq experiment define a list of direct targets of your knockdown gene), you can define them in the [target gene list](meta_data/target.txt). Or you can keep the existing link to the [signature gene list](meta_data/signature.txt). Please notice the name of the gene in these two file must be consistant with the names in the GTF file you used as "--sjdbGTFfile" option when you are doing STAR alignment. 


#### How to run it:

Once everything has been set up, you can run the pipeline in two stages:

1. Submit your jobs for each of your sample using the following command:
   
	```
	cut -f1 meta_data/group_info.txt | code/skipn 1 | xargs -n1 -I {} qsub -hard -pe threaded 8 -l tmp_free=190G -l tmp_token=24G -j Y -b Y -cwd -N {} -o {}.out code/By_Sample {}
	```
	
2. After all the samples are processed, you can run the following command to summarize your results:   
   
	```
        head -1 meta_data/group_info.txt | cut -f2- | xargs -n1 code/Summarize params pipeline/summarize
	```
3. The summary html report will be generated in the program directory and you need to unzip and open the html file for the report. 
