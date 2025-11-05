# GN\_01

Tool to perform a complete RNA-seq bioinformatic analysis from the raw reads obtained in the sequencing process to the differential expression analysis of the samples. It includes all the different steps divided in different folders. 

***Authors:** Nerea Barrio Cabezas and Gabriel Coll Rieffel*

# Installation

You just need to download the GN\_01 complete compressed folder and decompress it using the following command:

tar \-xzvf [GN01.tar.gz](http://GN01.tar.gz)

In order to execute all the scripts successfully, we created a conda environment named gabner, containing all the necessary dependencies and programs. It is located in the CBGP cluster. So, in order to activate this environment, you should type: 

conda activate /data/2025/grado-biotech/gabriel.coll/miniforge3/envs/gabner

Alternatively, in case you cannot access this environment or you do not have access to the CBGP cluster, here we provide a detailed list of all the dependencies needed in the RNA-seq analysis, the versions with which we worked and how to install them in a conda environment the same way we did. There is also a link to the corresponding github repository of each of them in case you want to look for more information.

- **FastQC (v0.12.1)**: tool to perform quality control checks on raw sequence data from a high-throughput sequencing analysis. ([https://github.com/s-andrews/FastQC](https://github.com/s-andrews/FastQC?tab=readme-ov-file))  
    
  mamba install \-c bioconda fastqc  
    
- **MultiQC (v1.28)**: Joins the separate results of FastQC into a single report with the information of all the samples (it can also be used after the alignment process to summarise the results). ([https://github.com/MultiQC/MultiQC](https://github.com/MultiQC/MultiQC))

	mamba install \-c bioconda multiqc

- **Fastp (v0.24.1)**: Preprocessing tool for quality filtering, trimming, and adapter removal of FASTQ files. ([https://github.com/OpenGene/fastp](https://github.com/OpenGene/fastp))

mamba install \-c bioconda fastp	

- **Hisat2 (v2.2.1)**: Fast alignment program for mapping sequencing reads to a reference genome taking into account the splicing process (eukaryote samples). ([https://github.com/DaehwanKimLab/hisat2](https://github.com/DaehwanKimLab/hisat2))

	mamba install \-c bioconda hisat2

- **Samtools (v1.21):** tool to convert and sort the sam files to bam files (this way they can be used later on by featureCounts and are stored in a binary, smaller file). ([https://github.com/samtools/samtools](https://github.com/samtools/samtools))  
    
  mamba install \-c bioconda samtools  
    
- **FeatureCounts (v2.1.1)**: Tool for efficiently counting reads mapped to genomic features such as genes or exons. ([https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html](https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html))  
    
  mamba install \-c bioconda subread


  
For the differential expression analysis, instead of working in Bash shell, the script is designed to be applied through the R programming language, so that the different dependencies and packages are installed during the execution of the script. There is no need for external installation processes. However, here we supply the list of all the packages used to perform the analysis and create the different plots: 

- CRAN packages: ggpubr, ggplot2, dplyr, tidyverse, gplots, RColorBrewer, BiocManager  
    
  install.packages(\<package\_name\>, dependencies \= TRUE)  
    
- Bioconductor packages: "biomaRt", "clusterProfiler", "DESeq2", "org.Hs.eg.db", "AnnotationDbi", "GO.db", "amap"  
    
  BiocManager::install(\<package\_name\>, ask \= FALSE, update \= TRUE)  
  


The two most important ones, in charge of performing the actual analysis (the rest are more related to the plots’ generation and other additional tasks) are: 

- **DESeq2 (v1.48.1):** Statistical tool for identifying differentially expressed genes from RNA-seq count data. ([https://github.com/thelovelab/DESeq2](https://github.com/thelovelab/DESeq2))  
- **clusterProfiler (v4.16.0):** R package for statistical analysis and visualization of functional enrichment of differentially expressed genes. ([https://github.com/YuLab-SMU/clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler))

# Description

GN\_01 is a folder containing all the necessary scripts to perform a complete RNA-seq bioinformatic analysis from the raw reads obtained from the sequencing process to the differential expression analysis of the samples. The structure is organised in the following folders:   
📂GN\_01:

- 📄README.md  
- 📂scripts  
  - 📂00\_raw\_data  
  - 📂01\_pre\_fastqc  
  - 📂02\_trimming\_filtering  
  - 📂03\_post\_fastqc  
  - 📂04\_reads\_alignment  
  - 📂05\_counts\_alignment  
  - 📂06\_differential\_analysis

In each folder, you can find the corresponding script (“0X\_script.sh”) to execute that step of the pipeline. In the moment the script is executed (except in the case of the last one, as it follows a different pattern), two more folders are created automatically: 

- 📂0X\_results: contains all the output files generated by the script  
- 📂0X\_ logs: stores the standard output (stdout) and standard error (stderr) created during the process, in order to improve the traceability and the communication with the user

Where “0X” makes reference to the step of the pipeline we are at and the number of the script we are executing. 

Here we provided a detailed description of the purpose of each of them, the different tools it uses, with their different options and arguments, and the output it generates.

## 📂00\_raw\_data

**Objective:** This folder is the only one that has two different scripts. 

- 00\_symlink\_script.sh: create symbolic links to the compressed FASTQ files downloaded with the sra-toolkit by our classmates. You need to have previously downloaded the FASTQ files in another path. Those files contain the raw reads from the sequencing process.  
    
- 00\_sra\_script.sh: download RNA-Seq sequence data from NCBI using the sra-tools utility with a file containing a list of the accession numbers of the SRA files (e.g., SRR1234567). This script was never tested, so it may contain execution errors. 

**Tools:** only the 00\_sra\_script.sh uses external tools. 

- **sra-tools** ([https://github.com/ncbi/sra-tools](https://github.com/ncbi/sra-tools)). It can be installed with “mamba install \-c bioconda sra-tools”. The version is not provided as we never got to install it.   
    
  Sra-tools commands used in this script:   
- Download SRA files from NCBI in .sra format in the specified output directory


  prefetch \<accession\_number\> \-O \<output\_directory\>


- Convert SRA files to FASTQ format and save them in the specified output directory


  fasterq-dump \<sra\_file\> \--outdir \<output\_directory\>

**Output:** 

- 00\_symlink\_script.sh  
  - 00\_results: one .[fastq.gz](http://fastq.gz) file for each sample  
  - 00\_logs: one folder 00\_logs\_sample\_name for each sample and, inside it, one file logs.out and logs.err for that specific sample. There is also a 00\_logs\_general folder with one file 00\_logs.out and one file 00\_logs.err which contain the stdout and stderr of the first part of the script, with all the argument checks, before starting processing the first sample.

- 00\_sra\_script.sh (never tested)  
  - 00\_results: one .[fastq.gz](http://fastq.gz) file for each sample  
  - 00\_logs: one folder 00\_logs\_sample\_name for each sample and, inside it, one file logs.out and logs.err for that specific sample, and also a 00\_logs\_general\_folder for the first part of the script. 

## 📂01\_pre\_fastqc

**Objective:** perform a quality control of raw reads from RNA-Seq sequence data, from NCBI, using FastQC and MultiQC as tools.

**Tools and commands used:** 

- **FastQC.** Performs the quality control of the raw reads and generates an html file with the results:

  fastqc \-o \<output\_directory\> \<input\_file\>


- **MultiQC.** Generates a report combining all the results delivered by FastQC:  
    
  multiqc \-o \<output\_directory\> \-n \<report\_name\> \<input\_directory\>


  
**Output:** 

- 01\_results:  
  - fastqc\_results: one .html file with the report and a .zip file with the reads of each sample.  
  - multiqc\_results: one multiqc\_report.html and a multiqc\_report\_data folder with all the detailed information of the quality checks carried out.  
- 01\_logs:  
  - 01\_fastqc\_logs: one folder 01\_logs\_sample\_name for each sample and, inside it, one file 01\_logs\_sample\_name.out and 01\_logs\_sample\_name.err for that specific sample.  
  - 01\_logs\_general: one file 01\_logs.out and one file 01\_logs.err which contain the stdout and stderr of the first part of the script, with all the argument checks, before starting processing the first sample.  
  - 01\_multiqc\_logs: one file 01\_multiqc\_logs.out and one file 01\_multiqc\_logs.err with the stdout and stderr of multiqc respectively.

## 📂02\_trimming\_filtering

**Objective:** trim the fastq.gz files of the RNA-seq using Fastp. By default it will detect and eliminate automatically the adapters for paired-end reads, do not evaluate the duplication rate, cut the bases of the sliding window if the mean quality of that window is lower than 20, and discard the reads shorter than 30 bases. The user can then add more options to the trimming process and modify the default values of the parameters.

**Tools and flags:** 

- **Fastp.** The following options of fastp will be used:   
    
  - \-i: input file with the forward read (fastq.gz)  
  - \-I: input file with the reverse read (fastq.gz)  
  - \-o: output file with the trimmed forward read (fastq.gz)  
  - \-O: output file with the trimmed reverse read (fastq.gz)  
  - \--failed\_out: output file with the failed reads (fastq.gz)  
  - \--detect\_adapter\_for\_pe: detect automatically adapters for paired-end reads and eliminate them  
  - \--dont\_eval\_duplication: do not evaluate the duplication rate (as it is a RNA-seq, we normally expect many duplications that do not mean something went wrong, so we will save time and resources if we do not evaluate them)  
  - \--trim\_front1: trim the front of the forward read in N number of bases (default: 0\)  
  - \--trim\_front2: trim the front of the reverse read in N number of bases (default: 0\)  
  - \--trim\_tail1: trim the tail of the forward read in N number of bases (default: 0\)  
  - \--trim\_tail2: trim the tail of the reverse read in N number of bases (default: 0\)  
  - \--trim\_poly\_x: trim poly-X sequences (where X is A, C, G or T) (useful for polyA in mRNA and polyG that appears in Illumina sequencing)  
  - \--cut\_front: cut the front of the read (useful as the extremes usually have a lower quality)  
  - \--cut\_tail: cut the tail of the read (useful as the extremes usually have a lower quality)  
  - \--cut\_mean\_quality: cut the bases of the sliding window if the mean quality of that window is lower than the specified value (20 is the default value)  
  - \--cut\_window\_size: size of the sliding window (4 is the default value)  
  - \--length\_required: minimum length of the read to be kept (30 is the default value), if the read is shorter than this value, it will be discarded  
  - \-j: output file with the json report (by default it will be the name of the sample without the extension followed by " fastp\_report.json")  
  - \-h: output file with the html report (by default it will be the name of the sample without the extension followed by " fastp\_report.html")  
  - \--thread: number of threads to use (useful to save time)  
  - \--report\_title: title that will appear in the html report (by default, it will be the name of the sample without the extension followed by " fastp report")


  
**Output:** 

- 02\_results: two .fastq.gz files with the trimmed reads (one for the forward read and another one for the reverse read as we assume we are working with paired-end reads), .html report and .json report for each sample  
- 02\_logs: one folder 02\_logs\_sample\_name for each sample and, inside it, one file 02\_logs\_sample\_name.out and 02\_logs\_sample\_name.err for that specific sample. There is also a 02\_logs\_general folder with one file 02\_logs.out and one file 02\_logs.err which contain the stdout and stderr of the first part of the script, with all the argument checks, before starting processing the first sample.

## 📂03\_post\_fastqc

**Objective:** perform a quality control of trimmed and processed reads from RNA-Seq sequence data, using FastQC and MultiQC as tools.

**Tools and commands:** 

- **FastQC.** Command used in this script:

  fastqc \-o \<output\_directory\> \<input\_file\>


- **MultiQC.** Command used in this script:  
    
  multiqc \-o \<output\_directory\> \-n \<report\_name\> \<input\_directory\>

**Output:**

- 03\_results:  
  - fastqc\_results: one .html file with the report and a .zip file with the reads of each sample.  
  - multiqc\_results: one multiqc\_report.html and a multiqc\_report\_data folder with all the detailed information of the quality checks carried out.  
- 03\_logs:  
  - 03\_fastqc\_logs: one folder 03\_logs\_sample\_name for each sample and, inside it, one file 03\_logs\_sample\_name.out and 03\_logs\_sample\_name.err for that specific sample.  
  - 03\_logs\_general: one file 03\_logs.out and one file 03\_logs.err which contain the stdout and stderr of the first part of the script, with all the argument checks, before starting processing the first sample.  
  - 03\_multiqc\_logs: one file 03\_multiqc\_logs.out and one file 03\_multiqc\_logs.err with the stdout and stderr of MultiQC respectively.

## 📂04\_reads\_alignment

**Objective:** perform an alignment of trimmed and processed reads from RNA-Seq sequencing data, using the tool HISAT2, and a later quality control of the results with MultiQC as tool.

**Tools and commands:** 

- **Hisat2.** Hisat2 commands used in this script: 

- Create the genome index:

		hisat2-build \<reference\_genome\> \<genome\_index\_prefix\>

- Perform the alignment:

  		hisat2 \-x \<genome\_index\> \-1 \<reads\_1.fq\> \-2 \<read2\_2.fq\> \-S output.sam \--summary-file \<summary\_text\_file.txt\>


- **Samtools.** Converts and sorts the sam files to bam files before passing them to MultiQC. Samtools commands used in this script:  
    
- Convert the sam files to bam files:


  samtools view \-b \<input.sam\> \-o \<output.bam\>


- Sort the bam files by genome position:


  samtools sort \<input.bam\> \-o \<output.sorted.bam\>


- Create the index for the sorted bam files:


  samtools index \<input.sorted.bam\>


- Generate a text file with the statistics of the alignment of the sorted bam files:


  samtools flagstat \<input.sorted.bam\> \> \<output.txt\>


- Generate a text file with the statistics of the sorted bam files by chromosomes:


  samtools idxstats \<input.sorted.bam\> \> \<output.txt\>


  

- **MultiQC.** Generates a report combining all the summary results of the alignment delivered by Hisat2.  
    
  multiqc \-o \<output\_directory\> \-n \<report\_name\> \<input\_directory\>

**Output:**

- 04\_results:  
  - 04\_alignment\_results: a genome\_index folder with the index created by the aligner to carry out the spliced alignment and five files for each sample:  
    - sample\_name.bam: binary file containing the aligned sequencing reads.  
    - sample\_name.bam.bai: index file for the BAM.  
    - sample\_name\_flagstat.txt: summary statistics of the alignment process (e.g., total reads, mapped reads, properly paired reads).  
    - sample\_name\_idxstats.txt: number of reads mapped to each chromosome.  
    - sample\_name\_summary.txt: detailed alignment report generated by the aligner.  
  - 04\_multiqc\_results: one multiqc\_report.html and a multiqc\_report\_data folder with all the detailed information of the quality checks carried out.  
- 04\_logs:  
  - 04\_alignment\_logs: one folder 04\_logs\_sample\_name for each sample and, inside it, one file 04\_logs\_sample\_name.out and 04\_logs\_sample\_name.err for that specific sample.  
  - 04\_logs\_general: one file 04\_logs.out and one file 04\_logs.err which contain the stdout and stderr of the first part of the script, with all the argument checks, before starting processing the first sample.  
  - 04\_multiqc\_logs: one file 04\_multiqc\_logs.out and one file 04\_multiqc\_logs.err with the stdout and stderr of MultiQC respectively.

## 📂05\_counts\_alignment

**Objective:** perform a quantification of reads per gene from RNA-Seq sequencing data, which has been previously processed (trimmed and aligned), using featureCounts as tool and the annotation file of the genome of interest.

**Tools and flags:** 

- **FeatureCounts.** Command used in the script:  
    
  featureCounts \-O \-p \-T \<threads\> \-s \<strand\> \-o \<output\_file\> \-a \<annotation\_file\> \<input\_BAM\_files\>


  
**Output:**

- 05\_results: one file results.tsv and another file results.tsv.summary. The first one is a matrix with a list of all the genes and the counts of reads assigned to each gene, representing together the counts of all the samples divided in different columns (each sample’s counts in a different column). The other one is a summary of the total amount of reads that were assigned and unassigned in each sample.  
- 05\_logs:   
  - 05\_logs\_featurecounts: one file 05\_logs\_featurecounts.out and one file 05\_logs\_featurecounts.err with the stdout and stderr of featureCounts, respectively.  
  - 05\_logs\_general: one file 05\_logs.out and one file 05\_logs.err which contain the stdout and stderr of the first part of the script, with all the argument checks, before starting processing the first sample.

## 📂06\_differential\_analysis

**Objective:** perform a differential expression analysis from the counts matrix generated by featureCounts and a later Gene Ontology enrichment analysis of the DEGs main functions. 

**Tools:** 

- **R:** the only tool you will need is R programming language. The rest of the packages, such as DESeq2 for the differential expression analysis and clusterProfiler for the GO enrichment analysis, are provided within the code. In case they are not already installed in your computer, they will be automatically installed and loaded.

**Output:**

The output will be stored in one location or another depending on which output\_directory you define at the beginning. However, regardless of where they are downloaded, the following files and folders will be created during the execution of the R script provided (DEseq2\_biotech\_2025\_raw\_modified.Rmd): 

- DEseq2\_biotech\_2025\_raw\_modified.html: interactive HTML report with the results of the differential expression analysis.  
- DEseq2\_biotech\_2025\_raw\_modified.pdf: static PDF version of the differential expression analysis report.  
- Functional\_Analysis\_GO\_BP.csv: results of the Gene Ontology (GO) enrichment analysis for Biological Processes.  
- plots folder with:  
  - barplot\_GO.png: bar plot of the top enriched GO terms.  
  - Dendrogram.png: hierarchical clustering of samples based on gene expression.  
  - Heatmap.png: heatmap of the most differentially expressed genes.  
  - Histogram\_with\_filter.png: histogram of log-transformed mean expression (logSumCounts) with filtering applied.  
  - Histogram\_without\_filter.png: histogram of log-transformed mean expression (logSumCounts) without filtering.  
  - MA.png: MA plot showing log-fold changes versus mean expression.  
  - PCA\_plot.png: Principal Component Analysis plot for sample clustering.  
  - volcano\_plot.png: volcano plot highlighting significant differentially expressed genes.  
- ConditionHighGlucose\_vs\_LowGlucose folder with:  
  - ConditionHighGlucose\_vs\_LowGlucose\_final\_results.csv: complete table of differential expression results for the comparison between High Glucose and Low Glucose conditions.

# Running

There is a detailed explanation in each of the scripts of how to use them and the different options and arguments that should be provided. Of course, you should follow the logical order of execution of the different scripts, starting with 01\_script and finishing with 05\_script in bash, plus an R execution of the code in DEseq2\_biotech\_2025\_raw\_modified.Rmd.

Despite that, you can look inside the scripts for more information, here we provide a list of all the options and arguments you should apply and also an usage example to make it clearer. All the examples have been made assuming you are located in the corresponding folder in each step of the pipeline (e.g. you are in folder 04\_reads\_alignment when doing the alignment). 

### **00\_symlink\_script**

**Input parameters and arguments:**

 \- i: Input directory where the FASTQ files are located  
 \- o: Output directory (optional, default: current directory)  
 \- h: Help message (optional, displays the usage of the script)  
 \- v: Version of the script (optional)

**Usage example:**

 ./00\_symlink\_script.sh \-i \</path/to/input/directory\> \-o \</path/to/output/directory\>

### **00\_sra\_script (never tested)**

**Input parameters and arguments:**

 \- a: File with the accessions number of the SRA files (e.g., SRR1234567) in each line  
 \- o: Output directory (optional, default: current directory)  
 \- h: Help message (optional, displays the usage of the script)  
 \- v: Version of the script (optional)

**Usage example:** 

./00\_sra\_script.sh \-a \<file\_with\_accession\_numbers\> \-o \</path/to/output/directory\>

### **01\_script**

**Input parameters and arguments:**

 \- i: Directory (path) with RNA-seq raw reads, in fastq format or fastq.gz  
 \- o: Output directory (optional, default: current directory)  
 \- h: Help message (optional, displays the usage of the script)  
 \- v: Version of the script (optional)

**Usage example:** 

./01\_script.sh \-i \<directory\_with\_raw\_reads\> \-o \</path/to/output/directory\>

### **02\_script**

**Input parameters and arguments:**

 \- i: Input directory where the FASTQ files are located  
 \- o: Output directory (optional, default: current directory)  
 \- d: Use this option if you want to store the failed reads in a separate file and specify the path and name of that file (e.g., ./failed\_reads.fastq.gz)  
 \- r: Number of bases to cut from the 5' end of the read forward (optional, default 0\)  
 \- R: Number of bases to cut from the 5' end of the read reverse (optional, default 0\)  
 \- e: Number of bases to cut from the 3' end of the read forward (optional, default 0\)  
 \- E: Number of bases to cut from the 3' end of the read reverse (optional, default 0\)  
 \- x: Use this option if you want to trim poly-X sequences   
 \- f: Use this option if you want to trim from 5' end  
 \- t: Use this option if you want to trim from 3' end  
 \- q: Quality threshold for trimming (optional, default: 20\)  
 \- w: Size of the sliding window (optional, default: 4\)  
 \- l: Minimum length of the read to be kept (optional, default: 30\)  
 \- T: Number of threads to use (optional, default: 1\)  
 \- h: Help message (optional, displays the usage of the script)  
 \- v: Version of the script (optional)

**Usage example:** 

./02\_script.sh \-i \</path/to/input/directory\> \-o \</path/to/output/directory\> \-d \<name\_of\_file\_to\_store\_failed\_reads.fastq.gz\> \-r \<number\_of\_bases\_to\_cut\_in\_read\_forward\_5'\> \-R \<number\_of\_bases\_to\_cut\_in\_read\_reverse\_5'\> \-e \<number\_of\_bases\_to\_cut\_in\_read\_forward\_3'\> \-E \<number\_of\_bases\_to\_cut\_in\_read\_reverse\_3'\> \-x \-f \-t \-q \<quality\_threshold\> \-w \<size\_of\_sliding\_window\> \-l \<minimum\_length\_of\_read\> \-T \<number\_of\_threads\>

### **03\_script**

**Input parameters and arguments:**

 \- i: Directory (path) with RNA-seq reads, in fastq format or fastq.gz  
 \- o: Output directory (optional, default: current directory)  
 \- h: Help message (optional, displays the usage of the script)  
 \- v: Version of the script (optional)

**Usage example:** 

./03\_script.sh \-i \<directory\_with\_reads\> \-o \</path/to/output/directory\>

### **04\_script**

**Input parameters and arguments:**

 \- i: Directory (path) with trimmed RNA-seq reads, in fastq.gz format  
 \- o: Output directory (optional, default: current directory)  
 \- g: Reference genome file (FASTA format)  
 \- a: Annotation file (GTF format)  
 \- h: Help message (optional, displays the usage of the script)  
 \- v: Version of the script (optional)

**Usage example:** 

./04\_script.sh \-i \<directory\_with\_reads\> \-o \</path/to/output/directory\> \-g \<reference\_genome\> \-a \<annotation\_file\>

### **05\_script**

**Input parameters and arguments:**

 \- i: Directory (path) with processed RNA-seq data, in BAM or SAM format  
 \- o: Output directory (optional, default: current directory)  
 \- a: Reference genome annotation file (GTF or GFF3 format)  
 \- T: Followed by the number of threads to use (optional, default: 1\)  
 \- O: Assigns reads to all their overlapping meta-features  
 \- p: Paired-end reads (optional, default: single-end reads)  
 \- s: Strandness of the data, followed by 0 (unstranded reads), 1 (stranded reads) or 2 (reversely stranded reads) (optional, default: 0\)  
 \- h: Help message (optional, displays the usage of the script)  
 \- v: Version of the script (optional)

**Usage example:** 

./05\_script.sh \-i \<directory\_with\_processed\_RNAseq\_data\> \-o \</path/to/output/directory\> \-a \<path/to/annotation/file\>

### **06\_script**

In this case, you’ll just have to run the whole script in R or RStudio specifying the correct output\_directory you want to store the results in and the input\_directory where the results.tsv file generated by featureCounts was kept. Previously, you will also have to create a metadata.csv file containing the particular samples ID and conditions used in your specific study (our metadata.csv file is provided within the 06\_differential\_analysis folder). 

# Technical details and reproducibility 

To perform the analysis, the reason why we chose fastp, hisat2 and featureCounts was mainly because they are faster and easier to use for beginner than other tools such as Trimmomatic for the trimming process, STAR for the spliced alignment and HTSeq-count for the quantification of the reads assigned to each gene. 

FastQC and MultiQC were chosen as they are considered the standard quality control tools for RNA-seq, which means they are pretty optimised for the task and provide very informative plots that could be harder to interpret in other tools such as fastp.

All the samples were processed individually to keep track of the different logs generated by each one of them. This way, if something went wrong during the execution process, we could very easily found the exact point and sample that was causing the trouble. The only one in which we did not follow this rule was the last one (05\_counts\_alignment) and just for a very simple reason: it was easier for DESeq2 to have the table with all the sample’s counts together rather than having to create it by merging the tables of the individual samples. This way, we could execute both scripts faster and in a simpler way. 

The purpose of this whole analysis was to study the RNA-seq raw reads from the following paper: 

Panigrahi G, Candia J, Dorsey TH, et al. (2023). *Diabetes-associated breast cancer is molecularly distinct and shows a DNA damage repair deficiency.* **JCI Insight**, 8(23): e170105. Publicado el 8 de diciembre de 2023\.  
doi: 10.1172/jci.insight.170105   
[https://pubmed.ncbi.nlm.nih.gov/37906280/](https://pubmed.ncbi.nlm.nih.gov/37906280/)

   
We were working with human cells obtained from breast cancer and comparing two conditions: high glucose vs low glucose. So the whole analysis focuses on the final objective of finding which biological processes are up or downregulated in breast cancer cell when grown in high glucose conditions vs low glucose conditions. However, we did not follow the exact same bioinformatic procedure that they did.

The reference human genome and annotation files used were obtained from Ensembl database: [https://www.ensembl.org/Homo\_sapiens/Info/Index](https://www.ensembl.org/Homo_sapiens/Info/Index)

The raw data from the paper were originally published in NCBI, and that is where we downloaded them from, under the ID of GSE236419.

Our data were paired-end, which we had to take into consideration while writing the scripts, that is why most of them have a specific option to inform whether the reads are or not paired-end. 

Finally, we end this explanation by leaving here the exact commands we used in the exact order they were used. This way, if anyone wants to reproduce the analysis and compare his/her results with ours, it can be easily done. We do not include the execution for the 06 step of the analysis because of the reasons already mentioned above.

Between one and another we always did a “cd” to the next folder in the analysis. We are not going to present the execution of the last step as it was carried out through R and not through the Bash shell. All the procedure was carried out in the CBGP cluster, specifically inside the path: `/data/2025/grado_biotech/gabriel.coll/GN_01`

**Execution carried out by the authors:** 

`./00_symlink_script.sh -i /data/2025/grado_biotech/marta.taboada/SRV_01/00-raw_data -o /data/2025/grado_biotech/gabriel.coll/GN_01/scripts/00_raw_data`

`./01_script.sh -i /data/2025/grado_biotech/gabriel.coll/GN_01/scripts/00_raw_data/00_results`

`./02_script.sh -r 10 -R 10 -x -f -t -q 20 -w 4 -l 30 -T 6 -i /data/2025/grado_biotech/gabriel.coll/GN_01/scripts/00_raw_data/00_results`

`./03_script.sh -i /data/2025/grado_biotech/gabriel.coll/GN_01/scripts/02_trimming_filtering/02_results`

`./04_script.sh -i /data/2025/grado_biotech/gabriel.coll/GN_01/scripts/02_trimming_filtering/02_results/ -g /data/2025/grado_biotech/gabriel.coll/GN_01/GenRef/reference_genome.fa -a /data/2025/grado_biotech/gabriel.coll/GN_01/GenRef/annotation_file.gtf`

`./05_script.sh -p -i /data/2025/grado_biotech/gabriel.coll/GN_01/scripts/04_reads_alignment/04_results/04_alignment_results -a /data/2025/grado_biotech/gabriel.coll/GN_01/GenRef/annotation_file.gtf`

