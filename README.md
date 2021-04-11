# RNA-seq workflow
## Content
- Tools used for the workflow
- Shell scripts and how to use them
- Results storage

# Tools for the workflow
## Anaconda
Anaconda is used for package management and installation. A Shell Script has been created to download Miniconda which will be used in this project.<br/>
The Shell Script downloads the lastest version of Miniconda directly from Anaconda repositories and will need some interaction for the installation to continue, accepting changes coming from Miniconda. <br/>
Furthermore the Sheel Script will add new channels to download the packages from and all packages needed for the analysis will be downloaded.
## R language and packages
For the last part of the analysis the R language will be neded just as some packages for Differential Expression analysis and plot generation.<br/>
To install R in Ubuntu based distros you can follow steps directly from this [LINK](https://cran.r-project.org/bin/linux/ubuntu/) which is from the CRAN project.<br/>
After installing r-base, you will need to install some other libs to compile R packages and dependencies for the R packages, following these commands (add sudo if needed):<br/>
```
apt install r-base-dev
```
```
apt install r-cran-xml r-cran-xml2 libxml2-dev libcurl4-openssl-dev libssl-dev
```
<br/>

Bioconductor will also be needed to download almost every package used, so install BiocManager following this [LINK](https://www.bioconductor.org/install/).<br/>
After installing all dependencies and BiocManager, start installing the packages needed for the R script, running these commands in the R command line interface in the same order as below:
```
install.packages(c("httr", "curl", "RCurl", "openssl", "XML"))
```
```
BiocManager::install(c("DESeq2","ggplot2","clusterProfiler","biomaRt","ReactomePA","ggsci","gage","dplyr","topGO","DOSE","org.Hs.eg.db","org.Mm.eg.db","pheatmap","genefilter","GO.db","KEGG.db","RColorBrewer")
```
## Guideline for the workflow
After installing everything, it's time to get ready and run the analysis.<br/>
But first of all, you will need to change some parameters to fit the analysis to your own RNA-seq experiment.<br/><br/>
After cloning or downloading the repository from Github, you will be able to see 4 Shell Scripts and a folder called DE_analysis containing R scripts and a metadata.txt file.<br/><br/>
### Follow these steps:
- Enter the directory of the downloaded repository and run the first Shell Script called minicondains if you need to install Miniconda into your computer:<br/>
```
bash minicondains.sh
```
- Then, after following the steps to install Minconda, you will need to restart your Terminal and run the second Shell Script called pkgminiconda to install all packages needed:<br/>
```
bash pkgminiconda.sh
```
- After installing all packages, create all the directories needed for result storage, index, annotation and files from your own experiment. Run:<br/>
```
bash directoryBuilder.sh
```
- You will be able to see that a whole new directory have been created to supply the analysis.<br/><br/>
Open the rnaseq_workflow folder and move your sample fastq files to the input folder.<br/><br/>
Now open the analysisRunner file and change the SAMPLES="" adding the file names of your own samples separated by single space and save, e.g.: SAMPLES="sample1 sample2 sample3".<br/><br/>
After that, open the DE_analysis folder and edit the metadata.txt file to your experiment type. An example has been given to help and guide.<br/><br/>
The next step is to run the analysisRunner file typing:<br/>
```
bash analysisRunner.sh
```
- After running the file, it will iterate over every .fastq file given in the SAMPLES variable and result in a final_counts.txt file.<br/>
Then, with this count table and the metadata prepared earlier you will be able to run the R script from source (either on command line or RStudio).<br/>
```
Rscript DErunner.R
```
- The R script will return 4 tables with normalized gene counts and the results with annotated genes. Also a heatmap from the expression value, volcano plot, KEGG enrichment and GO enrichment bar plot.

## Order of workflow and tool used
- FastQC: quality control of sample sequencing
- Trim Galore!: apply adapter and quality trimming to FastQ files
- HISAT2: alignment for mapping NGS reads to genomes
- featureCounts: read summarization program for counting reads from RNA or genomic DNA sequencing
- MultiQC: aggregation of results into a single report