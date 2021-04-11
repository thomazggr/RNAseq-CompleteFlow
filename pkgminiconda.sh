# Append new channel to download conda packages
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

# Download all packages needed to run RNA-seq workflow (star not included)
conda install -c anaconda git wget --yes
conda install -c bioconda fastqc trim-galore sortmerna subread multiqc hisat2 --yes

echo "All packages have been installed from Bioconda"
echo "Packages included: fastqc, trim galore, sortmerna, featurecounts, multiqc and hisat2"
