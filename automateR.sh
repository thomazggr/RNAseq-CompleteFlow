#!/bin/bash

sudo apt update -qq

sudo apt install --no-install-recommends software-properties-common dirmngr

sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

sudo apt install --no-install-recommends r-base r-base-dev libxml2-dev libcurl4-openssl-dev libssl-dev

echo "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) install.packages(\"BiocManager\")" | sudo R --no-save
echo "BiocManager::install(version = \"3.13\")" | sudo R --no-save

echo "install.packages(c(\"httr\", \"curl\", \"RCurl\", \"openssl\", \"XML\"))" | sudo R --no-save
echo "BiocManager::install(c( \"DESeq2\", \"ggplot2\", \"clusterProfiler\", \"AnnotationDbi\", \"ReactomePA\", \"gage\", \"ggsci\", \"dplyr\",\
\"DOSE\", \"org.Hs.eg.db\", \"org.Mm.eg.db\", \"pheatmap\", \"RColorBrewer\"))" | sudo R --no-save

