# Creating usage function to guide someone
usage()
{
cat << EOF
    Command line usage: bash analysisRunner.sh --D --O hsa --T n_threads --S sample_1 sample_2 sample_3 sample_n

    -D |   --download   (Optional) If called, indexes and annotation files will be downloaded
    -O |   --organism   (Required) hsa for homo sapies or mmu for mus musculus samples
    -T |   --threads    (Optional) If no number of threads to be used is passed, 4 is the default number
    -S |   --samples    (Required) Name of .fastq files to be used in the analysis. Each file has to be separated by a single space

EOF
}

# Declaring optional variables
THREADS=4
DOWNLOAD="N"

while [ "$1" != "" ]; do
    case $1 in
        -D | --download)
            shift
            DOWNLOAD="Y"
        ;;
        -O | --organism)
            shift
            ORG=$1
        ;;
        -T | --threads)
            shift
            THREADS=$1
        ;;
        -S | --samples)
            shift
            SAMPLES="$1"
        ;;
        -h | --help) usage
            exit
        ;;
        * ) usage
            exit 1
    esac
    shift
done

if [ -z $ORG ]; then
    echo "Organism is required. Please, provide either hsa or mmu with the -O or --organism flag in the command line when running analysisRunner.sh"
    exit
fi

if [ -z $SAMPLES ]; then
    echo "You need to pass your file name for the samples that will be used passing the -S or --samples flag when running analysisRunner.sh"
    exit
fi

if [ "$DOWNLOAD" -eq "Y" ]; then
    if [ "$ORG" -eq "hsa" ]; then
        # Downloads human genome annotation to the exact folder
        echo "Downloading .gtf annotation file, please wait..."
        echo
        wget -P annotation/ ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
        echo "Finished downloading annotation!"
        echo
        # Downloads indexed genome from Hisat2 AWS
        echo "Downloading indexed human genome from Hisat2, please wait..."
        echo
        wget -P index/ https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
        #mv genome indexes files to the index folder and remove left over folders from the download, changing the dir from the hisat2 command
        echo "Genome indexes has been downloaded!"
        echo
    elif [ "$ORG" -eq "mmu" ]; then
        # Downloads mus musculus genome annotation to the exact folder
        echo "Downloading .gtf annotation file, please wait..."
        echo
        wget -P annotation/ 
        echo "Finished downloading annotation!"
        echo
        # Downloads indexed MMU genome from Hisat2 AWS
        echo "Downloading indexed mus musculus genome from Hisat2, please wait..."
        echo
        wget -P index/ 
        #mv genome indexes files to the index folder and remove left over folders from the download, changing the dir from the hisat2 command
        echo "Genome indexes has been downloaded!"
        echo
    else
        echo "The $ORG organism option does not exist right now. If this is a valid organism please make a formal request to be added in the workflow sending an email to thomaz@vivaldi.net"
    fi

elif [ "$DOWNLOAD" -eq "N" ]; then
    echo "No files will be download. Using annotation and index files that have been places in the correct folders."
    echo "Following with the analysis..."
    echo
else
    echo "Something unexpected happened with the download option. Please try it again!"
fi

# Set your samples here! Set the name before file format .fastq, e.g.: sample_1.fastq, so it will be SAMPLES = "sample_1"
# If you have multiple samples, set with only 1 (one) space between names, as the next example: SAMPLES = "sample_1 sample_2 sample_3"
SAMPLES="sample1 sample2 sample3"

for SAMPLE in $SAMPLES; do
    # Starting fastqc
    echo "Starting FastQC for quality control, please wait..."
    echo
    fastqc -o results/1_initial_qc/ --noextract input/${SAMPLE}.fastq
    echo "FastQC has ended, you can check quality control at results/1_initial_qc/"
    echo

    # Starts trimmage using TrimGalore!
    echo "Starting file trimmage using TrimGalore!, please wait..."
    echo
    trim_galore --quality 20 --fastqc --gzip --length 25 --output_dir results/2_trimmed_output/ input/${SAMPLE}.fastq
    echo "Trimmage has ended!"
    echo

    # Starts alignment using hisat2 for low memory comsumption ~8gb for human genome an medium sized fastq file (~1.5gb)
    echo "Starting alignment using Hisat2, this can take a while, go grab a coffee!"
    echo
    echo "And also, please wait..."
    echo
    hisat2 -p 8 -x index/grch38_genome/grch38/genome -U results/2_trimmed_output/${SAMPLE}_trimmed.fq.gz -S results/4_aligned_sequences/${SAMPLE}.sam
    echo "Alignment has ended! HOORAY!"
    echo
done

# Go to dir where SAM files are
cd results/4_aligned_sequences

# Store list of SAM files as a variable
SAMLIST=$(ls -t ./*.sam | tr '\n' ' ')
echo $SAMLIST

# Final counts to be used in DE analysis
echo "Starting featureCounts for final counts, please wait..."
echo
featureCounts -a ../../annotation/* -o ../../results/5_final_counts/final_counts.txt -g 'gene_name' -T 4 $SAMLIST
echo "featureCounts has ended, counts will be used to DE analysis!"
echo
# Back to the main folder
cd ../../

# Run last analysis using MultiQC
echo "MultiQC has started, please wait..."
echo
multiqc results --outdir results/6_multiQC
echo "MultiQC has ended!"
echo
echo "All your results are stored in the results folder!"
echo