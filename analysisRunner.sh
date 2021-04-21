# Creating usage function to guide someone
usage()
{
cat << EOF
    Command line usage: bash analysisRunner.sh -D -O hsa -T n_threads -S sample_1,sample_2,sample_3,sample_n

    -D |   --download   (Optional) If called, indexes and annotation files will be downloaded
    -O |   --organism   (Required) hsa for homo sapies or mmu for mus musculus samples
    -T |   --threads    (Optional) If no number of threads to be used is passed, 4 is the default number
    -S |   --samples    (Required) Name of .fastq files to be used in the analysis. Each file has to be separated by comma ,

EOF
}

# Declaring optional variables
THREADS=4
DOWNLOAD="N"

# While loop to parse all arguments
while [ "$1" != "" ]; do
    case $1 in
        -D | --download)
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
        * )
    esac
    shift
done

# Check if organism are empty (mandatory argument)
if [ -z $ORG ]; then
    echo "Organism is required. Please, provide either hsa or mmu with the -O or --organism flag in the command line when running analysisRunner.sh"
    exit
fi

# Check if samples are empty (mandatory argument)
if [ -z $SAMPLES ]; then
    echo "You need to pass your file name for the samples that will be used passing the -S or --samples flag when running analysisRunner.sh"
    exit
fi

# Check if download argument has been passed, YES argument will download annotation and index for the organism selected (either hsa or mmu)
# If download has not been passed it will use index and annotation files inside each folder
if [ "$DOWNLOAD" = "Y" ]; then
    if [ "$ORG" = "hsa" ]; then
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
        cd grch38_genome/grch38/genome
        # Genome indexes files to the index folder and remove left over folders from the download, changing the dir from the hisat2 command
        mv * ../../../
        cd ../../../
        rm -r grch38_genome
        cd ../
        echo "Genome indexes has been downloaded!"
        echo
    elif [ "$ORG" = "mmu" ]; then
        # Downloads mus musculus genome annotation to the exact folder
        echo "Downloading .gtf annotation file, please wait..."
        echo
        wget -P annotation/ ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/gencode.vM26.annotation.gtf.gz
        echo "Finished downloading annotation!"
        echo
        # Downloads indexed MMU genome from Hisat2 AWS
        echo "Downloading indexed mus musculus genome from Hisat2, please wait..."
        echo
        wget -P index/ https://cloud.biohpc.swmed.edu/index.php/s/grcm38/download
        cd grcm38/
        # Genome indexes files to the index folder and remove left over folders from the download, changing the dir from the hisat2 command
        mv * ../
        cd ../
        rm -r grcm38
        cd ../
        echo "Genome indexes has been downloaded!"
        echo
    else
        echo "The $ORG organism option does not exist right now. If this is a valid organism please make a formal request to be added in the workflow sending an email to thomaz@vivaldi.net"
    fi

elif [ "$DOWNLOAD" = "N" ]; then
    echo "No files will be download. Using annotation and index files that have been places in the correct folders."
    echo "Following with the analysis..."
    echo
else
    echo "Something unexpected happened with the download option. Please try it again!"
fi

# All samples that has passed as argument will be used here to be passed to the for loop
SAMPLESREP=`echo $SAMPLES | tr ',' ' '`

# Tests that can be run
#echo $SAMPLESREP
#echo $THREADS
#echo $DOWNLOAD
#echo $ORG
#exit

for SAMPLE in $SAMPLESREP; do
    # Starting fastqc
    echo "Starting FastQC for quality control, please wait..."
    echo
    fastqc -o results/1_initial_qc/ --noextract input/${SAMPLE}.fastq.gz
    echo "FastQC has ended, you can check quality control at results/1_initial_qc/"
    echo

    # Starts trimmage using TrimGalore!
    echo "Starting file trimmage using TrimGalore!, please wait..."
    echo
    trim_galore --quality 20 --fastqc --length 25 --output_dir results/2_trimmed_output/ input/${SAMPLE}.fastq.gz
    echo "Trimmage has ended!"
    echo

    # Starts alignment using hisat2 for low memory comsumption ~8gb for human genome an medium sized fastq file (~1.5gb)
    echo "Starting alignment using Hisat2, this can take a while, go grab a coffee!"
    echo
    echo "And also, please wait..."
    echo
    hisat2 -p $THREADS -x index -U results/2_trimmed_output/${SAMPLE}_trimmed.fq.gz -S results/4_aligned_sequences/${SAMPLE}.sam
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