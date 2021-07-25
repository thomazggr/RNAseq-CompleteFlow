# Drawing the hello entrance from the script
hello()
{
cat << EOF
O       o O       o O       o O       o O       o O   
| O   o | | O   o | | O   o | | O   o | | O   o | | O  
| | O | | | | O | | | | O | | | | O | | | | O | | | | O 
| o   O | | o   O | | o   O | | o   O | | o   O | | o 
o       O o       O o       O o       O o       O o   
EOF
}

# Creating usage function to guide someone
usage()
{
cat << EOF
    Command line usage: bash analysisRunner.sh -D -P -O hsa -T n_threads -S sample_1,sample_2,sample_3,sample_n

    -D |   --download   (Optional) If called, indexes and annotation files will be downloaded
    -O |   --organism   (Required) hsa for homo sapies or mmu for mus musculus samples
    -T |   --threads    (Optional) If no number of threads to be used is passed, 4 is the default number
    -S |   --samples    (Required) Name of .fastq files to be used in the analysis. Each file has to be separated by comma ,
    -P |   --paired     (Optional) Pass -P or --paired if you are running paired-ended analysis. Sample's file name has to be as sample1_R1 sample1_R2 but only sample1 will be passed.

EOF
}

# Declaring optional variables
THREADS=4
DOWNLOAD="N"
PE="N"

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
        -P | --paired)
            shift
            PE="Y"
        ;;
        -h | --help)
            usage
            exit
        ;;
        * )
    esac
    shift
done

# Check if organism are empty (mandatory argument)
if [ -z $ORG ]; then
    echo -e "Organism is required. Please, provide either hsa or mmu with the -O or --organism flag in the command line when running analysisRunner.sh"
    exit
fi

# Check if samples are empty (mandatory argument)
if [ -z $SAMPLES ]; then
    echo -e "You need to pass your file name for the samples that will be used passing the -S or --samples flag when running analysisRunner.sh"
    exit
fi

# All samples that has passed as argument will be used here to be passed to the for loop
SAMPLESREP=`echo -e $SAMPLES | tr ',' ' '`
if [ "$PE" == "Y" ]; then # Paired ended samples
    for SAMPLE in $SAMPLESREP; do
        if [ ! -f input/${SAMPLE}_R1.fastq.gz ] && [ ! -f input/${SAMPLE}_R2.fastq.gz ]
        then
            echo -e "The file ${SAMPLE}.fastq.gz does not exist or is not placed in the input folder. Please check the name and the folder!"
            exit
        fi
    done
elif [ "$PE" == "N" ]; then # Single ended samples
    for SAMPLE in $SAMPLESREP; do
        if [ ! -f input/${SAMPLE}.fastq.gz ]
        then
            echo -e "The file ${SAMPLE}.fastq.gz does not exist or is not placed in the input folder. Please check the name and the folder!"
            exit
        fi
    done
fi
# Give the welcome package
STARTDT="$(date +'%Y-%m-%d_%H-%M')"
hello |& tee -a "loggersAT${STARTDT}.log"
echo -e "\n # # # # # # # # # # # # # # # # # # # #" |& tee -a "loggersAT${STARTDT}.log"
echo -e "\n ANALYSIS PARAMETERS TO BE USED" |& tee -a "loggersAT${STARTDT}.log"
echo -e "\n # # # # # # # # # # # # # # # # # # # #" |& tee -a "loggersAT${STARTDT}.log"
echo -e "Hello! Welcome to the RNA-seq workflow" |& tee -a "loggersAT${STARTDT}.log"
# Giving parameters away
echo -e "Start datetime: $STARTDT" |& tee -a "loggersAT${STARTDT}.log"
echo -e "Threads used in HISAT2: $THREADS" |& tee -a "loggersAT${STARTDT}.log"
echo -e "Paired-ended analysis? $PE" |& tee -a "loggersAT${STARTDT}.log" 
echo -e "Chosen organism: $ORG" |& tee -a "loggersAT${STARTDT}.log"
echo -e "Will download indexes and annotation files? $DOWNLOAD" |& tee -a "loggersAT${STARTDT}.log"
echo -e "Samples used: $SAMPLESREP" |& tee -a "loggersAT${STARTDT}.log"

# Check if download argument has been passed, YES argument will download annotation and index for the organism selected (either hsa or mmu)
# If download has not been passed it will use index and annotation files inside each folder
if [ "$DOWNLOAD" = "Y" ]; then
    echo -e "\n # # # # # # # # # # # # # # # # # # # #" |& tee -a "loggersAT${STARTDT}.log"
    echo -e "\n DOWNLOADING BOTH ANNOTATION AND INDEX FOR ORGANISM" |& tee -a "loggersAT${STARTDT}.log"
    echo -e "\n # # # # # # # # # # # # # # # # # # # # \n" |& tee -a "loggersAT${STARTDT}.log"
    if [ "$ORG" = "hsa" ]; then
        # Downloads human genome annotation to the exact folder
        echo -e "Downloading .gtf annotation file, please wait... \n" |& tee -a "loggersAT${STARTDT}.log"
        wget -P annotation/ ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
        echo -e "Finished downloading annotation! \n" |& tee -a "loggersAT${STARTDT}.log"
        # Downloads indexed genome from Hisat2 AWS
        echo -e "Downloading indexed human genome from Hisat2, please wait... \n" |& tee -a "loggersAT${STARTDT}.log"
        wget -P index/ https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
        cd grch38_genome/grch38/
        # Genome indexes files to the index folder and remove left over folders from the download, changing the dir from the hisat2 command
        mv * ../../
        cd ../../
        rm -r grch38_genome
        cd ../
        echo -e "Genome indexes has been downloaded! \n" |& tee -a "loggersAT${STARTDT}.log"
    elif [ "$ORG" = "mmu" ]; then
        # Downloads mus musculus genome annotation to the exact folder
        echo -e "Downloading .gtf annotation file, please wait... \n" |& tee -a "loggersAT${STARTDT}.log"
        wget -P annotation/ ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/gencode.vM26.annotation.gtf.gz
        echo -e "Finished downloading annotation! \n" |& tee -a "loggersAT${STARTDT}.log"
        # Downloads indexed MMU genome from Hisat2 AWS
        echo -e "Downloading indexed mus musculus genome from Hisat2, please wait... \n" |& tee -a "loggersAT${STARTDT}.log"
        wget -P index/ https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
        tar -xzf index/mm10_genome.tar.gz
        rm index/mm10_genome.tar.gz
        cd mm10
        # Genome indexes files to the index folder and remove left over folders from the download, changing the dir from the hisat2 command
        mv * ../index
        cd ..
        rm -r mm10/
        echo -e "Genome indexes has been downloaded! \n" |& tee -a "loggersAT${STARTDT}.log"
    else
        echo -e "The $ORG organism option does not exist right now. If this is a valid organism please make a formal request to be added in the workflow sending an email to thomaz@vivaldi.net \n" |& tee -a "loggersAT${STARTDT}.log"
    fi

elif [ "$DOWNLOAD" = "N" ]; then
    echo -e "No files will be download. Using annotation and index files that have been places in the correct folders. \nFollowing with the analysis..." |& tee -a "loggersAT${STARTDT}.log"
else
    echo -e "Something unexpected happened with the download option. Please try it again! \n" |& tee -a "loggersAT${STARTDT}.log"
fi

echo -e "\n # # # # # # # # # # # # # # # # # # # #" |& tee -a "loggersAT${STARTDT}.log"
echo -e "\n STARTING QC, TRIMMAGE AND ALIGNMENT FOR EVERY FILE" |& tee -a "loggersAT${STARTDT}.log"
echo -e "\n # # # # # # # # # # # # # # # # # # # # \n" |& tee -a "loggersAT${STARTDT}.log"
if [ "$PE" == "Y" ]; then # Paired ended analysis
    for SAMPLE in $SAMPLESREP; do
        echo -e "\n----------------------------------- FILE ${SAMPLE} -----------------------------------" |& tee -a "loggersAT${STARTDT}.log"
        # Starting fastqc
        echo -e "Starting FastQC for quality control, please wait... \n" |& tee -a "loggersAT${STARTDT}.log"
        fastqc -o results/1_initial_qc/ --noextract input/${SAMPLE}_R1.fastq.gz input/${SAMPLE}_R2.fastq.gz |& tee -a "loggersAT${STARTDT}.log"
        echo -e "FastQC has ended, you can check quality control at results/1_initial_qc/ \n" |& tee -a "loggersAT${STARTDT}.log"

        # Starts trimmage using TrimGalore!
        echo -e "Starting file trimmage using TrimGalore!, please wait... \n" |& tee -a "loggersAT${STARTDT}.log"
        trim_galore --quality 20 --fastqc --length 25 --paired --output_dir results/2_trimmed_output/ \
        input/${SAMPLE}_R1.fastq.gz input/${SAMPLE}_R2.fastq.gz |& tee -a "loggersAT${STARTDT}.log"
        echo -e "Trimmage has ended! \n" |& tee -a "loggersAT${STARTDT}.log"

        # Starts alignment using hisat2 for low memory comsumption ~8gb for human genome an medium sized fastq file (~1.5gb)
        echo -e "Starting alignment using Hisat2, this can take a while, go grab a coffee! \nAnd also, please wait..." |& tee -a "loggersAT${STARTDT}.log"
        hisat2 -p $THREADS -x index/genome -1 results/2_trimmed_output/${SAMPLE}_R1_val_1.fq.gz -2 results/2_trimmed_output/${SAMPLE}_R2_val_2.fq.gz \
        -S results/4_aligned_sequences/${SAMPLE}.sam |& tee -a "loggersAT${STARTDT}.log"
        echo -e "Alignment has ended! HOORAY! \n" |& tee -a "loggersAT${STARTDT}.log"
    done
elif [ "$PE" == "N" ]; then # Single ended analysis
    for SAMPLE in $SAMPLESREP; do
        echo -e "\n----------------------------------- FILE ${SAMPLE} -----------------------------------" |& tee -a "loggersAT${STARTDT}.log"
        # Starting fastqc
        echo -e "Starting FastQC for quality control, please wait... \n" |& tee -a "loggersAT${STARTDT}.log"
        fastqc -o results/1_initial_qc/ --noextract input/${SAMPLE}.fastq.gz |& tee -a "loggersAT${STARTDT}.log"
        echo -e "FastQC has ended, you can check quality control at results/1_initial_qc/ \n" |& tee -a "loggersAT${STARTDT}.log"

        # Starts trimmage using TrimGalore!
        echo -e "Starting file trimmage using TrimGalore!, please wait... \n" |& tee -a "loggersAT${STARTDT}.log"
        trim_galore --quality 20 --fastqc --length 25 --output_dir results/2_trimmed_output/ input/${SAMPLE}.fastq.gz |& tee -a "loggersAT${STARTDT}.log"
        echo -e "Trimmage has ended! \n" |& tee -a "loggersAT${STARTDT}.log"

        # Starts alignment using hisat2 for low memory comsumption ~8gb for human genome an medium sized fastq file (~1.5gb)
        echo -e "Starting alignment using Hisat2, this can take a while, go grab a coffee! \nAnd also, please wait..." |& tee -a "loggersAT${STARTDT}.log"
        hisat2 -p $THREADS -x index/genome -U results/2_trimmed_output/${SAMPLE}_trimmed.fq.gz \
        -S results/4_aligned_sequences/${SAMPLE}.sam |& tee -a "loggersAT${STARTDT}.log"
        echo -e "Alignment has ended! HOORAY! \n" |& tee -a "loggersAT${STARTDT}.log"
    done
fi


# Go to dir where SAM files are
cd results/4_aligned_sequences

# Store list of SAM files as a variable
SAMLIST=$(ls -t ./*.sam | tr '\n' ' ')
echo -e $SAMLIST |& tee -a "../../loggersAT${STARTDT}.log"

echo -e "\n # # # # # # # # # # # # # # # # # # # #" |& tee -a "loggersAT${STARTDT}.log"
echo -e "\n STARTING FEATURECOUNTS FOR ALL SAM FILES" |& tee -a "loggersAT${STARTDT}.log"
echo -e "\n # # # # # # # # # # # # # # # # # # # # \n" |& tee -a "loggersAT${STARTDT}.log"
# Final counts to be used in DE analysis
echo -e "Starting featureCounts for final counts, please wait... \n" |& tee -a "../../loggersAT${STARTDT}.log"
featureCounts -a ../../annotation/* -o ../../results/5_final_counts/final_counts.txt -g 'gene_name' -T 4 $SAMLIST |& tee -a "../../loggersAT${STARTDT}.log"
echo -e "featureCounts has ended, counts will be used to DE analysis! \n" |& tee -a "../../loggersAT${STARTDT}.log"
# Copy final counts file to the DE analysis folder
cp final_counts.txt ../../DE_analysis/
# Back to the main folder
cd ../../

echo -e "\n # # # # # # # # # # # # # # # # # # # #" |& tee -a "loggersAT${STARTDT}.log"
echo -e "\n STARTING MULTIQC SUMMARIZE RESULTS" |& tee -a "loggersAT${STARTDT}.log"
echo -e "\n # # # # # # # # # # # # # # # # # # # # \n" |& tee -a "loggersAT${STARTDT}.log"
# Run last analysis using MultiQC
echo -e "MultiQC has started, please wait... \n" |& tee -a "loggersAT${STARTDT}.log"
multiqc results --outdir results/6_multiQC |& tee -a "loggersAT${STARTDT}.log"
echo -e "MultiQC has ended! \nAll your results are stored in the results folder!" |& tee -a "loggersAT${STARTDT}.log"

ENDDT="$(date +'%Y-%m-%d_%H-%M')"
echo -e "Ending time: $ENDDT" |& tee -a "loggersAT${STARTDT}.log"