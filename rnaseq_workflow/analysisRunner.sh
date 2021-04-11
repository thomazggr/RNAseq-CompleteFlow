# Downloads human genome annotation to the exact folder
echo "Downloading .gtf annotation file, please wait..."
#wget -P annotation/ ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
echo "Finished downloading annotation!"

# Downloads indexed genome from Hisat2 AWS
echo "Downloading indexed human genome from Hisat2, please wait..."
#wget -P index/ https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
echo "Genome indexes has been downloaded!"

# Set your samples here! Set the name before file format .fastq, e.g.: sample_1.fastq, so it will be SAMPLES = "sample_1"
# If you have multiple samples, set with only 1 (one) space between names, as the next example: SAMPLES = "sample_1 sample_2 sample_3"
SAMPLES="case3"

for SAMPLE in $SAMPLES; do
    # Starting fastqc
    echo "Starting FastQC for quality control, please wait..."
    fastqc -o results/1_initial_qc/ --noextract input/${SAMPLE}.fastq
    echo "FastQC has ended, you can check quality control at results/1_initial_qc/"

    # Starts trimmage using TrimGalore!
    echo "Starting file trimmage using TrimGalore!, please wait..."
    trim_galore --quality 20 --fastqc --length 25 --output_dir results/2_trimmed_output/ input/${SAMPLE}.fastq
    echo "Trimmage has ended!"

    # Gunzip trimmed fq file to Hisat
    echo "Gunzipping your trimmed fq file, please wait..."
    gzip results/2_trimmed_output/${SAMPLE}_trimmed.fq
    echo "You've gunzipped your file, moving on >.>"

    # Starts alignment using hisat2 for low memory comsumption ~8gb for human genome an medium sized fastq file (~1.5gb)
    echo "Starting alignment using Hisat2, this can take a while, go grab a coffee!"
    echo "And also, please wait..."
    hisat2 -p 8 -x index/grch38_genome/grch38/genome -U results/2_trimmed_output/${SAMPLE}_trimmed.fq.gz -S results/4_aligned_sequences/${SAMPLE}.sam
    echo "Alignment has ended! HOORAY!"
done

# Go to dir where SAM files are
cd results/4_aligned_sequences

# Store list of SAM files as a variable
SAMLIST=$(ls -t ./*.sam | tr '\n' ' ')
echo $SAMLIST

# Final counts to be used in DE analysis
echo "Starting featureCounts for final counts, please wait..."
featureCounts -a ../../annotation/* -o ../../results/5_final_counts/final_counts.txt -g 'gene_name' -T 4 $SAMLIST
echo "featureCounts has ended, counts will be used to DE analysis!"

# Back to the main folder
cd ../../

# Run last analysis using MultiQC
echo "MultiQC has started, please wait..."
multiqc results --outdir results/6_multiQC
echo "MultiQC has ended!"
echo "All your results are stored in the results folder!"