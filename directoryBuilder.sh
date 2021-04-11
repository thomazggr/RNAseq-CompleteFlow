if [ -d rnaseq_workflow/ ]; then
    echo "Directory already exists =("
else
    mkdir -p rnaseq_workflow/{annotation,index,input,results/{1_initial_qc,2_trimmed_output,4_aligned_sequences,5_final_counts,6_multiQC}}
    echo "Directories have been created! You can now move your input data to the >input< folder and run the .sh for analysis =)"
fi

echo "Ending directoryBuilder, have fun and get some coffee!"
