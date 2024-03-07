# Create reference metadata for Drop-seq
sh $path_to_dropseq_tools/create_Drop-seq_reference_metadata.sh \
   -n $genome_annotation_name \
   -r $path_to_genome_fasta \
   -g $path_to_annotation_gtf \
   -d $path_to_dropseq_tools \
   -o $output_dir \
   -a $path_to_STAR \
   -b $path_to_bgzip \
   -i $path_to_samtools

# Convert FASTQ to SAM format
java -jar $path_to_picard FastqToSam \
    FASTQ=$path_to_fastq_R1 \
    FASTQ2=$path_to_fastq_R2 \
    OUTPUT=$unaligned_bam_output \
    SAMPLE_NAME=$sample_name

# Align reads and perform Drop-seq analysis
$path_to_dropseq_tools/Drop-seq_alignment.sh \
   -g $path_to_STAR \
   -r $path_to_reference_fasta \
   -d $path_to_dropseq_tools \
   -o $output_dir \
   -t $tmp_dir \
   -s $path_to_STAR \
   $unaligned_read_pairs_bam

# Digital expression
$path_to_dropseq_tools/DigitalExpression \
   -INPUT $unaligned_read_pairs_bam \
   -OUTPUT $output_matrix \
   -SUMMARY $output_summary \
   -MIN_NUM_GENES_PER_CELL 200 \
   -TMP_DIR $tmp_dir

