# Activate conda environment
conda activate $env_name

# Navigate to the clean data directory
cd $path_to_clean_data

# Trim adapters and perform quality control
trim_galore --fastqc --phred33 $path_to_raw_data/${sample_id}_1.fastq.gz -o $path_to_trimmed_data

# Align reads with STAR
STAR --runMode alignReads --runThreadN $num_threads \
     --genomeDir $path_to_genome_index \
     --readFilesIn $path_to_trimmed_data/${sample_id}_trimmed.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $path_to_alignment_output/${sample_id}/ \
     --quantMode TranscriptomeSAM GeneCounts \
     --outSAMtype BAM Unsorted SortedByCoordinate \
     --twopassMode Basic \
     --alignMatesGapMax 100000 \
     --alignIntronMax 100000 \
     --outReadsUnmapped None \
     --chimSegmentMin 12 \
     --chimJunctionOverhangMin 12 \
     --alignSJDBoverhangMin 10 \
     --chimSegmentReadGapMax 3 \
     --alignSJstitchMismatchNmax 5 -1 5 5 \
     --outFilterScoreMinOverLread 0.33 \
     --outSAMstrandField intronMotif \
     --chimOutJunctionFormat 1 \
     --chimOutType WithinBAM \
     --outSAMunmapped Within \
     --outSAMattributes NH HI AS nM NM MD jM jI XS

# Index the BAM file with samtools
samtools index $path_to_alignment_output/${sample_id}/${sample_id}gene_Aligned.sortedByCoord.out.bam

# Sort the BAM file by read name
java -jar $path_to_picard SortSam \
     INPUT=$path_to_alignment_output/${sample_id}/${sample_id}gene_Aligned.sortedByCoord.out.bam \
     OUTPUT=$path_to_alignment_output/${sample_id}/${sample_id}gene_Aligned.sortedByRName.out.bam \
     TMP_DIR=$path_to_tmp_dir \
     SORT_ORDER=queryname VALIDATION_STRINGENCY=SILENT

# Quantify features with featureCounts
featureCounts -T $num_threads -t exon -g gene_id \
              -a $path_to_annotation \
              -o $path_to_counts/${sample_id}read.count \
              $path_to_alignment_output/${sample_id}/${sample_id}gene_Aligned.sortedByRName.out.bam \
              1>featureCounts.log 2>&1

