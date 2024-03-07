# Define the array of sample IDs to process
declare -a SAMPLE_LIST=("Sample1" "Sample2" "Sample3") # Replace with actual sample IDs

# Define paths
RAW_FASTQ_DIR="/path/to/raw_fastq"
TRIMMED_FASTQ_DIR="/path/to/trimmed_fastq"
GENOME_DIR="/home/ychai/rnaseq/index_dir"
STAR_OUTPUT_DIR="/path/to/star_output"
RSEM_REF_DIR="/path/to/rsem_reference"
RSEM_OUTPUT_DIR="/path/to/rsem_output"

# read trimming
for SAMPLE_ID in ${SAMPLE_LIST[@]}
do 
    # Trim adapters and quality filter the reads
    trim_galore --fastqc --retain_unpaired --paired \
                ${RAW_FASTQ_DIR}/${SAMPLE_ID}_1.fastq.gz \
                ${RAW_FASTQ_DIR}/${SAMPLE_ID}_2.fastq.gz \
                -o ${TRIMMED_FASTQ_DIR}
done

#alignment
for SAMPLE_ID in ${SAMPLE_LIST[@]}
do 
    # Align reads to the genome with STAR
    STAR --runMode alignReads --runThreadN 1 \
         --genomeDir ${GENOME_DIR} \
         --outFilterMultimapNmax 20 \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverLmax 0.1 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --outFilterScoreMinOverLread 0.33 \
         --outFilterMatchNminOverLread 0.33 \
         --readFilesIn ${TRIMMED_FASTQ_DIR}/${SAMPLE_ID}_1_val_1.fq.gz \
                       ${TRIMMED_FASTQ_DIR}/${SAMPLE_ID}_2_val_2.fq.gz \
         --readFilesCommand zcat \
         --outFileNamePrefix ${STAR_OUTPUT_DIR}/${SAMPLE_ID} \
         --alignSoftClipAtReferenceEnds Yes \
         --quantMode TranscriptomeSAM GeneCounts \
         --outSAMtype BAM Unsorted SortedByCoordinate \
         --outSAMunmapped Within KeepPairs \
         --chimSegmentMin 15 \
         --chimJunctionOverhangMin 15 \
         --chimOutType WithinBAM SoftClip \
         --outSAMattributes NH HI AS nM NM MD jM jI XS \
         --outSAMattrRGline ID:rg1 SM:sm1
done

#quantification with rsem
rsem-prepare-reference --gtf ${RSEM_REF_DIR}/gencode.v38.annotation.gtf \
                       ${RSEM_REF_DIR}/hg38.fa \
                       ${RSEM_REF_DIR}/reference_name

for SAMPLE_ID in ${SAMPLE_LIST[@]}
do 
    # Calculate expression levels with RSEM
    rsem-calculate-expression --num-threads 12 \
                              --fragment-length-max 1000 \
                              --no-bam-output \
                              --paired-end \
                              --estimate-rspd \
                              --forward-prob 0.5 \
                              --bam \
                              ${STAR_OUTPUT_DIR}/${SAMPLE_ID}Aligned.toTranscriptome.out.bam \
                              ${RSEM_REF_DIR}/reference_name \
                              ${RSEM_OUTPUT_DIR}/${SAMPLE_ID}
done

