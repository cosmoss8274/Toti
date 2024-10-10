#!/bin/bash

mkdir -p ${BAM_DIR}

# Trimming
cat ${ACC_LIST_FILE} | while read SAMPLE_ID
do
    # i.e. for single-end data
    trim_galore -q 20 ${FASTQ_DIR}/${SAMPLE_ID}.fastq.gz --gzip -o ${TRIMMED_DIR} --cores 20 --fastqc

# Alignment
    bwa mem -t 20 -M ${GENOME_FA} ${TRIMMED_DIR}/${SAMPLE_ID}_trimmed.fq.gz |
    samtools sort -@ 32 -O BAM -o ${BAM_DIR}/${SAMPLE_ID}.sorted.bam
    samtools index ${BAM_DIR}/${SAMPLE_ID}.sorted.bam

# Deduplication
    gatk MarkDuplicates -I ${BAM_DIR}/${SAMPLE_ID}.sorted.bam --ADD_PG_TAG_TO_READS false \
    --REMOVE_SEQUENCING_DUPLICATES true \
    --CREATE_INDEX true \
    -O ${BAM_DIR}/${SAMPLE_ID}.rmdup.bam \
    -M ${BAM_DIR}/${SAMPLE_ID}.rmdup.matrix.txt

# Peak calling
    macs2 callpeak -c ${BAM_DIR}/Input-${SAMPLE_ID}.rmdup.bam \
    -t ${BAM_DIR}/${SAMPLE_ID}.rmdup.bam \
    -q 0.05 -f BAM -g hs -n ${SAMPLE_ID} \
    --outdir ${PEAK_DIR}
done

