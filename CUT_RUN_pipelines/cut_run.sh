#!/bin/bash

# Create output directories
mkdir -p ${PROJ_PATH}/{fastq,alignment/{sam,bam,bed,bedgraph},peakCalling/SEACR}

for SAMPLE in ${SAMPLE_LIST}; do
# QC analysis
    fastqc -o ${PROJ_PATH}/fastq -f fastq ${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz

# Merge technical replicates/lanes if necessary
# for SAMPLE in ${SAMPLE_LIST}; do
#     cat ${FASTQ_DIR}/${SAMPLE}/*_R1_*.fastq.gz > ${PROJ_PATH}/fastq/${SAMPLE}_R1.fastq.gz
#     cat ${FASTQ_DIR}/${SAMPLE}/*_R2_*.fastq.gz > ${PROJ_PATH}/fastq/${SAMPLE}_R2.fastq.gz
# done

# Alignment
    bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${CORES} \
    -x ${BOWTIE2_INDEX} -1 ${PROJ_PATH}/fastq/${SAMPLE}_R1.fastq.gz -2 ${PROJ_PATH}/fastq/${SAMPLE}_R2.fastq.gz \
    | samtools sort -@ ${CORES} -O bam -o ${PROJ_PATH}/bam/${SAMPLE}.sorted.bam
    samtools index ${PROJ_PATH}/bam/${SAMPLE}.sorted.bam

# Deduplication
    picard MarkDuplicates I=${PROJ_PATH}/bam/${SAMPLE}.sorted.bam O=${PROJ_PATH}/bam/${SAMPLE}.rmdup.bam \
    M=${PROJ_PATH}/bam/${SAMPLE}.rmdup.metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true

# Filtering and conversion to BEDPE
    samtools view -b -F 0x04 -q 20 ${PROJ_PATH}/bam/${SAMPLE}.rmdup.bam \
    | bedtools bamtobed -bedpe -i - \
    | awk '$1==$4 && $6-$2 < 1000 {print $0}' > ${PROJ_PATH}/bed/${SAMPLE}.clean.bed
    cut -f 1,2,6 ${PROJ_PATH}/bed/${SAMPLE}.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${PROJ_PATH}/bed/${SAMPLE}.fragments.bed

# BedGraph generation with genome coverage
    bedtools genomecov -bg -i ${PROJ_PATH}/bed/${SAMPLE}.fragments.bed -g ${CHROM_SIZE_FILE} > ${PROJ_PATH}/bedgraph/${SAMPLE}.fragments.bedgraph

# Peak calling with SEACR
    bash ${SEACR_SCRIPT} ${PROJ_PATH}/bedgraph/${SAMPLE}.fragments.bedgraph \
    0.01 non stringent ${PROJ_PATH}/peakCalling/SEACR/${SAMPLE}_seacr_top0.01.peaks
done
