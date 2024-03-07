#!/bin/bash

# Define variables for directory paths
FASTQ_DIR="path_to_raw_fastq_files"
TRIMMED_DIR="path_to_trimmed_fastq_files"
GENOME_DIR="path_to_genome_directory"
BAM_DIR="path_to_output_bam_files"
METHYL_DIR="path_to_methylation_output"

# Define variable for the accession list file
ACC_LIST_FILE="path_to_accession_list_file.txt"

# Quality control and adapter removal
# This step performs quality checking and trims adapters from the FASTQ files
cat ${ACC_LIST_FILE} | while read SAMPLE_ID
do
    trim_galore -q 20 --fastqc --gzip --rrbs ${FASTQ_DIR}/${SAMPLE_ID}.fastq.gz -o ${TRIMMED_DIR}/
done

# Mapping reads, for PBAT libraries, use the --pbat option. For other libraries, use default settings.
    bismark -p 32 --bam --bowtie2 --genome ${GENOME_DIR} ${TRIMMED_DIR}/${SAMPLE_ID}_trimmed.fq.gz
done

# Deduplication and sorting of BAM files
# RRBS libraries skip the deduplication step
    deduplicate_bismark --bam ${BAM_DIR}/${SAMPLE_ID}_trimmed_bismark_bt2.bam
    samtools sort -o ${BAM_DIR}/${SAMPLE_ID}.de.sorted.bam ${BAM_DIR}/${SAMPLE_ID}_trimmed_bismark_bt2.deduplicated.bam
    samtools index ${BAM_DIR}/${SAMPLE_ID}.de.sorted.bam
done

# Filtering out multiple mappings and low-quality reads
    samtools view -bF 4 -q 20 ${BAM_DIR}/${SAMPLE_ID}.de.sorted.bam > ${BAM_DIR}/${SAMPLE_ID}.filtered.bam
    samtools sort -o ${BAM_DIR}/${SAMPLE_ID}.filtered.sorted.bam ${BAM_DIR}/${SAMPLE_ID}.filtered.bam
    samtools index ${BAM_DIR}/${SAMPLE_ID}.filtered.sorted.bam

# quantifies methylation levels at CpG sites
    bismark_methylation_extractor --gzip --buffer_size 40G --comprehensive --bedgraph --cytosine_report --genome_folder ${GENOME_DIR} \
    --output ${METHYL_DIR} --parallel 20 --report ${BAM_DIR}/${SAMPLE_ID}.filtered.sorted.bam
    # Generation of a nucleotide coverage report
    bam2nuc --dir ${METHYL_DIR} --genome_folder ${GENOME_DIR} ${BAM_DIR}/${SAMPLE_ID}.filtered.sorted.bam
done

# Geerating Bismark reports for the processed samples
    bismark2report --dir ${METHYL_DIR} \
    --alignment_report ${BAM_DIR}/${SAMPLE_ID}_trimmed_bismark_bt2_SE_report.txt \
    --dedup_report ${BAM_DIR}/${SAMPLE_ID}.deduplication_report.txt \
    --splitting_report ${BAM_DIR}/${SAMPLE_ID}.splitting_report.txt \
    --mbias_report ${BAM_DIR}/${SAMPLE_ID}.M-bias.txt \
    --nucleotide_stats_report ${BAM_DIR}/${SAMPLE_ID}.nucleotide_stats.txt
done

