# Navigate to data download directory
cd $path_to_data_download

# Convert SRA to FASTQ and compress
fastq-dump --split-files --gzip $path_to_sra/$sra_id.sra -O $path_to_fastq_data

# Rename files to Cell Ranger compatible format
ls *_1.fastq.gz | cut -d"_" -f 1 > temp1
ls *_2.fastq.gz | cut -d"_" -f 1 > temp2
while read id; do
    mv ./${id}_1.fastq.gz ./${id}_S1_L001_R1_001.fastq.gz
done < temp1
while read id; do
    mv ./${id}_2.fastq.gz ./${id}_S1_L001_R2_001.fastq.gz
done < temp2
rm temp1 temp2

# Run Cell Ranger count
$path_to_cellranger count --id=$sample_id \
                          --fastqs=$path_to_fastq \
                          --sample=$sample_id \
                          --transcriptome=$path_to_transcriptome \
                          --localcores=$num_cores \
                          --localmem=$mem_gb
