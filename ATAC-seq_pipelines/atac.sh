# Step 1: Prefetch SRA files
cd $path_to_sra
cat srr_list.txt | while read line; do
    prefetch -X 60GB $line
done

# Step 2: Align reads with Bowtie2 and process with Samtools
cd $path_to_clean_data
cat sample_ids.txt | while read id; do
    bowtie2 --threads 16 -X 2000 --very-sensitive -x $path_to_reference/$genome_index -1 ${id}_1_val_1.fq.gz -2 ${id}_2_val_2.fq.gz | \
    samtools view -bS - > ${id}.bam
    samtools sort -o $path_to_bam/${id}.sorted.bam ${id}.bam
    cd $path_to_bam
    samtools index ${id}.sorted.bam
done
cp sample_ids.txt $path_to_bam/

# Step 3: Sort BAM files by read name
cd $path_to_bam
cat sample_ids.txt | while read id; do
    samtools sort -n -o $path_to_sorted/${id}.sn.bam ${id}.sorted.bam
    samtools index $path_to_sorted/${id}.sn.bam
done

# Step 4: Remove mitochondrial reads, filter, and sort
cd $path_to_bam
cat sample_ids.txt | while read id; do
    samtools view -h ${id}.sorted.bam | grep -v chrM | samtools view -bh - > ${id}.noMT.bam
    samtools sort ${id}.noMT.bam -o ${id}.noMT.sort
    samtools index ${id}.noMT.sort
    samtools view -bh -f 3 -q 10 ${id}.noMT.sort > ${id}.filt.noMT.bam
    samtools sort ${id}.filt.noMT.bam -o ${id}.filt.noMT.sort
    samtools index ${id}.filt.noMT.sort
done
cp sample_ids.txt $path_to_rmpcr/

# Step 5: PCR duplicate removal and conversion to BED
cd $path_to_rmpcr
cat sample_ids.txt | while read id; do
    picard MarkDuplicates I=$path_to_bam/${id}.filt.noMT.sort O=${id}.pcr.bam M=${id}.pcr.metrics.txt REMOVE_DUPLICATES=true
    samtools sort ${id}.pcr.bam -o ${id}.pcr.sort
    samtools index ${id}.pcr.sort
    samtools sort -n ${id}.pcr.sort -o ${id}.pcr.sortname
    samtools fixmate ${id}.pcr.sortname ${id}.pcr.fixed
    samtools view -bf 0x2 ${id}.pcr.fixed | bedtools bamtobed -i stdin -bedpe > ${id}.fixed.bed
done

# Step 6: Tn5 shift and minimal BEDPE conversion
cat sample_ids.txt | while read id; do
    sh $path_to_scripts/bedpeTn5shift.sh ${id}.fixed.bed > ${id}.tn5.bedpe
    sh $path_to_scripts/bedpeMinimalConvert.sh ${id}.tn5.bedpe > ${id}.minimal.bedpe
done

# Step 7: Call peaks with MACS2
cd $path_to_rmpcr
cat sample_ids.txt | while read id; do
    macs2 callpeak -t ${id}.minimal.bedpe -f BEDPE -g $genome_size --outdir $path_to_final -n final --nomodel --extsize 200 --shift 100 --qval 5e-2 -B --SPMR --call-summits --nolambda --keep-dup all
done

