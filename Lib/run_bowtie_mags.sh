#!/bin/bash

# Define arguments for ease of reading 
rgi_output=$1
reads1=$2
reads2=$3
sample=$4
outdir=$5

# Convert CSV from RGI to FASTA
awk -F '\t' 'NR>1 {gsub(/ /, "_", $9); print ">"$9"_"$3"-"$4"_"$5"\n"$18}' $rgi_output > $outdir/${sample}_arg_sequences.fasta

# Build Bowtie2 index
bowtie2-build $outdir/${sample}_arg_sequences.fasta $outdir/${sample}_arg_index
echo "finished building indexes"

# Align reads with Bowtie2
bowtie2 -x $outdir/${sample}_arg_index -1 $reads1 -2 $reads2 -S $outdir/${sample}_aligned.sam
echo "finished alignment"

# Convert to BAM files
samtools view -bS $outdir/${sample}_aligned.sam > $outdir/${sample}_aligned.bam
echo "finished converting to BAM"

# Sort BAM file
samtools sort $outdir/${sample}_aligned.bam -o $outdir/${sample}_aligned_sorted.bam
echo "finished sorting BAM"

# Index BAM file
samtools index $outdir/${sample}_aligned_sorted.bam
echo "finished indexing bam"

# Run idxstats
samtools idxstats $outdir/${sample}_aligned_sorted.bam > $outdir/${sample}_arg_counts.txt
echo "finished idxstats"

# Add headers to idxstats output
sed -i '1i ARG_Name\tSequence_Length\tMapped_Reads\tUnmapped_Reads' $outdir/${sample}_arg_counts.txt
echo "added headers to idxstats file"

# ./Lib/run_bowtie_mags.sh /dartfs/rc/lab/R/RossB/DartCF_infant_meta/Results/RGI/MAG/021021_01/021021_01.MAG.summary.txt /dartfs/rc/lab/R/RossB/DartCF_infant_meta/ReadData/021021_01/021021_01_R1.qc.humanDecontaminated.fastq.gz /dartfs/rc/lab/R/RossB/DartCF_infant_meta/ReadData/021021_01/021021_01_R2.qc.humanDecontaminated.fastq.gz 021021_01 /dartfs/rc/lab/R/RossB/DartCF_infant_meta/Results/RGI/MAG/021021_01/