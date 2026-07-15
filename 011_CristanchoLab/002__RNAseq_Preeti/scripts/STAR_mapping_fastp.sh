#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6   # Number of CPU cores per task


nthreads=6

x=(
    "ReN_Norm_Rep1"
    "ReN_Norm_Rep2"
    "ReN_Norm_Rep3"
    "ReN_Norm_Rep4"
    "ReN_Hypo_Rep1"
    "ReN_Hypo_Rep2"
    "ReN_Hypo_Rep3"
    "ReN_Hypo_Rep4"
    "PSC_Norm_Rep1"
    "PSC_Norm_Rep2"
    "PSC_Norm_Rep3"
    "PSC_Norm_Rep4"
    "PSC_Hypo_Rep1"
    "PSC_Hypo_Rep2"
    "PSC_Hypo_Rep3"
    "PSC_Hypo_Rep4"
    )

module load STAR/2.7.3a-GCC-9.3.0
module load SAMtools/1.16.1-GCC-11.3.0

for x in "${x[@]}"; do
	STAR --genomeDir ../../Master/meta/STAR_hg38/ \
		--runThreadN 6 \
		--readFilesCommand zcat \
		--readFilesIn output/fastp/${x}_1.fq.gz output/fastp/${x}_2.fq.gz \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix output/STAR/fastp/${x}_
    samtools index output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam
done

