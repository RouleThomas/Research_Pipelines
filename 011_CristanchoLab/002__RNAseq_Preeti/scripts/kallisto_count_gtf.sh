#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8


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
        
for x in "${x[@]}"; do
echo "Processing sample ${x}"
kallisto quant -i ../../Master/meta/kallisto/transcripts.idx \
    -o output/kallisto/${x}_quant \
    -b 100 output/fastp/${x}_1.fq.gz output/fastp/${x}_2.fq.gz \
    -g ../../Master/meta/gencode.v47.annotation.gtf -t 8 --rf-stranded --genomebam --chromosomes ../../Master/meta/GRCh38_chrom_sizes.tab
done







