#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


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
featureCounts -p -C -O -s 2 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v47.annotation.gtf \
	-o output/featurecounts/${x}.txt output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam
done



