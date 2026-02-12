#!/bin/bash
#SBATCH --mem=200G
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
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done





