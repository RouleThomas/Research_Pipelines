#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



multiBigwigSummary bins -b \
    output/bigwig/PSC_Norm_Rep1.bw \
    output/bigwig/PSC_Norm_Rep2.bw \
    output/bigwig/PSC_Norm_Rep3.bw \
    output/bigwig/PSC_Norm_Rep4.bw \
    output/bigwig/PSC_Hypo_Rep1.bw \
    output/bigwig/PSC_Hypo_Rep2.bw \
    output/bigwig/PSC_Hypo_Rep3.bw \
    output/bigwig/PSC_Hypo_Rep4.bw \
 -o output/bigwig/multiBigwigSummary_TPM_PSC.npz


