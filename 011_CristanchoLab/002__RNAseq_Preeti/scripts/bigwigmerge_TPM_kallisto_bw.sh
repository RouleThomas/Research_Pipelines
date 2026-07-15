#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00




# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig_kallisto/ReN_Norm_median.bedGraph \
    median output/bigwig_kallisto/ReN_Norm_Rep1.bw \
    output/bigwig_kallisto/ReN_Norm_Rep2.bw \
    output/bigwig_kallisto/ReN_Norm_Rep3.bw \
    output/bigwig_kallisto/ReN_Norm_Rep4.bw 
wiggletools write_bg output/bigwig_kallisto/ReN_Hypo_median.bedGraph \
    median output/bigwig_kallisto/ReN_Hypo_Rep1.bw \
    output/bigwig_kallisto/ReN_Hypo_Rep2.bw \
    output/bigwig_kallisto/ReN_Hypo_Rep3.bw \
    output/bigwig_kallisto/ReN_Hypo_Rep4.bw 
wiggletools write_bg output/bigwig_kallisto/PSC_Norm_median.bedGraph \
    median output/bigwig_kallisto/PSC_Norm_Rep1.bw \
    output/bigwig_kallisto/PSC_Norm_Rep2.bw \
    output/bigwig_kallisto/PSC_Norm_Rep3.bw \
    output/bigwig_kallisto/PSC_Norm_Rep4.bw 
wiggletools write_bg output/bigwig_kallisto/PSC_Hypo_median.bedGraph \
    median output/bigwig_kallisto/PSC_Hypo_Rep1.bw \
    output/bigwig_kallisto/PSC_Hypo_Rep2.bw \
    output/bigwig_kallisto/PSC_Hypo_Rep3.bw \
    output/bigwig_kallisto/PSC_Hypo_Rep4.bw 



# Sort the bedgraph 
bedtools sort -i output/bigwig_kallisto/ReN_Norm_median.bedGraph > output/bigwig_kallisto/ReN_Norm_median.sorted.bedGraph
bedtools sort -i output/bigwig_kallisto/ReN_Hypo_median.bedGraph > output/bigwig_kallisto/ReN_Hypo_median.sorted.bedGraph
bedtools sort -i output/bigwig_kallisto/PSC_Norm_median.bedGraph > output/bigwig_kallisto/PSC_Norm_median.sorted.bedGraph
bedtools sort -i output/bigwig_kallisto/PSC_Hypo_median.bedGraph > output/bigwig_kallisto/PSC_Hypo_median.sorted.bedGraph


# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig_kallisto/ReN_Norm_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_kallisto/ReN_Norm_median.bw

bedGraphToBigWig output/bigwig_kallisto/ReN_Hypo_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_kallisto/ReN_Hypo_median.bw

bedGraphToBigWig output/bigwig_kallisto/PSC_Norm_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_kallisto/PSC_Norm_median.bw

bedGraphToBigWig output/bigwig_kallisto/PSC_Hypo_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_kallisto/PSC_Hypo_median.bw