#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00




# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig/ReN_Norm_median.bedGraph \
    median output/bigwig/ReN_Norm_Rep1.bw \
    output/bigwig/ReN_Norm_Rep2.bw \
    output/bigwig/ReN_Norm_Rep3.bw \
    output/bigwig/ReN_Norm_Rep4.bw 
wiggletools write_bg output/bigwig/ReN_Hypo_median.bedGraph \
    median output/bigwig/ReN_Hypo_Rep1.bw \
    output/bigwig/ReN_Hypo_Rep2.bw \
    output/bigwig/ReN_Hypo_Rep3.bw \
    output/bigwig/ReN_Hypo_Rep4.bw 
wiggletools write_bg output/bigwig/PSC_Norm_median.bedGraph \
    median output/bigwig/PSC_Norm_Rep1.bw \
    output/bigwig/PSC_Norm_Rep2.bw \
    output/bigwig/PSC_Norm_Rep3.bw \
    output/bigwig/PSC_Norm_Rep4.bw 
wiggletools write_bg output/bigwig/PSC_Hypo_median.bedGraph \
    median output/bigwig/PSC_Hypo_Rep1.bw \
    output/bigwig/PSC_Hypo_Rep2.bw \
    output/bigwig/PSC_Hypo_Rep3.bw \
    output/bigwig/PSC_Hypo_Rep4.bw 



# Sort the bedgraph 
bedtools sort -i output/bigwig/ReN_Norm_median.bedGraph > output/bigwig/ReN_Norm_median.sorted.bedGraph
bedtools sort -i output/bigwig/ReN_Hypo_median.bedGraph > output/bigwig/ReN_Hypo_median.sorted.bedGraph
bedtools sort -i output/bigwig/PSC_Norm_median.bedGraph > output/bigwig/PSC_Norm_median.sorted.bedGraph
bedtools sort -i output/bigwig/PSC_Hypo_median.bedGraph > output/bigwig/PSC_Hypo_median.sorted.bedGraph


# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig/ReN_Norm_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/ReN_Norm_median.bw

bedGraphToBigWig output/bigwig/ReN_Hypo_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/ReN_Hypo_median.bw

bedGraphToBigWig output/bigwig/PSC_Norm_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/PSC_Norm_median.bw

bedGraphToBigWig output/bigwig/PSC_Hypo_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/PSC_Hypo_median.bw