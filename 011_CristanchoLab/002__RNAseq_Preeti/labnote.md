# Project

Two cell lines (ReNcell and hiPSCs) in two conditions (24h Normoxia vs 24h Hypoxia) in 4 Bio reps.


--> Directional mRNA library preparation (poly A enrichment), NovaSeq X Plus Series (PE150)


Key 	Sample name 	Belong experiment
A1	ReN_Nor_24h	Exp3 (Pooled)
A2	ReN_Nor_24h	Exp4 (Pooled)
A3	ReN_Nor_24h	Exp5 (Pooled)
A4	ReN_Nor_24h	Exp6 (Pooled)
A5	ReN_Hyp_24h	Exp3 (Pooled)
A6	ReN_Hyp_24h	Exp4 (Pooled)
A7	ReN_Hyp_24h	Exp5 (Pooled)
A8	ReN_Hyp_24h	Exp6 (Pooled)
B1	hiPSCs_Nor_24h 	Exp3 (Pooled)
B2	hiPSCs_Nor_24h 	Exp4 (Pooled)
B3	hiPSCs_Nor_24h 	Exp5 (Pooled)
B4	hiPSCs_Nor_24h 	Exp6 (Pooled)
B5	hiPSCs_Hyp_24h 	Exp3 (Pooled)
B6	hiPSCs_Hyp_24h 	Exp4 (Pooled)
B7	hiPSCs_Hyp_24h 	Exp5 (Pooled)
B8	hiPSCs_Hyp_24h 	Exp6 (Pooled)





# Data access

Access Cristancho lab folder from CHOP computer: `/mnt/isilon/cristancho_data/Preeti/nsc_rnaseq/downloaded_data/01.RawData` 
Sample name: `/mnt/isilon/cristancho_data/Preeti/nsc_rnaseq/nsc_rna_seq_sample details.xlsx` 






# Pipeline
- Copy data to my `/scr1` working env
- Rename files
- FastQC (fastqc)
- Trimming (fastp)
- Count with featureCounts
- DEGs with DESEQ2





# Cp / import data


```bash
# Data copied from Preeti folder to input_raw/

## cp all data to input folder
cp input_raw/*/*.fq.gz input/

#--> Manually rename all files 
```

--> All good, files succesfully renamed according to `/mnt/isilon/cristancho_data/Preeti/nsc_rnaseq/nsc_rna_seq_sample details.xlsx` with following nomenclature: `[PSC or ReN]_[Norm or Hypo]_[Rep1 - 4]_1.fq`


A1	ReN_Norm_Rep1   ReN_Nor_24h 
A2	ReN_Norm_Rep2   ReN_Nor_24h
A3	ReN_Norm_Rep3   ReN_Nor_24h
A4	ReN_Norm_Rep4   ReN_Nor_24h
A5	ReN_Hypo_Rep1   ReN_Hyp_24h
A6	ReN_Hypo_Rep2   ReN_Hyp_24h
A7	ReN_Hypo_Rep3   ReN_Hyp_24h
A8	ReN_Hypo_Rep4   ReN_Hyp_24h
B1	PSC_Norm_Rep1   hiPSCs_Nor_24h 
B2	PSC_Norm_Rep2   hiPSCs_Nor_24h 
B3	PSC_Norm_Rep3   hiPSCs_Nor_24h 
B4	PSC_Norm_Rep4   hiPSCs_Nor_24h 
B5	PSC_Hypo_Rep1   hiPSCs_Hyp_24h 
B6	PSC_Hypo_Rep2   hiPSCs_Hyp_24h 
B7	PSC_Hypo_Rep3   hiPSCs_Hyp_24h 
B8	PSC_Hypo_Rep4   hiPSCs_Hyp_24h 




# Fastp cleaning

```bash
sbatch scripts/fastp.sh # 60456895 ok
```

# Fastqc 

Let's run fastqc to check why so many reads assigned to no features in featurecounts...


```bash
# Raw fastq
sbatch scripts/fastqc_raw.sh # 60548567 xxx


# fastp-clean fastq
sbatch scripts/fastqc_fastp.sh # 60548570 xxx
```





# STAR mapping fastp trim

```bash
sbatch --dependency=afterany:60456895 scripts/STAR_mapping_fastp.sh # 60457138 ok


## Convert alignment to bigwig
conda activate deeptools
sbatch --dependency=afterany:60457138 scripts/STAR_TPM_bw.sh # 60457411 ok

## Calculate median
conda activate BedToBigwig
sbatch --dependency=afterany:60457411 scripts/bigwigmerge_STAR_TPM_bw.sh # 60458041 ok

```

--> All good





# deepTool bigwig QC - STAR mapping


## TPM bigwig


```bash
conda activate deeptools

###################################
# Include X chr - All samples #####################
###################################
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_STAR_TPM.sh # 60553000 ok

############################################
# Plot alls ###########
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_TPM.npz \
    --transpose \
    --ntop 0 \
    --labels ReN_Norm_Rep1 ReN_Norm_Rep2 ReN_Norm_Rep3 ReN_Norm_Rep4 ReN_Hypo_Rep1 ReN_Hypo_Rep2 ReN_Hypo_Rep3 ReN_Hypo_Rep4 PSC_Norm_Rep1 PSC_Norm_Rep2 PSC_Norm_Rep3 PSC_Norm_Rep4 PSC_Hypo_Rep1 PSC_Hypo_Rep2 PSC_Hypo_Rep3 PSC_Hypo_Rep4 \
    --colors blue blue blue blue red red red red blue blue blue blue red red red red \
    --markers 'o' 'o' 'o' 'o' 'o' 'o' 'o' 'o' 's' 's' 's' 's' 's' 's' 's' 's' \
    -o output/bigwig/multiBigwigSummary_TPM_plotPCA.pdf \
    --plotWidth 7

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_TPM.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels ReN_Norm_Rep1 ReN_Norm_Rep2 ReN_Norm_Rep3 ReN_Norm_Rep4 ReN_Hypo_Rep1 ReN_Hypo_Rep2 ReN_Hypo_Rep3 ReN_Hypo_Rep4 PSC_Norm_Rep1 PSC_Norm_Rep2 PSC_Norm_Rep3 PSC_Norm_Rep4 PSC_Hypo_Rep1 PSC_Hypo_Rep2 PSC_Hypo_Rep3 PSC_Hypo_Rep4 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_TPM_heatmap.pdf

#################################




###################################
# Include X chr - PSC samples #####################
###################################
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_STAR_TPM-PSC.sh # 60553025 ok

############################################
# Plot PSC ###########
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_TPM_PSC.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_Norm_Rep1 PSC_Norm_Rep2 PSC_Norm_Rep3 PSC_Norm_Rep4 PSC_Hypo_Rep1 PSC_Hypo_Rep2 PSC_Hypo_Rep3 PSC_Hypo_Rep4 \
    --colors blue blue blue blue red red red red \
    --markers 's' 'o' '>' 'x' 's' 'o' '>' 'x' \
    -o output/bigwig/multiBigwigSummary_TPM_PSC_plotPCA.pdf \
    --plotWidth 6 \
    --plotHeight 7

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_TPM_PSC.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_Norm_Rep1 PSC_Norm_Rep2 PSC_Norm_Rep3 PSC_Norm_Rep4 PSC_Hypo_Rep1 PSC_Hypo_Rep2 PSC_Hypo_Rep3 PSC_Hypo_Rep4 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_TPM_PSC_heatmap.pdf

#################################






###################################
# Include X chr - ReN samples #####################
###################################
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_STAR_TPM-ReN.sh # 60553031 ok
############################################
# Plot ReN ###########
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_ReN.npz \
    --transpose \
    --ntop 0 \
    --labels ReN_Norm_Rep1 ReN_Norm_Rep2 ReN_Norm_Rep3 ReN_Norm_Rep4 ReN_Hypo_Rep1 ReN_Hypo_Rep2 ReN_Hypo_Rep3 ReN_Hypo_Rep4 \
    --colors blue blue blue blue red red red red \
    --markers 's' 'o' '>' 'x' 's' 'o' '>' 'x' \
    -o output/bigwig/multiBigwigSummary_ReN_plotPCA.pdf \
    --plotWidth 6 \
    --plotHeight 7

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_ReN.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels ReN_Norm_Rep1 ReN_Norm_Rep2 ReN_Norm_Rep3 ReN_Norm_Rep4 ReN_Hypo_Rep1 ReN_Hypo_Rep2 ReN_Hypo_Rep3 ReN_Hypo_Rep4 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_ReN_heatmap.pdf

#################################



```

--> **PSC and ReN cluster well per conditions**; no bio rep batch effect (**PSC cluster a bit better than ReN**; would expect more DEGs in this cell type)






# Count with featureCounts


Data look stranded: To confirm with Preeti.

AMPD stringenT: featureCounts -p -C -O \
AMPD relax: featureCounts -p -C -O -M --fraction \


```bash
conda activate featurecounts

# slight test
## -s 2 for stranded
featureCounts -p -C -O -M --fraction -s 2 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v47.annotation.gtf \
	-o output/featurecounts/ReN_Norm_Rep1.txt output/STAR/fastp/ReN_Norm_Rep1_Aligned.sortedByCoord.out.bam
#--> 49% assigned

## -s 1 for stranded
featureCounts -p -C -O -M --fraction -s 1 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v47.annotation.gtf \
	-o output/featurecounts/ReN_Norm_Rep1.txt output/STAR/fastp/ReN_Norm_Rep1_Aligned.sortedByCoord.out.bam
#--> 16% assigned

## unstranded, not counting multimapped reads
featureCounts -p -C -O \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v47.annotation.gtf \
	-o output/featurecounts/ReN_Norm_Rep1.txt output/STAR/fastp/ReN_Norm_Rep1_Aligned.sortedByCoord.out.bam
#--> 50% assigned! mostly not assigned due to multimapping; seems to be stranded in the end

## unstranded, counting multimapped reads
featureCounts -p -C -O -M --fraction \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v47.annotation.gtf \
	-o output/featurecounts/ReN_Norm_Rep1.txt output/STAR/fastp/ReN_Norm_Rep1_Aligned.sortedByCoord.out.bam
#--> 60%! Still not great...; the rest fall into unassigned no features

## unstranded, counting multimapped reads, testing previous gene annotations
featureCounts -p -C -O -M --fraction \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI.gtf \
	-o output/featurecounts/ReN_Norm_Rep1.txt output/STAR/fastp/ReN_Norm_Rep1_Aligned.sortedByCoord.out.bam
#--> 60%! issue does not come from using the new gene annotation; the rest fall into unassigned no features

## count on gene (REMOVE -O), not counting multimapped reads 
featureCounts -p -C -s 2 -t gene -g gene_id \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v47.annotation.gtf \
	-o output/featurecounts/ReN_Norm_Rep1.txt output/STAR/fastp/ReN_Norm_Rep1_Aligned.sortedByCoord.out.bam
#--> 48%!

## count on gene (REMOVE -O), counting multimapped reads
featureCounts -p -C -M --fraction -s 2 -t gene -g gene_id \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v47.annotation.gtf \
	-o output/featurecounts/ReN_Norm_Rep1.txt output/STAR/fastp/ReN_Norm_Rep1_Aligned.sortedByCoord.out.bam
#--> 66%!


## -s 2 for stranded, and test basic gtf
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz
# gunzip gencode.v49.basic.annotation.gtf.gz

featureCounts -p -C -O -M --fraction -s 2 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v49.basic.annotation.gtf \
	-o output/featurecounts/ReN_Norm_Rep1.txt output/STAR/fastp/ReN_Norm_Rep1_Aligned.sortedByCoord.out.bam
#--> 49% assigned; same value so I can still use the comprehenssive GENCODE version



# all samples:
sbatch scripts/featurecounts_unstranded.sh # 60549123 ok
#--> 50-65% uniquely aligned reads - NO, data is stranded
sbatch scripts/featurecounts_multi_unstranded.sh # 60549134 ok
#--> 60-75% uniquely aligned reads - NO, data is stranded
sbatch scripts/featurecounts.sh # 60547513 ok
#--> 50-60% uniquely aligned reads - NO, many multimapped reads!

## The two options below are good:
sbatch scripts/featurecounts_multi.sh # 60547516 ok
#--> 55-65% uniquely aligned reads - YES: output/featurecounts_multi
sbatch scripts/featurecounts_multi_gene.sh # 61302605 ok
#--> 65-75% uniquely aligned reads - YES: output/featurecounts_multi_gene



```

test with `ReN_Norm_Rep1`
- ~50%, 16% alignment with stranded paramters `-s 2` and `-s 1`, respectively + looking at IGV bam, **files are stranded**
- Many mapping to *unassigned features*, and *multimapped reads*.. Let's generate two versions:
  - one counting multimapped reads (`*_multi`), and another one removing multimapped reads, improve a bit but not so much
- I noticed checking bigwig on IGV that many reads fall within introns, should explain the high *unassigned features*
  - For this let's do a version where I **count on gene, not exon**; **when counting on gene level do NOT use` -O`**; indeed this count reads in gene that overlap twice, but it will lead to artificial correlation and issue with DESEQ2 when counting on gene body...


--> Seems **count on exon, stranded, and count multimapped reads is the best approach**: `output/featurecounts_multi` (count on gene as also been generated at `output/featurecounts_multi_gene`)


--> Looking at IGV, some samples show reads within introns; seems occuring more in PSC than ReN; and seems not dependent on condition (happen for both Norm and Hypo on some replicates)



# Calculate TPM and RPKM


Use custom R script `RPKM_TPM_featurecounts.R` as follow:
```bash
conda activate deseq2
# output/featurecounts_multi
## Rscript scripts/RPKM_TPM_featurecounts.R INPUT OUTPUT_PREFIX
sbatch scripts/featurecounts_TPM.sh # 61302936 ok
## mv all output to output/tpm or rpkm folder
mv output/featurecounts_multi/*tpm* output/tpm_featurecounts_multi/
mv output/featurecounts_multi/*rpkm* output/rpkm_featurecounts_multi/


# output/featurecounts_multi_gene
## Rscript scripts/RPKM_TPM_featurecounts.R INPUT OUTPUT_PREFIX
sbatch scripts/featurecounts_multi_gene_TPM.sh # 61303337
## mv all output to output/tpm or rpkm folder
mv output/featurecounts_multi_gene/*tpm* output/tpm_featurecounts_multi_gene/
mv output/featurecounts_multi_gene/*rpkm* output/rpkm_featurecounts_multi_gene/
```

--> All good. 




**Display gene with TPM**:

```R
# library
library("DESeq2")
library("tidyverse")
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")
library("rtracklayer")

set.seed(42)



# import GTF for gene name
gtf <- import("../../Master/meta/gencode.v47.annotation.gtf")
## Extract geneId and geneSymbol
gene_table <- mcols(gtf) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct() %>%
  as_tibble()
## Rename columns
colnames(gene_table) <- c("Geneid", "geneSymbol")



# Plot with TPM
## import tpm
#### Generate TPM for ALL samples
#### collect all samples ID
samples <- c("PSC_Norm_Rep1", "PSC_Norm_Rep2" ,"PSC_Norm_Rep3" ,"PSC_Norm_Rep4" ,"PSC_Hypo_Rep1", "PSC_Hypo_Rep2", "PSC_Hypo_Rep3", "PSC_Hypo_Rep4", "ReN_Norm_Rep1", "ReN_Norm_Rep2" ,"ReN_Norm_Rep3" ,"ReN_Norm_Rep4" ,"ReN_Hypo_Rep1", "ReN_Hypo_Rep2", "ReN_Hypo_Rep3", "ReN_Hypo_Rep4")


## Make a loop for importing all tpm data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/tpm_featurecounts_multi/", sample, "_tpm.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    dplyr::select(Geneid, starts_with("output.STAR.")) %>%
    dplyr::rename(!!sample := starts_with("output.STAR."))
}

## Merge all dataframe into a single one and add gene names
tpm_all_sample <- purrr::reduce(sample_data, full_join, by = "Geneid") 

tpm_all_sample = tpm_all_sample %>%
  left_join(gene_table)



# write.table(  tpm_all_sample,  file = "output/tpm_featurecounts_multi/tpm_all_sample.tsv",  sep = "\t",  row.names = FALSE,  quote = FALSE)
### If need to import: read.table(  "output/tpm_featurecounts_multi/tpm_all_sample.tsv",  sep = "\t",  header = TRUE)

# plot some genes
tpm_all_sample_tidy <- tpm_all_sample %>%
  pivot_longer(
    cols = -c(Geneid, geneSymbol),
    names_to = "variable",
    values_to = "tpm"
  ) %>%
  separate(
    variable,
    into = c("tissue", "condition", "replicate"),
    sep = "_"
  ) %>%
  dplyr::rename(
    gene = Geneid
  )


# genes 
c("EFNB2", "EPHA4") # 


plot_data <- tpm_all_sample_tidy %>%
  unique() %>%
  filter(geneSymbol %in% c("EFNB2")) %>%
  group_by(geneSymbol, tissue, condition) %>%
  summarise(mean_log2tpm = mean(log2(tpm + 1)),
            se_log2tpm = sd(log2(tpm + 1)) / sqrt(n())) %>%
  ungroup()


plot_data$tissue <-
  factor(plot_data$tissue,
         c("PSC", "ReN"))
plot_data$condition <-
  factor(plot_data$condition,
         c("Norm", "Hypo"))


# Plot with ggpubr
c("EFNB2", "EPHA4", "MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "MCM7") # 
c("EFNB1", "EFNA5", "EPHB2", "EPHA3") # 



library("ggpubr")
my_comparisons <- list( c("Norm", "Hypo") )

pdf("output/tpm_featurecounts_multi/tpm-EPHB2.pdf", width=4, height=3)
tpm_all_sample_tidy %>%
  filter(geneSymbol %in% c("EPHB2") ) %>% 
  unique() %>%
  mutate(TPM = log2(tpm + 1) ) %>%
    ggboxplot(., x = "condition", y = "TPM",
                 fill = "condition",
                 palette = c("blue","red")) +
      # Add the statistical comparisons
      stat_compare_means(comparisons = my_comparisons, 
                        method = "t.test", 
                        aes(group = condition)) +
      theme_bw() +
      facet_wrap(~tissue) +
      ylab("log2(TPM + 1)") +
      ggtitle("EPHB2")
dev.off()






# dotplots with individual replicates
## Hypoxia marker genes 


pdf("output/tpm_featurecounts_multi/tpm-hypoxia_genes-dots_by_rep.pdf", width=5, height=15)
tpm_all_sample_tidy %>%
  filter(geneSymbol %in% c( "VEGFA",  "ADM",  "EGLN3",  "BNIP3",  "BNIP3L",  "CA9",  "CA12",  "LDHA",  "SLC2A1",  "ENO1",  "PFKP",  "HK2",  "ALDOA",  "PDK1")) %>%
  distinct(geneSymbol, tissue, condition, replicate, tpm, .keep_all = TRUE) %>%
  mutate(
    condition = factor(condition, levels = c("Norm", "Hypo")),
    replicate = factor(replicate),
    TPM = log2(tpm + 1)
  ) %>%
  ggplot(aes(x = condition, y = TPM, color = replicate)) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0),
    size = 2.2,
    alpha = 0.9
  ) +
  theme_bw() +
  facet_grid(geneSymbol ~ tissue, scales = "free_y") +
  ylab("log2(TPM + 1)") +
  xlab("") +
  guides(color = guide_legend(title = "Bio Rep"))

dev.off()






```


--> Check Hypoxia marker genes to check whether all replicate not responsive to Hypoxia in ReN; confirm no outlier Hypoxia-response.
  --> No replicate effect: ReN cells globally respond less to Hypoxia, all replicates behave similarly



# DEGs with deseq2 (featurecounts)

**IMPORTANT NOTE: Here it is advisable to REMOVE all genes from chromosome X and Y BEFORE doing the DEGs analysis (X chromosome re-activation occurs in some samples, notably these with more cell passage; in our case, the HET and KO)**
--> It is good to do this on the count matrix see [here](https://support.bioconductor.org/p/119932/)
### 'one-by-one' comparison
Comparison WT vs mutant:
- PSC Hypo vs Norm
- ReN Hypo vs Norm



### PSC Hypo vs Norm - without X/Y chr

```bash
conda activate deseq2
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")
library("rtracklayer")

set.seed(42)

# import GTF for gene name
gtf <- import("../../Master/meta/gencode.v47.annotation.gtf")
## Extract geneId and geneSymbol
gene_table <- mcols(gtf) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct() %>%
  as_tibble()
## Rename columns
colnames(gene_table) <- c("geneId", "geneSymbol")



# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("PSC_Norm_Rep1", "PSC_Norm_Rep2" ,"PSC_Norm_Rep3" ,"PSC_Norm_Rep4" ,"PSC_Hypo_Rep1", "PSC_Hypo_Rep2", "PSC_Hypo_Rep3", "PSC_Hypo_Rep4")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_multi/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/")) %>%
    dplyr::rename(!!sample := starts_with("output/STAR/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 

## Create colData file that describe all our samples
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("celltype", "condition", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = round(counts_all_matrix),
                              colData = coldata,
                              design= ~ condition)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$condition <- relevel(dds$condition, ref = "Norm")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="condition_Hypo_vs_Norm", type="apeglm")


# Add geneSymbol
res_tibble = as_tibble(rownames_to_column(as.data.frame(res), var = "geneId")) %>%
  left_join(gene_table)


# Identify DEGs and count them

## padj 0.05 FC 0.58 ##################################
res_df <- res_tibble %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.58 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.58 & res_df$padj == TRUE, na.rm = TRUE)



## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, 'Sky Blue',
    ifelse(res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.58)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.58)'

pdf("output/deseq2/plotVolcano_res_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Hypo vs Norm, PSC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.58,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) ) +
  annotate("text", x = 3, y = 140, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 140, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()






# Save as gene list for GO analysis:
### Complete table with geneSymbol
write.table(res_tibble, file = "output/deseq2/res_PSC_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/deseq2/upregulated_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/deseq2/downregulated_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


########################################
## Label hypoxia gene ####################
hypoxia_genes <- c(
  "VEGFA",
  "ADM",
  "EGLN3",
  "BNIP3",
  "BNIP3L",
  "CA9",
  "CA12",
  "LDHA",
  "SLC2A1",
  "ENO1",
  "PFKP",
  "HK2",
  "ALDOA",
  "PDK1"
)

pdf("output/deseq2/plotVolcano_res_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-hypoxia_genes.pdf",
    width = 7, height = 8)
EnhancedVolcano(
  res,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Hypo vs Norm, PSC',
  pCutoff = 5e-2,
  FCcutoff = 0.58,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  selectLab = hypoxia_genes,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "grey40",

  legendPosition = 'none'
) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 24)
  ) +
  annotate("text", x = 3, y = 140,
           label = paste(n_upregulated),
           hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 140,
           label = paste(n_downregulated),
           hjust = 0, size = 6, color = "darkred")
dev.off()






########################################
## Label Ephrin gene ####################
Ephrin_CellChat_genes <- c(
  "EFNB2",
  "EFNB1",
  "EFNA5",
  "EPHB2",
  "EPHA4",
  "EPHA3"
)

pdf("output/deseq2/plotVolcano_res_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-Ephrin_CellChat_genes.pdf",
    width = 7, height = 8)
EnhancedVolcano(
  res,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Hypo vs Norm, PSC',
  pCutoff = 5e-2,
  FCcutoff = 0.58,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  selectLab = Ephrin_CellChat_genes,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "grey40",

  legendPosition = 'none'
) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 24)
  ) +
  annotate("text", x = 3, y = 140,
           label = paste(n_upregulated),
           hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 140,
           label = paste(n_downregulated),
           hjust = 0, size = 6, color = "darkred")
dev.off()




```


--> EFNB1 on X chr...







### PSC Hypo vs Norm - including X/Y chr

```bash
conda activate deseq2
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")
library("rtracklayer")

set.seed(42)

# import GTF for gene name
gtf <- import("../../Master/meta/gencode.v47.annotation.gtf")
## Extract geneId and geneSymbol
gene_table <- mcols(gtf) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct() %>%
  as_tibble()
## Rename columns
colnames(gene_table) <- c("geneId", "geneSymbol")



# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("PSC_Norm_Rep1", "PSC_Norm_Rep2" ,"PSC_Norm_Rep3" ,"PSC_Norm_Rep4" ,"PSC_Hypo_Rep1", "PSC_Hypo_Rep2", "PSC_Hypo_Rep3", "PSC_Hypo_Rep4")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_multi/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/")) %>%
    dplyr::rename(!!sample := starts_with("output/STAR/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(dplyr::select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("celltype", "condition", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = round(counts_all_matrix),
                              colData = coldata,
                              design= ~ condition)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$condition <- relevel(dds$condition, ref = "Norm")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="condition_Hypo_vs_Norm", type="apeglm")


# Add geneSymbol
res_tibble = as_tibble(rownames_to_column(as.data.frame(res), var = "geneId")) %>%
  left_join(gene_table)


# Identify DEGs and count them

## padj 0.05 FC 0.58 ##################################
res_df <- res_tibble %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.58 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.58 & res_df$padj == TRUE, na.rm = TRUE)



## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, 'Sky Blue',
    ifelse(res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.58)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.58)'

pdf("output/deseq2/plotVolcano_resXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Hypo vs Norm, PSC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.58,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) ) +
  annotate("text", x = 3, y = 140, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 140, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()






# Save as gene list for GO analysis:
### Complete table with geneSymbol
write.table(res_tibble, file = "output/deseq2/resXinclude_PSC_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/deseq2/upregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/deseq2/downregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


########################################
## Label hypoxia gene ####################
hypoxia_genes <- c(
  "VEGFA",
  "ADM",
  "EGLN3",
  "BNIP3",
  "BNIP3L",
  "CA9",
  "CA12",
  "LDHA",
  "SLC2A1",
  "ENO1",
  "PFKP",
  "HK2",
  "ALDOA",
  "PDK1"
)

pdf("output/deseq2/plotVolcano_resXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-hypoxia_genes.pdf",
    width = 7, height = 8)
EnhancedVolcano(
  res,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Hypo vs Norm, PSC',
  pCutoff = 5e-2,
  FCcutoff = 0.58,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  selectLab = hypoxia_genes,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "grey40",

  legendPosition = 'none'
) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 24)
  ) +
  annotate("text", x = 3, y = 140,
           label = paste(n_upregulated),
           hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 140,
           label = paste(n_downregulated),
           hjust = 0, size = 6, color = "darkred")
dev.off()






########################################
## Label Ephrin gene ####################
Ephrin_CellChat_genes <- c(
  "EFNB2",
  "EFNB1",
  "EFNA5",
  "EPHB2",
  "EPHA4",
  "EPHA3"
)

pdf("output/deseq2/plotVolcano_resXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-Ephrin_CellChat_genes.pdf",
    width = 7, height = 8)
EnhancedVolcano(
  res,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Hypo vs Norm, PSC',
  pCutoff = 5e-2,
  FCcutoff = 0.58,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  selectLab = Ephrin_CellChat_genes,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "grey40",

  legendPosition = 'none'
) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 24)
  ) +
  annotate("text", x = 3, y = 140,
           label = paste(n_upregulated),
           hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 140,
           label = paste(n_downregulated),
           hjust = 0, size = 6, color = "darkred")
dev.off()


















## padj 0.00001 FC 0.58 ##################################
res_df <- res_tibble %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.00001, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 1 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -1 & res_df$padj == TRUE, na.rm = TRUE)



## Plot-volcano
# FILTER ON QVALUE 0.00001 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_tibble$log2FoldChange < -1 & res_tibble$padj < 0.00001, 'Sky Blue',
    ifelse(res_tibble$log2FoldChange > 1 & res_tibble$padj < 0.00001, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.00001; log2FC > 1)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.00001; log2FC < -1)'

pdf("output/deseq2/plotVolcano_resXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Hypo vs Norm, PSC',
  pCutoff = 0.00001,         #
  FCcutoff = 1,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) ) +
  annotate("text", x = 3, y = 140, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 140, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()






# Save as gene list for GO analysis:


### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 1 & res_tibble$padj < 0.00001, ]

#### Filter for down-regulated genes
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -1 & res_tibble$padj < 0.00001, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/deseq2/upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/deseq2/downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)




```


--> The default padj 0.05 and FC 0.58 treshold lead to too many DEGs; and some show weird profile with very high FC (weird) --> Lets prefer use padj 0.00001 (1e-5) and FC > 1




### ReN Hypo vs Norm - without X/Y chr

```bash
conda activate deseq2
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")
library("rtracklayer")


# import GTF for gene name
gtf <- import("../../Master/meta/gencode.v47.annotation.gtf")
## Extract geneId and geneSymbol
gene_table <- mcols(gtf) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct() %>%
  as_tibble()
## Rename columns
colnames(gene_table) <- c("geneId", "geneSymbol")



# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("ReN_Norm_Rep1", "ReN_Norm_Rep2" ,"ReN_Norm_Rep3" ,"ReN_Norm_Rep4" ,"ReN_Hypo_Rep1", "ReN_Hypo_Rep2", "ReN_Hypo_Rep3", "ReN_Hypo_Rep4")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_multi/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/")) %>%
    rename(!!sample := starts_with("output/STAR/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 

## Create colData file that describe all our samples
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("celltype", "condition", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = round(counts_all_matrix),
                              colData = coldata,
                              design= ~ condition)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$condition <- relevel(dds$condition, ref = "Norm")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="condition_Hypo_vs_Norm", type="apeglm")


# Add geneSymbol
res_tibble = as_tibble(rownames_to_column(as.data.frame(res), var = "geneId")) %>%
  left_join(gene_table)


# Identify DEGs and count them

## padj 0.05 FC 0.58 ##################################
res_df <- res_tibble %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.58 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.58 & res_df$padj == TRUE, na.rm = TRUE)



## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, 'Sky Blue',
    ifelse(res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.58)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.58)'

pdf("output/deseq2/plotVolcano_res_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Hypo vs Norm, ReN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.58,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) ) +
  annotate("text", x = 3, y = 140, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -2, y = 140, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()






# Save as gene list for GO analysis:
### Complete table with geneSymbol
write.table(res_tibble, file = "output/deseq2/res_ReN_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/deseq2/upregulated_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/deseq2/downregulated_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



## Label hypoxia gene
hypoxia_genes <- c(
  "VEGFA",
  "ADM",
  "EGLN3",
  "BNIP3",
  "BNIP3L",
  "CA9",
  "CA12",
  "LDHA",
  "SLC2A1",
  "ENO1",
  "PFKP",
  "HK2",
  "ALDOA",
  "PDK1"
)

pdf("output/deseq2/plotVolcano_res_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-hypoxia_genes.pdf",
    width = 7, height = 8)
EnhancedVolcano(
  res,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Hypo vs Norm, ReN',
  pCutoff = 5e-2,
  FCcutoff = 0.58,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  selectLab = hypoxia_genes,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "grey40",

  legendPosition = 'none'
) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 24)
  ) +
  annotate("text", x = 3, y = 140,
           label = paste(n_upregulated),
           hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -2, y = 140,
           label = paste(n_downregulated),
           hjust = 0, size = 6, color = "darkred")
dev.off()



```






### ReN Hypo vs Norm - including X/Y chr

```bash
conda activate deseq2
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")
library("rtracklayer")

set.seed(42)


# import GTF for gene name
gtf <- import("../../Master/meta/gencode.v47.annotation.gtf")
## Extract geneId and geneSymbol
gene_table <- mcols(gtf) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct() %>%
  as_tibble()
## Rename columns
colnames(gene_table) <- c("geneId", "geneSymbol")



# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("ReN_Norm_Rep1", "ReN_Norm_Rep2" ,"ReN_Norm_Rep3" ,"ReN_Norm_Rep4" ,"ReN_Hypo_Rep1", "ReN_Hypo_Rep2", "ReN_Hypo_Rep3", "ReN_Hypo_Rep4")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_multi/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/")) %>%
    rename(!!sample := starts_with("output/STAR/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")



# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(dplyr::select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("celltype", "condition", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = round(counts_all_matrix),
                              colData = coldata,
                              design= ~ condition)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$condition <- relevel(dds$condition, ref = "Norm")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="condition_Hypo_vs_Norm", type="apeglm")


# Add geneSymbol
res_tibble = as_tibble(rownames_to_column(as.data.frame(res), var = "geneId")) %>%
  left_join(gene_table)


# Identify DEGs and count them

## padj 0.05 FC 0.58 ##################################
res_df <- res_tibble %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.58 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.58 & res_df$padj == TRUE, na.rm = TRUE)



## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, 'Sky Blue',
    ifelse(res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.58)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.58)'

pdf("output/deseq2/plotVolcano_resXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Hypo vs Norm, ReN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.58,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) ) +
  annotate("text", x = 3, y = 140, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -2, y = 140, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()






# Save as gene list for GO analysis:
### Complete table with geneSymbol
write.table(res_tibble, file = "output/deseq2/resXinclude_ReN_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/deseq2/upregulatedXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/deseq2/downregulatedXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



## Label hypoxia gene
hypoxia_genes <- c(
  "VEGFA",
  "ADM",
  "EGLN3",
  "BNIP3",
  "BNIP3L",
  "CA9",
  "CA12",
  "LDHA",
  "SLC2A1",
  "ENO1",
  "PFKP",
  "HK2",
  "ALDOA",
  "PDK1"
)

pdf("output/deseq2/plotVolcano_resXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-hypoxia_genes.pdf",
    width = 7, height = 8)
EnhancedVolcano(
  res,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Hypo vs Norm, ReN',
  pCutoff = 5e-2,
  FCcutoff = 0.58,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  selectLab = hypoxia_genes,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "grey40",

  legendPosition = 'none'
) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 24)
  ) +
  annotate("text", x = 3, y = 140,
           label = paste(n_upregulated),
           hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -2, y = 140,
           label = paste(n_downregulated),
           hjust = 0, size = 6, color = "darkred")
dev.off()




########################################
## Label Ephrin gene ####################
Ephrin_CellChat_genes <- c(
  "EFNB2",
  "EFNB1",
  "EFNA5",
  "EPHB2",
  "EPHA4",
  "EPHA3"
)

pdf("output/deseq2/plotVolcano_resXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-Ephrin_CellChat_genes.pdf",
    width = 7, height = 8)
EnhancedVolcano(
  res,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Hypo vs Norm, ReN',
  pCutoff = 5e-2,
  FCcutoff = 0.58,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  selectLab = Ephrin_CellChat_genes,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "grey40",

  legendPosition = 'none'
) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 24)
  ) +
  annotate("text", x = 3, y = 140,
           label = paste(n_upregulated),
           hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 140,
           label = paste(n_downregulated),
           hjust = 0, size = 6, color = "darkred")
dev.off()




```






--> Much **more DEGs in PSC** (also samples looks more clean in their PCA)

--> Core **Hypoxia marker genes are induced in both conditions**!







# Functional analysis with enrichGO (single list of genes dotplot)


We will use clusterProfile package. Tutorial [here](https://hbctraining.github.io/DGE_workshop_salmon/lessons/functional_analysis_2019.html).

Let's do a test of the pipeline with genes from cluster4 amd cluster14 from the rlog counts. Our background list will be all genes tested for differential expression.

**IMPORTANT NOTE: When doing GO, do NOT set a universe (background list of genes) it perform better!**




```R
# packages
library("clusterProfiler")
library("pathview")
library("DOSE")
library("org.Hs.eg.db")
library("enrichplot")
library("rtracklayer")
library("tidyverse")

## Read GTF file
gtf_file <- "../../Master/meta/gencode.v47.annotation.gtf"
gtf_data <- import(gtf_file)

## Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name

## Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()


### GeneSymbol list of signif DEG qval 0.05 FC 0.58
output/deseq2/upregulated_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.txt
output/deseq2/downregulated_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.txt

output/deseq2/upregulated_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.txt
output/deseq2/downregulated_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.txt

output/deseq2/upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.txt
output/deseq2/downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.txt

### Genes with significant splicing changes
output/IsoformSwitchAnalyzeR_kallisto/PSC/significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-PSC.txt

output/IsoformSwitchAnalyzeR_kallisto/PSC/extractConsequenceEnrichment.txt
output/IsoformSwitchAnalyzeR_kallisto/ReN/extractConsequenceEnrichment.txt


############ PSC - UP ############

PSC_up = read_csv("output/deseq2/upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.txt", col_names = "gene_name")

ego <- enrichGO(gene = as.character(PSC_up$gene_name), 
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)
                
pdf("output/GO/dotplot_BP-upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf", width=7, height=7)
dotplot(ego, showCategory=20)
dev.off()

pdf("output/GO/dotplot_BP-upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf", width=5, height=4)
dotplot(ego, showCategory=10)
dev.off()

# Export gene GO output
ego_tbl_clean <- as.data.frame(ego) %>%
  separate(GeneRatio, into = c("GeneHits", "GeneSetSize"), sep = "/", convert = TRUE) %>%
  separate(BgRatio,   into = c("BgHits", "BgSetSize"),     sep = "/", convert = TRUE) %>%
  mutate(
    GeneRatioNumeric = GeneHits / GeneSetSize,
    BgRatioNumeric   = BgHits / BgSetSize
  ) %>%
  rename(
    GO_ID = ID,
    Term  = Description,
    PValue = pvalue,
    FDR = p.adjust,
    QValue = qvalue,
    Genes = geneID
  ) %>%
  select(
    GO_ID, Term, Count,
    GeneHits, GeneSetSize, GeneRatioNumeric,
    BgHits, BgSetSize, BgRatioNumeric,
    PValue, FDR, QValue, Genes
  )
write.table(
  ego_tbl_clean,
  file = "output/GO/BP-upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
###


entrez_genes <- as.character( mapIds(org.Hs.eg.db, as.character(PSC_up$gene_name), 'ENTREZID', 'SYMBOL') )

ekegg <- enrichKEGG(gene = entrez_genes, 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05)
                
pdf("output/GO/dotplot_KEGG-upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf", width=7, height=7)
dotplot(ekegg, showCategory=20)
dev.off()

pdf("output/GO/dotplot_KEGG-upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf", width=5, height=4)
dotplot(ekegg, showCategory=10)
dev.off()


# Export gene KEGG output
eekegg_tbl_clean <- as.data.frame(ekegg) %>%
  separate(GeneRatio, into = c("GeneHits", "GeneSetSize"), sep = "/", convert = TRUE) %>%
  separate(BgRatio,   into = c("BgHits", "BgSetSize"),     sep = "/", convert = TRUE) %>%
  mutate(
    GeneRatioNumeric = GeneHits / GeneSetSize,
    BgRatioNumeric   = BgHits / BgSetSize
  ) %>%
  rename(
    KEGG_ID = ID,
    Term  = Description,
    PValue = pvalue,
    FDR = p.adjust,
    QValue = qvalue,
    Genes = geneID
  ) %>%
  select(
    KEGG_ID, Term, Count,
    GeneHits, GeneSetSize, GeneRatioNumeric,
    BgHits, BgSetSize, BgRatioNumeric,
    PValue, FDR, QValue, Genes
  )
write.table(
  ego_tbl_clean,
  file = "output/GO/KEGG-upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
###





############ PSC - DOWN ############

PSC_down = read_csv("output/deseq2/downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.txt", col_names = "gene_name")

ego <- enrichGO(gene = as.character(PSC_down$gene_name), 
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)
                
pdf("output/GO/dotplot_BP-downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf", width=7, height=7)
dotplot(ego, showCategory=20)
dev.off()

pdf("output/GO/dotplot_BP-downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf", width=5, height=5)
dotplot(ego, showCategory=10)
dev.off()



# Export gene GO output
ego_tbl_clean <- as.data.frame(ego) %>%
  separate(GeneRatio, into = c("GeneHits", "GeneSetSize"), sep = "/", convert = TRUE) %>%
  separate(BgRatio,   into = c("BgHits", "BgSetSize"),     sep = "/", convert = TRUE) %>%
  mutate(
    GeneRatioNumeric = GeneHits / GeneSetSize,
    BgRatioNumeric   = BgHits / BgSetSize
  ) %>%
  rename(
    GO_ID = ID,
    Term  = Description,
    PValue = pvalue,
    FDR = p.adjust,
    QValue = qvalue,
    Genes = geneID
  ) %>%
  select(
    GO_ID, Term, Count,
    GeneHits, GeneSetSize, GeneRatioNumeric,
    BgHits, BgSetSize, BgRatioNumeric,
    PValue, FDR, QValue, Genes
  )
write.table(
  ego_tbl_clean,
  file = "output/GO/BP-downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
###





entrez_genes <- as.character( mapIds(org.Hs.eg.db, as.character(PSC_down$gene_name), 'ENTREZID', 'SYMBOL') )

ekegg <- enrichKEGG(gene = entrez_genes, 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05)
                
pdf("output/GO/dotplot_KEGG-downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf", width=7, height=7)
dotplot(ekegg, showCategory=20)
dev.off()

pdf("output/GO/dotplot_KEGG-downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf", width=5, height=4)
dotplot(ekegg, showCategory=10)
dev.off()



# Export gene KEGG output
eekegg_tbl_clean <- as.data.frame(ekegg) %>%
  separate(GeneRatio, into = c("GeneHits", "GeneSetSize"), sep = "/", convert = TRUE) %>%
  separate(BgRatio,   into = c("BgHits", "BgSetSize"),     sep = "/", convert = TRUE) %>%
  mutate(
    GeneRatioNumeric = GeneHits / GeneSetSize,
    BgRatioNumeric   = BgHits / BgSetSize
  ) %>%
  rename(
    KEGG_ID = ID,
    Term  = Description,
    PValue = pvalue,
    FDR = p.adjust,
    QValue = qvalue,
    Genes = geneID
  ) %>%
  select(
    KEGG_ID, Term, Count,
    GeneHits, GeneSetSize, GeneRatioNumeric,
    BgHits, BgSetSize, BgRatioNumeric,
    PValue, FDR, QValue, Genes
  )
write.table(
  ego_tbl_clean,
  file = "output/GO/KEGG-downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
###





############ PSC - SPLICING CHANGES ALL ############

PSC_splicing = read_csv("output/IsoformSwitchAnalyzeR_kallisto/PSC/significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-PSC.txt", col_names = "gene_name")

ego <- enrichGO(gene = as.character(PSC_splicing$gene_name), 
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)
                
pdf("output/GO/dotplot_BP-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-PSC-top20.pdf", width=7, height=7)
dotplot(ego, showCategory=20)
dev.off()

pdf("output/GO/dotplot_BP-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-PSC-top10.pdf", width=5, height=4)
dotplot(ego, showCategory=10)
dev.off()



# Export gene GO output
ego_tbl_clean <- as.data.frame(ego) %>%
  separate(GeneRatio, into = c("GeneHits", "GeneSetSize"), sep = "/", convert = TRUE) %>%
  separate(BgRatio,   into = c("BgHits", "BgSetSize"),     sep = "/", convert = TRUE) %>%
  mutate(
    GeneRatioNumeric = GeneHits / GeneSetSize,
    BgRatioNumeric   = BgHits / BgSetSize
  ) %>%
  rename(
    GO_ID = ID,
    Term  = Description,
    PValue = pvalue,
    FDR = p.adjust,
    QValue = qvalue,
    Genes = geneID
  ) %>%
  select(
    GO_ID, Term, Count,
    GeneHits, GeneSetSize, GeneRatioNumeric,
    BgHits, BgSetSize, BgRatioNumeric,
    PValue, FDR, QValue, Genes
  )
write.table(
  ego_tbl_clean,
  file = "output/GO/BP-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-PSC.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
###




entrez_genes <- as.character( mapIds(org.Hs.eg.db, as.character(PSC_splicing$gene_name), 'ENTREZID', 'SYMBOL') )

ekegg <- enrichKEGG(gene = entrez_genes, 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05)
                
pdf("output/GO/dotplot_KEGG-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-PSCi-top20.pdf", width=7, height=7)
dotplot(ekegg, showCategory=20)
dev.off()

pdf("output/GO/dotplot_KEGG-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-PSC-top10.pdf", width=5, height=4)
dotplot(ekegg, showCategory=10)
dev.off()




# Export gene KEGG output
eekegg_tbl_clean <- as.data.frame(ekegg) %>%
  separate(GeneRatio, into = c("GeneHits", "GeneSetSize"), sep = "/", convert = TRUE) %>%
  separate(BgRatio,   into = c("BgHits", "BgSetSize"),     sep = "/", convert = TRUE) %>%
  mutate(
    GeneRatioNumeric = GeneHits / GeneSetSize,
    BgRatioNumeric   = BgHits / BgSetSize
  ) %>%
  rename(
    KEGG_ID = ID,
    Term  = Description,
    PValue = pvalue,
    FDR = p.adjust,
    QValue = qvalue,
    Genes = geneID
  ) %>%
  select(
    KEGG_ID, Term, Count,
    GeneHits, GeneSetSize, GeneRatioNumeric,
    BgHits, BgSetSize, BgRatioNumeric,
    PValue, FDR, QValue, Genes
  )
write.table(
  ego_tbl_clean,
  file = "output/GO/KEGG-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-PSC.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
###













############ PSC - SPLICING Consequences ############

PSC_splicing = read_tsv("output/IsoformSwitchAnalyzeR_kallisto/PSC/extractConsequenceEnrichment.txt")


out_pdf <- "output/GO/dotplot_BP-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-extractConsequenceEnrichment-PSC-top10.pdf"
show_n <- 10
df <- PSC_splicing %>%
  filter(!is.na(switchConsequence), switchConsequence != "") %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  distinct(switchConsequence, gene_name)
categories <- sort(unique(df$switchConsequence))
print(categories)
pdf(out_pdf, width = 7, height = 5, onefile = TRUE)
for (cat in categories) {
  genes <- df %>%
    filter(switchConsequence == cat) %>%
    pull(gene_name) %>%
    unique()
  ego <- enrichGO(
    gene          = genes,
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  # handle empty results
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    p <- ggplot() +
      theme_void() +
      ggtitle(paste0("switchConsequence: ", cat, " (n genes = ", length(genes), ")")) +
      annotate("text", x = 0, y = 0,
               label = "No significant GO BP terms (p<0.05)", size = 5) +
      xlim(-1, 1) + ylim(-1, 1)
    print(p)

  } else {
    p <- dotplot(ego, showCategory = show_n) +
      ggtitle(paste0("GO BP — ", cat, " (n genes = ", length(genes), ")")) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }
}
dev.off()



# Export gene GO outputdf <- PSC_splicing %>%
df <- PSC_splicing %>%
  filter(!is.na(switchConsequence), switchConsequence != "") %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  distinct(switchConsequence, gene_name)

categories <- sort(unique(df$switchConsequence))

out_dir <- "output/GO/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# helper to make safe filenames
safe_name <- function(x) {
  x %>%
    gsub("[/\\\\]+", "_", .) %>%        # slashes
    gsub("[^A-Za-z0-9._-]+", "_", .) %>% # weird chars -> _
    gsub("_+", "_", .) %>%              # collapse __
    gsub("^_|_$", "", .)                # trim _
}
# ---- loop: run GO + save table ----
for (cat in categories) {
  genes <- df %>%
    filter(switchConsequence == cat) %>%
    pull(gene_name) %>%
    unique()
  ego <- enrichGO(
    gene          = genes,
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  out_file <- file.path(out_dir, paste0("BP-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-extractConsequenceEnrichment-", safe_name(cat),"-PSC" , ".tsv"))
  # if empty, still write an empty table with headers (nice for pipelines)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    empty <- data.frame(
      GO_ID = character(),
      Description = character(),
      GeneRatio = character(),
      BgRatio = character(),
      pvalue = numeric(),
      p.adjust = numeric(),
      qvalue = numeric(),
      geneID = character(),
      Count = integer()
    )
    write.table(empty, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
    next
  }
  ego_tbl <- as.data.frame(ego)
  # optional: add the category + gene list size as metadata columns
  ego_tbl <- ego_tbl %>%
    mutate(
      switchConsequence = cat,
      nGenesInput = length(genes)
    )
  write.table(ego_tbl, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
}
###





















out_pdf <- "output/GO/KEGG-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-extractConsequenceEnrichment-PSC-top10.pdf"
show_n <- 10

dir.create("output/GO", recursive = TRUE, showWarnings = FALSE)

df <- PSC_splicing %>%
  filter(!is.na(switchConsequence), switchConsequence != "") %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  distinct(switchConsequence, gene_name)

categories <- sort(unique(df$switchConsequence))

pdf(out_pdf, width = 7, height = 5, onefile = TRUE)
for (cat in categories) {
  genes_sym <- df %>%
    filter(switchConsequence == cat) %>%
    pull(gene_name) %>%
    unique()
  # map SYMBOL -> ENTREZ
  entrez <- mapIds(
    org.Hs.eg.db,
    keys     = genes_sym,
    column   = "ENTREZID",
    keytype  = "SYMBOL",
    multiVals = "first"
  )
  entrez_genes <- unique(na.omit(as.character(entrez)))
  # run KEGG (human)
  ekegg <- enrichKEGG(
    gene          = entrez_genes,
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
  if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
    p <- ggplot() +
      theme_void() +
      ggtitle(paste0("KEGG — ", cat,
                     " (n SYMBOL = ", length(genes_sym),
                     ", n ENTREZ = ", length(entrez_genes), ")")) +
      annotate("text", x = 0, y = 0,
               label = "No significant KEGG terms (p<0.05)", size = 5) +
      xlim(-1, 1) + ylim(-1, 1)
    print(p)
  } else {
    p <- dotplot(ekegg, showCategory = show_n) +
      ggtitle(paste0("KEGG — ", cat,
                     " (n SYMBOL = ", length(genes_sym),
                     ", n ENTREZ = ", length(entrez_genes), ")")) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }
}
dev.off()




# Export gene KEGG output
out_dir <- "output/GO"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_name <- function(x) {
  x %>%
    gsub("[/\\\\]+", "_", .) %>%          # slashes
    gsub("[^A-Za-z0-9._-]+", "_", .) %>%  # weird chars -> _
    gsub("_+", "_", .) %>%                # collapse __
    gsub("^_|_$", "", .)                  # trim _
}
df <- PSC_splicing %>%
  filter(!is.na(switchConsequence), switchConsequence != "") %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  distinct(switchConsequence, gene_name)

categories <- sort(unique(df$switchConsequence))

for (cat in categories) {
  genes_sym <- df %>%
    filter(switchConsequence == cat) %>%
    pull(gene_name) %>%
    unique()
  # map SYMBOL -> ENTREZ
  entrez <- mapIds(
    org.Hs.eg.db,
    keys      = genes_sym,
    column    = "ENTREZID",
    keytype   = "SYMBOL",
    multiVals = "first"
  )
  entrez_genes <- unique(na.omit(as.character(entrez)))
  ekegg <- enrichKEGG(
    gene          = entrez_genes,
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )

  out_file <- file.path(out_dir, paste0("KEGG-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-extractConsequenceEnrichment-", safe_name(cat),"-PSC" , ".tsv"))
  # write empty table if nothing significant (keeps pipeline consistent)
  if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
    empty <- data.frame(
      ID = character(),
      Description = character(),
      GeneRatio = character(),
      BgRatio = character(),
      pvalue = numeric(),
      p.adjust = numeric(),
      qvalue = numeric(),
      geneID = character(),
      Count = integer(),
      switchConsequence = character(),
      nGenesSymbol = integer(),
      nGenesEntrez = integer()
    )
    write.table(empty, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
    next
  }
  ekegg_tbl <- as.data.frame(ekegg) %>%
    mutate(
      switchConsequence = cat,
      nGenesSymbol = length(genes_sym),
      nGenesEntrez = length(entrez_genes)
    )
  write.table(ekegg_tbl, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
}
###

















############ ReN - UP ############

ReN_up = read_csv("output/deseq2/upregulated_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.txt", col_names = "gene_name")

ego <- enrichGO(gene = as.character(ReN_up$gene_name), 
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)
                
pdf("output/GO/dotplot_BP-upregulated_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-top20.pdf", width=7, height=7)
dotplot(ego, showCategory=20)
dev.off()

pdf("output/GO/dotplot_BP-upregulated_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-top10.pdf", width=5, height=4)
dotplot(ego, showCategory=10)
dev.off()



# Export gene GO output
ego_tbl_clean <- as.data.frame(ego) %>%
  separate(GeneRatio, into = c("GeneHits", "GeneSetSize"), sep = "/", convert = TRUE) %>%
  separate(BgRatio,   into = c("BgHits", "BgSetSize"),     sep = "/", convert = TRUE) %>%
  mutate(
    GeneRatioNumeric = GeneHits / GeneSetSize,
    BgRatioNumeric   = BgHits / BgSetSize
  ) %>%
  rename(
    GO_ID = ID,
    Term  = Description,
    PValue = pvalue,
    FDR = p.adjust,
    QValue = qvalue,
    Genes = geneID
  ) %>%
  select(
    GO_ID, Term, Count,
    GeneHits, GeneSetSize, GeneRatioNumeric,
    BgHits, BgSetSize, BgRatioNumeric,
    PValue, FDR, QValue, Genes
  )
write.table(
  ego_tbl_clean,
  file = "output/GO/BP-upregulatedresXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
###



entrez_genes <- as.character( mapIds(org.Hs.eg.db, as.character(ReN_up$gene_name), 'ENTREZID', 'SYMBOL') )

ekegg <- enrichKEGG(gene = entrez_genes, 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05)
                
pdf("output/GO/dotplot_KEGG-upregulated_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-top20.pdf", width=7, height=7)
dotplot(ekegg, showCategory=20)
dev.off()

pdf("output/GO/dotplot_KEGG-upregulated_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-top10.pdf", width=5, height=4)
dotplot(ekegg, showCategory=10)
dev.off()


# Export gene KEGG output
eekegg_tbl_clean <- as.data.frame(ekegg) %>%
  separate(GeneRatio, into = c("GeneHits", "GeneSetSize"), sep = "/", convert = TRUE) %>%
  separate(BgRatio,   into = c("BgHits", "BgSetSize"),     sep = "/", convert = TRUE) %>%
  mutate(
    GeneRatioNumeric = GeneHits / GeneSetSize,
    BgRatioNumeric   = BgHits / BgSetSize
  ) %>%
  rename(
    KEGG_ID = ID,
    Term  = Description,
    PValue = pvalue,
    FDR = p.adjust,
    QValue = qvalue,
    Genes = geneID
  ) %>%
  select(
    KEGG_ID, Term, Count,
    GeneHits, GeneSetSize, GeneRatioNumeric,
    BgHits, BgSetSize, BgRatioNumeric,
    PValue, FDR, QValue, Genes
  )
write.table(
  ego_tbl_clean,
  file = "output/GO/KEGG-upregulatedresXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
###



############ ReN - DOWN ############

ReN_down = read_csv("output/deseq2/downregulated_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.txt", col_names = "gene_name")

ego <- enrichGO(gene = as.character(ReN_down$gene_name), 
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)
                
pdf("output/GO/dotplot_BP-downregulated_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-top20.pdf", width=7, height=7)
dotplot(ego, showCategory=20)
dev.off()

pdf("output/GO/dotplot_BP-downregulated_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-top10.pdf", width=5, height=4)
dotplot(ego, showCategory=10)
dev.off()




# Export gene GO output
ego_tbl_clean <- as.data.frame(ego) %>%
  separate(GeneRatio, into = c("GeneHits", "GeneSetSize"), sep = "/", convert = TRUE) %>%
  separate(BgRatio,   into = c("BgHits", "BgSetSize"),     sep = "/", convert = TRUE) %>%
  mutate(
    GeneRatioNumeric = GeneHits / GeneSetSize,
    BgRatioNumeric   = BgHits / BgSetSize
  ) %>%
  rename(
    GO_ID = ID,
    Term  = Description,
    PValue = pvalue,
    FDR = p.adjust,
    QValue = qvalue,
    Genes = geneID
  ) %>%
  select(
    GO_ID, Term, Count,
    GeneHits, GeneSetSize, GeneRatioNumeric,
    BgHits, BgSetSize, BgRatioNumeric,
    PValue, FDR, QValue, Genes
  )
write.table(
  ego_tbl_clean,
  file = "output/GO/BP-downregulatedresXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
###




entrez_genes <- as.character( mapIds(org.Hs.eg.db, as.character(ReN_down$gene_name), 'ENTREZID', 'SYMBOL') )

ekegg <- enrichKEGG(gene = entrez_genes, 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05)
                
pdf("output/GO/dotplot_KEGG-downregulated_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-top20.pdf", width=7, height=7)
dotplot(ekegg, showCategory=20)
dev.off()

pdf("output/GO/dotplot_KEGG-downregulated_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-top10.pdf", width=5, height=4)
dotplot(ekegg, showCategory=10)
dev.off()




# Export gene KEGG output
eekegg_tbl_clean <- as.data.frame(ekegg) %>%
  separate(GeneRatio, into = c("GeneHits", "GeneSetSize"), sep = "/", convert = TRUE) %>%
  separate(BgRatio,   into = c("BgHits", "BgSetSize"),     sep = "/", convert = TRUE) %>%
  mutate(
    GeneRatioNumeric = GeneHits / GeneSetSize,
    BgRatioNumeric   = BgHits / BgSetSize
  ) %>%
  rename(
    KEGG_ID = ID,
    Term  = Description,
    PValue = pvalue,
    FDR = p.adjust,
    QValue = qvalue,
    Genes = geneID
  ) %>%
  select(
    KEGG_ID, Term, Count,
    GeneHits, GeneSetSize, GeneRatioNumeric,
    BgHits, BgSetSize, BgRatioNumeric,
    PValue, FDR, QValue, Genes
  )
write.table(
  ego_tbl_clean,
  file = "output/GO/KEGG-downregulatedresXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
###






############ ReN - SPLICING CHANGES ############

ReN_splicing = read_csv("output/IsoformSwitchAnalyzeR_kallisto/ReN/significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-ReN.txt", col_names = "gene_name")

ego <- enrichGO(gene = as.character(ReN_splicing$gene_name), 
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)
#--> Not signif           
pdf("output/GO/dotplot_BP-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-ReN-top20.pdf", width=7, height=7)
dotplot(ego, showCategory=20)
dev.off()

pdf("output/GO/dotplot_BP-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-ReN-top10.pdf", width=5, height=4)
dotplot(ego, showCategory=10)
dev.off()


entrez_genes <- as.character( mapIds(org.Hs.eg.db, as.character(ReN_splicing$gene_name), 'ENTREZID', 'SYMBOL') )

ekegg <- enrichKEGG(gene = entrez_genes, 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05)
                
pdf("output/GO/dotplot_KEGG-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-ReN-top20.pdf", width=7, height=7)
dotplot(ekegg, showCategory=20)
dev.off()

pdf("output/GO/dotplot_KEGG-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-ReN-top10.pdf", width=5, height=4)
dotplot(ekegg, showCategory=10)
dev.off()





```







# GSEA


Let's do GSEA analysis related to:
- [GOBP_EPHRIN_RECEPTOR_SIGNALING_PATHWAY](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_EPHRIN_RECEPTOR_SIGNALING_PATHWAY.html): Ephrin genes (many seen downregulated from scRNAseq; *cellchat genes downreg: Efnb2, Efnb1, Efna5; and receptors: Ephb2, Epha4, Epha3*)
- [REACTOME_EPHRIN_SIGNALING](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_EPHRIN_SIGNALING.html)


Two method tested; rank on:
- log2fc only
- log2fc and pvalue ranking (combined_score = log2FoldChange * -log10(pvalue) )

```R
# Packages
library("tidyverse")
library("clusterProfiler")
library("msigdbr") # BiocManager::install("msigdbr")
library("org.Mm.eg.db")
library("enrichplot") # for gseaplot2()
library("pheatmap")


set.seed(42)


####################
# FC only ##########
####################
# import DEGs
PSC = read.table("output/deseq2/res_PSC_Hypo_vs_Norm-featurecounts_multi.txt", header = TRUE, sep = "\t") %>%
  as_tibble() 
  
ReN = read.table("output/deseq2/res_ReN_Hypo_vs_Norm-featurecounts_multi.txt", header = TRUE, sep = "\t") %>%
  as_tibble() 

# import gene signature marker lists
EPHRIN = read.table("output/GSEA/geneList-GOBP_EPHRIN_RECEPTOR_SIGNALING_PATHWAY.txt", header = FALSE, sep = "\t") %>%
  as_tibble() %>%
  dplyr::rename("gene" = "V1") %>%
  add_column(term = "EPHRIN")
EPHRIN_REACTOME = read.table("output/GSEA/geneList-REACTOME_EPHRIN_SIGNALING.txt", header = FALSE, sep = "\t") %>%
  as_tibble() %>%
  dplyr::rename("gene" = "V1") %>%
  add_column(term = "EPHRIN_REACTOME")

# Order our DEG
## Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- ReN$log2FoldChange  ### CHAGNE HERE DATA!!!!!!!
names(lfc_vector) <- ReN$geneSymbol ### CHAGNE HERE DATA!!!!!!!
## We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

# run GSEA
## without pvalue
gsea_results <- GSEA(
  geneList = lfc_vector,
  minGSSize = 1,
  maxGSSize = 5000,
  pvalueCutoff = 1,
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = EPHRIN %>% dplyr::select(term,gene), # Need to be in that order...
)

gsea_result_df <- data.frame(gsea_results@result)
# Save output
readr::write_tsv(
  gsea_result_df,
  file.path("output/GSEA/gsea_results-ReN-FC-GOBP_EPHRIN_RECEPTOR_SIGNALING_PATHWAY.tsv"
  )
)

# plots
pdf("output/GSEA/gsea_results-ReN-FC-GOBP_EPHRIN_RECEPTOR_SIGNALING_PATHWAY.pdf", width=8, height=8)
enrichplot::gseaplot(
  gsea_results,
  geneSetID = "EPHRIN",
  title = "GOBP_EPHRIN_RECEPTOR_SIGNALING_PATHWAY",
  color.line = "#0d76ff"
)
dev.off()









####################
# FC+pvalue ########
####################
# import DEGs
PSC = read.table("output/deseq2/res_PSC_Hypo_vs_Norm-featurecounts_multi.txt", header = TRUE, sep = "\t") %>%
  as_tibble() %>%
  mutate(combined_score = log2FoldChange * -log10(pvalue))
  
ReN = read.table("output/deseq2/res_ReN_Hypo_vs_Norm-featurecounts_multi.txt", header = TRUE, sep = "\t") %>%
  as_tibble() %>%
  mutate(combined_score = log2FoldChange * -log10(pvalue))

# import gene signature marker lists
EPHRIN = read.table("output/GSEA/geneList-GOBP_EPHRIN_RECEPTOR_SIGNALING_PATHWAY.txt", header = FALSE, sep = "\t") %>%
  as_tibble() %>%
  dplyr::rename("gene" = "V1") %>%
  add_column(term = "EPHRIN")
EPHRIN_REACTOME = read.table("output/GSEA/geneList-REACTOME_EPHRIN_SIGNALING.txt", header = FALSE, sep = "\t") %>%
  as_tibble() %>%
  dplyr::rename("gene" = "V1") %>%
  add_column(term = "EPHRIN_REACTOME")

# Order our DEG
## Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- ReN$combined_score  ### CHAGNE HERE DATA!!!!!!!
names(lfc_vector) <- ReN$geneSymbol ### CHAGNE HERE DATA!!!!!!!
## We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

# run GSEA
## without pvalue
gsea_results <- GSEA(
  geneList = lfc_vector,
  minGSSize = 1,
  maxGSSize = 5000,
  pvalueCutoff = 1,
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = EPHRIN %>% dplyr::select(term,gene), # Need to be in that order...
)

gsea_result_df <- data.frame(gsea_results@result)
# Save output
readr::write_tsv(
  gsea_result_df,
  file.path("output/GSEA/gsea_results-ReN-FCpvalue-GOBP_EPHRIN_RECEPTOR_SIGNALING_PATHWAY.tsv"
  )
)

# plots
pdf("output/GSEA/gsea_results-ReN-FCpvalue-GOBP_EPHRIN_RECEPTOR_SIGNALING_PATHWAY.pdf", width=8, height=8)
enrichplot::gseaplot(
  gsea_results,
  geneSetID = "EPHRIN",
  title = "GOBP_EPHRIN_RECEPTOR_SIGNALING_PATHWAY",
  color.line = "#0d76ff"
)
dev.off()

```


--> **Nothing significant** for PSC, ReN for GO and REACTOME Ephrin-related terms; both gene ranking method tested





# Cell cycle

Let's use `cyclone` from `scran` to predict cell cycle (method for scRNAseq originally, but require a matrix of gene and expression)
--> Let's try it with log2tpm+1




```bash
conda activate scRNAseqV3
```



```R
# packages
library("scran")
library("tidyverse")

# import tpm
tpm_all_sample=  read.table(  "output/tpm_featurecounts_multi/tpm_all_sample.tsv",  sep = "\t",  header = TRUE) %>%
  as_tibble()

# Remove version suffix (ENSG... .12)
tpm_all_sample2 <- tpm_all_sample %>%
  mutate(Geneid = sub("\\..*$", "", Geneid)) %>%
  distinct(Geneid, .keep_all = TRUE) # keep first if duplicates


# Convert numeric columns to log2(TPM+1)
logtpm_tbl <- tpm_all_sample2 %>%
  mutate(across(-c(Geneid, geneSymbol), ~ log2(.x + 1)))

# Build expression matrix for cyclone: genes x samples, rownames = Ensembl IDs
logtpm_mat <- logtpm_tbl %>%
  column_to_rownames("Geneid") %>%
  select(-geneSymbol) %>%
  as.matrix()


#----------------------------
# 2) Build sample metadata by parsing column names
#   expected pattern: <tissue>_<condition>_<replicate>
#   e.g. PSC_Norm_Rep1, ReN_Hypo_Rep4
#----------------------------
sample_meta <- tibble(sample = colnames(logtpm_mat)) %>%
  separate(sample, into = c("tissue", "condition", "replicate"), sep = "_", remove = FALSE) %>%
  mutate(
    tissue = factor(tissue, levels = c("PSC", "ReN")),
    condition = factor(condition, levels = c("Norm", "Hypo")),
    replicate = factor(replicate, levels = paste0("Rep", 1:4))
  )




#----------------------------
# 3) Load HUMAN cyclone marker pairs
#   scran often ships 'human_cycle_markers.rds'
#   We'll check and stop with a clear message if not present.
#----------------------------
hs_path <- system.file("exdata", "human_cycle_markers.rds", package = "scran")
if (hs_path == "") {
  stop(
    "Could not find 'human_cycle_markers.rds' inside scran.\n",
    "Run: list.files(system.file('exdata', package='scran'), full.names=TRUE)\n",
    "If only mouse markers are present, tell me what files you see and I’ll adapt."
  )
}
hs.pairs <- readRDS(hs_path)

#----------------------------
# 4) Function to run cyclone on a subset of samples
#   Returns a tidy per-sample table with phase + scores.
#----------------------------
run_cyclone_subset <- function(mat_genes_x_samples, meta_subset, pairs) {
  # cyclone expects expression matrix genes x samples
  cyc <- cyclone(mat_genes_x_samples, pairs)

  # Convert results to tidy df
  out <- meta_subset %>%
    select(sample, tissue, condition, replicate) %>%
    mutate(
      phase = cyc$phases,
      G1 = cyc$scores[,"G1"],
      S  = cyc$scores[,"S"],
      G2M = cyc$scores[,"G2M"]
    )

  return(out)
}

#----------------------------
# 5) Run cyclone separately for each tissue x condition
#----------------------------
results_by_group <- sample_meta %>%
  group_split(tissue, condition) %>%
  map_dfr(function(meta_sub) {
    samples <- meta_sub$sample
    mat_sub <- logtpm_mat[, samples, drop = FALSE]
    run_cyclone_subset(mat_sub, meta_sub, hs.pairs)
  })

#----------------------------
# 6) Output
#----------------------------
results_by_group <- results_by_group %>%
  arrange(tissue, condition, replicate)

print(results_by_group)


#----------------------------
# 7) (Optional) quick visualization: scores per group
#----------------------------
# Long format for ggplot if you want:
scores_long <- results_by_group %>%
  pivot_longer(cols = c(G1, S, G2M), names_to = "phase_score", values_to = "score")

# Example quick summaries:
print(results_by_group %>% group_by(tissue, condition) %>%
        summarize(mean_G1 = mean(G1), mean_S = mean(S), mean_G2M = mean(G2M), .groups="drop"))
```


--> Result is weird as 0% cells in S phase; as tool design for scRNAseq not sure it is good to use it here.




# Splicing
## count with Kallisto
Follow instruction [here](https://pachterlab.github.io/kallisto/download.html)


```bash
conda activate kallisto


## run in sbatch
sbatch scripts/kallisto_count_gtf.sh # 62901198 ok

# Convert pseudoalignment to bigwig
conda activate deeptools

sbatch --dependency=afterany:62901198 scripts/TPM_kallisto_bw.sh # 62901271 ok


# Calculate median
conda activate BedToBigwig

sbatch --dependency=afterany:62901271 scripts/bigwigmerge_TPM_kallisto_bw.sh # 62901282 ok
```

- *NOTE: Added `--rf-stranded --genomebam` options for strandness and pseudobam alignemt generation*

--> ~70-80% pseudoalign reads



## Differential alternative mRNA splicing

### Run IsoformSwitchAnalyzeR prerequisets Part1 - webserver CPC2, PFAM, IUPRED2A, SIGNALP

After running `isoformSwitchAnalysisPart1()` from  `## IsoformSwitchAnalyzeR usage`


**Run these in Webserver** :
- coding potential with [CPC2](https://cpc2.gao-lab.org/), I put `output/IsoformSwitchAnalyzeR/isoformSwitchAnalyzeR_isoform_nt.fasta`;
  --> `IsoformSwitchAnalyzeR_kallisto/result_cpc2.txt`
    - PSC: Batch 260101661605011 
    - ReN: Batch 260101270241233 
- protein domain with [PFAM](https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan), manually on the cluster, see below:
  --> OK
- Prediction of Intrinsically Unstructured Proteins with [IUPred2](https://iupred2a.elte.hu/); V3 do not support multi FASTA file with Default: `IUPred2 long disorder (default)`;
  --> `IsoformSwitchAnalyzeR_kallisto/isoformSwitchAnalyzeR_isoform_AA_complete.result`
- Prediction of signal peptide with [SignalP 6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/) with option Eukarya/short output/fast; and save the *Prediction summary*
  - --> `IsoformSwitchAnalyzeR_kallisto/prediction_results_SignalIP6.txt` --> LOOK GOOD


**PFAM reformatting**:

```bash
# Load packages and modules
conda activate deseq2
module load HMMER


# Go to the software
cd ../../Master/software
cd pfam_scan

# Run it
./pfam_scan.py ../../../011_CristanchoLab/002__RNAseq_Preeti/output/IsoformSwitchAnalyzeR_kallisto/PSC/isoformSwitchAnalyzeR_isoform_AA_complete.fasta ../pfamdb/ -out ../../../011_CristanchoLab/002__RNAseq_Preeti/output/pfam/pfam_results_kallisto-PSC.txt
./pfam_scan.py ../../../011_CristanchoLab/002__RNAseq_Preeti/output/IsoformSwitchAnalyzeR_kallisto/ReN/isoformSwitchAnalyzeR_isoform_AA_complete.fasta ../pfamdb/ -out ../../../011_CristanchoLab/002__RNAseq_Preeti/output/pfam/pfam_results_kallisto-ReN.txt

# Re-format the output
cd ../../../011_CristanchoLab/002__RNAseq_Preeti
python3 scripts/reformat_pfam_kallisto-PSC.py
python3 scripts/reformat_pfam_kallisto-ReN.py

```



### Run IsoformSwitchAnalyzeR Part1 and part2




```bash
conda activate IsoformSwitchAnalyzeRv5
```

```R
# packages
library("IsoformSwitchAnalyzeR")
library("rtracklayer")


set.seed(42)


# Kallisto ####################################################
#########
## PSC ##
#########


# Importing the Data
salmonQuant <- importIsoformExpression(
  parentDir = "output/kallisto/",
  pattern = "/PSC_"
)

# metadata file

myDesign <- data.frame(
  sampleID  = colnames(salmonQuant$abundance)[-1],
  condition = sub("^PSC_(Hypo|Norm)_.*", "\\1", colnames(salmonQuant$abundance)[-1])
)
myDesign$condition <- factor(myDesign$condition, levels = c("Norm", "Hypo"))
myDesign <- myDesign[order(myDesign$condition, myDesign$sampleID), ]
rownames(myDesign) <- myDesign$sampleID



# 
aSwitchList <- importRdata(
    isoformCountMatrix   = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = "../../Master/meta/gencode.v47.chr_patch_hapl_scaff.annotation.gtf", # gencode.v47.annotation.gtf gencode.v47.chr_patch_hapl_scaff.annotation.gtf
    isoformNtFasta       = "../../Master/meta/salmon/Homo_sapiens.GRCh38.cdna.all.fa.gz",
    fixStringTieAnnotationProblem = TRUE,
    showProgress = FALSE
)
summary(aSwitchList)


SwitchList <- isoformSwitchAnalysisPart1(
    switchAnalyzeRlist   = aSwitchList,
    pathToOutput = 'output/IsoformSwitchAnalyzeR_kallisto/PSC',
    outputSequences      = TRUE, # change to TRUE whan analyzing your own data 
    prepareForWebServers = TRUE  # change to TRUE if you will use webservers for external sequence analysis
)

## SAVE IMAGE R SESSION
# save.image("output/IsoformSwitchAnalyzeR_kallisto/IsoformSwitchAnalyzeR_kallisto_PSC.RData")
# load("output/IsoformSwitchAnalyzeR_kallisto/IsoformSwitchAnalyzeR_kallisto_PSC.RData")
##


#--> Run in WebServer the CPC2, PFAM, IUPRED2A, SIGNALP

analysSwitchList <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = SwitchList, 
  n                         = 10,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = TRUE,
  pathToCPC2resultFile      = "output/IsoformSwitchAnalyzeR_kallisto/PSC/result_cpc2-PSC.txt",
  pathToPFAMresultFile      = "output/pfam/pfam_results_kallisto_reformat-PSC.txt",
  pathToIUPred2AresultFile  = "output/IsoformSwitchAnalyzeR_kallisto/PSC/isoformSwitchAnalyzeR_isoform_AA_complete-PSC.result",
  pathToSignalPresultFile   = "output/IsoformSwitchAnalyzeR_kallisto/PSC/prediction_results_SignalIP6-PSC.txt",
  outputPlots               = TRUE
)



## SAVE IMAGE R SESSION
# save.image("output/IsoformSwitchAnalyzeR_kallisto/IsoformSwitchAnalyzeR_kallisto_PSC.RData")
# load("output/IsoformSwitchAnalyzeR_kallisto/IsoformSwitchAnalyzeR_kallisto_PSC.RData")
##






# Genome-wide Summaries
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/PSC/extractSwitchOverlap.pdf', onefile = FALSE, height=6, width = 9)
extractSwitchOverlap(
    analysSwitchList,
    filterForConsequences=TRUE,
    plotIsoforms = FALSE
)
dev.off()


pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/PSC/extractConsequenceSummary.pdf', onefile = FALSE, height=5, width = 9)
extractConsequenceSummary(
    analysSwitchList,
    consequencesToAnalyze='all',
    plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
    asFractionTotal = FALSE      # enables analysis of fraction of significant features
)
dev.off()
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/PSC/extractConsequenceSummary_Genes.pdf', onefile = FALSE, height=5, width = 9)
extractConsequenceSummary(
    analysSwitchList,
    consequencesToAnalyze='all',
    plotGenes = TRUE,           # enables analysis of genes (instead of isoforms)
    asFractionTotal = FALSE      # enables analysis of fraction of significant features
)
dev.off()


# Consequence Enrichment Analysis
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/PSC/extractConsequenceEnrichment.pdf', onefile = FALSE, height=6, width = 8)
extractConsequenceEnrichment(
    analysSwitchList,
    consequencesToAnalyze='all',
    analysisOppositeConsequence = TRUE,
    localTheme = theme_bw(base_size = 14), # Increase font size in vignette
    returnResult = TRUE, 
)
dev.off()


# Extract consequence for each gene
res <- extractConsequenceEnrichment(
    analysSwitchList,
    consequencesToAnalyze='all',
    analysisOppositeConsequence = TRUE,
    localTheme = theme_bw(base_size = 14), # Increase font size in vignette
    returnResult = TRUE,
    returnSummary = FALSE # This option add the complete summary table
)
write.table(
  res,
  file = "output/IsoformSwitchAnalyzeR_kallisto/PSC/extractConsequenceEnrichment.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)





# Splicing Enrichment Analysis
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/PSC/extractSplicingEnrichment.pdf', onefile = FALSE, height=4, width = 8)
extractSplicingEnrichment(
    analysSwitchList,
    returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)
dev.off()


# Extract consequence for each gene
res <- extractSplicingEnrichment(
    analysSwitchList,
    returnResult = TRUE,
    returnSummary = FALSE
)
write.table(
  res,
  file = "output/IsoformSwitchAnalyzeR_kallisto/PSC/extractSplicingEnrichment.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)




#Overview Plots
## Volcano like plot
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/PSC/Overview_Plots.pdf', onefile = FALSE, height=3, width = 4)
ggplot(data=analysSwitchList$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
     geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
    facet_wrap( ~ condition_1) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
    theme_bw()
dev.off()

## count nb of isoforms:
significant_isoforms <- analysSwitchList$isoformFeatures %>%
  filter(abs(dIF) > 0.1 & isoform_switch_q_value < 0.05) %>%
  nrow()



## Import GTF to have gene name chromosome 
gtf <- import("../../Master/meta/gencode.v47.chr_patch_hapl_scaff.annotation.gtf")
## Convert to data frame
gtf_df <- as.data.frame(gtf)
## Keep only gene entries
gene_df <- gtf_df %>%
  filter(type == "gene") %>%
  dplyr::select(seqnames, gene_name) %>%
  distinct() %>%
  rename(chromosome = seqnames) %>%
  as_tibble()

# Save output list of genes
## Signif only AND with consequence
significant_isoforms__PSC_geneSymbol =  analysSwitchList$isoformFeatures %>%
  filter(abs(dIF) > 0.1 & isoform_switch_q_value < 0.05) %>%
  filter(switchConsequencesGene == TRUE) %>%
  dplyr::select(gene_name) %>% unique() %>% left_join(gene_df) %>% as_tibble() %>%   dplyr::select(gene_name)
write.table(
  significant_isoforms__PSC_geneSymbol,
  file = "output/IsoformSwitchAnalyzeR_kallisto/PSC/significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-PSC.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

### Switch vs Gene changes:
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/PSC/Overview_Plots2.pdf', onefile = FALSE, height=3, width = 4)
ggplot(data=analysSwitchList$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
    geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) + 
    facet_wrap(~ condition_1) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    geom_hline(yintercept = 0, linetype='dashed') +
    geom_vline(xintercept = 0, linetype='dashed') +
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='Gene log2 fold change', y='dIF') +
    theme_bw()
dev.off()





## Generate plot for a gene - All panels
c("EFNB2", "EFNB1", "EFNA5", "EPHB2", "EPHA4", "EPHA3") # CellChat downreg Ephrin


pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/PSC/switchPlot-NormHypo-EPHA3.pdf', onefile = FALSE, height=6, width = 9)
switchPlot(analysSwitchList, gene= "EPHA3", condition1= "Norm", condition2= "Hypo", reverseMinus = FALSE)
dev.off()


## Only gene and isoform expression
pdf("output/IsoformSwitchAnalyzeR_kallisto/PSC/switchPlotGeneExp-NormHypo-EFNB2.pdf", width=3, height=5)
switchPlotGeneExp(
  switchAnalyzeRlist = analysSwitchList,
  gene = "EFNB2",
  condition1 = "Norm",
  condition2 = "Hypo"
)  +
  scale_fill_manual(values = c("Norm" = "blue", "Hypo" = "red")) +
  scale_color_manual(values = c("Norm" = "blue", "Hypo" = "red"))
dev.off()



## SAVE IMAGE R SESSION
# save.image("output/IsoformSwitchAnalyzeR_kallisto/IsoformSwitchAnalyzeR_kallisto_PSC.RData")
# load("output/IsoformSwitchAnalyzeR_kallisto/IsoformSwitchAnalyzeR_kallisto_PSC.RData")
##






#########
## ReN ##
#########


# Importing the Data
salmonQuant <- importIsoformExpression(
  parentDir = "output/kallisto/",
  pattern = "/ReN_"
)

# metadata file
myDesign <- data.frame(
  sampleID  = colnames(salmonQuant$abundance)[-1],
  condition = sub("^ReN_(Hypo|Norm)_.*", "\\1", colnames(salmonQuant$abundance)[-1])
)
myDesign$condition <- factor(myDesign$condition, levels = c("Norm", "Hypo"))
myDesign <- myDesign[order(myDesign$condition, myDesign$sampleID), ]
rownames(myDesign) <- myDesign$sampleID

myDesign
# 
aSwitchList <- importRdata(
    isoformCountMatrix   = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = "../../Master/meta/gencode.v47.chr_patch_hapl_scaff.annotation.gtf", # gencode.v47.annotation.gtf gencode.v47.chr_patch_hapl_scaff.annotation.gtf
    isoformNtFasta       = "../../Master/meta/salmon/Homo_sapiens.GRCh38.cdna.all.fa.gz",
    fixStringTieAnnotationProblem = TRUE,
    showProgress = FALSE
)
summary(aSwitchList)




SwitchList <- isoformSwitchAnalysisPart1(
    switchAnalyzeRlist   = aSwitchList,
    pathToOutput = 'output/IsoformSwitchAnalyzeR_kallisto/ReN',
    outputSequences      = TRUE, # change to TRUE whan analyzing your own data 
    prepareForWebServers = TRUE  # change to TRUE if you will use webservers for external sequence analysis
)

## SAVE IMAGE R SESSION
# save.image("output/IsoformSwitchAnalyzeR_kallisto/IsoformSwitchAnalyzeR_kallisto_ReN.RData")
# load("output/IsoformSwitchAnalyzeR_kallisto/IsoformSwitchAnalyzeR_kallisto_ReN.RData")
##


#--> Run in WebServer the CPC2, PFAM, IUPRED2A, SIGNALP

analysSwitchList <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = SwitchList, 
  n                         = 10,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = TRUE,
  pathToCPC2resultFile      = "output/IsoformSwitchAnalyzeR_kallisto/ReN/result_cpc2-ReN.txt",
  pathToPFAMresultFile      = "output/pfam/pfam_results_kallisto_reformat-ReN.txt",
  pathToIUPred2AresultFile  = "output/IsoformSwitchAnalyzeR_kallisto/ReN/isoformSwitchAnalyzeR_isoform_AA_complete-ReN.result",
  pathToSignalPresultFile   = "output/IsoformSwitchAnalyzeR_kallisto/ReN/prediction_results_SignalIP6-ReN.txt",
  outputPlots               = TRUE
)

## SAVE IMAGE R SESSION
# save.image("output/IsoformSwitchAnalyzeR_kallisto/IsoformSwitchAnalyzeR_kallisto_ReN.RData")
# load("output/IsoformSwitchAnalyzeR_kallisto/IsoformSwitchAnalyzeR_kallisto_ReN.RData")
##






# Genome-wide Summaries
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/ReN/extractSwitchOverlap.pdf', onefile = FALSE, height=6, width = 9)
extractSwitchOverlap(
    analysSwitchList,
    filterForConsequences=TRUE,
    plotIsoforms = FALSE
)
dev.off()


pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/ReN/extractConsequenceSummary.pdf', onefile = FALSE, height=5, width = 9)
extractConsequenceSummary(
    analysSwitchList,
    consequencesToAnalyze='all',
    plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
    asFractionTotal = FALSE      # enables analysis of fraction of significant features
)
dev.off()
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/ReN/extractConsequenceSummary_Genes.pdf', onefile = FALSE, height=5, width = 9)
extractConsequenceSummary(
    analysSwitchList,
    consequencesToAnalyze='all',
    plotGenes = TRUE,           # enables analysis of genes (instead of isoforms)
    asFractionTotal = FALSE      # enables analysis of fraction of significant features
)
dev.off()


# Consequence Enrichment Analysis
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/ReN/extractConsequenceEnrichment.pdf', onefile = FALSE, height=4, width = 8)
extractConsequenceEnrichment(
    analysSwitchList,
    consequencesToAnalyze='all',
    analysisOppositeConsequence = TRUE,
    localTheme = theme_bw(base_size = 14), # Increase font size in vignette
    returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)
dev.off()



# Extract consequence for each gene
res <- extractConsequenceEnrichment(
    analysSwitchList,
    consequencesToAnalyze='all',
    analysisOppositeConsequence = TRUE,
    localTheme = theme_bw(base_size = 14), # Increase font size in vignette
    returnResult = TRUE,
    returnSummary = FALSE # This option add the complete summary table
)
write.table(
  res,
  file = "output/IsoformSwitchAnalyzeR_kallisto/ReN/extractConsequenceEnrichment.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)





# Splicing Enrichment Analysis
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/ReN/extractSplicingEnrichment.pdf', onefile = FALSE, height=4, width = 8)
extractSplicingEnrichment(
    analysSwitchList,
    returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)
dev.off()


# Extract consequence for each gene
res <- extractSplicingEnrichment(
    analysSwitchList,
    returnResult = TRUE,
    returnSummary = FALSE
)
write.table(
  res,
  file = "output/IsoformSwitchAnalyzeR_kallisto/ReN/extractSplicingEnrichment.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)





#Overview Plots
## Volcano like plot
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/ReN/Overview_Plots.pdf', onefile = FALSE, height=3, width = 4)
ggplot(data=analysSwitchList$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
     geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
    facet_wrap( ~ condition_1) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
    theme_bw()
dev.off()

## count nb of isoforms:
significant_isoforms <- analysSwitchList$isoformFeatures %>%
  filter(abs(dIF) > 0.1 & isoform_switch_q_value < 0.05) %>%
  nrow()



## Import GTF to have gene name chromosome 
gtf <- import("../../Master/meta/gencode.v47.chr_patch_hapl_scaff.annotation.gtf")
## Convert to data frame
gtf_df <- as.data.frame(gtf)
## Keep only gene entries
gene_df <- gtf_df %>%
  filter(type == "gene") %>%
  dplyr::select(seqnames, gene_name) %>%
  distinct() %>%
  rename(chromosome = seqnames) %>%
  as_tibble()

# Save output list of genes
## Signif only AND with consequence
significant_isoforms__ReN_geneSymbol =  analysSwitchList$isoformFeatures %>%
  filter(abs(dIF) > 0.1 & isoform_switch_q_value < 0.05) %>%
  filter(switchConsequencesGene == TRUE) %>%
  dplyr::select(gene_name) %>% unique() %>% left_join(gene_df) %>% as_tibble() %>%   dplyr::select(gene_name)
write.table(
  significant_isoforms__ReN_geneSymbol,
  file = "output/IsoformSwitchAnalyzeR_kallisto/ReN/significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-ReN.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

### Switch vs Gene changes:
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/ReN/Overview_Plots2.pdf', onefile = FALSE, height=3, width = 4)
ggplot(data=analysSwitchList$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
    geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) + 
    facet_wrap(~ condition_1) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    geom_hline(yintercept = 0, linetype='dashed') +
    geom_vline(xintercept = 0, linetype='dashed') +
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='Gene log2 fold change', y='dIF') +
    theme_bw()
dev.off()





## Generate plot for a gene - All panels
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/ReN/switchPlot-NormHypo-GAB1.pdf', onefile = FALSE, height=6, width = 9)
switchPlot(analysSwitchList, gene= "GAB1", condition1= "Norm", condition2= "Hypo", reverseMinus = FALSE)
dev.off()


## Only gene and isoform expression
pdf("output/IsoformSwitchAnalyzeR_kallisto/ReN/switchPlotGeneExp-NormHypo-RRM2.pdf", width=3, height=5)
switchPlotGeneExp(
  switchAnalyzeRlist = analysSwitchList,
  gene = "RRM2",
  condition1 = "Norm",
  condition2 = "Hypo"
)  +
  scale_fill_manual(values = c("Norm" = "blue", "Hypo" = "red")) +
  scale_color_manual(values = c("Norm" = "blue", "Hypo" = "red"))
dev.off()



## SAVE IMAGE R SESSION
# save.image("output/IsoformSwitchAnalyzeR_kallisto/IsoformSwitchAnalyzeR_kallisto_ReN.RData")
# load("output/IsoformSwitchAnalyzeR_kallisto/IsoformSwitchAnalyzeR_kallisto_ReN.RData")
##


```


--> There should be ways to investigate the list of isoforms more carefully but I m just gonna do a GO..
  --> **Gene list of interest**: Gene with isforom switch and a functional consequence (loss of coding potential... etc..) `output/IsoformSwitchAnalyzeR_kallisto/significant_isoforms_dIF01qval05switchConsequencesGeneTRUE__[OEKO or KO]_geneSymbol.txt`
    --> Check in `001*/018*` whether these genes are bound with EZH1 (Venn overlap using consens peak does ot show high overlap)


--> **IsoformSwitchAnalyzeR evaluates changes within each isoform, so the X chromosome doesn’t bias the overall analysis**. I can simply exclude isoform switches detected on chrX.



More info [here](https://bioconductor.org/packages/devel/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html)

