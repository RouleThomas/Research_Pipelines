--> This folder is the **clean + reproducible** version of: `001_CristanchoLab/002__RNAseq_Preeti` (original exploration)

# Data Overview

- **Data type**: Bulk RNA-seq
- **Library**: Directional mRNA (poly(A) enrichment)
- **Platform**: Illumina NovaSeq X Plus (paired-end, 150 bp)
- **Factors**:
  - **Cell type**: ReN, PSC
  - **Condition**: Normoxia (24h), Hypoxia (24h)
  - **Replicates**: 4 biological replicates per group (total n = 16)



# Data processing

## Sample Renaming

| Key | Original sample name | New sample name |
| --: | -------------------- | --------------- |
|  A1 | ReN_Nor_24h_Exp3     | ReN_Norm_Rep1   |
|  A2 | ReN_Nor_24h_Exp4     | ReN_Norm_Rep2   |
|  A3 | ReN_Nor_24h_Exp5     | ReN_Norm_Rep3   |
|  A4 | ReN_Nor_24h_Exp6     | ReN_Norm_Rep4   |
|  A5 | ReN_Hyp_24h_Exp3     | ReN_Hypo_Rep1   |
|  A6 | ReN_Hyp_24h_Exp4     | ReN_Hypo_Rep2   |
|  A7 | ReN_Hyp_24h_Exp5     | ReN_Hypo_Rep3   |
|  A8 | ReN_Hyp_24h_Exp6     | ReN_Hypo_Rep4   |
|  B1 | hiPSCs_Nor_24h_Exp3  | PSC_Norm_Rep1   |
|  B2 | hiPSCs_Nor_24h_Exp4  | PSC_Norm_Rep2   |
|  B3 | hiPSCs_Nor_24h_Exp5  | PSC_Norm_Rep3   |
|  B4 | hiPSCs_Nor_24h_Exp6  | PSC_Norm_Rep4   |
|  B5 | hiPSCs_Hyp_24h_Exp3  | PSC_Hypo_Rep1   |
|  B6 | hiPSCs_Hyp_24h_Exp4  | PSC_Hypo_Rep2   |
|  B7 | hiPSCs_Hyp_24h_Exp5  | PSC_Hypo_Rep3   |
|  B8 | hiPSCs_Hyp_24h_Exp6  | PSC_Hypo_Rep4   |



## Read pre-processing & Quality filtering

Reads were trimmed using **fastp (VERSION)** with default adapter detection and quality filtering parameters.

The following command was run for each (paired-end) sample:
```bash
fastp \
  -i <sample>_1.fq.gz \
  -I <sample>_2.fq.gz \
  -o <sample>_1.fq.gz \
  -O <sample>_2.fq.gz \
  -j <sample>.json \
  -h <sample>.html
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/fastp.sh`


## Read mapping

Reads were aligned to the human reference genome (**hg38**, *GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta*) using **STAR** (v2.7.3a).  
Resulting BAM files were sorted by coordinate and indexed using **SAMtools** (v1.16.1).

The following command was run for each (paired-end) sample:
```bash
STAR \
  --genomeDir <genome_index> \
  --readFilesIn <sample>_1.fq.gz <sample>_2.fq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate

samtools index <sample>_Aligned.sortedByCoord.out.bam
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/STAR_mapping_fastp.sh`


## Gene-level quantification

Gene-level read counts were generated using **featureCounts** (v2.0.1) with the **GENCODE v47** gene annotation.

The following command was run for each sample:
```bash
featureCounts -p -C -O -M --fraction -s 2
	-a <annotation.gtf> \
	-o <sample>.txt <sample>_Aligned.sortedByCoord.out.bam
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/featurecounts_multi.sh`

**OUTPUT FILES**:
- Output of featureCounts: `output/featurecounts_multi/[sample name].txt`

## Isoform-level quantification


Isoform-level read counts were generated using **kallisto** (v0.44.0) with the **GENCODE v47** gene annotation.


The following command was run for each sample:
```bash
kallisto quant \
  -i <transcripts.idx> \
  -o <sample>_quant \
  -b 100 \
  -g <annotation.gtf> \
  --rf-stranded \
  --genomebam \
  --chromosomes GRCh38_chrom_sizes.tab \
  <sample>_1.fq.gz <sample>_2.fq.gz
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/kallisto_count_gtf.sh`

**OUTPUT FILES**:
- Output of kallisto: `output/kallisto/[sample name]_quant/abundance.tsv`


## Transcript abundance estimation (TPM/RPKM)

Gene-level TPM and RPKM values were computed from the featureCounts output using a custom R script `RPKM_TPM_featurecounts.R`. For each input featureCounts table, the script converts raw counts to TPM/RPKM using the feature length column (Length) and outputs two tab-separated files: `*_tpm.txt` and `*_rpkm.txt`.

The following command was run for each sample:
```bash
Rscript scripts/RPKM_TPM_featurecounts.R <featureCounts_output>.txt <output_prefix>
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/featurecounts_TPM.sh`

**OUTPUT FILES**:
- TPM count per sample: `output/tpm_featurecounts_multi/[sample name].txt`


## Coverage track & QC
### Generation of bigwig coverage file
Genome-wide coverage tracks (.bigWig) were generated from coordinate-sorted BAM files using **bamCoverage** (v3.5.1). Coverage was normalized using BPM (bins per million mapped reads) with single–base resolution.

The following command was run for each sample:
```bash
bamCoverage \
  --bam <sample>_Aligned.sortedByCoord.out.bam \
  --outFileName <sample>.bw \
  --outFileFormat bigwig \
  --normalizeUsing BPM \
  --binSize 1
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/bigwigmerge_STAR_TPM_bw.sh `

For each biological condition, replicate bigWig files were aggregated by computing the median coverage signal across using **wiggletools** (v1.0.0). Median coverage tracks were first generated in bedGraph format, subsequently sorted using **bedtools** (v2.30.0), and finally converted back to bigWig format using **bedGraphToBigWig** (v4).

```bash
# Compute median coverage across replicates
wiggletools write_bg <condition>_median.bedGraph \
  median <rep1>.bw <rep2>.bw <rep3>.bw <rep4>.bw
# Sort bedGraph
bedtools sort -i <condition>_median.bedGraph > <condition>_median.sorted.bedGraph
# Convert bedGraph to bigWig
bedGraphToBigWig \
  <condition>_median.sorted.bedGraph \
  GRCh38_chrom_sizes.tab \
  <condition>_median.bw
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/bigwigmerge_STAR_TPM_bw.sh`

**OUTPUT FILES**:
- Bigwig coverage files: `output/bigwig/[ReN or PSC]_[Norm or Hypo]_[median or replicate].bw`

### QC - PCA

BPM-normalized genome-wide coverage tracks (.bigWig) were used for Principal component analysis (PCA) using **multiBigwigSummary** (v3.5.1) which compile a matrix of per-sample coverage values, then vizualized with **plotPCA** (v3.5.1). 

The following command was run for PSC, and ReN samples:
```bash
# Compute matrix
multiBigwigSummary bins \
    -b output/bigwig/PSC_Norm_Rep1.bw \
       output/bigwig/PSC_Norm_Rep2.bw \
       output/bigwig/PSC_Norm_Rep3.bw \
       output/bigwig/PSC_Norm_Rep4.bw \
       output/bigwig/PSC_Hypo_Rep1.bw \
       output/bigwig/PSC_Hypo_Rep2.bw \
       output/bigwig/PSC_Hypo_Rep3.bw \
       output/bigwig/PSC_Hypo_Rep4.bw \
    -o output/bigwig/multiBigwigSummary_TPM_PSC.npz

# Generate PCA plot
plotPCA \
    -in output/bigwig/multiBigwigSummary_TPM_PSC.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_Norm_Rep1 PSC_Norm_Rep2 PSC_Norm_Rep3 PSC_Norm_Rep4 \
             PSC_Hypo_Rep1 PSC_Hypo_Rep2 PSC_Hypo_Rep3 PSC_Hypo_Rep4 \
    --colors blue blue blue blue red red red red \
    --markers 's' 'o' '>' 'x' 's' 'o' '>' 'x' \
    -o output/bigwig/multiBigwigSummary_TPM_PSC_plotPCA.pdf \
    --plotWidth 6 \
    --plotHeight 7
```

**OUTPUT FILES**:
- PCA plot of PSC samples: `output/bigwig/multiBigwigSummary_TPM_PSC_plotPCA.pdf`
- PCA plot of ReN samples: `output/bigwig/multiBigwigSummary_ReN_plotPCA.pdf`


# Data analysis

## Gene expression (TPM)

Let's display transcript abundance in TPM for some genes of interest in **R** (v4.2.2):

```R
# Load packages
library("tidyverse") # v2.0.0
library("rtracklayer") # v1.58.0
library("ggpubr") # v0.6.0

set.seed(42) # set seed for reproducibility

# Import gene annotation file to convert gene ID to gene symbol
gtf <- import("../../Master/meta/gencode.v47.annotation.gtf")
gene_table <- mcols(gtf) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct() %>%
  as_tibble()
colnames(gene_table) <- c("Geneid", "geneSymbol")


# Import TPM counts for each sample and merge them in a unique table
## Collect samples IDs
samples <- c("PSC_Norm_Rep1", "PSC_Norm_Rep2" ,"PSC_Norm_Rep3" ,"PSC_Norm_Rep4" ,"PSC_Hypo_Rep1", "PSC_Hypo_Rep2", "PSC_Hypo_Rep3", "PSC_Hypo_Rep4", "ReN_Norm_Rep1", "ReN_Norm_Rep2" ,"ReN_Norm_Rep3" ,"ReN_Norm_Rep4" ,"ReN_Hypo_Rep1", "ReN_Hypo_Rep2", "ReN_Hypo_Rep3", "ReN_Hypo_Rep4")

## Loop over samples and import TPM
sample_data <- list()
for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/tpm_featurecounts_multi/", sample, "_tpm.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    dplyr::select(Geneid, starts_with("output.STAR.")) %>%
    dplyr::rename(!!sample := starts_with("output.STAR."))
}

## Merge all samples into a single table and add gene symbol
tpm_all_sample <- purrr::reduce(sample_data, full_join, by = "Geneid")  %>%
  left_join(gene_table)
###### --- Export --- ######
write.table(tpm_all_sample,  file = "output/tpm_featurecounts_multi/tpm_all_sample.tsv",  sep = "\t",  row.names = FALSE,  quote = FALSE)
###### -------------- ######


# Re-shape TPM counts table
## Shape to long format and extract sample metada from the column names
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

## Plot example: EPHB2 (log2(TPM+1)), Norm vs Hypo, for each tissue
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


## Replicate-level TPM expression of hypoxia marker genes
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

**OUTPUT FILES**:
- Gene-level TPM expression matrix for all samples: `output/tpm_featurecounts_multi/tpm_all_sample.tsv`
- Gene-specific TPM boxplots: `output/tpm_featurecounts_multi/tpm-[GENE OF INTEREST]`
- Replicate-level TPM expression of hypoxia marker genes: `output/tpm_featurecounts_multi/tpm-hypoxia_genes-dots_by_rep.pdf`


**CONCLUSION**:

- Due to ReN cells showing far less DEGs than PSC; we checked whether some replicates were more or less sensitive to the Hypoxia condition by checking hypoxia marker genes and label replicates: **No replicate effect detected; all Hypoxia marker genes respond similarly between replicates**.
  - **ReN is less responsive that PSC to Hypoxia**
  

## Differential Gene Expression 

Differentially expressed genes (DEGs) between normoxia and hypoxia in PSC and ReN cells in **R** (v4.2.2).

### PSC

```R
# Load packages
library("tidyverse") # v2.0.0
library("rtracklayer") # v1.58.0
library("DESeq2") # v1.38.3
library("EnhancedVolcano") # v1.16.0
library("apeglm") # v1.20.0

set.seed(42) # set seed for reproducibility

# Import gene annotation file to convert gene ID to gene symbol
gtf <- import("../../Master/meta/gencode.v47.annotation.gtf")
gene_table <- mcols(gtf) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct() %>%
  as_tibble()
colnames(gene_table) <- c("Geneid", "geneSymbol")

# Import featureCounts counts for each sample and merge them in a unique table
## Collect samples IDs
samples <- c("PSC_Norm_Rep1", "PSC_Norm_Rep2" ,"PSC_Norm_Rep3" ,"PSC_Norm_Rep4" ,"PSC_Hypo_Rep1", "PSC_Hypo_Rep2", "PSC_Hypo_Rep3", "PSC_Hypo_Rep4")

## Loop over samples and import counts
sample_data <- list()
for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_multi/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/")) %>%
    dplyr::rename(!!sample := starts_with("output/STAR/"))
}

## Merge all samples into a single table 
counts_all <- purrr::reduce(sample_data, full_join, by = "Geneid")

# Prepare table for DESEQ2
## Convert dataframe to matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
counts_all_matrix = make_matrix(dplyr::select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Construct sample metadata information
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("celltype", "condition", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## Convert sample metadata to matrix format
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

# Run DESEQ2
## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = round(counts_all_matrix),
                              colData = coldata,
                              design= ~ condition)

## Filter out lowly expressed genes (minimum total count ≥ 5)
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Set reference (control) condition
dds$condition <- relevel(dds$condition, ref = "Norm")

## Run DESeq2 model fitting
dds <- DESeq(dds)

## Extract results with shrunken log2 fold changes
res <- lfcShrink(dds, coef = "condition_Hypo_vs_Norm", type = "apeglm")

## Add gene symbols to DESeq2 results
res_tibble = as_tibble(rownames_to_column(as.data.frame(res), var = "Geneid")) %>%
  left_join(gene_table)

# Identify DEGs
############################################################
####### adjusted p-value < 0.05 and |log2FC| > 0.58 ########
############################################################
res_df <- res_tibble %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.58 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.58 & res_df$padj == TRUE, na.rm = TRUE)

### Volcano plot of differential expression results
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


### Volcano plot of differential expression results with hypoxia marker genes labelled
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

### Volcano plot of differential expression results with ephrin genes labelled
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


############################################################
####### adjusted p-value < 1e-5 and |log2FC| > 1 ########
############################################################

res_df <- res_tibble %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.00001, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 1 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -1 & res_df$padj == TRUE, na.rm = TRUE)

### Volcano plot of differential expression results
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


###### --- Export --- ######
# DESEQ2 result table
write.table(res_tibble, file = "output/deseq2/resXinclude_PSC_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# Signif up and down -regulated genes
## adjusted p-value < 0.05 and |log2FC| > 0.58
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, ]
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, ]
write.table(upregulated$geneSymbol, file = "output/deseq2/upregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/deseq2/downregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
## adjusted p-value < 1e-5 and |log2FC| > 1
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 1 & res_tibble$padj < 0.00001, ]
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -1 & res_tibble$padj < 0.00001, ]
write.table(upregulated$geneSymbol, file = "output/deseq2/upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/deseq2/downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# R session info
loaded_pkgs <- sessionInfo()$otherPkgs
pkg_versions <- data.frame(
  package = names(loaded_pkgs),
  version = sapply(loaded_pkgs, function(x) x$Version)
)
write.table(
  pkg_versions,
  file = "output/deseq2/RsessionInfo-DESEQ2_DEG.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

###### -------------- ######
```

**OUTPUT FILES:**
- Complete DESEQ2 output table: `output/deseq2/resXinclude_PSC_Hypo_vs_Norm-featurecounts_multi.txt`
- Volcano plot of DEGs *adjusted p-value < 0.05 and |log2FC| > 0.58* : `output/deseq2/plotVolcano_resXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.pdf`
- Volcano plot of DEGs with hypoxia marker genes labelled: `output/deseq2/plotVolcano_resXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-hypoxia_genes.pdf`
- Volcano plot of DEGs with ephrin genes labelled: `output/deseq2/plotVolcano_resXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-Ephrin_CellChat_genes.pdf`
- Volcano plot of DEGs *adjusted p-value < 1e-5 and |log2FC| > 1*: `output/deseq2/plotVolcano_resXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.pdf`
- List of up regulated genes *adjusted p-value < 0.05 and |log2FC| > 0.58*: `output/deseq2/upregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.txt`
- List of up regulated genes *adjusted p-value < 1e-5 and |log2FC| > 1*: `output/deseq2/upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.txt`
- List of down regulated genes *adjusted p-value < 0.05 and |log2FC| > 0.58*: `output/deseq2/downregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.txt`
- List of down regulated genes *adjusted  p-value < 1e-5 and |log2FC| > 1*: `output/deseq2/downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.txt`
- R packages version: `output/deseq2/RsessionInfo-DESEQ2_DEG.txt`





### ReN

```R
# Load packages
library("tidyverse") # v2.0.0
library("rtracklayer") # v1.58.0
library("DESeq2") # v1.38.3
library("EnhancedVolcano") # v1.16.0
library("apeglm") # v1.20.0

set.seed(42) # set seed for reproducibility

# Import gene annotation file to convert gene ID to gene symbol
gtf <- import("../../Master/meta/gencode.v47.annotation.gtf")
gene_table <- mcols(gtf) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct() %>%
  as_tibble()
colnames(gene_table) <- c("Geneid", "geneSymbol")

# Import featureCounts counts for each sample and merge them in a unique table
## Collect samples IDs
samples <- c("ReN_Norm_Rep1", "ReN_Norm_Rep2" ,"ReN_Norm_Rep3" ,"ReN_Norm_Rep4" ,"ReN_Hypo_Rep1", "ReN_Hypo_Rep2", "ReN_Hypo_Rep3", "ReN_Hypo_Rep4")

## Loop over samples and import counts
sample_data <- list()
for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_multi/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/")) %>%
    dplyr::rename(!!sample := starts_with("output/STAR/"))
}

## Merge all samples into a single table 
counts_all <- purrr::reduce(sample_data, full_join, by = "Geneid")

# Prepare table for DESEQ2
## Convert dataframe to matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
counts_all_matrix = make_matrix(dplyr::select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Construct sample metadata information
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("celltype", "condition", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## Convert sample metadata to matrix format
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

# Run DESEQ2
## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = round(counts_all_matrix),
                              colData = coldata,
                              design= ~ condition)

## Filter out lowly expressed genes (minimum total count ≥ 5)
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Set reference (control) condition
dds$condition <- relevel(dds$condition, ref = "Norm")

## Run DESeq2 model fitting
dds <- DESeq(dds)

## Extract results with shrunken log2 fold changes
res <- lfcShrink(dds, coef = "condition_Hypo_vs_Norm", type = "apeglm")

## Add gene symbols to DESeq2 results
res_tibble = as_tibble(rownames_to_column(as.data.frame(res), var = "Geneid")) %>%
  left_join(gene_table)

# Identify DEGs
############################################################
####### adjusted p-value < 0.05 and |log2FC| > 0.58 ########
############################################################
res_df <- res_tibble %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.58 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.58 & res_df$padj == TRUE, na.rm = TRUE)

### Volcano plot of differential expression results
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


### Volcano plot of differential expression results with hypoxia marker genes labelled
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

### Volcano plot of differential expression results with ephrin genes labelled
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




###### --- Export --- ######
# DESEQ2 result table
write.table(res_tibble, file = "output/deseq2/resXinclude_ReN_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# Signif up and down -regulated genes
## adjusted p-value < 0.05 and |log2FC| > 0.58
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, ]
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, ]
write.table(upregulated$geneSymbol, file = "output/deseq2/upregulatedXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/deseq2/downregulatedXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# R session info
loaded_pkgs <- sessionInfo()$otherPkgs
pkg_versions <- data.frame(
  package = names(loaded_pkgs),
  version = sapply(loaded_pkgs, function(x) x$Version)
)
write.table(
  pkg_versions,
  file = "output/deseq2/RsessionInfo-DESEQ2_DEG.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

###### -------------- ######
```


**OUTPUT FILES:**
- Complete DESEQ2 output table: `output/deseq2/resXinclude_PSC_Hypo_vs_Norm-featurecounts_multi.txt`
- Volcano plot of DEGs *adjusted p-value < 0.05 and |log2FC| > 0.58* : `output/deseq2/plotVolcano_resXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.pdf`
- Volcano plot of DEGs with hypoxia marker genes labelled: `output/deseq2/plotVolcano_resXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-hypoxia_genes.pdf`
- Volcano plot of DEGs with ephrin genes labelled: `output/deseq2/plotVolcano_resXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-Ephrin_CellChat_genes.pdf`
- List of up regulated genes *adjusted p-value < 0.05 and |log2FC| > 0.58*: `output/deseq2/upregulatedXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.txt`
- List of down regulated genes *adjusted p-value < 0.05 and |log2FC| > 0.58*: `output/deseq2/downregulatedXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.txt`
- R packages version: `output/deseq2/RsessionInfo-DESEQ2_DEG.txt`





## Alternative splicing 


Identify, Annotate and Visualize Alternative Splicing and Isoform Switches with Functional Consequences using [IsoformSwitchAnalyzeR](https://github.com/kvittingseerup/IsoformSwitchAnalyzeR).

The ***IsoformSwitchAnalyzeR workflow consists of two main steps***:
- **isoformSwitchAnalysisPart1**: Identification of significant isoform switches, ORF prediction, and extraction of nucleotide and amino acid sequences for downstream analyses (generate files: `isoformSwitchAnalyzeR_isoform_nt.fasta` and `isoformSwitchAnalyzeR_isoform_AA_complete.fasta`)
- **isoformSwitchAnalysisPart2**: Functional annotation of isoform switches, alternative splicing analysis, and visualization of predicted functional consequences.


***IMPORTANT***: Running `isoformSwitchAnalysisPart2()` requires results from several external sequence analysis tools, which must be run after Part 1 and provided as input files:
- `isoformSwitchAnalyzeR_isoform_nt.fasta` is used for **[CPC2](https://cpc2.gao-lab.org/) (coding potential)**
- `isoformSwitchAnalyzeR_isoform_AA_complete.fasta` is used for **PFAM (protein domain), [IUPred2](https://iupred2a.elte.hu/) (intrinsically disordered region) and [SignalP 6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/) (signal peptide)**
  - Use web server for [CPC2](https://cpc2.gao-lab.org/), [IUPred2](https://iupred2a.elte.hu/) and [SignalP 6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/)
  - Perform analysis on the cluster for PFAM (see below)

PFAM analysis with **HMMER** (v1.0) (After running `isoformSwitchAnalysisPart1()`). The following command was run:
```bash
# Run PFAM scan
./pfam_scan.py <isoformSwitchAnalyzeR_isoform_nt.fasta> ../pfamdb/ -out <isoformSwitchAnalyzeR_isoform_nt.pfamresult>

# Reformat PFAM output for IsoformSwitchAnalyzeR using a custom script
python3 scripts/reformat_pfam_kallisto-PSC.py
```
--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/reformat_pfam_kallisto-PSC.py`


See below **representative example of complete analysis for PSC:**

```R
# Load packages
library("IsoformSwitchAnalyzeR") # v2.2.0
library("rtracklayer") # v1.62.0


set.seed(42) # set seed for reproducibility


# Import gene annotation file to convert gene ID to gene symbol
gtf <- import("../../Master/meta/gencode.v47.chr_patch_hapl_scaff.annotation.gtf")
gtf_df <- as.data.frame(gtf)
gene_df <- gtf_df %>%
  filter(type == "gene") %>%
  dplyr::select(seqnames, gene_name) %>%
  distinct() %>%
  rename(chromosome = seqnames) %>%
  as_tibble()

# Import isoform-level quantification data from kallisto
salmonQuant <- importIsoformExpression(
  parentDir = "output/kallisto/",  # Indicate kallisto output folder
  pattern = "/PSC_"                # Indicate sample pattern
)

# Process count for isoformSwitchAnalysisPart1
## Construct sample metadata 
myDesign <- data.frame(
  sampleID  = colnames(salmonQuant$abundance)[-1],
  condition = sub("^PSC_(Hypo|Norm)_.*", "\\1", colnames(salmonQuant$abundance)[-1])
)
myDesign$condition <- factor(myDesign$condition, levels = c("Norm", "Hypo"))
myDesign <- myDesign[order(myDesign$condition, myDesign$sampleID), ]
rownames(myDesign) <- myDesign$sampleID

## Create IsoformSwitchAnalyzeR input object
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

# Run IsoformSwitchAnalyzeR core analysis (part 1)
SwitchList <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist    = aSwitchList,
  pathToOutput          = 'output/IsoformSwitchAnalyzeR_kallisto/PSC',
  outputSequences       = TRUE,  # required for downstream sequence-based analyses
  prepareForWebServers  = TRUE   # enables compatibility with external web tools
)

#--> Recommended to SAVE image session or .rds object of `SwitchList`

#--> Run in WebServer the CPC2, PFAM, IUPRED2A, SIGNALP

# Run IsoformSwitchAnalyzeR core analysis (part 2)
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

# Data exploration
## Volcano plot showing number of significant isoform switch
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


## Consequence of the isoform switch (ie. ORF loss, transcript is non-coding, ...)
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/PSC/extractConsequenceEnrichment.pdf', onefile = FALSE, height=6, width = 8)
extractConsequenceEnrichment(
    analysSwitchList,
    consequencesToAnalyze='all',
    analysisOppositeConsequence = TRUE,
    localTheme = theme_bw(base_size = 14), # Increase font size in vignette
    returnResult = TRUE, 
)
dev.off()

## Enrichment of Alternative Splicing event (ie. Exon Skipping, Intron Retention, ...)
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/PSC/extractSplicingEnrichment.pdf', onefile = FALSE, height=4, width = 8)
extractSplicingEnrichment(
    analysSwitchList,
    returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)
dev.off()

## Display isoform switch of a gene of interest (ONLY WORK ON SIGNFICANT GENE!)
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/PSC/switchPlot-NormHypo-EPHA2.pdf', onefile = FALSE, height=6, width = 9)
switchPlot(analysSwitchList, gene= "EPHA3", condition1= "Norm", condition2= "Hypo", reverseMinus = FALSE)
dev.off()


###### --- Export --- ######

# Consequence of the isoform switch (ie. ORF loss, transcript is non-coding,...)
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

# Enrichment of Alternative Splicing event (ie. Exon Skipping, Intron Retention, ...)
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

# Genes with significant isofom switch with a functional consequence
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

# R session info
loaded_pkgs <- sessionInfo()$otherPkgs
pkg_versions <- data.frame(
  package = names(loaded_pkgs),
  version = sapply(loaded_pkgs, function(x) x$Version)
)
write.table(
  pkg_versions,
  file = "output/IsoformSwitchAnalyzeR_kallisto/RsessionInfo-IsoformSwitchAnalyzeR.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

###### -------------- ######
```

**OUTPUT FILES**:
- R packages version: `output/IsoformSwitchAnalyzeR_kallisto/RsessionInfo-IsoformSwitchAnalyzeR.txt`
- **PSC**
  - Output of `isoformSwitchAnalysisPart1()`:
    - AA sequences: `output/IsoformSwitchAnalyzeR_kallisto/PSC/isoformSwitchAnalyzer_isoform_AA_complete.fasta` (**and subset version**)
    - Nucleotide sequences: `output/IsoformSwitchAnalyzeR_kallisto/PSC/isoformSwitchAnalyzer_isoform_nt.fasta`
  - Volcano plot showing number of significant isoform switch: `output/IsoformSwitchAnalyzeR_kallisto/PSC/Overview_Plots.pdf`
  - Consequence of the isoform switch (ie. ORF loss, transcript is non-coding,...): 
    - plot: `output/IsoformSwitchAnalyzeR_kallisto/PSC/extractConsequenceEnrichment.pdf`
    - consequence for each gene: `output/IsoformSwitchAnalyzeR_kallisto/PSC/extractConsequenceEnrichment.txt`
  - Enrichment of Alternative Splicing event (ie. Exon Skipping, Intron Retention, ...): 
    - plot: `output/IsoformSwitchAnalyzeR_kallisto/PSC/extractSplicingEnrichment.pdf`
    - AS event for each gene: `output/IsoformSwitchAnalyzeR_kallisto/PSC/extractSplicingEnrichment.txt`
  - Genes with significant isofom switch with a functional consequence: `output/IsoformSwitchAnalyzeR_kallisto/PSC/significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-PSC.txt`
  - Plot of isoform switch of a gene of interest: `output/IsoformSwitchAnalyzeR_kallisto/PSC/switchPlot-NormHypo-[GENE OF INTEREST].pdf`
- **ReN**
  - Output of `isoformSwitchAnalysisPart1()`:
    - AA sequences: `output/IsoformSwitchAnalyzeR_kallisto/ReN/isoformSwitchAnalyzer_isoform_AA_complete.fasta` (**and subset version**)
    - Nucleotide sequences: `output/IsoformSwitchAnalyzeR_kallisto/ReN/isoformSwitchAnalyzer_isoform_nt.fasta`
  - Volcano plot showing number of significant isoform switch: `output/IsoformSwitchAnalyzeR_kallisto/ReN/Overview_Plots.pdf`
  - Consequence of the isoform switch (ie. ORF loss, transcript is non-coding,...): 
    - plot: `output/IsoformSwitchAnalyzeR_kallisto/ReN/extractConsequenceEnrichment.pdf`
    - consequence for each gene: `output/IsoformSwitchAnalyzeR_kallisto/ReN/extractConsequenceEnrichment.txt`
  - Enrichment of Alternative Splicing event (ie. Exon Skipping, Intron Retention, ...): 
    - plot: `output/IsoformSwitchAnalyzeR_kallisto/ReN/extractSplicingEnrichment.pdf`
    - AS event for each gene: `output/IsoformSwitchAnalyzeR_kallisto/ReN/extractSplicingEnrichment.txt`
  - Genes with significant isofom switch with a functional consequence: `output/IsoformSwitchAnalyzeR_kallisto/ReN/significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-ReN.txt`
  - Plot of isoform switch of a gene of interest: `output/IsoformSwitchAnalyzeR_kallisto/ReN/switchPlot-NormHypo-[GENE OF INTEREST].pdf`


## Functional enrichment 

Functional enrichment analyses (GO, KEGG) on differentially expressed genes and genes exhibiting significant isoform switching using the **clusterProfiler** package in **R** (v4.2.2).

List of genes to be analyzed:
- Up and Down -regulated genes in PSC *adjusted p-value < 0.05 and |log2FC| > 0.58*
- Up and Down -regulated genes in PSC *adjusted p-value < 1e-5 and |log2FC| > 1*
- Up and Down -regulated genes in ReN *adjusted p-value < 0.05 and |log2FC| > 0.58*
- Genes with significant isoform switching in PSC 
- Genes with significant isoform switching in PSC by consequence (ie. for genes leading to transcript becoming non-coding, IDR loss, loss of functional domain...)

--> One representative example with Up-regulated genes in PSC *adjusted p-value < 1e-5 and |log2FC| > 1* provided

```R
# Load packages
library("clusterProfiler") # v4.6.2
library("pathview") # v1.38.0
library("DOSE") # v3.24.2
library("org.Hs.eg.db") # v3.16.0
library("enrichplot") # v1.18.4
library("rtracklayer") # v1.58.0
library("tidyverse") # v2.0.0

set.seed(42) # set seed for reproducibility

# Import gene annotation file to convert gene ID to gene symbol
gtf_file <- "../../Master/meta/gencode.v47.annotation.gtf"
gtf_data <- import(gtf_file)
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()


#########################################################################
###### Example with PSC --> adjusted p-value < 1e-5 and log2FC > 1 ######
#########################################################################

# Import list of up-regulated genes in PSC
PSC_up = read_csv("output/deseq2/upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.txt", col_names = "gene_name")

# Run GO BP enrichment
ego <- enrichGO(gene = as.character(PSC_up$gene_name), 
                keyType = "SYMBOL",     
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)

## Plot top 20 term                
pdf("output/GO/dotplot_BP-upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf", width=7, height=7)
dotplot(ego, showCategory=20)
dev.off()
## Plot top 10 term                
pdf("output/GO/dotplot_BP-upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf", width=5, height=4)
dotplot(ego, showCategory=10)
dev.off()

###### --- Export --- ######
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
###### -------------- ######





# Run KEGG pathway enrichment
entrez_genes <- as.character( mapIds(org.Hs.eg.db, as.character(PSC_up$gene_name), 'ENTREZID', 'SYMBOL') )
ekegg <- enrichKEGG(gene = entrez_genes, 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05)           
## Plot top 20 term                
pdf("output/GO/dotplot_KEGG-upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf", width=7, height=7)
dotplot(ekegg, showCategory=20)
dev.off()
## Plot top 10 term                
pdf("output/GO/dotplot_KEGG-upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf", width=5, height=4)
dotplot(ekegg, showCategory=10)
dev.off()

###### --- Export --- ######
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
###### -------------- ######
```


**OUTPUT FILES**:
- **Differentially Expressed Genes (DEGs)** - **PSC**
  - **GO Biological Processes**
    - adjusted p-value < 0.05 and |log2FC| > 0.58
      - Top 20 GO terms for up-regulated genes: `output/GO/dotplot_BP-upregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf`
      - Top 10 GO terms for up-regulated genes:  `output/GO/dotplot_BP-upregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf`
      - Top 20 GO terms for down-regulated genes: `output/GO/dotplot_BP-downregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf`
      - Top 10 GO terms for down-regulated genes: `output/GO/dotplot_BP-downregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf`
      - Output table with GO term and associated genes for up-regulated genes: `output/GO/BP-upregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.tsv`
      - Output table with GO term and associated genes for down-regulated genes: `output/GO/BP-downregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.tsv`
    - adjusted p-value < 1e-5 and |log2FC| > 1
      - Top 20 GO terms for up-regulated genes: `output/GO/dotplot_BP-upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf`
      - Top 10 GO terms for up-regulated genes: `output/GO/dotplot_BP-upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf`
      - Top 20 GO terms for down-regulated genes: `output/GO/dotplot_BP-downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf`
      - Top 10 GO terms for down-regulated genes: `output/GO/dotplot_BP-downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf`
      - Output table with GO term and associated genes for up-regulated genes: `output/GO/BP-upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.tsv`
      - Output table with GO term and associated genes for down-regulated genes: `output/GO/BP-downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.tsv`
  - **KEGG Pathways**
    - adjusted p-value < 0.05 and |log2FC| > 0.58
      - Top 20 KEGG terms for up-regulated genes: `output/GO/dotplot_KEGG-upregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf`
      - Top 20 KEGG terms for up-regulated genes: `output/GO/dotplot_KEGG-upregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf`
      - Top 20 KEGG terms for down-regulated genes: `output/GO/dotplot_KEGG-downregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf`
      - Top 20 KEGG terms for down-regulated genes: `output/GO/dotplot_KEGG-downregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf`
      - Output table with KEGG term and associated genes for up-regulated genes: `output/GO/KEGG-upregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.tsv`
      - Output table with KEGG term and associated genes for down-regulated genes: `output/GO/KEGG-downregulatedresXinclude_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi.tsv`
    - adjusted p-value < 1e-5 and |log2FC| > 1
      - Top 20 KEGG terms for up-regulated genes: `output/GO/dotplot_KEGG-upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf`
      - Top 20 KEGG terms for up-regulated genes: `output/GO/dotplot_KEGG-upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf`
      - Top 20 KEGG terms for down-regulated genes: `output/GO/dotplot_KEGG-downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf`
      - Top 20 KEGG terms for down-regulated genes: `output/GO/dotplot_KEGG-downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf`
      - Output table with KEGG term and associated genes for up-regulated genes: `output/GO/KEGG-upregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.tsv`
      - Output table with KEGG term and associated genes for down-regulated genes: `output/GO/KEGG-downregulatedresXinclude_q00001fc1_PSC_Hypo_vs_Norm-featurecounts_multi.tsv`

- **Differentially Expressed Genes (DEGs)** - **ReN**
  - **GO Biological Processes**
    - adjusted p-value < 0.05 and |log2FC| > 0.58
      - Top 20 GO terms for up-regulated genes: `output/GO/dotplot_BP-upregulatedresXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-top20.pdf`
      - Top 10 GO terms for up-regulated genes:  `output/GO/dotplot_BP-upregulatedresXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-top10.pdf`
      - Top 20 GO terms for down-regulated genes: `output/GO/dotplot_BP-downregulatedresXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-top20.pdf`
      - Top 10 GO terms for down-regulated genes: `output/GO/dotplot_BP-downregulatedresXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-top10.pdf`
      - Output table with GO term and associated genes for up-regulated genes: `output/GO/BP-upregulatedresXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.tsv`
      - Output table with GO term and associated genes for down-regulated genes: `output/GO/BP-downregulatedresXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.tsv`
  - **KEGG Pathways**
    - adjusted p-value < 0.05 and |log2FC| > 0.58
      - Top 20 KEGG terms for up-regulated genes: `output/GO/dotplot_KEGG-upregulatedresXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-top20.pdf`
      - Top 10 KEGG terms for up-regulated genes:  `output/GO/dotplot_KEGG-upregulatedresXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-top10.pdf`
      - Top 20 KEGG terms for down-regulated genes: `output/GO/dotplot_KEGG-downregulatedresXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-top20.pdf`
      - Top 10 KEGG terms for down-regulated genes: `output/GO/dotplot_KEGG-downregulatedresXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi-top10.pdf`
      - Output table with KEGG term and associated genes for up-regulated genes: `output/GO/KEGG-upregulatedresXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.tsv`
      - Output table with KEGG term and associated genes for down-regulated genes: `output/GO/KEGG-downregulatedresXinclude_q05fc058_ReN_Hypo_vs_Norm-featurecounts_multi.tsv`

- **Genes with significant isoform switch** - **PSC**
  - ***All genes*** 
    - **GO Biological Processes**
      - Top 20 GO terms: `output/GO/dotplot_BP-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-PSC-top20.pdf`
      - Top 10 GO terms: `output/GO/dotplot_BP-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-PSC-top10.pdf`
      - Output table with GO term and associated genes: `output/GO/BP-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-PSC.tsv`
    - **KEGG**
      - Top 20 KEGG terms: `output/GO/dotplot_KEGG-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-PSC-top20.pdf`
      - Top 10 KEGG terms: `output/GO/dotplot_KEGG-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-PSC-top10.pdf`
      - Output table with KEGG term and associated genes: `output/GO/KEGG-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-PSC.tsv`
  - ***Functional enrichment per isoform switch consequence:***
    - **GO Biological Processes**
      - Top 10 GO terms; **one page for each category**: `output/GO/dotplot_BP-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-extractConsequenceEnrichment-PSC-top10.pdf`
      - Output table with GO term and associated genes: `BP-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-extractConsequenceEnrichment-[CONSEQUENCE ID]-PSC.tsv`
    - **KEGG**
      - Top 10 KEGG terms; **one page for each category**: `output/GO/KEGG-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-extractConsequenceEnrichment-PSC-top10.pdf`
      - Output table with KEGG term and associated genes: `KEGG-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-extractConsequenceEnrichment-[CONSEQUENCE ID]-PSC.tsv`



- **Genes with significant isoform switch** - **ReN**
  - ***All genes*** 
    - **GO Biological Processes** --> Not signif.
    - **KEGG**
      - Top 10 KEGG terms: `output/GO/dotplot_KEGG-significant_isoforms_dIF01qval05switchConsequencesGeneTRUE_geneSymbol-ReN-top10.pdf`
  - ***Functional enrichment per isoform switch consequence:*** --> NA: no significant consequence enrichment in ReN

