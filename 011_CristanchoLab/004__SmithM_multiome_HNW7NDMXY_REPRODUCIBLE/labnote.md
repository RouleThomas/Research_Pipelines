--> This folder is the **clean + reproducible** version of: `001_CristanchoLab/001__BDRhapsody_v1` (original exploration)

# Data Overview

- **Data type**: 10X multiome
- **Factors**:
  - **Organism**: Mice cortex
  - **Condition**: Normoxia, Hypoxia
  - **Replicates**: 10 biological replicates per group (total n = 20)



# Data processing

## Data import

Analysis started with the counts (output of `cellranger-arc count`) at `\cristancho_data\rawdata\CristanchoA_SmithM_multiome_HNW7NDMXY_240226_H3VKGDSXC_240227\CristanchoA_SmithM_multiome_HNW7NDMXY_240226_H3VKGDSXC_240227\count-outs`



## Data analysis




### Doublet detection - RNA assay 

Identify doublet/multiplet cells with [scrublet](https://github.com/swolock/scrublet). With data being multiome I made a custom script for scrublet to work (ie. as design for non-multiome dataset).

```bash
# Command example (ie. repeated for each sample individually)
python3 scripts/scrublet_doublets_histogram-10xMultiome.py AnaCristancho-20240206-1957-H-M_Multiome/outs/filtered_feature_bc_matrix output/doublets_10xMultiome/AnaCristancho-20240206-1957-H-M_Multiome.tsv 0.25
```



### ambiant RNA contamination correction - RNA assay 

Let's correct our RNA count matrix for ambient RNA contamination with *soupX*.

Below code ran individually for each sample:
```R
# install.packages('SoupX')
library("SoupX")
library("Seurat")
library("tidyverse")

# soupX decontamination
## Decontaminate one channel of 10X data mapped with cellranger
sc = load10X('AnaCristancho-20240213-2264-H-F_Multiome/outs')   # CHANGE FILE NAME HERE

## Assess % of conta
pdf("output/soupX/autoEstCont-AnaCristancho-20240213-2264-H-F_Multiome.pdf", width=10, height=10)   # CHANGE FILE NAME HEREs
sc = autoEstCont(sc) 
dev.off()
## Generate the corrected matrix
out = adjustCounts(sc)
## Save the matrix
save(out, file = "output/soupX/AnaCristancho-20240213-2264-H-F_Multiome.RData") # CHANGE FILE NAME HERE
```


### QC filtering - RNA assay 

QC plots (ie. nCount_RNA, nFeature_RNA...) generated for each sample individually; outlier cells and doublet removed for each sample: specific QC treshold for each sample.


```R
set.seed(42)

# packages
library("SoupX")
library("Seurat")
library("tidyverse")
library("dplyr")
library("Seurat")
library("patchwork")
library("sctransform")
library("glmGamPoi")
library("celldex")
library("SingleR")
library("gprofiler2")


################################################
############# Import RNA ########################
################################################

## Load the matrix and Create SEURAT object
samples <- list(
  AnaCristancho_20240206_1957_H_M = "output/soupX/AnaCristancho-20240206-1957-H-M_Multiome.RData",
  AnaCristancho_20240206_1966_N_M = "output/soupX/AnaCristancho-20240206-1966-N-M_Multiome.RData",
  AnaCristancho_20240206_1987_H_F = "output/soupX/AnaCristancho-20240206-1987-H-F_Multiome.RData",
  AnaCristancho_20240206_2027_N_F = "output/soupX/AnaCristancho-20240206-2027-N-F_Multiome.RData",
  
  AnaCristancho_20240207_1966_N_F = "output/soupX/AnaCristancho-20240207-1966-N-F_Multiome.RData",
  AnaCristancho_20240207_2011_H_M = "output/soupX/AnaCristancho-20240207-2011-H-M_Multiome.RData",
  AnaCristancho_20240207_2023_H_F = "output/soupX/AnaCristancho-20240207-2023-H-F_Multiome.RData",
  AnaCristancho_20240207_2028_N_M = "output/soupX/AnaCristancho-20240207-2028-N-M_Multiome.RData",
  
  AnaCristancho_20240212_1965_N_F = "output/soupX/AnaCristancho-20240212-1965-N-F_Multiome.RData",
  AnaCristancho_20240212_2009_H_M = "output/soupX/AnaCristancho-20240212-2009-H-M_Multiome.RData",
  AnaCristancho_20240212_2012_H_F = "output/soupX/AnaCristancho-20240212-2012-H-F_Multiome.RData",
  AnaCristancho_20240212_2027_N_M = "output/soupX/AnaCristancho-20240212-2027-N-M_Multiome.RData",
  
  AnaCristancho_20240213_2022_H_M = "output/soupX/AnaCristancho-20240213-2022-H-M_Multiome.RData",
  AnaCristancho_20240213_2113_N_M = "output/soupX/AnaCristancho-20240213-2113-N-M_Multiome.RData",
  AnaCristancho_20240213_2264_H_F = "output/soupX/AnaCristancho-20240213-2264-H-F_Multiome.RData",
  
  AnaCristancho_20240214_1965_N_M = "output/soupX/AnaCristancho-20240214-1965-N-M_Multiome.RData",
  AnaCristancho_20240214_2023_H_M = "output/soupX/AnaCristancho-20240214-2023-H-M_Multiome.RData",
  AnaCristancho_20240214_2028_N_F = "output/soupX/AnaCristancho-20240214-2028-N-F_Multiome.RData",
  AnaCristancho_20240214_2113_N_F = "output/soupX/AnaCristancho-20240214-2113-N-F_Multiome.RData",
  AnaCristancho_20240214_2139_H_F = "output/soupX/AnaCristancho-20240214-2139-H-F_Multiome.RData"
)


seurat_objects <- list()

for (sample_name in names(samples)) {
  load(samples[[sample_name]])
  seurat_objects[[sample_name]] <- CreateSeuratObject(counts = out, project = sample_name)
}

## Function to assign Seurat objects to variables (unlist the list)
assign_seurat_objects <- function(seurat_objects_list) {
  for (sample_name in names(seurat_objects_list)) {
    assign(sample_name, seurat_objects_list[[sample_name]], envir = .GlobalEnv)
  }
}
assign_seurat_objects(seurat_objects) # This apply the function

# Add mitochondrial and ribosomal content
add_quality_control <- function(seurat_object) {
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
  seurat_object[["percent.rb"]] <- PercentageFeatureSet(seurat_object, pattern = "^Rp[sl]")
  return(seurat_object)
}
seurat_objects <- lapply(seurat_objects, add_quality_control) 
assign_seurat_objects(seurat_objects) 


# Add doublet information
add_doublet_information <- function(sample_name, seurat_object) {
  # Convert R-safe sample name back to file name format
  file_sample_name <- gsub("_", "-", sample_name)
  doublet_file <- paste0("output/doublets_10xMultiome/", file_sample_name, "_Multiome.tsv")
  doublets <- read.table(doublet_file, header = FALSE, row.names = 1)
  colnames(doublets) <- c("Doublet_score", "Is_doublet")
  seurat_object <- AddMetaData(seurat_object, doublets)
  seurat_object$Doublet_score <- as.numeric(seurat_object$Doublet_score)
  return(seurat_object)
}
for (sample_name in names(seurat_objects)) {
  seurat_objects[[sample_name]] <- add_doublet_information(
    sample_name,
    seurat_objects[[sample_name]]
  )
}
assign_seurat_objects(seurat_objects)


###################################
# AnaCristancho_20240206_1957_H_M #
###################################

## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240206_1957_H_M-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240206_1957_H_M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240206_1957_H_M-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240206_1957_H_M, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240206_1957_H_M-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240206_1957_H_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,20000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240206_1957_H_M-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240206_1957_H_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,2500)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 250) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 6000) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 400) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 25000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 2) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 3) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240206_1957_H_M <- apply_qc(AnaCristancho_20240206_1957_H_M)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240206_1957_H_M$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240206_1957_H_M"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240206_1957_H_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240206_1957_H_M <- subset(
  AnaCristancho_20240206_1957_H_M,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240206_1957_H_M) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240206_1957_H_M <- NormalizeData(
  AnaCristancho_20240206_1957_H_M,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240206_1957_H_M)
AnaCristancho_20240206_1957_H_M <- ScaleData(
  AnaCristancho_20240206_1957_H_M,
  features = all.genes
)
AnaCristancho_20240206_1957_H_M <- CellCycleScoring(
  AnaCristancho_20240206_1957_H_M,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240206_1957_H_M$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240206_1957_H_M"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240206_1957_H_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240206_1957_H_M,
  file = "output/seurat/AnaCristancho_20240206_1957_H_M-QCPass_vKeepBottom.rds"
)
#################################################################



###################################
# AnaCristancho_20240206_1966_N_M #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240206_1966_N_M-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240206_1966_N_M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240206_1966_N_M-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240206_1966_N_M, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240206_1966_N_M-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240206_1966_N_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,2500)
dev.off()


##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 300) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 7500) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 500) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 40000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 1.5) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 3) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240206_1966_N_M <- apply_qc(AnaCristancho_20240206_1966_N_M)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240206_1966_N_M$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240206_1966_N_M"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240206_1966_N_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240206_1966_N_M <- subset(
  AnaCristancho_20240206_1966_N_M,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240206_1966_N_M) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240206_1966_N_M <- NormalizeData(
  AnaCristancho_20240206_1966_N_M,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240206_1966_N_M)
AnaCristancho_20240206_1966_N_M <- ScaleData(
  AnaCristancho_20240206_1966_N_M,
  features = all.genes
)
AnaCristancho_20240206_1966_N_M <- CellCycleScoring(
  AnaCristancho_20240206_1966_N_M,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240206_1966_N_M$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240206_1966_N_M"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240206_1966_N_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240206_1966_N_M,
  file = "output/seurat/AnaCristancho_20240206_1966_N_M-QCPass_vKeepBottom.rds"
)
#################################################################


###################################
# AnaCristancho_20240206_1987_H_F #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240206_1987_H_F-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240206_1987_H_F, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240206_1987_H_F-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240206_1987_H_F, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240206_1987_H_F-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240206_1987_H_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 250) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 5000) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 200) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 25000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 1.5) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 4) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240206_1987_H_F <- apply_qc(AnaCristancho_20240206_1987_H_F)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240206_1987_H_F$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240206_1987_H_F"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240206_1987_H_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240206_1987_H_F <- subset(
  AnaCristancho_20240206_1987_H_F,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240206_1987_H_F) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240206_1987_H_F <- NormalizeData(
  AnaCristancho_20240206_1987_H_F,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240206_1987_H_F)
AnaCristancho_20240206_1987_H_F <- ScaleData(
  AnaCristancho_20240206_1987_H_F,
  features = all.genes
)
AnaCristancho_20240206_1987_H_F <- CellCycleScoring(
  AnaCristancho_20240206_1987_H_F,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240206_1987_H_F$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240206_1987_H_F"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240206_1987_H_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240206_1987_H_F,
  file = "output/seurat/AnaCristancho_20240206_1987_H_F-QCPass_vKeepBottom.rds"
)
#################################################################



###################################
# AnaCristancho_20240206_2027_N_F #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240206_2027_N_F-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240206_2027_N_F, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240206_2027_N_F-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240206_2027_N_F, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240206_2027_N_F-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240206_2027_N_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,50000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 200) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 5000) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 200) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 20000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 2) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 3) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240206_2027_N_F <- apply_qc(AnaCristancho_20240206_2027_N_F)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240206_2027_N_F$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240206_2027_N_F"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240206_2027_N_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240206_2027_N_F <- subset(
  AnaCristancho_20240206_2027_N_F,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240206_2027_N_F) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240206_2027_N_F <- NormalizeData(
  AnaCristancho_20240206_2027_N_F,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240206_2027_N_F)
AnaCristancho_20240206_2027_N_F <- ScaleData(
  AnaCristancho_20240206_2027_N_F,
  features = all.genes
)
AnaCristancho_20240206_2027_N_F <- CellCycleScoring(
  AnaCristancho_20240206_2027_N_F,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240206_2027_N_F$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240206_2027_N_F"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240206_2027_N_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240206_2027_N_F,
  file = "output/seurat/AnaCristancho_20240206_2027_N_F-QCPass_vKeepBottom.rds"
)
#################################################################



###################################
# AnaCristancho_20240207_1966_N_F #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240207_1966_N_F-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240207_1966_N_F, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240207_1966_N_F-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240207_1966_N_F, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240207_1966_N_F-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240207_1966_N_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,25000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240207_1966_N_F-nCount_RNA1.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240207_1966_N_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 250) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 6000) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 250) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 20000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 2.5) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 4) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240207_1966_N_F <- apply_qc(AnaCristancho_20240207_1966_N_F)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240207_1966_N_F$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240207_1966_N_F"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240207_1966_N_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240207_1966_N_F <- subset(
  AnaCristancho_20240207_1966_N_F,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240207_1966_N_F) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240207_1966_N_F <- NormalizeData(
  AnaCristancho_20240207_1966_N_F,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240207_1966_N_F)
AnaCristancho_20240207_1966_N_F <- ScaleData(
  AnaCristancho_20240207_1966_N_F,
  features = all.genes
)
AnaCristancho_20240207_1966_N_F <- CellCycleScoring(
  AnaCristancho_20240207_1966_N_F,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240207_1966_N_F$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240207_1966_N_F"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240207_1966_N_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240207_1966_N_F,
  file = "output/seurat/AnaCristancho_20240207_1966_N_F-QCPass_vKeepBottom.rds"
)
#################################################################



###################################
# AnaCristancho_20240207_2011_H_M #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240207_2011_H_M-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240207_2011_H_M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240207_2011_H_M-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240207_2011_H_M, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240207_2011_H_M-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240207_2011_H_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,25000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240207_2011_H_M-nCount_RNA1.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240207_2011_H_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 250) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 7000) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 250) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 30000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 2) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 4) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240207_2011_H_M <- apply_qc(AnaCristancho_20240207_2011_H_M)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240207_2011_H_M$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240207_2011_H_M"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240207_2011_H_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240207_2011_H_M <- subset(
  AnaCristancho_20240207_2011_H_M,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240207_2011_H_M) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240207_2011_H_M <- NormalizeData(
  AnaCristancho_20240207_2011_H_M,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240207_2011_H_M)
AnaCristancho_20240207_2011_H_M <- ScaleData(
  AnaCristancho_20240207_2011_H_M,
  features = all.genes
)
AnaCristancho_20240207_2011_H_M <- CellCycleScoring(
  AnaCristancho_20240207_2011_H_M,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240207_2011_H_M$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240207_2011_H_M"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240207_2011_H_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240207_2011_H_M,
  file = "output/seurat/AnaCristancho_20240207_2011_H_M-QCPass_vKeepBottom.rds"
)
#################################################################




###################################
# AnaCristancho_20240207_2023_H_F #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240207_2023_H_F-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240207_2023_H_F, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240207_2023_H_F-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240207_2023_H_F, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240207_2023_H_F-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240207_2023_H_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,25000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240207_2023_H_F-nCount_RNA1.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240207_2023_H_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 200) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 6000) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 200) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 20000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 2) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 4) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240207_2023_H_F <- apply_qc(AnaCristancho_20240207_2023_H_F)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240207_2023_H_F$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240207_2023_H_F"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240207_2023_H_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240207_2023_H_F <- subset(
  AnaCristancho_20240207_2023_H_F,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240207_2023_H_F) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240207_2023_H_F <- NormalizeData(
  AnaCristancho_20240207_2023_H_F,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240207_2023_H_F)
AnaCristancho_20240207_2023_H_F <- ScaleData(
  AnaCristancho_20240207_2023_H_F,
  features = all.genes
)
AnaCristancho_20240207_2023_H_F <- CellCycleScoring(
  AnaCristancho_20240207_2023_H_F,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240207_2023_H_F$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240207_2023_H_F"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240207_2023_H_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240207_2023_H_F,
  file = "output/seurat/AnaCristancho_20240207_2023_H_F-QCPass_vKeepBottom.rds"
)
#################################################################


###################################
# AnaCristancho_20240207_2028_N_M #
###################################

## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240207_2028_N_M-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240207_2028_N_M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240207_2028_N_M-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240207_2028_N_M, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240207_2028_N_M-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240207_2028_N_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,25000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240207_2028_N_M-nCount_RNA1.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240207_2028_N_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 250) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 6000) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 250) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 25000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 3) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 4) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240207_2028_N_M <- apply_qc(AnaCristancho_20240207_2028_N_M)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240207_2028_N_M$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240207_2028_N_M"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240207_2028_N_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240207_2028_N_M <- subset(
  AnaCristancho_20240207_2028_N_M,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240207_2028_N_M) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240207_2028_N_M <- NormalizeData(
  AnaCristancho_20240207_2028_N_M,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240207_2028_N_M)
AnaCristancho_20240207_2028_N_M <- ScaleData(
  AnaCristancho_20240207_2028_N_M,
  features = all.genes
)
AnaCristancho_20240207_2028_N_M <- CellCycleScoring(
  AnaCristancho_20240207_2028_N_M,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240207_2028_N_M$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240207_2028_N_M"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240207_2028_N_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240207_2028_N_M,
  file = "output/seurat/AnaCristancho_20240207_2028_N_M-QCPass_vKeepBottom.rds"
)
#################################################################




###################################
# AnaCristancho_20240212_1965_N_F #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240212_1965_N_F-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240212_1965_N_F, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240212_1965_N_F-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240212_1965_N_F, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240212_1965_N_F-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240212_1965_N_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,25000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240212_1965_N_F-nCount_RNA1.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240212_1965_N_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 200) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 4000) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 250) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 12000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 2) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 5) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240212_1965_N_F <- apply_qc(AnaCristancho_20240212_1965_N_F)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240212_1965_N_F$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240212_1965_N_F"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240212_1965_N_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240212_1965_N_F <- subset(
  AnaCristancho_20240212_1965_N_F,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240212_1965_N_F) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240212_1965_N_F <- NormalizeData(
  AnaCristancho_20240212_1965_N_F,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240212_1965_N_F)
AnaCristancho_20240212_1965_N_F <- ScaleData(
  AnaCristancho_20240212_1965_N_F,
  features = all.genes
)
AnaCristancho_20240212_1965_N_F <- CellCycleScoring(
  AnaCristancho_20240212_1965_N_F,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240212_1965_N_F$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240212_1965_N_F"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240212_1965_N_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240212_1965_N_F,
  file = "output/seurat/AnaCristancho_20240212_1965_N_F-QCPass_vKeepBottom.rds"
)
#################################################################



###################################
# AnaCristancho_20240212_2009_H_M #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240212_2009_H_M-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240212_2009_H_M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240212_2009_H_M-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240212_2009_H_M, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240212_2009_H_M-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240212_2009_H_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,25000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240212_2009_H_M-nCount_RNA1.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240212_2009_H_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 250) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 5000) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 250) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 30000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 2) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 5) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240212_2009_H_M <- apply_qc(AnaCristancho_20240212_2009_H_M)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240212_2009_H_M$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240212_2009_H_M"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240212_2009_H_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240212_2009_H_M <- subset(
  AnaCristancho_20240212_2009_H_M,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240212_2009_H_M) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240212_2009_H_M <- NormalizeData(
  AnaCristancho_20240212_2009_H_M,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240212_2009_H_M)
AnaCristancho_20240212_2009_H_M <- ScaleData(
  AnaCristancho_20240212_2009_H_M,
  features = all.genes
)
AnaCristancho_20240212_2009_H_M <- CellCycleScoring(
  AnaCristancho_20240212_2009_H_M,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240212_2009_H_M$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240212_2009_H_M"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240212_2009_H_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240212_2009_H_M,
  file = "output/seurat/AnaCristancho_20240212_2009_H_M-QCPass_vKeepBottom.rds"
)
#################################################################



###################################
# AnaCristancho_20240212_2012_H_F #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240212_2012_H_F-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240212_2012_H_F, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240212_2012_H_F-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240212_2012_H_F, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240212_2012_H_F-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240212_2012_H_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,25000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240212_2012_H_F-nCount_RNA1.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240212_2012_H_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 250) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 7000) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 200) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 40000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 2) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 4) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240212_2012_H_F <- apply_qc(AnaCristancho_20240212_2012_H_F)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240212_2012_H_F$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240212_2012_H_F"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240212_2012_H_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240212_2012_H_F <- subset(
  AnaCristancho_20240212_2012_H_F,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240212_2012_H_F) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240212_2012_H_F <- NormalizeData(
  AnaCristancho_20240212_2012_H_F,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240212_2012_H_F)
AnaCristancho_20240212_2012_H_F <- ScaleData(
  AnaCristancho_20240212_2012_H_F,
  features = all.genes
)
AnaCristancho_20240212_2012_H_F <- CellCycleScoring(
  AnaCristancho_20240212_2012_H_F,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240212_2012_H_F$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240212_2012_H_F"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240212_2012_H_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240212_2012_H_F,
  file = "output/seurat/AnaCristancho_20240212_2012_H_F-QCPass_vKeepBottom.rds"
)
#################################################################



###################################
# AnaCristancho_20240212_2027_N_M #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240212_2027_N_M-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240212_2027_N_M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240212_2027_N_M-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240212_2027_N_M, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240212_2027_N_M-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240212_2027_N_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,25000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240212_2027_N_M-nCount_RNA1.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240212_2027_N_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 200) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 7500) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 200) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 25000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 4) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 4) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240212_2027_N_M <- apply_qc(AnaCristancho_20240212_2027_N_M)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240212_2027_N_M$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240212_2027_N_M"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240212_2027_N_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240212_2027_N_M <- subset(
  AnaCristancho_20240212_2027_N_M,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240212_2027_N_M) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240212_2027_N_M <- NormalizeData(
  AnaCristancho_20240212_2027_N_M,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240212_2027_N_M)
AnaCristancho_20240212_2027_N_M <- ScaleData(
  AnaCristancho_20240212_2027_N_M,
  features = all.genes
)
AnaCristancho_20240212_2027_N_M <- CellCycleScoring(
  AnaCristancho_20240212_2027_N_M,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240212_2027_N_M$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240212_2027_N_M"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240212_2027_N_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240212_2027_N_M,
  file = "output/seurat/AnaCristancho_20240212_2027_N_M-QCPass_vKeepBottom.rds"
)
#################################################################



###################################
# AnaCristancho_20240213_2022_H_M #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240213_2022_H_M-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240213_2022_H_M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240213_2022_H_M-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240213_2022_H_M, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240213_2022_H_M-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240213_2022_H_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,25000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240213_2022_H_M-nCount_RNA1.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240213_2022_H_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 200) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 6000) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 1600) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 25000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 2) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 4) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240213_2022_H_M <- apply_qc(AnaCristancho_20240213_2022_H_M)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240213_2022_H_M$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240213_2022_H_M"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240213_2022_H_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240213_2022_H_M <- subset(
  AnaCristancho_20240213_2022_H_M,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240213_2022_H_M) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240213_2022_H_M <- NormalizeData(
  AnaCristancho_20240213_2022_H_M,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240213_2022_H_M)
AnaCristancho_20240213_2022_H_M <- ScaleData(
  AnaCristancho_20240213_2022_H_M,
  features = all.genes
)
AnaCristancho_20240213_2022_H_M <- CellCycleScoring(
  AnaCristancho_20240213_2022_H_M,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240213_2022_H_M$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240213_2022_H_M"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240213_2022_H_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240213_2022_H_M,
  file = "output/seurat/AnaCristancho_20240213_2022_H_M-QCPass_vKeepBottom.rds"
)
#################################################################



###################################
# AnaCristancho_20240213_2113_N_M #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240213_2113_N_M-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240213_2113_N_M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240213_2113_N_M-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240213_2113_N_M, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240213_2113_N_M-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240213_2113_N_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,25000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240213_2113_N_M-nCount_RNA1.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240213_2113_N_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 250) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 7500) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 300) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 50000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 2) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 5) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240213_2113_N_M <- apply_qc(AnaCristancho_20240213_2113_N_M)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240213_2113_N_M$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240213_2113_N_M"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240213_2113_N_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240213_2113_N_M <- subset(
  AnaCristancho_20240213_2113_N_M,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240213_2113_N_M) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240213_2113_N_M <- NormalizeData(
  AnaCristancho_20240213_2113_N_M,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240213_2113_N_M)
AnaCristancho_20240213_2113_N_M <- ScaleData(
  AnaCristancho_20240213_2113_N_M,
  features = all.genes
)
AnaCristancho_20240213_2113_N_M <- CellCycleScoring(
  AnaCristancho_20240213_2113_N_M,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240213_2113_N_M$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240213_2113_N_M"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240213_2113_N_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240213_2113_N_M,
  file = "output/seurat/AnaCristancho_20240213_2113_N_M-QCPass_vKeepBottom.rds"
)
#################################################################


###################################
# AnaCristancho_20240213_2264_H_F #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240213_2264_H_F-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240213_2264_H_F, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240213_2264_H_F-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240213_2264_H_F, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240213_2264_H_F-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240213_2264_H_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,25000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240213_2264_H_F-nCount_RNA1.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240213_2264_H_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 250) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 6000) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 250) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 20000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 2) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 5) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240213_2264_H_F <- apply_qc(AnaCristancho_20240213_2264_H_F)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240213_2264_H_F$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240213_2264_H_F"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240213_2264_H_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240213_2264_H_F <- subset(
  AnaCristancho_20240213_2264_H_F,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240213_2264_H_F) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240213_2264_H_F <- NormalizeData(
  AnaCristancho_20240213_2264_H_F,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240213_2264_H_F)
AnaCristancho_20240213_2264_H_F <- ScaleData(
  AnaCristancho_20240213_2264_H_F,
  features = all.genes
)
AnaCristancho_20240213_2264_H_F <- CellCycleScoring(
  AnaCristancho_20240213_2264_H_F,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240213_2264_H_F$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240213_2264_H_F"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240213_2264_H_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240213_2264_H_F,
  file = "output/seurat/AnaCristancho_20240213_2264_H_F-QCPass_vKeepBottom.rds"
)
#################################################################



###################################
# AnaCristancho_20240214_1965_N_M #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_1965_N_M-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_1965_N_M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_1965_N_M-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_1965_N_M, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_1965_N_M-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_1965_N_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,50000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_1965_N_M-nCount_RNA1.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_1965_N_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 250) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 7500) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 250) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 50000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 2) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 4) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240214_1965_N_M <- apply_qc(AnaCristancho_20240214_1965_N_M)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240214_1965_N_M$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240214_1965_N_M"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240214_1965_N_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240214_1965_N_M <- subset(
  AnaCristancho_20240214_1965_N_M,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240214_1965_N_M) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240214_1965_N_M <- NormalizeData(
  AnaCristancho_20240214_1965_N_M,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240214_1965_N_M)
AnaCristancho_20240214_1965_N_M <- ScaleData(
  AnaCristancho_20240214_1965_N_M,
  features = all.genes
)
AnaCristancho_20240214_1965_N_M <- CellCycleScoring(
  AnaCristancho_20240214_1965_N_M,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240214_1965_N_M$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240214_1965_N_M"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240214_1965_N_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240214_1965_N_M,
  file = "output/seurat/AnaCristancho_20240214_1965_N_M-QCPass_vKeepBottom.rds"
)
#################################################################



###################################
# AnaCristancho_20240214_2023_H_M #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_2023_H_M-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_2023_H_M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_2023_H_M-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_2023_H_M, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_2023_H_M-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_2023_H_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,50000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_2023_H_M-nCount_RNA1.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_2023_H_M, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 250) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 5000) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 300) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 20000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 3) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 3.5) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240214_2023_H_M <- apply_qc(AnaCristancho_20240214_2023_H_M)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240214_2023_H_M$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240214_2023_H_M"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240214_2023_H_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240214_2023_H_M <- subset(
  AnaCristancho_20240214_2023_H_M,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240214_2023_H_M) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240214_2023_H_M <- NormalizeData(
  AnaCristancho_20240214_2023_H_M,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240214_2023_H_M)
AnaCristancho_20240214_2023_H_M <- ScaleData(
  AnaCristancho_20240214_2023_H_M,
  features = all.genes
)
AnaCristancho_20240214_2023_H_M <- CellCycleScoring(
  AnaCristancho_20240214_2023_H_M,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240214_2023_H_M$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240214_2023_H_M"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240214_2023_H_M.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240214_2023_H_M,
  file = "output/seurat/AnaCristancho_20240214_2023_H_M-QCPass_vKeepBottom.rds"
)
#################################################################




###################################
# AnaCristancho_20240214_2028_N_F #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_2028_N_F-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_2028_N_F, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_2028_N_F-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_2028_N_F, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_2028_N_F-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_2028_N_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,50000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_2028_N_F-nCount_RNA1.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_2028_N_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 300) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 7500) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 300) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 50000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 1) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 5) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240214_2028_N_F <- apply_qc(AnaCristancho_20240214_2028_N_F)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240214_2028_N_F$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240214_2028_N_F"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240214_2028_N_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240214_2028_N_F <- subset(
  AnaCristancho_20240214_2028_N_F,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240214_2028_N_F) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240214_2028_N_F <- NormalizeData(
  AnaCristancho_20240214_2028_N_F,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240214_2028_N_F)
AnaCristancho_20240214_2028_N_F <- ScaleData(
  AnaCristancho_20240214_2028_N_F,
  features = all.genes
)
AnaCristancho_20240214_2028_N_F <- CellCycleScoring(
  AnaCristancho_20240214_2028_N_F,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240214_2028_N_F$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240214_2028_N_F"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240214_2028_N_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240214_2028_N_F,
  file = "output/seurat/AnaCristancho_20240214_2028_N_F-QCPass_vKeepBottom.rds"
)
#################################################################



###################################
# AnaCristancho_20240214_2113_N_F #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_2113_N_F-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_2113_N_F, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_2113_N_F-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_2113_N_F, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_2113_N_F-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_2113_N_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,50000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_2113_N_F-nCount_RNA1.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_2113_N_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 300) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 7500) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 400) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 45000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 1) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 6) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240214_2113_N_F <- apply_qc(AnaCristancho_20240214_2113_N_F)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240214_2113_N_F$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240214_2113_N_F"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240214_2113_N_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240214_2113_N_F <- subset(
  AnaCristancho_20240214_2113_N_F,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240214_2113_N_F) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240214_2113_N_F <- NormalizeData(
  AnaCristancho_20240214_2113_N_F,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240214_2113_N_F)
AnaCristancho_20240214_2113_N_F <- ScaleData(
  AnaCristancho_20240214_2113_N_F,
  features = all.genes
)
AnaCristancho_20240214_2113_N_F <- CellCycleScoring(
  AnaCristancho_20240214_2113_N_F,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240214_2113_N_F$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240214_2113_N_F"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240214_2113_N_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240214_2113_N_F,
  file = "output/seurat/AnaCristancho_20240214_2113_N_F-QCPass_vKeepBottom.rds"
)
#################################################################


###################################
# AnaCristancho_20240214_2139_H_F #
###################################
## QC plot
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_2139_H_F-all.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_2139_H_F, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_2139_H_F-nFeature_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_2139_H_F, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,4000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_2139_H_F-nCount_RNA.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_2139_H_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,50000)
dev.off()
pdf("output/seurat/VlnPlot_QC-AnaCristancho_20240214_2139_H_F-nCount_RNA1.pdf", width = 15, height = 6)
VlnPlot(AnaCristancho_20240214_2139_H_F, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()

##### QC filtering keeping multiple QC reasons
apply_qc <- function(seurat_object) {
  meta <- seurat_object@meta.data
  qc_reasons <- rep("Pass", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    reasons <- c()
    if (meta$Is_doublet[i] == "True") reasons <- c(reasons, "Doublet")
    if (meta$nFeature_RNA[i] < 200) reasons <- c(reasons, "Low_nFeature")
    if (meta$nFeature_RNA[i] > 5000) reasons <- c(reasons, "High_nFeatureRNA")
    if (meta$nCount_RNA[i] < 200) reasons <- c(reasons, "Low_nCountRNA")
    if (meta$nCount_RNA[i] > 20000) reasons <- c(reasons, "High_nCountRNA")
    if (meta$percent.mt[i] > 2.5) reasons <- c(reasons, "High_MT")
    if (meta$percent.rb[i] > 3.5) reasons <- c(reasons, "High_RB")
    if (length(reasons) > 0) {
      qc_reasons[i] <- paste(reasons, collapse = ",")
    }
  }
  meta$QC <- qc_reasons
  seurat_object@meta.data <- meta
  return(seurat_object)
}
AnaCristancho_20240214_2139_H_F <- apply_qc(AnaCristancho_20240214_2139_H_F)
#### Write QC summary for one sample
qc_summary <- table(AnaCristancho_20240214_2139_H_F$QC)
qc_summary_df <- as.data.frame(qc_summary)
qc_summary_df$Sample <- "AnaCristancho_20240214_2139_H_F"
write.table(
  qc_summary_df,
  file = "output/seurat/QC_summary_vKeepBottom-AnaCristancho_20240214_2139_H_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
## subset Seurat object to keep cells that pass QC
AnaCristancho_20240214_2139_H_F <- subset(
  AnaCristancho_20240214_2139_H_F,
  subset = QC == "Pass"
)
# Normalize, scale data, and run cell cycle scoring
DefaultAssay(AnaCristancho_20240214_2139_H_F) <- "RNA"
set.seed(42)
## Load gene markers of cell cycle
mmus_s <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
mmus_g2m <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name
# Process one Seurat object
AnaCristancho_20240214_2139_H_F <- NormalizeData(
  AnaCristancho_20240214_2139_H_F,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
all.genes <- rownames(AnaCristancho_20240214_2139_H_F)
AnaCristancho_20240214_2139_H_F <- ScaleData(
  AnaCristancho_20240214_2139_H_F,
  features = all.genes
)
AnaCristancho_20240214_2139_H_F <- CellCycleScoring(
  AnaCristancho_20240214_2139_H_F,
  s.features = mmus_s,
  g2m.features = mmus_g2m
)
# Write cell cycle phase summary
phase_summary <- table(AnaCristancho_20240214_2139_H_F$Phase)
phase_summary_df <- as.data.frame(phase_summary)
phase_summary_df$Sample <- "AnaCristancho_20240214_2139_H_F"
write.table(
  phase_summary_df,
  file = "output/seurat/CellCyclePhase_vKeepBottom-AnaCristancho_20240214_2139_H_F.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
############ SAVE sample ########################################
saveRDS(
  AnaCristancho_20240214_2139_H_F,
  file = "output/seurat/AnaCristancho_20240214_2139_H_F-QCPass_vKeepBottom.rds"
)
#################################################################
```

Each sample is saved as an individual seurat RNA clean seurat object.



### QC filtering - ATAC assay

Add fragment file and clean ATAC assay sample per sample.

```R
set.seed(42)

# library
library("Signac")
library("Seurat")
library("tidyverse")
library("EnsDb.Mmusculus.v79") # mm10
library("harmony")
use_python("~/anaconda3/envs/SignacV5/bin/python") 



# Load RNA seurat object
AnaCristancho_20240206_1957_H_M = readRDS("output/seurat/AnaCristancho_20240206_1957_H_M-QCPass_vKeepBottom.rds")
AnaCristancho_20240206_1966_N_M = readRDS("output/seurat/AnaCristancho_20240206_1966_N_M-QCPass_vKeepBottom.rds")
AnaCristancho_20240206_1987_H_F = readRDS("output/seurat/AnaCristancho_20240206_1987_H_F-QCPass_vKeepBottom.rds")
AnaCristancho_20240206_2027_N_F = readRDS("output/seurat/AnaCristancho_20240206_2027_N_F-QCPass_vKeepBottom.rds")
AnaCristancho_20240207_1966_N_F = readRDS("output/seurat/AnaCristancho_20240207_1966_N_F-QCPass_vKeepBottom.rds")

AnaCristancho_20240207_2011_H_M = readRDS("output/seurat/AnaCristancho_20240207_2011_H_M-QCPass_vKeepBottom.rds")
AnaCristancho_20240207_2023_H_F = readRDS("output/seurat/AnaCristancho_20240207_2023_H_F-QCPass_vKeepBottom.rds")
AnaCristancho_20240207_2028_N_M = readRDS("output/seurat/AnaCristancho_20240207_2028_N_M-QCPass_vKeepBottom.rds")
AnaCristancho_20240212_1965_N_F = readRDS("output/seurat/AnaCristancho_20240212_1965_N_F-QCPass_vKeepBottom.rds")
AnaCristancho_20240212_2009_H_M = readRDS("output/seurat/AnaCristancho_20240212_2009_H_M-QCPass_vKeepBottom.rds")

AnaCristancho_20240212_2012_H_F = readRDS("output/seurat/AnaCristancho_20240212_2012_H_F-QCPass_vKeepBottom.rds")
AnaCristancho_20240212_2027_N_M = readRDS("output/seurat/AnaCristancho_20240212_2027_N_M-QCPass_vKeepBottom.rds")
AnaCristancho_20240213_2022_H_M = readRDS("output/seurat/AnaCristancho_20240213_2022_H_M-QCPass_vKeepBottom.rds")
AnaCristancho_20240213_2113_N_M = readRDS("output/seurat/AnaCristancho_20240213_2113_N_M-QCPass_vKeepBottom.rds")
AnaCristancho_20240213_2264_H_F = readRDS("output/seurat/AnaCristancho_20240213_2264_H_F-QCPass_vKeepBottom.rds")

AnaCristancho_20240214_1965_N_M = readRDS("output/seurat/AnaCristancho_20240214_1965_N_M-QCPass_vKeepBottom.rds")
AnaCristancho_20240214_2023_H_M = readRDS("output/seurat/AnaCristancho_20240214_2023_H_M-QCPass_vKeepBottom.rds")
AnaCristancho_20240214_2028_N_F = readRDS("output/seurat/AnaCristancho_20240214_2028_N_F-QCPass_vKeepBottom.rds")
AnaCristancho_20240214_2113_N_F = readRDS("output/seurat/AnaCristancho_20240214_2113_N_F-QCPass_vKeepBottom.rds")
AnaCristancho_20240214_2139_H_F = readRDS("output/seurat/AnaCristancho_20240214_2139_H_F-QCPass_vKeepBottom.rds")


########################
#### Add ATAC ######
########################


##############################################
# AnaCristancho_20240206_1957_H_M ############
##############################################
AnaCristancho_20240206_1957_H_M_h5 <- Read10X_h5("AnaCristancho-20240206-1957-H-M_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 


# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240206_1957_H_M$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240206_1957_H_M_h5$Peaks # Here we load the h5 file

rna_names<-colnames(AnaCristancho_20240206_1957_H_M$RNA)
atac_names<-colnames(AnaCristancho_20240206_1957_H_M_h5$Peaks)

intersect <- intersect(atac_names, rna_names)

intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240206-1957-H-M_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240206_1957_H_M[["ATAC"]] <- chrom_assay



##############################################
# AnaCristancho_20240206_1966_N_M ############
##############################################
AnaCristancho_20240206_1966_N_M_h5 <- Read10X_h5("AnaCristancho-20240206-1966-N-M_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240206_1966_N_M$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240206_1966_N_M_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240206_1966_N_M$RNA)
atac_names<-colnames(AnaCristancho_20240206_1966_N_M_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240206-1966-N-M_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240206_1966_N_M[["ATAC"]] <- chrom_assay




##############################################
# AnaCristancho_20240206_1987_H_F ############
##############################################
AnaCristancho_20240206_1987_H_F_h5 <- Read10X_h5("AnaCristancho-20240206-1987-H-F_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240206_1987_H_F$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240206_1987_H_F_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240206_1987_H_F$RNA)
atac_names<-colnames(AnaCristancho_20240206_1987_H_F_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240206-1987-H-F_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240206_1987_H_F[["ATAC"]] <- chrom_assay



##############################################
# AnaCristancho_20240206_2027_N_F ############
##############################################
AnaCristancho_20240206_2027_N_F_h5 <- Read10X_h5("AnaCristancho-20240206-2027-N-F_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240206_2027_N_F$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240206_2027_N_F_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240206_2027_N_F$RNA)
atac_names<-colnames(AnaCristancho_20240206_2027_N_F_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240206-2027-N-F_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240206_2027_N_F[["ATAC"]] <- chrom_assay



##############################################
# AnaCristancho_20240207_1966_N_F ############
##############################################
AnaCristancho_20240207_1966_N_F_h5 <- Read10X_h5("AnaCristancho-20240207-1966-N-F_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240207_1966_N_F$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240207_1966_N_F_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240207_1966_N_F$RNA)
atac_names<-colnames(AnaCristancho_20240207_1966_N_F_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240207-1966-N-F_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240207_1966_N_F[["ATAC"]] <- chrom_assay




##############################################
# AnaCristancho_20240207_2011_H_M ############
##############################################
AnaCristancho_20240207_2011_H_M_h5 <- Read10X_h5("AnaCristancho-20240207-2011-H-M_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240207_2011_H_M$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240207_2011_H_M_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240207_2011_H_M$RNA)
atac_names<-colnames(AnaCristancho_20240207_2011_H_M_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240207-2011-H-M_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240207_2011_H_M[["ATAC"]] <- chrom_assay




##############################################
# AnaCristancho_20240207_2023_H_F ############
##############################################
AnaCristancho_20240207_2023_H_F_h5 <- Read10X_h5("AnaCristancho-20240207-2023-H-F_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240207_2023_H_F$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240207_2023_H_F_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240207_2023_H_F$RNA)
atac_names<-colnames(AnaCristancho_20240207_2023_H_F_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240207-2023-H-F_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240207_2023_H_F[["ATAC"]] <- chrom_assay




##############################################
# AnaCristancho_20240207_2028_N_M ############
##############################################
AnaCristancho_20240207_2028_N_M_h5 <- Read10X_h5("AnaCristancho-20240207-2028-N-M_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240207_2028_N_M$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240207_2028_N_M_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240207_2028_N_M$RNA)
atac_names<-colnames(AnaCristancho_20240207_2028_N_M_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240207-2028-N-M_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240207_2028_N_M[["ATAC"]] <- chrom_assay




##############################################
# AnaCristancho_20240212_1965_N_F ############
##############################################
AnaCristancho_20240212_1965_N_F_h5 <- Read10X_h5("AnaCristancho-20240212-1965-N-F_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240212_1965_N_F$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240212_1965_N_F_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240212_1965_N_F$RNA)
atac_names<-colnames(AnaCristancho_20240212_1965_N_F_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240212-1965-N-F_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240212_1965_N_F[["ATAC"]] <- chrom_assay




##############################################
# AnaCristancho_20240212_2009_H_M ############
##############################################
AnaCristancho_20240212_2009_H_M_h5 <- Read10X_h5("AnaCristancho-20240212-2009-H-M_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240212_2009_H_M$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240212_2009_H_M_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240212_2009_H_M$RNA)
atac_names<-colnames(AnaCristancho_20240212_2009_H_M_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240212-2009-H-M_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240212_2009_H_M[["ATAC"]] <- chrom_assay





##############################################
# AnaCristancho_20240212_2012_H_F ############
##############################################
AnaCristancho_20240212_2012_H_F_h5 <- Read10X_h5("AnaCristancho-20240212-2012-H-F_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240212_2012_H_F$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240212_2012_H_F_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240212_2012_H_F$RNA)
atac_names<-colnames(AnaCristancho_20240212_2012_H_F_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240212-2012-H-F_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240212_2012_H_F[["ATAC"]] <- chrom_assay




##############################################
# AnaCristancho_20240212_2027_N_M ############
##############################################
AnaCristancho_20240212_2027_N_M_h5 <- Read10X_h5("AnaCristancho-20240212-2027-N-M_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240212_2027_N_M$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240212_2027_N_M_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240212_2027_N_M$RNA)
atac_names<-colnames(AnaCristancho_20240212_2027_N_M_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240212-2027-N-M_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240212_2027_N_M[["ATAC"]] <- chrom_assay



##############################################
# AnaCristancho_20240213_2022_H_M ############
##############################################
AnaCristancho_20240213_2022_H_M_h5 <- Read10X_h5("AnaCristancho-20240213-2022-H-M_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240213_2022_H_M$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240213_2022_H_M_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240213_2022_H_M$RNA)
atac_names<-colnames(AnaCristancho_20240213_2022_H_M_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240213-2022-H-M_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240213_2022_H_M[["ATAC"]] <- chrom_assay



##############################################
# AnaCristancho_20240213_2113_N_M ############
##############################################
AnaCristancho_20240213_2113_N_M_h5 <- Read10X_h5("AnaCristancho-20240213-2113-N-M_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240213_2113_N_M$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240213_2113_N_M_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240213_2113_N_M$RNA)
atac_names<-colnames(AnaCristancho_20240213_2113_N_M_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240213-2113-N-M_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240213_2113_N_M[["ATAC"]] <- chrom_assay




##############################################
# AnaCristancho_20240213_2264_H_F ############
##############################################
AnaCristancho_20240213_2264_H_F_h5 <- Read10X_h5("AnaCristancho-20240213-2264-H-F_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240213_2264_H_F$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240213_2264_H_F_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240213_2264_H_F$RNA)
atac_names<-colnames(AnaCristancho_20240213_2264_H_F_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240213-2264-H-F_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240213_2264_H_F[["ATAC"]] <- chrom_assay



##############################################
# AnaCristancho_20240214_1965_N_M ############
##############################################
AnaCristancho_20240214_1965_N_M_h5 <- Read10X_h5("AnaCristancho-20240214-1965-N-M_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240214_1965_N_M$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240214_1965_N_M_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240214_1965_N_M$RNA)
atac_names<-colnames(AnaCristancho_20240214_1965_N_M_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240214-1965-N-M_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240214_1965_N_M[["ATAC"]] <- chrom_assay



##############################################
# AnaCristancho_20240214_2023_H_M ############
##############################################
AnaCristancho_20240214_2023_H_M_h5 <- Read10X_h5("AnaCristancho-20240214-2023-H-M_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240214_2023_H_M$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240214_2023_H_M_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240214_2023_H_M$RNA)
atac_names<-colnames(AnaCristancho_20240214_2023_H_M_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240214-2023-H-M_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240214_2023_H_M[["ATAC"]] <- chrom_assay



##############################################
# AnaCristancho_20240214_2028_N_F ############
##############################################
AnaCristancho_20240214_2028_N_F_h5 <- Read10X_h5("AnaCristancho-20240214-2028-N-F_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240214_2028_N_F$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240214_2028_N_F_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240214_2028_N_F$RNA)
atac_names<-colnames(AnaCristancho_20240214_2028_N_F_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240214-2028-N-F_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240214_2028_N_F[["ATAC"]] <- chrom_assay



##############################################
# AnaCristancho_20240214_2113_N_F ############
##############################################
AnaCristancho_20240214_2113_N_F_h5 <- Read10X_h5("AnaCristancho-20240214-2113-N-F_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240214_2113_N_F$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240214_2113_N_F_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240214_2113_N_F$RNA)
atac_names<-colnames(AnaCristancho_20240214_2113_N_F_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240214-2113-N-F_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240214_2113_N_F[["ATAC"]] <- chrom_assay



##############################################
# AnaCristancho_20240214_2139_H_F ############
##############################################
AnaCristancho_20240214_2139_H_F_h5 <- Read10X_h5("AnaCristancho-20240214-2139-H-F_Multiome/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 

# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- AnaCristancho_20240214_2139_H_F$RNA # Here we load the soupX scrublet QC clean corrected seurat
atac_counts <- AnaCristancho_20240214_2139_H_F_h5$Peaks # Here we load the h5 file
rna_names<-colnames(AnaCristancho_20240214_2139_H_F$RNA)
atac_names<-colnames(AnaCristancho_20240214_2139_H_F_h5$Peaks)
intersect <- intersect(atac_names, rna_names)
intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "AnaCristancho-20240214-2139-H-F_Multiome/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
AnaCristancho_20240214_2139_H_F[["ATAC"]] <- chrom_assay




##################
## QC ATAC #########
##################

##############################################
# AnaCristancho_20240206_1957_H_M ############
##############################################
DefaultAssay(AnaCristancho_20240206_1957_H_M) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240206_1957_H_M <- NucleosomeSignal(object = AnaCristancho_20240206_1957_H_M)
## compute TSS enrichment score per cell
AnaCristancho_20240206_1957_H_M <- TSSEnrichment(object = AnaCristancho_20240206_1957_H_M, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240206_1957_H_M.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240206_1957_H_M, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240206_1957_H_M.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240206_1957_H_M,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240206_1957_H_M-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240206_1957_H_M,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 10000)
dev.off()
AnaCristancho_20240206_1957_H_M$high.tss <- ifelse(AnaCristancho_20240206_1957_H_M$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240206_1957_H_M.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240206_1957_H_M, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240206_1957_H_M$nucleosome_group <- ifelse(AnaCristancho_20240206_1957_H_M$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240206_1957_H_M.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240206_1957_H_M, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).

## subset cells that pass QC
AnaCristancho_20240206_1957_H_M <- subset(
  x = AnaCristancho_20240206_1957_H_M,
  subset = nCount_ATAC > 2000 &
    nCount_ATAC < 20000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 7.5 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.25
)
AnaCristancho_20240206_1957_H_M
saveRDS(AnaCristancho_20240206_1957_H_M, file = "output/seurat/AnaCristancho_20240206_1957_H_M-QCPass_RNA_ATAC_vKeepBottom.rds") 



##############################################
# AnaCristancho_20240206_1966_N_M ############
##############################################
DefaultAssay(AnaCristancho_20240206_1966_N_M) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240206_1966_N_M <- NucleosomeSignal(object = AnaCristancho_20240206_1966_N_M)
## compute TSS enrichment score per cell
AnaCristancho_20240206_1966_N_M <- TSSEnrichment(object = AnaCristancho_20240206_1966_N_M, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240206_1966_N_M.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240206_1966_N_M, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240206_1966_N_M.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240206_1966_N_M,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240206_1966_N_M-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240206_1966_N_M,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 10000)
dev.off()
AnaCristancho_20240206_1966_N_M$high.tss <- ifelse(AnaCristancho_20240206_1966_N_M$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240206_1966_N_M.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240206_1966_N_M, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240206_1966_N_M$nucleosome_group <- ifelse(AnaCristancho_20240206_1966_N_M$nucleosome_signal > 1.5, 'NS > 1.5', 'NS < 1.5') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240206_1966_N_M.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240206_1966_N_M, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240206_1966_N_M <- subset(
  x = AnaCristancho_20240206_1966_N_M,
  subset = nCount_ATAC > 5000 &
    nCount_ATAC < 50000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 10 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.5
)
AnaCristancho_20240206_1966_N_M

saveRDS(AnaCristancho_20240206_1966_N_M, file = "output/seurat/AnaCristancho_20240206_1966_N_M-QCPass_RNA_ATAC_vKeepBottom.rds") 


##############################################
# AnaCristancho_20240206_1987_H_F ############
##############################################
DefaultAssay(AnaCristancho_20240206_1987_H_F) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240206_1987_H_F <- NucleosomeSignal(object = AnaCristancho_20240206_1987_H_F)
## compute TSS enrichment score per cell
AnaCristancho_20240206_1987_H_F <- TSSEnrichment(object = AnaCristancho_20240206_1987_H_F, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240206_1987_H_F.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240206_1987_H_F, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240206_1987_H_F.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240206_1987_H_F,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240206_1987_H_F-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240206_1987_H_F,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 10000)
dev.off()
AnaCristancho_20240206_1987_H_F$high.tss <- ifelse(AnaCristancho_20240206_1987_H_F$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240206_1987_H_F.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240206_1987_H_F, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240206_1987_H_F$nucleosome_group <- ifelse(AnaCristancho_20240206_1987_H_F$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240206_1987_H_F.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240206_1987_H_F, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240206_1987_H_F <- subset(
  x = AnaCristancho_20240206_1987_H_F,
  subset = nCount_ATAC > 2500 &
    nCount_ATAC < 30000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 12 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.25
)
AnaCristancho_20240206_1987_H_F

saveRDS(AnaCristancho_20240206_1987_H_F, file = "output/seurat/AnaCristancho_20240206_1987_H_F-QCPass_RNA_ATAC_vKeepBottom.rds") 




##############################################
# AnaCristancho_20240206_2027_N_F ############
##############################################
DefaultAssay(AnaCristancho_20240206_2027_N_F) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240206_2027_N_F <- NucleosomeSignal(object = AnaCristancho_20240206_2027_N_F)
## compute TSS enrichment score per cell
AnaCristancho_20240206_2027_N_F <- TSSEnrichment(object = AnaCristancho_20240206_2027_N_F, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240206_2027_N_F.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240206_2027_N_F, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240206_2027_N_F.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240206_2027_N_F,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240206_2027_N_F-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240206_2027_N_F,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 10000)
dev.off()
AnaCristancho_20240206_2027_N_F$high.tss <- ifelse(AnaCristancho_20240206_2027_N_F$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240206_2027_N_F.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240206_2027_N_F, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240206_2027_N_F$nucleosome_group <- ifelse(AnaCristancho_20240206_2027_N_F$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240206_2027_N_F.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240206_2027_N_F, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240206_2027_N_F <- subset(
  x = AnaCristancho_20240206_2027_N_F,
  subset = nCount_ATAC > 2500 &
    nCount_ATAC < 30000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 12 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.25
)
AnaCristancho_20240206_2027_N_F

saveRDS(AnaCristancho_20240206_2027_N_F, file = "output/seurat/AnaCristancho_20240206_2027_N_F-QCPass_RNA_ATAC_vKeepBottom.rds") 





##############################################
# AnaCristancho_20240207_1966_N_F ############
##############################################
DefaultAssay(AnaCristancho_20240207_1966_N_F) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240207_1966_N_F <- NucleosomeSignal(object = AnaCristancho_20240207_1966_N_F)
## compute TSS enrichment score per cell
AnaCristancho_20240207_1966_N_F <- TSSEnrichment(object = AnaCristancho_20240207_1966_N_F, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240207_1966_N_F.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240207_1966_N_F, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240207_1966_N_F.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240207_1966_N_F,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240207_1966_N_F-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240207_1966_N_F,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 10000)
dev.off()
AnaCristancho_20240207_1966_N_F$high.tss <- ifelse(AnaCristancho_20240207_1966_N_F$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240207_1966_N_F.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240207_1966_N_F, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240207_1966_N_F$nucleosome_group <- ifelse(AnaCristancho_20240207_1966_N_F$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240207_1966_N_F.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240207_1966_N_F, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240207_1966_N_F <- subset(
  x = AnaCristancho_20240207_1966_N_F,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 30000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 12 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.25
)
AnaCristancho_20240207_1966_N_F

saveRDS(AnaCristancho_20240207_1966_N_F, file = "output/seurat/AnaCristancho_20240207_1966_N_F-QCPass_RNA_ATAC_vKeepBottom.rds") 






##############################################
# AnaCristancho_20240207_2011_H_M ############
##############################################
DefaultAssay(AnaCristancho_20240207_2011_H_M) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240207_2011_H_M <- NucleosomeSignal(object = AnaCristancho_20240207_2011_H_M)
## compute TSS enrichment score per cell
AnaCristancho_20240207_2011_H_M <- TSSEnrichment(object = AnaCristancho_20240207_2011_H_M, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240207_2011_H_M.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240207_2011_H_M, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240207_2011_H_M.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240207_2011_H_M,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240207_2011_H_M-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240207_2011_H_M,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 10000)
dev.off()
AnaCristancho_20240207_2011_H_M$high.tss <- ifelse(AnaCristancho_20240207_2011_H_M$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240207_2011_H_M.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240207_2011_H_M, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240207_2011_H_M$nucleosome_group <- ifelse(AnaCristancho_20240207_2011_H_M$nucleosome_signal > 1.50, 'NS > 1.50', 'NS < 1.50') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240207_2011_H_M.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240207_2011_H_M, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240207_2011_H_M <- subset(
  x = AnaCristancho_20240207_2011_H_M,
  subset = nCount_ATAC > 3000 &
    nCount_ATAC < 50000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 10 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.5
)
AnaCristancho_20240207_2011_H_M

saveRDS(AnaCristancho_20240207_2011_H_M, file = "output/seurat/AnaCristancho_20240207_2011_H_M-QCPass_RNA_ATAC_vKeepBottom.rds") 




##############################################
# AnaCristancho_20240207_2023_H_F ############
##############################################
DefaultAssay(AnaCristancho_20240207_2023_H_F) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240207_2023_H_F <- NucleosomeSignal(object = AnaCristancho_20240207_2023_H_F)
## compute TSS enrichment score per cell
AnaCristancho_20240207_2023_H_F <- TSSEnrichment(object = AnaCristancho_20240207_2023_H_F, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240207_2023_H_F.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240207_2023_H_F, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240207_2023_H_F.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240207_2023_H_F,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240207_2023_H_F-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240207_2023_H_F,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 2500)
dev.off()
AnaCristancho_20240207_2023_H_F$high.tss <- ifelse(AnaCristancho_20240207_2023_H_F$TSS.enrichment > 1.5, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240207_2023_H_F.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240207_2023_H_F, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240207_2023_H_F$nucleosome_group <- ifelse(AnaCristancho_20240207_2023_H_F$nucleosome_signal > 1.50, 'NS > 1.50', 'NS < 1.50') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240207_2023_H_F.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240207_2023_H_F, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240207_2023_H_F <- subset(
  x = AnaCristancho_20240207_2023_H_F,
  subset = nCount_ATAC > 250 &
    nCount_ATAC < 20000  &
    TSS.enrichment > 1.5 &
    TSS.enrichment < 12 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.5
)
AnaCristancho_20240207_2023_H_F

saveRDS(AnaCristancho_20240207_2023_H_F, file = "output/seurat/AnaCristancho_20240207_2023_H_F-QCPass_RNA_ATAC_vKeepBottom.rds") 






##############################################
# AnaCristancho_20240207_2028_N_M ############
##############################################
DefaultAssay(AnaCristancho_20240207_2028_N_M) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240207_2028_N_M <- NucleosomeSignal(object = AnaCristancho_20240207_2028_N_M)
## compute TSS enrichment score per cell
AnaCristancho_20240207_2028_N_M <- TSSEnrichment(object = AnaCristancho_20240207_2028_N_M, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240207_2028_N_M.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240207_2028_N_M, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240207_2028_N_M.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240207_2028_N_M,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240207_2028_N_M-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240207_2028_N_M,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 5000)
dev.off()
AnaCristancho_20240207_2028_N_M$high.tss <- ifelse(AnaCristancho_20240207_2028_N_M$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240207_2028_N_M.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240207_2028_N_M, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240207_2028_N_M$nucleosome_group <- ifelse(AnaCristancho_20240207_2028_N_M$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240207_2028_N_M.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240207_2028_N_M, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240207_2028_N_M <- subset(
  x = AnaCristancho_20240207_2028_N_M,
  subset = nCount_ATAC > 2000 &
    nCount_ATAC < 30000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 12 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.5
)
AnaCristancho_20240207_2028_N_M

saveRDS(AnaCristancho_20240207_2028_N_M, file = "output/seurat/AnaCristancho_20240207_2028_N_M-QCPass_RNA_ATAC_vKeepBottom.rds") 





##############################################
# AnaCristancho_20240212_1965_N_F ############
##############################################
DefaultAssay(AnaCristancho_20240212_1965_N_F) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240212_1965_N_F <- NucleosomeSignal(object = AnaCristancho_20240212_1965_N_F)
## compute TSS enrichment score per cell
AnaCristancho_20240212_1965_N_F <- TSSEnrichment(object = AnaCristancho_20240212_1965_N_F, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240212_1965_N_F.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240212_1965_N_F, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240212_1965_N_F.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240212_1965_N_F,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240212_1965_N_F-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240212_1965_N_F,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 5000)
dev.off()
AnaCristancho_20240212_1965_N_F$high.tss <- ifelse(AnaCristancho_20240212_1965_N_F$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240212_1965_N_F.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240212_1965_N_F, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240212_1965_N_F$nucleosome_group <- ifelse(AnaCristancho_20240212_1965_N_F$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240212_1965_N_F.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240212_1965_N_F, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240212_1965_N_F <- subset(
  x = AnaCristancho_20240212_1965_N_F,
  subset = nCount_ATAC > 2000 &
    nCount_ATAC < 25000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 12 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.25
)
AnaCristancho_20240212_1965_N_F

saveRDS(AnaCristancho_20240212_1965_N_F, file = "output/seurat/AnaCristancho_20240212_1965_N_F-QCPass_RNA_ATAC_vKeepBottom.rds") 







##############################################
# AnaCristancho_20240212_2009_H_M ############
##############################################
DefaultAssay(AnaCristancho_20240212_2009_H_M) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240212_2009_H_M <- NucleosomeSignal(object = AnaCristancho_20240212_2009_H_M)
## compute TSS enrichment score per cell
AnaCristancho_20240212_2009_H_M <- TSSEnrichment(object = AnaCristancho_20240212_2009_H_M, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240212_2009_H_M.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240212_2009_H_M, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240212_2009_H_M.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240212_2009_H_M,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240212_2009_H_M-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240212_2009_H_M,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 5000)
dev.off()
AnaCristancho_20240212_2009_H_M$high.tss <- ifelse(AnaCristancho_20240212_2009_H_M$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240212_2009_H_M.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240212_2009_H_M, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240212_2009_H_M$nucleosome_group <- ifelse(AnaCristancho_20240212_2009_H_M$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240212_2009_H_M.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240212_2009_H_M, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240212_2009_H_M <- subset(
  x = AnaCristancho_20240212_2009_H_M,
  subset = nCount_ATAC > 2000 &
    nCount_ATAC < 20000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 12 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.25
)
AnaCristancho_20240212_2009_H_M

saveRDS(AnaCristancho_20240212_2009_H_M, file = "output/seurat/AnaCristancho_20240212_2009_H_M-QCPass_RNA_ATAC_vKeepBottom.rds") 






##############################################
# AnaCristancho_20240212_2012_H_F ############
##############################################
DefaultAssay(AnaCristancho_20240212_2012_H_F) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240212_2012_H_F <- NucleosomeSignal(object = AnaCristancho_20240212_2012_H_F)
## compute TSS enrichment score per cell
AnaCristancho_20240212_2012_H_F <- TSSEnrichment(object = AnaCristancho_20240212_2012_H_F, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240212_2012_H_F.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240212_2012_H_F, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240212_2012_H_F.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240212_2012_H_F,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240212_2012_H_F-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240212_2012_H_F,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 5000)
dev.off()
AnaCristancho_20240212_2012_H_F$high.tss <- ifelse(AnaCristancho_20240212_2012_H_F$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240212_2012_H_F.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240212_2012_H_F, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240212_2012_H_F$nucleosome_group <- ifelse(AnaCristancho_20240212_2012_H_F$nucleosome_signal > 1.35, 'NS > 1.35', 'NS < 1.35') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240212_2012_H_F.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240212_2012_H_F, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240212_2012_H_F <- subset(
  x = AnaCristancho_20240212_2012_H_F,
  subset = nCount_ATAC > 3000 &
    nCount_ATAC < 30000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 12 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.35
)
AnaCristancho_20240212_2012_H_F

saveRDS(AnaCristancho_20240212_2012_H_F, file = "output/seurat/AnaCristancho_20240212_2012_H_F-QCPass_RNA_ATAC_vKeepBottom.rds") 







##############################################
# AnaCristancho_20240212_2027_N_M ############
##############################################
DefaultAssay(AnaCristancho_20240212_2027_N_M) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240212_2027_N_M <- NucleosomeSignal(object = AnaCristancho_20240212_2027_N_M)
## compute TSS enrichment score per cell
AnaCristancho_20240212_2027_N_M <- TSSEnrichment(object = AnaCristancho_20240212_2027_N_M, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240212_2027_N_M.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240212_2027_N_M, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240212_2027_N_M.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240212_2027_N_M,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240212_2027_N_M-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240212_2027_N_M,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 5000)
dev.off()
AnaCristancho_20240212_2027_N_M$high.tss <- ifelse(AnaCristancho_20240212_2027_N_M$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240212_2027_N_M.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240212_2027_N_M, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240212_2027_N_M$nucleosome_group <- ifelse(AnaCristancho_20240212_2027_N_M$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240212_2027_N_M.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240212_2027_N_M, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240212_2027_N_M <- subset(
  x = AnaCristancho_20240212_2027_N_M,
  subset = nCount_ATAC > 2500 &
    nCount_ATAC < 25000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 12 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.25
)
AnaCristancho_20240212_2027_N_M

saveRDS(AnaCristancho_20240212_2027_N_M, file = "output/seurat/AnaCristancho_20240212_2027_N_M-QCPass_RNA_ATAC_vKeepBottom.rds") 





##############################################
# AnaCristancho_20240213_2022_H_M ############
##############################################
DefaultAssay(AnaCristancho_20240213_2022_H_M) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240213_2022_H_M <- NucleosomeSignal(object = AnaCristancho_20240213_2022_H_M)
## compute TSS enrichment score per cell
AnaCristancho_20240213_2022_H_M <- TSSEnrichment(object = AnaCristancho_20240213_2022_H_M, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240213_2022_H_M.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240213_2022_H_M, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240213_2022_H_M.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240213_2022_H_M,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240213_2022_H_M-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240213_2022_H_M,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 5000)
dev.off()
AnaCristancho_20240213_2022_H_M$high.tss <- ifelse(AnaCristancho_20240213_2022_H_M$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240213_2022_H_M.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240213_2022_H_M, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240213_2022_H_M$nucleosome_group <- ifelse(AnaCristancho_20240213_2022_H_M$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240213_2022_H_M.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240213_2022_H_M, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240213_2022_H_M <- subset(
  x = AnaCristancho_20240213_2022_H_M,
  subset = nCount_ATAC > 3000 &
    nCount_ATAC < 30000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 15 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.25
)
AnaCristancho_20240213_2022_H_M

saveRDS(AnaCristancho_20240213_2022_H_M, file = "output/seurat/AnaCristancho_20240213_2022_H_M-QCPass_RNA_ATAC_vKeepBottom.rds") 




##############################################
# AnaCristancho_20240213_2113_N_M ############
##############################################
DefaultAssay(AnaCristancho_20240213_2113_N_M) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240213_2113_N_M <- NucleosomeSignal(object = AnaCristancho_20240213_2113_N_M)
## compute TSS enrichment score per cell
AnaCristancho_20240213_2113_N_M <- TSSEnrichment(object = AnaCristancho_20240213_2113_N_M, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240213_2113_N_M.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240213_2113_N_M, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240213_2113_N_M.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240213_2113_N_M,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240213_2113_N_M-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240213_2113_N_M,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 5000)
dev.off()
AnaCristancho_20240213_2113_N_M$high.tss <- ifelse(AnaCristancho_20240213_2113_N_M$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240213_2113_N_M.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240213_2113_N_M, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240213_2113_N_M$nucleosome_group <- ifelse(AnaCristancho_20240213_2113_N_M$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240213_2113_N_M.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240213_2113_N_M, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240213_2113_N_M <- subset(
  x = AnaCristancho_20240213_2113_N_M,
  subset = nCount_ATAC > 3000 &
    nCount_ATAC < 50000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 15 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.25
)
AnaCristancho_20240213_2113_N_M

saveRDS(AnaCristancho_20240213_2113_N_M, file = "output/seurat/AnaCristancho_20240213_2113_N_M-QCPass_RNA_ATAC_vKeepBottom.rds") 






##############################################
# AnaCristancho_20240213_2264_H_F ############
##############################################
DefaultAssay(AnaCristancho_20240213_2264_H_F) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240213_2264_H_F <- NucleosomeSignal(object = AnaCristancho_20240213_2264_H_F)
## compute TSS enrichment score per cell
AnaCristancho_20240213_2264_H_F <- TSSEnrichment(object = AnaCristancho_20240213_2264_H_F, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240213_2264_H_F.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240213_2264_H_F, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240213_2264_H_F.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240213_2264_H_F,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240213_2264_H_F-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240213_2264_H_F,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 5000)
dev.off()
AnaCristancho_20240213_2264_H_F$high.tss <- ifelse(AnaCristancho_20240213_2264_H_F$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240213_2264_H_F.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240213_2264_H_F, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240213_2264_H_F$nucleosome_group <- ifelse(AnaCristancho_20240213_2264_H_F$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240213_2264_H_F.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240213_2264_H_F, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240213_2264_H_F <- subset(
  x = AnaCristancho_20240213_2264_H_F,
  subset = nCount_ATAC > 2000 &
    nCount_ATAC < 20000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 15 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.25
)
AnaCristancho_20240213_2264_H_F

saveRDS(AnaCristancho_20240213_2264_H_F, file = "output/seurat/AnaCristancho_20240213_2264_H_F-QCPass_RNA_ATAC_vKeepBottom.rds") 



##############################################
# AnaCristancho_20240214_1965_N_M ############
##############################################
DefaultAssay(AnaCristancho_20240214_1965_N_M) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240214_1965_N_M <- NucleosomeSignal(object = AnaCristancho_20240214_1965_N_M)
## compute TSS enrichment score per cell
AnaCristancho_20240214_1965_N_M <- TSSEnrichment(object = AnaCristancho_20240214_1965_N_M, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240214_1965_N_M.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240214_1965_N_M, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240214_1965_N_M.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240214_1965_N_M,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240214_1965_N_M-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240214_1965_N_M,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 5000)
dev.off()
AnaCristancho_20240214_1965_N_M$high.tss <- ifelse(AnaCristancho_20240214_1965_N_M$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240214_1965_N_M.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240214_1965_N_M, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240214_1965_N_M$nucleosome_group <- ifelse(AnaCristancho_20240214_1965_N_M$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240214_1965_N_M.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240214_1965_N_M, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240214_1965_N_M <- subset(
  x = AnaCristancho_20240214_1965_N_M,
  subset = nCount_ATAC > 3000 &
    nCount_ATAC < 30000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 12 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.25
)
AnaCristancho_20240214_1965_N_M

saveRDS(AnaCristancho_20240214_1965_N_M, file = "output/seurat/AnaCristancho_20240214_1965_N_M-QCPass_RNA_ATAC_vKeepBottom.rds") 




##############################################
# AnaCristancho_20240214_2023_H_M ############
##############################################
DefaultAssay(AnaCristancho_20240214_2023_H_M) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240214_2023_H_M <- NucleosomeSignal(object = AnaCristancho_20240214_2023_H_M)
## compute TSS enrichment score per cell
AnaCristancho_20240214_2023_H_M <- TSSEnrichment(object = AnaCristancho_20240214_2023_H_M, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240214_2023_H_M.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240214_2023_H_M, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240214_2023_H_M.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240214_2023_H_M,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240214_2023_H_M-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240214_2023_H_M,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 5000)
dev.off()
AnaCristancho_20240214_2023_H_M$high.tss <- ifelse(AnaCristancho_20240214_2023_H_M$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240214_2023_H_M.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240214_2023_H_M, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240214_2023_H_M$nucleosome_group <- ifelse(AnaCristancho_20240214_2023_H_M$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240214_2023_H_M.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240214_2023_H_M, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240214_2023_H_M <- subset(
  x = AnaCristancho_20240214_2023_H_M,
  subset = nCount_ATAC > 3000 &
    nCount_ATAC < 30000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 15 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.25
)
AnaCristancho_20240214_2023_H_M

saveRDS(AnaCristancho_20240214_2023_H_M, file = "output/seurat/AnaCristancho_20240214_2023_H_M-QCPass_RNA_ATAC_vKeepBottom.rds") 




##############################################
# AnaCristancho_20240214_2028_N_F ############
##############################################
DefaultAssay(AnaCristancho_20240214_2028_N_F) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240214_2028_N_F <- NucleosomeSignal(object = AnaCristancho_20240214_2028_N_F)
## compute TSS enrichment score per cell
AnaCristancho_20240214_2028_N_F <- TSSEnrichment(object = AnaCristancho_20240214_2028_N_F, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240214_2028_N_F.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240214_2028_N_F, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240214_2028_N_F.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240214_2028_N_F,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240214_2028_N_F-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240214_2028_N_F,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 5000)
dev.off()
AnaCristancho_20240214_2028_N_F$high.tss <- ifelse(AnaCristancho_20240214_2028_N_F$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240214_2028_N_F.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240214_2028_N_F, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240214_2028_N_F$nucleosome_group <- ifelse(AnaCristancho_20240214_2028_N_F$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240214_2028_N_F.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240214_2028_N_F, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240214_2028_N_F <- subset(
  x = AnaCristancho_20240214_2028_N_F,
  subset = nCount_ATAC > 3000 &
    nCount_ATAC < 60000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 10 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.25
)
AnaCristancho_20240214_2028_N_F

saveRDS(AnaCristancho_20240214_2028_N_F, file = "output/seurat/AnaCristancho_20240214_2028_N_F-QCPass_RNA_ATAC_vKeepBottom.rds") 





##############################################
# AnaCristancho_20240214_2113_N_F ############
##############################################
DefaultAssay(AnaCristancho_20240214_2113_N_F) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240214_2113_N_F <- NucleosomeSignal(object = AnaCristancho_20240214_2113_N_F)
## compute TSS enrichment score per cell
AnaCristancho_20240214_2113_N_F <- TSSEnrichment(object = AnaCristancho_20240214_2113_N_F, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240214_2113_N_F.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240214_2113_N_F, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240214_2113_N_F.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240214_2113_N_F,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240214_2113_N_F-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240214_2113_N_F,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 5000)
dev.off()
AnaCristancho_20240214_2113_N_F$high.tss <- ifelse(AnaCristancho_20240214_2113_N_F$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240214_2113_N_F.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240214_2113_N_F, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240214_2113_N_F$nucleosome_group <- ifelse(AnaCristancho_20240214_2113_N_F$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240214_2113_N_F.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240214_2113_N_F, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240214_2113_N_F <- subset(
  x = AnaCristancho_20240214_2113_N_F,
  subset = nCount_ATAC > 2000 &
    nCount_ATAC < 40000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 15 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.25
)
AnaCristancho_20240214_2113_N_F

saveRDS(AnaCristancho_20240214_2113_N_F, file = "output/seurat/AnaCristancho_20240214_2113_N_F-QCPass_RNA_ATAC_vKeepBottom.rds") 





##############################################
# AnaCristancho_20240214_2139_H_F ############
##############################################
DefaultAssay(AnaCristancho_20240214_2139_H_F) <- "ATAC" # 

# QC metrics
## compute nucleosome signal score per cell
AnaCristancho_20240214_2139_H_F <- NucleosomeSignal(object = AnaCristancho_20240214_2139_H_F)
## compute TSS enrichment score per cell
AnaCristancho_20240214_2139_H_F <- TSSEnrichment(object = AnaCristancho_20240214_2139_H_F, fast = FALSE)

### QC plots
pdf("output/seurat/DensityScatter_vKeepBottom-AnaCristancho_20240214_2139_H_F.pdf", width=5, height=5)
DensityScatter(AnaCristancho_20240214_2139_H_F, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240214_2139_H_F.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240214_2139_H_F,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot_QC_vKeepBottom-AnaCristancho_20240214_2139_H_F-nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = AnaCristancho_20240214_2139_H_F,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 5000)
dev.off()
AnaCristancho_20240214_2139_H_F$high.tss <- ifelse(AnaCristancho_20240214_2139_H_F$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot_QC_vKeepBottom-AnaCristancho_20240214_2139_H_F.pdf", width=12, height=6)
TSSPlot(AnaCristancho_20240214_2139_H_F, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

AnaCristancho_20240214_2139_H_F$nucleosome_group <- ifelse(AnaCristancho_20240214_2139_H_F$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/seurat/FragmentHistogram_QC_vKeepBottom-AnaCristancho_20240214_2139_H_F.pdf", width=12, height=6)
FragmentHistogram(object = AnaCristancho_20240214_2139_H_F, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).


## subset cells that pass QC
AnaCristancho_20240214_2139_H_F <- subset(
  x = AnaCristancho_20240214_2139_H_F,
  subset = nCount_ATAC > 2500 &
    nCount_ATAC < 25000  &
    TSS.enrichment > 2 &
    TSS.enrichment < 15 &
    nucleosome_signal > 0.1 &
    nucleosome_signal < 1.25
)
AnaCristancho_20240214_2139_H_F

saveRDS(AnaCristancho_20240214_2139_H_F, file = "output/seurat/AnaCristancho_20240214_2139_H_F-QCPass_RNA_ATAC_vKeepBottom.rds") 
```





### Sample integration with SCT


```R
# LOAD all clean seurat object
AnaCristancho_20240206_1957_H_M = readRDS("output/seurat/AnaCristancho_20240206_1957_H_M-QCPass_RNA_ATAC_vKeepBottom.rds")
AnaCristancho_20240206_1966_N_M = readRDS("output/seurat/AnaCristancho_20240206_1966_N_M-QCPass_RNA_ATAC_vKeepBottom.rds")
AnaCristancho_20240206_1987_H_F = readRDS("output/seurat/AnaCristancho_20240206_1987_H_F-QCPass_RNA_ATAC_vKeepBottom.rds")
AnaCristancho_20240206_2027_N_F = readRDS("output/seurat/AnaCristancho_20240206_2027_N_F-QCPass_RNA_ATAC_vKeepBottom.rds")
AnaCristancho_20240207_1966_N_F = readRDS("output/seurat/AnaCristancho_20240207_1966_N_F-QCPass_RNA_ATAC_vKeepBottom.rds")

AnaCristancho_20240207_2011_H_M = readRDS("output/seurat/AnaCristancho_20240207_2011_H_M-QCPass_RNA_ATAC_vKeepBottom.rds")
AnaCristancho_20240207_2023_H_F = readRDS("output/seurat/AnaCristancho_20240207_2023_H_F-QCPass_RNA_ATAC_vKeepBottom.rds")
AnaCristancho_20240207_2028_N_M = readRDS("output/seurat/AnaCristancho_20240207_2028_N_M-QCPass_RNA_ATAC_vKeepBottom.rds")
AnaCristancho_20240212_1965_N_F = readRDS("output/seurat/AnaCristancho_20240212_1965_N_F-QCPass_RNA_ATAC_vKeepBottom.rds")
AnaCristancho_20240212_2009_H_M = readRDS("output/seurat/AnaCristancho_20240212_2009_H_M-QCPass_RNA_ATAC_vKeepBottom.rds")

AnaCristancho_20240212_2012_H_F = readRDS("output/seurat/AnaCristancho_20240212_2012_H_F-QCPass_RNA_ATAC_vKeepBottom.rds")
AnaCristancho_20240212_2027_N_M = readRDS("output/seurat/AnaCristancho_20240212_2027_N_M-QCPass_RNA_ATAC_vKeepBottom.rds")
AnaCristancho_20240213_2022_H_M = readRDS("output/seurat/AnaCristancho_20240213_2022_H_M-QCPass_RNA_ATAC_vKeepBottom.rds")
AnaCristancho_20240213_2113_N_M = readRDS("output/seurat/AnaCristancho_20240213_2113_N_M-QCPass_RNA_ATAC_vKeepBottom.rds")
AnaCristancho_20240213_2264_H_F = readRDS("output/seurat/AnaCristancho_20240213_2264_H_F-QCPass_RNA_ATAC_vKeepBottom.rds")

AnaCristancho_20240214_1965_N_M = readRDS("output/seurat/AnaCristancho_20240214_1965_N_M-QCPass_RNA_ATAC_vKeepBottom.rds")
AnaCristancho_20240214_2023_H_M = readRDS("output/seurat/AnaCristancho_20240214_2023_H_M-QCPass_RNA_ATAC_vKeepBottom.rds")
AnaCristancho_20240214_2028_N_F = readRDS("output/seurat/AnaCristancho_20240214_2028_N_F-QCPass_RNA_ATAC_vKeepBottom.rds")
AnaCristancho_20240214_2113_N_F = readRDS("output/seurat/AnaCristancho_20240214_2113_N_F-QCPass_RNA_ATAC_vKeepBottom.rds")
AnaCristancho_20240214_2139_H_F = readRDS("output/seurat/AnaCristancho_20240214_2139_H_F-QCPass_RNA_ATAC_vKeepBottom.rds")

# Integrate All samples (Run in Slurm)
DefaultAssay(AnaCristancho_20240206_1957_H_M) <- "RNA"
DefaultAssay(AnaCristancho_20240206_1966_N_M) <- "RNA"
DefaultAssay(AnaCristancho_20240206_1987_H_F) <- "RNA"
DefaultAssay(AnaCristancho_20240206_2027_N_F) <- "RNA"
DefaultAssay(AnaCristancho_20240207_2011_H_M) <- "RNA"

DefaultAssay(AnaCristancho_20240207_1966_N_F) <- "RNA"
DefaultAssay(AnaCristancho_20240207_2023_H_F) <- "RNA"
DefaultAssay(AnaCristancho_20240207_2028_N_M) <- "RNA"
DefaultAssay(AnaCristancho_20240212_1965_N_F) <- "RNA"
DefaultAssay(AnaCristancho_20240212_2009_H_M) <- "RNA"

DefaultAssay(AnaCristancho_20240212_2012_H_F) <- "RNA"
DefaultAssay(AnaCristancho_20240212_2027_N_M) <- "RNA"
DefaultAssay(AnaCristancho_20240213_2022_H_M) <- "RNA"
DefaultAssay(AnaCristancho_20240213_2113_N_M) <- "RNA"
DefaultAssay(AnaCristancho_20240213_2264_H_F) <- "RNA"

DefaultAssay(AnaCristancho_20240214_1965_N_M) <- "RNA"
DefaultAssay(AnaCristancho_20240214_2023_H_M) <- "RNA"
DefaultAssay(AnaCristancho_20240214_2028_N_F) <- "RNA"
DefaultAssay(AnaCristancho_20240214_2113_N_F) <- "RNA"
DefaultAssay(AnaCristancho_20240214_2139_H_F) <- "RNA"

options(future.globals.maxSize = 8000 * 1024^2) # options(future.globals.maxSize = 8000 * 1024^2) [TO RUN to avoid error: Error in getGlobalsAndPackages(expr, envir = envir, globals = globals) :]

AnaCristancho_20240206_1957_H_M <- SCTransform(AnaCristancho_20240206_1957_H_M, method = "glmGamPoi", ncells = 5568, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)
AnaCristancho_20240206_1966_N_M <- SCTransform(AnaCristancho_20240206_1966_N_M, method = "glmGamPoi", ncells = 6188, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)
AnaCristancho_20240206_1987_H_F <- SCTransform(AnaCristancho_20240206_1987_H_F, method = "glmGamPoi", ncells = 8774, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)
AnaCristancho_20240206_2027_N_F <- SCTransform(AnaCristancho_20240206_2027_N_F, method = "glmGamPoi", ncells = 9007, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)
AnaCristancho_20240207_2011_H_M <- SCTransform(AnaCristancho_20240207_2011_H_M, method = "glmGamPoi", ncells = 6968, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)

AnaCristancho_20240207_1966_N_F <- SCTransform(AnaCristancho_20240207_1966_N_F, method = "glmGamPoi", ncells = 14386, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000) 
AnaCristancho_20240207_2023_H_F <- SCTransform(AnaCristancho_20240207_2023_H_F, method = "glmGamPoi", ncells = 14078, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)
AnaCristancho_20240207_2028_N_M <- SCTransform(AnaCristancho_20240207_2028_N_M, method = "glmGamPoi", ncells = 9462, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)
AnaCristancho_20240212_1965_N_F <- SCTransform(AnaCristancho_20240212_1965_N_F, method = "glmGamPoi", ncells = 7240, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)
AnaCristancho_20240212_2009_H_M <- SCTransform(AnaCristancho_20240212_2009_H_M, method = "glmGamPoi", ncells = 6679, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)

AnaCristancho_20240212_2012_H_F <- SCTransform(AnaCristancho_20240212_2012_H_F, method = "glmGamPoi", ncells = 5562, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)
AnaCristancho_20240212_2027_N_M <- SCTransform(AnaCristancho_20240212_2027_N_M, method = "glmGamPoi", ncells = 8631, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)
AnaCristancho_20240213_2022_H_M <- SCTransform(AnaCristancho_20240213_2022_H_M, method = "glmGamPoi", ncells = 5319, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)
AnaCristancho_20240213_2113_N_M <- SCTransform(AnaCristancho_20240213_2113_N_M, method = "glmGamPoi", ncells = 7998, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)
AnaCristancho_20240213_2264_H_F <- SCTransform(AnaCristancho_20240213_2264_H_F, method = "glmGamPoi", ncells = 6763, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)

AnaCristancho_20240214_1965_N_M <- SCTransform(AnaCristancho_20240214_1965_N_M, method = "glmGamPoi", ncells = 3531, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)
AnaCristancho_20240214_2023_H_M <- SCTransform(AnaCristancho_20240214_2023_H_M, method = "glmGamPoi", ncells = 8013, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)
AnaCristancho_20240214_2028_N_F <- SCTransform(AnaCristancho_20240214_2028_N_F, method = "glmGamPoi", ncells = 4246, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)
AnaCristancho_20240214_2113_N_F <- SCTransform(AnaCristancho_20240214_2113_N_F, method = "glmGamPoi", ncells = 8648, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)
AnaCristancho_20240214_2139_H_F <- SCTransform(AnaCristancho_20240214_2139_H_F, method = "glmGamPoi", ncells = 7614, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), verbose = TRUE, variable.features.n = 2000)


# Data integration (check active assay is 'SCT')
srat.list <- list(AnaCristancho_20240206_1957_H_M = AnaCristancho_20240206_1957_H_M, AnaCristancho_20240206_1966_N_M = AnaCristancho_20240206_1966_N_M, AnaCristancho_20240206_1987_H_F = AnaCristancho_20240206_1987_H_F, AnaCristancho_20240206_2027_N_F = AnaCristancho_20240206_2027_N_F, AnaCristancho_20240207_2011_H_M = AnaCristancho_20240207_2011_H_M, AnaCristancho_20240207_1966_N_F = AnaCristancho_20240207_1966_N_F, AnaCristancho_20240207_2023_H_F = AnaCristancho_20240207_2023_H_F, AnaCristancho_20240207_2028_N_M = AnaCristancho_20240207_2028_N_M, AnaCristancho_20240212_1965_N_F = AnaCristancho_20240212_1965_N_F, AnaCristancho_20240212_2009_H_M = AnaCristancho_20240212_2009_H_M, AnaCristancho_20240212_2012_H_F = AnaCristancho_20240212_2012_H_F, AnaCristancho_20240212_2027_N_M = AnaCristancho_20240212_2027_N_M, AnaCristancho_20240213_2022_H_M = AnaCristancho_20240213_2022_H_M, AnaCristancho_20240213_2113_N_M = AnaCristancho_20240213_2113_N_M, AnaCristancho_20240213_2264_H_F = AnaCristancho_20240213_2264_H_F, AnaCristancho_20240214_1965_N_M = AnaCristancho_20240214_1965_N_M,AnaCristancho_20240214_2023_H_M = AnaCristancho_20240214_2023_H_M,AnaCristancho_20240214_2028_N_F = AnaCristancho_20240214_2028_N_F,AnaCristancho_20240214_2113_N_F = AnaCristancho_20240214_2113_N_F,AnaCristancho_20240214_2139_H_F = AnaCristancho_20240214_2139_H_F)

features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 2000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)

srat.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
SmithM_multiome.sct <- IntegrateData(anchorset = srat.anchors, normalization.method = "SCT")




### SET sex and condition in my seurat object ###
#### get orig.ident as character
orig <- as.character(SmithM_multiome.sct$orig.ident)
#### extract final pattern like "N_F" or "H_M" (allows "_" or "-")
m <- regmatches(orig, regexec("([NH])[_-]([FM])$", orig))
#### group captures (NA if no match)
cond_init <- sapply(m, function(x) if(length(x) >= 3) toupper(x[2]) else NA_character_)
sex_init  <- sapply(m, function(x) if(length(x) >= 3) toupper(x[3]) else NA_character_)
#### map to full labels
SmithM_multiome.sct$condition <- ifelse(is.na(cond_init), NA,
  ifelse(cond_init == "N", "Normoxia",
    ifelse(cond_init == "H", "Hypoxia", cond_init)))
SmithM_multiome.sct$sex <- ifelse(is.na(sex_init), NA,
  ifelse(sex_init == "F", "Female",
    ifelse(sex_init == "M", "Male", sex_init)))
#### make factors with sensible order
SmithM_multiome.sct$condition <- factor(SmithM_multiome.sct$condition, levels = c("Normoxia", "Hypoxia"))
SmithM_multiome.sct$sex <- factor(SmithM_multiome.sct$sex, levels = c("Female", "Male"))
####################################################




DefaultAssay(SmithM_multiome.sct) <- "integrated"

SmithM_multiome.sct <- RunPCA(SmithM_multiome.sct, verbose = FALSE, npcs = 50)
SmithM_multiome.sct <- RunUMAP(SmithM_multiome.sct, reduction = "pca", dims = 1:50, verbose = FALSE)
SmithM_multiome.sct <- FindNeighbors(SmithM_multiome.sct, reduction = "pca", k.param = 50, dims = 1:50)
SmithM_multiome.sct <- FindClusters(SmithM_multiome.sct, resolution = 0.27, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


pdf("output/seurat/SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-seurat_clusters.pdf", width = 7, height = 6)
DimPlot(SmithM_multiome.sct, group.by = "seurat_clusters", reduction = "umap") + ggtitle("mtrbNotRegress_vKeepBottom_dim50kparam50res025algo4")
dev.off()

pdf("output/seurat/SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-seurat_clustersLABEL.pdf", width = 7, height = 6)
DimPlot(SmithM_multiome.sct, group.by = "seurat_clusters", reduction = "umap",  label = TRUE, repel = TRUE) + ggtitle("mtrbNotRegress_vKeepBottom_dim50kparam50res025algo4")
dev.off()



pdf("output/seurat/DimPlot-SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-condition.pdf", width = 7, height = 6)
DimPlot(SmithM_multiome.sct, group.by = "condition", reduction = "umap",
        cols = c("blue","red")) + ggtitle("condition")
dev.off()
pdf("output/seurat/DimPlot-SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-SPLITcondition.pdf", width = 12, height = 6)
DimPlot(SmithM_multiome.sct, split.by = "condition", reduction = "umap") + ggtitle("condition")
dev.off()
pdf("output/seurat/DimPlot-SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-sex.pdf", width = 7, height = 6)
DimPlot(SmithM_multiome.sct, group.by = "sex", reduction = "umap") + ggtitle("sex")
dev.off()
pdf("output/seurat/DimPlot-SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-SPLITsex.pdf", width = 12, height = 6)
DimPlot(SmithM_multiome.sct, split.by = "sex", reduction = "umap") + ggtitle("sex")
dev.off()


###########################################################################
# saveRDS(SmithM_multiome.sct, file = "output/seurat/SmithM_multiome.sct.mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4.rds")
SmithM_multiome.sct <- readRDS("output/seurat/SmithM_multiome.sct.mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4.rds") 
###########################################################################
```





### Cell type marker RNA assay

Identify marker genes in each cell types

```R
Idents(SmithM_multiome.sct) <- "seurat_clusters"


# Normalize and scale RNA assay prior DGEs
DefaultAssay(SmithM_multiome.sct) <- "RNA"
SmithM_multiome.sct <- NormalizeData(SmithM_multiome.sct, normalization.method = "LogNormalize", scale.factor = 10000) 
all.genes <- rownames(SmithM_multiome.sct)
SmithM_multiome.sct <- ScaleData(SmithM_multiome.sct, features = all.genes) 

all_markers <- FindAllMarkers(SmithM_multiome.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_markers, file = "output/seurat/FindAllMarkers-onlypospct025fc025-SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4.txt", sep = "\t", quote = FALSE, row.names = TRUE)
```



### Azimuth cell type annotation


```R
library("Seurat")
library("Signac")
library("Azimuth")

# Load seurat object
SmithM_multiome.sct <- readRDS("output/seurat/SmithM_multiome.sct.mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4.rds") 

# Run Azimuth
DefaultAssay(SmithM_multiome.sct) <- "RNA"

SmithM_multiome.sct <- RunAzimuth(SmithM_multiome.sct, reference = "mousecortexref")


pdf("output/seurat/RunAzimuth_mousecortexref-predicted.class-SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-seurat_clusters.pdf", width = 7, height = 6)
DimPlot(SmithM_multiome.sct, group.by = "predicted.class", label = FALSE, label.size = 3, reduction = "umap")
dev.off()
pdf("output/seurat/RunAzimuth_mousecortexref-predicted.subclass-SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-seurat_clusters.pdf", width = 9, height = 6)
DimPlot(SmithM_multiome.sct, group.by = "predicted.subclass", label = TRUE, label.size = 5, reduction = "umap")
dev.off()
pdf("output/seurat/RunAzimuth_mousecortexref-predicted.cluster-SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-seurat_clusters.pdf", width = 9, height = 6)
DimPlot(SmithM_multiome.sct, group.by = "predicted.cluster", label = TRUE, label.size = 5, reduction = "umap")
dev.off()


###########################################################################
# saveRDS(SmithM_multiome.sct, file = "output/seurat/SmithM_multiome.sct.mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4.Azimuth_mousecortexref.rds") 
SmithM_multiome.sct <- readRDS(file = "output/seurat/SmithM_multiome.sct.mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4.Azimuth_mousecortexref.rds")
###########################################################################
```




### ATAC integration and harmony cleaning

```R
##########################################
# Data is not clean so we run harmony ####
##########################################

###  pre-processing and dimensional reductio
DefaultAssay(SmithM_multiome.sct) <- "ATAC"
SmithM_multiome.sct <- RunTFIDF(SmithM_multiome.sct)
SmithM_multiome.sct <- FindTopFeatures(SmithM_multiome.sct, min.cutoff = 'q0') # usually q0 is used, but here ATAC UMAP is fragmented per orig.ident...
SmithM_multiome.sct <- RunSVD(SmithM_multiome.sct)
SmithM_multiome.sct <- RunUMAP(SmithM_multiome.sct, reduction = 'lsi', dims = 2:40, reduction.name = "umap.atac", reduction.key = "atacUMAP_") # We exclude the first dimension as this is typically correlated with sequencing depth

# Harmony
SmithM_multiome.sct <- RunHarmony(
  object = SmithM_multiome.sct,
  group.by.vars = "orig.ident",
  reduction.use = "lsi",
  dims.use = 2:40,
  project.dim = FALSE,
  reduction.save = "harmony_lsi"
)

SmithM_multiome.sct <- RunUMAP(SmithM_multiome.sct, dims = 1:39, reduction = 'harmony_lsi', reduction.name = "umap.atac_harm") # harmony dim 2:40 save as new dims into 1:39
SmithM_multiome.sct <- FindMultiModalNeighbors(SmithM_multiome.sct, reduction.list = list("integrated_dr", "harmony_lsi"), dims.list = list(1:40, 1:39)) # harmony dim 2:40 save as new dims into 1:39
SmithM_multiome.sct <- RunUMAP(SmithM_multiome.sct, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
SmithM_multiome.sct <- FindClusters(SmithM_multiome.sct, graph.name = "wsnn", algorithm = 3, verbose = TRUE)




# plot gene expression, ATAC-seq, or WNN analysis
p1 <- DimPlot(SmithM_multiome.sct, reduction = "umap", group.by = "predicted.subclass", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(SmithM_multiome.sct, reduction = "umap.atac_harm", group.by = "predicted.subclass", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(SmithM_multiome.sct, reduction = "wnn.umap", group.by = "predicted.subclass", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

pdf("output/seurat/UMAP_RNAATACWNN-SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-ATACdim240q0harmony.pdf", width=12, height=5)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

SmithM_multiome.sct$condition <- factor(SmithM_multiome.sct$condition, levels = c("Normoxia", "Hypoxia")) # Reorder untreated 1st

pdf("output/seurat/UMAP_WNN-SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-RNAATACWNN-ATACdim240harmony.pdf", width=5, height=5)
p = DimPlot(SmithM_multiome.sct, reduction = "wnn.umap", group.by = "predicted.subclass", label = FALSE, label.size = 3, repel = TRUE) 
LabelClusters(p, id = "predicted.subclass", fontface = "bold", color = "black", size = 3)
dev.off()
pdf("output/seurat/UMAP_WNN-SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-RNAATACWNN-ATACdim240harmony.pdf", width=5, height=5)
DimPlot(SmithM_multiome.sct, reduction = "wnn.umap", group.by = "condition", label = FALSE, cols = c("blue", "red")) 
dev.off()
pdf("output/seurat/UMAP_WNN-SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-RNAATACWNN-ATACdim240harmony-orig.ident.pdf", width=10, height=5)
DimPlot(SmithM_multiome.sct, reduction = "wnn.umap", group.by = "orig.ident", label = FALSE) 
dev.off()


# Calculate gene activity (count ATAC peak within gene and promoter)
GeneActivity = GeneActivity(
  SmithM_multiome.sct,
  assay = "ATAC",
  features = NULL,
  extend.upstream = 2000,
  extend.downstream = 0,
  biotypes = NULL,
  max.width = NULL,
  process_n = 2000,
  gene.id = FALSE,
  verbose = TRUE
)
## Add gene activity as new assay

SmithM_multiome.sct[["GeneActivity"]] <- CreateAssayObject(counts = GeneActivity)

# SAVE ##########################################################################################
# saveRDS(SmithM_multiome.sct, file = "output/seurat/SmithM_multiome.sct.mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4.Azimuth_mousecortexref.harmony1.GeneActivity.rds")
SmithM_multiome.sct <- readRDS(file = "output/seurat/SmithM_multiome.sct.mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4.Azimuth_mousecortexref.harmony1.GeneActivity.rds")
##########################################################################################
```




### Cluster naming based on Azimuth

Azimuth previsouly provided an annotation confidence score for each cluster, we now rename each of our individual seurat_clusters

```R
############################################
# Re-annotate cluster based on Azimuth #####
############################################
## Extract seurat_clusters from the Azimuth object
az_clusters <- SmithM_multiome.Azimuth.sct$seurat_clusters

## Make sure names are cell names
names(az_clusters) <- colnames(SmithM_multiome.Azimuth.sct)

## Check overlap between objects
length(intersect(colnames(SmithM_multiome.sct), names(az_clusters)))
ncol(SmithM_multiome.sct)

SmithM_multiome.sct$seurat_clusters1 <- az_clusters[colnames(SmithM_multiome.sct)]


meta <- SmithM_multiome.sct@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    seurat_clusters1 = as.character(seurat_clusters1),
    predicted.subclass = as.character(predicted.subclass)
  )
# Count overlap between Seurat clusters and Azimuth subclass
overlap <- meta %>%
  count(seurat_clusters1, predicted.subclass, name = "n") %>%
  group_by(seurat_clusters1) %>%
  mutate(
    cluster_total = sum(n),
    pct_of_cluster = 100 * n / cluster_total
  ) %>%
  ungroup() %>%
  group_by(predicted.subclass) %>%
  mutate(
    subclass_total = sum(n),
    pct_of_subclass = 100 * n / subclass_total
  ) %>%
  ungroup() %>%
  arrange(as.numeric(seurat_clusters1), desc(n))

cluster_annot_map <- c(
  "1"  = "L23_IT1",
  "2"  = "L5_IT1",
  "3"  = "L2356_IT",
  "4"  = "L6_CT1",
  "5"  = "L23_IT2",
  "6"  = "Astro",
  "7"  = "Oligo1",
  "8"  = "Micro_PVM",
  "9"  = "L5_ET",
  "10" = "OPC",
  "11" = "Vip_Lamp5",
  "12" = "Pvalb",
  "13" = "L23_IT3",
  "14" = "L6_IT1",
  "15" = "L23_IT4",
  "16" = "Sst",
  "17" = "L56_NP",
  "18" = "L5_IT2",
  "19" = "Oligo2",
  "20" = "Endo",
  "21" = "L6b",
  "22" = "L6_IT2",
  "23" = "L6_CT2"
)

SmithM_multiome.sct$cluster.annot <- unname(
  cluster_annot_map[as.character(SmithM_multiome.sct$seurat_clusters1)]
)

# plot gene expression, ATAC-seq, or WNN analysis
p1 <- DimPlot(SmithM_multiome.sct, reduction = "umap", group.by = "cluster.annot", label = TRUE, label.size = 3.25, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(SmithM_multiome.sct, reduction = "umap.atac_harm", group.by = "cluster.annot", label = TRUE, label.size =3.25, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(SmithM_multiome.sct, reduction = "wnn.umap", group.by = "cluster.annot", label = TRUE, label.size = 3.25, repel = TRUE) + ggtitle("WNN")

pdf("output/seurat/UMAP_RNAATACWNN-SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-ATACdim240q0harmonyclusterannot.pdf", width=12, height=5)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()


pdf("output/seurat/UMAP_WNN-SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-RNAATACWNN-ATACdim240harmonyclusterannot.pdf", width=5, height=5)
DimPlot(SmithM_multiome.sct, reduction = "wnn.umap", group.by = "condition", label = FALSE, cols = c("blue", "red")) 
dev.off()
```




### DGEs analysis - RNA assay


Identify DEGs (Normoxia vs Hypoxia) in each cell types with *Wilcox* and *MAST* methods


```R
####################################
# WILCOX DEG - cluster.annot #####
####################################

# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
DefaultAssay(SmithM_multiome.sct) <- "RNA"


SmithM_multiome.sct <- NormalizeData(SmithM_multiome.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(SmithM_multiome.sct)
SmithM_multiome.sct <- ScaleData(SmithM_multiome.sct, features = all.genes) 

SmithM_multiome.sct$annot.stim <- paste(SmithM_multiome.sct$cluster.annot, SmithM_multiome.sct$condition,
    sep = "-")
Idents(SmithM_multiome.sct) <- "annot.stim"


### Here is to automatize the process cluster per cluster:
#### Define the cluster pairs for comparison
cluster_pairs <- list(
  c("L23_IT1-Normoxia",   "L23_IT1-Hypoxia",   "L23_IT1"),
  c("L5_IT1-Normoxia",    "L5_IT1-Hypoxia",    "L5_IT1"),
  c("L2356_IT-Normoxia",   "L2356_IT-Hypoxia",   "L2356_IT"),
  c("L6_CT1-Normoxia",   "L6_CT1-Hypoxia",   "L6_CT1"),
  c("L23_IT2-Normoxia",   "L23_IT2-Hypoxia",   "L23_IT2"),
  c("Astro-Normoxia",   "Astro-Hypoxia",   "Astro"),
  c("Oligo1-Normoxia",   "Oligo1-Hypoxia",   "Oligo1"),
  c("Micro_PVM-Normoxia",   "Micro_PVM-Hypoxia",   "Micro_PVM"),
  c("L5_ET-Normoxia",     "L5_ET-Hypoxia",     "L5_ET"),
  c("OPC-Normoxia",     "OPC-Hypoxia",     "OPC"),
  c("Vip_Lamp5-Normoxia",     "Vip_Lamp5-Hypoxia",     "Vip_Lamp5"),
  c("Pvalb-Normoxia",     "Pvalb-Hypoxia",     "Pvalb"),
  c("L23_IT3-Normoxia",     "L23_IT3-Hypoxia",     "L23_IT3"),
  c("L6_IT1-Normoxia",      "L6_IT1-Hypoxia",      "L6_IT1"),
  c("L23_IT4-Normoxia",  "L23_IT4-Hypoxia",  "L23_IT4"),
  c("Sst-Normoxia",     "Sst-Hypoxia",     "Sst"),
  c("L56_NP-Normoxia",     "L56_NP-Hypoxia",     "L56_NP"),
  c("L5_IT2-Normoxia",     "L5_IT2-Hypoxia",     "L5_IT2"),
  c("Oligo2-Normoxia",      "Oligo2-Hypoxia",      "Oligo2"),
  c("Endo-Normoxia",      "Endo-Hypoxia",      "Endo"),
  c("L6b-Normoxia",     "L6b-Hypoxia",     "L6b"),
  c("L6_IT2-Normoxia",      "L6_IT2-Hypoxia",      "L6_IT2"),
  c("L6_CT2-Normoxia",      "L6_CT2-Hypoxia",      "L6_CT2")
)

## Function to run FindMarkers and return a tibble with cluster and gene information
run_find_markers <- function(pair) {
  control <- pair[1]      # Normoxia
  treatment <- pair[2]    # Hypoxia
  subclass_name <- pair[3]
  
  markers <- FindMarkers(
    object = SmithM_multiome.sct,
    ident.1 = treatment,
    ident.2 = control,
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    min.diff.pct = -Inf,
    assay = "RNA"
  ) %>%
    rownames_to_column(var = "gene") %>%
    mutate(subclass = subclass_name) %>%
    relocate(subclass, gene)
  return(markers)
}

FindMarkers_NormoxiavsHypoxia_clusterannot <- map_dfr(
  cluster_pairs,
  run_find_markers
) %>%
  as_tibble()

## save output
write.table(FindMarkers_NormoxiavsHypoxia_clusterannot, file = "output/seurat/FindMarkers-SmithM_multiome-NormoxiavsHypoxia_subclass-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-clusterannot.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)






####################################
# MAST DEG - cluster.annot #####
####################################
# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
DefaultAssay(SmithM_multiome.sct) <- "RNA"


SmithM_multiome.sct <- NormalizeData(SmithM_multiome.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(SmithM_multiome.sct)
SmithM_multiome.sct <- ScaleData(SmithM_multiome.sct, features = all.genes) 

SmithM_multiome.sct$annot.stim <- paste(SmithM_multiome.sct$cluster.annot, SmithM_multiome.sct$condition,
    sep = "-")
Idents(SmithM_multiome.sct) <- "annot.stim"

### Here is to automatize the process cluster per cluster:
#### Define the cluster pairs for comparison
cluster_pairs <- list(
  c("L23_IT1-Normoxia",   "L23_IT1-Hypoxia",   "L23_IT1"),
  c("L5_IT1-Normoxia",    "L5_IT1-Hypoxia",    "L5_IT1"),
  c("L2356_IT-Normoxia",   "L2356_IT-Hypoxia",   "L2356_IT"),
  c("L6_CT1-Normoxia",   "L6_CT1-Hypoxia",   "L6_CT1"),
  c("L23_IT2-Normoxia",   "L23_IT2-Hypoxia",   "L23_IT2"),
  c("Astro-Normoxia",   "Astro-Hypoxia",   "Astro"),
  c("Oligo1-Normoxia",   "Oligo1-Hypoxia",   "Oligo1"),
  c("Micro_PVM-Normoxia",   "Micro_PVM-Hypoxia",   "Micro_PVM"),
  c("L5_ET-Normoxia",     "L5_ET-Hypoxia",     "L5_ET"),
  c("OPC-Normoxia",     "OPC-Hypoxia",     "OPC"),
  c("Vip_Lamp5-Normoxia",     "Vip_Lamp5-Hypoxia",     "Vip_Lamp5"),
  c("Pvalb-Normoxia",     "Pvalb-Hypoxia",     "Pvalb"),
  c("L23_IT3-Normoxia",     "L23_IT3-Hypoxia",     "L23_IT3"),
  c("L6_IT1-Normoxia",      "L6_IT1-Hypoxia",      "L6_IT1"),
  c("L23_IT4-Normoxia",  "L23_IT4-Hypoxia",  "L23_IT4"),
  c("Sst-Normoxia",     "Sst-Hypoxia",     "Sst"),
  c("L56_NP-Normoxia",     "L56_NP-Hypoxia",     "L56_NP"),
  c("L5_IT2-Normoxia",     "L5_IT2-Hypoxia",     "L5_IT2"),
  c("Oligo2-Normoxia",      "Oligo2-Hypoxia",      "Oligo2"),
  c("Endo-Normoxia",      "Endo-Hypoxia",      "Endo"),
  c("L6b-Normoxia",     "L6b-Hypoxia",     "L6b"),
  c("L6_IT2-Normoxia",      "L6_IT2-Hypoxia",      "L6_IT2"),
  c("L6_CT2-Normoxia",      "L6_CT2-Hypoxia",      "L6_CT2")
)

## Function to run FindMarkers and return a tibble with cluster and gene information
run_find_markers <- function(pair) {
  control <- pair[1]      # Normoxia
  treatment <- pair[2]    # Hypoxia
  subclass_name <- pair[3]
  
  markers <- FindMarkers(
    object = SmithM_multiome.sct,
    ident.1 = treatment,
    ident.2 = control,
    test.use = "MAST",
    logfc.threshold = 0,
    min.pct = 0,
    min.diff.pct = -Inf,
    assay = "RNA", # Specify the RNA assay (default for raw counts)
    slot = "data") %>% # Use lognorm data for MAST
    rownames_to_column(var = "gene") %>%
    mutate(subclass = subclass_name) %>%
    relocate(subclass, gene)
  return(markers)
}

FindMarkers_MAST_NormoxiavsHypoxia_clusterannot <- map_dfr(
  cluster_pairs,
  run_find_markers
) %>%
  as_tibble()

## save output
write.table(FindMarkers_MAST_NormoxiavsHypoxia_clusterannot, file = "output/seurat/FindMarkers_MAST-SmithM_multiome-NormoxiavsHypoxia_subclass-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-clusterannot.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# SAVE ##########################################################################################
# saveRDS(SmithM_multiome.sct, file = "output/seurat/SmithM_multiome.sct.mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4.Azimuth_mousecortexref.harmony.clusterannot.rds")
SmithM_multiome.sct <- readRDS(file = "output/seurat/SmithM_multiome.sct.mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4.Azimuth_mousecortexref.harmony.clusterannot.rds")
##########################################################################################
```


### DARs analysis with Signac - ATAC assay

DAR Normoxia vs Hypoxia in each cluster

```R
# Identify diff access regions (DAR) - LR stat method + latent.vars ("nCount_ATAC")
## PRIOR Lets switch to ATAC assay and normalize and scale before doing the DEGs
DefaultAssay(SmithM_multiome.sct) <- "ATAC"

#--> FOR ATAC do NOT NormalizeData() and ScaleData()

SmithM_multiome.sct$annot.stim <- paste(SmithM_multiome.sct$cluster.annot, SmithM_multiome.sct$condition,
    sep = "-")
Idents(SmithM_multiome.sct) <- "annot.stim"


### Here is to automatize the process clsuter per cluster:
#### Define the cluster pairs for comparison
cluster_pairs <- list(
  c("L23_IT1-Normoxia",   "L23_IT1-Hypoxia",   "L23_IT1"),
  c("L5_IT1-Normoxia",    "L5_IT1-Hypoxia",    "L5_IT1"),
  c("L2356_IT-Normoxia",   "L2356_IT-Hypoxia",   "L2356_IT"),
  c("L6_CT1-Normoxia",   "L6_CT1-Hypoxia",   "L6_CT1"),
  c("L23_IT2-Normoxia",   "L23_IT2-Hypoxia",   "L23_IT2"),
  c("Astro-Normoxia",   "Astro-Hypoxia",   "Astro"),
  c("Oligo1-Normoxia",   "Oligo1-Hypoxia",   "Oligo1"),
  c("Micro_PVM-Normoxia",   "Micro_PVM-Hypoxia",   "Micro_PVM"),
  c("L5_ET-Normoxia",     "L5_ET-Hypoxia",     "L5_ET"),
  c("OPC-Normoxia",     "OPC-Hypoxia",     "OPC"),
  c("Vip_Lamp5-Normoxia",     "Vip_Lamp5-Hypoxia",     "Vip_Lamp5"),
  c("Pvalb-Normoxia",     "Pvalb-Hypoxia",     "Pvalb"),
  c("L23_IT3-Normoxia",     "L23_IT3-Hypoxia",     "L23_IT3"),
  c("L6_IT1-Normoxia",      "L6_IT1-Hypoxia",      "L6_IT1"),
  c("L23_IT4-Normoxia",  "L23_IT4-Hypoxia",  "L23_IT4"),
  c("Sst-Normoxia",     "Sst-Hypoxia",     "Sst"),
  c("L56_NP-Normoxia",     "L56_NP-Hypoxia",     "L56_NP"),
  c("L5_IT2-Normoxia",     "L5_IT2-Hypoxia",     "L5_IT2"),
  c("Oligo2-Normoxia",      "Oligo2-Hypoxia",      "Oligo2"),
  c("Endo-Normoxia",      "Endo-Hypoxia",      "Endo"),
  c("L6b-Normoxia",     "L6b-Hypoxia",     "L6b"),
  c("L6_IT2-Normoxia",      "L6_IT2-Hypoxia",      "L6_IT2"),
  c("L6_CT2-Normoxia",      "L6_CT2-Hypoxia",      "L6_CT2")
)

## Function to run FindMarkers and return a tibble with cluster and gene information
run_find_markers <- function(pair) {
  control <- pair[1]      # Normoxia
  treatment <- pair[2]    # Hypoxia
  subclass_name <- pair[3]
  
  markers <- FindMarkers(
    object = SmithM_multiome.sct,
    ident.1 = treatment,
    ident.2 = control,
    test.use = "LR",
    min.pct = 0.1, # 
    assay = "ATAC",
  latent.vars = c("nCount_ATAC" ) 
  ) %>%
    rownames_to_column(var = "query_region") %>%
    mutate(subclass = subclass_name) %>%
    relocate(subclass, query_region)
  return(markers)
}

FindMarkers_ATAC_NormoxiavsHypoxia_clusterannot <- map_dfr(
  cluster_pairs,
  run_find_markers
) %>%
  as_tibble()

## save output
write.table(FindMarkers_ATAC_NormoxiavsHypoxia_clusterannot, file = "output/seurat/FindMarkers_ATAC-SmithM_multiome-NormoxiavsHypoxia_subclass-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4-clusterannot.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

## Find closest gene to DAR
DefaultAssay(SmithM_multiome.sct) <- "ATAC"

DAR_closestGene <- ClosestFeature(SmithM_multiome.sct, regions = FindMarkers_ATAC_NormoxiavsHypoxia_clusterannot$query_region)
DAR_genes = FindMarkers_ATAC_NormoxiavsHypoxia_clusterannot %>% 
  left_join(DAR_closestGene) %>%
  as_tibble() %>%
  unique()
 
## save output
write.table(DAR_genes, file = "output/seurat/DAR_genes-SmithM_multiome-LR_latentnCount_ATAC-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
```






### TSS plot

Plot of ATAC signal around gene TSS

```R
##################################
## TSSplot - All genes #################
##################################
DefaultAssay(SmithM_multiome.sct) <- "ATAC"

SmithM_multiome.sct = TSSEnrichment(SmithM_multiome.sct, assay = "ATAC", fast = FALSE, verbose = TRUE)

# all genes - all cluster
pdf("output/seurat/TSSPlot-SmithM_multiome-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4.pdf", width = 6, height = 5)
TSSPlot(
  SmithM_multiome.sct,
  assay = "ATAC",
  group.by = "condition"
) +
  facet_null() +
  scale_color_manual(values = c(
    "Normoxia" = "#2b05ff",      # 
    "Hypoxia" = "#ff0000"   # 
  )) +
  theme_bw()
dev.off()


# One plot per cell type ### 
clusters_all <- unique(SmithM_multiome.sct$cluster.annot)

for (clust in clusters_all) {
  obj_sub <- subset(
    SmithM_multiome.sct,
    subset = cluster.annot == clust
  )
  p <- TSSPlot(
    obj_sub,
    assay = "ATAC",
    group.by = "condition"
  ) +
    facet_null() +
    scale_color_manual(values = c(
      "Normoxia" = "#2b05ff",      # 
      "Hypoxia" = "#ff0000"   # 
    )) +
    theme_bw() +
    ggtitle(clust) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  file_name <- paste0(
    "output/seurat/TSSPlot-SmithM_multiome-",
    gsub("[ /]", "_", clust),
    "-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4.pdf"
  )
  pdf(file_name, width = 6, height = 5)
  print(p)
  dev.off()
}


##################################
## TSSplot - DAR genes #################
##################################
DefaultAssay(SmithM_multiome.sct) <- "ATAC"

## 1. Extract significant DAR genes: Open + Close together
DAR_sig <- DAR_genes %>%
  dplyr::filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 0.25
  ) %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      avg_log2FC > 0.25  ~ "Open",
      avg_log2FC < -0.25 ~ "Close",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Direction))


all_DAR_genes <- DAR_sig %>%
  dplyr::pull(gene_name) %>%
  unique() %>%
  na.omit() %>%
  as.character()


## 2. Get TSS positions for only those DAR genes
annotations <- Annotation(SmithM_multiome.sct[["ATAC"]])
tss_all <- GetTSSPositions(ranges = annotations)

# Keep only TSS where gene_name is one of your DAR genes
tss_DAR_all <- tss_all[
  toupper(mcols(tss_all)$gene_name) %in% toupper(all_DAR_genes)
]

# Check result
length(tss_DAR_all)
tss_DAR_all


## 3. Global TSSPlot using only Open + Close DAR genes
SmithM_multiome.sct_DAR_TSS <- TSSEnrichment(
  object = SmithM_multiome.sct,
  assay = "ATAC",
  tss.positions = tss_DAR_all,
  fast = FALSE,
  process_n = 2000,
  verbose = TRUE
)
p_global_DAR_TSS <- TSSPlot(
  SmithM_multiome.sct_DAR_TSS,
  assay = "ATAC",
  group.by = "condition"
) +
  facet_null() +
  scale_color_manual(values = c(
    "Normoxia" = "#2b05ff",
    "Hypoxia"  = "#ff0000"
  )) +
  theme_bw() +
  ggtitle("All cell types - DAR genes only") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

pdf(
  "output/seurat/TSSPlot-SmithM_multiome-DARgenes_OpenClose_ALLcelltypes-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4.pdf",
  width = 6,
  height = 5
)
print(p_global_DAR_TSS)
dev.off()


## One TSSPlot per cell type
## using DAR genes specific to that cell type


safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

clusters_all <- unique(SmithM_multiome.sct$cluster.annot)

for (clust in clusters_all) {
  
  message("Processing: ", clust)
    ## 1. Get DAR genes specific to this cell type
  DAR_genes_clust <- DAR_sig %>%
    dplyr::filter(subclass == clust) %>%
    dplyr::pull(gene_name) %>%
    unique() %>%
    na.omit() %>%
    as.character()
  
  if (length(DAR_genes_clust) == 0) {
    message("  Skipping ", clust, ": no DAR genes.")
    next
  }
  ## 2. Isolate TSSs of those DAR genes only  
  tss_DAR_clust <- tss_all[
    toupper(mcols(tss_all)$gene_name) %in% toupper(DAR_genes_clust)
  ]
  
  tss_DAR_clust <- unique(tss_DAR_clust)
  
  if (length(tss_DAR_clust) == 0) {
    message("  Skipping ", clust, ": no matching TSS found.")
    next
  }
  
  message("  DAR genes: ", length(DAR_genes_clust))
  message("  TSSs found: ", length(tss_DAR_clust))
    ## 3. Subset object to this cell type
  obj_sub <- subset(
    SmithM_multiome.sct,
    subset = cluster.annot == clust
  )
  ## 4. Recalculate TSS enrichment using only this cell type's DAR-gene TSSs
  obj_sub <- TSSEnrichment(
    object = obj_sub,
    assay = "ATAC",
    tss.positions = tss_DAR_clust,
    fast = FALSE,
    process_n = 2000,
    verbose = TRUE
  )
  ## 5. Plot Normoxia vs Hypoxia  
  p <- TSSPlot(
    obj_sub,
    assay = "ATAC",
    group.by = "condition"
  ) +
    facet_null() +
    scale_color_manual(values = c(
      "Normoxia" = "#2b05ff",
      "Hypoxia"  = "#ff0000"
    )) +
    theme_bw() +
    ggtitle(paste0(clust, " - cell_type_specificDARgenes")) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  file_name <- paste0(
    "output/seurat/TSSPlot-SmithM_multiome-DARgenes_OpenClose_CELLTYPE_SPECIFIC-",
    safe_name(clust),
    "-mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4.pdf"
  )
  
  pdf(file_name, width = 6, height = 5)
  print(p)
  dev.off()
}
```




### SCPA - functional analysis - RNA assay

GO and pathway enrichment analysis with [SCPA](https://jackbibby1.github.io/SCPA/)


```R
library("Seurat")
library("Signac")
library("tidyverse")
library("SCPA")
library("msigdbr")
library("dplyr")
library("ggplot2")
library("stringr")


SmithM_multiome.sct <- readRDS(file = "output/seurat/SmithM_multiome.sct.mtrbNotRegress_vKeepBottom_dim50kparam50res027algo4.Azimuth_mousecortexref.harmony.clusterannot.rds")

DefaultAssay(SmithM_multiome.sct) <- "RNA" # 

## list available pathways
msigdbr_collections(db_species = "MM")
## import Pathways
pathways <- msigdbr(
  db_species = "MM",
  species = "Mus musculus",
  collection = "M5",
  subcollection = "GO:BP",
) %>%
  format_pathways()

pathways <- msigdbr(
  db_species = "MM",
  species = "Mus musculus",
  collection = "M2"
) %>%
  format_pathways()
  
names(pathways) <- sapply(pathways, function(x) x$Pathway[1]) # just to name the list, so easier to visualise

# Code to save output for each cell type comparison
clusters = c(
    # Glutamatergic / excitatory neurons
  "L23_IT1",
  "L23_IT2",
  "L23_IT3",
  "L23_IT4",
  "L2356_IT",
  "L5_IT1",
  "L5_IT2",
  "L5_ET",
  "L56_NP",
  "L6_IT1",
  "L6_IT2",
  "L6_CT1",
  "L6_CT2",
  "L6b",
    # GABAergic / inhibitory neurons
  "Vip_Lamp5",
  "Sst",
  "Pvalb",
    # Non-neuronal
  "Astro",
  "Oligo1",
  "Oligo2",
  "OPC",
  "Micro_PVM",
  "Endo"
)


### Loop through each value
for (cluster in clusters) {
  #### Extract data for Normoxia and Hypoxia based on current value
  Normoxia <- seurat_extract(SmithM_multiome.sct,
                       meta1 = "condition", value_meta1 = "Normoxia",
                       meta2 = "cluster.annot", value_meta2 = cluster)

  Hypoxia <- seurat_extract(SmithM_multiome.sct,
                           meta1 = "condition", value_meta1 = "Hypoxia",
                           meta2 = "cluster.annot", value_meta2 = cluster)

  ##### Compare pathways
  Hypoxia_Normoxia <- compare_pathways(samples = list(Hypoxia, Normoxia),  # list(population1,population2) FC = population 2 - population 1; here FC >1 = more in Hypoxia ; <0 = less in hypoxia
                                pathways = pathways,
                                parallel = TRUE, cores = 8)

  ##### Write to file using the current value in the filename
  output_filename <- paste0("output/Pathway/SCPA_M5_GOBP-SmithM_multiome-dim50kparam50res027algo4clusterannot-", cluster, ".txt")       # CHANGE HERE PATHWAYS NAME 
  write.table(Hypoxia_Normoxia, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}

####################
# DATA ANALYSIS  ###
####################
# Cell types in the same order used for the SCPA analysis
clusters <- c(
  # Glutamatergic / excitatory neurons
  "L23_IT1",
  "L23_IT2",
  "L23_IT3",
  "L23_IT4",
  "L2356_IT",
  "L5_IT1",
  "L5_IT2",
  "L5_ET",
  "L56_NP",
  "L6_IT1",
  "L6_IT2",
  "L6_CT1",
  "L6_CT2",
  "L6b",

  # GABAergic / inhibitory neurons
  "Vip_Lamp5",
  "Sst",
  "Pvalb",

  # Non-neuronal
  "Astro",
  "Oligo1",
  "Oligo2",
  "OPC",
  "Micro_PVM",
  "Endo"
)


#############################################
# Import all SCPA results into one tibble ###
#############################################

#  SCPA_M5_GOBP
scpa_results <- map_dfr(clusters, function(cluster) {

  input_filename <- file.path(
    "output/Pathway",
    paste0(
      "SCPA_M5_GOBP-SmithM_multiome-",  # SCPA_M5_GOBP
      "dim50kparam50res027algo4clusterannot-",
      cluster,
      ".txt"
    )
  )

  if (!file.exists(input_filename)) {
    warning("File not found for cluster: ", cluster)
    return(NULL)
  }

  read_tsv(
    input_filename,
    show_col_types = FALSE
  ) %>%
    mutate(cluster = cluster, .before = 1)
})

# Inspect results
head(scpa_results)
# dotplot summary
# Keep significant pathways
scpa_sig <- scpa_results %>%
  filter(
    adjPval < 0.05,
    is.finite(FC),
    is.finite(qval)
  ) %>%
  mutate(
    cluster = factor(cluster, levels = clusters),

    # Remove the common GOBP_ prefix to make labels shorter
    Pathway_short = str_remove(Pathway, "^GOBP_")
  )

# Order pathways first by the cell type in which they have the
# strongest qval, and then by qval
pathway_order <- scpa_sig %>%
  group_by(Pathway_short) %>%
  slice_max(
    order_by = qval,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(cluster, qval) %>%
  pull(Pathway_short)

scpa_sig <- scpa_sig %>%
  mutate(
    Pathway_short = factor(
      Pathway_short,
      levels = pathway_order
    )
  )
pdf("output/Pathway/dotplot-SCPA_M5_GOBP-sig-allClusters.pdf", width = 7, height = 8)
ggplot(
  scpa_sig,
  aes(
    x = cluster,
    y = Pathway_short
  )
) +
  geom_point(
    aes(
      size = qval,
      fill = FC
    ),
    shape = 21,
    color = "black",
    stroke = 0.25,
    alpha = 0.9
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "FC"
  ) +
  scale_size_continuous(
    range = c(2.5, 8),
    name = "SCPA qval"
  ) +
  labs(
    title = "Significant SCPA pathways across cell types",
    subtitle = "Only pathways with adjusted P-value < 0.05",
    x = "Cell type",
    y = "GO biological process"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 7
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      linewidth = 0.2,
      color = "grey90"
    ),
    legend.position = "right"
  )
dev.off()



#  SCPA_M2
scpa_results <- map_dfr(clusters, function(cluster) {

  input_filename <- file.path(
    "output/Pathway",
    paste0(
      "SCPA_M2-SmithM_multiome-",  # 
      "dim50kparam50res027algo4clusterannot-",
      cluster,
      ".txt"
    )
  )

  if (!file.exists(input_filename)) {
    warning("File not found for cluster: ", cluster)
    return(NULL)
  }

  read_tsv(
    input_filename,
    show_col_types = FALSE
  ) %>%
    mutate(cluster = cluster, .before = 1)
})

# Inspect results
head(scpa_results)

# Keep significant pathways
scpa_sig <- scpa_results %>%
  filter(
    adjPval < 0.05,
    is.finite(FC),
    is.finite(qval)
  ) %>%
  mutate(
    cluster = factor(cluster, levels = clusters),

    # Remove the common GOBP_ prefix to make labels shorter
    Pathway_short = str_remove(Pathway, "^GOBP_")
  )

# Order pathways first by the cell type in which they have the
# strongest qval, and then by qval
pathway_order <- scpa_sig %>%
  group_by(Pathway_short) %>%
  slice_max(
    order_by = qval,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(cluster, qval) %>%
  pull(Pathway_short)

scpa_sig <- scpa_sig %>%
  mutate(
    Pathway_short = factor(
      Pathway_short,
      levels = pathway_order
    )
  )

pdf("output/Pathway/dotplot-SCPA_M2-sig-allClusters.pdf", width = 8, height = 8)
ggplot(
  scpa_sig,
  aes(
    x = cluster,
    y = Pathway_short
  )
) +
  geom_point(
    aes(
      size = qval,
      fill = FC
    ),
    shape = 21,
    color = "black",
    stroke = 0.25,
    alpha = 0.9
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "FC"
  ) +
  scale_size_continuous(
    range = c(2.5, 8),
    name = "SCPA qval"
  ) +
  labs(
    title = "Significant SCPA pathways across cell types",
    subtitle = "Only pathways with adjusted P-value < 0.05",
    x = "Cell type",
    y = "GO biological process"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 7
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      linewidth = 0.2,
      color = "grey90"
    ),
    legend.position = "right"
  )
dev.off()
```




