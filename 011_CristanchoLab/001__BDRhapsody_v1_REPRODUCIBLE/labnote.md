--> This folder is the **clean + reproducible** version of: `001_CristanchoLab/001__BDRhapsody_v1` (original exploration)

# Data Overview

- **Data type**: BD Rhapsody multiome - Whole Transcriptome Analysis (WTA) RNA-seq and ATAC-seq
- **Factors**:
  - **Organism**: Mice
  - **Condition**: Normoxia (24h), Hypoxia (24h)
  - **Replicates**: 2 biological replicates per group (total n = 4)



# Data processing

## Data import

Data downloaded from [Seven Bridge - BD Rhapsody](https://igor.sbgenomics.com/home)

- Seurat object : `/mnt/isilon/cristancho_data/Thomas/011_CristanchoLab/001__BDRhapsody_v1/output/seurat/ATACMultiomewithST-CC-additionalSMK_Seurat.rds`
- ATAC fragment file: `/mnt/isilon/cristancho_data/Thomas/011_CristanchoLab/001__BDRhapsody_v1/output/seurat/ATACMultiomewithST-CC-additionalSMK_ATAC_Fragments.bed.gz`


## Data analysis

*NOTE: BD Rhpasody data already identified doublet so no need to run scrublet*

### QC filtering - RNA assay 

```R
set.seed(42)

# library/packages
library("Signac")
library("Seurat")
library("tidyverse")
library("reticulate") 
library("metap") 
use_python("~/anaconda3/envs/SignacV5/bin/python") 


# load seurat object
ATACMultiomewithST_SMK <- readRDS(file = "output/seurat/ATACMultiomewithST-CC-additionalSMK_Seurat.rds")

# Add mit and rib genes information
ATACMultiomewithST_SMK[["percent.mt"]] <- PercentageFeatureSet(ATACMultiomewithST_SMK, pattern = "^mt-")
ATACMultiomewithST_SMK[["percent.rb"]] <- PercentageFeatureSet(ATACMultiomewithST_SMK, pattern = "^Rp[sl]")


# QC plots
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK-nFeature_RNA.pdf", width = 10, height = 6)
VlnPlot(ATACMultiomewithST_SMK, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK-nCount_RNA.pdf", width = 10, height = 6)
VlnPlot(ATACMultiomewithST_SMK, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK-percent.mt.pdf", width = 10, height = 6)
VlnPlot(ATACMultiomewithST_SMK, features = c("percent.mt"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK-percent.rb.pdf", width = 10, height = 6)
VlnPlot(ATACMultiomewithST_SMK, features = c("percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
#--> Manual inspection of QC plots to identify QC treshold


# Subset cells passing QC
ATACMultiomewithST_SMK_QCv1 <- subset(
  ATACMultiomewithST_SMK,
  subset =
    nFeature_RNA > 250 &
    nFeature_RNA < 2000 &
    nCount_RNA > 250 &
    nCount_RNA < 5000 &
    percent.mt < 7.5 &
    percent.rb < 5 &
    Sample_Name != "Multiplet" &
    Sample_Name != "Undetermined"
)
#--> 27,173 to 21,685 cells


```

After manual inspection of QC plots, we remove the outlier cells and multiplet; we go from 27,173 (raw) to 21,685 (RNA clean) cells .


### QC filtering - ATAC assay


```R


# Path to fragments file
fragpath <- "output/seurat/ATACMultiomewithST-CC-additionalSMK_ATAC_Fragments.bed.gz"


# Add fragment file to seurat object
ATACMultiomewithST_SMK_QCv1[['peaks']]@fragments <- list(
  CreateFragmentObject(
    path = fragpath,
    cellnames = colnames(ATACMultiomewithST_SMK_QCv1),
    validate.fragments = TRUE
  )
)

DefaultAssay(ATACMultiomewithST_SMK_QCv1) <- "peaks" # 

# compute nucleosome signal score per cell
ATACMultiomewithST_SMK_QCv1 <- NucleosomeSignal(object = ATACMultiomewithST_SMK_QCv1)
# compute TSS enrichment score per cell
ATACMultiomewithST_SMK_QCv1 <- TSSEnrichment(object = ATACMultiomewithST_SMK_QCv1, fast = FALSE)



# QC plots
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK_V2QCv1-nCount_peaks.pdf", width=12, height=6)
VlnPlot(
  object = ATACMultiomewithST_SMK_QCv1,
  features = c('nCount_peaks'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK_V2QCv1-TSSenrichment.pdf", width=12, height=6)
VlnPlot(
  object = ATACMultiomewithST_SMK_QCv1,
  features = c('TSS.enrichment'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK_V2QCv1-nucleosome_signal.pdf", width=12, height=6)
VlnPlot(
  object = ATACMultiomewithST_SMK_QCv1,
  features = c('nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
ATACMultiomewithST_SMK_QCv1$high.tss <- ifelse(ATACMultiomewithST_SMK_QCv1$TSS.enrichment > 2, 'High', 'Low') 
pdf("output/seurat/TSSPlot-ATACMultiomewithST_SMK_QCv1-TSSenrichment2.pdf", width=12, height=6)
TSSPlot(ATACMultiomewithST_SMK_QCv1, group.by = 'high.tss') + NoLegend()
dev.off()


# Subset cells passing QC
ATACMultiomewithST_SMK_QCv2 <- subset(
  x = ATACMultiomewithST_SMK_QCv1,
  subset = nCount_peaks > 100 &
    nCount_peaks < 25000 &
    TSS.enrichment > 2 &
    TSS.enrichment < 11 &
    nucleosome_signal > 0.11 &
    nucleosome_signal < 0.28
)
#--> 21,685 to 21,333 cells
```

After manual inspection of QC plots, we remove the outlier cells and multiplet; we go from 21,685 (RNA clean cells) to 21,333 cells.


### Sample integration with SCT

Integrate each individual samples with SCT



```R
# Individualize sample
Normoxia_F <- subset(ATACMultiomewithST_SMK_QCv2, subset = Sample_Name == "Normoxia_F") # 4488 cells
Normoxia_M <- subset(ATACMultiomewithST_SMK_QCv2, subset = Sample_Name == "Normoxia_M") # 5178 cells

Hypoxia_F <- subset(ATACMultiomewithST_SMK_QCv2, subset = Sample_Name == "Hypoxia_F") # 6881 cells
Hypoxia_M <- subset(ATACMultiomewithST_SMK_QCv2, subset = Sample_Name == "Hypoxia_M") # 4786 cells


# Add sex/condition as metrics
Normoxia_F$sex <- "Female"
Normoxia_M$sex <- "Male"
Hypoxia_F$sex <- "Female"
Hypoxia_M$sex <- "Male"
Normoxia_F$condition <- "Normoxia"
Normoxia_M$condition <- "Normoxia"
Hypoxia_F$condition <- "Hypoxia"
Hypoxia_M$condition <- "Hypoxia"

DefaultAssay(Normoxia_F) <- "RNA"
DefaultAssay(Normoxia_M) <- "RNA"
DefaultAssay(Hypoxia_F) <- "RNA"
DefaultAssay(Hypoxia_M) <- "RNA"


# SCTransform each sample individually
Normoxia_F <- SCTransform(Normoxia_F, method = "glmGamPoi", ncells = 4488, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "nFeature_RNA")) 
Normoxia_M <- SCTransform(Normoxia_M, method = "glmGamPoi", ncells = 5178, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "nFeature_RNA")) 
Hypoxia_F <- SCTransform(Hypoxia_F, method = "glmGamPoi", ncells = 6881, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "nFeature_RNA")) 
Hypoxia_M <- SCTransform(Hypoxia_M, method = "glmGamPoi", ncells = 4786, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "nFeature_RNA")) 


# SCT integration
srat.list <- list(Normoxia_F = Normoxia_F, Normoxia_M = Normoxia_M, Hypoxia_F = Hypoxia_F, Hypoxia_M = Hypoxia_M)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
ATACMultiomewithST_SMK_V2QCv2.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)

ATACMultiomewithST_SMK_V2QCv2.sct <- IntegrateData(anchorset = ATACMultiomewithST_SMK_V2QCv2.anchors, normalization.method = "SCT")


DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- "integrated"

ATACMultiomewithST_SMK_V2QCv2.sct <- RunPCA(ATACMultiomewithST_SMK_V2QCv2.sct, verbose = FALSE, npcs = 30)
ATACMultiomewithST_SMK_V2QCv2.sct <- RunUMAP(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
ATACMultiomewithST_SMK_V2QCv2.sct <- FindNeighbors(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "pca", k.param = 30, dims = 1:30)
ATACMultiomewithST_SMK_V2QCv2.sct <- FindClusters(ATACMultiomewithST_SMK_V2QCv2.sct, resolution = 0.4, verbose = FALSE, algorithm = 4) # method = "igraph" needed for large nb of cells


ATACMultiomewithST_SMK_V2QCv2.sct$condition <- factor(ATACMultiomewithST_SMK_V2QCv2.sct$condition, levels = c("Normoxia", "Hypoxia")) # 

# UMAP plots
pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.pdf", width=7, height=6)
DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap", label=TRUE)
dev.off()
pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30-groupCondition.pdf", width=7, height=6)
DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap", label=FALSE, group.by= "condition" )
dev.off()


###########################################################################
# saveRDS(ATACMultiomewithST_SMK_V2QCv2.sct, file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.rds") 
ATACMultiomewithST_SMK_V2QCv2.sct <- readRDS(file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.rds")
###########################################################################
```

Seurat object at this stage can be accessed at: `output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.rds`



### Cluster naming

Cluster were annotated based on below marker genes expression:
- Cluster1= Glut1= Neurod2, Neurod6, Rbfox1
- Cluster2= Glut2= Neurod2, Neurod6, Rbfox1
- Cluster3= Glut3= Neurod2, Neurod6, Rbfox1
- Cluster4= Glut4= Neurod2, Neurod6, Rbfox1
- Cluster5= NSC= Mki67, Top2a, Dlx6os1
- Cluster6= RG= Fabp7, Slc1a3, Tnc
- Cluster7= In1= Gad2, Nrxn3, Lhx6
- Cluster8= In2= Gad2, Nrxn3, Lhx6
- Cluster9= In3= Gad2, Nrxn3, Lhx6
- Cluster10= Glut5= Neurod2, Neurod6, Rbfox1
- Cluster11= Glut6= Neurod2, Neurod6, Rbfox1
- Cluster12= In4= Gad2, Nrxn3, Lhx6
- Cluster13= Astro= Aqp4
- Cluster14= NPC= Eomes
- Cluster15= OPC= Olig1, Olig2, Pdgfra
- Cluster16= Neuron= Reln, Slc17a6
- Cluster17= Endo= Col4a1, Col4a2
- Cluster18= Mg= C1qa, C1qb



```R
Idents(ATACMultiomewithST_SMK_V2QCv2.sct) <- "seurat_clusters"

# Cell type name in order of seurat_cluster
new.cluster.ids <- c(
  "Glut1",
  "Glut2",
  "Glut3",
  "Glut4",
  "NSC",
  "RG",
  "In1",
  "In2",
  "In3",
  "Glut5",
  "Glut6",
  "In4",
  "Astro",
  "NPC",
  "OPC",
  "Neuron",
  "Endo",
  "Mg"
)

names(new.cluster.ids) <- levels(ATACMultiomewithST_SMK_V2QCv2.sct)
ATACMultiomewithST_SMK_V2QCv2.sct <- RenameIdents(ATACMultiomewithST_SMK_V2QCv2.sct, new.cluster.ids)
ATACMultiomewithST_SMK_V2QCv2.sct$cluster.annot <- Idents(ATACMultiomewithST_SMK_V2QCv2.sct) # create a new slot in my seurat object


# UMAP with update cell type name
pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-label.pdf", width=15, height=6)
DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 6)
dev.off()
pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelNoSplit.pdf", width=7, height=6)
DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap",  label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 5)
dev.off()



# Gene expression marker in dotplot
all_markers <- c(
  "Aqp4",
  "Col4a1", "Col4a2",
  "Neurod2", "Neurod6", "Rbfox1",
  "Gad2", "Nrxn3", "Lhx6",
  "C1qa", "C1qb",
  "Reln", "Slc17a6",
  "Eomes",
  "Mki67", "Top2a", "Dlx6os1",
  "Olig1", "Olig2", "Pdgfra",
  "Fabp7", "Slc1a3", "Tnc"
)



levels(ATACMultiomewithST_SMK_V2QCv2.sct) <- c(
  "Astro",
  "Endo",
  "Glut1",
  "Glut2",
  "Glut3",
  "Glut4",
  "Glut5",
  "Glut6",
  "In1",
  "In2",
  "In3",
  "In4",
  "Mg",
  "Neuron",
  "NPC",
  "NSC",
  "OPC",
  "RG"
)


pdf("output/seurat/DotPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-all_markers.pdf", width=8, height=4.5)
DotPlot(ATACMultiomewithST_SMK_V2QCv2.sct, assay = "SCT", features = all_markers, cols = c("#1A00FF", "#FF0000")) + RotatedAxis() + scale_y_discrete(limits = rev)
dev.off()
```


### Cell type marker RNA assay

Identify marker genes in each cell types

```R
Idents(ATACMultiomewithST_SMK_V2QCv2.sct) <- "cluster.annot"

# Normalize and scale RNA assay prior DGEs
DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- "RNA"
ATACMultiomewithST_SMK_V2QCv2.sct <- NormalizeData(ATACMultiomewithST_SMK_V2QCv2.sct, normalization.method = "LogNormalize", scale.factor = 10000) 
all.genes <- rownames(ATACMultiomewithST_SMK_V2QCv2.sct)
ATACMultiomewithST_SMK_V2QCv2.sct <- ScaleData(ATACMultiomewithST_SMK_V2QCv2.sct, features = all.genes) 

all_markers <- FindAllMarkers(ATACMultiomewithST_SMK_V2QCv2.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_markers, file = "output/seurat/FindAllMarkers-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-all_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)
```





### DGEs analysis - RNA assay


Identify DEGs (Normoxia vs Hypoxia) in each cell types



```R
DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- "RNA"

# Normalize and scale RNA assay prior DGEs
ATACMultiomewithST_SMK_V2QCv2.sct <- NormalizeData(ATACMultiomewithST_SMK_V2QCv2.sct, normalization.method = "LogNormalize", scale.factor = 10000) # 
all.genes <- rownames(ATACMultiomewithST_SMK_V2QCv2.sct)
ATACMultiomewithST_SMK_V2QCv2.sct <- ScaleData(ATACMultiomewithST_SMK_V2QCv2.sct, features = all.genes) 


ATACMultiomewithST_SMK_V2QCv2.sct$celltype.stim <- paste(ATACMultiomewithST_SMK_V2QCv2.sct$cluster.annot, ATACMultiomewithST_SMK_V2QCv2.sct$condition,
    sep = "-")
Idents(ATACMultiomewithST_SMK_V2QCv2.sct) <- "celltype.stim"

cell_types <- c(  "Astro",
  "Endo",
  "Glut1",
  "Glut2",
  "Glut3",
  "Glut4",
  "Glut5",
  "Glut6",
  "In1",
  "In2",
  "In3",
  "In4",
  "Mg",
  "Neuron",
  "NPC",
  "NSC",
  "OPC",
  "RG")

for (cell_type in cell_types) {
  response_name <- paste(cell_type, "Hypoxia-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.response", sep = ".")
  ident_1 <- paste0(cell_type, "-Hypoxia")
  ident_2 <- paste0(cell_type, "-Normoxia")

  file_name <- paste0("output/seurat/", cell_type,
                      "-Hypoxia-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.txt")

  message("Running DE for cluster ", cell_type, " ...")

  tryCatch(
    {
      response <- FindMarkers(
        ATACMultiomewithST_SMK_V2QCv2.sct,
        assay = "RNA",
        ident.1 = ident_1,
        ident.2 = ident_2,
        verbose = FALSE
      )
      
      print(head(response, 15))
      write.table(response, file = file_name, sep = "\t", quote = FALSE, row.names = TRUE)
    },
    error = function(e) {
      message("⚠️  Skipping cluster ", cell_type, ": ", conditionMessage(e))
    }
  )
}


###########################################################################
# saveRDS(ATACMultiomewithST_SMK_V2QCv2.sct, file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1.rds") 
ATACMultiomewithST_SMK_V2QCv2.sct <- readRDS(file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1.rds")
###########################################################################
```

DEGs save as output/seurat/[CLUSTERNAME]-Hypoxia-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.txt


Seurat object at this stage can be accessed at: `output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1.rds`





### WNN integration

Generation of WNN UMAP (UMAP taking into account RNA and ATAC assays)


```R
# pre-processing and dimensional reduction
DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- "peaks"
ATACMultiomewithST_SMK_V2QCv2.sct <- RunTFIDF(ATACMultiomewithST_SMK_V2QCv2.sct)
ATACMultiomewithST_SMK_V2QCv2.sct <- FindTopFeatures(ATACMultiomewithST_SMK_V2QCv2.sct, min.cutoff = 'q0')
ATACMultiomewithST_SMK_V2QCv2.sct <- RunSVD(ATACMultiomewithST_SMK_V2QCv2.sct)
ATACMultiomewithST_SMK_V2QCv2.sct <- RunUMAP(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = 'lsi', dims = 2:40, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
#--> We exclude the first dimension as this is typically correlated with sequencing depth


ATACMultiomewithST_SMK_V2QCv2.sct <- FindMultiModalNeighbors(ATACMultiomewithST_SMK_V2QCv2.sct, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:40))
ATACMultiomewithST_SMK_V2QCv2.sct <- RunUMAP(ATACMultiomewithST_SMK_V2QCv2.sct, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
ATACMultiomewithST_SMK_V2QCv2.sct <- FindClusters(ATACMultiomewithST_SMK_V2QCv2.sct, graph.name = "wsnn", algorithm = 3, verbose = TRUE)

# QC plot to check which dims is corr with seq depth
DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- "peaks"
pdf("output/seurat/DepthCor-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04_ATACdim240-ATACdim240.pdf", width=10, height=5)
DepthCor(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "lsi", n = 40)
dev.off()


# plot gene expression, ATAC-seq, or WNN analysis
p1 <- DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap", group.by = "cluster.annot", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap.atac", group.by = "cluster.annot", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "wnn.umap", group.by = "cluster.annot", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("WNN")
pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04_ATACdim240-RNAATACWNN.pdf", width=14, height=5)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()


# plot separating condition
ATACMultiomewithST_SMK_V2QCv2.sct$condition <- factor(ATACMultiomewithST_SMK_V2QCv2.sct$condition, levels = c("Normoxia", "Hypoxia")) # Reorder untreated 1st
pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04_ATACdim240-WNN.pdf", width=5, height=5)
p = DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "wnn.umap", group.by = "cluster.annot", label = FALSE, label.size = 4, repel = TRUE) + ggtitle("Cell type") + NoLegend()
LabelClusters(p, id = "cluster.annot", fontface = "bold", color = "black", size = 3)
dev.off()

pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04_ATACdim240-RNAcondition.pdf", width=5, height=5)
DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap", group.by = "condition", label = FALSE, cols = c("blue", "red")) + ggtitle("RNA")  + NoLegend()
dev.off()
pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04_ATACdim240-ATACcondition.pdf", width=5, height=5)
DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap.atac", group.by = "condition", label = FALSE, cols = c("blue", "red")) + ggtitle("ATAC")  + NoLegend()
dev.off()
pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04_ATACdim240-WNNcondition.pdf", width=5, height=5)
DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "wnn.umap", group.by = "condition", label = FALSE, cols = c("blue", "red")) + ggtitle("WNN")  + NoLegend()
dev.off()
```



### Calculate promoter activity (GeneActivity)


Count ATAC reads within gene promoters (TSS up 2kb)

```R
GeneActivity = GeneActivity(
  ATACMultiomewithST_SMK_V2QCv2.sct,
  assay = "peaks",
  features = NULL,
  extend.upstream = 2000,
  extend.downstream = 0,
  biotypes = NULL,
  max.width = NULL,
  process_n = 2000,
  gene.id = FALSE,
  verbose = TRUE
)

# Add gene activity as a new assay
ATACMultiomewithST_SMK_V2QCv2.sct[["GeneActivity"]] <- CreateAssayObject(counts = GeneActivity)


###########################################################################
# saveRDS(ATACMultiomewithST_SMK_V2QCv2.sct, file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1GeneActivity.rds") 
ATACMultiomewithST_SMK_V2QCv2.sct <- readRDS(file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1GeneActivity.rds")
###########################################################################


# Marker genes GeneActivity
all_markers <- c(
  "Aqp4",
  "Col4a1", "Col4a2",
  "Neurod2", "Neurod6", "Rbfox1",
  "Gad2", "Nrxn3", "Lhx6",
  "C1qa", "C1qb",
  "Reln", "Slc17a6",
  "Eomes",
  "Mki67", "Top2a", "Dlx6os1",
  "Olig1", "Olig2", "Pdgfra",
  "Fabp7", "Slc1a3", "Tnc"
)


markers_ATAC_Cristancho <- c(
  "Neurod6", "Tle4", "Zfpm2", "Satb2", "Dok5", "Sema3c", "Sema3a",  "Reln", "Ldb2", "Gad2", "Lhx6", "Sst", "Isl1", "Drd2", "Eomes", "Mki67", "Top2a", "Tnc", "Olig2", "Slc1a3", "Fabp7", "Col4a1", "Cldn5", "C1qb", "C1qa"
)

Idents(ATACMultiomewithST_SMK_V2QCv2.sct) <- "cluster.annot"

levels(ATACMultiomewithST_SMK_V2QCv2.sct) <- c(
  "Astro",
  "Endo",
  "Glut1",
  "Glut2",
  "Glut3",
  "Glut4",
  "Glut5",
  "Glut6",
  "In1",
  "In2",
  "In3",
  "In4",
  "Mg",
  "Neuron",
  "NPC",
  "NSC",
  "OPC",
  "RG"
)


pdf("output/seurat/DotPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-all_markers-GeneActivity.pdf", width=8, height=4.5)
DotPlot(ATACMultiomewithST_SMK_V2QCv2.sct, assay = "GeneActivity", features = all_markers, cols = c("#0A7B83", "#E78F0B")) + RotatedAxis() + scale_y_discrete(limits = rev)
dev.off()
```

Seurat object at this stage can be accessed at: `output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1GeneActivity.rds`




### DARs analysis with Signac - ATAC assay


```R
# Identify diff access regions (DAR)
DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- 'peaks'



Idents(ATACMultiomewithST_SMK_V2QCv2.sct) <- "celltype.stim"

### Here is to automatize the process clsuter per cluster:
#### Define the cluster pairs for comparison
cluster_pairs <- list(
  c("Astro-Normoxia",   "Astro-Hypoxia",   "Astro"),
  c("Endo-Normoxia",   "Endo-Hypoxia",    "Endo"),
  c("Glut1-Normoxia",  "Glut1-Hypoxia",   "Glut1"),
  c("Glut2-Normoxia",  "Glut2-Hypoxia",   "Glut2"),
  c("Glut3-Normoxia",  "Glut3-Hypoxia",   "Glut3"),
  c("Glut4-Normoxia",  "Glut4-Hypoxia",   "Glut4"),
  c("Glut5-Normoxia",  "Glut5-Hypoxia",   "Glut5"),
  c("Glut6-Normoxia",  "Glut6-Hypoxia",   "Glut6"),
  c("In1-Normoxia",    "In1-Hypoxia",     "In1"),
  c("In2-Normoxia",    "In2-Hypoxia",     "In2"),
  c("In3-Normoxia",    "In3-Hypoxia",     "In3"),
  c("In4-Normoxia",    "In4-Hypoxia",     "In4"),
  c("Mg-Normoxia",     "Mg-Hypoxia",      "Mg"),
  c("Neuron-Normoxia", "Neuron-Hypoxia",  "Neuron"),
  c("NPC-Normoxia",    "NPC-Hypoxia",     "NPC"),
  c("NSC-Normoxia",    "NSC-Hypoxia",     "NSC"),
  c("OPC-Normoxia",    "OPC-Hypoxia",     "OPC"),
  c("RG-Normoxia",     "RG-Hypoxia",      "RG")
)


run_find_markers <- function(pair) {
  ident_1 <- pair[2] 
  ident_2 <- pair[1]
  cluster_name <- pair[3]
  
  # Run FindMarkers for each pair of clusters
  markers <- FindMarkers(
    object = ATACMultiomewithST_SMK_V2QCv2.sct,
    ident.1 = ident_1,
    ident.2 = ident_2,    verbose = TRUE,
    test.use = "LR",
    min.pct = 0.1, # 
    assay = "peaks",
  latent.vars = c("nCount_peaks" )) %>%
  rownames_to_column(var = "query_region") %>%  # Add gene names as a column
  add_column(cluster = cluster_name)    # Add cluster name as a new column
  
  return(markers)
}

# Combine all FindMarkers results into one tibble
all_markers_tibble <- map_dfr(cluster_pairs, run_find_markers) %>%
  as_tibble()
# save output
write.table(all_markers_tibble, file = "output/seurat/DAR_peaks-LR_latentnCount_peaks-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE) # new version with LR test anmd latent vars "nCount_peaks", "TSS.enrichment" 
```

### Find closest gene to DARs



```R
# Find closest gene to DAR
# Download genome file mm39 and gene version v106 - work with BD Rhapsody
library("AnnotationHub")
library("ensembldb")

ah <- AnnotationHub()

# Query EnsDb mouse + Ensembl release 106 (GRCm39/mm39-era)
hits <- query(ah, c("EnsDb", "Mus musculus", "Ensembl", "106"))
hits
EnsDb.Mmusculus.v106 <- ah[["AH100674"]]  
EnsDb.Mmusculus.v106
library("BSgenome.Mmusculus.UCSC.mm39") # BiocManager::install("BSgenome.Mmusculus.UCSC.mm39")

# Add genomic information to Seurat object
# Convert rownames of ATAC counts to GRanges
grange.counts <- StringToGRanges(rownames(ATACMultiomewithST_SMK_V2QCv2.sct[["peaks"]]), sep = c(":", "-"))
# Filter for standard chromosomes
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- ATACMultiomewithST_SMK_V2QCv2.sct[as.vector(grange.use), ]
# Get annotations for the mouse genome
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v106) # EnsDb.Mmusculus.v79 if using mm10
# Adjust chromosome naming style
seqlevelsStyle(annotations) <- 'UCSC'

# Identify closest gene to DAR
DAR_closestGene <- ClosestFeature(ATACMultiomewithST_SMK_V2QCv2.sct, regions = all_markers_tibble$query_region, annotation = annotations)
DAR_genes = all_markers_tibble %>% 
  left_join(DAR_closestGene) %>%
  as_tibble() %>%
  unique()


write.table(DAR_genes, file = "output/seurat/DAR_genes-LR_latentnCount_peaks-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
```



### Link peak to genes



```R
# Link peaks to genes - Find peaks that are correlated with the expression of nearby genes
ATACMultiomewithST_SMK_V2QCv2.sct <- RegionStats(
  object = ATACMultiomewithST_SMK_V2QCv2.sct,
  genome = BSgenome.Mmusculus.UCSC.mm39,
  assay = "peaks",
  verbose = TRUE
)

# Run  LinkPeaks
ATACMultiomewithST_SMK_V2QCv2.sct = LinkPeaks(
  ATACMultiomewithST_SMK_V2QCv2.sct,
  peak.assay = "peaks",
  expression.assay = "RNA",
  peak.slot = "counts",
  expression.slot = "data",
  method = "pearson",
  gene.coords = NULL,
  distance = 5e+05,
  min.distance = NULL,
  min.cells = 10,
  genes.use = NULL,
  n_sample = 200,
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
  gene.id = FALSE,
  verbose = TRUE
)


###########################################################################
# saveRDS(ATACMultiomewithST_SMK_V2QCv2.sct, file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1GeneActivityLinkPeaks.rds") 
ATACMultiomewithST_SMK_V2QCv2.sct <- readRDS(file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1GeneActivityLinkPeaks.rds")
###########################################################################



Links = as_tibble(Links(ATACMultiomewithST_SMK_V2QCv2.sct))

### Adjust the pvalue and select positive corr
Links$adjusted_pvalue <- p.adjust(Links$pvalue, method = "BH")
write.table(Links, file = c("output/seurat/LinkPeaks-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.txt"),sep="\t", quote=FALSE, row.names=FALSE)

# Isolate signif genes
Links_signif = Links %>%
  dplyr::filter(adjusted_pvalue < 0.05, score >0) %>%  
  dplyr::select(gene) %>% # 1,942 Link Signif
  unique()

Links %>%
  dplyr::filter(adjusted_pvalue < 0.05, score >0) %>%  
  dplyr::select(peak) %>% # 3,740 Link peak Signif
  unique()

# Normalize and scale GeneActivity
DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- "GeneActivity"

ATACMultiomewithST_SMK_V2QCv2.sct <- NormalizeData(ATACMultiomewithST_SMK_V2QCv2.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(ATACMultiomewithST_SMK_V2QCv2.sct)
ATACMultiomewithST_SMK_V2QCv2.sct <- ScaleData(ATACMultiomewithST_SMK_V2QCv2.sct, features = all.genes) # zero-centres and scales it



# Find all markers 
Idents(ATACMultiomewithST_SMK_V2QCv2.sct) <- "cluster.annot"
DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- "RNA"

Links_markers <- FindAllMarkers(ATACMultiomewithST_SMK_V2QCv2.sct, features = Links_signif$gene, assay = "RNA", only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
# Identify in which cluster the Links_markers gene is highly express
Links_markers_pval= as_tibble(Links_markers) %>%
  group_by(gene) %>%
  dplyr::filter(p_val == min(p_val)) %>%
  dplyr::select(gene, cluster)



# plot heatmap
pdf("output/seurat/DoHeatmap-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1GeneActivityLinkPeaks-LinksPadj05Score0-SCTscaledata.pdf", width=8, height=4)
DoHeatmap(ATACMultiomewithST_SMK_V2QCv2.sct, assay = "SCT", slot= "scale.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5, draw.lines = FALSE, label = FALSE)
dev.off()
pdf("output/seurat/DoHeatmap-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1GeneActivityLinkPeaks-LinksPadj05Score0-GeneActivityscaledata1.pdf", width=8, height=4)
DoHeatmap(ATACMultiomewithST_SMK_V2QCv2.sct, assay = "GeneActivity", slot= "scale.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5, draw.lines = FALSE, label = FALSE, disp.max = 1.25, disp.min = -2)
dev.off()
```


### TSS plot


```R
# TSSplot 
DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- "peaks"
ATACMultiomewithST_SMK_V2QCv2.sct = TSSEnrichment(ATACMultiomewithST_SMK_V2QCv2.sct, assay = "peaks", fast = FALSE, verbose = TRUE)

# all genes - all cluster
pdf("output/seurat/TSSPlot-ATACMultiomewithST_SMK_V2QCv2-all.pdf", width = 6, height = 5)
TSSPlot(
  ATACMultiomewithST_SMK_V2QCv2.sct,
  assay = "peaks",
  group.by = "condition"
) +
  facet_null() +
  scale_color_manual(values = c(
    "Normoxia" = "#3700ff",      # cyan
    "Hypoxia" = "#f10e0e"   # red
  )) +
  theme_bw()
dev.off()


# One plot per cell type ### 
clusters_all <- unique(ATACMultiomewithST_SMK_V2QCv2.sct$cluster.annot)

for (clust in clusters_all) {
  obj_sub <- subset(
    ATACMultiomewithST_SMK_V2QCv2.sct,
    subset = cluster.annot == clust
  )
  p <- TSSPlot(
    obj_sub,
    assay = "peaks",
    group.by = "condition"
  ) +
    facet_null() +
    scale_color_manual(values = c(
    "Normoxia" = "#3700ff",      # cyan
    "Hypoxia" = "#f10e0e"   # red
    )) +
    theme_bw() +
    ggtitle(clust) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  file_name <- paste0(
    "output/seurat/TSSPlot_",
    gsub("[ /]", "_", clust),
    "_ATACMultiomewithST_SMK_V2QCv2.pdf"
  )
  pdf(file_name, width = 6, height = 5)
  print(p)
  dev.off()
}
```






### DARs analysis with archR - ATAC assay


```R
library("ArchR")
library("Seurat")
library("Signac")
library("SingleCellExperiment")
library("BSgenome.Mmusculus.UCSC.mm10")

set.seed(42)


# load
ATACMultiomewithST_SMK_V2QCv2.sct <- readRDS(file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1GeneActivityLinkPeaks.rds")

Idents(ATACMultiomewithST_SMK_V2QCv2.sct) <- "Sample_Name"
levels = levels(ATACMultiomewithST_SMK_V2QCv2.sct@active.ident)

addArchRGenome("mm10")

# Load fragment
files <- c(
  SMK = "output/seurat/ATACMultiomewithST-CC-additionalSMK_ATAC_Fragments.bed.gz"
)

ArrowFiles <- createArrowFiles(
  inputFiles = files,
  sampleNames = names(files),
  promoterRegion = c(2000, 100),
  minTSS = 1,
  minFrags = 100,
  maxFrags = 1e5,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  excludeChr = c("chrM"),
  force = TRUE
)

# Create ArchR project
SMK_proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "output/archr",
  copyArrows = TRUE
)

SMK_proj



##########################################################
# Add cell type annotation from seurat to archR object ###
##########################################################
ATACMultiomewithST_SMK_V2QCv2.sct


## Seurat object
seu <- ATACMultiomewithST_SMK_V2QCv2.sct

## Metadata columns to transfer
cols_to_transfer <- c(
  "cluster.annot",
  "Sample_Name",
  "sex",
  "condition",
  "seurat_clusters"
)

## Extract Seurat metadata
meta <- seu@meta.data[, cols_to_transfer, drop = FALSE]

## Convert factors to characters
meta[] <- lapply(meta, function(x) {
  if (is.factor(x)) as.character(x) else x
})
meta$seurat_clusters <- as.character(meta$seurat_clusters)


## Make Seurat cell names match ArchR cell names
meta$ArchR_cell <- paste0("SMK#", rownames(meta))

archr_cells <- getCellNames(SMK_proj)

cat("Overlap after adding SMK# prefix:\n")
cat(sum(meta$ArchR_cell %in% archr_cells), "Seurat cells found in ArchR\n")
cat(length(archr_cells), "total ArchR cells\n")
cat(nrow(meta), "total Seurat cells\n")


## Keep only cells present in ArchR
meta_keep <- meta[meta$ArchR_cell %in% archr_cells, ]

cat("Final cells transferred:", nrow(meta_keep), "\n")


rownames(meta_keep) <- meta_keep$ArchR_cell

for (col in cols_to_transfer) {
  
  SMK_proj <- addCellColData(
    ArchRProj = SMK_proj,
    data = meta_keep[[col]],
    cells = rownames(meta_keep),
    name = col,
    force = TRUE
  )
  
  cat("Added:", col, "\n")
}


table(SMK_proj$cluster.annot)
table(SMK_proj$Sample_Name)
table(SMK_proj$sex)
table(SMK_proj$condition)
table(SMK_proj$seurat_clusters)


#################################################
# Only keep archr cells overlapping with seurat
#################################################
archr_cells <- getCellNames(SMK_proj)

## Seurat cells are like "45285702"
## ArchR cells are like "SMK#45285702"
seurat_cells_archr_style <- paste0("SMK#", colnames(seu))

overlap_cells <- intersect(archr_cells, seurat_cells_archr_style)

cat("Overlap cells:", length(overlap_cells), "\n")
cat("ArchR cells:", length(archr_cells), "\n")
cat("Seurat cells:", ncol(seu), "\n")


SMK_proj_multiome <- subsetArchRProject(
  ArchRProj = SMK_proj,
  cells = overlap_cells,
  outputDirectory = "output/archr",
  force = TRUE
)

SMK_proj_multiome = loadArchRProject(path = "output/archr")


#################################################
# Get RNA raw counts from Seurat
#################################################

DefaultAssay(seu_sub) <- "RNA"

rna_counts <- GetAssayData(
  seu_sub,
  assay = "RNA",
  layer = "counts"
)

## Rename columns to ArchR style
colnames(rna_counts) <- paste0("SMK#", colnames(rna_counts))

## Keep/reorder cells exactly like ArchR
rna_counts <- rna_counts[, getCellNames(SMK_proj_multiome)]

# get gene annotation from archr
geneAnno <- getGeneAnnotation(SMK_proj_multiome)
genes_gr <- geneAnno$genes

mcols(genes_gr) |> head()



gene_symbols <- mcols(genes_gr)$symbol

## Remove duplicated gene symbols in ArchR annotation
keep_gene_anno <- !is.na(gene_symbols) & !duplicated(gene_symbols)

genes_gr <- genes_gr[keep_gene_anno]
gene_symbols <- gene_symbols[keep_gene_anno]

## Match Seurat RNA genes to ArchR gene annotation
common_genes <- intersect(rownames(rna_counts), gene_symbols)

cat("Common genes:", length(common_genes), "\n")
cat("Seurat RNA genes:", nrow(rna_counts), "\n")
cat("ArchR annotated genes:", length(gene_symbols), "\n")



rna_counts2 <- rna_counts[common_genes, ]

genes_gr2 <- genes_gr[match(common_genes, gene_symbols)]

## Make sure rownames and rowRanges names match
names(genes_gr2) <- common_genes

seRNA <- SingleCellExperiment(
  assays = list(counts = rna_counts2),
  rowRanges = genes_gr2
)

rownames(seRNA) <- common_genes
colnames(seRNA) <- colnames(rna_counts2)


SMK_proj_multiome <- addGeneExpressionMatrix(
  input = SMK_proj_multiome,
  seRNA = seRNA,
  strictMatch = TRUE,
  force = TRUE
)

getAvailableMatrices(SMK_proj_multiome)

#--> GeneExpressionMatrix is here

#################################################
# Add PeakMatrix in ArchR #####
#################################################
## run LSI
SMK_proj_multiome <- addIterativeLSI(
  ArchRProj = SMK_proj_multiome,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 4,
  clusterParams = list(
    resolution = c(0.8),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 2:40,  # 2:40 as in the seurat
  force = TRUE
)


SMK_proj_multiome$cluster_condition <- paste0(
  SMK_proj_multiome$cluster.annot,
  "_",
  SMK_proj_multiome$condition
)


## pseudo-bulk coverages
SMK_proj_multiome <- addGroupCoverages(
  ArchRProj = SMK_proj_multiome,
  groupBy = "cluster_condition",
  minCells = 30,
  maxCells = 100000,
  maxReplicates = 8,
  force = TRUE
)




# Call macs2 from another conda env
pathToMacs2 <- "/home/roulet/anaconda3/envs/macs2/bin/macs2"

file.exists(pathToMacs2)
system2(pathToMacs2, "--version", stdout = TRUE, stderr = TRUE)



pdf("output/archr/addReproduciblePeakSet.pdf", width=10, height=10)
SMK_proj_multiome <- addReproduciblePeakSet(
  ArchRProj = SMK_proj_multiome,
  groupBy = "cluster_condition",
  pathToMacs2 = pathToMacs2,
  force = TRUE,
  maxPeaks = 10*10^7
)
dev.off()
## add the peak matrix
SMK_proj_multiome <- addPeakMatrix(
  ArchRProj = SMK_proj_multiome,
  force = TRUE
)

getAvailableMatrices(SMK_proj_multiome)





##################
# Save ArchR project
SMK_proj_multiome <- saveArchRProject(
  ArchRProj = SMK_proj_multiome,
  outputDirectory = "output/archr",
  overwrite = TRUE,
  load = TRUE,
  dropCells = FALSE,
  logFile = createLogFile("output/archr/saveArchRProject"),
  threads = getArchRThreads()
)

# LOAD ArchR project
SMK_proj_multiome = loadArchRProject(path = "output/archr")
####################




#################################################
# Add co-accessibility and  peak to gene linkage #####
#################################################


SMK_proj_multiome = addPeak2GeneLinks(ArchRProj = SMK_proj_multiome, useMatrix = "GeneExpressionMatrix", reducedDims= "IterativeLSI", dimsToUse= 1:30)


p2g = getPeak2GeneLinks(ArchRProj = SMK_proj_multiome, corCutOff = 0.45, resolution = 1, returnLoops= TRUE)


pdf("output/archr/plotPeak2GeneHeatmap-corCutOff045res1groupBycluster_condition.pdf", width=5, height=5)
plotPeak2GeneHeatmap(ArchRProj = SMK_proj_multiome, groupBy= "cluster_condition")
dev.off()


h_sub= plotPeak2GeneHeatmap(ArchRProj = SMK_proj_multiome, groupBy= "cluster_condition", k=6)
h_sub1= plotPeak2GeneHeatmap(ArchRProj = SMK_proj_multiome, groupBy= "cluster.annot", k=6)
h_sub2= plotPeak2GeneHeatmap(ArchRProj = SMK_proj_multiome, groupBy= "condition", k=6)
h_sub_mat = plotPeak2GeneHeatmap(ArchRProj = SMK_proj_multiome, groupBy= "cluster_condition", k=6, returnMatrices= TRUE)

ATAC_p2g = h_sub@ht_list[["ATAC Z-Scores\n2383 P2GLinks"]]
original_order = ATAC_p2g@column_order


pdf("output/archr/ATAC_p2g-corCutOff045res1groupBycluster_condition_k6.pdf", width=5, height=5)
ATAC_p2g
dev.off()


Knames = h_sub_mat@listData[["ATAC"]]@listData[["colData"]]@rownames
cellnames = h_sub_mat@listData[["ATAC"]]@listData[["colData"]]@listData[["groupBy"]]


df <- data.frame(
  Knames = Knames,
  cellnames = cellnames,
  stringsAsFactors = FALSE
)

df[c("cluster.annot", "condition")] <- str_split_fixed(df$cellnames, "_", 2)

## Check
head(df)
table(df$condition)
table(df$cluster.annot)


## Make condition order explicit
df$condition <- factor(df$condition, levels = c("Normoxia", "Hypoxia"))

## Sort by condition, then cell type
df_sorted <- df[order(df$condition, df$cluster.annot), ]

## This is the order to apply to the heatmap columns
reorder <- as.numeric(rownames(df_sorted))

ATAC_p2g@column_order <- reorder


df_ann <- df[, c("cluster.annot", "condition")]

df_ann$condition <- as.character(df_ann$condition)
df_ann$celltype <- as.character(df_ann$cluster.annot)
df_ann$cluster.annot <- NULL





pdf("output/archr/ATAC_p2g-corCutOff045res1groupBycluster_condition_k6-column_order.pdf", width=5, height=5)
ATAC_p2g
dev.off()


colors <- list(
    condition = c('Normoxia' = "tomato", 'Hypoxia' = 'deepskyblue3'))

colors <- list(
  condition = c(
    "Normoxia" = "tomato",
    "Hypoxia"  = "deepskyblue3"
  ),
  
  celltype = c(
    "Astro"  = "blue",
    "Endo"   = "goldenrod3",
    
    # Glut clusters = green gradient
    "Glut1"  = "#00441B",
    "Glut2"  = "#006D2C",
    "Glut3"  = "#238B45",
    "Glut4"  = "#41AB5D",
    "Glut5"  = "#74C476",
    
    # Inhibitory clusters = red gradient
    "In1"    = "#FCAE91",
    "In2"    = "#FB6A4A",
    "In3"    = "#CB181D",
    
    "Mg"     = "brown4",
    "Neuron" = "royalblue4",
    "NSC"    = "lightblue2",
    "OPC"    = "purple",
    "RG"     = "yellow"
  )
)

colAnn <- HeatmapAnnotation(
  df = df_ann,
  which = "column",
  na_col = "white",
  col = colors,
  simple_anno_size = unit(0.6, "cm"),
  gap = unit(1, "mm"),
  annotation_legend_param = list(
    condition = list(
      nrow = 2,
      title = "condition",
      title_position = "topcenter",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 10)
    ),
    celltype = list(
      nrow = 16,
      title = "celltype",
      title_position = "topcenter",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 10)
    )
  )
)


pdf("output/archr/ATAC_p2g-corCutOff045res1groupBycluster_condition_k6-column_order.pdf", width=5, height=5)
ATAC_p2g@top_annotation<-colAnn
draw(ATAC_p2g, heatmap_legend_side = "bottom")
dev.off()


#################
# heatmap FDR ###
#################
mat_atac <- ATAC_p2g@matrix

## Get row clusters from the heatmap k = 6
pdf("output/archr/ATAC_p2g-corCutOff045res1groupBycluster_condition_k6-ht_tmp.pdf", width=5, height=5)
ht_tmp <- draw(ATAC_p2g)
dev.off()

row_groups <- row_order(ht_tmp)

## If ComplexHeatmap returns nested list, simplify
if (length(row_groups) == 1 && is.list(row_groups[[1]])) {
  row_groups <- row_groups[[1]]
}

names(row_groups) <- paste0("cluster", seq_along(row_groups))

cluster_mats <- lapply(row_groups, function(ii) {
  mat_atac[ii, , drop = FALSE]
})


Knames <- h_sub_mat@listData[["ATAC"]]@listData[["colData"]]@rownames
cellnames <- h_sub_mat@listData[["ATAC"]]@listData[["colData"]]@listData[["groupBy"]]

df2 <- data.frame(
  Knames = Knames,
  cellnames = cellnames,
  stringsAsFactors = FALSE
)

df2[c("celltype", "condition")] <- stringr::str_split_fixed(df2$cellnames, "_", 2)

df2$condition <- factor(df2$condition, levels = c("Normoxia", "Hypoxia"))

celltype_order <- c(
  "Astro", "Endo",
  paste0("Glut", 1:5),
  paste0("In", 1:3),
  "Mg", "Neuron", "NSC", "OPC", "RG"
)

celltype_order <- celltype_order[celltype_order %in% unique(df2$celltype)]

df2$celltype <- factor(df2$celltype, levels = celltype_order)

## Check matrix columns match metadata rows
stopifnot(ncol(mat_atac) == nrow(df2))



for (cl in names(cluster_mats)) {
  df2[[cl]] <- apply(cluster_mats[[cl]], 2, median, na.rm = TRUE)
}

################################################
## Permutation test: Normoxia vs Hypoxia per cell type and cluster
################################################
library("purrr")
library("perm") # install.packages("perm")
library("tidyr")


s <- split(df2, df2$celltype)

df_result <- data.frame(
  celltype = character(),
  cluster = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (ct in names(s)) {
  
  for (cl in names(cluster_mats)) {
    
    tmp <- s[[ct]]
    
    ## skip if cluster column missing
    if (!cl %in% colnames(tmp)) next
    
    ## keep only non-NA values
    tmp <- tmp %>%
      dplyr::filter(
        !is.na(.data[[cl]]),
        !is.na(condition)
      )
    
    ## need both conditions
    if (length(unique(tmp$condition)) < 2) next
    
    ## need at least 2 total values
    if (nrow(tmp) < 2) next
    
    message("Testing ", ct, " / ", cl)
    
    result <- permKS(tmp[[cl]], tmp$condition)
    
    df_result1 <- data.frame(
      celltype = ct,
      cluster = cl,
      p_value = result[["p.value"]]
    )
    
    df_result <- rbind(df_result, df_result1)
  }
}


df_result$fdr <- p.adjust(df_result$p_value, method = "BH")




  
write.csv(
  df_result,
  "output/archr/ATAC_p2g_density_permKS_FDR_results.csv",
  row.names = FALSE
)


################################################
## Make FDR matrix for heatmap
################################################

matrix_fdr <- df_result %>%
  mutate(
    cluster = gsub("cluster", "", cluster),
    celltype = factor(celltype, levels = celltype_order)
  ) %>%
  dplyr::select(celltype, cluster, fdr) %>%
  pivot_wider(names_from = cluster, values_from = fdr) %>%
  arrange(celltype) %>%
  as.data.frame()

rownames(matrix_fdr) <- matrix_fdr$celltype
matrix_fdr$celltype <- NULL
matrix_fdr <- as.matrix(matrix_fdr)



col_fdr <- colorRamp2(
  c(0, 0.01, 0.05, 0.1, 0.5, 1),
  c("#034E7B", "#2987BB", "#8CB3D4", "#D7D6E9", "#FFF7FB", "#fafcfa")
)



pdf("output/archr/ATAC_p2g_density_permKS_FDR_heatmap.pdf", width = 4, height = 8)
Heatmap(
  matrix_fdr,
  name = "FDR",
  col = col_fdr,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_side = "top",
  column_names_rot = 0,
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title = "Cluster",
  column_title_gp = gpar(fontface = "bold"),
  na_col = "grey90"
)

dev.off()







##################
# Save ArchR project
SMK_proj_multiome <- saveArchRProject(
  ArchRProj = SMK_proj_multiome,
  outputDirectory = "output/archr",
  overwrite = TRUE,
  load = TRUE,
  dropCells = FALSE,
  logFile = createLogFile("output/archr/saveArchRProject"),
  threads = getArchRThreads()
)

# LOAD ArchR project
SMK_proj_multiome = loadArchRProject(path = "output/archr")
####################






###############################################
## Distribution of peaks and counts ###########
###############################################

peak_set<-getPeakSet(SMK_proj_multiome)

names<-as.data.frame(peak_set@ranges@NAMES)
names<-as.data.frame(str_split_fixed(names$`peak_set@ranges@NAMES`, "_", 2))

celltype<-names$V1
peak_set@elementMetadata@listData[["Celltype"]]<-celltype
condition<-names$V2

peak_set@elementMetadata@listData[["Condition"]]<-condition
peak_set



df<-data.frame(names = peak_set@ranges@NAMES, celltype = celltype, condition = condition, location = peak_set@elementMetadata@listData[["peakType"]])
table(df$condition, df$location)


pdf("output/archr/ggplot-PeakDistribution.pdf", width = 4, height = 2)
distribution<-ggplot(df, aes(x=condition, fill = location)) +
  geom_bar(position="fill") + theme_classic()+coord_flip()+ggtitle("Distribution of Peaks")+scale_y_continuous(expand = c(0,0))
distribution 
dev.off()

###############################################
## DA testing ###########
###############################################

clusters <- sort(unique(SMK_proj_multiome$cluster.annot))

DA_peak_list <- list()

for (cl in clusters) {
  
  use_group <- paste0(cl, "_Hypoxia") # Positive FC = more in hypoxia
  bgd_group <- paste0(cl, "_Normoxia")
  
  available_groups <- unique(SMK_proj_multiome$cluster_condition)
  
  if (!all(c(use_group, bgd_group) %in% available_groups)) {
    message("Skipping ", cl, ": missing Hypoxia or Normoxia")
    next
  }
  
  n_use <- sum(SMK_proj_multiome$cluster_condition == use_group)
  n_bgd <- sum(SMK_proj_multiome$cluster_condition == bgd_group)
  
  if (n_use < 10 || n_bgd < 10) {
    message("Skipping ", cl, ": too few cells. Hypoxia=", n_use, ", Normoxia=", n_bgd)
    next
  }
  
  message("Running DA for ", cl, ": ", use_group, " vs ", bgd_group)
  
  markerTest <- getMarkerFeatures(
    ArchRProj = SMK_proj_multiome,
    useMatrix = "PeakMatrix",
    groupBy = "cluster_condition",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = use_group,
    bgdGroups = bgd_group
  )
  
  marker_df <- getMarkers(
    markerTest,
    cutOff = "FDR <= 1"
  )[[use_group]]
  
  marker_df <- as.data.frame(marker_df)
  marker_df$cluster.annot <- cl
  marker_df$comparison <- paste0("Hypoxia_vs_Normoxia_", cl)
  
  DA_peak_list[[cl]] <- marker_df
}


DA_peaks_all <- do.call(rbind, DA_peak_list)
head(DA_peaks_all)

write.csv(
  DA_peaks_all,
  "output/archr/DA_peaks-Hypoxia_vs_Normoxia-by_cluster_all.csv",
  row.names = FALSE
)



## Filter DARs
DA_filt <- DA_peaks_all %>%
  as.data.frame() %>%
  dplyr::filter(
    FDR < 0.05,
    abs(Log2FC) >= 1
  ) %>%
  mutate(
    direction = ifelse(Log2FC > 0, "up", "down")
  )

## Since your comparison is Hypoxia_vs_Normoxia:
## up   = more accessible in Hypoxia
## down = more accessible in Normoxia
table(DA_filt$cluster.annot, DA_filt$direction)




## count DARs per cell type and direction
library("tidyr")
df_count <- DA_filt %>%
  count(cluster.annot, direction, name = "Count") %>%
  complete(
    cluster.annot,
    direction,
    fill = list(Count = 0)
  )

## plot
pdf("output/archr/df_count-DAR_cluster.pdf", width = 5, height = 4)
ggplot(df_count, aes(fill = direction, y = Count, x = cluster.annot)) +
  geom_bar(position = "stack", stat = "identity") +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limits = rev(levels(factor(df_count$cluster.annot)))) +
  scale_fill_manual(
    values = c(
      "up" = "darkseagreen",
      "down" = "goldenrod3"
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black", size = 12),
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_blank()
  ) +
  ylim(0,4000)
dev.off()
```










