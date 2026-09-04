##### SUPPLEMENTAL FIGURE 1 ######
### Kaplan-meier overall survival analysis 

## Curve Generation
#load packages
library(survival)
library(survminer)
library(ggplot2)

#load NSCLC patient datasets downloaded from cBioPortal after subsetting for KRAS and TP53 mutation status
data1 <- read.csv("K_PROCESSED.csv")
data2 <- read.csv("KP_HET_PROCESSED.csv")
data3 <- read.csv("KP_NULL_PROCESSED.csv")

# Combine datasets
combined_data <- rbind(data1, data2, data3)

# Cap overall survival time at 100 months
combined_data$Time_capped <- pmin(combined_data$Time, 100)
combined_data$Status_capped <- ifelse(
  combined_data$Time > 100,
  0,
  combined_data$Status
)

# Create Kaplan-Meier survival object and fit
surv_obj <- Surv(
  time = combined_data$Time_capped,
  event = combined_data$Status_capped
)

surv_fit <- survfit(
  surv_obj ~ Strata,
  data = combined_data
)

# Define Strata factor and patient group labels
combined_data$Strata <- factor(
  combined_data$Strata,
  levels = c("K", "KP_HET", "KP_NULL"),
  labels = c(
    "KRAS",
    "KRAS + single-allele TP53 loss (KP+/-)",
    "KRAS + complete TP53 loss (KP)"
  )
)

##### 1A #####
# Kaplan-Meier survival plot
res <- ggsurvplot(
  surv_fit,
  data = combined_data,
  xlab = "Overall Survival Time (months)",
  ylab = "Survival Probability",
  legend.title = "Patient Mutational Status",
  legend.labs = c(
    "KRAS",
    "KRAS + single-allele TP53 loss (KP+/-)",
    "KRAS + complete TP53 loss (KP)"
  ),
  palette = c("green3", "dodgerblue3", "firebrick3"),
  surv.median.line = "hv",
  pval = TRUE,
  pval.method = TRUE,
  risk.table = FALSE,
  censor.shape = 124,
  censor.size = 2,
  xlim = c(0, 100),
  break.time.by = 25,
  ggtheme = theme_classic(
    base_size = 14,
    base_family = "Arial"
  )
)

# Figure formatting
res$plot <- res$plot +
  theme(
    legend.title = element_text(
      size = 35,
      face = "bold"
    ),
    legend.text = element_text(
      size = 35,
      face = "bold"
    ),
    axis.title = element_text(
      size = 35,
      face = "bold"
    ),
    axis.text = element_text(
      size = 35,
      face = "bold"
    ),
    plot.title = element_text(
      size = 35,
      face = "bold",
      hjust = 0.5
    ),
    plot.caption = element_text(
      size = 35,
      face = "bold"
    )
  )

print(res)


## Survival Statistics to report in the paper
# Median survival time and 95% confidence intervals
summary(surv_fit)$table[
  ,
  c("median", "0.95LCL", "0.95UCL")
]

# Full survival summary
summary(surv_fit)

# Survival probabilities and numbers at risk/events over time
as.data.frame(
  summary(surv_fit)[
    c("time", "n.risk", "n.event", "surv")
  ]
)

# Overall log-rank test using full survival follow-up
survdiff(
  Surv(Time, Status) ~ Strata,
  data = combined_data
)

# Pairwise log-rank comparisons using full survival follow-up
pairwise_survdiff(
  Surv(Time, Status) ~ Strata,
  data = combined_data
)

# Number of patients, events, median survival, and 95% confidence intervals
summary(surv_fit)$table[
  ,
  c(
    "records",
    "events",
    "median",
    "0.95LCL",
    "0.95UCL"
  )
]

# Print Kaplan-Meier fit
surv_fit


##### SUPPLEMENTAL FIGURE 2 #####
#### scRNA-seq preprocessing and lineage identification

### Processed scRNA-seq feature-barcode matrices are available from the published NCBI GEO accession GSEXXXXXX
## Representative preprocessing is shown for sample K_1 
# The same workflow was applied independently to all six biological samples --> K_1, K_2, K_3, KP_1, KP_2, and KP_3 (KP here refers to KP+/-)

# Load packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(scCustomize)
library(SCpubr)
library(flexdashboard)

## Representative preprocessing: K_1
# Download the following processed files for sample K_1 from GSE343294:
  # K_1_barcodes.tsv.gz
  # K_1_features.tsv.gz
  # K_1_matrix.mtx.gz

# Place the files in:
  # data/GSE343294/K_1/

# The files should be renamed to the standard 10x Genomics filenames:
  # barcodes.tsv.gz
  # features.tsv.gz
  # matrix.mtx.gz

# Load the Cell Ranger filtered feature-barcode matrix
data <- Read10X(data.dir = "data/GSE343294/K_1/")

# Create unfiltered object for calculation of the retained cell fraction
K_1_Raw <- CreateSeuratObject(counts = data, min.cells = 0, min.features = 0)

# Create Seurat object for downstream analysis
K_1 <- CreateSeuratObject(counts = data, min.cells = 3,min.features = 200)

# Calculate mitochondrial transcript percentage
K_1[["percent.mt"]] <- PercentageFeatureSet(K_1, pattern = "^mt-")

# Normalize expression and identify variable features
K_1 <- NormalizeData(K_1)
K_1 <- FindVariableFeatures(K_1)

## Quality-control assessment
VlnPlot(K_1,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  )
)

FeatureScatter(
  K_1,
  feature1 = "nCount_RNA", feature2 = "nFeature_RNA"
) +
  ggtitle("Relationship between nCount_RNA and nFeature_RNA") +
  theme_bw()

summary(K_1$nFeature_RNA)
summary(K_1$nCount_RNA)

quantile(
  K_1$nFeature_RNA,
  probs = c(0.95, 0.99, 0.995, 0.999)
)

## Cell filtering and dimensional reduction
# Filter cells based on the QC thresholds used for K_1
K_1 <- subset(
  K_1,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 7500 &
    percent.mt < 15
)

K_1 <- ScaleData(K_1) %>% RunPCA(npcs = 30)

ElbowPlot(K_1, ndims = 30)

K_1 <- FindNeighbors(K_1)
K_1 <- FindClusters(K_1)
K_1 <- RunUMAP(K_1, dims = 1:30)

DimPlot(K_1, label = TRUE)

## Initial lineage-marker assessment
lineage_markers <- c(
  "Cdh1", "Epcam", "Nkx2-1",                                   # epithelial / lung epithelial
  "Hnf4a", "Gkn2", "Krt14", "Krt8", "Krt18", "Sox9",           # epithelial tumor-state / differentiation markers
  "Rtkn2", "Lamp3",                                            # alveolar epithelial markers
  "Scgb1a1", "Krt5", "Foxj1", "Sox2",                          # airway epithelial markers
  "Icam2", "Cdh5", "Pecam1",                                   # endothelial
  "Pdgfra", "Pdgfrb", "Col1a1",                                # mesenchymal / stromal
  "Ptprc",                                                     # pan-immune
  "Cd2", "Cd3e", "Cd3d", "Cd3g",                               # T lymphocytes
  "Cd19", "Cd79a", "Ms4a1",                                    # B lymphocytes
  "Itgax", "Spp1", "Trem2", "Cd64", "Cd14", "Cx3cr1",          # myeloid / monocyte-macrophage-associated
  "Mzb1",                                                      # plasma cells
  "Ccr7",                                                      # activated/migratory immune cells
  "Ly6g",                                                      # neutrophils
  "Itgb3",                                                     # megakaryocyte/platelet-associated and other stromal/myeloid populations
  "Mki67"                                                      # proliferation
)

DotPlot(K_1, features = lineage_markers) + coord_flip()
VlnPlot(K_1, features = c("Cdh1", "Epcam", "Icam2", "Cdh5", "Ptprc","Krt8", "Col1a1", "Pdgfra", "Pdgfrb" ))

## Additional QC metrics
K_1[["percent.rb"]] <- PercentageFeatureSet(K_1, pattern = "^Rp[SL][[:digit:]]|^Rplp[[:digit:]]|^Rpsa")

VlnPlot_scCustom(K_1,
  features = c(
    "nFeature_RNA",
    "percent.mt",
    "percent.rb"
  )
)

K_1[["umi_gene_ratio"]] <- log10(K_1@meta.data$nCount_RNA / K_1@meta.data$nFeature_RNA)

VlnPlot(K_1, features = "umi_gene_ratio")

## Cluster marker identification
Idents(K_1) <- "seurat_clusters"

Markers_K_1 <- FindAllMarkers(
  K_1,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(Markers_K_1, "Markers_K_1.csv", row.names = FALSE)

# Remove ribosomal, mitochondrial, predicted, and other uninformative genes
markers_clean <- Markers_K_1 %>%
  dplyr::filter(
    !grepl(
      "^(Rps|Rpl|Gm|mt-)",
      gene,
      ignore.case = TRUE
    )
  )

Top100_Markers_K_1 <- markers_clean %>%
  dplyr::filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.5
  ) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = avg_log2FC,
    n = 100,
    with_ties = FALSE
  ) %>%
  dplyr::arrange(
    desc(avg_log2FC)
  )

write.csv(Top100_Markers_K_1, "Top100_Markers_K_1_Filtered.csv", row.names = FALSE)

## Iterative removal of low-quality / uninformative clusters
# Remove initial cluster 0
Idents(K_1) <- "seurat_clusters"

K_1_filtered <- subset(K_1,idents = c(1:19))
K_1_filtered <- ScaleData(K_1_filtered)
K_1_filtered <- RunPCA(K_1_filtered)
K_1_filtered <- FindNeighbors(K_1_filtered,dims = 1:30)
K_1_filtered <- FindClusters(K_1_filtered)
K_1_filtered <- RunUMAP(K_1_filtered,dims = 1:30)

DimPlot(K_1_filtered,label = TRUE)

DotPlot(K_1_filtered, features = lineage_markers) + coord_flip()

## Final filtering and reclustering
# Cluster 17 was identified as low-quality / empty droplet and removed after examining gene expression
DEG_CLUSTER_17 <- FindMarkers(
  K_1_filtered,
  ident.1 = 17,
  ident.2 = c(0:16,18:23),
  min.pct = 0.5,
  logfc.threshold = 0.1,
  min.cells.feature = 40
)

write.csv(DEG_CLUSTER_17,file = "DEG_CLUSTER_17.csv")

# Remaining clusters were retained for lineage annotation.
Idents(K_1_filtered) <- "seurat_clusters"

K_1_filtered_X2 <- subset(K_1_filtered, idents = c(0:16, 18:23))
K_1_filtered_X2 <- ScaleData(K_1_filtered_X2)
K_1_filtered_X2 <- RunPCA(K_1_filtered_X2)
K_1_filtered_X2 <- FindNeighbors(K_1_filtered_X2, dims = 1:30)
K_1_filtered_X2 <- FindClusters(K_1_filtered_X2)
K_1_filtered_X2 <- RunUMAP(K_1_filtered_X2, dims = 1:30)

DimPlot(K_1_filtered_X2, label = TRUE)
DotPlot(K_1_filtered_X2, features = lineage_markers) +coord_flip()

## Broad lineage annotation
Idents(K_1_filtered_X2) <- "seurat_clusters"

K_1_filtered_X2$cell_lineage <- K_1_filtered_X2@active.ident

K_1_filtered_X2$cell_lineage <- plyr::mapvalues(
  x = K_1_filtered_X2$cell_lineage,
  from = c(0:23),
  to = c(
    "Immune",
    "Immune",
    "Endothelial",
    "Endothelial",
    "Immune",
    "Epithelial",
    "Endothelial",
    "Epithelial",
    "Epithelial",
    "Immune",
    "Mesenchymal",
    "Immune",
    "Immune",
    "Mesenchymal",
    "Endothelial",
    "Epithelial",
    "Immune",
    "Proliferating Immune",
    "Immune",
    "Endothelial",
    "Endothelial",
    "Immune",
    "Epithelial",
    "Immune"
  )
)

Idents(K_1_filtered_X2) <- "cell_lineage"
DimPlot(K_1_filtered_X2,repel = TRUE,label = TRUE)

## Percentage of cells retained after QC
# Calculate the percentage of cells retained after filtering
QC_ratio <- round(
  ncol(K_1_filtered_X2) / ncol(K_1_Raw) * 100,
  2
)

QC_ratio

# Visualize the retained-cell percentage
flexdashboard::gauge(
  QC_ratio,
  min = 0,
  max = 100,
  symbol = "%",
  gaugeSectors(
    success = c(50, 100),
    danger = c(0, 50)
  )
)

## Merge all six independently processed biological samples
merged_K_KP <- merge(
  K_1_filtered_X2,
  y = c(
    K_2_filtered_X2,
    K_3_filtered_X2,
    KP_1_filtered,
    KP_2_filtered,
    KP_3_filtered
  ),
  add.cell.ids = c(
    "K_1",
    "K_2",
    "K_3",
    "KP_1",
    "KP_2",
    "KP_3"
  )
)

merged_K_KP <- FindVariableFeatures(merged_K_KP, nfeatures = 10000)

merged_K_KP <- ScaleData(merged_K_KP)
merged_K_KP <- RunPCA(merged_K_KP)
merged_K_KP <- FindNeighbors(merged_K_KP, dims = 1:30)
merged_K_KP <- FindClusters(merged_K_KP)
merged_K_KP <- RunUMAP(merged_K_KP,dims = 1:30)

## Assign genotype labels
Idents(merged_K_KP) <- "Sample_ID"

merged_K_KP$Genotype_ID <- merged_K_KP@active.ident

merged_K_KP$Genotype_ID <- plyr::mapvalues(
  x = merged_K_KP$Genotype_ID,
  from = c(
    "K_1",
    "K_2",
    "K_3",
    "KP_1",
    "KP_2",
    "KP_3"
  ),
  to = c(
    "K",
    "K",
    "K",
    "KP+/-",
    "KP+/-",
    "KP+/-"
  )
)

#

##### 2B #####
## Canonical lineage-marker expression

Idents(merged_K_KP) <- "cell_lineage"

DotPlot(
  merged_K_KP,
  features = c(
    "Plvap",
    "Icam2",
    "Cdh5",
    "Ptprc",
    "Cdh1",
    "Epcam",
    "Krt8",
    "Krt18",
    "Col1a1",
    "Col3a1",
    "Pdgfrb",
    "Pdgfra",
    "Mki67"
  )
) +
  scale_color_gradient2(
    low = "#0066FF",
    mid = "#FFFFFF",
    high = "#990000",
    midpoint = 0
  ) +
  RotatedAxis()


## Lineage color palette
lineage_cols <- c(
  "Epithelial" = "#33A02C",
  "Endothelial" = "#4D4D4D",
  "Immune" = "#D95F02",
  "Mesenchymal" = "#7570B3",
  "Proliferating Immune" = "#C71585"
)


##### 2A #####
## UMAP visualization of broad lineages by genotype

p_umap <- DimPlot(
  merged_K_KP,
  reduction = "umap",
  group.by = "cell_lineage",
  split.by = "Genotype_ID",
  cols = lineage_cols,
  pt.size = 1
)

p_umap


##### 2C #####
## Relative lineage composition by genotype

p <- SCpubr::do_BarPlot(
  merged_K_KP,
  group.by = "cell_lineage",
  split.by = "Genotype_ID",
  position = "fill",
  colors.use = lineage_cols,
  flip = FALSE
)

p_bar <- p +
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = lineage_cols,
    drop = FALSE
  ) +
  theme_classic(
    base_size = 14
  ) +
  theme(
    axis.line = element_line(
      color = "black",
      size = 0.6
    ),
    axis.ticks = element_line(
      color = "black",
      size = 0.5
    ),
    axis.text = element_text(
      color = "black"
    ),
    axis.title = element_text(
      face = "bold",
      color = "black"
    ),
    legend.position = "bottom"
  ) +
  ylab("Percentage of cells (%)") +
  xlab("Genotype")

p_bar


##### SUPPLEMENTAL FIGURE 3 #####
### Identification of normal and malignant epithelial populations

## Load packages
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scCustomize)
library(EnhancedVolcano)
library(UCell)
library(viridis)


## Isolate epithelial cells from K and KP+/- scRNA-seq dataset

Idents(merged_K_KP) <- "cell_lineage"
K_KP_Epi <- subset(merged_K_KP, idents = "Epithelial")

K_KP_Epi <- FindVariableFeatures(K_KP_Epi, nfeatures = 10000)
K_KP_Epi <- ScaleData(K_KP_Epi)
K_KP_Epi <- RunPCA(K_KP_Epi)
K_KP_Epi <- FindNeighbors(K_KP_Epi, dims = 1:30, reduction = "pca")
K_KP_Epi <- FindClusters(K_KP_Epi, resolution = 1.5, cluster.name = "unintegrated_clusters")
K_KP_Epi <- RunUMAP(K_KP_Epi, dims = 1:30, reduction = "pca")

Idents(K_KP_Epi) <- "seurat_clusters"

DimPlot_scCustom(K_KP_Epi, label = TRUE)


## Assess residual non-epithelial lineage markers
DotPlot_scCustom(K_KP_Epi, features = c("Cdh1","Epcam","Krt8","Krt18","Sox9","Plvap","Icam2","Cdh5","Lyve1","Vwf","Ptprc","Cd3e","Cd3d","Cd2","Ms4a1","Cd79a","Cd19","Cd14","Cd68","Cd163","Nkg7","Ctsw","Aif1","Axl","Col1a1","Col3a1","Pdgfrb","Pdgfra","Fbln1","Lum","Mki67")) + coord_flip()

## Iterative assessment of ambiguous cluster 24 (potentially mesothelial)
VlnPlot_scCustom(K_KP_Epi, features = c("Ptprc","Cdh5","Runx2","Pecam1","Wt1","Msln","Krt19"))

K_KP_Epi <- JoinLayers(K_KP_Epi)

Cluster_Mesothelial_vs_rest <- FindMarkers(K_KP_Epi, ident.1 = 24, ident.2 = c(0:23,25,26), min.pct = 0.01, logfc.threshold = 0.01, pseudocount.use = 0.001)

boring.genes <- grep("^Rps", rownames(K_KP_Epi@assays$RNA), perl = TRUE, value = TRUE)
ultra.boring.genes <- grep("^Gm", rownames(K_KP_Epi@assays$RNA), perl = TRUE, value = TRUE)
very.boring.genes <- grep("^mt", rownames(K_KP_Epi@assays$RNA), perl = TRUE, value = TRUE)
hated.genes <- grep("^Rpl", rownames(K_KP_Epi@assays$RNA), perl = TRUE, value = TRUE)

Cluster_Mesothelial_vs_rest <- Cluster_Mesothelial_vs_rest[!rownames(Cluster_Mesothelial_vs_rest) %in% boring.genes, ]
Cluster_Mesothelial_vs_rest <- Cluster_Mesothelial_vs_rest[!rownames(Cluster_Mesothelial_vs_rest) %in% ultra.boring.genes, ]
Cluster_Mesothelial_vs_rest <- Cluster_Mesothelial_vs_rest[!rownames(Cluster_Mesothelial_vs_rest) %in% hated.genes, ]
Cluster_Mesothelial_vs_rest <- Cluster_Mesothelial_vs_rest[!rownames(Cluster_Mesothelial_vs_rest) %in% very.boring.genes, ]

EnhancedVolcano(Cluster_Mesothelial_vs_rest, lab = rownames(Cluster_Mesothelial_vs_rest), x = "avg_log2FC", y = "p_val_adj", FCcutoff = 0.5, pCutoff = 10e-5, xlim = c(-max(abs(Cluster_Mesothelial_vs_rest$avg_log2FC)) - 0.1, max(abs(Cluster_Mesothelial_vs_rest$avg_log2FC)) + 0.1), gridlines.major = FALSE, gridlines.minor = FALSE, xlab = "ln(fold change)", colAlpha = 1, col = c("black","green","blue","red"), labSize = 5, pointSize = 1, arrowheads = TRUE, widthConnectors = 0.2, max.overlaps = 15, drawConnectors = TRUE)

write.csv(Cluster_Mesothelial_vs_rest, "Cluster_Mesothelial_vs_rest.csv")

FeaturePlot_scCustom(K_KP_Epi, features = c("Trp63","Krt5"), order = TRUE)

# Cluster 24 showed basal/squamous-associated epithelial features and was retained.

## Remove residual endothelial and immune clusters and recluster
# Cluster 21: endothelial contamination
# Cluster 23: immune contamination

Idents(K_KP_Epi) <- "seurat_clusters"

K_KP_Epi_Filtered <- subset(K_KP_Epi, idents = c(0:20,22,24:26))

K_KP_Epi_Filtered <- FindVariableFeatures(K_KP_Epi_Filtered, nfeatures = 10000)
K_KP_Epi_Filtered <- ScaleData(K_KP_Epi_Filtered)
K_KP_Epi_Filtered <- RunPCA(K_KP_Epi_Filtered)
K_KP_Epi_Filtered <- FindNeighbors(K_KP_Epi_Filtered, dims = 1:30, reduction = "pca")
K_KP_Epi_Filtered <- FindClusters(K_KP_Epi_Filtered, resolution = 1.5, cluster.name = "unintegrated_clusters")
K_KP_Epi_Filtered <- RunUMAP(K_KP_Epi_Filtered, dims = 1:30, reduction = "pca")

Idents(K_KP_Epi_Filtered) <- "seurat_clusters"
DimPlot_scCustom(K_KP_Epi_Filtered, label = TRUE, split.by = "Genotype_ID", pt.size = 2, label.size = 6)

##### 3A #####
## Epithelial-cell UMAP showing final unsupervised clusters

no_green_26_palette <- c("#F6224F","#D95F02","#FA0087","#AA0DFE","#C71585","#7570B3","#FEAF16","#A87E63","#708090","#4D4D4D","#2ED9FF","#F8A19F","#FFFF00","#FF00FF","#6A3D9A","#A1CAF1","#0090FF","#FFD700","#1F78B4","#E31A3C","#B3A2C7","#C49A6C","#9970AB","#FDB462","#D9D9D9","#BC80BD")

K_KP_Epi_Filtered$seurat_clusters <- factor(K_KP_Epi_Filtered$seurat_clusters, levels = 0:25)

p_umap_final <- DimPlot(K_KP_Epi_Filtered, group.by = "seurat_clusters", cols = no_green_26_palette, label = FALSE, pt.size = 1, repel = TRUE) + theme_void()

p_umap_final

ggsave("SuppFig3A_Epithelial_clusters.pdf", p_umap_final, width = 5.2, height = 4)


## Reference-based epithelial annotation using LungMAP
# Load LungMAP mouse epithelial reference atlas

LungMAP_Mouse_Atlas_Epithelial <- readRDS("LungMAP_Mouse_Atlas_Epithelial.rds")

# Set LungMAP cell-type identities
Idents(LungMAP_Mouse_Atlas_Epithelial) <- "celltype_level2"

# Visualize reference annotations
DimPlot(LungMAP_Mouse_Atlas_Epithelial, label = TRUE)

## Transfer LungMAP epithelial cell-type labels to K and KP+/- epithelial cells

ref <- LungMAP_Mouse_Atlas_Epithelial

anchors <- FindTransferAnchors(
  reference = ref,
  query = K_KP_Epi_Filtered,
  dims = 1:30,
  reference.reduction = "pca",
  features = VariableFeatures(K_KP_Epi_Filtered)
)

cell_type_predictions <- TransferData(
  anchorset = anchors,
  refdata = ref$celltype_level2,
  dims = 1:30
)

K_KP_Epi_Filtered$LungMAP_annotation <- cell_type_predictions$predicted.id

## Inspect transferred LungMAP annotations

DimPlot(
  K_KP_Epi_Filtered,
  group.by = "LungMAP_annotation",
  split.by = "Sample_ID",
  label = TRUE
)

Idents(K_KP_Epi_Filtered) <- "LungMAP_annotation"

SCpubr::do_BarPlot(
  K_KP_Epi_Filtered,
  group.by = "LungMAP_annotation",
  split.by = "Genotype_ID",
  position = "fill",
  flip = FALSE
)


##### 3B #####
## LungMAP annotation of K and KP+/- epithelial cells

# Fix plotting order of LungMAP annotations
lungmap_levels <- c(
  "AT1",
  "AT1/AT2",
  "AT2",
  "Basal",
  "Ciliated",
  "PNEC",
  "Secretory"
)

K_KP_Epi_Filtered$LungMAP_annotation <- factor(
  K_KP_Epi_Filtered$LungMAP_annotation,
  levels = lungmap_levels
)

# Colors used for LungMAP annotation
lungmap_cols_new <- c(
  "AT1"       = "#FF4500",
  "AT1/AT2"   = "#E6E6FA",
  "AT2"       = "#87CEEB",
  "Basal"     = "#4B0082",
  "Ciliated"  = "#FFD700",
  "PNEC"      = "#555555",
  "Secretory" = "#F08080"
)

DimPlot(
  K_KP_Epi_Filtered,
  group.by = "LungMAP_annotation",
  cols = lungmap_cols_new,
  label = FALSE,
  repel = TRUE,
  pt.size = 1
)


## Normal epithelial-cell signature scoring using UCell
# Gene signatures representing normal lung epithelial populations were derived from published lung epithelial scRNA-seq datasets

Basal_Cell_signature <- read.table("Basal_Cell_signature.txt")
Ciliated_Cell_Signature <- read.table("Ciliated_Cell_Signature.txt")
Club_Cell_Signature <- read.table("Club_Cell_Signature.txt")
AT1_Cell_Signature <- read.table("AT1_metagene_score_Dalia.txt")
AT2_Cell_Signature <- read.table("AT2_metagene_score_Dalia.txt")
AT2_15wk_Signature <- read.table("Dalia_AT2_Signature_15wk.txt")

AT2_Cell_Signature <- list(AT2_Cell_Signature = AT2_Cell_Signature$V1)
gene.sets <- c(Basal_Cell_signature,Ciliated_Cell_Signature,Club_Cell_Signature,AT1_Cell_Signature,AT2_Cell_Signature,AT2_15wk_Signature)

K_KP_Epi_Filtered <- JoinLayers(K_KP_Epi_Filtered)
K_KP_Epi_Filtered <- AddModuleScore_UCell(K_KP_Epi_Filtered, features = gene.sets)

FeaturePlot_scCustom(K_KP_Epi_Filtered, features = c("Basal_Cell_signature_UCell","Club_Cell_Signature_UCell","Ciliated_Cell_Signature_UCell","AT1_Cell_Signature_UCell","AT2_Cell_Signature_UCell","AT2_15wk_Signature_UCell"), colors_use = viridis_inferno_dark_high, pt.size = 0.5, label = TRUE)


## Summarize UCell signature scores by epithelial cell cluster

obj <- K_KP_Epi_Filtered
Idents(obj) <- "seurat_clusters"

sig_cols <- grep("_UCell$", colnames(obj@meta.data), value = TRUE)

scores_mean <- obj@meta.data |> mutate(cluster = Idents(obj)) |> group_by(cluster) |> summarise(across(all_of(sig_cols), ~mean(.x, na.rm = TRUE)), .groups = "drop")
scores_median <- obj@meta.data |> mutate(cluster = Idents(obj)) |> group_by(cluster) |> summarise(across(all_of(sig_cols), ~median(.x, na.rm = TRUE)), .groups = "drop")

write.csv(scores_mean, "UCell_scores_by_cluster_mean.csv", row.names = FALSE)
write.csv(scores_median, "UCell_scores_by_cluster_median.csv", row.names = FALSE)

## Mean UCell scores across clusters

scores_mean_long <- obj@meta.data %>% mutate(cluster = Idents(obj)) %>% group_by(cluster) %>% summarise(across(all_of(sig_cols), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>% pivot_longer(-cluster, names_to = "signature", values_to = "mean_score") %>% mutate(signature = sub("_UCell$", "", signature))
ggplot(scores_mean_long, aes(x = cluster, y = mean_score, fill = signature)) + geom_col(width = 0.75) + facet_wrap(~signature, scales = "free_y", ncol = 3) + labs(x = "Cluster", y = "Mean UCell", fill = "Signature") + theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# Clusters with the strongest enrichment for normal epithelial signatures:
# Cluster 3  = AT2
# Cluster 13 = Ciliated
# Cluster 17 = AT1
# Cluster 18 = Club


## UCell signature score distributions across epithelial clusters

df <- obj@meta.data |> mutate(cluster = Idents(obj)) |> dplyr::select(cluster, any_of("Sample_ID"), any_of("Genotype_ID"), all_of(sig_cols)) |> pivot_longer(cols = all_of(sig_cols), names_to = "signature", values_to = "score") |> mutate(signature = sub("_UCell$", "", signature))

## Example: Ciliated-cell signature
# The same plotting workflow was used for AT2, AT1 and Club signatures

sig <- "Ciliated_Cell_Signature"

df1 <- df |> filter(signature == sig)
top_val <- df1 |> group_by(cluster) |> summarise(med = mean(score, na.rm = TRUE)) |> summarise(top = max(med)) |> pull(top)
top_clusters <- df1 |> group_by(cluster) |> summarise(med = mean(score, na.rm = TRUE)) |> filter(med == top_val) |> pull(cluster)
df1 <- df1 |> mutate(top_flag = ifelse(cluster %in% top_clusters, "Top", "Other"))

##### 3C-F #####
ggplot(df1, aes(x = cluster, y = score, fill = top_flag, color = top_flag)) + 
    geom_boxplot(outlier.shape = NA, width = 0.7) + 
    geom_jitter(width = 0.2, alpha = 0.25, size = 0.5) + 
    scale_fill_manual(values = c(Other = "#262626", Top = "#4B0082")) + 
    scale_color_manual(values = c(Other = "#262626", Top = "#4B0082")) + 
    labs(title = paste(sig, "UCell score"), x = "Cluster", y = "UCell score") + 
    theme_minimal(base_size = 20, base_family = "Arial") + 
    theme(panel.grid = element_blank(), axis.line = element_line(color = "black", linewidth = 0.5), axis.text = element_text(color = "black", family = "Arial"), axis.title = element_text(color = "black", family = "Arial"), plot.title = element_text(color = "black", family = "Arial", face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")


## Final annotation of K and KP+/- epithelial cell types
# Final epithelial identities were assigned using concordant evidence from UCell normal epithelial signatures, LungMAP label transfer, and marker expression
Idents(K_KP_Epi_Filtered) <- "seurat_clusters"

K_KP_Epi_Filtered$epithelial_cell_type <- K_KP_Epi_Filtered@active.ident
K_KP_Epi_Filtered$epithelial_cell_type <- plyr::mapvalues(x = K_KP_Epi_Filtered$epithelial_cell_type,
                                                          from = c(0:25),
                                                          to   =  c("Malignant","Malignant","Malignant","AT2","Malignant","Malignant","Malignant","Malignant","Malignant","Malignant","Malignant","Malignant","Malignant","Airway_Ciliated","Malignant","Malignant","Malignant","AT1","Airway_Club","Malignant","Malignant","Malignant","Malignant","Malignant","Malignant","Malignant"))

Idents(K_KP_Epi_Filtered) <- 'epithelial_cell_type'
DimPlot(K_KP_Epi_Filtered, label = F, pt.size = 1)

##### 3G #####
## Final normal-versus-malignant epithelial annotation

lungmap_cols_vivid_green <- c(
  "AT1"             = "#1B4F3A",
  "AT2"             = "#DFFF00",
  "Malignant"       = "#7FC97F",
  "Airway_Ciliated" = "#808000",
  "Airway_Club"     = "#00FF7F"
)

DimPlot(
  K_KP_Epi_Filtered,
  reduction = "umap",
  group.by = "epithelial_cell_type",
  cols = lungmap_cols_vivid_green,
  label = TRUE,
  repel = TRUE,
  pt.size = 2
)

## Cell numbers by genotype and epithelial cell type
table(
  K_KP_Epi_Filtered$Genotype_ID,
  K_KP_Epi_Filtered$epithelial_cell_type
)

##### 3H #####
## Marker-expression validation of final epithelial annotations

DotPlot(
  K_KP_Epi_Filtered,
  features = c(
    # AT2
    "Sftpc","Sftpa1","Sftpa2","Sftpb","Abca3","Etv5","Napsa","Slc34a2",
    "Lamp3","Scd1","Kcnj2","Chia1","Lyz2",
    
    # Ciliated
    "Foxj1","Pifo","Tubb4b","Dnah5","Dnah9","Tekt1","Tekt2",
    "Tmem212","Fam183b","Rsph1","Rsph4a","Dynlrb2","Tspan1",
    
    # AT1
    "Ager","Aqp5","Clic5","Emp2","Cavin2","Selenof","Cav1","Cav2",
    "Pdpn","Hopx","Scnn1a",
    
    # Club / secretory
    "Scgb1a1","Cldn10","Scgb3a2","Fmo3","Pon1","Cyp2f2","Muc5b"
  )
) +
  scale_color_gradient2(
    low = "#0066FF",
    mid = "#FFFFFF",
    high = "#990000",
    midpoint = 0
  ) +
  RotatedAxis()


##### FIGURE 1 & SUPPLEMENTAL FIGURE 4 #####
#### Malignant epithelial state  characterization

## Load packages
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scCustomize)
library(harmony)
library(MAST)
library(RColorBrewer)
library(msigdbr)
library(fgsea)
library(qvalue)
library(purrr)
library(stringr)
library(ggrepel)


## Isolate malignant epithelial cells
K_KP_Epi_Filtered <- readRDS("K_KP_Epi_Filtered.rds")

Idents(K_KP_Epi_Filtered) <- "epithelial_cell_type"

K_KP_Malignant <- subset(K_KP_Epi_Filtered, idents = "Malignant")

K_KP_Malignant <- FindVariableFeatures(K_KP_Malignant, nfeatures = 10000)
K_KP_Malignant <- ScaleData(K_KP_Malignant)
K_KP_Malignant <- RunPCA(K_KP_Malignant)
K_KP_Malignant <- FindNeighbors(K_KP_Malignant, dims = 1:30, reduction = "pca")
K_KP_Malignant <- FindClusters(K_KP_Malignant, resolution = 0.5, cluster.name = "unintegrated_clusters")
K_KP_Malignant <- RunUMAP(K_KP_Malignant, dims = 1:30, reduction = "pca")

Idents(K_KP_Malignant) <- "seurat_clusters"
DimPlot_scCustom(K_KP_Malignant, label = TRUE, split.by = "Genotype_ID")


## Harmony integration of malignant epithelial cells across biological samples
K_KP_Malignant[["RNA"]] <- split(K_KP_Malignant[["RNA"]], f = K_KP_Malignant$Sample_ID)

DefaultAssay(K_KP_Malignant) <- "RNA"

K_KP_Malignant_Harmony <- IntegrateLayers(
  object = K_KP_Malignant,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)

K_KP_Malignant_Harmony <- FindNeighbors(K_KP_Malignant_Harmony, reduction = "harmony", dims = 1:30)
K_KP_Malignant_Harmony <- FindClusters(K_KP_Malignant_Harmony, resolution = 1, cluster.name = "harmony_clusters")
K_KP_Malignant_Harmony <- RunUMAP(K_KP_Malignant_Harmony, reduction = "harmony", dims = 1:30)

Idents(K_KP_Malignant_Harmony) <- "harmony_clusters"

DimPlot(K_KP_Malignant_Harmony, label = TRUE, split.by = "Genotype_ID")


## Initial malignant-cluster marker identification
K_KP_Malignant <- JoinLayers(K_KP_Malignant)

Idents(K_KP_Malignant) <- "seurat_clusters"

Markers_K_KP_Malignant <- FindAllMarkers(
  K_KP_Malignant,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "MAST",
  latent.vars = "Genotype_ID"
)

Top100_Markers_K_KP_Malignant <- Markers_K_KP_Malignant %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 100, with_ties = FALSE) %>%
  dplyr::arrange(desc(avg_log2FC))

write.csv(Top100_Markers_K_KP_Malignant, "Top100_Markers_K_KP_Malignant.csv", row.names = FALSE)

# Remove ribosomal, mitochondrial, predicted, and other uninformative genes
markers_clean <- Markers_K_KP_Malignant %>%
  dplyr::filter(!grepl("^(Rps|Rpl|Gm|mt-)", gene, ignore.case = TRUE))

Top100_Markers_K_KP_Malignant_Filtered <- markers_clean %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 100, with_ties = FALSE) %>%
  dplyr::arrange(desc(avg_log2FC))

write.csv(Top100_Markers_K_KP_Malignant_Filtered, "Top100_Markers_K_KP_Malignant_Filtered.csv", row.names = FALSE)


## Iterative malignant-cell filtering and final Harmony reclustering
# Normal epithelial contaminants and low-quality/uninformative clusters were removed iteratively based on cluster marker expression. The retained cells were subsequently reclustered and reintegrated with Harmony.

## Final refinement of malignant epithelial states
# Load malignant epithelial object following initial iterative filtering, Harmony integration, and preliminary cluster annotation.

K_KP_Malignant_Harmony_Filtered_Harmony <- readRDS("K_KP_Malignant_Harmony_Filtered_Harmony.rds")

# Inspect preliminary malignant-state annotations
unique(K_KP_Malignant_Harmony_Filtered_Harmony$clustering_annotation)

## Recode preliminary malignant-state annotations to final manuscript labels
old_col <- "clustering_annotation"

# Copy preliminary annotations into a new metadata column
K_KP_Malignant_Harmony_Filtered_Harmony$cell_type_v2 <-
  K_KP_Malignant_Harmony_Filtered_Harmony[[old_col]][,1]

# Map preliminary annotations to final malignant-state names
recode_map <- c(
  "Neuroepithelial AT2-like" = "Ncoa1+ AT2-like",
  "GI-like" = "Basal-like",
  "Hepatic-Gastric AT2-like" = "Fetub+ Gastric AT2-like",
  "Neuroepithelial" = "Fras1+ ECM_producing",
  "Cryab+ AT1-like" = "Cryab+ AT1-like",
  "Invasive Neuroepithelial" = "PI3K/Akt3+ Ductal-like",
  "Transitional AT2-like" = "Transitional AT2-like",
  "Hepatic-Gastric AT1/AT2" = "Fetub+ Gastric AT2-like",
  "EMT" = "EMT",
  "AT1-like" = "AT1-like",
  "Mixed" = "To_Remove",
  "AT2-like" = "AT2-like",
  "Myc+ Gastric" = "Mdk+ AT2-like",
  "Myc+ Neuroepithelial" = "Nfat5+ neovascularzation",
  "High Cycling" = "High-cycling"
)

old_vals <- K_KP_Malignant_Harmony_Filtered_Harmony$cell_type_v2

K_KP_Malignant_Harmony_Filtered_Harmony$cell_type_v2 <- ifelse(
  old_vals %in% names(recode_map),
  recode_map[old_vals],
  old_vals
)

Idents(K_KP_Malignant_Harmony_Filtered_Harmony) <- "cell_type_v2"

DimPlot(
  K_KP_Malignant_Harmony_Filtered_Harmony,
  group.by = "cell_type_v2",
  label = TRUE
) + NoLegend()


## Remove final low-quality / uninformative malignant population
# The preliminary "Mixed" population could not be confidently annotated and was excluded prior to final malignant-state analyses.

K_KP_Malignant_Harmony_Filtered_Harmony_Filtered <- subset(
  K_KP_Malignant_Harmony_Filtered_Harmony,
  idents = c(
    "Ncoa1+ AT2-like",
    "Basal-like",
    "Fetub+ Gastric AT2-like",
    "Fras1+ ECM_producing",
    "Cryab+ AT1-like",
    "PI3K/Akt3+ Ductal-like",
    "Transitional AT2-like",
    "EMT",
    "AT1-like",
    "AT2-like",
    "Mdk+ AT2-like",
    "Nfat5+ neovascularzation",
    "High-cycling"
  )
)

## Reclustering after removal of the Mixed population

K_KP_Malignant_Harmony_Filtered_Harmony_Filtered <- FindVariableFeatures(K_KP_Malignant_Harmony_Filtered_Harmony_Filtered, selection.method = "vst", nfeatures = 10000)
K_KP_Malignant_Harmony_Filtered_Harmony_Filtered <- ScaleData(K_KP_Malignant_Harmony_Filtered_Harmony_Filtered)
K_KP_Malignant_Harmony_Filtered_Harmony_Filtered <- RunPCA(K_KP_Malignant_Harmony_Filtered_Harmony_Filtered)
K_KP_Malignant_Harmony_Filtered_Harmony_Filtered <- FindNeighbors(K_KP_Malignant_Harmony_Filtered_Harmony_Filtered, dims = 1:30, reduction = "pca")
K_KP_Malignant_Harmony_Filtered_Harmony_Filtered <- FindClusters(K_KP_Malignant_Harmony_Filtered_Harmony_Filtered,resolution = 0.5,cluster.name = "res_0.5")
K_KP_Malignant_Harmony_Filtered_Harmony_Filtered <- RunUMAP(K_KP_Malignant_Harmony_Filtered_Harmony_Filtered, dims = 1:30, reduction = "pca")

Idents(K_KP_Malignant_Harmony_Filtered_Harmony_Filtered) <- "cell_type_v2"

DimPlot(
  K_KP_Malignant_Harmony_Filtered_Harmony_Filtered,
  reduction = "umap",
  split.by = "Genotype_ID",
  label = TRUE,
  label.size = 6,
  pt.size = 1.5
)

## Final Harmony integration across biological samples
# Split RNA assay by biological sample prior to Harmony integration
K_KP_Malignant_Harmony_Filtered_Harmony_Filtered[["RNA"]] <- split(K_KP_Malignant_Harmony_Filtered_Harmony_Filtered[["RNA"]], f = K_KP_Malignant_Harmony_Filtered_Harmony_Filtered$Sample_ID)

obj <- K_KP_Malignant_Harmony_Filtered_Harmony_Filtered
obj <- FindVariableFeatures(obj, assay = "RNA", selection.method = "vst", nfeatures = 20000)
obj <- ScaleData(obj, assay = "RNA")
obj <- RunPCA(obj, assay = "RNA", npcs = 50)

DefaultAssay(obj) <- "RNA"
obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj, resolution = 1, cluster.name = "harmony_clusters")
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30)


##### Fig 1D #####
DimPlot(
  obj,
  group.by = "Genotype_ID",
  cols = c("K" = "#FF0000", "KP+/-" = "#0000FF"),
  pt.size = 1,
  order = c("K", "KP+/-")
)

## Differential expression between KP+/- and K malignant cells
## Gene list used for Metascape/STRING analyses

Idents(obj) <- "Genotype_ID"

obj <- JoinLayers(obj)

DEG_KP_Vs_K <- FindMarkers(
  obj,
  ident.1 = "KP+/-",
  ident.2 = "K",
  min.pct = 0.01,
  logfc.threshold = 0.01,
  pseudocount.use = 0.001,
  test.use = "MAST"
)

all_genes <- rownames(obj[["RNA"]])

boring.to.remove <- unique(c(
  grep("^Rps", all_genes, value = TRUE),
  grep("^Rpl", all_genes, value = TRUE),
  grep("^Gm", all_genes, value = TRUE),
  grep("^mt", all_genes, value = TRUE)
))

DEG_KP_Vs_K <- DEG_KP_Vs_K[!rownames(DEG_KP_Vs_K) %in% boring.to.remove, ]

KP_top_genes <- DEG_KP_Vs_K %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 1.5) %>%
  rownames()

length(KP_top_genes)

write.table(KP_top_genes, "KP_UP_signatures_2.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Background list supplied to Metascape
all_detected_genes <- rownames(obj[["RNA"]])

write.table(all_detected_genes, "Full_Background_List.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

### KP+/- upregulated genes and the complete detected-gene background were uploaded to Metascape. 
## Metascape enrichment outputs were visualized in Cytoscape for Supplementary Fig. 4A; the enrichment bar plot shown in Supplementary Fig 4A is from Metascape outputs. 
# Supplementary figure 4B was generated from the Metascape output.

## Marker expression defining malignant epithelial states
Idents(obj) <- "cell_type_v2"

Markers_Malignant_States <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.4,
  logfc.threshold = 0.4,
  pseudocount.use = 0.001,
  test.use = "MAST",
  latent.vars = "Genotype_ID"
)

markers_clean <- Markers_Malignant_States %>%
  dplyr::filter(!grepl("^(Rps|Rpl|Gm|mt-)", gene, ignore.case = TRUE))

Top250_Markers_Malignant_States <- markers_clean %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 250, with_ties = FALSE) %>%
  dplyr::arrange(desc(avg_log2FC))

write.csv(Top250_Markers_Malignant_States, "Top250_Markers_Malignant_States.csv", row.names = FALSE)

##### Fig 1E #####
DotPlot(
  obj,
  features = c("Fras1","Ncoa1","Rtkn2","Fetub","Lamp3","Msln","Cryab","Sftpc","Akt3","Mki67","Nfat5","Mdk")
) +
  RotatedAxis() +
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))

##### Fig 1F #####
## Final malignant epithelial-state annotation
distinct_tumor_palette <- c(
  "Fras1+ ECM_remodeling"     = "#F8A19F",
  "Ncoa1+ AT2-like"           = "#0000FF",
  "Basal-like"                = "#FA0087",
  "AT1-like"                  = "#18FF2D",
  "Fetub+ Gastric AT2-like"   = "#2ED9FF",
  "Transitional AT2-like"     = "#F6222E",
  "Hypoxic"                   = "#FFFF00",
  "Cryab+ AT1-like"           = "#00A08B",
  "AT2-like"                  = "#FEAF16",
  "PI3K/Akt3+ Ductal-like"    = "#AA0DFE",
  "High-cycling"              = "#A87E63",
  "Nfat5+ neovascularization" = "#708090",
  "Mdk+ AT2-like"             = "#00008B"
)

p_final <- DimPlot(
  obj,
  group.by = "cell_type_v2",
  split.by = "Genotype_ID",
  cols = distinct_tumor_palette,
  label = TRUE,
  repel = TRUE,
  pt.size = 1.5
) + NoLegend() + NoAxes()

p_final

##### Fig 1G #####
## Relative abundance of malignant states by genotype
Idents(obj) <- "cell_type_v2"

p_bar <- SCpubr::do_BarPlot(
  obj,
  group.by = "cell_type_v2",
  split.by = "Genotype_ID",
  position = "fill",
  colors.use = distinct_tumor_palette,
  flip = FALSE
)

p_bar <- p_bar +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_fill_manual(values = distinct_tumor_palette, drop = FALSE) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    legend.position = "bottom"
  ) +
  ylab("Percentage of tumor states") +
  xlab("Genotype")

p_bar

##### Fig 1H #####
## Expression of lung- and gastrointestinal-lineage transcription factors
FeaturePlot_scCustom(obj, features = c("Nkx2-1","Hnf4a"), repel = TRUE)

##### Supp Fig 4D #####
## Number of detected features across malignant epithelial states
Idents(obj) <- "cell_type_v2"

p_vln <- VlnPlot(
  obj,
  features = "nFeature_RNA",
  group.by = "cell_type_v2",
  pt.size = 0,
  cols = distinct_tumor_palette
)

p_vln$layers[[1]]$aes_params$width <- 0.6

p_vln <- p_vln +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    legend.position = "none",
    plot.margin = margin(5.5, 5.5, 25, 5.5)
  ) +
  xlab("Tumor state") +
  ylab("nFeature_RNA")

p_vln

##### Fig 1I #####
## Z-scored malignant-state marker heatmap

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(grid)

## Marker genes used to distinguish malignant epithelial states
my_genes <- c(
  "Hnf1aos1","Cyp2c68","Cubn","Adgrg7",                         # Fras1+ ECM
  "B830012L14Rik","Adam22","St6galnac3","Col23a1",              # Ncoa1+ AT2-like
  "Mif","Nme2","Ppia","Naca",                                   # Basal-like
  "Sema3a","Galnt18","Rtkn2","Col4a4",                          # AT1-like
  "Hp","Apob","Muc1","Sphk1","Igfbp6",                          # Fetub+ Gastric AT2-like
  "Prmt8","Pde7b","Mitf","Itga9",                               # Transitional AT2-like
  "Cldn4","Procr","Gkn2","Tpm2",                                # EMT
  "Spock2","Cavin2","Sparc","Hopx",                             # Cryab+ AT1-like
  "Cxcl15","H2-Aa","Lcn2","Cd74",                               # AT2-like
  "Srgap3","Ampd3","Pfkp","Cyp2s1",                             # PI3K/Akt3+ Ductal-like
  "Aurkb","Top2a","Mki67","Ccna2",                              # High-cycling
  "Glt1d1","Prdm5","Vav3","Mtmr3",                              # Nfat5+ neovascularization
  "Klk8","C3","Marcksl1","Tspan3"                               # Mdk+ AT2-like
)

## Calculate average expression by malignant epithelial state
Idents(obj) <- "cell_type_v2"

avg_exp <- AverageExpression(
  obj,
  features = my_genes,
  group.by = "ident",
  slot = "data"
)$RNA


## Z-score each gene across malignant epithelial states
exp_scaled <- t(
  scale(
    t(
      as.matrix(avg_exp)
    )
  )
)

# Preserve the manually defined marker-gene order
exp_scaled <- exp_scaled[
  my_genes,
  ,
  drop = FALSE
]



## Transpose matrix for final heatmap orientation
# Rows = malignant epithelial states
# Columns = marker genes
exp_flipped <- t(exp_scaled)

## Figure 1I heatmap
pdf(
  "Figure1I_Malignant_State_Heatmap.pdf",
  width = 4.5,
  height = 1.5
)

Heatmap(
  exp_flipped,
  name = "Z-score",
  
  col = colorRamp2(
    c(-2, 0, 2),
    c("#0066FF", "#FFFFFF", "#990000")
  ),
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  # Gene labels
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 3),
  column_names_rot = 45,
  column_names_side = "bottom",
  column_names_centered = FALSE,
  
  # Malignant-state labels
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(
    fontsize = 5,
    fontface = "bold"
  ),
  
  # Heatmap legend
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6),
    labels_gp = gpar(fontsize = 5),
    direction = "horizontal",
    legend_width = unit(2, "cm")
  )
)

dev.off()

##### Figure 1J-L; Supp Fig 4F, H, and J #####
### GSEA of malignant epithelial states
library(Seurat)
library(dplyr)
library(tidyr)
library(msigdbr)
library(fgsea)
library(qvalue)
library(ggplot2)
library(purrr)
library(stringr)

## Generate ranked marker lists for each malignant epithelial state
Idents(obj) <- "cell_type_v2"

# Join RNA layers prior to differential-expression analysis
obj <- JoinLayers(obj)

# Identify genes enriched or depleted in each malignant state. MAST was used with genotype included as a latent variable.
all_markers_0.4 <- FindAllMarkers(
  object = obj,
  assay = "RNA",
  only.pos = FALSE,
  min.pct = 0.4,
  logfc.threshold = 0.4,
  test.use = "MAST",
  latent.vars = "Genotype_ID",
  pseudocount.use = 0.001
)

## Remove ribosomal, mitochondrial, and predicted genes from ranked lists
all_genes <- rownames(obj[["RNA"]])

boring.genes <- grep("^Rps", all_genes, perl = TRUE, value = TRUE)
hated.genes <- grep("^Rpl", all_genes, perl = TRUE, value = TRUE)
ultra.boring.genes <- grep("^Gm", all_genes, perl = TRUE, value = TRUE)
very.boring.genes <- grep("^mt", all_genes, perl = TRUE, value = TRUE)

boring.to.remove <- unique(
  c(
    boring.genes,
    hated.genes,
    ultra.boring.genes,
    very.boring.genes
  )
)


## Load MSigDB gene-set collections
# Hallmark
mouse_H <- msigdbr(
  species = "Mus musculus",
  category = "H"
)

# KEGG
mouse_KEGG <- msigdbr(
  species = "Mus musculus",
  category = "C2",
  subcategory = "CP:KEGG"
)

# WikiPathways
mouse_WIKI <- msigdbr(
  species = "Mus musculus",
  category = "C2",
  subcategory = "CP:WIKIPATHWAYS"
)

# BioCarta
mouse_BIOCARTA <- msigdbr(
  species = "Mus musculus",
  category = "C2",
  subcategory = "CP:BIOCARTA"
)

# Gene Ontology: Biological Process
mouse_GOBP <- msigdbr(
  species = "Mus musculus",
  category = "C5",
  subcategory = "GO:BP"
)

# Gene Ontology: Molecular Function
mouse_GOMF <- msigdbr(
  species = "Mus musculus",
  category = "C5",
  subcategory = "GO:MF"
)

# Gene Ontology: Cellular Component
mouse_GOCC <- msigdbr(
  species = "Mus musculus",
  category = "C5",
  subcategory = "GO:CC"
)

# C8: Cell Type Signatures
mouse_C8 <- msigdbr(
  species = "Mus musculus",
  category = "C8"
)

# C4: Cancer Modules
mouse_C4_CM <- msigdbr(
  species = "Mus musculus",
  category = "C4",
  subcategory = "CM"
)


## Combine MSigDB collections
gene_sets <- rbind(
  mouse_H,
  mouse_KEGG,
  mouse_WIKI,
  mouse_BIOCARTA,
  mouse_GOBP,
  mouse_GOMF,
  mouse_GOCC,
  mouse_C8,
  mouse_C4_CM
)

# Convert MSigDB table to a named list required by fgsea
gene_sets_list <- split(
  x = gene_sets$gene_symbol,
  f = gene_sets$gs_name
)

## Function to perform GSEA for each malignant epithelial state
run_fgsea_for_cluster <- function(marker_df, cluster_id) {
  
  # Remove ribosomal, mitochondrial, and predicted genes
  marker_df <- marker_df[
    !marker_df$gene %in% boring.to.remove,
  ]
  
  # Skip clusters with too few genes after filtering
  if (nrow(marker_df) < 2) {
    warning(
      paste0(
        "Cluster ",
        cluster_id,
        ": not enough genes after filtering. Skipping."
      )
    )
    return(NULL)
  }
  
  # Rank genes by average log2 fold change
  marker_df <- marker_df %>%
    dplyr::arrange(
      dplyr::desc(avg_log2FC)
    )
  
  gene_list <- marker_df$avg_log2FC
  names(gene_list) <- marker_df$gene
  
  # Remove missing values
  gene_list <- gene_list[
    !is.na(gene_list)
  ]
  
  # For duplicated gene symbols, retain the maximum log2 fold change
  gene_list <- tapply(
    gene_list,
    names(gene_list),
    max
  )
  
  gene_list <- sort(
    gene_list,
    decreasing = TRUE
  )
  
  # Run positive-enrichment GSEA
  fgsea_res <- fgsea(
    pathways = gene_sets_list,
    stats = gene_list,
    scoreType = "pos"
  )
  
  # Skip if no pathways are returned
  if (nrow(fgsea_res) == 0) {
    warning(
      paste0(
        "Cluster ",
        cluster_id,
        ": fgsea returned no pathways."
      )
    )
    return(NULL)
  }
  
  # Store leading-edge information
  leadingEdge_list <- fgsea_res$leadingEdge
  
  fgsea_res$GeneSet_Size <- fgsea_res$size
  fgsea_res$Overlap_N <- lengths(leadingEdge_list)
  
  fgsea_res$Overlap_Genes <- vapply(
    leadingEdge_list,
    function(x) paste(x, collapse = ";"),
    ""
  )
  
  
  ## Calculate q-values
  calculate_qvalue <- function(pvals) {
    
    if (length(pvals) == 0) {
      return(numeric(0))
    }
    
    qobj <- tryCatch(
      qvalue(
        pvals,
        lambda = 0.05,
        pi0.method = "bootstrap"
      ),
      error = function(e) NULL
    )
    
    if (inherits(qobj, "qvalue")) {
      return(qobj$qvalues)
    }
    
    return(rep(NA_real_, length(pvals)))
  }
  
  fgsea_res$qval <- calculate_qvalue(
    fgsea_res$pval
  )
  
  # Add malignant-state identity
  fgsea_res$cluster <- cluster_id
  
  
  ## Convert output to a writeable data frame
  
  fgsea_df <- as.data.frame(
    fgsea_res
  )
  
  fgsea_df <- fgsea_df %>%
    mutate(
      across(
        where(is.list),
        ~sapply(
          .,
          function(x) paste(x, collapse = ";")
        )
      )
    )
  
  
  ## Save cluster-specific GSEA results
  out_csv <- paste0(
    "GSEA_cluster_",
    cluster_id,
    ".csv"
  )
  
  write.csv(
    fgsea_df,
    file = out_csv,
    row.names = FALSE
  )
  
  
  ## Generate cluster-specific NES bar plot
  p <- ggplot(
    fgsea_df,
    aes(
      x = reorder(pathway, NES),
      y = NES
    )
  ) +
    geom_col(
      aes(
        fill = padj < 0.05
      )
    ) +
    coord_flip() +
    theme_minimal() +
    labs(
      title = paste0(
        "GSEA NES - ",
        cluster_id
      ),
      x = "Pathway",
      y = "Normalized Enrichment Score (NES)"
    )
  
  ggsave(
    filename = paste0(
      "GSEA_cluster_",
      cluster_id,
      "_NESplot.pdf"
    ),
    plot = p,
    width = 8,
    height = 10,
    units = "in"
  )
  
  return(fgsea_df)
}


## Run GSEA independently for all malignant epithelial states
cluster_list <- split(
  all_markers_0.4,
  all_markers_0.4$cluster
)

fgsea_by_cluster <- lapply(
  names(cluster_list),
  function(clust_id) {
    
    message(
      "Running GSEA for ",
      clust_id,
      " ..."
    )
    
    run_fgsea_for_cluster(
      marker_df = cluster_list[[clust_id]],
      cluster_id = clust_id
    )
  }
)

# Combine all cluster-specific results
fgsea_all <- bind_rows(
  fgsea_by_cluster
)

write.csv(
  fgsea_all,
  file = "GSEA_ALL_CLUSTERS.csv",
  row.names = FALSE
)

## Identify significantly enriched pathways
fgsea_sig <- fgsea_all %>%
  dplyr::filter(
    pval < 0.05 |
      qval < 0.25
  )


## Generate NES matrix across malignant epithelial states
NES_mat <- fgsea_sig %>%
  dplyr::select(
    pathway,
    cluster,
    NES
  ) %>%
  tidyr::pivot_wider(
    names_from = cluster,
    values_from = NES
  )

write.csv(
  NES_mat,
  file = "GSEA_ALL_CLUSTERS_NES_matrix.csv",
  row.names = FALSE
)


## Selected pathways from the corresponding fgsea output were visualized by normalized enrichment score (NES), with bar fill representing adjusted P value (padj)
# Representative GSEA barplot ---> KP+/- vs K Hypoxic GSEA
library(ggplot2)
library(dplyr)
library(stringr)

# Load selected GSEA results for the KP+/- vs K Hypoxic comparison
df <- read.csv(
  "KP_vs_K_Hypoxic_GSEA.csv",
  check.names = FALSE
)

# Remove pathways lacking NES or adjusted P values
df_clean <- df %>%
  filter(!is.na(NES), !is.na(padj)) %>%
  mutate(
    Pathway = str_wrap(as.character(pathway), width = 40),
    GeneCount = Size,
    logPval = -log10(padj)
  )

# Order pathways by NES
ord <- order(df_clean$NES)
df_clean$Pathway <- factor(
  df_clean$Pathway,
  levels = df_clean$Pathway[ord]
)

# Select pathways displayed in the manuscript panel
plot_data <- df_clean %>%
  top_n(15, wt = abs(NES)) %>%
  mutate(
    pathway_short = gsub(
      "GOMF_|MURARO_|HALLMARK_|GOBP_|AIZARANI_|WP_|DESCARTES_|GOCC_|LAKE_|KEGG_|BIOCARTA_|MANNO_",
      "",
      pathway
    ),
    pathway_short = gsub("_", " ", pathway_short),
    pathway_wrapped = str_wrap(pathway_short, width = 25)
  )

# Plot normalized enrichment scores
p_gsea <- ggplot(
  plot_data,
  aes(
    x = reorder(pathway_wrapped, NES),
    y = NES
  )
) +
  geom_bar(
    stat = "identity",
    aes(fill = padj),
    color = "black",
    linewidth = 0.2,
    width = 0.7
  ) +
  coord_flip() +
  scale_fill_gradient(
    low = "#FFFF00",
    high = "#4B0082"
  ) +
  theme_classic(base_size = 7) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(
      color = "black",
      size = 6,
      lineheight = 0.8
    ),
    axis.text.x = element_text(
      size = 6
    ),
    axis.title.x = element_text(
      size = 7,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.key.height = unit(
      0.1,
      "cm"
    ),
    legend.title = element_text(
      size = 6
    ),
    legend.text = element_text(
      size = 5
    ),
    plot.margin = margin(
      5,
      10,
      5,
      5
    )
  ) +
  labs(
    y = "Normalized Enrichment Score (NES)",
    fill = "padj"
  )

p_gsea


##### Supp Fig 4E, I, and K #####
### Genotype-specific volcano plots
# Example: KP+/- versus K Basal-like 

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(MAST)

## Create combined genotype / malignant-state identities
Idents(obj) <- "cell_type_v2"

obj$genotype_cell_type <- paste(
  obj$Genotype_ID,
  obj$cell_type_v2,
  sep = "_"
)

Idents(obj) <- "genotype_cell_type"

# Join RNA layers prior to differential-expression analysis
obj <- JoinLayers(obj)

## Representative comparison: KP+/- versus K Basal-like cells
cluster_KP_vs_K_Basal_like <- FindMarkers(
  obj,
  ident.1 = "KP+/-_Basal-like",
  ident.2 = "K_Basal-like",
  min.pct = 0.25,
  logfc.threshold = 0.25,
  pseudocount.use = 0.001,
  test.use = "MAST"
)


## Remove ribosomal, mitochondrial, and predicted genes
all_genes <- rownames(obj[["RNA"]])

genes_to_remove <- unique(c(
  grep("^Rps", all_genes, perl = TRUE, value = TRUE),
  grep("^Rpl", all_genes, perl = TRUE, value = TRUE),
  grep("^Gm",  all_genes, perl = TRUE, value = TRUE),
  grep("^mt",  all_genes, perl = TRUE, value = TRUE)
))

cluster_KP_vs_K_Basal_like <- cluster_KP_vs_K_Basal_like[
  !rownames(cluster_KP_vs_K_Basal_like) %in% genes_to_remove,
]

# Save differential-expression results
write.csv(
  cluster_KP_vs_K_Basal_like,
  file = "KP_vs_K_Basal_like_MAST_markers.csv",
  row.names = TRUE
)


## KP+/- versus K Basal-like volcano plot
# Genes highlighted in the manuscript panel
genes_to_label <- c(
  "Cdkn2a",
  "Ldha",
  "Pgk1",
  "Pkm",
  "Ybx3",
  "Fam49b",
  "Sox4",
  "Dctpp1",
  "Prmt1",
  "Pgam1",
  "Retnla",
  "Sftpa1",
  "Lyz2",
  "Cxcl15",
  "H2-Aa",
  "Chil1",
  "H2-Ab1",
  "Ly6k",
  "Lamp3"
)

# Convert differential-expression results to data frame
df <- as.data.frame(cluster_KP_vs_K_Basal_like)

# Retain highlighted genes present in the results
genes_to_label <- intersect(
  genes_to_label,
  rownames(df)
)

# Calculate -log10 adjusted P value.
# Replace adjusted P values of zero with a small finite value for plotting.
nz <- min(
  df$p_val_adj[df$p_val_adj > 0],
  na.rm = TRUE
)

df$yplot <- -log10(
  ifelse(
    df$p_val_adj > 0,
    df$p_val_adj,
    nz * 0.1
  )
)

# Highlight selected genes; all other genes are shown in grey
df$color <- ifelse(
  rownames(df) %in% genes_to_label,
  "#FA0087",
  "grey70"
)

# Prepare labels for highlighted genes
label_df <- data.frame(
  x = df[genes_to_label, "avg_log2FC"],
  y = df[genes_to_label, "yplot"],
  label = genes_to_label,
  stringsAsFactors = FALSE
)

p_volcano_basal <- ggplot(
  df,
  aes(
    x = avg_log2FC,
    y = yplot
  )
) +
  
  # Fold-change thresholds
  geom_vline(
    xintercept = c(-0.5, 0.5),
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey30"
  ) +
  
  # Adjusted P-value threshold
  geom_hline(
    yintercept = -log10(1e-5),
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey30"
  ) +
  
  # Genes
  geom_point(
    aes(color = color),
    size = 2.5,
    alpha = 1,
    show.legend = FALSE
  ) +
  
  scale_color_identity() +
  
  # Labels for selected genes
  geom_text_repel(
    data = label_df,
    mapping = aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    color = "black",
    fontface = "italic",
    size = 5,
    segment.color = NA,
    nudge_x = 0.10,
    hjust = 0,
    box.padding = 0.1,
    point.padding = 0.1,
    max.overlaps = Inf
  ) +
  
  labs(
    x = "log2[fold change (KP+/- Basal-like vs K Basal-like)]",
    y = "-log10 adjusted p-value"
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    text = element_text(
      family = "Arial",
      color = "black"
    ),
    axis.text = element_text(
      color = "black"
    )
  )

p_volcano_basal


##### Supp Fig 4G #####
library(Seurat)
library(ggplot2)
library(RColorBrewer)

# 1. Lock the Manuscript Order (1-13)
manuscript_order <- c(
  "Nfat5+ neovascularization", "Fetub+ Gastric AT2-like", "Fras1+ ECM_remodeling",
  "Transitional AT2-like", "PI3K/Akt3+ Ductal-like", "Ncoa1+ AT2-like",
  "Cryab+ AT1-like", "Mdk+ AT2-like", "High-cycling", "Basal-like",
  "AT1-like", "AT2-like", "Hypoxic"
)
Idents(obj) <- "cell_type_v2"
levels(obj) <- manuscript_order

# 2. Use the 10 most repeating genes from your lists
top_repeating_genes <- c(
  "Il15",   # Rank 1
  "Akt3",   # Rank 2
  "Pik3r3", # Rank 3
  "Map2k1", # Rank 4
  "Vav3",   # Rank 5
  "Dclk1",  # Rank 6
  "Raf1",   # Rank 7
  "Myo5a",  # Rank 8
  "Prlr",   # Rank 9
  "Slc7a5"  # Rank 10
)

# 3. Generate the DotPlot
p_dot_top <- DotPlot(obj, features = top_repeating_genes, dot.scale = 8) + 
  RotatedAxis() +
  # Using the Diverging Red-Blue palette for contrast
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  theme_classic(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(face = "italic", size = 9, color = "black"),
    axis.text.y = element_text(size = 8, face = "bold", color = "black"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    panel.grid.major = element_line(color = "grey95") # Adds subtle guide lines
  ) +
  labs(color = "Avg. Exp.", size = "% Exp.")

p_dot_top


##### FIGURE 2 & SUPPLEMENTAL FIGURE 5 #####
## Reciprocal comparison of K/KP+/- and published KP-/- malignant states

# Load packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(scCustomize)
library(dittoSeq)
library(scales)
library(ggalluvial)
library(UCell)
library(ComplexHeatmap)
library(circlize)
library(grid)


## Load published KP-/- reference dataset
KP_Null <- readRDS("KP_Null.rds")

Idents(KP_Null) <- "Cluster.Name"

DimPlot(
  KP_Null,
  group.by = "Cluster.Name",
  label = TRUE,
  label.size = 4,
  repel = TRUE
)


## Published KP-/- tumor-state color palette
tj_kp_palette <- c(
  "Mesenchymal-1" = "#5D4037",
  "Mesenchymal-1 (Met)" = "#FF8C00",
  "Mesenchymal-2" = "#4B0082",
  "Mesenchymal-2 (Met)" = "#00FF7F",
  "AT2-like" = "#D3D3D3",
  "AT1-like" = "#555555",
  "Endoderm-like" = "#E6BEFF",
  "Gastric-like" = "#808000",
  "Late Gastric" = "#008080",
  "Early gastric" = "#FFD700",
  "High plasticity" = "#FA8072",
  "Lung progenitor-like" = "#4169E1",
  "Early EMT-1" = "#8B0000",
  "Early EMT-2" = "#C71585",
  "Pre-EMT" = "#FF1493"
)

DimPlot(
  KP_Null,
  group.by = "Cluster.Name",
  cols = tj_kp_palette,
  label = TRUE,
  repel = TRUE,
  pt.size = 1.2
) +
  NoLegend() +
  NoAxes()


## Preprocess published KP-/- reference for label transfer
KP_Null <- PercentageFeatureSet(
  KP_Null,
  pattern = "^mt-",
  col.name = "percent.mt"
)

# Remove obsolete assay-associated metadata fields
KP_Null[["originalexp"]] <- NULL
KP_Null[["nCount_originalexp"]] <- NULL
KP_Null[["nFeature_originalexp"]] <- NULL

# Normalize and generate PCA/UMAP embeddings required for label transfer
KP_Null_processed <- NormalizeData(KP_Null)
KP_Null_processed <- FindVariableFeatures(KP_Null_processed, selection.method = "vst")
KP_Null_processed <- ScaleData(KP_Null_processed)
KP_Null_processed <- RunPCA(KP_Null_processed)

KP_Null_processed <- FindNeighbors(
  KP_Null_processed,
  dims = 1:30,
  reduction = "pca"
)

KP_Null_processed <- FindClusters(
  KP_Null_processed,
  resolution = 1,
  cluster.name = "unintegrated_clusters"
)

#to generate refUMAP
KP_Null_processed <- RunUMAP(
  KP_Null_processed,
  dims = 1:30,
  reduction = "pca",
  return.model = TRUE
)


##### Fig 2A, Supp Fig 5A #####
## Transfer published KP-/- tumor-state labels to K and KP+/- malignant cells
anchors <- FindTransferAnchors(
  reference = KP_Null_processed,
  query = obj,
  dims = 1:30,
  reference.reduction = "pca"
)

cell_type_predictions <- TransferData(
  anchorset = anchors,
  refdata = KP_Null_processed$Cluster.Name,
  dims = 1:30
)

obj$TJ_annotation <- cell_type_predictions$predicted.id

obj <- AddMetaData(
  obj,
  metadata = cell_type_predictions
)

DimPlot(
  obj,
  group.by = "TJ_annotation",
  split.by = "Genotype_ID",
  cols = tj_kp_palette,
  label = FALSE,
  pt.size = 0.5,
  repel = TRUE
) +
  NoAxes()


## Prepare datasets for integrated comparison
# Rename KP-/- reference annotation column for consistency
colnames(KP_Null_processed@meta.data)[
  colnames(KP_Null_processed@meta.data) == "Cluster.Name"
] <- "TJ_annotation"

# Assign genotype label to published KP-/- cells
KP_Null_processed$Genotype_ID <- "KP-/-"

# Ensure K and KP+/- labels are consistent
obj$Genotype_ID <- factor(
  as.character(obj$Genotype_ID),
  levels = c("K", "KP+/-")
)


## Merge K, KP+/-, and KP-/- malignant cells
K_KP_KPNull_Malignant <- merge(obj,y = KP_Null_processed, add.cell.ids = c("K_KP_Het", "KP_Null"))
K_KP_KPNull_Malignant <- FindVariableFeatures(K_KP_KPNull_Malignant)
K_KP_KPNull_Malignant <- ScaleData(K_KP_KPNull_Malignant)
K_KP_KPNull_Malignant <- RunPCA(K_KP_KPNull_Malignant)
K_KP_KPNull_Malignant <- FindNeighbors(K_KP_KPNull_Malignant, dims = 1:30, reduction = "pca")
K_KP_KPNull_Malignant <- FindClusters(K_KP_KPNull_Malignant, resolution = 0.5, cluster.name = "unintegrated_clusters")
K_KP_KPNull_Malignant <- RunUMAP(K_KP_KPNull_Malignant, dims = 1:30, reduction = "pca")


## Harmony integration across genotypes
K_KP_KPNull_Malignant[["RNA"]] <- split(K_KP_KPNull_Malignant[["RNA"]],f = K_KP_KPNull_Malignant$Genotype_ID)

DefaultAssay(K_KP_KPNull_Malignant) <- "RNA"

K_KP_KPNull_Malignant_Harmony <- IntegrateLayers(object = K_KP_KPNull_Malignant,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)

K_KP_KPNull_Malignant_Harmony <- FindNeighbors(K_KP_KPNull_Malignant_Harmony, reduction = "harmony", dims = 1:30)
K_KP_KPNull_Malignant_Harmony <- FindClusters(K_KP_KPNull_Malignant_Harmony, resolution = 1.5, cluster.name = "harmony_clusters")
K_KP_KPNull_Malignant_Harmony <- RunUMAP(K_KP_KPNull_Malignant_Harmony, reduction = "harmony", dims = 1:30)

DimPlot(
  K_KP_KPNull_Malignant_Harmony,
  group.by = "TJ_annotation",
  split.by = "Genotype_ID",
  label = TRUE
)


##### Supp Fig 5B #####
K_KP_KPNull_Malignant_Harmony$Genotype_ID <- factor(
  as.character(K_KP_KPNull_Malignant_Harmony$Genotype_ID),
  levels = c("K", "KP+/-", "KP-/-")
)

Idents(K_KP_KPNull_Malignant_Harmony) <- "Genotype_ID"

cells_ordered <- c(
  colnames(K_KP_KPNull_Malignant_Harmony)[
    K_KP_KPNull_Malignant_Harmony$Genotype_ID == "K"
  ],
  colnames(K_KP_KPNull_Malignant_Harmony)[
    K_KP_KPNull_Malignant_Harmony$Genotype_ID == "KP+/-"
  ],
  colnames(K_KP_KPNull_Malignant_Harmony)[
    K_KP_KPNull_Malignant_Harmony$Genotype_ID == "KP-/-"
  ]
)

dittoBarPlot(
  object = K_KP_KPNull_Malignant_Harmony,
  var = "TJ_annotation",
  group.by = "Genotype_ID",
  color.panel = tj_kp_palette,
  cells.use = cells_ordered
) +
  scale_y_continuous(
    labels = label_percent(),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    legend.position = "bottom"
  ) +
  ylab("Percentage of Tumor States") +
  xlab("Genotype")


##### Fig 2B #####
## Correspondence between published KP-/- states and K + KP+/- malignant states

# Generate cell-count table used for alluvial/Sankey visualization
df_sankey <- as.data.frame(
  table(
    TJ_annotation = obj$TJ_annotation,
    cell_type_v2 = obj$cell_type_v2
  )
)

colnames(df_sankey)[3] <- "Freq"


## K/KP+/- malignant-state palette
distinct_tumor_palette <- c(
  "Fras1+ ECM_producing" = "#F8A19F",
  "Ncoa1+ AT2-like" = "#0000FF",
  "Basal-like" = "#FA0087",
  "AT1-like" = "#18FF2D",
  "Fetub+ Gastric AT2-like" = "#2ED9FF",
  "Transitional AT2-like" = "#F6222E",
  "EMT" = "#FFFF00",
  "Cryab+ AT1-like" = "#00A08B",
  "AT2-like" = "#FEAF16",
  "PI3K/Akt3+ Ductal-like" = "#AA0DFE",
  "High-cycling" = "#A87E63",
  "Nfat5+ neovascularzation" = "#708090",
  "Mdk+ AT2-like" = "#00008B"
)

ggplot(
  df_sankey,
  aes(
    axis1 = TJ_annotation,
    axis2 = cell_type_v2,
    y = Freq
  )
) +
  geom_alluvium(
    aes(fill = cell_type_v2),
    width = 0.25,
    alpha = 0.8
  ) +
  geom_stratum(
    width = 0.25,
    fill = "grey95",
    color = "black"
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 2.2,
    min.y = 50
  ) +
  scale_fill_manual(
    values = distinct_tumor_palette
  ) +
  scale_x_discrete(
    limits = c(
      "KP-/- Reference",
      "K/KP+/- Final"
    ),
    expand = c(0.1, 0.1)
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 8)
  )


##### Fig 2C #####
## Score K and KP+/- malignant-state signatures in published KP-/- tumor states
# tumor_signatures contains gene signatures derived from the malignant states identified in the K and KP+/- dataset generated in this study

KP_Null_processed <- AddModuleScore_UCell(
  KP_Null_processed,
  features = tumor_signatures,
  name = NULL
)

scores_df <- KP_Null_processed@meta.data %>%
  group_by(TJ_annotation) %>%
  summarise(
    across(
      all_of(names(tumor_signatures)),
      ~mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  as.data.frame()

rownames(scores_df) <- scores_df$TJ_annotation
scores_df$TJ_annotation <- NULL


## Z-score tumor-signature enrichment across published KP-/- states
# Rows = K and KP+/- malignant-state signatures; Columns = published KP-/- tumor states

plot_matrix <- t(scale(scores_df))

## Define final manuscript ordering
column_order <- c(
  "AT1-like",
  "AT2-like",
  "Lung progenitor-like",
  "Early gastric",
  "Gastric-like",
  "Late Gastric",
  "Endoderm-like",
  "High plasticity",
  "Pre-EMT",
  "Early EMT-1",
  "Early EMT-2",
  "Mesenchymal-2",
  "Mesenchymal-2 (Met)",
  "Mesenchymal-1",
  "Mesenchymal-1 (Met)"
)

row_order <- c(
  "AT1-like",
  "Cryab+ AT1-like",
  "AT2-like",
  "Mdk+ AT2-like",
  "Ncoa1+ AT2-like",
  "Transitional AT2-like",
  "Fetub+ Gastric AT2-like",
  "Fras1+ ECM_remodeling",
  "Nfat5+ Neovascularization",
  "PI3K/Akt+ Ductal-like",
  "Hypoxic",
  "Basal-like",
  "High_cycling"
)

plot_matrix_ordered <- plot_matrix[
  row_order,
  column_order,
  drop = FALSE
]

## K and KP+/- tumor-signature heatmap across KP-/- reference states
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#0066FF", "#FFFFFF", "#990000")
)

Heatmap(
  plot_matrix_ordered,
  name = "Z-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(
    fontsize = 10,
    fontfamily = "sans",
    col = "black"
  ),
  column_names_gp = gpar(
    fontsize = 10,
    fontfamily = "sans",
    col = "black"
  ),
  column_names_rot = 45,
  border = "black",
  rect_gp = gpar(
    col = "white",
    lwd = 0.5
  ),
  heatmap_legend_param = list(
    title_gp = gpar(
      fontsize = 10,
      fontfamily = "sans",
      fontface = "bold"
    )
  )
)


##### Fig 2D-F, Supp Fig 5C-D #####
## Reciprocal label transfer: K and KP+/- reference -> published KP-/- query

# The reciprocal analysis was performed using the same transfer, projection, correspondence, and signature-scoring workflow shown above, with the reference and query datasets exchanged

# Reference:
obj_LT <- RunUMAP(obj, dims = 1:30, 
                  reduction = "pca",return.model = T)
reference <- obj_LT

# Query:
query <- KP_Null_processed

# Reciprocal label transfer:
anchors_reciprocal <- FindTransferAnchors(reference = reference, query = query, 
                                    dims = 1:20, reference.reduction = "pca")

cell_type_predictions <- TransferData(anchorset = anchors_reciprocal, refdata = obj_LT$cell_type_v2, dims = 1:20)

query$MV_annotation <-cell_type_predictions$predicted.id


## Projection, Sankey correspondence, relative state composition, and reciprocal signature-scoring analyses were generated using the same procedures above


##### FIGURE 3 & SUPPLEMENTAL FIGURE 6 #####
## Representative inferCNV workflow: The complete workflow below is shown for the published KP-/- malignant epithelial dataset
# The same inferCNV workflow and parameters were applied independently to K and KP+/- malignant cells using the corresponding genotype-specific malignant-cell input objects

## Load packages
library(Seurat)
library(dplyr)
library(infercnv)
library(readr)
library(ggplot2)
library(ggpubr)

## Prepare normal aged Tabula Muris T-cell reference population
# Load the aged Tabula Muris/Tabula Senis lung reference object --> T cells from 18- and 21-month-old mice were used as the diploid reference population for inferCNV

Tabula_Muris_Lung_18m_21m <- readRDS(
  "Tabula_Muris_Lung_18m_21m.rds"
)

# Subset immune cells
Idents(Tabula_Muris_Lung_18m_21m) <- "cell_lineage"

Tabula_Muris_Lung_18m_21m_Imm <- subset(
  Tabula_Muris_Lung_18m_21m,
  idents = "Immune"
)

# Inspect immune-cell annotations
unique(
  Tabula_Muris_Lung_18m_21m_Imm$cell_type
)

## Subset Tabula Muris T-cell populations
# T-cell populations used for the normal diploid reference:
  # "T cell"
  # "CD4-positive, alpha-beta T cell"
  # "regulatory T cell"
  # "CD8-positive, alpha-beta T cell"

Idents(Tabula_Muris_Lung_18m_21m_Imm) <- "cell_type"

TM_T_cells <- subset(
  Tabula_Muris_Lung_18m_21m_Imm,
  idents = c(
    "T cell",
    "CD4-positive, alpha-beta T cell",
    "regulatory T cell",
    "CD8-positive, alpha-beta T cell"
  )
)

## Assign common inferCNV reference annotation
# All selected Tabula Muris T-cell populations were pooled as a single normal diploid reference population for inferCNV.

TM_T_cells$inferCNV_annotation <- "Tabula_Muris_T"

Idents(TM_T_cells) <- "inferCNV_annotation"

table(TM_T_cells$inferCNV_annotation)

DimPlot(
  TM_T_cells,
  group.by = "inferCNV_annotation",
  label = TRUE
)

## Randomly subsample published KP-/- malignant cells
set.seed(1)

N <- 5000

cells_keep <- sample(
  colnames(KP_Null),
  size = min(N, ncol(KP_Null)),
  replace = FALSE
)

KP_Null_subsampled <- subset(
  KP_Null,
  cells = cells_keep
)

table(KP_Null_subsampled$Cluster.Name)

# Subsampled and full-dataset inferCNV analyses produced concordant CNV patterns; the subsampled dataset was used for the workflow shown here.


## Prepare KP-/- malignant-cell annotations
colnames(KP_Null_subsampled@meta.data)[
  colnames(KP_Null_subsampled@meta.data) == "Cluster.Name"
] <- "inferCNV_annotation"

Idents(KP_Null_subsampled) <- "inferCNV_annotation"

DimPlot(KP_Null_subsampled, group.by = "inferCNV_annotation")

## Merge KP-/- malignant cells with diploid T-cell reference
Annotated_TM_T_KP_Null_Subsampled <- merge(KP_Null_subsampled, y = TM_T_cells, add.cell.ids = c("KP_Null_subsampled","TM_T_cells"))
Annotated_TM_T_KP_Null_Subsampled <- ScaleData(Annotated_TM_T_KP_Null_Subsampled)
Annotated_TM_T_KP_Null_Subsampled <- RunPCA(Annotated_TM_T_KP_Null_Subsampled)
Annotated_TM_T_KP_Null_Subsampled <- FindNeighbors(Annotated_TM_T_KP_Null_Subsampled,dims = 1:30)
Annotated_TM_T_KP_Null_Subsampled <- FindClusters(Annotated_TM_T_KP_Null_Subsampled, resolution = 0.8)
Annotated_TM_T_KP_Null_Subsampled <- RunUMAP(Annotated_TM_T_KP_Null_Subsampled,dims = 1:30)

Idents(Annotated_TM_T_KP_Null_Subsampled) <- "inferCNV_annotation"
DimPlot(Annotated_TM_T_KP_Null_Subsampled,group.by = "inferCNV_annotation")


# Join RNA layers before exporting the inferCNV input matrix
Annotated_TM_T_KP_Null_Subsampled <- JoinLayers(Annotated_TM_T_KP_Null_Subsampled)

## Generate inferCNV input files
raw <- Annotated_TM_T_KP_Null_Subsampled

# Raw count matrix
exp.rawdata <- as.matrix(
  raw[["RNA"]]$counts
)

write.table(
  exp.rawdata,
  "input_count_matrix.txt",
  row.names = TRUE,
  quote = FALSE,
  sep = "\t"
)

# Cell annotation file
raw_metadata <- raw[[]]

anno <- raw_metadata %>%
  dplyr::select(
    inferCNV_annotation
  )

colnames(anno) <- NULL

write.table(
  anno,
  "input_cell_annotation.txt",
  row.names = TRUE,
  quote = FALSE,
  sep = "\t"
)


## Gene-position annotation for inferCNV
# Genomic positions were obtained from the inferCNV/Trinity CTAT
  # mouse GRCm39 gene-position annotation.

full.anno <- read.delim(
  "mouse_gencode.GRCm39.vM32.basic.annotation.by_gene_name.infercnv_positions.txt",
  header = FALSE,
  stringsAsFactors = TRUE
)

colnames(full.anno) <- c(
  "gene_name",
  "chr",
  "start",
  "end"
)

rownames(full.anno) <- full.anno$gene_name
colnames(full.anno) <- NULL

write.table(
  full.anno,
  "input_gene_annotation.txt",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

## Create inferCNV object
numcores <- 7

initial_infercnv_object <- CreateInfercnvObject(
  raw_counts_matrix = "input_count_matrix.txt",
  annotations_file = "input_cell_annotation.txt",
  delim = "\t",
  gene_order_file = "input_gene_annotation.txt",
  ref_group_names = c(
    "Tabula_Muris_T"
  )
)

write_rds(
  initial_infercnv_object,
  "initial_infercnv_object.rds"
)

## Run inferCNV
start_time <- Sys.time()

infercnv_output <- infercnv::run(
  initial_infercnv_object,
  num_threads = numcores - 1,
  out_dir = "infercnv_output",
  cutoff = 0.1,
  window_length = 101,
  max_centered_threshold = 3,
  cluster_by_groups = TRUE,
  plot_steps = FALSE,
  HMM = TRUE,
  denoise = TRUE,
  sd_amplifier = 1.3,
  write_expr_matrix = FALSE,
  analysis_mode = "samples",
  png_res = 500
)

end_time <- Sys.time()

print(start_time)
print(end_time)
print(end_time - start_time)

write_rds(
  infercnv_output,
  "infercnv_output.rds"
)

## Export inferCNV expression matrix
infercnvobs <- as.data.frame(
  infercnv_output@expr.data
)

write.table(
  infercnvobs,
  "infercnvobs.txt",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)


##### Fig 3A-D #####
## Plot inferred copy-number profiles
plot_cnv(
  infercnv_output,
  out_dir = ".",
  output_filename = "KP_Null_inferCNV",
  output_format = "png",
  png_res = 500
)

##### Supp Fig 6A-C #####
## Genome-wide absolute inferred CNV score
# The same calculation and visualization workflow was applied independently to K and KP+/- cells using their corresponding inferCNV outputs and genotype-specific Seurat objects

## The representative workflow below is shown for the published KP-/- cells:
# Load final inferCNV output
infercnv_obj <- readRDS(
  "KP_Null_infercnv_output.rds"
)

# Extract final inferCNV expression-intensity matrix.
# inferCNV values are centered around 1, where:
  #   1   = inferred diploid state
  #   > 1 = inferred copy-number gain
  #   < 1 = inferred copy-number loss

expr_data <- infercnv_obj@expr.data

## Calculate genome-wide absolute CNV score for each cell
# Calculate the absolute deviation from the inferred diploid baseline
abs_deviation_matrix <- abs(expr_data - 1)

# Sum absolute deviations across all genomic features for each cell
absolute_cnv_score <- colSums(abs_deviation_matrix)


## Match inferCNV cell names to the original KP-/- Seurat object
# Remove the prefix introduced when the malignant and reference populations were merged prior to inferCNV analysis.
names(absolute_cnv_score) <- gsub(
  "^KP_Null_subsampled_",
  "",
  names(absolute_cnv_score)
)

# Identify cells present in both the inferCNV output and Seurat object
common_cells <- intersect(
  names(absolute_cnv_score),
  colnames(KP_Null_subsampled)
)

message(
  "Found ",
  length(common_cells),
  " matching cells."
)


## Add absolute CNV score to Seurat metadata
KP_Null_subsampled <- AddMetaData(
  KP_Null_subsampled,
  metadata = absolute_cnv_score[common_cells],
  col.name = "abs_cnv_score"
)

summary(KP_Null_subsampled$abs_cnv_score)

## Generate visualization-specific CNV score
# Values above 800 were clipped to 800 for visualization only;  underlying absolute CNV scores were not modified.

KP_Null_subsampled$abs_cnv_visual <- pmin(
  KP_Null_subsampled$abs_cnv_score,
  800
)

## Representative absolute-CNV UMAP
FeaturePlot(
  KP_Null_subsampled,
  features = "abs_cnv_visual",
  pt.size = 1.2,
  label = TRUE
) +
  scale_colour_gradientn(
    colors = c(
      "#0000FF",
      "lightgrey",
      "#FFCC00",
      "#FF0000",
      "#800000"
    ),
    values = scales::rescale(
      c(0, 300, 400, 600, 800)
    ),
    name = "CNV Score"
  ) +
  ggtitle(
    "Absolute CNV Magnitude (|x - 1|)"
  ) +
  theme_classic()


##### Fig 3E #####
## Inferred Trp53 copy-number score across K and KP+/- malignant epithelial states
# InferCNV was run jointly on K and KP+/- malignant cells using the same workflow and parameters described above --> The inferred value at the Trp53 locus was extracted for each cell and mapped back to the final Seurat object (obj)

## Extract inferred Trp53 copy-number value
# Remove prefixes introduced during construction of the inferCNV input object
colnames(expr_data) <- gsub(
  "^K_KPHet_Tumor_|^TM_T_cells_",
  "",
  colnames(expr_data)
)

# Identify cells shared between the inferCNV output and Seurat object
common_cells <- intersect(
  colnames(expr_data),
  colnames(obj)
)

# Extract the inferred Trp53 value for matched malignant cells
Trp53_row <- expr_data[
  "Trp53",
  common_cells
]

# Add inferred Trp53 copy-number value to Seurat metadata
obj <- AddMetaData(
  obj,
  metadata = Trp53_row,
  col.name = "Trp53_genomic_val"
)


### Trp53 inferred copy-number score across malignant epithelial states
Idents(obj) <- "cell_type_v2"

VlnPlot(
  obj,
  features = "Trp53_genomic_val",
  group.by = "cell_type_v2",
  cols = distinct_tumor_palette,
  pt.size = 0
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    x = "Tumor state",
    y = "Inferred Trp53 copy-number score"
  )

##### Supp Fig 6D #####
## Compare inferred Trp53 copy-number scores between K and KP+/- cells within each malignant epithelial state

plot_data <- FetchData(
  obj,
  vars = c(
    "Trp53_genomic_val",
    "cell_type_v2",
    "Genotype_ID"
  )
)

# Order tumor states by median inferred Trp53 value
plot_data$cell_type_v2 <- reorder(
  plot_data$cell_type_v2,
  plot_data$Trp53_genomic_val,
  median
)

## Calculate mean Trp53 score and standard error by genotype/state
summary_df <- plot_data %>%
  group_by(
    cell_type_v2,
    Genotype_ID
  ) %>%
  summarise(
    mean_val = mean(
      Trp53_genomic_val,
      na.rm = TRUE
    ),
    se_val = sd(
      Trp53_genomic_val,
      na.rm = TRUE
    ) / sqrt(
      sum(!is.na(Trp53_genomic_val))
    ),
    .groups = "drop"
  )

## Plot genotype-specific Trp53 copy-number scores
ggplot(
  summary_df,
  aes(
    x = cell_type_v2,
    y = mean_val,
    fill = Genotype_ID
  )
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(0.8),
    width = 0.7
  ) +
  geom_errorbar(
    aes(
      ymin = mean_val - se_val,
      ymax = mean_val + se_val
    ),
    width = 0.2,
    position = position_dodge(0.8)
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    alpha = 0.5
  ) +
  stat_compare_means(
    data = plot_data,
    aes(
      x = cell_type_v2,
      y = Trp53_genomic_val,
      group = Genotype_ID
    ),
    label = "p.signif",
    method = "wilcox.test",
    inherit.aes = FALSE,
    label.y = 1.02,
    size = 6
  ) +
  scale_fill_manual(
    values = c(
      "K" = "#FF4B4B",
      "KP+/-" = "#4B4BFF"
    )
  ) +
  coord_cartesian(
    ylim = c(0.85, 1.05)
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    x = "Tumor state",
    y = "Mean inferred Trp53 copy-number score (+/- SE)"
  )

##### Supp Fig 6E-F #####
### Inferred Kras copy-number score
## Kras copy-number values were extracted and analyzed using the same workflow shown above for Trp53

# Kras was analyzed identically to Trp53 by substituting:
Kras_row <- expr_data["Kras", common_cells]

obj <- AddMetaData(
  obj,
  metadata = Kras_row,
  col.name = "Kras_genomic_val"
)

## Subsequent visualization and statistical analysis were performed as shown above for Trp53, using "Kras_genomic_val" in place of "Trp53_genomic_val"
# The corresponding "Kras" row of the inferCNV expression matrix was extracted for each cell and mapped back to the Seurat object and the same plotting was then applied

##### FIGURE 4 & SUPPLEMENTARY FIGURE 7 #####

## CytoTRACE analysis of malignant epithelial states
#The complete workflow below is shown for the combined K and KP+/- malignant epithelial dataset. The same workflow was applied independently to the published KP-/- malignant epithelial dataset.

library(Seurat)
library(CytoTRACE)
library(dplyr)
library(ggplot2)
library(ggpubr)


# CytoTRACE was run using non-normalized RNA counts
counts <- as.matrix(
  GetAssayData(
    obj,
    assay = "RNA",
    slot = "counts"
  )
)

results <- CytoTRACE(
  counts,
  ncores = 8
)


## Add CytoTRACE scores to obj metadata
obj$CytoTRACE_score <- results$CytoTRACE[
  match(
    colnames(obj),
    names(results$CytoTRACE)
  )
]

head(
  obj@meta.data[
    ,
    c(
      "cell_type_v2",
      "Genotype_ID",
      "CytoTRACE_score"
    )
  ]
)


##### Fig 4A #####
## Compare CytoTRACE scores between K and KP+/- cells within tumor states
ggplot(
  obj@meta.data,
  aes(
    x = cell_type_v2,
    y = CytoTRACE_score,
    fill = Genotype_ID
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.8,
    color = "black",
    linewidth = 0.6
  ) +
  stat_compare_means(
    aes(group = Genotype_ID),
    label = "p.signif",
    method = "wilcox.test",
    label.y = 1.02
  ) +
  scale_fill_manual(
    values = c(
      "K" = "#FF0000",
      "KP+/-" = "#0000FF"
    )
  ) +
  ylim(0, 1.1) +
  theme_classic() +
  theme(
    axis.text = element_text(
      color = "black",
      size = 12
    ),
    axis.title = element_text(
      color = "black",
      face = "bold",
      size = 14
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.8
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.text = element_text(
      color = "black"
    ),
    legend.title = element_text(
      color = "black",
      face = "bold"
    )
  ) +
  labs(
    x = "Tumor state",
    y = "CytoTRACE Score"
  )

## Summary of genotype-specific statistical comparisons
stats_results <- obj@meta.data %>%
  group_by(cell_type_v2) %>%
  filter(
    n_distinct(Genotype_ID) == 2
  ) %>%
  summarise(
    p_val = wilcox.test(
      CytoTRACE_score ~ Genotype_ID
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      p_val < 0.0001 ~ "****",
      p_val < 0.001  ~ "***",
      p_val < 0.01   ~ "**",
      p_val < 0.05   ~ "*",
      TRUE           ~ "ns"
    )
  )

print(stats_results)


##### Fig 4B #####

## CytoTRACE scores across published KP-/- malignant epithelial states
# The same CytoTRACE workflow shown above was applied independently to the published KP-/- malignant epithelial object using its raw RNA count matrix

##### Fig 4C ####
## Relationship between copy-number burden and plasticity across tumor states
# Mean absolute inferCNV copy number burden and mean CytoTRACE plasticity score were calculated for each K/KP+/- malignant epithelial state and a linear regression was fit across tumor-state means.


library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)


## Extract cell-level CytoTRACE and inferCNV scores
plot_df <- FetchData(
  obj,
  vars = c(
    "cell_type_v2",
    "CytoTRACE_score",
    "abs_cnv_score"
  )
)


## Calculate mean CNV burden and CytoTRACE score for each tumor state
cluster_summary <- plot_df %>%
  group_by(cell_type_v2) %>%
  summarise(
    mean_CNV = mean(
      abs_cnv_score,
      na.rm = TRUE
    ),
    mean_CytoTRACE = mean(
      CytoTRACE_score,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


## Plot relationship between tumor-state CNV burden and plasticity
ggplot(
  cluster_summary,
  aes(
    x = mean_CNV,
    y = mean_CytoTRACE
  )
) +
  
  # Linear regression with 95% confidence interval
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    linetype = "dashed",
    color = "black",
    fill = "grey70",
    alpha = 0.5,
    linewidth = 0.8,
    se = TRUE
  ) +
  
  # Tumor-state means
  geom_point(
    size = 3,
    color = "darkred"
  ) +
  
  # Tumor-state labels
  geom_text_repel(
    aes(
      label = cell_type_v2
    ),
    size = 3,
    color = "black",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = NA,
    show.legend = FALSE
  ) +
  
  theme_classic() +
  
  theme(
    axis.text = element_text(
      color = "black"
    ),
    axis.title = element_text(
      color = "black",
      face = "bold"
    )
  ) +
  
  labs(
    x = "Copy Number Burden Score",
    y = "Plasticity Score"
  )


##### Fig 4D  #####
## Inferred transcription-factor activity using decoupleR / CollecTRI
# TF activities were inferred from the final K and KP+/- malignant seurat  object using the CollecTRI transcription-factor regulatory network and the weighted-mean (wmean) method implemented in decoupleR

library(Seurat)
library(decoupleR)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)
library(scCustomize)


# Join RNA layers prior to activity inference
obj <- JoinLayers(obj)

# Use final malignant-state annotations
Idents(obj) <- "cell_type_v2"


## Load CollecTRI regulatory network

# The CollecTRI regulatory network is publicly available from OmniPath. During the original analysis, online retrieval through R was unsuccessful; therefore, the publicly available CollecTRI interaction table was downloaded and read locally.

net_raw <- readr::read_tsv(
  "Collectri_raw.txt",
  show_col_types = FALSE
)

## Prepare CollecTRI mode-of-regulation information

net_processed <- net_raw %>%
  mutate(
    mor = case_when(
      is_stimulation == "True" | is_stimulation == 1 ~ 1,
      is_inhibition  == "True" | is_inhibition  == 1 ~ -1,
      TRUE ~ 1
    )
  ) %>%
  mutate(
    source_id = gsub(
      "COMPLEX:",
      "",
      source
    ),
    target_id = gsub(
      "COMPLEX:",
      "",
      target
    )
  ) %>%
  mutate(
    # For multi-protein identifiers, retain the first listed component
    source_id = sub(
      "_.*",
      "",
      source_id
    ),
    target_id = sub(
      "_.*",
      "",
      target_id
    )
  )


## CollecTRI identifier mapping and mouse gene-symbol formatting

# Extract unique UniProt identifiers from the CollecTRI network
all_uniprots <- unique(
  c(
    net_processed$source_id,
    net_processed$target_id
  )
)

# Map UniProt identifiers to gene symbols using the local annotation database
id_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = all_uniprots,
  columns = "SYMBOL",
  keytype = "UNIPROT"
)

# Format mapped symbols to mouse-style capitalization
id_map <- id_map %>%
  filter(
    !is.na(SYMBOL)
  ) %>%
  mutate(
    mouse_symbol = paste0(
      substr(
        SYMBOL,
        1,
        1
      ),
      tolower(
        substr(
          SYMBOL,
          2,
          nchar(SYMBOL)
        )
      )
    )
  ) %>%
  dplyr::select(
    UNIPROT,
    mouse_symbol
  ) %>%
  distinct(
    UNIPROT,
    .keep_all = TRUE
  )

# Create UniProt-to-symbol lookup vector
lookup <- setNames(
  id_map$mouse_symbol,
  id_map$UNIPROT
)

# Apply mapped gene symbols to source and target identifiers
net_mapped <- net_processed %>%
  mutate(
    source_mouse = ifelse(
      source_id %in% names(lookup),
      lookup[source_id],
      source_id
    ),
    target_mouse = ifelse(
      target_id %in% names(lookup),
      lookup[target_id],
      target_id
    )
  )


## Final CollecTRI network cleanup

# Resolve Trp53 naming and collapse duplicate regulator-target interactions
net <- net_mapped %>%
  mutate(
    source_mouse = if_else(
      source_mouse %in% c(
        "Tp53",
        "P04637"
      ),
      "Trp53",
      source_mouse
    ),
    target_mouse = if_else(
      target_mouse %in% c(
        "Tp53",
        "P04637"
      ),
      "Trp53",
      target_mouse
    )
  ) %>%
  dplyr::select(
    source = source_mouse,
    target = target_mouse,
    mor
  ) %>%
  group_by(
    source,
    target
  ) %>%
  summarise(
    mor = mean(
      mor,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  distinct()



## Infer transcription-factor activities using decoupleR
activities <- run_wmean(
  mat = GetAssayData(
    obj,
    assay = "RNA",
    layer = "data"
  ),
  net = net,
  .source = "source",
  .target = "target",
  .mor = "mor",
  times = 100,
  minsize = 5
)

## Convert inferred TF activities to TF-by-cell matrix

# Collapse any duplicate TF/cell entries by their mean inferred activity
activities_wide <- activities %>%
  group_by(
    source,
    condition
  ) %>%
  summarise(
    score = mean(
      score,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols = source,
    names_from = condition,
    values_from = score
  ) %>%
  column_to_rownames(
    "source"
  )

tf_matrix <- as.matrix(
  activities_wide
)

tf_matrix[
  is.na(tf_matrix)
] <- 0

storage.mode(
  tf_matrix
) <- "double"


## Add inferred TF-activity matrix to Seurat object
# TF activity scores are continuous inferred values rather than sequencing counts, and are therefore stored in the assay data layer.

obj[["tfsulm"]] <- CreateAssayObject(
  data = tf_matrix
)

DefaultAssay(obj) <- "tfsulm"

Idents(obj) <- "cell_type_v2"

# Top variable inferred TF activities across K andKP+/- malignant states
n_tfs <- 25

# Extract cell-level inferred TF activity and calculate mean activity within each malignant epithelial state
activity_long <- t(
  as.matrix(
    GetAssayData(
      obj,
      assay = "tfsulm",
      layer = "data"
    )
  )
) %>%
  as.data.frame() %>%
  mutate(
    cell_type = Idents(obj)
  ) %>%
  pivot_longer(
    cols = -cell_type,
    names_to = "source",
    values_to = "score"
  ) %>%
  group_by(
    cell_type,
    source
  ) %>%
  summarise(
    mean_activity = mean(
      score,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Rank TFs by variability in mean inferred activity across tumor states
top_tfs <- activity_long %>%
  group_by(source) %>%
  summarise(
    variability = sd(
      mean_activity,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(variability)
  ) %>%
  slice_head(
    n = n_tfs
  ) %>%
  pull(source)

# Generate tumor-state x TF activity matrix
heatmap_mat <- activity_long %>%
  filter(
    source %in% top_tfs
  ) %>%
  pivot_wider(
    id_cols = cell_type,
    names_from = source,
    values_from = mean_activity
  ) %>%
  column_to_rownames(
    "cell_type"
  ) %>%
  as.matrix()

# Diverging blue-white-red palette
palette_cols <- colorRampPalette(
  rev(
    brewer.pal(
      n = 11,
      name = "RdBu"
    )
  )
)(100)

# Z-score each TF across malignant states for visualization
pheatmap(
  mat = heatmap_mat,
  color = palette_cols,
  border_color = "white",
  scale = "column",
  main = "Top variable inferred TF activity by K - KP+/- tumor states"
)


##### Supp Fig 7A #####
## Myc inferred TF activity and RNA expression

# Myc inferred TF activity
DefaultAssay(obj) <- "tfsulm"

p_myc_activity <- FeaturePlot_scCustom(
  obj,
  features = "Myc",
  reduction = "umap"
) +
  ggtitle(
    "Myc Activity"
  )

p_myc_activity


# Myc expression
DefaultAssay(obj) <- "RNA"

p_myc_expression <- FeaturePlot_scCustom(
  obj,
  features = "Myc",
  reduction = "umap"
) +
  ggtitle(
    "Myc Expression"
  )

p_myc_expression


##### Fig 4E-G #####
## Monocle3 pseudotime analysis of genotype-specific malignant trajectories
# Evolutionary trajectories were inferred independently for K, KP+/-, and published KP-/- malignant epithelial cells using Monocle3

## The complete workflow below is shown for K tumor cells
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(igraph)
library(ggplot2)
library(ggrepel)
library(grid)


## Subset K malignant epithelial cells

DefaultAssay(obj) <- "RNA"
Idents(obj) <- "Genotype_ID"

obj_K <- subset(
  obj,
  idents = "K"
)

obj_K <- JoinLayers(
  obj_K,
  assay = "RNA"
)

# Convert Seurat object to Monocle3 cell_data_set
cds <- as.cell_data_set(
  obj_K
)

cds <- estimate_size_factors(
  cds
)

cds <- cluster_cells(
  cds
)

cds <- learn_graph(
  cds,
  use_partition = FALSE
)

cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(
  obj_K
)

# AT1-like cells were selected interactively as the root population
cds <- order_cells(
  cds
)

# Add final malignant-state annotations
colData(cds)$label <- obj_K$cell_type_v2


## K malignant-cell pseudotime trajectory
p <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  cell_size = 1,
  label_branch_points = FALSE,
  label_leaves = FALSE
) +
  NoAxes()

# Calculate tumor-state label positions
data_df <- data.frame(
  colData(cds),
  reducedDims(cds)$UMAP
)

colnames(data_df)[
  (ncol(data_df) - 1):ncol(data_df)
] <- c(
  "U1",
  "U2"
)

label_pos <- data_df %>%
  group_by(label) %>%
  summarise(
    U1 = median(U1),
    U2 = median(U2),
    .groups = "drop"
  )

final_plot_K <- p +
  geom_text_repel(
    data = label_pos,
    aes(
      x = U1,
      y = U2,
      label = label
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 3.5,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    min.segment.length = 0
  )

final_plot_K

# Directional arrows shown in the final assembled figure were added during figure preparation and were not generated by the Monocle3 analysis code


##### Supp Fig 7B-D #####
## Preparation of K/KP+/- malignant epithelial data for Palantir analysis

# Palantir trajectory analysis was performed in Python/Jupyter on the joint K and KP+/- malignant seurat object.
# The R code below exports the final Seurat object to AnnData format while preserving normalized RNA expression, cell metadata, Harmony embeddings, and UMAP coordinates.

library(reticulate)
library(Seurat)
library(anndata)
library(Matrix)

# Use the Python environment containing Palantir and AnnData dependencies
use_condaenv(
  "palantir_env",
  required = TRUE
)


# Join Seurat v5 RNA layers prior to matrix export
obj_joined <- JoinLayers(
  obj,
  assay = "RNA"
)

# Extract normalized RNA expression matrix
# AnnData requires cells as rows and genes as columns
expression_matrix <- t(
  GetAssayData(
    obj_joined,
    assay = "RNA",
    layer = "data"
  )
)

# Extract cell metadata
cell_metadata <- obj_joined@meta.data

# Preserve dimensional reductions used for downstream Palantir analyses
harmony_embeddings <- Embeddings(
  obj_joined,
  reduction = "harmony"
)

umap_embeddings <- Embeddings(
  obj_joined,
  reduction = "umap"
)

# Construct AnnData object
adata <- AnnData(
  X = expression_matrix,
  obs = cell_metadata,
  obsm = list(
    X_harmony = harmony_embeddings,
    X_umap = umap_embeddings
  )
)


### Palantir trajectory analysis
## The AnnData object generated above was used for downstream Palantir pseudotime, trajectory, and gene-expression trend analyses in Python/Jupyter
## The complete notebook used to generate Supplementary Figure 7B-D is provided separately as: Palantir_Run_K_KPHet.ipynb

##### FIGURE 5 & SUPPLEMENTARY FIGURE 8 #####
### Rank-based K/KP+/- tumor-state signature scoring in KRAS-mutant TCGA-LUAD

## KRAS-mutant TCGA-LUAD tumor samples were identified from GDC somatic mutation calls
# Corresponding GDC STAR-count RNA-seq files were isolated, TPM values were assembled into a gene-by-sample expression matrix, transformed as log2(TPM + 1), and scored using singscore.

library(data.table)
library(tidyverse)
library(fs)
library(singscore)

## Identify KRAS-mutant TCGA-LUAD cases from GDC somatic mutation calls

maf_path <- "cohortMAF.2026-06-22.maf"

mutations <- fread(
  maf_path,
  skip = "Hugo_Symbol"
)

kras_mutants <- mutations %>%
  filter(
    Hugo_Symbol == "KRAS"
  ) %>%
  select(
    Hugo_Symbol,
    Variant_Classification,
    HGVSp_Short,
    Tumor_Sample_Barcode
  ) %>%
  distinct()

# Patient-level TCGA identifiers
kras_patient_barcodes <- unique(
  kras_mutants$Tumor_Sample_Barcode
)

kras_patients_short <- unique(
  substr(
    kras_patient_barcodes,
    1,
    12
  )
)

## Isolate GDC RNA-seq files from KRAS-mutant tumor samples

rna_root_dir <- "GDC_TCGA_LUAD_RNA"
output_dir <- file.path(
  rna_root_dir,
  "rna_counts_files"
)

rna_sample_sheet_path <- "gdc_sample_sheet_RNA.tsv"

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

rna_sample_sheet <- fread(
  rna_sample_sheet_path
)

# Locate GDC STAR-count files
all_rna_files <- list.files(
  path = rna_root_dir,
  pattern = "_star_gene_counts\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

copied_counter <- 0

walk(
  all_rna_files,
  function(file_path) {
    
    folder_uuid <- basename(
      dirname(file_path)
    )
    
    matched_row <- rna_sample_sheet[
      rna_sample_sheet$`File ID` == folder_uuid,
    ]
    
    if (nrow(matched_row) > 0) {
      
      barcode <- matched_row$`Sample ID`[1]
      barcode_short <- substr(barcode, 1, 12)
      vial_code <- substr(barcode, 14, 15)
      
      # Retain KRAS-mutant tumor samples only.
      # TCGA sample type codes 01-09 denote tumor samples.
      if (
        barcode_short %in% kras_patients_short &&
        as.numeric(vial_code) < 10
      ) {
        
        dest_name <- paste0(
          barcode,
          "_star_gene_counts.tsv"
        )
        
        dest_path <- file.path(
          output_dir,
          dest_name
        )
        
        file_copy(
          file_path,
          dest_path,
          overwrite = TRUE
        )
        
        copied_counter <<- copied_counter + 1
      }
    }
  }
)

message(
  "KRAS-mutant tumor RNA files isolated: ",
  copied_counter
)

## Build KRAS-mutant TCGA-LUAD TPM expression matrix

rna_files <- list.files(
  path = output_dir,
  pattern = "\\.tsv$",
  full.names = TRUE
)

master_expression_list <- map(
  rna_files,
  function(file_path) {
    
    # GDC STAR-count files contain four metadata rows before the table
    dt <- fread(
      file_path,
      skip = 4
    )
    
    setDT(dt)
    
    patient_barcode <- str_remove(
      basename(file_path),
      "_star_gene_counts\\.tsv"
    )
    
    # GDC STAR-count columns used in the original analysis:
    # V1 = Ensembl gene ID
    # V2 = gene symbol
    # V7 = TPM unstranded
    selected_data <- dt[
      ,
      .(
        gene_id = V1,
        gene_symbol = V2,
        expression_value = V7
      )
    ]
    
    colnames(selected_data)[3] <- patient_barcode
    
    setkey(
      selected_data,
      gene_id
    )
    
    selected_data
  }
)

## Merge samples into a single gene-by-sample matrix

master_expression_matrix <- purrr::reduce(
  master_expression_list,
  function(df1, df2) {
    
    merge(
      df1,
      df2,
      by = c(
        "gene_id",
        "gene_symbol"
      ),
      all = TRUE
    )
  }
)

# Remove STAR alignment-summary rows
master_expression_matrix <- master_expression_matrix[
  !startsWith(
    gene_id,
    "N_"
  )
]

# Remove entries lacking a gene symbol
master_expression_matrix <- master_expression_matrix[
  gene_symbol != "" &
    !is.na(gene_symbol)
]

# Rename gene-symbol column
setnames(
  master_expression_matrix,
  "gene_symbol",
  "Hybridization REF"
)

patient_cols <- setdiff(
  colnames(master_expression_matrix),
  c(
    "gene_id",
    "Hybridization REF"
  )
)

setcolorder(
  master_expression_matrix,
  c(
    "Hybridization REF",
    patient_cols
  )
)

# Ensembl IDs were used only for merging and were removed afterward
master_expression_matrix[
  ,
  gene_id := NULL
]

## Transform TPM expression
patient_cols <- setdiff(
  colnames(master_expression_matrix),
  "Hybridization REF"
)

master_expression_matrix[
  ,
  (patient_cols) := lapply(
    .SD,
    function(x) log2(x + 1)
  ),
  .SDcols = patient_cols
]

# Collapse duplicate gene symbols by their mean log2(TPM + 1) expression
master_expression_matrix <- master_expression_matrix[
  ,
  lapply(
    .SD,
    mean,
    na.rm = TRUE
  ),
  by = "Hybridization REF"
]

# Save final expression matrix
fwrite(
  master_expression_matrix,
  file = "aggregated_KRAS_tumor_RNAseq_TPM_unstranded_FINAL.txt",
  sep = "\t",
  row.names = FALSE
)


## Convert expression table to numeric matrix for singscore

expr_mat <- as.matrix(
  master_expression_matrix[
    ,
    -1,
    with = FALSE
  ]
)

rownames(expr_mat) <- master_expression_matrix$`Hybridization REF`

storage.mode(expr_mat) <- "numeric"

# Confirm unique gene symbols
sum(
  duplicated(
    rownames(expr_mat)
  )
)


## Human ortholog signatures corresponding to K and KP+/- malignant states

Fras1_ECM_remodeling <- c(
  "HNF1A-AS1", "CYP2C8", "CYP2C9", "CUBN", "UNC13C",
  "ADGRG7", "EYA4", "PRTG", "KLHL3", "TMEFF2", "PCSK5"
)

Ncoa1_AT2_like <- c(
  "C17orf82", "SLC24A4", "ZNF830", "ADAM22",
  "ST6GALNAC3", "COL23A1", "MEG3"
)

Basal_like <- c(
  "EEF1B2", "GAPDH", "RACK1", "MIF", "NME2",
  "UBA52", "PPIA", "NACA", "TMSB10"
)

AT1_like <- c(
  "SH3RF3", "SEMA3A", "GALNT18", "COL4A4", "RNF150",
  "NEBL", "RADIL", "COL4A3", "SHROOM3", "RTKN2"
)

Fetub_Gastric_AT2_like <- c(
  "HP", "HSPA2", "CYP2F1", "APOB", "FZD8", "MYCT1",
  "MUC1", "RASL11A", "SPHK1", "IGFBP6", "CYP1B1",
  "CFI", "CXCL5", "FETUB", "IGFBP5", "MUC5B",
  "ALB", "HSD11B2", "NAPSA", "TF"
)

Transitional_AT2_like <- c(
  "MECOM", "CACNA1C", "PRMT8", "PDE7B", "PPP1R14C",
  "BMP1", "MITF", "ROR1", "ITGA9", "ADAM19"
)

Hypoxic <- c(
  "CLDN4", "MSLN", "PROCR", "R3HDML", "GZMA",
  "GKN2", "CTSE", "LY6G6C", "AREG", "TPM2",
  "GJB4", "TNFRSF23"
)

Cryab_AT1_like <- c(
  "PAKAP", "SPOCK2", "SEC14L3", "IGFBP2", "CRYAB",
  "CAVIN2", "SPARC", "HOPX", "IFIT3", "TPPP3",
  "PMP22", "CYP4B1", "PRDX6", "HIST1H2BC",
  "CAV1", "CLIC3", "SERPINB9"
)

AT2_like <- c(
  "CXCL17", "HLA-DRA", "LCN2", "CD74"
)

PI3K_Akt_Ductal <- c(
  "SRGAP3", "ITGA2", "AMPD3", "PFKP",
  "CDH13", "CYP2S1", "ZNF330"
)

High_Cycling <- c(
  "DEPDC1", "CDC25C", "ANKLE1", "ASPM", "KIF18B",
  "AURKB", "SGOL1", "TOP2A", "KIFC1", "CCNB1",
  "NEIL3", "MKI67"
)

Nfat_Neovascularization <- c(
  "DST", "C15ORF48", "TRERF1", "GLT1D1", "PRDM5",
  "VAV3", "PAM", "PPM1H", "PPP2R3A", "NAV3",
  "SNHG11", "MTUS1", "MTMR3", "SFMBT2", "TFAP2B",
  "GRIP1", "PTP4A3", "PVT1", "SLC45A4"
)

Mdk_AT2_like <- c(
  "A2ML1", "FABP5", "GSTA4", "PCBD1", "NPM1",
  "CDK4", "NHP2", "PRDX2", "CSTB", "KLK8",
  "PLET1", "C3", "C11ORF96", "MARCKSL1", "TSPAN3"
)

all_signatures <- list(
  Fras1_ECM_producing = Fras1_ECM_producing,
  Ncoa1_AT2_like = Ncoa1_AT2_like,
  Basal_like = Basal_like,
  Cryab_AT1_like = Cryab_AT1_like,
  AT1_like = AT1_like,
  Fetub_Gastric_AT2_like = Fetub_Gastric_AT2_like,
  Transitional_AT2_like = Transitional_AT2_like,
  EMT_like = EMT_like,
  AT2_like = AT2_like,
  PI3K_Akt_Ductal = PI3K_Akt_Ductal,
  High_Cycling = High_Cycling,
  Nfat_Neovascularization = Nfat_Neovascularization,
  Mdk_AT2_like = Mdk_AT2_like
)


## Rank-based singscore calculation
ranked_expr <- rankGenes(
  expr_mat
)

singscore_list <- lapply(
  names(all_signatures),
  function(sig) {
    
    genes <- all_signatures[[sig]]
    
    genes_present <- intersect(
      genes,
      rownames(expr_mat)
    )
    
    message(
      sig,
      ": ",
      length(genes_present),
      "/",
      length(genes),
      " genes detected"
    )
    
    if (length(genes_present) < 1) {
      return(
        rep(
          NA,
          ncol(expr_mat)
        )
      )
    }
    
    scored <- simpleScore(
      ranked_expr,
      upSet = genes_present
    )
    
    scored$TotalScore
  }
)

singscore_matrix <- do.call(
  rbind,
  singscore_list
)

rownames(singscore_matrix) <- names(
  all_signatures
)

colnames(singscore_matrix) <- colnames(
  expr_mat
)

score_matrix <- singscore_matrix


## Save signature scores

singscore_dt <- as.data.table(
  singscore_matrix,
  keep.rownames = "Signature"
)

fwrite(
  singscore_dt,
  file = "KRAS_patients_singscore_signature_scores.txt",
  sep = "\t"
)


##### Fig 5A #####
## Distribution of malignant-state singscores across KRAS-mutant TCGA-LUAD

distinct_tumor_palette <- c(
  "Fras1_ECM_producing"     = "#F8A19F",
  "Ncoa1_AT2_like"          = "#0000FF",
  "Basal_like"              = "#FA0087",
  "Cryab_AT1_like"          = "#00A08B",
  "AT1_like"                = "#18FF2D",
  "Fetub_Gastric_AT2_like"  = "#2ED9FF",
  "Transitional_AT2_like"   = "#F6222E",
  "EMT_like"                = "#FFFF00",
  "AT2_like"                = "#FEAF16",
  "PI3K_Akt_Ductal"         = "#AA0DFE",
  "High_Cycling"            = "#A87E63",
  "Nfat_Neovascularization" = "#708090",
  "Mdk_AT2_like"            = "#00008B"
)

long_plot_df <- as.data.frame(
  score_matrix
) %>%
  rownames_to_column(
    "Signature"
  ) %>%
  pivot_longer(
    -Signature,
    names_to = "Sample",
    values_to = "Singscore"
  ) %>%
  mutate(
    Signature = factor(
      Signature,
      levels = names(
        distinct_tumor_palette
      )
    )
  )

ggplot(
  long_plot_df,
  aes(
    x = Signature,
    y = Singscore,
    fill = Signature
  )
) +
  geom_violin(
    alpha = 0.5,
    scale = "width",
    color = "gray30",
    linewidth = 0.4
  ) +
  geom_jitter(
    width = 0.1,
    size = 0.4,
    alpha = 0.15,
    color = "gray40"
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    color = "black",
    alpha = 0.8,
    fill = "white",
    linewidth = 0.15
  ) +
  scale_fill_manual(
    values = distinct_tumor_palette
  ) +
  coord_flip() +
  theme_bw(
    base_size = 12
  ) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(
      face = "bold",
      color = "black"
    ),
    axis.text.x = element_text(
      color = "black"
    ),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "",
    y = "Enrichment Score"
  )


#### Cox proportional-hazards analysis of K and KP+/- tumor-state signatures
### KRAS-mutant TCGA-LUAD singscores were merged with TCGA clinical data
## Cox proportional-hazards models were run for each signature independently
# Hazard ratios represent the effect of a 1-SD increase in signature score
# Models were adjusted for age at diagnosis and stratified by pathological stage group (Early Stage vs Late Stage).

library(TCGAbiolinks)
library(data.table)
library(tidyverse)
library(survival)
library(broom)
library(ggplot2)
library(survminer)

## Download and prepare TCGA-LUAD clinical data

all_luad_clinical <- GDCquery_clinic(
  project = "TCGA-LUAD",
  type = "clinical"
)

setDT(all_luad_clinical)

clinical_extracted <- all_luad_clinical[
  ,
  .(
    patient_id = toupper(submitter_id),
    
    Vital_Status = tolower(vital_status),
    
    OS_event = ifelse(
      tolower(vital_status) == "dead",
      1,
      0
    ),
    
    days_to_death = as.numeric(days_to_death),
    
    days_to_followup = as.numeric(
      days_to_last_follow_up
    ),
    
    Gender = factor(
      sex_at_birth
    ),
    
    Age_at_Diagnosis = as.numeric(
      age_at_index
    ),
    
    Pathologic_Stage = factor(
      ajcc_pathologic_stage
    )
  )]


## Calculate overall-survival time

clinical_extracted[
  ,
  OS_time := ifelse(
    OS_event == 1,
    days_to_death,
    days_to_followup
  )
]

# Exclude patients lacking survival time or with <30 days follow-up
clinical_extracted <- clinical_extracted[
  !is.na(OS_time) &
    OS_time >= 30
]


## Define pathological stage groups

# Remove substage designations to obtain Stage I-IV
clinical_extracted$Main_Stage <- gsub(
  "[ABCD1234]",
  "",
  clinical_extracted$Pathologic_Stage
)

clinical_extracted$Main_Stage <- trimws(
  clinical_extracted$Main_Stage
)

clinical_extracted$Main_Stage <- factor(
  clinical_extracted$Main_Stage,
  levels = c(
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IV"
  )
)

# Remove patients lacking pathological-stage information
clinical_extracted <- clinical_extracted[
  !is.na(Main_Stage)
]

# Group Stage I-II as Early Stage and Stage III-IV as Late Stage
clinical_extracted$Stage_Group <- ifelse(
  clinical_extracted$Main_Stage %in% c(
    "Stage I",
    "Stage II"
  ),
  "Early Stage",
  "Late Stage"
)

clinical_extracted$Stage_Group <- factor(
  clinical_extracted$Stage_Group,
  levels = c(
    "Early Stage",
    "Late Stage"
  )
)

message(
  "Total patients remaining for analysis: ",
  nrow(clinical_extracted)
)

print(
  table(
    clinical_extracted$Stage_Group
  )
)


## Convert singscore matrix to patient-level dataframe

sig_cols <- rownames(
  singscore_matrix
)

score_df <- as.data.frame(
  t(
    singscore_matrix
  )
)

score_df$sample_barcode <- rownames(
  score_df
)

score_df$patient_id <- substr(
  score_df$sample_barcode,
  1,
  12
)

# Prefer TCGA 01A tumor aliquots when duplicate 01A/01B samples exist
score_df <- score_df %>%
  mutate(
    sample_priority = ifelse(
      grepl(
        "-01A",
        sample_barcode
      ),
      1,
      2
    )
  ) %>%
  arrange(
    patient_id,
    sample_priority,
    sample_barcode
  ) %>%
  distinct(
    patient_id,
    .keep_all = TRUE
  ) %>%
  select(
    -sample_priority
  )

# Confirm one tumor sample per patient
sum(
  duplicated(
    score_df$patient_id
  )
)

## Merge clinical and singscore data

master_modeling_df <- merge(
  score_df,
  clinical_extracted,
  by = "patient_id"
)

master_modeling_df <- master_modeling_df %>%
  filter(
    !is.na(OS_time),
    !is.na(OS_event),
    !is.na(Age_at_Diagnosis),
    !is.na(Stage_Group)
  )

message(
  "Final integrated cohort size: ",
  nrow(master_modeling_df)
)

print(
  table(
    master_modeling_df$OS_event
  )
)

print(
  table(
    master_modeling_df$Stage_Group
  )
)

fwrite(
  as.data.table(
    master_modeling_df
  ),
  "KRAS_integrated_clinical_singscore_signatures.txt",
  sep = "\t"
)


### Cox proportional-hazards analysis
## Hazard ratios are calculated per 1-SD increase in signature score.
# Models are adjusted for age and stratified by stage group.

model_df <- master_modeling_df %>%
  filter(
    !is.na(OS_time),
    !is.na(OS_event),
    !is.na(Age_at_Diagnosis),
    !is.na(Stage_Group)
  )

# Z-score continuous signature scores
model_df[
  ,
  sig_cols
] <- scale(
  model_df[
    ,
    sig_cols
  ]
)


## Fit Cox model independently for each tumor-state signature

cox_results <- map_dfr(
  sig_cols,
  function(sig) {
    
    fml <- as.formula(
      paste0(
        "Surv(OS_time, OS_event) ~ ",
        sig,
        " + Age_at_Diagnosis + strata(Stage_Group)"
      )
    )
    
    fit <- coxph(
      fml,
      data = model_df
    )
    
    tidy_fit <- tidy(
      fit,
      exponentiate = TRUE,
      conf.int = TRUE
    ) %>%
      filter(
        term == sig
      )
    
    # Global proportional-hazards test
    ph_test <- tryCatch(
      {
        cox.zph(fit)$table[
          "GLOBAL",
          "p"
        ]
      },
      error = function(e) NA
    )
    
    tibble(
      Signature = sig,
      HR = tidy_fit$estimate,
      CI_low = tidy_fit$conf.low,
      CI_high = tidy_fit$conf.high,
      Cox_p = tidy_fit$p.value,
      PH_p = ph_test
    )
  }
) %>%
  mutate(
    FDR = p.adjust(
      Cox_p,
      method = "BH"
    ),
    
    Cox_label = paste0(
      "HR = ",
      round(HR, 2),
      " [",
      round(CI_low, 2),
      "-",
      round(CI_high, 2),
      "]",
      ", Cox p = ",
      signif(Cox_p, 3)
    )
  ) %>%
  arrange(
    Cox_p
  )

cox_results

## Identify nominally significant signatures

sig_cox <- cox_results %>%
  filter(
    Cox_p < 0.05
  ) %>%
  pull(
    Signature
  )

sig_cox


##### Fig 5B #####
## Forest plot of signature-associated hazard ratios

forest_df <- cox_results %>%
  filter(
    !is.na(HR),
    !is.na(CI_low),
    !is.na(CI_high),
    CI_low > 0,
    CI_high > 0
  ) %>%
  arrange(
    HR
  ) %>%
  mutate(
    Signature = factor(
      Signature,
      levels = Signature
    ),
    
    Sig_status = case_when(
      FDR < 0.05 ~ "FDR < 0.05",
      Cox_p < 0.05 ~ "Nominal p < 0.05",
      TRUE ~ "NS"
    ),
    
    # Visual precision proxy used to scale point size
    CI_Width = CI_high - CI_low,
    Plot_Weight = 1 / CI_Width
  )


forest_plot <- ggplot(
  forest_df,
  aes(
    x = HR,
    y = Signature
  )
) +
  
  # Null hazard ratio
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "gray50"
  ) +
  
  # 95% confidence intervals
  geom_errorbarh(
    aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    height = 0.2,
    linewidth = 0.7,
    color = "black"
  ) +
  
  # Hazard-ratio estimates
  geom_point(
    aes(
      color = Sig_status,
      size = Plot_Weight
    ),
    shape = 15
  ) +
  
  scale_x_log10(
    breaks = c(
      0.25,
      0.33,
      0.50,
      0.75,
      0.90,
      1.00,
      1.10,
      1.25,
      1.50,
      2.00
    ),
    labels = c(
      "0.25",
      "0.33",
      "0.50",
      "0.75",
      "0.90",
      "1.00",
      "1.10",
      "1.25",
      "1.50",
      "2.00"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "FDR < 0.05" = "firebrick",
      "Nominal p < 0.05" = "firebrick",
      "NS" = "black"
    )
  ) +
  
  scale_size_continuous(
    range = c(
      2,
      5.5
    )
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    
    axis.text.y = element_text(
      color = "black"
    ),
    
    legend.position = "right"
  ) +
  
  labs(
    title = "Age-adjusted, stage-stratified Cox analysis",
    subtitle = "HR per 1 SD increase in signature score",
    x = "Hazard ratio",
    y = "",
    color = "Significance",
    size = "Relative precision"
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        size = 3.5
      )
    )
  )

print(
  forest_plot
)

#### Four-group Kaplan-Meier survival analysis
### Patients were divided by pathological stage group (Early vs Late) and by median signature score (High vs Low), generating four survival groups:
  # Early Stage / Low Signature
  #   Early Stage / High Signature
  #   Late Stage / Low Signature
  #   Late Stage / High Signature


library(ggplot2)
library(survminer)
library(survival)
library(dplyr)

## Helper function for four-group log-rank test
get_logrank_p <- function(df, group_col) {
  
  fit <- survdiff(
    as.formula(
      paste0(
        "Surv(OS_time, OS_event) ~ ",
        group_col
      )
    ),
    data = df
  )
  
  pchisq(
    fit$chisq,
    length(fit$n) - 1,
    lower.tail = FALSE
  )
}


## Function to perform and plot four-group Kaplan-Meier analysis
plot_4way_km <- function(
    sig,
    df = model_df,
    outdir = "."
) {
  
  # Retrieve Cox-model result for this signature
  cox_row <- cox_results %>%
    filter(
      Signature == sig
    )
  
  # Prepare stage x signature groups
  tmp <- df %>%
    filter(
      !is.na(OS_time),
      !is.na(OS_event),
      !is.na(Stage_Group),
      !is.na(.data[[sig]])
    ) %>%
    mutate(
      Stage_Group = trimws(
        as.character(Stage_Group)
      ),
      
      # Median split of signature score
      Sig_Group = ifelse(
        .data[[sig]] >= median(
          .data[[sig]],
          na.rm = TRUE
        ),
        "High",
        "Low"
      ),
      
      KM4_Group = paste(
        Stage_Group,
        Sig_Group,
        sep = "_"
      )
    )
  
  # Fix plotting order
  tmp$KM4_Group <- factor(
    tmp$KM4_Group,
    levels = c(
      "Early Stage_Low",
      "Early Stage_High",
      "Late Stage_Low",
      "Late Stage_High"
    )
  )
  
  tmp <- tmp %>%
    filter(
      !is.na(KM4_Group)
    )
  
  
  
  ## Four-way log-rank test
  logrank_p <- get_logrank_p(
    tmp,
    "KM4_Group"
  )
  
  
  
  ## Kaplan-Meier fit
  surv_fit <- survfit(
    Surv(
      OS_time,
      OS_event
    ) ~ KM4_Group,
    data = tmp
  )
  
  
  
  ## Plot
  
  label_text <- paste0(
    sig,
    "\nCox p = ",
    signif(
      cox_row$Cox_p,
      3
    ),
    "; 4-way log-rank p = ",
    signif(
      logrank_p,
      3
    )
  )
  
  p <- ggsurvplot(
    surv_fit,
    data = tmp,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    censor = TRUE,
    title = label_text,
    xlab = "Overall survival time",
    ylab = "Overall survival probability",
    legend.title = "",
    risk.table.height = 0.25,
    ggtheme = theme_classic(
      base_size = 13
    )
  )
  
  # Save Kaplan-Meier curve and risk table
  pdf(
    file.path(
      outdir,
      paste0(
        sig,
        "_4way_KM.pdf"
      )
    ),
    width = 7,
    height = 7
  )
  
  print(p)
  
  dev.off()
  
  return(
    list(
      plot = p,
      logrank_p = logrank_p
    )
  )
}


##### Fig 5C #####
## Mdk+ AT2-like signature
Mdk_KM <- plot_4way_km(
  "Mdk_AT2_like"
)

Mdk_KM$plot
Mdk_KM$logrank_p


##### Supp Fig 8A #####
## AT2-like signature

AT2_KM <- plot_4way_km(
  "AT2_like"
)

AT2_KM$plot
AT2_KM$logrank_p


#### K and KP+/- tumor-state signature expression in the human scRNA NSCLC atlas 
### The processed tumor-cell object from Salcher et al. (2022) was used to evaluate AT2-like and Mdk+ AT2-like signature enrichment across major NSCLC tumor subtypes.

library(Seurat)
library(dplyr)
library(ggplot2)


## Load processed Salcher et al. NSCLC tumor-cell object
NSCLC_Atlas_2022_Tumor <- readRDS(
  "NSCLC_Atlas_2022_Tumor.rds"
)

# The object was previously normalized and tumor cells were subsetted from the global NSCLC atlas
# Gene identifiers were converted from Ensembl IDs to gene symbols before downstream signature analysis

Idents(NSCLC_Atlas_2022_Tumor) <- "cell_type_tumor"

unique(
  NSCLC_Atlas_2022_Tumor$cell_type_tumor
)


## Group malignant-cell annotations by pathological subtype
NSCLC_Atlas_2022_Tumor@meta.data <- NSCLC_Atlas_2022_Tumor@meta.data %>%
  mutate(
    pathology = case_when(
      grepl("LUAD", cell_type_tumor)  ~ "LUAD",
      grepl("LUSC", cell_type_tumor)  ~ "LUSC",
      grepl("NSCLC", cell_type_tumor) ~ "NSCLC-NOS",
      TRUE                            ~ "Other"
    )
  )

# Confirm pathological group assignments
table(
  NSCLC_Atlas_2022_Tumor$pathology
)

## Pathology color palette
pathology_custom_colors <- c(
  "LUAD"      = "#6B4450",
  "LUSC"      = "#A3C6C0",
  "NSCLC-NOS" = "#1D2A44"
)


##### Fig 5D #####

Idents(NSCLC_Atlas_2022_Tumor) <- "pathology"

DimPlot(
  NSCLC_Atlas_2022_Tumor,
  group.by = "pathology",
  cols = pathology_custom_colors,
  label = FALSE,
  pt.size = 0.5
)


VlnPlot(
  NSCLC_Atlas_2022_Tumor,
  features = c(
    "AT2_like_UCell",
    "Mdk_AT2_like_UCell"
  ),
  group.by = "pathology",
  cols = pathology_custom_colors,
  pt.size = 0,
  ncol = 1
)

##### Supp Fig 8B #####
#### K/KP+/- malignant-state signature enrichment across human NSCLC tumor states
### Human ortholog signatures derived from the K and KP+/- malignant epithelial states were scored in the Salcher et al. (2022) NSCLC tumor-cell atlas using UCell. Cell-level signature scores were Z-scored and then averaged within each annotated human NSCLC tumor state for visualization.

library(Seurat)
library(UCell)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(pheatmap)

## Define NSCLC tumor-state color palette
state_pathology_palette <- c(
  "Tumor cells LUAD"         = "#FF1493",
  "Tumor cells LUAD mitotic" = "#800020",
  "Tumor cells LUAD EMT"     = "#FF00FF",
  "Tumor cells LUAD NE"      = "#F497A9",
  "Tumor cells LUAD MSLN"    = "#B22222",
  "Tumor cells LUSC"         = "#A3C6C0",
  "Tumor cells LUSC mitotic" = "#4A736E",
  "Tumor cells NSCLC mixed"  = "#1D2A44"
)


## Visualize annotated tumor states
Idents(NSCLC_Atlas_2022_Tumor) <- "cell_type_tumor"

DimPlot(
  NSCLC_Atlas_2022_Tumor,
  group.by = "cell_type_tumor",
  cols = state_pathology_palette,
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5
)


## Score K/KP+/- malignant-state signatures using UCell

NSCLC_Atlas_2022_Tumor <- AddModuleScore_UCell(
  NSCLC_Atlas_2022_Tumor,
  features = all_signatures
)

# Metadata columns generated by UCell
original_names <- names(all_signatures)

sig_names <- paste0(
  original_names,
  "_UCell"
)

# Map UCell metadata-column names back to clean signature names
clean_labels <- setNames(
  original_names,
  sig_names
)


## Extract cell-level UCell scores

plot_data <- FetchData(
  NSCLC_Atlas_2022_Tumor,
  vars = c(
    sig_names,
    "cell_type_tumor"
  )
)


## Z-score each signature across all tumor cells

plot_data_scaled <- plot_data %>%
  tibble::rownames_to_column(
    var = "Cell_ID"
  ) %>%
  mutate(
    across(
      all_of(sig_names),
      ~ as.vector(scale(.))
    )
  ) %>%
  arrange(
    cell_type_tumor
  ) %>%
  tibble::column_to_rownames(
    var = "Cell_ID"
  )


## Average Z-scored signature activity by annotated NSCLC tumor state
plot_data_avg <- plot_data_scaled %>%
  group_by(
    cell_type_tumor
  ) %>%
  summarise(
    across(
      all_of(sig_names),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  as.data.frame()

rownames(plot_data_avg) <- plot_data_avg$cell_type_tumor

plot_data_avg$cell_type_tumor <- NULL


## Prepare signature-by-tumor-state matrix
heatmap_matrix_avg <- t(
  as.matrix(
    plot_data_avg
  )
)

# Replace UCell metadata names with manuscript signature names
rownames(heatmap_matrix_avg) <- clean_labels[
  rownames(heatmap_matrix_avg)
]


### Heatmap of average signature enrichment across human NSCLC tumor states
pheatmap(
  heatmap_matrix_avg,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = TRUE,
  color = colorRampPalette(
    c(
      "#2166ac",
      "#f7f7f7",
      "#b2182b"
    )
  )(100),
  breaks = seq(
    -1.5,
    1.5,
    length.out = 100
  ),
  angle_col = 45,
  main = "Average signature activity across NSCLC tumor states"
)


####  Preparation of KRAS-mutant TCGA-LUAD copy-number segmentation data
### KRAS-mutant TCGA-LUAD cases were identified from GDC somatic mutation calls
## Corresponding tumor copy-number segmentation files were then isolated using the GDC sample sheet and TCGA sample barcodes.

library(data.table)
library(tidyverse)
library(fs)


## Identify KRAS-mutant TCGA-LUAD cases
maf_path <- "cohortMAF.2026-06-22.maf"

mutations <- fread(
  maf_path,
  skip = "Hugo_Symbol"
)

kras_mutants <- mutations %>%
  filter(
    Hugo_Symbol == "KRAS"
  ) %>%
  select(
    Hugo_Symbol,
    Variant_Classification,
    HGVSp_Short,
    Tumor_Sample_Barcode
  ) %>%
  distinct()

# TCGA barcodes containing KRAS mutations
kras_patient_barcodes <- unique(
  kras_mutants$Tumor_Sample_Barcode
)

# Reduce to patient-level TCGA identifier (TCGA-XX-XXXX)
kras_patients_short <- unique(
  substr(
    kras_patient_barcodes,
    1,
    12
  )
)

message(
  "KRAS-mutant TCGA-LUAD patients identified: ",
  length(kras_patients_short)
)


## Load GDC copy-number segmentation sample sheet
wgs_root_dir <- "GDC_TCGA_LUAD_WGS/WGS"

sample_sheet_path <- file.path(
  wgs_root_dir,
  "gdc_sample_sheet.2026-06-23.tsv"
)

sample_sheet <- fread(
  sample_sheet_path
)


## Locate GDC copy-number segmentation files
all_seg_files <- list.files(
  path = wgs_root_dir,
  pattern = "_grch38\\.seg\\.v2\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

message(
  "Total segmentation files identified: ",
  length(all_seg_files)
)


## Isolate KRAS-mutant tumor segmentation files
segment_folder <- file.path(
  wgs_root_dir,
  "segment_files"
)

if (!dir.exists(segment_folder)) {
  dir.create(
    segment_folder,
    recursive = TRUE
  )
}

copied_counter <- 0

walk(
  all_seg_files,
  function(file_path) {
    
    # GDC File ID corresponds to the parent directory of each downloaded file
    folder_uuid <- basename(
      dirname(file_path)
    )
    
    matched_row <- sample_sheet[
      sample_sheet$`File ID` == folder_uuid,
    ]
    
    if (nrow(matched_row) > 0) {
      
      barcode <- matched_row$`Sample ID`[1]
      
      # Patient-level TCGA identifier
      barcode_short <- substr(
        barcode,
        1,
        12
      )
      
      # TCGA sample-type code
      vial_code <- substr(
        barcode,
        14,
        15
      )
      
      # Retain KRAS-mutant tumor samples only.
      # TCGA sample-type codes 01-09 correspond to tumor samples.
      if (
        barcode_short %in% kras_patients_short &&
        as.numeric(vial_code) < 10
      ) {
        
        dest_name <- paste0(
          barcode,
          "_grch38.seg.v2.txt"
        )
        
        dest_path <- file.path(
          segment_folder,
          dest_name
        )
        
        file_copy(
          file_path,
          dest_path,
          overwrite = TRUE
        )
        
        copied_counter <<- copied_counter + 1
      }
    }
  }
)

message(
  "KRAS-mutant tumor segmentation files isolated: ",
  copied_counter
)


## Merge KRAS-mutant tumor segmentation files
mutant_files <- list.files(
  path = segment_folder,
  pattern = "_grch38\\.seg\\.v2\\.txt$",
  full.names = TRUE
)

message(
  "Merging ",
  length(mutant_files),
  " KRAS-mutant tumor segmentation files."
)

master_wgs_matrix <- map_dfr(
  mutant_files,
  function(file_path) {
    
    dt <- fread(
      file_path
    )
    
    # Recover human-readable TCGA sample barcode from filename
    sample_barcode <- str_remove(
      basename(file_path),
      "_grch38\\.seg\\.v2\\.txt$"
    )
    
    # Remove GDC internal aliquot identifier if present
    if ("GDC_Aliquot" %in% colnames(dt)) {
      dt[
        ,
        GDC_Aliquot := NULL
      ]
    }
    
    # Add TCGA sample barcode as the first column
    clean_dt <- data.table(
      Sample = sample_barcode,
      dt
    )
    
    clean_dt
  }
)

fwrite(
  master_wgs_matrix,
  file = output_master_file,
  sep = "\t",
  row.names = FALSE
)

print(
  head(
    master_wgs_matrix
  )
)


#### Genome-wide CNA frequency profiles
### Human CNA frequencies were calculated from KRAS-mutant TCGA-LUAD tumor copy-number segmentation data
## Mouse CNA frequencies were calculated from inferCNV profiles of KP+/- malignant cells.

## load packages
library(CNTools)
library(GenomicRanges)
library(dplyr)
library(ggplot2)


### Human KRAS-mutant TCGA-LUAD CNA frequency profile
seg <- master_wgs_matrix


## Human chromosome lengths and cumulative genomic coordinates
chr_lengths <- c(
  "1" = 248956422,
  "2" = 242193529,
  "3" = 198295559,
  "4" = 190214555,
  "5" = 181538259,
  "6" = 170805979,
  "7" = 159345973,
  "8" = 145138636,
  "9" = 138394717,
  "10" = 133797422,
  "11" = 135086622,
  "12" = 133275309,
  "13" = 114364328,
  "14" = 107043718,
  "15" = 101991189,
  "16" = 90338345,
  "17" = 83257441,
  "18" = 80373285,
  "19" = 58617616,
  "20" = 64444167,
  "21" = 46709983,
  "22" = 50818468
)

chr_df <- data.frame(
  chr = names(chr_lengths),
  chr_length = as.numeric(chr_lengths)
)

chr_df$cum_start <- c(
  0,
  cumsum(chr_df$chr_length)[-nrow(chr_df)]
)

chr_df$cum_end <- chr_df$cum_start + chr_df$chr_length


## Initial 1-Mb bin-based human CNA frequency calculation
seg$bin <- floor(
  seg$Start / 1e6
)

bins <- tileGenome(
  seqlengths = chr_lengths,
  tilewidth = 1e6,
  cut.last.tile.in.chrom = TRUE
)

seg_gr <- GRanges(
  seqnames = as.character(seg$Chromosome),
  ranges = IRanges(
    seg$Start,
    seg$End
  ),
  score = seg$Segment_Mean,
  sample = seg$Sample
)

hits <- findOverlaps(
  bins,
  seg_gr
)

ov <- data.frame(
  bin = queryHits(hits),
  sample = mcols(seg_gr)$sample[
    subjectHits(hits)
  ],
  score = mcols(seg_gr)$score[
    subjectHits(hits)
  ]
)

gain_threshold <- 0.2
loss_threshold <- -0.2

gain_freq <- ov %>%
  group_by(bin) %>%
  summarise(
    freq_gain = mean(
      score > gain_threshold
    )
  )

loss_freq <- ov %>%
  group_by(bin) %>%
  summarise(
    freq_loss = mean(
      score < loss_threshold
    )
  )

freq_df <- full_join(
  gain_freq,
  loss_freq,
  by = "bin"
)

freq_df$freq_gain[
  is.na(freq_df$freq_gain)
] <- 0

freq_df$freq_loss[
  is.na(freq_df$freq_loss)
] <- 0


## Initial human CNA frequency visualization
ggplot(
  freq_df,
  aes(x = bin)
) +
  geom_area(
    aes(y = freq_gain),
    fill = "firebrick",
    alpha = 0.8
  ) +
  geom_area(
    aes(y = -freq_loss),
    fill = "steelblue",
    alpha = 0.8
  ) +
  theme_classic() +
  labs(
    x = "Genome Position (1 Mb bins)",
    y = "Fraction of Samples",
    title = "TCGA CNV Frequency Plot"
  )


## Human chromosome-concatenated CNA frequency calculation
human_bins <- seg %>%
  group_by(
    Chromosome,
    bin
  ) %>%
  summarise(
    gain = mean(
      Segment_Mean > 0.2
    ),
    loss = mean(
      Segment_Mean < -0.2
    ),
    .groups = "drop"
  )

human_bins <- merge(
  human_bins,
  chr_df,
  by.x = "Chromosome",
  by.y = "chr"
)

human_bins$genome_pos <- (
  human_bins$cum_start +
    human_bins$bin * 1e6
)

##### Supp Fig 8C #####
## Final human CNA frequency plot

ggplot(
  human_bins,
  aes(genome_pos)
) +
  geom_area(
    aes(y = gain),
    fill = "firebrick",
    alpha = 0.8
  ) +
  geom_area(
    aes(y = -loss),
    fill = "steelblue",
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = chr_df$cum_start,
    colour = "grey70",
    linewidth = 0.3
  ) +
  scale_x_continuous(
    breaks = (
      chr_df$cum_start +
        chr_df$chr_length / 2
    ),
    labels = chr_df$chr
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  coord_cartesian(
    ylim = c(-1, 1)
  ) +
  theme_classic() +
  labs(
    x = "Human Chromosome",
    y = "Fraction of Samples",
    title = "TCGA LUAD CNV Frequency (GISTIC-style)"
  )


### Mouse KP+/- inferCNV frequency profile
Mousegenmatrix <- read.table(
  "KPHetinfercnvobs.txt",
  sep = "\t",
  header = TRUE,
  row.names = 1
)

Mousegenedf <- read.table(
  "Geneorderforinfercnv.txt"
)

colnames(Mousegenedf) <- c(
  "gene",
  "chr",
  "start",
  "end"
)


## Calculate gene-level mouse gain/loss frequencies
gain_threshold <- 1.1
loss_threshold <- 0.9

freq_df <- data.frame(
  gene = rownames(Mousegenmatrix),
  
  gain_freq = rowMeans(
    Mousegenmatrix > gain_threshold
  ),
  
  loss_freq = rowMeans(
    Mousegenmatrix < loss_threshold
  )
)

freq_df <- merge(
  freq_df,
  Mousegenedf,
  by = "gene"
)

head(freq_df)

freq_df$bin <- floor(
  freq_df$start / 1e6
)

mouse_bins <- freq_df %>%
  group_by(
    chr,
    bin
  ) %>%
  summarise(
    gain = mean(gain_freq),
    loss = mean(loss_freq),
    .groups = "drop"
  )


## Preliminary chromosome-track visualization
ggplot(
  mouse_bins,
  aes(x = bin)
) +
  geom_area(
    aes(y = gain),
    fill = "red"
  ) +
  geom_area(
    aes(y = -loss),
    fill = "blue"
  ) +
  facet_grid(
    chr ~ .,
    scales = "free_x",
    space = "free_x"
  ) +
  theme_classic() +
  labs(
    x = "Genomic position",
    y = "Frequency"
  )


## Mouse chromosome lengths and cumulative genomic coordinates
mouse_chr_lengths <- c(
  chr1  = 195471971,
  chr2  = 182113224,
  chr3  = 160039680,
  chr4  = 156508116,
  chr5  = 151834684,
  chr6  = 149736546,
  chr7  = 145441459,
  chr8  = 129401213,
  chr9  = 124595110,
  chr10 = 130694993,
  chr11 = 122082543,
  chr12 = 120129022,
  chr13 = 120421639,
  chr14 = 124902244,
  chr15 = 104043685,
  chr16 = 98207768,
  chr17 = 94987271,
  chr18 = 90702639,
  chr19 = 61431566
)

chr_df <- data.frame(
  chr = names(mouse_chr_lengths),
  chr_length = as.numeric(mouse_chr_lengths)
)

chr_df$cum_start <- c(
  0,
  cumsum(chr_df$chr_length)[-nrow(chr_df)]
)

chr_df$cum_end <- (
  chr_df$cum_start +
    chr_df$chr_length
)

chr_df


## Convert mouse bins to concatenated genome coordinates
mouse_bins <- merge(
  mouse_bins,
  chr_df,
  by = "chr"
)

mouse_bins$genome_pos <- (
  mouse_bins$cum_start +
    mouse_bins$bin * 1e6
)

chr_breaks <- chr_df$cum_start
chr_labels <- chr_df$chr


##### Supp Fig 8D #####
ggplot(
  mouse_bins,
  aes(genome_pos)
) +
  geom_area(
    aes(y = gain),
    fill = "firebrick",
    alpha = 0.8
  ) +
  geom_area(
    aes(y = -loss),
    fill = "steelblue",
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = chr_df$cum_start,
    colour = "grey70",
    linewidth = 0.3
  ) +
  scale_x_continuous(
    breaks = (
      chr_df$cum_start +
        chr_df$chr_length / 2
    ),
    labels = chr_df$chr
  ) +
  theme_classic() +
  labs(
    x = "Mouse Chromosome",
    y = "Frequency",
    title = "Mouse inferCNV Frequency Plot"
  ) +
  coord_cartesian(
    ylim = c(-1, 1)
  ) +
  geom_hline(
    yintercept = 0,
    colour = "black",
    linewidth = 0.4
  )



#### Overlay of human KRAS-mutant LUAD and mouse KP+/- CNA frequency profiles
## Human and mouse CNA-frequency summaries generated above were labeled by species, normalized to a common maximum CNA frequency, combined into a
#  shared table, and visualized on their chromosome-concatenated genomic coordinate scales

## Annotate species and calculate relative genomic positions
human_bins$species <- "Human"
mouse_bins$species <- "Mouse"

human_bins$scaled_pos <- human_bins$genome_pos / max(human_bins$genome_pos)
mouse_bins$scaled_pos <- mouse_bins$genome_pos / max(mouse_bins$genome_pos)



## Normalize CNA frequencies across human and mouse datasets
global_max <- max(
  c(
    human_bins$gain,
    human_bins$loss,
    mouse_bins$gain,
    mouse_bins$loss
  ),
  na.rm = TRUE
)

human_bins$gain <- human_bins$gain / global_max
human_bins$loss <- human_bins$loss / global_max

mouse_bins$gain <- mouse_bins$gain / global_max
mouse_bins$loss <- mouse_bins$loss / global_max



## Prepare chromosome-boundary coordinates and harmonize table structure
human_breaks <- chr_df$cum_start
mouse_breaks <- mouse_bins$cum_start

colnames(mouse_bins)[
  colnames(mouse_bins) == "chr"
] <- "Chromosome"

all_bins <- rbind(
  human_bins,
  mouse_bins
)



## Convert combined CNA-frequency table to long format
all_long <- all_bins |>
  pivot_longer(
    cols = c(
      gain,
      loss
    ),
    names_to = "event",
    values_to = "freq"
  )



## Define chromosome-boundary positions for comparative plotting
human_vlines <- chr_df$cum_start
mouse_vlines <- mouse_bins$cum_start

## Generate event labels and plotting-color annotations

all_long$event1 <- paste(
  all_long$species,
  all_long$event,
  sep = "_"
)

all_long$color <- with(
  all_long,
  ifelse(
    species == "Human" & event == "gain",
    "red",
    ifelse(
      species == "Human" & event == "loss",
      "blue",
      ifelse(
        species == "Mouse" & event == "gain",
        "darkred",
        "darkblue"
      )
    )
  )
)


##### Fig 5E #####
## Comparative human-mouse CNA frequency plot

ggplot(
  all_long,
  aes(
    x = genome_pos,
    y = freq
  )
) +
  
  # Gains
  geom_area(
    data = subset(
      all_long,
      event == "gain"
    ),
    aes(
      fill = interaction(
        species,
        event
      )
    ),
    alpha = 0.5
  ) +
  
  # Losses
  geom_area(
    data = subset(
      all_long,
      event == "loss"
    ),
    aes(
      y = -freq,
      fill = interaction(
        species,
        event
      )
    ),
    alpha = 0.5
  ) +
  
  # Human chromosome boundaries
  geom_vline(
    xintercept = human_vlines,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.3
  ) +
  
  # Mouse chromosome boundaries
  geom_vline(
    xintercept = mouse_vlines,
    linetype = "dotted",
    color = "grey80",
    linewidth = 0.3
  ) +
  
  # Manual colors
  scale_fill_manual(
    values = c(
      "Human.gain" = "red",
      "Human.loss" = "blue",
      "Mouse.gain" = "green",
      "Mouse.loss" = "yellow"
    )
  ) +
  
  coord_cartesian(
    ylim = c(
      -1,
      1
    )
  ) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  
  theme_classic() +
  
  labs(
    x = "Genome Position (Human + Mouse concatenated scale)",
    y = "Normalized CNV Frequency",
    title = "Human vs Mouse CNV Frequency (GISTIC-style comparative plot)"
  )



#### Cross-species synteny analysis of human and mouse CNA landscapes
### Human KRAS-mutant TCGA-LUAD copy-number segmentation data and mouse KP+/- inferCNV profiles were mapped to gene-level genomic coordinates
## One-to-one human-mouse orthologs were identified using biomaRt and used to visualize cross-species correspondence of selected genomic loci

# load packages
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(GenomicRanges)
library(dplyr)
library(biomaRt)
library(circlize)
library(ComplexHeatmap)


## Human gene coordinates
genes_human <- genes(
  EnsDb.Hsapiens.v86,
  return.type = "GRanges"
)

genes_human <- data.frame(
  gene = genes_human$gene_name,
  chr = as.character(seqnames(genes_human)),
  start = start(genes_human),
  end = end(genes_human)
)

genes_human <- genes_human %>%
  group_by(gene) %>%
  summarise(
    chr = first(chr),
    start = min(start),
    end = max(end),
    .groups = "drop"
  )


## Map human CNA segmentation data to genes
seg_gr <- GRanges(
  seqnames = seg$Chromosome,
  ranges = IRanges(
    seg$Start,
    seg$End
  ),
  score = seg$Segment_Mean,
  sample = seg$Sample
)

gene_gr <- GRanges(
  seqnames = genes_human$chr,
  ranges = IRanges(
    genes_human$start,
    genes_human$end
  ),
  gene = genes_human$gene
)

hits <- findOverlaps(
  gene_gr,
  seg_gr
)

human_gene_cnv <- data.frame(
  gene = mcols(gene_gr)$gene[queryHits(hits)],
  sample = mcols(seg_gr)$sample[subjectHits(hits)],
  cnv = mcols(seg_gr)$score[subjectHits(hits)]
)


## Human gene-level CNA gain/loss frequencies
human_freq <- human_gene_cnv %>%
  group_by(
    gene,
    sample
  ) %>%
  summarise(
    cnv = mean(cnv),
    .groups = "drop"
  ) %>%
  mutate(
    gain = cnv > 0.2,
    loss = cnv < -0.2
  ) %>%
  group_by(gene) %>%
  summarise(
    gain_freq = mean(gain),
    loss_freq = mean(loss),
    .groups = "drop"
  )


## Mouse gene coordinates
genes <- genes(
  EnsDb.Mmusculus.v79,
  return.type = "GRanges"
)

mouse_gene_coords <- data.frame(
  gene = genes$gene_name,
  chr = as.character(seqnames(genes)),
  start = start(genes),
  end = end(genes)
)

mouse_gene_coords <- mouse_gene_coords %>%
  group_by(gene) %>%
  summarise(
    chr = first(chr),
    start = min(start),
    end = max(end),
    .groups = "drop"
  )


## Mouse gene-level inferCNV gain/loss frequencies
mouse_gene_cnv <- data.frame(
  gene = rownames(Mousegenmatrix),
  gain_freq = rowMeans(
    Mousegenmatrix > gain_threshold,
    na.rm = TRUE
  ),
  loss_freq = rowMeans(
    Mousegenmatrix < loss_threshold,
    na.rm = TRUE
  )
)


## Retrieve one-to-one human-mouse orthologs
human <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

mouse <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl"
)

orthologs <- getBM(
  attributes = c(
    "hgnc_symbol",
    "mmusculus_homolog_associated_gene_name",
    "mmusculus_homolog_orthology_type"
  ),
  mart = human
)

orthologs <- orthologs %>%
  dplyr::filter(
    mmusculus_homolog_associated_gene_name != "",
    mmusculus_homolog_orthology_type == "ortholog_one2one"
  ) %>%
  dplyr::select(
    human_gene = hgnc_symbol,
    mouse_gene = mmusculus_homolog_associated_gene_name
  )

colnames(orthologs) <- c(
  "human_gene",
  "mouse_gene"
)


## Merge orthologs with mouse CNA information
mouse_syn <- orthologs %>%
  merge(
    mouse_gene_cnv,
    by.x = "mouse_gene",
    by.y = "gene"
  ) %>%
  merge(
    mouse_gene_coords,
    by.x = "mouse_gene",
    by.y = "gene"
  )

# Retain canonical mouse chromosomes for inspection
mouse_gene_coords1 <- mouse_gene_coords %>%
  dplyr::filter(
    chr %in% c(
      as.character(1:19),
      "X",
      "Y"
    )
  )


## Preliminary mouse circos visualization
circos.clear()

mouse_chr_lengths <- tapply(
  mouse_gene_coords$end,
  mouse_gene_coords$chr,
  max
)

mouse_chr_df <- data.frame(
  chr = names(mouse_chr_lengths),
  length = as.numeric(mouse_chr_lengths)
)

circos.initialize(
  factors = mouse_chr_df$chr,
  xlim = cbind(
    0,
    mouse_chr_df$length
  )
)

circos.trackPlotRegion(
  ylim = c(-1, 1),
  panel.fun = function(x, y) {
    
    chr <- CELL_META$sector.index
    
    df <- mouse_syn %>%
      dplyr::filter(
        chr == !!chr
      )
    
    circos.points(
      df$start,
      df$gain_freq,
      col = "red",
      pch = 16,
      cex = 0.3
    )
    
    circos.points(
      df$start,
      -df$loss_freq,
      col = "blue",
      pch = 16,
      cex = 0.3
    )
  },
  track.height = 0.2
)


## Merge orthologs with human CNA information
human_syn <- orthologs %>%
  merge(
    human_freq,
    by.x = "human_gene",
    by.y = "gene"
  ) %>%
  merge(
    genes_human,
    by.x = "human_gene",
    by.y = "gene"
  )


## Summarize human and mouse ortholog coordinates and CNA frequencies
human_syn1 <- human_syn %>%
  dplyr::rename(
    gene = 1
  ) %>%
  dplyr::mutate(
    gene = as.character(gene)
  ) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    chr.human = first(chr),
    pos.human = mean(start),
    gain_freq.human = mean(gain_freq),
    loss_freq.human = mean(loss_freq),
    .groups = "drop"
  )

mouse_syn1 <- mouse_syn %>%
  dplyr::rename(
    gene = 1
  ) %>%
  dplyr::mutate(
    gene = as.character(gene)
  ) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    chr.mouse = first(chr),
    pos.mouse = mean(start),
    gain_freq.mouse = mean(gain_freq),
    loss_freq.mouse = mean(loss_freq),
    .groups = "drop"
  )

orthologs1 <- orthologs %>%
  dplyr::rename(
    human_gene = 1,
    mouse_gene = 2
  )

syn_h <- orthologs1 %>%
  dplyr::inner_join(
    human_syn1,
    by = c(
      "human_gene" = "gene"
    )
  )

Finalsyn <- syn_h %>%
  dplyr::inner_join(
    mouse_syn1,
    by = c(
      "mouse_gene" = "gene"
    )
  )


## Preliminary combined human-mouse circos visualization
circos.clear()

human_chr <- unique(
  Finalsyn$chr.human
)

mouse_chr <- unique(
  Finalsyn$chr.mouse
)

chr_df <- Finalsyn %>%
  dplyr::group_by(
    chr.mouse
  ) %>%
  dplyr::summarise(
    chr_len = max(pos.mouse),
    .groups = "drop"
  )

circos.initialize(
  factors = chr_df$chr.mouse,
  xlim = cbind(
    0,
    chr_df$chr_len
  )
)

# Chromosome labels
circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.05,
  bg.border = NA,
  panel.fun = function(x, y) {
    
    chr <- CELL_META$sector.index
    
    circos.text(
      x = CELL_META$xcenter,
      y = 0.5,
      labels = chr,
      facing = "bending.inside",
      cex = 0.8,
      col = "black"
    )
  }
)

# Mouse CNA track
circos.trackPlotRegion(
  ylim = c(-1, 1),
  track.height = 0.15,
  panel.fun = function(x, y) {
    
    chr <- CELL_META$sector.index
    
    df <- Finalsyn %>%
      dplyr::filter(
        chr.mouse == chr
      )
    
    circos.points(
      df$pos.mouse,
      df$gain_freq.mouse,
      col = "red",
      pch = 16,
      cex = 0.3
    )
    
    circos.points(
      df$pos.mouse,
      -df$loss_freq.mouse,
      col = "blue",
      pch = 16,
      cex = 0.3
    )
  }
)

# Human CNA track
circos.trackPlotRegion(
  ylim = c(-1, 1),
  track.height = 0.15,
  panel.fun = function(x, y) {
    
    chr <- CELL_META$sector.index
    
    df <- Finalsyn %>%
      dplyr::filter(
        chr.human == chr
      )
    
    circos.points(
      df$pos.human,
      df$gain_freq.human,
      col = "darkred",
      pch = 16,
      cex = 0.3
    )
    
    circos.points(
      df$pos.human,
      -df$loss_freq.human,
      col = "darkblue",
      pch = 16,
      cex = 0.3
    )
  }
)

Finalsyn$lwd <- 1 + 4 * abs(
  Finalsyn$gain_freq.human -
    Finalsyn$gain_freq.mouse
)

for (i in seq_len(nrow(Finalsyn))) {
  
  circos.link(
    sector.index1 = Finalsyn$chr.human[i],
    point1 = Finalsyn$pos.human[i],
    sector.index2 = Finalsyn$chr.mouse[i],
    point2 = Finalsyn$pos.mouse[i],
    col = rgb(
      0,
      0,
      0,
      0.25
    ),
    lwd = Finalsyn$lwd[i],
    border = NA
  )
}

## Retain valid shared genomic coordinates
valid_chr <- get.all.sector.index()

syn <- Finalsyn %>%
  dplyr::filter(
    chr.human %in% valid_chr,
    chr.mouse %in% valid_chr
  )

syn <- syn %>%
  dplyr::filter(
    pos.human > 0,
    pos.mouse > 0,
    pos.human < 2.5e8,
    pos.mouse < 2.0e8
  )

syn$lwd <- 1 + 4 * abs(
  syn$gain_freq.human -
    syn$gain_freq.mouse
)


## Build separate human and mouse chromosome sectors

# Prefix human and mouse chromosome names to prevent duplicate sector IDs
syn1 <- syn %>%
  mutate(
    chr.human = paste0(
      "H_",
      chr.human
    ),
    chr.mouse = paste0(
      "M_",
      chr.mouse
    )
  )

human_chr <- sort(
  unique(
    syn1$chr.human
  )
)

mouse_chr <- sort(
  unique(
    syn1$chr.mouse
  )
)

circos_order <- c(
  human_chr,
  mouse_chr
)

## Prepare human and mouse CNA-bin tracks for circos plotting

colnames(human_bins)[
  colnames(human_bins) == "Chromosome"
] <- "chr"

colnames(mouse_bins)[
  colnames(mouse_bins) == "Chromosome"
] <- "chr"

mouse_bins$chr <- gsub(
  "chr",
  "",
  mouse_bins$chr
)

human_bins1 <- human_bins %>%
  mutate(
    chr = paste0(
      "H_",
      chr
    )
  )

mouse_bins1 <- mouse_bins %>%
  mutate(
    chr = paste0(
      "M_",
      chr
    )
  )


## Determine chromosome lengths and plotting order

human_chr <- sort(
  unique(
    human_bins1$chr
  )
)

mouse_chr <- sort(
  unique(
    mouse_bins1$chr
  )
)

human_len <- human_bins1 %>%
  group_by(chr) %>%
  summarise(
    chr_len = max(
      end,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

mouse_len <- mouse_bins1 %>%
  group_by(chr) %>%
  summarise(
    chr_len = max(
      end,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

circos_order <- c(
  human_chr,
  mouse_chr
)

chr_lengths <- bind_rows(
  human_len,
  mouse_len
)

chr_lengths <- chr_lengths[
  match(
    circos_order,
    chr_lengths$chr
  ),
]

# Sanity checks
stopifnot(
  all(
    chr_lengths$chr == circos_order
  )
)

stopifnot(
  anyDuplicated(
    chr_lengths$chr
  ) == 0
)


## Convert CNA frequency bins to genomic peaks
bin_size <- 5e6

human_peak <- human_bins1 %>%
  group_by(
    chr,
    start
  ) %>%
  summarise(
    cnv_score = mean(
      gain - loss,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    end = start + bin_size
  ) %>%
  dplyr::filter(
    end > start
  )

mouse_peak <- mouse_bins1 %>%
  group_by(
    chr,
    start
  ) %>%
  summarise(
    cnv_score = mean(
      gain - loss,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    end = start + bin_size
  ) %>%
  dplyr::filter(
    end > start
  )


## Scale CNA profiles
scale_range <- function(x) {
  
  (
    x - min(
      x,
      na.rm = TRUE
    )
  ) /
    (
      max(
        x,
        na.rm = TRUE
      ) -
        min(
          x,
          na.rm = TRUE
        )
    )
}

human_peak$cnv_scaled <- scale_range(
  human_peak$cnv_score
)

mouse_peak$cnv_scaled <- scale_range(
  mouse_peak$cnv_score
)


## Final circos-layout sanity checks
stopifnot(
  anyDuplicated(
    chr_lengths$chr
  ) == 0,
  anyDuplicated(
    circos_order
  ) == 0,
  all(
    circos_order %in%
      chr_lengths$chr
  )
)

chr_lengths <- chr_lengths[
  match(
    circos_order,
    chr_lengths$chr
  ),
]


## Initialize final human-mouse circos layout

circos.clear()

n <- length(
  chr_lengths$chr
)

gap.after <- rep(
  1,
  n
)

# Large gap separating human and mouse chromosome sectors
gap.after[
  length(human_chr)
] <- 20

circos.par(
  cell.padding = c(
    0,
    0,
    0,
    0
  ),
  track.margin = c(
    0.01,
    0.01
  ),
  start.degree = 90,
  gap.after = gap.after
)

circos.initialize(
  factors = chr_lengths$chr,
  xlim = cbind(
    0,
    chr_lengths$chr_len
  )
)


## Chromosome labels
circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.08,
  bg.border = NA,
  panel.fun = function(x, y) {
    
    chr <- CELL_META$sector.index
    
    circos.text(
      x = CELL_META$xcenter,
      y = 0.5,
      labels = chr,
      facing = "bending.inside",
      cex = 0.6,
      niceFacing = TRUE
    )
  }
)


## CNA gain/loss tracks

# Human
circos.genomicTrack(
  human_bins1[
    ,
    c(
      "chr",
      "start",
      "end",
      "gain",
      "loss"
    )
  ],
  panel.fun = function(region, value, ...) {
    
    col <- ifelse(
      value[[1]] > value[[2]],
      "#d73027",
      "#4575b4"
    )
    
    circos.genomicRect(
      region,
      value,
      col = col,
      border = NA
    )
  },
  track.height = 0.08
)

# Mouse
circos.genomicTrack(
  mouse_bins1[
    ,
    c(
      "chr",
      "start",
      "end",
      "gain",
      "loss"
    )
  ],
  panel.fun = function(region, value, ...) {
    
    col <- ifelse(
      value[[1]] > value[[2]],
      "#d73027",
      "#4575b4"
    )
    
    circos.genomicRect(
      region,
      value,
      col = col,
      border = NA
    )
  },
  track.height = 0.08
)


## CNA frequency-bar tracks
circos.genomicTrack(
  human_peak[
    ,
    c(
      "chr",
      "start",
      "end",
      "cnv_score"
    )
  ],
  panel.fun = function(region, value, ...) {
    
    circos.genomicLines(
      region,
      value,
      type = "h",
      col = ifelse(
        value[[1]] > 0,
        "#d73027",
        "#4575b4"
      ),
      lwd = 2
    )
  },
  track.height = 0.15
)

circos.genomicTrack(
  mouse_peak[
    ,
    c(
      "chr",
      "start",
      "end",
      "cnv_score"
    )
  ],
  panel.fun = function(region, value, ...) {
    
    circos.genomicLines(
      region,
      value,
      type = "h",
      col = ifelse(
        value[[1]] > 0,
        "#d73027",
        "#4575b4"
      ),
      lwd = 2
    )
  },
  track.height = 0.15
)


##### Supp Fig 8E #####
## Cross-species CNA correspondence at canonical KRAS-associated loci

driver_genes <- c(
  "TP53",
  "KRAS",
  "STK11",
  "KEAP1"
)

syn_driver1 <- syn1 %>%
  dplyr::filter(
    human_gene %in% driver_genes
  )


## Plot selected human-mouse KRAS-associated syntenic pair

for (i in seq_len(nrow(syn_driver1))) {
  
  circos.link(
    syn_driver1$chr.human[i],
    syn_driver1$pos.human[i],
    syn_driver1$chr.mouse[i],
    syn_driver1$pos.mouse[i],
    col = adjustcolor(
      "gold",
      alpha.f = 0.6
    ),
    lwd = 2,
    border = NA
  )
}


## Add human and mouse driver-gene labels

circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.06,
  bg.border = NA,
  panel.fun = function(x, y) {
    
    sector <- CELL_META$sector.index
    
    # Human labels
    if (
      sector %in%
      syn_driver1$chr.human
    ) {
      
      idx <- which(
        syn_driver1$chr.human ==
          sector
      )
      
      circos.text(
        x = syn_driver1$pos.human[idx],
        y = 0.5,
        labels = syn_driver1$human_gene[idx],
        cex = 0.65,
        col = "black",
        facing = "inside",
        niceFacing = TRUE
      )
    }
    
    # Mouse labels
    if (
      sector %in%
      syn_driver1$chr.mouse
    ) {
      
      idx <- which(
        syn_driver1$chr.mouse ==
          sector
      )
      
      circos.text(
        x = syn_driver1$pos.mouse[idx],
        y = 0.5,
        labels = syn_driver1$mouse_gene[idx],
        cex = 0.55,
        col = "grey30",
        facing = "inside",
        niceFacing = TRUE
      )
    }
  }
)

##### Fig 5F #####
## Cross-species CNA correspondence at KP+/- malignant-state marker loci

kp_marker_genes <- c(
  "ST6GALNAC3",
  "CUBN",
  "AKT3",
  "PRMT8",
  "NCOA1",
  "FRAS1",
  "PDE7B",
  "EYA4",
  "MYC",
  "CDKN2A"
)

syn_kp_markers <- syn1 %>%
  dplyr::filter(
    human_gene %in%
      kp_marker_genes
  )


## Plot selected KP+/- marker-gene syntenic pairs
for (i in seq_len(nrow(syn_kp_markers))) {
  
  circos.link(
    syn_kp_markers$chr.human[i],
    syn_kp_markers$pos.human[i],
    syn_kp_markers$chr.mouse[i],
    syn_kp_markers$pos.mouse[i],
    col = adjustcolor(
      "gold",
      alpha.f = 0.6
    ),
    lwd = 2,
    border = NA
  )
}


## Add human and mouse KP+/- marker-gene labels

circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.06,
  bg.border = NA,
  panel.fun = function(x, y) {
    
    sector <- CELL_META$sector.index
    
    # Human labels
    if (
      sector %in%
      syn_kp_markers$chr.human
    ) {
      
      idx <- which(
        syn_kp_markers$chr.human ==
          sector
      )
      
      circos.text(
        x = syn_kp_markers$pos.human[idx],
        y = 0.5,
        labels = syn_kp_markers$human_gene[idx],
        cex = 0.65,
        col = "black",
        facing = "inside",
        niceFacing = TRUE
      )
    }
    
    # Mouse ortholog labels
    if (
      sector %in%
      syn_kp_markers$chr.mouse
    ) {
      
      idx <- which(
        syn_kp_markers$chr.mouse ==
          sector
      )
      
      circos.text(
        x = syn_kp_markers$pos.mouse[idx],
        y = 0.5,
        labels = syn_kp_markers$mouse_gene[idx],
        cex = 0.55,
        col = "grey30",
        facing = "inside",
        niceFacing = TRUE
      )
    }
  }
)


### Manual final figure formatting
## Ribbon colors in the final manuscript panels were manually adjusted during figure preparation after inspection of the corresponding human and mouse CNA profiles

