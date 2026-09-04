# Trp53-heterozygosity-lung-cancer

This repository contains the custom R and Python code used for computational analyses and figure generation in the accompanying manuscript.

The analysis code is organized according to the corresponding main and supplementary figure panels. Where the same analytical workflow was applied independently to multiple genotypes or datasets, a representative complete workflow is provided and the corresponding applications are indicated in the code.

---

## Repository contents

### `Source_Code.R`

Primary R script containing the computational analyses used throughout the manuscript.

The script is organized by manuscript figure and includes code for:

- single-cell RNA-seq preprocessing, quality filtering, dimensionality reduction, clustering, integration, and annotation;
- differential gene-expression and pathway-enrichment analyses;
- malignant tumor-state annotation and signature generation;
- UCell tumor-state signature scoring;
- reciprocal label-transfer analyses between K,KP+/- and published KP tumor states;
- inferCNV copy-number inference and copy-number burden analyses;
- CytoTRACE analysis of malignant-cell differentiation potential/plasticity;
- transcription-factor activity inference using decoupleR and the CollecTRI regulatory network;
- Monocle3 pseudotime analysis;
- preparation of the joint K,KP+/- malignant-cell dataset for downstream Palantir analysis;
- analysis of TCGA-LUAD clinical, transcriptomic, somatic mutation, and copy-number data;
- rank-based tumor-signature scoring and survival analyses;
- analysis of a published human NSCLC single-cell RNA-seq atlas; and
- comparative human–mouse CNA and ortholog/synteny analyses.

Section headers within the script identify the associated manuscript figure panels.

### `Palantir_Run_K_KPHet.ipynb`

Python/Jupyter notebook containing the Palantir trajectory analyses performed on the joint K and KP+/- malignant epithelial dataset.

The notebook contains the code used for:

- Supplementary Figure 7B, 7C and 7D
  
The corresponding R code used to prepare the AnnData object from the final Seurat object is included in `Source_Code.R`.

---

## Software and package versions

The following software/package versions are reported for the corresponding analyses in the manuscript Methods:

### Primary R / single-cell analysis

- Cell Ranger 7.0.1
- Seurat 5.1.0
- Harmony 1.2.1
- UCell 2.10.1
- fgsea 1.32.0
- msigdbr 7.5.1
- inferCNV 1.22.0
- CytoTRACE 0.3.1
- decoupleR 2.12.0
- Monocle3 1.4.27

Additional R packages used by `Source_Code.R` include:

`SeuratWrappers`, `MAST`, `SCpubr`, `scCustomize`, `ComplexHeatmap`,
`circlize`, `pheatmap`, `RColorBrewer`, `ggplot2`, `ggpubr`, `ggrepel`,
`ggalluvial`, `dittoSeq`, `TCGAbiolinks`, `singscore`, `survival`,
`survminer`, `GenomicRanges`, `biomaRt`, `EnsDb.Hsapiens.v86`,
`EnsDb.Mmusculus.v79`, `data.table`, `dplyr`, `tidyr`, `readr`,
`purrr`, and associated dependencies.

### Python / Palantir analysis

The Palantir analyses were performed in a dedicated Python environment using the packages required by the accompanying notebook, including:

- `palantir`
- `scanpy`
- `anndata`
- `pandas`
- `numpy`
- `matplotlib`
- `magic-impute`

Exact Python-environment package versions can be obtained from the analysis environment if required.

---

## Input data

The analyses use both study-generated and publicly available datasets.

### Study-generated data

Study-generated single-cell RNA-sequencing data are described in the manuscript and its Data Availability statement.

**Accession information is intentionally not reproduced in this peer-review repository.**

### Public datasets and resources

Publicly available datasets/resources used in the analyses include:

- TCGA Lung Adenocarcinoma (TCGA-LUAD) clinical, transcriptomic, somatic mutation, and copy-number data obtained through the NCI Genomic Data Commons;
- Tabula Muris lung single-cell data used to obtain age-matched normal T-cell references for inferCNV;
- a published mouse KP lung cancer single-cell RNA-seq dataset;
- the LungMAP normal mouse lung epithelial reference atlas;
- a published human NSCLC single-cell RNA-seq atlas;
- MSigDB gene-set collections accessed using `msigdbr`; and
- the CollecTRI transcription-factor regulatory network.

For the transcription-factor activity analysis, the publicly available CollecTRI interaction table was downloaded from OmniPath and read locally as `Collectri_raw.txt`, as documented in `Source_Code.R`.

---

## Running the analyses

### R analyses

Open `Source_Code.R` in R or RStudio.

The script is divided into sections corresponding to individual main and supplementary figure panels. The RStudio Document Outline can be used to navigate among figure-specific sections.

Local file paths in the script reflect the directory structure used for the original analyses and should be replaced with the corresponding paths on the user's system.

Several analyses use intermediate Seurat/RDS objects or publicly available reference datasets generated or processed during preceding sections of the workflow. Comments within the script indicate the relevant input objects and representative workflows.

Where a common analysis was applied independently to multiple genotypes or datasets, the complete workflow is shown for a representative dataset and the equivalent applications are indicated in the corresponding figure section.

### Palantir analyses

Open `Palantir_Run_K_KPHet.ipynb` in Jupyter Notebook or JupyterLab.

The notebook contains the Python code used for the Palantir analyses underlying Supplementary Figure 7B-D and may also be viewed directly on GitHub without execution.

---

## Reproducibility notes

The code in this repository reflects the analytical workflows used to generate the results reported in the manuscript.

Some final graphical formatting and annotation steps were performed during manuscript figure assembly and are explicitly noted in the source code where applicable.

Examples include selected graphical annotations or manually formatted elements that do not alter the underlying quantitative analyses.

The manuscript Methods and figure legends provide additional information on statistical tests, thresholds, data-processing parameters, and publicly available datasets used in each analysis.

---

## Data availability during peer review

This repository is provided for confidential peer review of the computational analyses associated with the manuscript.

Study-generated sequencing data and their accession information are described in the manuscript submission materials and will be made available in accordance with journal data-sharing requirements.

---

## Contact

Questions regarding the code or analyses should be directed to the corresponding author(s) of the manuscript.
