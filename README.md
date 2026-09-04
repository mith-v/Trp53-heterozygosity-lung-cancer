# Trp53-heterozygosity-lung-cancer

This repository contains the R and Python code used for computational analyses and figure generation in the accompanying manuscript.

## Scope of submitted code

No standalone software package or custom mathematical algorithm was developed for this study. The submitted R and Python code comprises analysis scripts that implement and adapt established computational methods and publicly available software packages for the datasets analyzed in the manuscript.

## Repository contents

Code is organized according to the corresponding main and supplementary figure panels. Where the same workflow was applied independently to multiple genotypes or datasets, a representative complete workflow is provided.

### `Source_Code.R`

Primary R script containing analyses for:

- single-cell RNA-seq preprocessing, integration, clustering, and annotation;
- differential expression and pathway enrichment;
- tumor-state annotation, signature generation, and UCell scoring;
- reciprocal label transfer between K, KP+/- and published KP tumor states;
- inferCNV copy-number inference and CNV burden analyses;
- CytoTRACE differentiation/plasticity analysis;
- transcription-factor activity inference using decoupleR/CollecTRI;
- Monocle3 pseudotime analysis;
- preparation of data for Palantir analysis;
- TCGA-LUAD transcriptomic, clinical/survival, mutation, and copy-number analyses;
- analysis of a published human NSCLC single-cell atlas; and
- comparative human–mouse CNA and ortholog/synteny analyses.

Figure-specific section headers identify the corresponding manuscript panels.

### `Palantir_Run_K_KPHet.ipynb`

Jupyter notebook containing the Palantir analyses underlying Supplementary Figure 7B-D. The R code used to generate the AnnData input object is included in `Source_Code.R`.

## Software

Major software/packages used include:

**R:** Seurat 5.1.0, Harmony 1.2.1, UCell 2.10.1, fgsea 1.32.0, msigdbr 7.5.1, inferCNV 1.22.0, CytoTRACE, decoupleR 2.12.0, and Monocle3 1.4.27.

**Python:** Palantir, Scanpy, AnnData, pandas, NumPy, matplotlib, and MAGIC.

Additional package dependencies are loaded within the corresponding sections of `Source_Code.R`.

## Input data

The analyses use study-generated data and publicly available datasets described in the manuscript Methods and Data Availability statement.

Public resources include TCGA-LUAD data from the NCI Genomic Data Commons, Tabula Muris, published mouse KP lung cancer and human NSCLC single-cell RNA-seq datasets, LungMAP, MSigDB, and the CollecTRI regulatory network.

## Running the analyses

Open `Source_Code.R` in R/RStudio and navigate to the relevant figure-specific section. Local file paths reflect the directory structure used for the original analyses and should be replaced as appropriate.

Open `Palantir_Run_K_KPHet.ipynb` in Jupyter Notebook/JupyterLab for the Palantir analyses.

Some analyses require intermediate RDS/Seurat objects or reference datasets generated or processed in preceding sections of the workflow.

## Reproducibility

The code reflects the analytical workflows used to generate the results reported in the manuscript. Some final graphical formatting and annotations were performed during figure assembly and are noted in the source code where applicable; these do not alter the underlying quantitative analyses.

Additional analysis parameters, statistical tests, and dataset information are provided in the manuscript Methods and figure legends.

