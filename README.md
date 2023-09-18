# Breast Cancer Gene Expression Biomarker Analysis

## Table of Contents

- [Introduction](#introduction)
- [Phase 1: Differential Expression Analysis of Breast Cancer Dataset GSE25055 using limma in R: Comparing Molecular Subtypes and Grading System](#phase-1-differential-expression-analysis-of-breast-cancer-dataset-gse25055-using-limma-in-r-comparing-molecular-subtypes-and-grading-system)
- [Phase 2: Differential Expression Analysis for Biomarker Discovery in Breast Cancer Grade 1 vs 3](#phase-2-differential-expression-analysis-for-biomarker-discovery-in-breast-cancer-grade-1-vs-3)
- [Phase 3: Integration of Gene Expression Data through Meta-Analysis for Robust Biomarker Discovery in Breast Cancer Grade 1 vs 3](#phase-3-integration-of-gene-expression-data-through-meta-analysis-for-robust-biomarker-discovery-in-breast-cancer-grade-1-vs-3)

- [Research Question](#research-question)
- [Research Hypothesis](#research-hypothesis)
- [Code Pages](#code-pages)
- [Dependencies](#dependencies)
- [How to Run the Code](#how-to-run-the-code)

------------------------------------------------------------------------------------------------------
## Introduction

This project seeks to identify potential biomarkers for breast cancer by conducting differential gene expression analysis across multiple datasets and employing the Meta-Analysis of Random-Effects Models. The goal is to uncover robust gene signatures that could be pivotal for early diagnosis, treatment, and management of the disease. For a more detailed background and objectives, please refer to [Project Summary](https://mohammadrezamohajeri.github.io/Breast-Cancer-Gene-Expression-Biomarker-Analysis/Pages/Project_Summary.html).

## Phase 1: Differential Expression Analysis of Breast Cancer Dataset (GSE25055) using limma in R: Comparing Molecular Subtypes and Grading System

* **For details, please refer to the**: [Phase 1 description](https://mohammadrezamohajeri.github.io/Breast-Cancer-Gene-Expression-Biomarker-Analysis/Pages/About_The_Project1.html)
  

## Phase 2: Differential Expression Analysis for Biomarker Discovery in Breast Cancer (Grade 1 vs 3)

* **For details, please refer to the** [Phase 2 description](https://mohammadrezamohajeri.github.io/Breast-Cancer-Gene-Expression-Biomarker-Analysis/Pages/About_The_Project2.html).

## Phase 3: Integration of Gene Expression Data through Meta-Analysis for Robust Biomarker Discovery in Breast Cancer (Grade 1 vs 3)

* **For details, please refer to the** [Phase 3 description](https://mohammadrezamohajeri.github.io/Breast-Cancer-Gene-Expression-Biomarker-Analysis/Pages/About_The_Project3.html).

------------------------------------------------------------------------------------------------------
## Research Question

### Title: Evolution of the Research Question

For details about the research question, please visit [Research Question](https://mohammadrezamohajeri.github.io/Breast-Cancer-Gene-Expression-Biomarker-Analysis/Pages/Research_Question.html).

------------------------------------------------------------------------------------------------------
## Research Hypothesis

### Title: Research Hypothesis Derived from the Research Question in This Study

For details about the research hypothesis, please visit [Research Hypothesis](https://mohammadrezamohajeri.github.io/Breast-Cancer-Gene-Expression-Biomarker-Analysis/Pages/Research_Hypothesis.html).

------------------------------------------------------------------------------------------------------
## Code Pages

For the source code and notebooks related to each phase of this project, please refer to the following pages:

1. [Code for Phase 1](https://mohammadrezamohajeri.github.io/Breast-Cancer-Gene-Expression-Biomarker-Analysis/Pages/R_Code_1Dataset_Code_Page.html) - *1-Dataset Analysis - limma*
2. [Code for Phase 2](https://mohammadrezamohajeri.github.io/Breast-Cancer-Gene-Expression-Biomarker-Analysis/Pages/R_Code_4Datasets_FC_FDR_Code_Page.html) - *Differential Expression Analysis for Biomarker Discovery in Breast Cancer (G1/G3)*
3. [Code for Phase 3](https://mohammadrezamohajeri.github.io/Breast-Cancer-Gene-Expression-Biomarker-Analysis/Pages/R_Code_4Datasets_Meta_Analysis_Code_Page.html) - *Integration of Gene Expression Data through Meta-Analysis for Robust Biomarker Discovery in Breast Cancer (G1/G3) (4 datasets)*

------------------------------------------------------------------------------------------------------
## Dependencies

### General Requirements
* R programming language
* RStudio (optional, but recommended for easier R script development)

### Phase 1: Differential Expression Analysis of Breast Cancer Dataset (GSE25055) using limma in R: Comparing Molecular Subtypes and Grading System
To run the Phase 1 analysis, the following R packages are required:

1. `require(GEOquery)` - For downloading and processing GEO datasets
2. `require(limma)` - For performing differential expression analysis
3. `require(tidyverse)` - For data manipulation and visualization
4. `require(plotly)` - For interactive plots

```
# R code to install Phase 1 packages
install.packages(c("GEOquery", "limma", "tidyverse", "plotly"))
```

### Phase 2: Differential Expression Analysis for Biomarker Discovery in Breast Cancer (Grade 1 vs 3)
To run the Phase 2 analysis, the following R packages are required:

1. `require(GEOquery)` - For accessing GEO data
2. `require(limma)` - For differential expression analysis
3. `require(tidyverse)` - For data manipulation and visualization
4. `require(plotly)` - For interactive plots
5. `require(ggvenn)` - For creating Venn diagrams using ggplot2

```
# R code to install Phase 2 packages
install.packages(c("GEOquery", "limma", "tidyverse", "plotly", "ggvenn"))
``` 
### Phase 3: Integration of Gene Expression Data through Meta-Analysis for Robust Biomarker Discovery in Breast Cancer (Grade 1 vs 3)
To run the Phase 3 analysis, the following R packages are required:

1. `require(tidyverse)` - For comprehensive data manipulation and visualization, including dplyr, ggplot2, and ...
2. `require(GEOquery)` - For accessing and retrieving gene expression data from the Gene Expression Omnibus (GEO) database
3. `require(reshape2)` - For reshaping data from wide to long format and vice versa, useful for data preprocessing
4. `require(caret)` - For ML and predictive modeling, provides a unified interface for algorithms and performance evaluation
5. `require(GeneMeta)` - For meta-analysis of gene expression data, allowing combining results from multiple studies
6. `require(limma)` - For analyzing microarray data using linear models, including differential expression analysis
7. `require(rgl)` - For creating interactive 3D visualizations, useful for exploring complex data structures
8. `require(sva)` -For surrogate variable analysis, to identify & adjust for unwanted variation in high-dimensional datasets
9. `require(plotly)` - For creating interactive and dynamic plots, including scatter plots, bar plots, and more
10. `require(ggvenn)` - For creating Venn diagrams using ggplot2, useful for visualizing set relationships
11. `require(ggpubr)` - For creating publication-ready plots using ggplot2, providing additional customization options
12. `require(ComplexHeatmap)` - For creating complex heatmaps, allowing vis. of multiple layers of data and annotations
13. `require(RColorBrewer)` - For providing color palettes suitable for data visualization
14. `require(msigdbr)` - For accessing and using the Molecular Signatures Database (MSigDB) for gene set enrichment analysis
15. `require(fgsea)` - For performing gene set enrichment analysis (GSEA) using fast gene set testing algorithms
16. `require(locfit)` - For local regression modeling, useful for fitting flexible curves to data
17. `require(vsn)` - For variance stabilization and normalization, specifically for microarray data analysis

```
# R code to install Phase 3 packages
install.packages(c("tidyverse", "GEOquery", "reshape2", "caret", "GeneMeta", "limma", "rgl", "sva", "plotly", "ggvenn", "ggpubr", "ComplexHeatmap", "RColorBrewer", "msigdbr", "fgsea", "locfit", "vsn"))
``` 
------------------------------------------------------------------------------------------------------
## How to Run the Code

Instructions for running the analysis, including any parameters that should be set.
