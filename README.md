# Breast-Cancer-Gene-Expression-Biomarker-Analysis
This meta-analysis project aimed to identify potential biomarkers for breast cancer diagnosis and prognosis through a comprehensive gene expression profiling analysis. The study utilized 468 grade 1 and 3 breast cancer samples from four microarray breast cancer GEO datasets, all of which were on the same platform, technology, and biological specimens. Integration was achieved using REM meta-analysis with the GeneMeta Bioconductor R package. Common differentially expressed genes (DEGs) were identified among all four datasets, and meta-combined DEGs were obtained and subjected to functional analysis.

[My website (Gene Expression Profiling)](https://mohammadrezamohajeri.github.io/Breast-Cancer-Gene-Expression-Biomarker-Analysis/index.html)


# Header 1
## Header 2
### Header 3
#### Header 4
##### Header 5
###### Header 6

* Item 1
* Item 2
* Item 3

- Another item 1
- Another item 2
- Another item 3

`counts <- data.frame(read_tsv("countMatrix.tsv"))`

1. First item
2. Second item
3. Third item

For inline code, wrap the code in backticks (`). 


```
 require(readr)
  require(tidyr)
  require(gridExtra)
  require(reshape2)
  require(viridis)
  require(ggplot2)
  require(DESeq2)
  require(biomaRt)
  require(genefilter)
  require(org.Hs.eg.db)
  require(ComplexHeatmap)
  require(clusterProfiler)
  require(readr)
  require(knitr)
# ===========================================
# Importing_Counts_Matrix
# ===========================================
# Set the random seed for reproducible results
set.seed(1234)
# Reading and preparing files:
counts <- data.frame(read_tsv("countMatrix.tsv"))
# Set the row names of the 'counts' data frame to be the first column (usually gene IDs or similar identifiers)
rownames(counts) <- counts[,1]
``` 






