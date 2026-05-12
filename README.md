# Transcriptomic Analysis of MYC Target Genes as a Proxy for MYC Transcriptional Activity

Supplementary code for the book chapter:

> Arnal Segura, M. *Transcriptomic Analysis of MYC Target Genes as a Proxy for MYC Transcriptional Activity*. In: **The MYC Gene: Methods and Protocols**, Third Edition.

## Contents

| File | Description |
| --- | --- |
| `MYC_pipeline_MDAMB231.Rmd` | Annotated R Markdown notebook reproducing the bulk RNA-seq branch of the chapter workflow on the public GEO series GSE309250 (MDA-MB-231 cells, 30 µM Omomyc vs. vehicle, 120 h). Renders to a self-contained HTML report. |
| `MYC_pipeline_MDAMB231.R` | Equivalent plain R script with the same code and section comments. |

Both files implement the same pipeline; choose one or the other depending on whether an HTML report is desired.

## Pipeline overview

1. Import of the count matrix from the GEO supplementary FTP.
2. Filtering of low-count genes and TMM normalization with `edgeR`.
3. Differential expression analysis with `limma` + `voom`.
4. MYC transcriptional activity inferred from three complementary approaches:
   * Gene Set Enrichment Analysis (GSEA) with `clusterProfiler`, using MYC-related gene sets from the MSigDB Hallmark and C2 collections retrieved with `msigdbr`.
   * Gene Set Variation Analysis (GSVA) with the `GSVA` package on the log<sub>2</sub>-TMM normalized matrix.
   * Transcription-factor activity inference with the Univariate Linear Model (ULM) implemented in `decoupleR`, using the CollecTRI MYC regulon.

## Software requirements

* R ≥ 4.4.2
* Bioconductor packages: `edgeR`, `limma`, `clusterProfiler`, `enrichplot`, `GSEABase`, `GSVA`, `SummarizedExperiment`, `ComplexHeatmap`, `decoupleR`
* CRAN packages: `data.table`, `dplyr`, `msigdbr` (v7.5.1; if v9.0.0 or later is used the `category` argument must be renamed to `collection`), `fastcluster`

A reproducible record of the package versions used to test the script is printed by the final `sessionInfo()` call.

## How to run

From an interactive R session:

```r
rmarkdown::render("MYC_pipeline_MDAMB231.Rmd")
# or
source("MYC_pipeline_MDAMB231.R")
```

From a terminal:

```bash
Rscript -e 'rmarkdown::render("MYC_pipeline_MDAMB231.Rmd")'
Rscript MYC_pipeline_MDAMB231.R
```

## Data

The example uses the public GEO series **GSE309250**. The count matrix is downloaded on the fly from:

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309250/suppl/GSE309250_counts.txt.gz
```

No locally stored data file is required.

## Contact

M. Arnal Segura — Peptomyc S.L., Barcelona, Spain.
