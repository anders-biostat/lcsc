
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lcsc (linked charts for single cells)

<!-- badges: start -->

<!-- badges: end -->

This package allows you to classify single cells based on
nearest-neighbor smoothing, not rely on unsupervised learning methods.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("anders-biostat/lc-sc")
```

## Example

We start of with the standard workflow consisting of quality control,
normalization, and linear/non-linear dimensional reduction using
`Seurat`.

``` r
library(dplyr)
library(Seurat)

pbmc.data <- Read10X(data.dir = "inst/extdata/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:30)
```

Extract the necessary data from the `Seurat` object.

``` r
# From the Seurat object
counts <- GetAssayData(pbmc, "counts")
meta_data <- pbmc[[]]
pc_space <- Embeddings(pbmc, "pca")
embedding <- Embeddings(pbmc, "umap")

# How is the sample column named?
s = "orig.ident"
cells = tibble::tibble(
  id = rownames(meta_data),
  sample = meta_data[[s]]
)
```

Start the `linked charts` application after generating a nearest
neighborhood graph per sample.

``` r
library(lcsc)

# Build nearest neighborhood graph per sample using nn
k = 50
nn <- run_nn(cells, pc_space, k=k, dim=30)

lc_vis(cells, counts, pc_space, embedding, nn, k=k)
```

The application would look like this after selecting macrophages using
the smoothed expression of CD68.

![](man/figures/pbmc_lcsc_application.png)
