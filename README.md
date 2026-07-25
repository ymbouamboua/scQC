
<!-- README.md is generated from README.Rmd. Please edit README.Rmd only. -->

<div align="center">

<img src="man/figures/scQC_logo.svg" width="300" alt="scQC logo"/>

</div>

<!-- badges: start -->

![R](https://img.shields.io/badge/R-%3E%3D4.3-blue) [![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<!-- badges: end -->

## Overview

**scQC** is an R package for reproducible quality control of single-cell
RNA-seq datasets.

It accepts 10x Genomics outputs, matrices, and Seurat objects, computes
standard QC metrics, filters low-quality cells, optionally detects
doublets, and returns a QC-filtered Seurat object together with
reproducible summaries.

## Features

- Read 10x Genomics directories and HDF5 files.
- Accept Seurat objects, dense matrices, and sparse matrices.
- Calculate mitochondrial and ribosomal percentages.
- Calculate genes-per-UMI and dropout metrics.
- Apply MAD-based, fixed-threshold, or no filtering.
- Optionally detect doublets with `scDblFinder`.
- Process multiple samples.
- Export QC summaries and reproducible outputs.
- Return either a Seurat object or a complete QC result object.

## Installation

You can install the development version from GitHub with:

``` r
# install.packages("pak")
pak::pak("ymbouamboua/scQC")
```

Alternatively, install the package from a local clone:

``` r
devtools::install()
```

## Load scQC

``` r
library(scQC)
```

## Supported inputs

| Input type                      | Supported |
|:--------------------------------|:---------:|
| 10x Genomics directory          |    Yes    |
| `filtered_feature_bc_matrix.h5` |    Yes    |
| `raw_feature_bc_matrix.h5`      |    Yes    |
| `dgCMatrix`                     |    Yes    |
| Base R matrix                   |    Yes    |
| Seurat object                   |    Yes    |
| List of supported inputs        |    Yes    |

## Quick start

### From a 10x Genomics directory

``` r
library(scQC)

result <- run_qc(
  x = "path/to/filtered_feature_bc_matrix",
  sample_id = "Sample1",
  species = "human"
)
```

Inspect the filtered Seurat object:

``` r
result$object
```

Inspect the QC summary:

``` r
result$summary
```

### From a Seurat object

``` r
result <- run_qc(
  x = seurat_object,
  sample_id = "Sample1",
  species = "human"
)
```

### From a count matrix

``` r
result <- run_qc(
  x = counts,
  sample_id = "Sample1",
  species = "human"
)
```

The matrix must contain genes in rows and cells in columns.

## Reading input separately

Use `read_sc_input()` when you only want to convert an input into a
Seurat object.

``` r
object <- read_sc_input(
  x = "path/to/filtered_feature_bc_matrix",
  sample_id = "Sample1"
)
```

The function supports:

- Seurat objects;
- dense count matrices;
- sparse count matrices;
- 10x Genomics directories;
- 10x Genomics HDF5 files.

For multimodal 10x data, the gene-expression feature type is selected
automatically when available.

## Quality-control workflow

``` text
Input
  |
  v
Read and validate input
  |
  v
Calculate QC metrics
  |
  v
Determine filtering thresholds
  |
  v
Remove low-quality cells
  |
  v
Optional doublet detection
  |
  v
Write summaries and outputs
  |
  v
Return QC result
```

## Main QC metrics

`scQC` can calculate the following cell-level metrics:

- total UMI count;
- number of detected genes;
- mitochondrial transcript percentage;
- ribosomal transcript percentage;
- genes-per-UMI ratio;
- log10 genes-per-UMI;
- dropout fraction;
- doublet score and classification when `scDblFinder` is enabled.

## Filtering methods

### MAD-based filtering

MAD-based filtering calculates data-driven thresholds from the
distribution of QC metrics.

``` r
result <- run_qc(
  x = seurat_object,
  sample_id = "Sample1",
  species = "human",
  method = "mad",
  mad_n = 5
)
```

### Fixed-threshold filtering

Fixed filtering uses explicitly defined thresholds.

``` r
result <- run_qc(
  x = seurat_object,
  sample_id = "Sample1",
  species = "human",
  method = "fixed",
  min_feat = 200,
  min_umi = 500,
  max_mito = 10,
  max_ribo = 50
)
```

### No cell filtering

Use `method = "none"` to calculate metrics without removing cells based
on QC thresholds.

``` r
result <- run_qc(
  x = seurat_object,
  sample_id = "Sample1",
  species = "human",
  method = "none",
  rm_dbl = FALSE
)
```

## Doublet detection

Doublet detection can be enabled with:

``` r
result <- run_qc(
  x = seurat_object,
  sample_id = "Sample1",
  species = "human",
  rm_dbl = TRUE
)
```

`scQC` uses `scDblFinder` when doublet detection is requested.

Doublet detection should generally be performed independently for each
biological sample or capture.

## Multiple samples

A named list can be supplied to `run_qc()`.

``` r
inputs <- list(
  Sample1 = "data/Sample1/filtered_feature_bc_matrix",
  Sample2 = "data/Sample2/filtered_feature_bc_matrix",
  Sample3 = sample3_seurat
)

result <- run_qc(
  x = inputs,
  species = "human"
)
```

Naming the list is recommended because the names are used as sample
identifiers when possible.

## Returning only the Seurat object

By default, `run_qc()` returns a complete result object.

To return only the filtered Seurat object:

``` r
qc_object <- run_qc(
  x = seurat_object,
  sample_id = "Sample1",
  species = "human",
  return = "object"
)
```

## Result structure

A complete result may contain:

``` r
names(result)
```

Typical components include:

- `object`: filtered Seurat object;
- `objects`: sample-level Seurat objects for multi-sample runs;
- `summary`: sample-level QC summary;
- `thresholds`: thresholds applied to each metric;
- `doublets`: doublet-detection summary;
- `removed_cells`: cells removed during QC;
- `failures`: inputs that could not be processed;
- `parameters`: parameters used for the run;
- `files`: paths to generated output files.

## Example QC run

``` r
result <- run_qc(
  x = "data/Sample1/filtered_feature_bc_matrix",
  sample_id = "Sample1",
  species = "human",
  method = "mad",
  min_feat = 200,
  min_umi = 500,
  mad_n = 5,
  max_mito = 10,
  calc_ribo = TRUE,
  max_ribo = 50,
  calc_drop = TRUE,
  rm_dbl = TRUE,
  save_object = TRUE,
  outdir = "results/qc/Sample1",
  return = "result",
  verbose = TRUE
)
```

## Output files

Depending on the selected options, `scQC` can write:

- filtered Seurat objects;
- QC summary tables;
- threshold tables;
- removed-cell tables;
- doublet summaries;
- run parameters;
- processing logs or failure summaries.

Open the function documentation with:

``` r
?run_qc
?read_sc_input
```

## Development status

`scQC` is currently under active development.

The interface, defaults, and returned result structure may change before
the first stable release.

## Contributing

Issues, bug reports, and pull requests are welcome.

When reporting a problem, please include:

- the input type;
- the command used;
- the complete error message;
- the output of `sessionInfo()`.

## Citation

A formal citation will be added in a future release.

Until then, please cite the GitHub repository and the package version
used in your analysis.

## License

This project is licensed under the MIT License.
