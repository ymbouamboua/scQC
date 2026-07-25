#' Read single-cell count data
#'
#' Converts supported single-cell inputs into a Seurat object.
#'
#' Supported inputs include:
#'
#' - a Seurat object;
#' - a feature-by-cell count matrix;
#' - a 10x Genomics matrix directory;
#' - a 10x Genomics HDF5 file;
#' - a Cell Ranger `outs` directory.
#'
#' @param x A Seurat object, count matrix, 10x directory, HDF5 file,
#'   or Cell Ranger output directory.
#' @param sample_id Sample identifier used when constructing a Seurat object.
#' @param assay Name of the assay to create.
#' @param min_cells Minimum number of cells in which a feature must be detected.
#' @param min_features Minimum number of detected features required when
#'   initially constructing the object. For QC workflows, zero is recommended
#'   so that filtering occurs later in `run_qc()`.
#' @param feature_type Feature type to select when a 10x input contains
#'   multiple modalities. Defaults to `"Gene Expression"`.
#' @param gene_column Column of `features.tsv` used as feature names.
#' @param strip_suffix Remove the terminal barcode suffix when reading
#'   matrix-format 10x data.
#' @param verbose Display input detection messages.
#'
#' @return A Seurat object containing raw counts.
#'
#' @export
read_sc_input <- function(
    x,
    sample_id = NULL,
    assay = "RNA",
    min_cells = 0L,
    min_features = 0L,
    feature_type = "Gene Expression",
    gene_column = 2L,
    strip_suffix = FALSE,
    verbose = TRUE
) {
  .qc_require_package("SeuratObject")
  .qc_require_package("Matrix")

  if (inherits(x, "Seurat")) {
    .qc_message(verbose, "[INPUT] Seurat object detected.")

    return(
      .validate_seurat_input(
        object = x,
        assay = assay
      )
    )
  }

  if (.is_count_matrix(x)) {
    .qc_message(verbose, "[INPUT] Count matrix detected.")

    counts <- .validate_count_matrix(x)

    return(
      .matrix_to_seurat(
        counts = counts,
        sample_id = sample_id,
        assay = assay,
        min_cells = min_cells,
        min_features = min_features
      )
    )
  }

  if (is.character(x) && length(x) == 1L) {
    return(
      .path_to_seurat(
        path = x,
        sample_id = sample_id,
        assay = assay,
        min_cells = min_cells,
        min_features = min_features,
        feature_type = feature_type,
        gene_column = gene_column,
        strip_suffix = strip_suffix,
        verbose = verbose
      )
    )
  }

  stop(
    paste(
      "`x` must be a Seurat object, count matrix,",
      "10x matrix directory, 10x HDF5 file,",
      "or Cell Ranger output directory."
    ),
    call. = FALSE
  )
}
