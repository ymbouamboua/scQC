.is_count_matrix <- function(x) {
  is.matrix(x) ||
    inherits(x, "Matrix") ||
    inherits(x, "dgCMatrix")
}


.validate_count_matrix <- function(counts) {
  if (length(dim(counts)) != 2L) {
    stop(
      "The count input must be a two-dimensional matrix.",
      call. = FALSE
    )
  }

  if (nrow(counts) == 0L || ncol(counts) == 0L) {
    stop(
      "The count matrix must contain at least one feature and one cell.",
      call. = FALSE
    )
  }

  if (is.null(rownames(counts))) {
    stop(
      "The count matrix must have feature names as row names.",
      call. = FALSE
    )
  }

  if (is.null(colnames(counts))) {
    stop(
      "The count matrix must have cell barcodes as column names.",
      call. = FALSE
    )
  }

  if (anyDuplicated(rownames(counts))) {
    warning(
      "Duplicate feature names were made unique.",
      call. = FALSE
    )

    rownames(counts) <- make.unique(rownames(counts))
  }

  if (anyDuplicated(colnames(counts))) {
    stop(
      "Cell names must be unique.",
      call. = FALSE
    )
  }

  if (anyNA(counts)) {
    stop(
      "The count matrix contains missing values.",
      call. = FALSE
    )
  }

  if (any(counts < 0)) {
    stop(
      "The count matrix contains negative values.",
      call. = FALSE
    )
  }

  if (!inherits(counts, "dgCMatrix")) {
    counts <- methods::as(
      Matrix::Matrix(counts, sparse = TRUE),
      "dgCMatrix"
    )
  }

  counts
}


.matrix_to_seurat <- function(
    counts,
    sample_id,
    assay,
    min_cells,
    min_features
) {
  sample_id <- sample_id %||% "Sample1"

  object <- SeuratObject::CreateSeuratObject(
    counts = counts,
    assay = assay,
    project = sample_id,
    min.cells = min_cells,
    min.features = min_features
  )

  object$sample <- sample_id
  object$orig.ident <- sample_id

  object
}


.validate_seurat_input <- function(object, assay) {
  if (ncol(object) == 0L) {
    stop(
      "The Seurat object contains no cells.",
      call. = FALSE
    )
  }

  available_assays <- names(object@assays)

  if (!assay %in% available_assays) {
    if (length(available_assays) == 1L) {
      assay <- available_assays[[1L]]

      warning(
        "Assay `RNA` was not found. Using assay `",
        assay,
        "`.",
        call. = FALSE
      )
    } else {
      stop(
        "Assay `",
        assay,
        "` was not found. Available assays: ",
        paste(available_assays, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  SeuratObject::DefaultAssay(object) <- assay

  object
}



.path_to_seurat <- function(
    path,
    sample_id,
    assay,
    min_cells,
    min_features,
    feature_type,
    gene_column,
    strip_suffix,
    verbose
) {
  .qc_require_package("Seurat")

  path <- path.expand(path)

  if (!file.exists(path)) {
    stop(
      "Input path does not exist: ",
      path,
      call. = FALSE
    )
  }

  resolved <- .resolve_10x_path(path)

  .qc_message(
    verbose,
    "[INPUT] Detected ",
    resolved$type,
    ": ",
    resolved$path
  )

  counts <- switch(
    resolved$type,

    "10x_h5" = Seurat::Read10X_h5(
      filename = resolved$path,
      use.names = TRUE,
      unique.features = TRUE
    ),

    "10x_directory" = Seurat::Read10X(
      data.dir = resolved$path,
      gene.column = gene_column,
      unique.features = TRUE,
      strip.suffix = strip_suffix
    ),

    stop(
      "Unsupported resolved input type.",
      call. = FALSE
    )
  )

  counts <- .select_feature_type(
    counts = counts,
    feature_type = feature_type
  )

  counts <- .validate_count_matrix(counts)

  if (is.null(sample_id)) {
    sample_id <- .infer_sample_id(
      original_path = path,
      resolved_path = resolved$path
    )
  }

  .matrix_to_seurat(
    counts = counts,
    sample_id = sample_id,
    assay = assay,
    min_cells = min_cells,
    min_features = min_features
  )
}


.resolve_10x_path <- function(path) {
  if (file.info(path)$isdir) {
    normalized_path <- normalizePath(
      path,
      winslash = "/",
      mustWork = TRUE
    )

    direct_h5 <- c(
      file.path(
        normalized_path,
        "filtered_feature_bc_matrix.h5"
      ),
      file.path(
        normalized_path,
        "raw_feature_bc_matrix.h5"
      )
    )

    direct_matrix_dirs <- c(
      normalized_path,
      file.path(
        normalized_path,
        "filtered_feature_bc_matrix"
      ),
      file.path(
        normalized_path,
        "raw_feature_bc_matrix"
      )
    )

    outs_h5 <- c(
      file.path(
        normalized_path,
        "outs",
        "filtered_feature_bc_matrix.h5"
      ),
      file.path(
        normalized_path,
        "outs",
        "raw_feature_bc_matrix.h5"
      )
    )

    outs_matrix_dirs <- c(
      file.path(
        normalized_path,
        "outs",
        "filtered_feature_bc_matrix"
      ),
      file.path(
        normalized_path,
        "outs",
        "raw_feature_bc_matrix"
      )
    )

    h5_candidates <- c(
      direct_h5,
      outs_h5
    )

    existing_h5 <- h5_candidates[
      file.exists(h5_candidates)
    ]

    if (length(existing_h5) > 0L) {
      return(
        list(
          type = "10x_h5",
          path = existing_h5[[1L]]
        )
      )
    }

    matrix_candidates <- c(
      direct_matrix_dirs,
      outs_matrix_dirs
    )

    valid_matrix_dirs <- matrix_candidates[
      vapply(
        matrix_candidates,
        .is_10x_matrix_directory,
        logical(1)
      )
    ]

    if (length(valid_matrix_dirs) > 0L) {
      return(
        list(
          type = "10x_directory",
          path = valid_matrix_dirs[[1L]]
        )
      )
    }

    stop(
      paste0(
        "No compatible 10x matrix was found under: ",
        normalized_path,
        "\nExpected either a filtered/raw feature matrix HDF5 file ",
        "or a directory containing matrix.mtx, features.tsv, ",
        "and barcodes.tsv."
      ),
      call. = FALSE
    )
  }

  extension <- tolower(tools::file_ext(path))

  if (extension %in% c("h5", "hdf5")) {
    return(
      list(
        type = "10x_h5",
        path = normalizePath(
          path,
          winslash = "/",
          mustWork = TRUE
        )
      )
    )
  }

  stop(
    "Unsupported input file: ",
    path,
    call. = FALSE
  )
}


.is_10x_matrix_directory <- function(path) {
  if (!dir.exists(path)) {
    return(FALSE)
  }

  matrix_files <- c(
    "matrix.mtx",
    "matrix.mtx.gz"
  )

  feature_files <- c(
    "features.tsv",
    "features.tsv.gz",
    "genes.tsv",
    "genes.tsv.gz"
  )

  barcode_files <- c(
    "barcodes.tsv",
    "barcodes.tsv.gz"
  )

  has_matrix <- any(
    file.exists(
      file.path(path, matrix_files)
    )
  )

  has_features <- any(
    file.exists(
      file.path(path, feature_files)
    )
  )

  has_barcodes <- any(
    file.exists(
      file.path(path, barcode_files)
    )
  )

  has_matrix && has_features && has_barcodes
}


.select_feature_type <- function(
    counts,
    feature_type = "Gene Expression"
) {
  if (!is.list(counts)) {
    return(counts)
  }

  if (length(counts) == 0L) {
    stop(
      "The 10x input did not contain any count matrices.",
      call. = FALSE
    )
  }

  available <- names(counts)

  if (
    !is.null(available) &&
    feature_type %in% available
  ) {
    return(counts[[feature_type]])
  }

  if (length(counts) == 1L) {
    return(counts[[1L]])
  }

  stop(
    "The 10x input contains multiple feature types: ",
    paste(available, collapse = ", "),
    ". Specify one using `feature_type`.",
    call. = FALSE
  )
}


.qc_normalize_inputs <- function(
    x,
    sample_id = NULL
) {
  is_single_input <- (
    inherits(x, "Seurat") ||
      is.matrix(x) ||
      inherits(x, "Matrix") ||
      (
        is.character(x) &&
          length(x) == 1L
      )
  )

  if (is_single_input) {
    name <- sample_id

    if (is.null(name) || !nzchar(name)) {
      name <- .qc_infer_input_name(x)
    }

    result <- list(x)
    names(result) <- name

    return(result)
  }

  if (!is.list(x) || length(x) == 0L) {
    stop(
      paste(
        "`x` must be a supported single input",
        "or a non-empty list of supported inputs."
      ),
      call. = FALSE
    )
  }

  if (!is.null(sample_id)) {
    if (length(sample_id) != length(x)) {
      stop(
        "For multiple inputs, `sample_id` must have one value per input.",
        call. = FALSE
      )
    }

    names(x) <- sample_id
  }

  if (is.null(names(x))) {
    names(x) <- paste0(
      "Sample_",
      seq_along(x)
    )
  }

  missing_names <- (
    is.na(names(x)) |
      !nzchar(names(x))
  )

  names(x)[missing_names] <- paste0(
    "Sample_",
    which(missing_names)
  )

  names(x) <- make.unique(names(x))

  x
}


.infer_sample_id <- function(
    original_path,
    resolved_path
) {
  original_path <- normalizePath(
    original_path,
    winslash = "/",
    mustWork = TRUE
  )

  resolved_path <- normalizePath(
    resolved_path,
    winslash = "/",
    mustWork = TRUE
  )

  container_names <- c(
    "outs",
    "filtered_feature_bc_matrix",
    "raw_feature_bc_matrix"
  )

  if (dir.exists(original_path)) {
    current_path <- original_path
    sample_id <- basename(current_path)

    while (
      sample_id %in% container_names &&
      dirname(current_path) != current_path
    ) {
      current_path <- dirname(current_path)
      sample_id <- basename(current_path)
    }

    return(sample_id)
  }

  filename <- basename(resolved_path)

  if (
    filename %in% c(
      "filtered_feature_bc_matrix.h5",
      "raw_feature_bc_matrix.h5"
    )
  ) {
    current_path <- dirname(resolved_path)
    sample_id <- basename(current_path)

    while (
      sample_id %in% container_names &&
      dirname(current_path) != current_path
    ) {
      current_path <- dirname(current_path)
      sample_id <- basename(current_path)
    }

    return(sample_id)
  }

  sample_id <- sub(
    "_(?:filtered|raw)_feature_bc_matrix\\.h5$",
    "",
    filename,
    ignore.case = TRUE
  )

  sample_id
}


.qc_infer_input_name <- function(
    x,
    fallback = "Sample_1"
) {
  if (inherits(x, "Seurat")) {
    metadata <- x[[]]

    if ("orig.ident" %in% colnames(metadata)) {
      values <- unique(
        as.character(metadata$orig.ident)
      )

      values <- values[
        !is.na(values) &
          nzchar(values)
      ]

      if (length(values) == 1L) {
        return(values)
      }
    }

    return(fallback)
  }

  if (
    is.character(x) &&
    length(x) == 1L &&
    file.exists(x)
  ) {
    resolved_path <- tryCatch(
      .resolve_10x_path(x),
      error = function(e) NULL
    )

    if (is.null(resolved_path)) {
      return(fallback)
    }

    return(
      .infer_sample_id(
        original_path = x,
        resolved_path = resolved_path
      )
    )
  }

  fallback
}
