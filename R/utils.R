`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


.qc_message <- function(verbose, ...) {
  if (isTRUE(verbose)) {
    message(...)
  }

  invisible(NULL)
}


.qc_require_package <- function(package) {
  if (!requireNamespace(
    package,
    quietly = TRUE
  )) {
    stop(
      "Package `",
      package,
      "` is required for this operation.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.qc_bind_rows <- function(x) {
  x <- x[
    !vapply(
      x,
      is.null,
      logical(1)
    )
  ]

  if (length(x) == 0L) {
    return(data.frame())
  }

  result <- do.call(
    rbind,
    x
  )

  rownames(result) <- NULL
  result
}


.qc_collect_failures <- function(results) {
  if (length(results) == 0L) {
    return(
      data.frame(
        sample = character(),
        stage = character(),
        error = character(),
        stringsAsFactors = FALSE
      )
    )
  }

  do.call(
    rbind,
    lapply(
      results,
      function(result) {
        data.frame(
          sample = result$sample,
          stage = result$stage,
          error = result$error,
          stringsAsFactors = FALSE
        )
      }
    )
  )
}


.qc_apply_parallel <- function(
    indices,
    worker,
    n_cores,
    verbose
) {
  n_cores <- min(
    as.integer(n_cores),
    length(indices)
  )

  .qc_message(
    verbose,
    "[QC] Using ",
    n_cores,
    " worker(s)."
  )

  if (.Platform$OS.type == "unix") {
    return(
      parallel::mclapply(
        indices,
        worker,
        mc.cores = n_cores,
        mc.preschedule = FALSE
      )
    )
  }

  warning(
    paste(
      "Parallel processing is disabled on this platform",
      "to avoid copying large Seurat objects to PSOCK workers."
    ),
    call. = FALSE
  )

  lapply(indices, worker)
}



.qc_resolve_assay <- function(
    object,
    assay
) {
  available <- names(object@assays)

  if (assay %in% available) {
    return(assay)
  }

  if (length(available) == 1L) {
    warning(
      "Assay `",
      assay,
      "` was not found. Using `",
      available[[1L]],
      "`.",
      call. = FALSE
    )

    return(available[[1L]])
  }

  stop(
    "Assay `",
    assay,
    "` was not found. Available assays: ",
    paste(available, collapse = ", "),
    "."
  )
}


.qc_resolve_sample <- function(
    object,
    sample_id,
    sample_col
) {
  if (!is.null(sample_id) &&
      length(sample_id) == 1L &&
      nzchar(sample_id)) {
    return(sample_id)
  }

  metadata <- object[[]]

  if (sample_col %in% colnames(metadata)) {
    values <- unique(
      as.character(metadata[[sample_col]])
    )

    values <- values[
      !is.na(values) &
        nzchar(values)
    ]

    if (length(values) == 1L) {
      return(values)
    }
  }

  "Sample_1"
}


.qc_gene_patterns <- function(
    species,
    mito_pat,
    ribo_pat
) {
  default_mito <- switch(
    species,
    human = "^MT-",
    mouse = "^mt-"
  )

  default_ribo <- switch(
    species,
    human = "^RP[LS]",
    mouse = "^Rp[ls]"
  )

  list(
    mito = if (is.null(mito_pat)) {
      default_mito
    } else {
      mito_pat
    },
    ribo = if (is.null(ribo_pat)) {
      default_ribo
    } else {
      ribo_pat
    }
  )
}



.qc_removal_reason <- function(filter_table) {
  reasons <- rep(
    "retained",
    nrow(filter_table)
  )

  reason_columns <- c(
    fail_min_feat = "low_features",
    fail_max_feat = "high_features",
    fail_min_umi = "low_umi",
    fail_max_umi = "high_umi",
    fail_mito = "high_mito",
    fail_ribo = "high_ribo",
    fail_drop = "high_dropout",
    fail_g2u = "low_genes_per_umi"
  )

  for (i in seq_len(nrow(filter_table))) {
    current <- names(reason_columns)[
      vapply(
        names(reason_columns),
        function(column) {
          isTRUE(filter_table[[column]][i])
        },
        logical(1)
      )
    ]

    current_reasons <- unname(
      reason_columns[current]
    )

    if (isTRUE(filter_table$doublet[i])) {
      current_reasons <- c(
        current_reasons,
        "doublet"
      )
    }

    if (length(current_reasons) > 0L) {
      reasons[i] <- paste(
        current_reasons,
        collapse = ";"
      )
    }
  }

  reasons
}




.qc_combine_objects <- function(
    objects,
    merge,
    join_layers,
    assay,
    verbose
) {
  if (length(objects) == 1L) {
    return(objects[[1L]])
  }

  if (!merge) {
    return(objects)
  }

  .qc_message(
    verbose,
    "[QC] Merging ",
    length(objects),
    " filtered objects."
  )

  add_ids <- names(objects)

  merged <- merge(
    x = objects[[1L]],
    y = objects[-1L],
    add.cell.ids = add_ids,
    merge.data = TRUE
  )

  if (
    join_layers &&
    "JoinLayers" %in%
    getNamespaceExports("SeuratObject")
  ) {
    merged <- tryCatch(
      SeuratObject::JoinLayers(
        object = merged,
        assay = assay
      ),
      error = function(e) {
        warning(
          "Objects were merged, but layers were not joined: ",
          conditionMessage(e),
          call. = FALSE
        )

        merged
      }
    )
  }

  merged
}



.qc_write_results <- function(
    objects,
    summary,
    thresholds,
    doublets,
    removed_cells,
    failures,
    parameters,
    outdir,
    save_object
) {
  if (is.null(outdir)) {
    return(list())
  }

  files <- list()

  files$summary <- file.path(
    outdir,
    "QC_summary.tsv"
  )

  files$thresholds <- file.path(
    outdir,
    "QC_thresholds.tsv"
  )

  files$doublets <- file.path(
    outdir,
    "Doublet_summary.tsv"
  )

  files$removed_cells <- file.path(
    outdir,
    "QC_removed_cells.tsv"
  )

  files$parameters <- file.path(
    outdir,
    "QC_parameters.rds"
  )

  utils::write.table(
    summary,
    files$summary,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  utils::write.table(
    thresholds,
    files$thresholds,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  utils::write.table(
    doublets,
    files$doublets,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  utils::write.table(
    removed_cells,
    files$removed_cells,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  saveRDS(
    parameters,
    files$parameters
  )

  if (nrow(failures) > 0L) {
    files$failures <- file.path(
      outdir,
      "QC_failures.tsv"
    )

    utils::write.table(
      failures,
      files$failures,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  if (save_object) {
    object_dir <- file.path(
      outdir,
      "objects"
    )

    dir.create(
      object_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    files$objects <- character(
      length(objects)
    )

    for (i in seq_along(objects)) {
      safe_name <- gsub(
        "[^A-Za-z0-9._-]",
        "_",
        names(objects)[i]
      )

      object_file <- file.path(
        object_dir,
        paste0(
          safe_name,
          "_qc.rds"
        )
      )

      saveRDS(
        objects[[i]],
        object_file
      )

      files$objects[i] <- object_file
    }

    names(files$objects) <- names(objects)
  }

  files
}


.qc_validate_arguments <- function(
    min_feat,
    min_umi,
    mad_n,
    max_mito,
    max_ribo,
    max_drop,
    min_g2u,
    dbl_score_thr,
    method,
    fixed_thr,
    n_cores
) {
  nonnegative <- list(
    min_feat = min_feat,
    min_umi = min_umi,
    mad_n = mad_n,
    max_mito = max_mito,
    max_ribo = max_ribo
  )

  for (name in names(nonnegative)) {
    value <- nonnegative[[name]]

    if (
      length(value) != 1L ||
      !is.numeric(value) ||
      !is.finite(value) ||
      value < 0
    ) {
      stop(
        "`",
        name,
        "` must be one finite non-negative number.",
        call. = FALSE
      )
    }
  }

  proportions <- list(
    max_drop = max_drop,
    dbl_score_thr = dbl_score_thr
  )

  for (name in names(proportions)) {
    value <- proportions[[name]]

    if (
      length(value) != 1L ||
      !is.numeric(value) ||
      !is.finite(value) ||
      value < 0 ||
      value > 1
    ) {
      stop(
        "`",
        name,
        "` must be between 0 and 1.",
        call. = FALSE
      )
    }
  }

  if (
    length(min_g2u) != 1L ||
    !is.numeric(min_g2u) ||
    !is.finite(min_g2u)
  ) {
    stop(
      "`min_g2u` must be one finite numeric value.",
      call. = FALSE
    )
  }

  if (
    length(n_cores) != 1L ||
    !is.numeric(n_cores) ||
    !is.finite(n_cores) ||
    n_cores < 1
  ) {
    stop(
      "`n_cores` must be a positive integer.",
      call. = FALSE
    )
  }

  if (method == "fixed") {
    required <- c(
      "max_feat",
      "max_umi"
    )

    if (
      !is.list(fixed_thr) ||
      !all(required %in% names(fixed_thr))
    ) {
      stop(
        paste(
          "`fixed_thr` must be a named list containing",
          "`max_feat` and `max_umi`."
        ),
        call. = FALSE
      )
    }

    values <- unlist(
      fixed_thr[required],
      use.names = FALSE
    )

    if (
      any(!is.finite(values)) ||
      any(values <= 0)
    ) {
      stop(
        "Fixed upper thresholds must be positive finite values.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}



#' Save a lightweight QC result
#'
#' Removes Seurat objects and removed-cell records before saving an
#' `scqc_result`.
#'
#' @param x An object returned by [run_qc()].
#' @param file Output file name.
#' @param outdir Output directory.
#'
#' @return Invisibly returns the saved file path.
#'
#' @export
save_qc <- function(x, file, outdir = ".") {
  if (!inherits(x, "scqc_result")) {
    stop("`x` must be an `scqc_result`.", call. = FALSE)
  }

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  x[c("objects", "object", "removed_cells")] <- NULL
  path <- file.path(outdir, file)

  saveRDS(x, path)
  invisible(path)
}

