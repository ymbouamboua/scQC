.qc_process_object <- function(
    object,
    sample_id,
    assay,
    layer,
    min_feat,
    min_umi,
    mad_n,
    max_mito,
    calc_ribo,
    max_ribo,
    calc_drop,
    max_drop,
    log_g2u,
    min_g2u,
    rm_dbl,
    dbl_method,
    dbl_score_thr,
    method,
    fixed_thr,
    mito_pat,
    ribo_pat,
    sample_col,
    seed,
    verbose,
    decontaminate=FALSE, ambient_method="decontX", ambient_raw=NULL,
    ambient_batch=NULL, ambient_clusters=NULL, ambient_args=list(),
    keep_raw=FALSE, ambient_assay="RNAdecont", feature_type="Gene Expression"
) {
  started <- Sys.time()
  current_stage <- "quality_control"

  tryCatch(
    {
      assay <- .qc_resolve_assay(
        object = object,
        assay = assay
      )

      SeuratObject::DefaultAssay(object) <- assay

      sample <- .qc_resolve_sample(
        object = object,
        sample_id = sample_id,
        sample_col = sample_col
      )

      called_cells <- colnames(object)
      batches <- clusters <- NULL
      if (decontaminate) {
        batches <- .qc_ambient_groups(object,ambient_batch)
        if (!is.null(ambient_clusters)) clusters <- .qc_ambient_groups(object,ambient_clusters)
      }
      object$sample <- sample

      input_cells <- ncol(object)

      if (input_cells == 0L) {
        stop("The input object contains no cells.")
      }

      .qc_message(
        verbose,
        "[QC] ",
        sample,
        ": ",
        input_cells,
        " input cells."
      )

      object <- .qc_add_metrics(
        object = object,
        assay = assay,
        layer = layer,
        mito_pat = mito_pat,
        ribo_pat = ribo_pat,
        calc_ribo = calc_ribo,
        calc_drop = calc_drop,
        log_g2u = log_g2u
      )

      thresholds <- .qc_calculate_thresholds(
        object = object,
        method = method,
        mad_n = mad_n,
        fixed_thr = fixed_thr
      )

      filter_table <- .qc_filter_table(
        object = object,
        min_feat = min_feat,
        max_feat = thresholds$max_feat,
        min_umi = min_umi,
        max_umi = thresholds$max_umi,
        max_mito = max_mito,
        calc_ribo = calc_ribo,
        max_ribo = max_ribo,
        calc_drop = calc_drop,
        max_drop = max_drop,
        log_g2u = log_g2u,
        min_g2u = min_g2u
      )

      keep_cells <- filter_table$cell[
        filter_table$pass_qc
      ]

      if (length(keep_cells) == 0L) {
        stop(
          "No cells passed the requested QC thresholds."
        )
      }

      object <- object[, keep_cells]

      after_qc_cells <- ncol(object)

      doublet_result <- list(
        object = object,
        detected = 0L,
        removed = 0L,
        method = "not_run",
        threshold = NA_real_,
        doublet_cells = character()
      )

      if (rm_dbl) {
        doublet_result <- .qc_detect_doublets(
          object = object,
          assay = assay,
          dbl_method = dbl_method,
          dbl_score_thr = dbl_score_thr,
          seed = seed
        )

        object <- doublet_result$object
      }

      after_doublet_cells <- ncol(object)
      corrected_mito_cells <- character()

      ambient_result <- list(summary=data.frame(),cells=data.frame())
      if (decontaminate) {
        current_stage <- "ambient"
        ambient_result <- .qc_decontaminate(
          object,assay,layer,sample,ambient_method,ambient_raw,batches,clusters,
          called_cells,keep_raw,ambient_assay,ambient_args,feature_type,seed,verbose)
        object <- ambient_result$object
        corrected_assay <- if (keep_raw) ambient_assay else assay
        object <- .qc_add_metrics(object,corrected_assay,"counts",mito_pat,ribo_pat,
          calc_ribo,calc_drop,FALSE)
        if (log_g2u) {
          cm <- .qc_get_layer(object,corrected_assay,"counts")
          object$log10_genes_per_umi <- log10(pmax(Matrix::colSums(cm>0)/pmax(Matrix::colSums(cm),1),.Machine$double.xmin))
        }

        # Ambient correction can change the mitochondrial fraction. Enforce
        # max_mito on corrected counts too, so the returned object and its
        # refreshed pct_mt metric remain consistent with the requested limit.
        corrected_mito_cells <- colnames(object)[
          is.na(object$pct_mt) | object$pct_mt > max_mito
        ]

        if (length(corrected_mito_cells) > 0L) {
          object <- object[, setdiff(colnames(object), corrected_mito_cells)]
        }

        if (ncol(object) == 0L) {
          stop("No cells met `max_mito` after ambient RNA correction.")
        }

        current_stage <- "quality_control"
      }

      final_cells <- ncol(object)

      filter_table$doublet <- (
        filter_table$cell %in%
          doublet_result$doublet_cells
      )

      filter_table$fail_corrected_mito <- (
        filter_table$cell %in% corrected_mito_cells
      )

      filter_table$retained <- (
        filter_table$cell %in%
          colnames(object)
      )

      filter_table$sample <- sample

      filter_table$removal_reason <- .qc_removal_reason(
        filter_table
      )

      elapsed <- as.numeric(
        difftime(
          Sys.time(),
          started,
          units = "secs"
        )
      )

      summary <- data.frame(
        sample = sample,
        input_cells = input_cells,
        after_qc_cells = after_qc_cells,
        after_doublet_cells = after_doublet_cells,
        final_cells = final_cells,
        qc_removed = input_cells - after_qc_cells,
        doublets_removed = after_qc_cells - after_doublet_cells,
        ambient_mito_removed = length(corrected_mito_cells),
        total_removed = input_cells - final_cells,
        qc_removed_pct = round(
          100 * (input_cells - after_qc_cells) /
            input_cells,
          3
        ),
        doublets_removed_pct = round(
          100 * (after_qc_cells - after_doublet_cells) /
            input_cells,
          3
        ),
        ambient_mito_removed_pct = round(
          100 * length(corrected_mito_cells) /
            input_cells,
          3
        ),
        total_removed_pct = round(
          100 * (input_cells - final_cells) /
            input_cells,
          3
        ),
        elapsed_seconds = round(elapsed, 3),
        stringsAsFactors = FALSE
      )

      threshold_table <- data.frame(
        sample = sample,
        method = method,
        mad_n = if (identical(method, "MAD")) {
          as.numeric(mad_n)
        } else {
          NA_real_
        },
        min_feat = min_feat,
        max_feat = thresholds$max_feat,
        min_umi = min_umi,
        max_umi = thresholds$max_umi,
        max_mito = max_mito,
        max_ribo = if (calc_ribo) {
          max_ribo
        } else {
          NA_real_
        },
        max_drop = if (calc_drop) {
          max_drop
        } else {
          NA_real_
        },
        min_log10_g2u = if (log_g2u) {
          min_g2u
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )

      doublet_table <- data.frame(
        sample = sample,
        method = doublet_result$method,
        detected = doublet_result$detected,
        removed = doublet_result$removed,
        score_threshold = doublet_result$threshold,
        stringsAsFactors = FALSE
      )

      .qc_message(
        verbose,
        "[QC] ",
        sample,
        ": retained ",
        final_cells,
        "/",
        input_cells,
        " cells."
      )

      list(
        status = "success",
        sample = sample,
        object = object,
        summary = summary,
        thresholds = threshold_table,
        doublets = doublet_table,
        ambient = ambient_result$summary,
        ambient_cells = ambient_result$cells,
        removed_cells = filter_table[
          !filter_table$retained,
          ,
          drop = FALSE
        ]
      )
    },
    error = function(e) {
      list(
        status = "failed",
        sample = sample_id,
        stage = current_stage,
        error = conditionMessage(e)
      )
    }
  )
}



.qc_add_metrics <- function(
    object,
    assay,
    layer,
    mito_pat,
    ribo_pat,
    calc_ribo,
    calc_drop,
    log_g2u
) {
  object$pct_mt <- Seurat::PercentageFeatureSet(
    object = object,
    pattern = mito_pat,
    assay = assay
  )

  if (calc_ribo) {
    object$pct_ribo <- Seurat::PercentageFeatureSet(
      object = object,
      pattern = ribo_pat,
      assay = assay
    )
  }

  if (calc_drop) {
    counts <- .qc_get_layer(
      object = object,
      assay = assay,
      layer = layer
    )

    detected_features <- Matrix::colSums(
      counts > 0
    )

    object$dropout <- 1 -
      detected_features / nrow(counts)

    rm(counts, detected_features)
  }

  if (log_g2u) {
    metadata <- object[[]]

    genes_per_umi <- metadata$nFeature_RNA /
      pmax(metadata$nCount_RNA, 1)

    object$log10_genes_per_umi <- log10(
      pmax(
        genes_per_umi,
        .Machine$double.xmin
      )
    )
  }

  object
}



.qc_get_layer <- function(
    object,
    assay,
    layer
) {
  counts <- tryCatch(
    SeuratObject::LayerData(
      object = object,
      assay = assay,
      layer = layer
    ),
    error = function(e) {
      tryCatch(
        SeuratObject::GetAssayData(
          object = object,
          assay = assay,
          layer = layer
        ),
        error = function(e2) {
          stop(
            "Unable to retrieve layer `",
            layer,
            "` from assay `",
            assay,
            "`: ",
            conditionMessage(e2)
          )
        }
      )
    }
  )

  if (nrow(counts) == 0L ||
      ncol(counts) == 0L) {
    stop("The selected count layer is empty.")
  }

  counts
}


.qc_calculate_thresholds <- function(
    object,
    method,
    mad_n,
    fixed_thr
) {
  metadata <- object[[]]

  if (method == "MAD") {
    max_feat <- .qc_upper_mad(
      metadata$nFeature_RNA,
      mad_n
    )

    max_umi <- .qc_upper_mad(
      metadata$nCount_RNA,
      mad_n
    )
  } else if (method == "fixed") {
    max_feat <- fixed_thr$max_feat
    max_umi <- fixed_thr$max_umi
  } else {
    max_feat <- Inf
    max_umi <- Inf
  }

  list(
    max_feat = max_feat,
    max_umi = max_umi
  )
}


.qc_upper_mad <- function(x, mad_n) {
  x <- x[
    is.finite(x) &
      !is.na(x)
  ]

  if (length(x) == 0L) {
    return(Inf)
  }

  centre <- stats::median(
    x,
    na.rm = TRUE
  )

  dispersion <- stats::mad(
    x,
    center = centre,
    constant = 1.4826,
    na.rm = TRUE
  )

  if (!is.finite(dispersion) ||
      dispersion <= 0) {
    return(Inf)
  }

  centre + mad_n * dispersion
}



.qc_filter_table <- function(
    object,
    min_feat,
    max_feat,
    min_umi,
    max_umi,
    max_mito,
    calc_ribo,
    max_ribo,
    calc_drop,
    max_drop,
    log_g2u,
    min_g2u
) {
  metadata <- object[[]]

  result <- data.frame(
    cell = rownames(metadata),
    nFeature_RNA = metadata$nFeature_RNA,
    nCount_RNA = metadata$nCount_RNA,
    pct_mt = metadata$pct_mt,
    stringsAsFactors = FALSE
  )

  result$fail_min_feat <- (
    is.na(result$nFeature_RNA) |
      result$nFeature_RNA < min_feat
  )

  result$fail_max_feat <- (
    is.na(result$nFeature_RNA) |
      result$nFeature_RNA > max_feat
  )

  result$fail_min_umi <- (
    is.na(result$nCount_RNA) |
      result$nCount_RNA < min_umi
  )

  result$fail_max_umi <- (
    is.na(result$nCount_RNA) |
      result$nCount_RNA > max_umi
  )

  result$fail_mito <- (
    is.na(result$pct_mt) |
      result$pct_mt > max_mito
  )

  if (calc_ribo) {
    result$pct_ribo <- metadata$pct_ribo

    result$fail_ribo <- (
      is.na(result$pct_ribo) |
        result$pct_ribo > max_ribo
    )
  } else {
    result$fail_ribo <- FALSE
  }

  if (calc_drop) {
    result$dropout <- metadata$dropout

    result$fail_drop <- (
      is.na(result$dropout) |
        result$dropout > max_drop
    )
  } else {
    result$fail_drop <- FALSE
  }

  if (log_g2u) {
    result$log10_genes_per_umi <-
      metadata$log10_genes_per_umi

    result$fail_g2u <- (
      is.na(result$log10_genes_per_umi) |
        result$log10_genes_per_umi < min_g2u
    )
  } else {
    result$fail_g2u <- FALSE
  }

  failure_columns <- grep(
    "^fail_",
    colnames(result),
    value = TRUE
  )

  result$pass_qc <- !apply(
    result[, failure_columns, drop = FALSE],
    1,
    any
  )

  result
}



.qc_detect_doublets <- function(
    object,
    assay,
    dbl_method,
    dbl_score_thr,
    seed
) {
  .qc_require_package("Seurat")
  .qc_require_package("scDblFinder")
  .qc_require_package("SingleCellExperiment")
  .qc_require_package("SummarizedExperiment")

  if (ncol(object) < 50L) {
    warning(
      paste(
        "Fewer than 50 cells remain.",
        "Doublet detection was skipped."
      ),
      call. = FALSE
    )

    return(
      list(
        object = object,
        detected = 0L,
        removed = 0L,
        method = "skipped_small_object",
        threshold = NA_real_,
        doublet_cells = character()
      )
    )
  }

  set.seed(seed)

  sce <- Seurat::as.SingleCellExperiment(
    object,
    assay = assay
  )

  sce <- scDblFinder::scDblFinder(sce)

  col_data <- SummarizedExperiment::colData(sce)

  score <- col_data$scDblFinder.score
  class <- as.character(
    col_data$scDblFinder.class
  )

  if (is.null(score) || is.null(class)) {
    stop(
      "scDblFinder did not return its expected score and class columns."
    )
  }

  object$scDblFinder_score <- score
  object$scDblFinder_class <- class

  if (dbl_method == "class") {
    is_doublet <- class == "doublet"
    threshold <- NA_real_
  } else {
    is_doublet <- score > dbl_score_thr
    threshold <- dbl_score_thr
  }

  is_doublet[is.na(is_doublet)] <- FALSE

  doublet_cells <- colnames(object)[
    is_doublet
  ]

  object <- object[
    ,
    !is_doublet,
    drop = FALSE
  ]

  list(
    object = object,
    detected = length(doublet_cells),
    removed = length(doublet_cells),
    method = dbl_method,
    threshold = threshold,
    doublet_cells = doublet_cells
  )
}



