#' Run quality control on single-cell RNA-seq data
#'
#' Reads and quality-controls a 10x matrix, count matrix, Seurat object,
#' or a named list of supported inputs.
#'
#' @param x A Seurat object, feature-by-cell count matrix, 10x directory,
#'   10x HDF5 file, Cell Ranger output directory, or named list of these.
#' @param sample_id Optional sample identifier. For multiple inputs, this can
#'   be a character vector matching the length of `x`.
#' @param assay Assay containing raw counts.
#' @param layer Layer containing raw counts.
#' @param feature_type Feature type selected from multimodal 10x inputs.
#' @param species Species used to define mitochondrial and ribosomal patterns.
#' @param min_feat Minimum detected features per cell.
#' @param min_umi Minimum UMIs per cell.
#' @param mad_n Number of MADs used for adaptive upper thresholds.
#' @param max_mito Maximum mitochondrial percentage.
#' @param calc_ribo Calculate and filter ribosomal percentage.
#' @param max_ribo Maximum ribosomal percentage.
#' @param calc_drop Calculate and filter dropout fraction.
#' @param max_drop Maximum dropout fraction.
#' @param log_g2u Calculate and filter log10 genes-per-UMI.
#' @param min_g2u Minimum log10 genes-per-UMI value.
#' @param rm_dbl Detect and remove doublets using scDblFinder.
#' @param dbl_method Use scDblFinder classification or a manual score cutoff.
#' @param dbl_score_thr Score cutoff used when `dbl_method = "score"`.
#' @param method Upper threshold method: `"MAD"`, `"fixed"`, or `"none"`.
#' @param fixed_thr Named list containing `max_feat` and `max_umi`.
#' @param mito_pat Optional mitochondrial feature regular expression.
#' @param ribo_pat Optional ribosomal feature regular expression.
#' @param sample_col Metadata column containing the sample identifier.
#' @param outdir Output directory. Set to `NULL` to disable file writing.
#' @param save_object Save filtered Seurat objects as RDS files.
#' @param merge Merge multiple successfully processed objects.
#' @param join_layers Join Seurat v5 layers after merging.
#' @param parallel Process multiple inputs in parallel.
#' @param n_cores Number of parallel workers.
#' @param seed Random seed used by scDblFinder.
#' @param return Return `"result"` or only the processed `"object"`.
#' @param verbose Print progress messages.
#'
#' @return A structured `scqc_result`, a Seurat object, or a list of objects.
#'
#' @export
run_qc <- function(
    x,
    sample_id = NULL,
    assay = "RNA",
    layer = "counts",
    feature_type = "Gene Expression",
    species = c("human", "mouse"),
    min_feat = 200,
    min_umi = 500,
    mad_n = 5,
    max_mito = 5,
    calc_ribo = FALSE,
    max_ribo = 3,
    calc_drop = FALSE,
    max_drop = 0.95,
    log_g2u = FALSE,
    min_g2u = 0.80,
    rm_dbl = FALSE,
    dbl_method = c("class", "score"),
    dbl_score_thr = 0.5,
    method = c("MAD", "fixed", "none"),
    fixed_thr = list(
      max_feat = 6000,
      max_umi = 20000
    ),
    mito_pat = NULL,
    ribo_pat = NULL,
    sample_col = "orig.ident",
    outdir = "QC",
    save_object = TRUE,
    merge = FALSE,
    join_layers = TRUE,
    parallel = FALSE,
    n_cores = 2L,
    seed = 1234L,
    return = c("result", "object"),
    verbose = TRUE
) {
  species <- match.arg(species)
  method <- match.arg(method)
  dbl_method <- match.arg(dbl_method)
  return <- match.arg(return)

  .qc_validate_arguments(
    min_feat = min_feat,
    min_umi = min_umi,
    mad_n = mad_n,
    max_mito = max_mito,
    max_ribo = max_ribo,
    max_drop = max_drop,
    min_g2u = min_g2u,
    dbl_score_thr = dbl_score_thr,
    method = method,
    fixed_thr = fixed_thr,
    n_cores = n_cores
  )

  patterns <- .qc_gene_patterns(
    species = species,
    mito_pat = mito_pat,
    ribo_pat = ribo_pat
  )

  inputs <- .qc_normalize_inputs(
    x = x,
    sample_id = sample_id
  )

  if (!is.null(outdir)) {
    dir.create(
      outdir,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }

  .qc_message(
    verbose,
    "[QC] Starting quality control for ",
    length(inputs),
    " input(s)."
  )

  worker <- function(i) {
    input_name <- names(inputs)[i]

    object <- tryCatch(
      read_sc_input(
        x = inputs[[i]],
        sample_id = input_name,
        assay = assay,
        min_cells = 0,
        min_features = 0,
        feature_type = feature_type,
        verbose = verbose && !parallel
      ),
      error = function(e) {
        return(
          list(
            status = "failed",
            sample = input_name,
            stage = "input",
            error = conditionMessage(e)
          )
        )
      }
    )

    if (is.list(object) &&
        identical(object$status, "failed")) {
      return(object)
    }

    .qc_process_object(
      object = object,
      sample_id = input_name,
      assay = assay,
      layer = layer,
      min_feat = min_feat,
      min_umi = min_umi,
      mad_n = mad_n,
      max_mito = max_mito,
      calc_ribo = calc_ribo,
      max_ribo = max_ribo,
      calc_drop = calc_drop,
      max_drop = max_drop,
      log_g2u = log_g2u,
      min_g2u = min_g2u,
      rm_dbl = rm_dbl,
      dbl_method = dbl_method,
      dbl_score_thr = dbl_score_thr,
      method = method,
      fixed_thr = fixed_thr,
      mito_pat = patterns$mito,
      ribo_pat = patterns$ribo,
      sample_col = sample_col,
      seed = seed + i - 1L,
      verbose = verbose && !parallel
    )
  }

  if (parallel && length(inputs) > 1L) {
    results <- .qc_apply_parallel(
      indices = seq_along(inputs),
      worker = worker,
      n_cores = n_cores,
      verbose = verbose
    )
  } else {
    results <- vector("list", length(inputs))

    pb <- NULL

    if (verbose && length(inputs) > 1L) {
      pb <- utils::txtProgressBar(
        min = 0,
        max = length(inputs),
        style = 3
      )
    }

    for (i in seq_along(inputs)) {
      results[[i]] <- worker(i)

      if (!is.null(pb)) {
        utils::setTxtProgressBar(pb, i)
      }
    }

    if (!is.null(pb)) {
      close(pb)
    }
  }

  success <- vapply(
    results,
    function(z) identical(z$status, "success"),
    logical(1)
  )

  failures <- .qc_collect_failures(results[!success])

  if (!any(success)) {
    msg <- "Quality control failed for every input."

    if (nrow(failures) > 0L) {
      details <- paste0(
        failures$sample,
        " [",
        failures$stage,
        "]: ",
        failures$error,
        collapse = "\n"
      )

      msg <- paste(msg, details, sep = "\n")
    }

    stop(msg, call. = FALSE)
  }

  successful <- results[success]

  objects <- lapply(successful, `[[`, "object")

  names(objects) <- vapply(
    successful,
    `[[`,
    character(1),
    "sample"
  )

  qc_summary <- .qc_bind_rows(
    lapply(successful, `[[`, "summary")
  )

  threshold_summary <- .qc_bind_rows(
    lapply(successful, `[[`, "thresholds")
  )

  doublet_summary <- .qc_bind_rows(
    lapply(successful, `[[`, "doublets")
  )

  removed_cells <- .qc_bind_rows(
    lapply(successful, `[[`, "removed_cells")
  )

  output_object <- .qc_combine_objects(
    objects = objects,
    merge = merge,
    join_layers = join_layers,
    assay = assay,
    verbose = verbose
  )

  parameters <- list(
    assay = assay,
    layer = layer,
    feature_type = feature_type,
    species = species,
    min_feat = min_feat,
    min_umi = min_umi,
    mad_n = mad_n,
    max_mito = max_mito,
    calc_ribo = calc_ribo,
    max_ribo = max_ribo,
    calc_drop = calc_drop,
    max_drop = max_drop,
    log_g2u = log_g2u,
    min_g2u = min_g2u,
    rm_dbl = rm_dbl,
    dbl_method = dbl_method,
    dbl_score_thr = dbl_score_thr,
    method = method,
    fixed_thr = fixed_thr,
    mito_pat = patterns$mito,
    ribo_pat = patterns$ribo,
    sample_col = sample_col,
    merge = merge,
    seed = seed
  )

  files <- .qc_write_results(
    objects = objects,
    summary = qc_summary,
    thresholds = threshold_summary,
    doublets = doublet_summary,
    removed_cells = removed_cells,
    failures = failures,
    parameters = parameters,
    outdir = outdir,
    save_object = save_object
  )

  result <- structure(
    list(
      object = output_object,
      objects = objects,
      summary = qc_summary,
      thresholds = threshold_summary,
      doublets = doublet_summary,
      removed_cells = removed_cells,
      failures = failures,
      parameters = parameters,
      files = files
    ),
    class = "scqc_result"
  )

  .qc_message(
    verbose,
    "[QC] Completed. Retained ",
    sum(qc_summary$final_cells),
    " of ",
    sum(qc_summary$input_cells),
    " cells."
  )

  if (return == "object") {
    return(output_object)
  }

  result
}
