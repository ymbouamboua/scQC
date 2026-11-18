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
#' @param seed Random seed used by scDblFinder and ambient correction.
#' @param return Return `"result"` or only the processed `"object"`.
#' @param verbose Print progress messages.
#'
#' @param decontaminate Enable ambient RNA correction after raw-count QC and
#'   doublet removal. Default `FALSE` preserves the existing workflow.
#' @param ambient_method Correction method: `"decontX"` (default) or `"SoupX"`.
#'   The standalone decontX package is preferred; compatible celda is supported.
#' @param ambient_raw Optional unfiltered droplet matrix or 10x path. For several
#'   inputs, use a list named like the inputs. For several libraries inside one
#'   object, use a nested list named by library. Required for SoupX, optional for
#'   decontX. All original called-cell barcodes are excluded from background.
#' @param ambient_batch Optional metadata column identifying capture libraries
#'   or channels, not biological donors pooled within the same channel. Without
#'   it, each input is treated as one library. Correction is separate per library.
#' @param ambient_clusters Optional metadata column of preliminary cell clusters.
#'   Required for SoupX; decontX can infer clusters when omitted.
#' @param ambient_args Named backend arguments. For decontX these are additional
#'   decontX arguments. For SoupX they are autoEstCont arguments, or
#'   `list(contamination_fraction=0.1)` for an explicitly supplied fraction.
#' @param keep_raw Keep raw and corrected assays in the returned object. Default
#'   `FALSE` keeps corrected counts only in the source assay. Keep raw inputs on
#'   disk before assigning the result over your only copy. The supplied input
#'   object is never modified in place. `TRUE` retains raw counts and adds the
#'   corrected assay, increasing memory and RDS size.
#' @param ambient_assay Name for the corrected assay when `keep_raw=TRUE`.
#'   Default `"RNAdecont"`; existing assays are not overwritten.
#' @details QC thresholds and doublet detection use raw counts. Corrected counts
#'   are not rounded and can be fractional. The corrected assay becomes default.
#'   In compact mode normalized/scaled data and cached reductions, graphs,
#'   neighbors and commands are removed because they predate correction. Rerun
#'   normalization and downstream analyses. Raw QC metrics are kept with a
#'   `qc_raw_` prefix, and corrected QC metrics are refreshed without refiltering.
#'   Per-cell contamination and before/after summaries are stored in
#'   `misc$scqc_ambient`, returned as `ambient` / `ambient_cells`, and exported
#'   when `outdir` is set. Merged results retain per-library provenance.
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
    verbose = TRUE,
    decontaminate = FALSE,
    ambient_method = c("decontX", "SoupX"),
    ambient_raw = NULL,
    ambient_batch = NULL,
    ambient_clusters = NULL,
    ambient_args = list(),
    keep_raw = FALSE,
    ambient_assay = "RNA"
) {
  species <- match.arg(species)
  method <- match.arg(method)
  dbl_method <- match.arg(dbl_method)
  return <- match.arg(return)
  ambient_method <- match.arg(ambient_method)
  .qc_ambient_validate(decontaminate, ambient_method, keep_raw,
    ambient_assay, assay, ambient_batch, ambient_clusters, ambient_args)

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

  raw_inputs <- if (decontaminate) .qc_ambient_input_raw(ambient_raw,names(inputs)) else
    stats::setNames(rep(list(NULL),length(inputs)),names(inputs))

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
      verbose = verbose && !parallel,
      decontaminate = decontaminate,
      ambient_method = ambient_method,
      ambient_raw = raw_inputs[[input_name]],
      ambient_batch = ambient_batch,
      ambient_clusters = ambient_clusters,
      ambient_args = ambient_args,
      keep_raw = keep_raw,
      ambient_assay = ambient_assay,
      feature_type = feature_type
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

  ambient_summary <- .qc_bind_rows(lapply(successful, `[[`, "ambient"))
  ambient_cells <- .qc_bind_rows(lapply(successful, `[[`, "ambient_cells"))

  output_object <- .qc_combine_objects(
    objects = objects,
    merge = merge,
    join_layers = join_layers,
    assay = assay,
    verbose = verbose
  )

  if (decontaminate && merge && length(objects)>1L) {
    if (join_layers && keep_raw) output_object <- SeuratObject::JoinLayers(output_object,assay=ambient_assay)
    output_object@misc$scqc_ambient <- list(
      method=ambient_method,keep_raw=keep_raw,
      output_assay=if (keep_raw) ambient_assay else assay,
      summary=ambient_summary,
      inputs=stats::setNames(lapply(objects,function(z) z@misc$scqc_ambient),names(objects)),
      cell_names="Per-input original barcodes; merged IDs are recorded in merged_cell.")
    if (nrow(ambient_cells)) {
      ambient_cells$merged_cell <- paste0(ambient_cells$sample,"_",ambient_cells$cell)
      output_object@misc$scqc_ambient$cells <- ambient_cells
    }
  }

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
    seed = seed,
    decontaminate = decontaminate, ambient_method = ambient_method,
    ambient_batch = ambient_batch, ambient_clusters = ambient_clusters,
    ambient_args = ambient_args, ambient_raw_supplied = !is.null(ambient_raw),
    keep_raw = keep_raw, ambient_assay = ambient_assay
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

  if (decontaminate && !is.null(outdir)) {
    files$ambient <- file.path(outdir,"Ambient_summary.tsv")
    files$ambient_cells <- file.path(outdir,"Ambient_cells.tsv")
    utils::write.table(ambient_summary,files$ambient,sep="\t",quote=FALSE,row.names=FALSE)
    utils::write.table(ambient_cells,files$ambient_cells,sep="\t",quote=FALSE,row.names=FALSE)
    files$ambient_umap <- tryCatch(
      .qc_save_ambient_umaps(ambient_cells,outdir),
      error=function(e) {
        warning("Could not export DecontX UMAPs: ",conditionMessage(e),call.=FALSE)
        character()
      }
    )
  }

  result <- structure(
    list(
      object = output_object,
      objects = objects,
      summary = qc_summary,
      thresholds = threshold_summary,
      doublets = doublet_summary,
      ambient = ambient_summary,
      ambient_cells = ambient_cells,
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
