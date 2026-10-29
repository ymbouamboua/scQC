# Ambient correction is opt-in. Library boundaries must never be pooled.
.qc_ambient_validate <- function(decontaminate, ambient_method, keep_raw,
                                ambient_assay, assay, ambient_batch, ambient_clusters, ambient_args) {
  for (nm in c("decontaminate","keep_raw")) {
    value <- get(nm)
    if (!is.logical(value) || length(value)!=1L || is.na(value))
      stop("`",nm,"` must be TRUE or FALSE.",call.=FALSE)
  }
  if (!decontaminate) return(invisible(NULL))
  string <- function(x) is.character(x) && length(x)==1L && !is.na(x) && nzchar(x)
  if (!string(ambient_assay) || (keep_raw && ambient_assay==assay))
    stop("Choose a distinct `ambient_assay` when keeping raw counts.",call.=FALSE)
  for (nm in c("ambient_batch","ambient_clusters")) {
    value <- get(nm)
    if (!is.null(value) && !string(value)) stop("`",nm,"` must be a metadata column name.",call.=FALSE)
  }
  if (!is.list(ambient_args) || (length(ambient_args) && (is.null(names(ambient_args)) ||
      anyNA(names(ambient_args)) || any(!nzchar(names(ambient_args))) || anyDuplicated(names(ambient_args)))))
    stop("`ambient_args` must be a uniquely named list.",call.=FALSE)
  reserved <- c("x","z","batch","background","bgBatch","assayName","bgAssayName","seed","verbose","sc","doPlot")
  if (any(names(ambient_args) %in% reserved))
    stop("Use the dedicated run_qc arguments instead of reserved `ambient_args`: ",
         paste(intersect(names(ambient_args),reserved),collapse=", "),call.=FALSE)
  if (ambient_method=="SoupX" && is.null(ambient_clusters))
    stop("SoupX requires `ambient_clusters` naming a metadata clustering column.",call.=FALSE)
  invisible(NULL)
}

.qc_ambient_counts <- function(x) {
  x <- .validate_count_matrix(x)
  if (any(!is.finite(x@x)) || any(abs(x@x-round(x@x))>1e-8))
    stop("Ambient correction requires finite, non-negative raw integer counts.",call.=FALSE)
  if (anyNA(rownames(x)) || anyNA(colnames(x)) || any(!nzchar(rownames(x))) || any(!nzchar(colnames(x))))
    stop("Ambient matrices require nonempty feature and barcode names.",call.=FALSE)
  x
}

.qc_ambient_backend <- function(method) {
  if (method=="SoupX") {
    .qc_require_package("SoupX")
    return("SoupX")
  }
  if (requireNamespace("decontX",quietly=TRUE)) return("decontX")
  if (requireNamespace("celda",quietly=TRUE)) return("celda")
  stop("Install Bioconductor package 'decontX' (or a compatible 'celda') for ambient correction.",call.=FALSE)
}

.qc_read_ambient_raw <- function(x,feature_type) {
  if (.is_count_matrix(x)) return(.qc_ambient_counts(x))
  if (!is.character(x) || length(x)!=1L || is.na(x) || !file.exists(path.expand(x)))
    stop("`ambient_raw` must be a raw-droplet matrix or an existing 10x path.",call.=FALSE)
  x <- path.expand(x)
  if (dir.exists(x)) {
    choices <- c(file.path(x,"raw_feature_bc_matrix.h5"),file.path(x,"raw_feature_bc_matrix"),
                 file.path(x,"outs","raw_feature_bc_matrix.h5"),file.path(x,"outs","raw_feature_bc_matrix"))
    if (any(file.exists(choices))) x <- choices[file.exists(choices)][1L]
  }
  resolved <- .resolve_10x_path(x)
  if (grepl("filtered_feature",resolved$path,fixed=TRUE))
    stop("Supply unfiltered raw droplets, not a filtered_feature matrix.",call.=FALSE)
  raw <- if (resolved$type=="10x_h5") Seurat::Read10X_h5(resolved$path,use.names=TRUE,unique.features=TRUE) else
    Seurat::Read10X(resolved$path,gene.column=2,unique.features=TRUE,strip.suffix=FALSE)
  .qc_ambient_counts(.select_feature_type(raw,feature_type))
}

.qc_ambient_input_raw <- function(raw,input_names) {
  if (is.null(raw)) return(stats::setNames(rep(list(NULL),length(input_names)),input_names))
  if (is.list(raw) && !is.null(names(raw)) && !anyDuplicated(names(raw)) && setequal(names(raw),input_names)) return(raw[input_names])
  if (length(input_names)==1L) return(stats::setNames(list(raw),input_names))
  if (!is.list(raw) || is.null(names(raw)) || anyDuplicated(names(raw)) ||
      !setequal(names(raw),input_names))
    stop("For multiple inputs, `ambient_raw` must be a list named exactly like the inputs.",call.=FALSE)
  raw[input_names]
}

.qc_ambient_groups <- function(object,column) {
  cells <- colnames(object)
  if (is.null(column)) return(stats::setNames(rep("library",length(cells)),cells))
  if (!column %in% names(object[[]])) stop("Ambient library/cluster column not found: ",column,call.=FALSE)
  groups <- as.character(object[[]][[column]])
  if (anyNA(groups) || any(!nzchar(groups))) stop("Ambient library/cluster IDs must not be missing or empty.",call.=FALSE)
  stats::setNames(groups,cells)
}

.qc_ambient_fit <- function(counts,background,clusters,method,backend,args,seed,verbose) {
  if (method=="decontX") {
    fun <- getExportedValue(backend,"decontX")
    fit <- do.call(fun,c(list(x=counts,z=if (is.null(clusters)) NULL else as.integer(factor(clusters)),
      background=background,seed=seed,verbose=verbose),args))
    return(list(counts=fit$decontXcounts,contamination=fit$contamination))
  }
  if (is.null(background)) stop("SoupX requires unfiltered `ambient_raw` droplets.",call.=FALSE)
  # Keep original raw droplet counts for retained cells; empty droplets are
  # supplied separately after removing every called-cell barcode.
  tod <- cbind(background,counts)
  sc <- SoupX::SoupChannel(tod=tod,toc=counts,calcSoupProfile=TRUE)
  sc <- SoupX::setClusters(sc,stats::setNames(clusters,colnames(counts)))
  rho <- args$contamination_fraction
  args$contamination_fraction <- NULL
  if (!is.null(rho)) {
    if (!is.numeric(rho) || length(rho)!=1L || !is.finite(rho) || rho<0 || rho>=1)
      stop("SoupX contamination_fraction must be in [0, 1).",call.=FALSE)
    if (length(args)) stop("SoupX auto-estimation arguments cannot accompany contamination_fraction.",call.=FALSE)
    sc <- SoupX::setContaminationFraction(sc,rho,forceAccept=FALSE)
  } else sc <- do.call(SoupX::autoEstCont,c(list(sc=sc,doPlot=FALSE,verbose=verbose),args))
  corrected <- SoupX::adjustCounts(sc,roundToInt=FALSE,verbose=as.integer(verbose))
  list(counts=corrected,contamination=sc$metaData[colnames(counts),"rho"])
}

.qc_decontaminate <- function(object,assay,layer,sample_id,method,raw,batches,clusters,
                             called_cells,keep_raw,ambient_assay,args,feature_type,seed,verbose) {
  if (identical(object@misc$scqc_ambient$output_assay,assay))
    stop("This assay is already ambient-corrected. Start from raw counts to avoid correcting twice.",call.=FALSE)
  backend <- .qc_ambient_backend(method)
  # Reject split or ambiguous layers rather than accidentally correcting only
  # the first library. Batch boundaries can be specified explicitly after joining.
  layers <- SeuratObject::Layers(object[[assay]],search=NA)
  selected <- layers[layers==layer | startsWith(layers,paste0(layer,"."))]
  if (length(selected)!=1L) stop("Select one joined raw-count layer for ambient correction; process libraries separately or supply ambient_batch.",call.=FALSE)
  counts <- .qc_ambient_counts(.qc_get_layer(object,assay,selected))
  if (!identical(colnames(counts),colnames(object))) stop("Selected counts must cover all current cells in object order.",call.=FALSE)
  if (keep_raw && ambient_assay %in% names(object@assays))
    stop("Ambient output assay already exists: ",ambient_assay,call.=FALSE)
  b <- batches[colnames(counts)]
  if (anyNA(b)) stop("Cannot align library labels to the retained cells.",call.=FALSE)
  libs <- unique(b)
  raw_by_library <- if (length(libs)==1L) {
    if (is.list(raw) && identical(names(raw), unname(libs))) raw else list(raw)
  } else {
    if (is.null(raw)) rep(list(NULL),length(libs)) else {
      if (!is.list(raw) || is.null(names(raw)) || anyDuplicated(names(raw)) || !setequal(names(raw),libs))
        stop("For several ambient libraries, supply `ambient_raw` as a list named by library.",call.=FALSE)
      raw[libs]
    }
  }
  names(raw_by_library) <- libs
  # Preserve the global RNG state even for stochastic backend algorithms.
  had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if (had_seed) old_seed <- get(".Random.seed",envir=.GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed",old_seed,envir=.GlobalEnv) else
    if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv))
  corrected_parts <- diagnostics <- summaries <- list()
  for (i in seq_along(libs)) {
    lib <- libs[i];cells <- colnames(counts)[b==lib]
    mat <- counts[,cells,drop=FALSE]
    if (any(Matrix::colSums(mat)==0)) stop("Zero-count cells cannot be decontaminated; increase the raw QC thresholds.",call.=FALSE)
    bg <- if (is.null(raw_by_library[[lib]])) NULL else .qc_read_ambient_raw(raw_by_library[[lib]],feature_type)
    if (!is.null(bg)) {
      if (!all(rownames(mat) %in% rownames(bg))) stop("Raw droplets lack features present in the filtered matrix; use matching feature names.",call.=FALSE)
      bg <- bg[rownames(mat),!colnames(bg) %in% called_cells,drop=FALSE]
      bg <- bg[,Matrix::colSums(bg)>0,drop=FALSE]
      if (!ncol(bg)) stop("No nonempty background droplets remain after excluding all called cells.",call.=FALSE)
    }
    z <- if (is.null(clusters)) NULL else clusters[cells]
    if (!is.null(z) && (anyNA(z) || length(unique(z))<2L))
      stop("Ambient clusters must cover the cells and contain at least two groups per library.",call.=FALSE)
    set.seed(seed+i-1L)
    .qc_message(verbose,"[AMBIENT] ",sample_id," / ",lib,": ",method)
    fit <- .qc_ambient_fit(mat,bg,z,method,backend,args,seed+i-1L,verbose)
    clean <- fit$counts
    if (is.null(clean) || !identical(dim(clean),dim(mat)) || is.null(rownames(clean)) || is.null(colnames(clean)) ||
        !setequal(rownames(clean),rownames(mat)) || !setequal(colnames(clean),cells))
      stop("Ambient backend returned incompatible count dimensions or cell/feature names.",call.=FALSE)
    clean <- methods::as(clean[rownames(mat),cells,drop=FALSE],"dgCMatrix")
    if (any(!is.finite(clean@x)) || any(clean@x<0)) stop("Ambient backend returned invalid corrected counts.",call.=FALSE)
    rho <- as.numeric(fit$contamination)
    if (length(rho)!=length(cells) || any(!is.finite(rho)) || any(rho<0 | rho>1))
      stop("Ambient backend returned invalid contamination estimates.",call.=FALSE)
    before <- Matrix::colSums(mat);after <- Matrix::colSums(clean)
    if (any(after>before+1e-6*pmax(1,before))) stop("Ambient correction unexpectedly increased cell totals.",call.=FALSE)
    diagnostics[[lib]] <- data.frame(sample=sample_id,library=lib,cell=cells,method=method,
      contamination=rho,raw_umi=as.numeric(before),corrected_umi=as.numeric(after),
      removed_fraction=as.numeric(1-after/before),raw_features=as.integer(Matrix::colSums(mat>0)),
      corrected_features=as.integer(Matrix::colSums(clean>0)),stringsAsFactors=FALSE)
    summaries[[lib]] <- data.frame(sample=sample_id,library=lib,method=method,backend=backend,
      backend_version=as.character(utils::packageVersion(backend)),n_cells=length(cells),
      n_background=if (is.null(bg)) 0L else ncol(bg),
      raw_umi=sum(before),corrected_umi=sum(after),removed_fraction=1-sum(after)/sum(before),
      median_contamination=stats::median(rho),stringsAsFactors=FALSE)
    corrected_parts[[lib]] <- clean
  }
  clean <- do.call(cbind,corrected_parts)[,colnames(counts),drop=FALSE]
  diagnostic <- do.call(rbind,diagnostics);diagnostic <- diagnostic[match(colnames(object),diagnostic$cell),,drop=FALSE]
  rownames(diagnostic) <- NULL
  original_md <- object[[]]
  raw_fields <- intersect(c(paste0("nCount_",assay),paste0("nFeature_",assay),"percent_mito","percent_ribo","dropout","log10_genes_per_umi"),names(original_md))
  for (nm in raw_fields) object[[paste0("qc_raw_",nm)]] <- original_md[[nm]]
  target <- if (keep_raw) ambient_assay else assay
  replacement <- if (inherits(object[[assay]],"Assay5")) SeuratObject::CreateAssay5Object(counts=clean) else SeuratObject::CreateAssayObject(counts=clean)
  # Retain gene annotations, but not normalized data, scaled data or variable features.
  fm <- object[[assay]][[]]
  if (ncol(fm)) replacement <- SeuratObject::AddMetaData(replacement,metadata=fm[rownames(clean),,drop=FALSE])
  object[[target]] <- replacement
  invalidated <- list(reductions=character(),graphs=character(),neighbors=character())
  if (!keep_raw) {
    invalidated <- list(reductions=names(object@reductions),graphs=names(object@graphs),neighbors=names(object@neighbors))
    object@reductions <- list();object@graphs <- list();object@neighbors <- list();object@commands <- list()
  }
  object$ambient_contamination <- diagnostic$contamination
  SeuratObject::DefaultAssay(object) <- target
  # Corrected metadata is refreshed by the caller using the same QC patterns.
  object@misc$scqc_ambient <- list(method=method,backend=backend,
    input_assay=assay,input_layer=selected,output_assay=target,keep_raw=keep_raw,
    library_ids=libs,arguments=args,seed=seed,summary=do.call(rbind,summaries),
    cells=diagnostic,invalidated=invalidated,
    note="QC and doublet filtering used raw counts. Corrected values may be fractional; normalize/recompute downstream analyses.")
  list(object=object,summary=do.call(rbind,summaries),cells=diagnostic)
}
