#' Compare quality-control metrics across samples
#'
#' @description
#' Compare scRNA-seq and snRNA-seq quality-control metrics across samples
#' and experimental groups. Supports Seurat objects, cellqc results, and
#' named lists of cellqc results. Generates publication-ready QC figures.
#'
#' @details
#' The function provides:
#' \itemize{
#'   \item QC distributions (violin and boxplots).
#'   \item Cell counts and QC retention per sample.
#'   \item Heatmap of standardized median QC metrics.
#'   \item UMI versus detected genes scatter plots.
#'   \item Combined QC dashboard.
#' }
#'
#' QC retention is reported only when pass/fail classifications are available.
#' No filtering or statistical testing is performed.
#'
#' @param object Seurat object, cellqc result, or named list of cellqc results.
#' @param sample.col Metadata column containing sample identities.
#' @param group.col Optional metadata column containing experimental groups.
#' @param sample.group Optional named vector mapping samples to groups.
#' @param assay Assay used to calculate missing count-based QC metrics.
#' @param species Species: "human" or "mouse".
#' @param metrics QC metrics to visualize.
#' @param color.by Color figures by "group" or "sample".
#' @param group.cols Optional named color vector for experimental groups.
#' @param sample.cols Optional named color vector for individual samples.
#' @param max.points Maximum number of cells plotted per sample.
#' @param seed Random seed for reproducible subsampling.
#' @param outdir Output directory. NULL disables file export.
#' @param verbose Print processing information.
#'
#' @return A list containing:
#' \describe{
#'   \item{qc}{Per-cell QC metadata.}
#'   \item{summary}{Per-sample QC summary.}
#'   \item{plots}{Individual ggplot objects.}
#'   \item{dashboard}{Combined patchwork dashboard.}
#'   \item{parameters}{Analysis parameters.}
#' }
#'
#' @export
qc_compare <- function(
    object,
    sample.col="sample",
    group.col="group",
    sample.group=NULL,
    assay="RNA",
    species=c("human","mouse"),
    metrics=c("nCount_RNA","nFeature_RNA","percent.mt","complexity"),
    color.by=c("group","sample"),
    group.cols=NULL,
    sample.cols=NULL,
    max.points=2000L,
    seed=1234L,
    outdir="QC/comparison",
    verbose=TRUE
) {

  # ------------------------------------------------------------------------- #
  # Setup
  # ------------------------------------------------------------------------- #
  species <- match.arg(species)
  color.by <- match.arg(color.by)

  for(pkg in c("ggplot2","dplyr","tidyr","patchwork"))
    if(!requireNamespace(pkg,quietly=TRUE))
      stop("Required package: ",pkg,call.=FALSE)

  if(length(max.points)!=1L || !is.finite(max.points) || max.points<1)
    stop("`max.points` must be a positive integer.",call.=FALSE)

  is_cellqc <- function(x)
    is.list(x) && is.data.frame(x$qc) &&
    all(c("nCount_RNA","nFeature_RNA") %in% names(x$qc))

  # ------------------------------------------------------------------------- #
  # Input — Seurat or cellqc results
  # ------------------------------------------------------------------------- #
  if(inherits(object,"Seurat")) {

    md <- object[[]]

    if(!sample.col %in% names(md))
      stop("Missing sample column: ",sample.col,call.=FALSE)

    qc <- data.frame(
      cell=rownames(md),
      sample=as.character(md[[sample.col]]),
      stringsAsFactors=FALSE
    )

    if(!is.null(group.col) && group.col %in% names(md))
      qc$group <- as.character(md[[group.col]])

    cols <- intersect(
      c(metrics,"nCount_RNA","nFeature_RNA","percent.mt",
        "percent.ribo","percent.hb","complexity",
        "pass","status","cellqc_status"),
      names(md)
    )

    for(nm in cols) qc[[nm]] <- md[[nm]]

    if(!all(c("nCount_RNA","nFeature_RNA") %in% names(qc)))
      stop("Missing nCount_RNA or nFeature_RNA.",call.=FALSE)

    if("percent.mt" %in% metrics && !"percent.mt" %in% names(qc)) {

      if(!requireNamespace("SeuratObject",quietly=TRUE) ||
         !requireNamespace("Matrix",quietly=TRUE))
        stop("SeuratObject and Matrix are required.",call.=FALSE)

      mat <- SeuratObject::LayerData(
        object,assay=assay,layer="counts"
      )

      mt <- grepl(
        if(species=="human") "^MT-" else "^mt-",
        rownames(mat),
        ignore.case=species=="mouse"
      )

      if(any(mt)) {
        totals <- Matrix::colSums(mat)
        values <- 100*Matrix::colSums(mat[mt,,drop=FALSE])/pmax(totals,1)
        qc$percent.mt <- unname(values[qc$cell])
      } else {
        warning("No mitochondrial genes detected.")
      }
    }

  } else if(is_cellqc(object)) {

    qc <- object$qc

  } else if(is.list(object) && length(object)>0L &&
            all(vapply(object,is_cellqc,logical(1)))) {

    ids <- names(object)

    if(is.null(ids) || anyNA(ids) || any(!nzchar(ids)))
      ids <- paste0("Sample",seq_along(object))

    qc <- dplyr::bind_rows(lapply(seq_along(object),function(i) {
      d <- object[[i]]$qc
      d$sample <- ids[i]
      d
    }))

  } else {
    stop("Provide a Seurat object or cellqc result(s).",call.=FALSE)
  }

  # ------------------------------------------------------------------------- #
  # Metadata
  # ------------------------------------------------------------------------- #
  if(!"cell" %in% names(qc))
    qc$cell <- seq_len(nrow(qc))

  if(!"sample" %in% names(qc))
    stop("Missing sample column.",call.=FALSE)

  qc$sample <- as.character(qc$sample)
  qc <- qc[!is.na(qc$sample) & nzchar(qc$sample),,drop=FALSE]

  if(!nrow(qc))
    stop("No QC records available.",call.=FALSE)

  if(!is.null(sample.group)) {

    if(is.null(names(sample.group)))
      stop("`sample.group` must be named.",call.=FALSE)

    qc$group <- unname(as.character(sample.group[qc$sample]))
  }

  if(!"group" %in% names(qc))
    qc$group <- "All"

  qc$group <- as.character(qc$group)
  qc$group[is.na(qc$group) | !nzchar(qc$group)] <- "Unknown"

  if(!"complexity" %in% names(qc))
    qc$complexity <- qc$nFeature_RNA/pmax(qc$nCount_RNA,1)

  if(!"pass" %in% names(qc)) {

    status <- if("status" %in% names(qc)) qc$status else
      if("cellqc_status" %in% names(qc)) qc$cellqc_status else NULL

    qc$pass <- if(is.null(status)) rep(NA,nrow(qc)) else
      ifelse(as.character(status)=="Pass",TRUE,
             ifelse(as.character(status)=="Flagged",FALSE,NA))
  }

  qc$pass <- as.logical(qc$pass)

  metrics <- unique(intersect(metrics,names(qc)))
  metrics <- metrics[vapply(qc[metrics],is.numeric,logical(1))]

  if(!length(metrics))
    stop("No numeric QC metrics available.",call.=FALSE)

  qc$sample <- factor(qc$sample,levels=unique(qc$sample))
  qc$group <- factor(qc$group,levels=unique(qc$group))

  # ------------------------------------------------------------------------- #
  # Colors — group or sample
  # ------------------------------------------------------------------------- #
  make_cols <- function(x,cols=NULL) {

    defaults <- setNames(
      grDevices::hcl.colors(length(x),"Dark 3"),
      x
    )

    if(is.null(cols)) return(defaults)

    if(is.null(names(cols))) {
      if(length(cols)!=length(x))
        stop("Unnamed palette must match category count.",call.=FALSE)
      return(setNames(cols,x))
    }

    missing <- setdiff(x,names(cols))

    if(length(missing))
      stop("Missing colors: ",paste(missing,collapse=", "),call.=FALSE)

    cols[x]
  }

  group.cols <- make_cols(levels(qc$group),group.cols)
  sample.cols <- make_cols(levels(qc$sample),sample.cols)

  color.col <- if(color.by=="group") "group" else "sample"
  color.values <- if(color.by=="group") group.cols else sample.cols
  color.title <- if(color.by=="group") "Group" else "Sample"

  fill.scale <- ggplot2::scale_fill_manual(
    values=color.values,name=color.title,drop=FALSE
  )

  color.scale <- ggplot2::scale_colour_manual(
    values=color.values,name=color.title,drop=FALSE
  )

  # ------------------------------------------------------------------------- #
  # Summary
  # ------------------------------------------------------------------------- #
  med <- function(x)
    if(all(is.na(x))) NA_real_ else stats::median(x,na.rm=TRUE)

  summary <- dplyr::summarise(
    dplyr::group_by(qc,.data$sample,.data$group),
    cells=dplyr::n(),
    median_UMI=med(.data$nCount_RNA),
    median_genes=med(.data$nFeature_RNA),
    median_mt=if("percent.mt" %in% names(qc))
      med(.data$percent.mt) else NA_real_,
    median_complexity=med(.data$complexity),
    pct_pass=if(all(is.na(.data$pass))) NA_real_
    else 100*mean(.data$pass,na.rm=TRUE),
    .groups="drop"
  )

  long <- tidyr::pivot_longer(
    qc,
    cols=dplyr::all_of(metrics),
    names_to="metric",
    values_to="value"
  )

  long$metric <- factor(long$metric,levels=metrics)

  # ------------------------------------------------------------------------- #
  # Plot theme
  # ------------------------------------------------------------------------- #
  th <- ggplot2::theme_bw(base_size=11)+
    ggplot2::theme(
      strip.text=ggplot2::element_text(face="bold"),
      panel.grid.minor=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_text(angle=45,hjust=1),
      legend.position="bottom"
    )

  # ------------------------------------------------------------------------- #
  # 1. QC distributions
  # ------------------------------------------------------------------------- #
  p_dist <- ggplot2::ggplot(
    long,
    ggplot2::aes(
      x=.data$sample,
      y=.data$value,
      fill=.data[[color.col]]
    )
  )+
    ggplot2::geom_violin(
      scale="width",trim=TRUE,alpha=.65,na.rm=TRUE
    )+
    ggplot2::geom_boxplot(
      width=.13,outlier.shape=NA,fill="white",alpha=.75,na.rm=TRUE
    )+
    ggplot2::facet_wrap(~metric,scales="free_y",ncol=2)+
    fill.scale+
    ggplot2::labs(
      title="QC distributions by sample",
      x=NULL,y=NULL
    )+th

  # ------------------------------------------------------------------------- #
  # 2. Cell counts
  # ------------------------------------------------------------------------- #
  p_cells <- ggplot2::ggplot(
    summary,
    ggplot2::aes(
      x=.data$sample,
      y=.data$cells,
      fill=.data[[color.col]]
    )
  )+
    ggplot2::geom_col(width=.7)+
    fill.scale+
    ggplot2::labs(
      title="Barcodes per sample",
      x=NULL,y="Barcodes"
    )+th

  # ------------------------------------------------------------------------- #
  # 3. QC retention
  # ------------------------------------------------------------------------- #
  p_ret <- NULL

  if(any(!is.na(summary$pct_pass))) {

    p_ret <- ggplot2::ggplot(
      summary,
      ggplot2::aes(
        x=.data$sample,
        y=.data$pct_pass,
        fill=.data[[color.col]]
      )
    )+
      ggplot2::geom_col(width=.7,na.rm=TRUE)+
      ggplot2::coord_cartesian(ylim=c(0,100))+
      fill.scale+
      ggplot2::labs(
        title="QC retention",
        x=NULL,y="Passing barcodes (%)"
      )+th
  }

  # ------------------------------------------------------------------------- #
  # 4. QC heatmap
  # ------------------------------------------------------------------------- #
  hm <- dplyr::summarise(
    dplyr::group_by(long,.data$sample,.data$metric),
    value=med(.data$value),
    .groups="drop"
  )

  hm <- dplyr::mutate(
    dplyr::group_by(hm,.data$metric),
    z={
      v <- .data$value
      s <- stats::sd(v,na.rm=TRUE)
      if(!is.finite(s) || s==0)
        rep(0,length(v))
      else
        (v-mean(v,na.rm=TRUE))/s
    }
  ) |>
    dplyr::ungroup()

  p_heat <- ggplot2::ggplot(
    hm,
    ggplot2::aes(
      x=.data$metric,
      y=.data$sample,
      fill=.data$z
    )
  )+
    ggplot2::geom_tile(colour="white")+
    ggplot2::scale_fill_gradient2(
      low="#377EB8",
      mid="white",
      high="#E41A1C",
      midpoint=0,
      na.value="grey90"
    )+
    ggplot2::labs(
      title="Median QC metrics (Z-score)",
      x=NULL,y=NULL,fill="Z-score"
    )+th

  # ------------------------------------------------------------------------- #
  # 5. UMI vs detected genes
  # ------------------------------------------------------------------------- #
  set.seed(seed)

  idx <- unlist(
    lapply(split(seq_len(nrow(qc)),qc$sample),function(i)
      if(length(i)>max.points)
        sample(i,max.points)
      else i
    ),
    use.names=FALSE
  )

  d <- qc[idx,,drop=FALSE]
  d <- d[d$nCount_RNA>0 & d$nFeature_RNA>0,,drop=FALSE]

  p_scatter <- ggplot2::ggplot(
    d,
    ggplot2::aes(
      x=.data$nCount_RNA,
      y=.data$nFeature_RNA,
      colour=.data[[color.col]]
    )
  )+
    ggplot2::geom_point(size=.5,alpha=.5)+
    ggplot2::scale_x_log10()+
    ggplot2::scale_y_log10()+
    ggplot2::facet_wrap(~sample)+
    color.scale+
    ggplot2::labs(
      title="UMI vs detected genes",
      x="UMI (log10)",
      y="Genes (log10)"
    )+th

  # ------------------------------------------------------------------------- #
  # Combine figures
  # ------------------------------------------------------------------------- #
  plots <- list(
    distributions=p_dist,
    cells=p_cells,
    retention=p_ret,
    heatmap=p_heat,
    scatter=p_scatter
  )

  plots <- Filter(Negate(is.null),plots)

  dashboard <- patchwork::wrap_plots(plots,ncol=2)+
    patchwork::plot_annotation(
      title="scQC — sample quality comparison"
    )

  # ------------------------------------------------------------------------- #
  # Parameters
  # ------------------------------------------------------------------------- #
  params <- list(
    sample.col=sample.col,
    group.col=group.col,
    metrics=metrics,
    species=species,
    color.by=color.by,
    group.cols=group.cols,
    sample.cols=sample.cols,
    max.points=max.points,
    seed=seed
  )

  # ------------------------------------------------------------------------- #
  # Export
  # ------------------------------------------------------------------------- #
  if(!is.null(outdir)) {

    dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

    utils::write.csv(
      qc,file.path(outdir,"qc_compare_cells.csv"),row.names=FALSE
    )

    utils::write.csv(
      summary,file.path(outdir,"qc_compare_summary.csv"),row.names=FALSE
    )

    ggplot2::ggsave(
      file.path(outdir,"qc_compare_dashboard.pdf"),
      dashboard,width=14,height=12
    )

    for(nm in names(plots))
      ggplot2::ggsave(
        file.path(outdir,paste0("qc_compare_",nm,".pdf")),
        plots[[nm]],width=10,height=7
      )

    saveRDS(
      list(summary=summary,parameters=params),
      file.path(outdir,"qc_compare_results.rds")
    )
  }

  if(verbose)
    message(
      "[QC_COMPARE] ",nrow(summary),
      " samples | ",nrow(qc)," barcodes"
    )

  structure(
    list(
      qc=qc,
      summary=summary,
      plots=plots,
      dashboard=dashboard,
      parameters=params
    ),
    class=c("qc_compare_result","list")
  )
}
