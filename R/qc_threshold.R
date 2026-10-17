#' Calculate and visualize sample-specific QC thresholds
#'
#' @description
#' Estimate quality-control thresholds independently for each sample using
#' fixed cutoffs, median absolute deviation (MAD), or empirical quantiles.
#' Supports Seurat objects and cellqc results.
#'
#' @details
#' For MAD-based thresholds, UMI and gene counts are transformed using
#' log10(x + 1). Lower thresholds identify low-quality cells, whereas upper
#' thresholds identify unusually high-count cells requiring further review.
#'
#' The function is diagnostic only and does not filter the input object.
#' Mitochondrial thresholds can be capped using `max.mt`.
#'
#' @param object Seurat object, cellqc result, or named list of cellqc results.
#' @param sample.col Metadata column containing sample identities.
#' @param group.col Optional metadata column containing experimental groups.
#' @param sample.group Optional named vector mapping samples to groups.
#' @param assay Assay used to retrieve count-based QC metrics.
#' @param species Species: "human" or "mouse".
#' @param type Data type: "snRNA" or "scRNA".
#' @param method Threshold method: "MAD", "fixed", or "quantile".
#' @param metrics QC metrics used for threshold estimation.
#' @param mad.n Number of MADs used for adaptive thresholds.
#' @param quantiles Lower and upper quantiles for adaptive thresholds.
#' @param min.features Minimum number of detected genes.
#' @param min.counts Minimum number of UMI counts.
#' @param max.features Optional maximum number of detected genes.
#' @param max.counts Optional maximum number of UMI counts.
#' @param max.mt Maximum mitochondrial percentage; NULL disables the cap.
#' @param cap.fixed Apply fixed limits to adaptive thresholds.
#' @param upper.counts Evaluate upper thresholds for UMI and genes.
#' @param color.by Color figures by "group" or "sample".
#' @param group.cols Optional named group colors.
#' @param sample.cols Optional named sample colors.
#' @param max.points Maximum plotted cells per sample.
#' @param seed Random seed for reproducible subsampling.
#' @param outdir Output directory; NULL disables exports.
#' @param verbose Print processing information.
#'
#' @return List containing:
#' \describe{
#'   \item{qc}{Per-cell QC metrics and classifications.}
#'   \item{thresholds}{Sample-specific thresholds.}
#'   \item{summary}{QC retention per sample.}
#'   \item{plots}{Individual ggplot objects.}
#'   \item{dashboard}{Combined diagnostic dashboard.}
#'   \item{parameters}{Analysis parameters.}
#' }
#'
#' @export
qc_threshold <- function(
    object,
    sample.col="sample",
    group.col="group",
    sample.group=NULL,
    assay="RNA",
    species=c("human","mouse"),
    type=c("snRNA","scRNA"),
    method=c("MAD","fixed","quantile"),
    metrics=c("nCount_RNA","nFeature_RNA","percent.mt"),
    mad.n=5,
    quantiles=c(.01,.99),
    min.features=200,
    min.counts=500,
    max.features=NULL,
    max.counts=NULL,
    max.mt=5,
    cap.fixed=TRUE,
    upper.counts=TRUE,
    color.by=c("group","sample"),
    group.cols=NULL,
    sample.cols=NULL,
    max.points=2000L,
    seed=1234L,
    outdir="QC/threshold",
    verbose=TRUE
) {
  
  # ------------------------------------------------------------------------- #
  # Setup
  # ------------------------------------------------------------------------- #
  species <- match.arg(species)
  type <- match.arg(type)
  method <- match.arg(method)
  color.by <- match.arg(color.by)
  
  for(pkg in c("ggplot2","dplyr","tidyr","patchwork"))
    if(!requireNamespace(pkg,quietly=TRUE))
      stop("Required package: ",pkg,call.=FALSE)
  
  if(length(mad.n)!=1L || !is.finite(mad.n) || mad.n<=0)
    stop("`mad.n` must be positive.",call.=FALSE)
  
  if(length(quantiles)!=2L || anyNA(quantiles) ||
     any(quantiles<=0 | quantiles>=1) || quantiles[1]>=quantiles[2])
    stop("`quantiles` must contain two increasing probabilities in (0,1).",
         call.=FALSE)
  
  if(length(max.points)!=1L || !is.finite(max.points) || max.points<1)
    stop("`max.points` must be positive.",call.=FALSE)
  
  if(!is.null(max.mt) && (length(max.mt)!=1L ||
                          !is.finite(max.mt) || max.mt<0))
    stop("`max.mt` must be NULL or a non-negative number.",call.=FALSE)
  
  is_cellqc <- function(x)
    is.list(x) && is.data.frame(x$qc) &&
    all(c("nCount_RNA","nFeature_RNA") %in% names(x$qc))
  
  med <- function(x)
    if(all(is.na(x))) NA_real_ else stats::median(x,na.rm=TRUE)
  
  # ------------------------------------------------------------------------- #
  # Input
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
        "percent.ribo","percent.hb","complexity"),
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
  
  qc$sample <- factor(qc$sample,levels=unique(qc$sample))
  qc$group <- factor(qc$group,levels=unique(qc$group))
  
  metrics <- unique(intersect(metrics,names(qc)))
  metrics <- metrics[vapply(qc[metrics],is.numeric,logical(1))]
  
  if(!length(metrics))
    stop("No numeric QC metrics available.",call.=FALSE)
  
  # ------------------------------------------------------------------------- #
  # Colors
  # ------------------------------------------------------------------------- #
  make_cols <- function(x,cols=NULL) {
    
    defaults <- setNames(
      grDevices::hcl.colors(length(x),"Dark 3"),x
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
  
  fill.scale <- ggplot2::scale_fill_manual(
    values=color.values,drop=FALSE
  )
  
  # ------------------------------------------------------------------------- #
  # Threshold estimation
  # ------------------------------------------------------------------------- #
  estimate <- function(x,metric) {
    
    x <- x[is.finite(x)]
    
    if(!length(x))
      return(c(lower=NA_real_,upper=NA_real_))
    
    is.count <- metric %in% c("nCount_RNA","nFeature_RNA")
    is.mt <- metric=="percent.mt"
    
    fixed.low <- if(metric=="nCount_RNA") min.counts else
      if(metric=="nFeature_RNA") min.features else NA_real_
    
    fixed.high <- if(metric=="nCount_RNA") max.counts else
      if(metric=="nFeature_RNA") max.features else
        if(is.mt) max.mt else NULL
    
    if(is.null(fixed.high)) fixed.high <- Inf
    
    if(method=="fixed") {
      
      low <- fixed.low
      high <- fixed.high
      
    } else {
      
      y <- if(is.count) log10(pmax(x,0)+1) else x
      
      if(method=="MAD") {
        
        center <- stats::median(y)
        deviation <- stats::mad(y,constant=1)
        
        low <- center-mad.n*deviation
        high <- center+mad.n*deviation
        
      } else {
        
        q <- stats::quantile(
          y,probs=quantiles,na.rm=TRUE,names=FALSE
        )
        
        low <- q[1]
        high <- q[2]
      }
      
      if(is.count) {
        low <- pmax(0,10^low-1)
        high <- pmax(0,10^high-1)
      }
      
      if(cap.fixed) {
        if(is.finite(fixed.low)) low <- max(low,fixed.low)
        if(is.finite(fixed.high)) high <- min(high,fixed.high)
      }
    }
    
    if(is.mt) low <- NA_real_
    
    if(is.count && !upper.counts)
      high <- Inf
    
    c(lower=low,upper=high)
  }
  
  thresholds <- dplyr::bind_rows(
    lapply(levels(qc$sample),function(s) {
      
      d <- qc[qc$sample==s,,drop=FALSE]
      
      dplyr::bind_rows(lapply(metrics,function(m) {
        
        limits <- estimate(d[[m]],m)
        
        data.frame(
          sample=s,
          metric=m,
          lower=unname(limits["lower"]),
          upper=unname(limits["upper"]),
          method=method
        )
      }))
    })
  )
  
  thresholds$sample <- factor(
    thresholds$sample,levels=levels(qc$sample)
  )
  
  thresholds$metric <- factor(
    thresholds$metric,levels=metrics
  )
  
  # ------------------------------------------------------------------------- #
  # Per-cell QC classification
  # ------------------------------------------------------------------------- #
  qc$pass <- TRUE
  qc$fail_reason <- ""
  
  for(m in metrics) {
    
    limits <- thresholds[thresholds$metric==m,,drop=FALSE]
    
    lower <- limits$lower[match(qc$sample,limits$sample)]
    upper <- limits$upper[match(qc$sample,limits$sample)]
    
    value <- qc[[m]]
    
    fail.low <- !is.na(lower) & !is.na(value) & value<lower
    fail.high <- !is.na(upper) & !is.na(value) & value>upper
    fail <- fail.low | fail.high
    
    qc$pass[fail] <- FALSE
    
    qc$fail_reason[fail.low] <- paste0(
      qc$fail_reason[fail.low],m,":low;"
    )
    
    qc$fail_reason[fail.high] <- paste0(
      qc$fail_reason[fail.high],m,":high;"
    )
  }
  
  qc$status <- ifelse(qc$pass,"Pass","Flagged")
  qc$fail_reason <- sub(";$","",qc$fail_reason)
  qc$fail_reason[qc$pass] <- "Pass"
  
  # ------------------------------------------------------------------------- #
  # Summary
  # ------------------------------------------------------------------------- #
  summary <- dplyr::summarise(
    dplyr::group_by(qc,.data$sample,.data$group),
    cells=dplyr::n(),
    passed=sum(.data$pass),
    flagged=sum(!.data$pass),
    pct_pass=100*mean(.data$pass),
    median_UMI=med(.data$nCount_RNA),
    median_genes=med(.data$nFeature_RNA),
    median_mt=if("percent.mt" %in% names(qc))
      med(.data$percent.mt) else NA_real_,
    .groups="drop"
  )
  
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
  long <- tidyr::pivot_longer(
    qc,
    cols=dplyr::all_of(metrics),
    names_to="metric",
    values_to="value"
  )
  
  long$metric <- factor(long$metric,levels=metrics)
  
  limits <- tidyr::pivot_longer(
    thresholds,
    cols=c("lower","upper"),
    names_to="bound",
    values_to="threshold"
  )
  
  limits <- limits[is.finite(limits$threshold),,drop=FALSE]
  
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
      width=.13,outlier.shape=NA,fill="white",
      alpha=.75,na.rm=TRUE
    )+
    ggplot2::geom_point(
      data=limits,
      ggplot2::aes(
        x=.data$sample,
        y=.data$threshold,
        shape=.data$bound
      ),
      inherit.aes=FALSE,
      size=2.5,
      colour="black"
    )+
    ggplot2::scale_shape_manual(
      values=c(lower=25,upper=24)
    )+
    ggplot2::facet_wrap(
      ~metric,scales="free_y",ncol=2
    )+
    fill.scale+
    ggplot2::labs(
      title="QC distributions and thresholds",
      x=NULL,y=NULL,shape="Threshold"
    )+th
  
  # ------------------------------------------------------------------------- #
  # 2. QC retention
  # ------------------------------------------------------------------------- #
  p_ret <- ggplot2::ggplot(
    summary,
    ggplot2::aes(
      x=.data$sample,
      y=.data$pct_pass,
      fill=.data[[color.col]]
    )
  )+
    ggplot2::geom_col(width=.7)+
    ggplot2::coord_cartesian(ylim=c(0,100))+
    fill.scale+
    ggplot2::labs(
      title="QC retention by sample",
      x=NULL,y="Passing cells (%)"
    )+th
  
  # ------------------------------------------------------------------------- #
  # 3. Exclusion reasons
  # ------------------------------------------------------------------------- #
  reasons <- qc[!qc$pass,c("sample","fail_reason"),drop=FALSE]
  
  if(nrow(reasons)) {
    
    reasons <- tidyr::separate_rows(
      reasons,.data$fail_reason,sep=";"
    )
    
    reasons <- dplyr::count(
      reasons,.data$sample,.data$fail_reason,
      name="cells"
    )
    
    p_excl <- ggplot2::ggplot(
      reasons,
      ggplot2::aes(
        x=.data$sample,
        y=.data$cells,
        fill=.data$fail_reason
      )
    )+
      ggplot2::geom_col()+
      ggplot2::labs(
        title="QC exclusion criteria",
        x=NULL,y="Flagged cells",
        fill="Criterion"
      )+th
    
  } else {
    
    p_excl <- NULL
  }
  
  # ------------------------------------------------------------------------- #
  # 4. UMI vs genes — QC status
  # ------------------------------------------------------------------------- #
  set.seed(seed)
  
  idx <- unlist(
    lapply(split(seq_len(nrow(qc)),qc$sample),function(i)
      if(length(i)>max.points) sample(i,max.points) else i
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
      colour=.data$status
    )
  )+
    ggplot2::geom_point(size=.45,alpha=.4)+
    ggplot2::scale_colour_manual(
      values=c(Pass="#377EB8",Flagged="#E41A1C")
    )+
    ggplot2::scale_x_log10()+
    ggplot2::scale_y_log10()+
    ggplot2::facet_wrap(~sample)+
    ggplot2::labs(
      title="UMI vs genes — QC classification",
      x="UMI (log10)",
      y="Genes (log10)",
      colour="QC status"
    )+th
  
  # ------------------------------------------------------------------------- #
  # 5. Threshold heatmap
  # ------------------------------------------------------------------------- #
  heat <- thresholds |>
    dplyr::mutate(
      label=paste0(
        ifelse(is.finite(.data$lower),
               round(.data$lower,1),"-"),
        " / ",
        ifelse(is.finite(.data$upper),
               round(.data$upper,1),"-")
      )
    )
  
  p_heat <- ggplot2::ggplot(
    heat,
    ggplot2::aes(
      x=.data$metric,
      y=.data$sample,
      fill=.data$upper
    )
  )+
    ggplot2::geom_tile(colour="white")+
    ggplot2::geom_text(
      ggplot2::aes(label=.data$label),size=3
    )+
    ggplot2::scale_fill_viridis_c(
      option="viridis",na.value="grey90"
    )+
    ggplot2::labs(
      title="Sample-specific QC thresholds",
      subtitle="Lower / upper",
      x=NULL,y=NULL,fill="Upper"
    )+th
  
  # ------------------------------------------------------------------------- #
  # Dashboard
  # ------------------------------------------------------------------------- #
  plots <- list(
    distributions=p_dist,
    retention=p_ret,
    exclusion=p_excl,
    scatter=p_scatter,
    thresholds=p_heat
  )
  
  plots <- Filter(Negate(is.null),plots)
  
  dashboard <- patchwork::wrap_plots(plots,ncol=2)+
    patchwork::plot_annotation(
      title=paste0("scQC — QC thresholds (",method,")")
    )
  
  # ------------------------------------------------------------------------- #
  # Parameters
  # ------------------------------------------------------------------------- #
  params <- list(
    sample.col=sample.col,
    group.col=group.col,
    species=species,
    type=type,
    method=method,
    metrics=metrics,
    mad.n=mad.n,
    quantiles=quantiles,
    min.features=min.features,
    min.counts=min.counts,
    max.features=max.features,
    max.counts=max.counts,
    max.mt=max.mt,
    cap.fixed=cap.fixed,
    upper.counts=upper.counts,
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
      qc,file.path(outdir,"qc_threshold_cells.csv"),
      row.names=FALSE
    )
    
    utils::write.csv(
      thresholds,file.path(outdir,"qc_threshold_thresholds.csv"),
      row.names=FALSE
    )
    
    utils::write.csv(
      summary,file.path(outdir,"qc_threshold_summary.csv"),
      row.names=FALSE
    )
    
    ggplot2::ggsave(
      file.path(outdir,"qc_threshold_dashboard.pdf"),
      dashboard,width=14,height=12
    )
    
    for(nm in names(plots))
      ggplot2::ggsave(
        file.path(outdir,paste0("qc_threshold_",nm,".pdf")),
        plots[[nm]],width=10,height=7
      )
    
    saveRDS(
      list(
        thresholds=thresholds,
        summary=summary,
        parameters=params
      ),
      file.path(outdir,"qc_threshold_results.rds")
    )
  }
  
  if(verbose)
    message(
      "[QC_THRESHOLD] ",method,
      " | ",nrow(summary)," samples",
      " | ",sum(qc$pass),"/",nrow(qc)," passed (",
      round(100*mean(qc$pass),1),"%)"
    )
  
  structure(
    list(
      qc=qc,
      thresholds=thresholds,
      summary=summary,
      plots=plots,
      dashboard=dashboard,
      parameters=params
    ),
    class=c("qc_threshold_result","list")
  )
}