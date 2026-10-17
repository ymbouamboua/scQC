#' Diagnose Cell Ranger scRNA-seq or snRNA-seq count quality
#'
#' Diagnostic only: never removes barcodes. Uses existing scQC input and
#' count-layer helpers. A barcode-rank plot from filtered 10x data cannot
#' identify empty droplets.
#'
#' @param x Seurat object, counts matrix, Cell Ranger directory or H5 file.
#' @param sample_id Sample identifier.
#' @param assay Assay name.
#' @param layer Raw-count layer.
#' @param species Human or mouse.
#' @param type scRNA or snRNA (controls default mitochondrial cutoff only).
#' @param min_feat,min_umi Lower thresholds for detected genes and UMIs.
#' @param max_feat,max_umi Optional upper thresholds.
#' @param max_mito Optional mitochondrial percentage threshold.
#' @param min_complexity Optional lower detected-genes / UMI threshold.
#' @param max_points Maximum number of points in scatter plots.
#' @param outdir Output directory; NULL disables writing.
#' @param seed Random seed for plotting subsample.
#' @param verbose Print status.
#' @return A list with object, qc, summary, plots, dashboard and parameters.
#' @export
cellqc <- function(x, sample_id="Sample", assay="RNA", layer="counts",
                   species=c("human","mouse"), type=c("scRNA","snRNA"),
                   min_feat=200, min_umi=500, max_feat=Inf, max_umi=Inf,
                   max_mito=NULL, min_complexity=NULL, max_points=50000L,
                   outdir="QC", seed=1234L, verbose=TRUE) {
  species <- match.arg(species)
  type <- match.arg(type)
  if(is.null(max_mito)) max_mito <- if(type=="snRNA") 5 else 20
  for(pkg in c("Seurat","SeuratObject","Matrix","ggplot2","patchwork"))
    if(!requireNamespace(pkg,quietly=TRUE)) stop("Package required: ",pkg,call.=FALSE)
  if(!is.numeric(max_points) || length(max_points)!=1L || is.na(max_points) || max_points<1)
    stop("`max_points` must be a positive number.",call.=FALSE)
  obj <- read_sc_input(x,sample_id=sample_id,assay=assay,verbose=verbose)
  mat <- .qc_get_layer(obj,assay=assay,layer=layer)
  if(any(mat@x < 0)) stop("Counts cannot be negative.",call.=FALSE)
  umi <- Matrix::colSums(mat)
  genes <- Matrix::colSums(mat>0)
  rn <- rownames(mat)
  mt <- grepl(if(species=="human") "^MT-" else "^mt-",rn)
  rb <- grepl(if(species=="human") "^RP[LS]" else "^Rp[ls]",rn)
  hb <- grepl(if(species=="human") "^HB[ABDEGQMZ]" else "^Hb[abdegmqz]",rn)
  pct <- function(ix) if(any(ix)) 100*as.numeric(Matrix::colSums(mat[ix,,drop=FALSE]))/pmax(umi,1) else rep(NA_real_,length(umi))
  qc <- data.frame(cell=colnames(mat),sample=sample_id,nCount_RNA=as.numeric(umi),
                   nFeature_RNA=as.numeric(genes),percent.mt=pct(mt),
                   percent.ribo=pct(rb),percent.hb=pct(hb),
                   complexity=as.numeric(genes)/pmax(umi,1),stringsAsFactors=FALSE)
  qc$log10_UMI <- log10(qc$nCount_RNA+1)
  qc$log10_genes <- log10(qc$nFeature_RNA+1)
  qc$low.features <- qc$nFeature_RNA<min_feat
  qc$high.features <- qc$nFeature_RNA>max_feat
  qc$low.counts <- qc$nCount_RNA<min_umi
  qc$high.counts <- qc$nCount_RNA>max_umi
  qc$high.mt <- !is.na(qc$percent.mt) & qc$percent.mt>max_mito
  qc$low.complexity <- if(is.null(min_complexity)) rep(FALSE,nrow(qc)) else qc$complexity<min_complexity
  qc$pass <- !Reduce(`|`,qc[c("low.features","high.features","low.counts","high.counts","high.mt","low.complexity")])
  qc$status <- factor(ifelse(qc$pass,"Pass","Flagged"),levels=c("Flagged","Pass"))
  med <- function(v) if(all(is.na(v))) NA_real_ else median(v,na.rm=TRUE)
  summary <- data.frame(sample=sample_id,barcodes=nrow(qc),features=nrow(mat),
                        median_UMI=med(qc$nCount_RNA),median_genes=med(qc$nFeature_RNA),
                        median_mt=med(qc$percent.mt),median_ribo=med(qc$percent.ribo),
                        median_complexity=med(qc$complexity),pass=sum(qc$pass),
                        flagged=sum(!qc$pass),pct_pass=100*mean(qc$pass))
  obj <- SeuratObject::AddMetaData(obj,metadata=setNames(qc$percent.mt,qc$cell),col.name="percent.mt")
  obj <- SeuratObject::AddMetaData(obj,metadata=setNames(as.character(qc$status),qc$cell),col.name="cellqc_status")
  set.seed(seed)
  d <- qc[if(nrow(qc)>max_points) sample.int(nrow(qc),max_points) else seq_len(nrow(qc)),,drop=FALSE]
  th <- ggplot2::theme_bw(base_size=11)+ggplot2::theme(plot.title=ggplot2::element_text(face="bold"),panel.grid.minor=ggplot2::element_blank(),legend.position="bottom")
  cols <- c(Flagged="#E41A1C",Pass="#377EB8")
  histplot <- function(col,title,xlab,cut=NULL,fill="#377EB8") {
    p <- ggplot2::ggplot(qc,ggplot2::aes(x=.data[[col]]))+ggplot2::geom_histogram(bins=65,fill=fill)+ggplot2::labs(title=title,x=xlab,y="Barcodes")+th
    if(!is.null(cut) && is.finite(cut)) p <- p+ggplot2::geom_vline(xintercept=cut,linetype="dashed")
    p
  }
  scatter <- function(y,title,ylab,cut=NULL) {
    p <- ggplot2::ggplot(d,ggplot2::aes(x=.data$log10_UMI,y=.data[[y]],colour=.data$status))+
      ggplot2::geom_point(size=.5,alpha=.4)+ggplot2::scale_colour_manual(values=cols)+
      ggplot2::labs(title=title,x="log10(UMI + 1)",y=ylab,colour=NULL)+th
    if(!is.null(cut) && is.finite(cut)) p <- p+ggplot2::geom_hline(yintercept=cut,linetype="dashed")
    p
  }
  rank <- sort(qc$nCount_RNA[qc$nCount_RNA>0],decreasing=TRUE)
  ranks <- data.frame(rank=seq_along(rank),umi=rank)
  plots <- list(
    umi=histplot("log10_UMI","UMI distribution","log10(UMI + 1)",log10(min_umi+1)),
    genes=histplot("log10_genes","Detected genes","log10(Genes + 1)",log10(min_feat+1),"#4DAF4A"),
    mitochondrial=histplot("percent.mt","Mitochondrial fraction","Mitochondrial (%)",max_mito,"#E41A1C"),
    complexity=scatter("log10_genes","UMI vs detected genes","log10(Genes + 1)"),
    mt_scatter=scatter("percent.mt","Mitochondrial QC","Mitochondrial (%)",max_mito),
    retention=ggplot2::ggplot(qc,ggplot2::aes(x=.data$status,fill=.data$status))+
      ggplot2::geom_bar(width=.65)+ggplot2::scale_fill_manual(values=cols)+
      ggplot2::labs(title="QC retention",x=NULL,y="Barcodes")+th+ggplot2::theme(legend.position="none"),
    barcode_rank=ggplot2::ggplot(ranks,ggplot2::aes(x=.data$rank,y=.data$umi))+
      ggplot2::geom_line(colour="#377EB8",linewidth=.5)+ggplot2::scale_x_log10()+
      ggplot2::scale_y_log10()+ggplot2::labs(title="Barcode rank (input barcodes)",x="Barcode rank",y="UMI")+th
  )
  dashboard <- patchwork::wrap_plots(plots,ncol=2)+patchwork::plot_annotation(
    title=paste0("CellQC — ",sample_id),
    subtitle=sprintf("%s barcodes | %.1f%% pass QC",format(nrow(qc),big.mark=","),summary$pct_pass))
  params <- list(species=species,type=type,min_feat=min_feat,min_umi=min_umi,
                 max_feat=max_feat,max_umi=max_umi,max_mito=max_mito,
                 min_complexity=min_complexity)
  if(!is.null(outdir)) {
    dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
    utils::write.csv(qc,file.path(outdir,"cellqc_metadata.csv"),row.names=FALSE)
    utils::write.csv(summary,file.path(outdir,"cellqc_summary.csv"),row.names=FALSE)
    ggplot2::ggsave(file.path(outdir,"cellqc_dashboard.pdf"),dashboard,width=12,height=14)
    for(nm in names(plots)) ggplot2::ggsave(file.path(outdir,paste0("cellqc_",nm,".pdf")),plots[[nm]],width=6,height=5)
    saveRDS(list(qc=qc,summary=summary,parameters=params),file.path(outdir,"cellqc_results.rds"))
  }
  structure(list(object=obj,qc=qc,summary=summary,plots=plots,dashboard=dashboard,
                 parameters=params),class=c("cellqc_result","list"))
}
