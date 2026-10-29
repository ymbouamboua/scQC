ambient_fixture <- function() {
  set.seed(9)
  m <- matrix(rpois(20*60,1),20,dimnames=list(c("MT-ND1",paste0("G",2:20)),paste0("cell",1:60)))
  m[2:6,1:30] <- m[2:6,1:30]+10
  m[7:11,31:60] <- m[7:11,31:60]+10
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(m,sparse=TRUE))
  object$cluster <- rep(c("A","B"),each=30)
  object$library <- rep(c("L1","L2"),30)
  object
}
ambient_qc <- function(object,...) run_qc(object,sample_id="s",min_feat=0,min_umi=0,
  max_mito=100,method="none",outdir=NULL,save_object=FALSE,verbose=FALSE,...)

test_that("decontamination is opt-in, defaults to decontX and compact mode preserves input", {
  testthat::with_mocked_bindings({
  obj <- ambient_fixture();raw <- .qc_get_layer(obj,"RNA","counts")
  original <- obj
  no <- ambient_qc(obj)
  expect_false("ambient_contamination" %in% names(no$object[[]]))
  expect_equal(.qc_get_layer(no$object,"RNA","counts"),raw)
  yes <- ambient_qc(obj,decontaminate=TRUE)
  expect_identical(names(yes$object@assays),"RNA")
  expect_equal(.qc_get_layer(yes$object,"RNA","counts"),raw*.8)
  expect_equal(unname(yes$object$nCount_RNA),as.numeric(Matrix::colSums(raw*.8)))
  expect_equal(unname(yes$object$qc_raw_nCount_RNA),as.numeric(Matrix::colSums(raw)))
  expect_equal(unname(yes$object$ambient_contamination),rep(.2,ncol(obj)))
  expect_identical(yes$parameters$ambient_method,"decontX")
  expect_identical(obj,original)
  expect_equal(yes$ambient$removed_fraction,.2)
  expect_false(any(vapply(yes$object@misc$scqc_ambient,inherits,logical(1),"Matrix")))
  expect_equal(yes$summary$final_cells,60)
  expect_error(ambient_qc(yes$object,decontaminate=TRUE),"already ambient-corrected")
  },
    .qc_ambient_backend=function(method) "Matrix",
    .qc_ambient_fit=function(counts,...) list(counts=counts*.8,contamination=rep(.2,ncol(counts))),
    .package="scQC")
})

test_that("optional raw assay is retained and stale compact analyses are removed", {
  testthat::with_mocked_bindings({
  obj <- ambient_fixture()
  xy <- matrix(seq_len(120),60,dimnames=list(colnames(obj),c("UMAP_1","UMAP_2")))
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(xy,key="UMAP_",assay="RNA")
  obj <- Seurat::NormalizeData(obj,verbose=FALSE)
  raw <- .qc_get_layer(obj,"RNA","counts")
  both <- ambient_qc(obj,decontaminate=TRUE,keep_raw=TRUE)
  expect_setequal(names(both$object@assays),c("RNA","RNAdecont"))
  expect_equal(.qc_get_layer(both$object,"RNA","counts"),raw)
  expect_equal(.qc_get_layer(both$object,"RNAdecont","counts"),raw*.8)
  expect_identical(SeuratObject::DefaultAssay(both$object),"RNAdecont")
  expect_true("umap" %in% names(both$object@reductions))
  compact <- ambient_qc(obj,decontaminate=TRUE)
  expect_length(compact$object@reductions,0)
  expect_false("data" %in% SeuratObject::Layers(compact$object[["RNA"]],search=NA))
  expect_true("umap" %in% compact$object@misc$scqc_ambient$invalidated$reductions)
  expect_error(ambient_qc(both$object,decontaminate=TRUE,keep_raw=TRUE),"already exists")
  },
    .qc_ambient_backend=function(method) "Matrix",
    .qc_ambient_fit=function(counts,...) list(counts=counts*.8,contamination=rep(.2,ncol(counts))),
    .package="scQC")
})

test_that("libraries are corrected separately and background excludes all original called cells", {
  obj <- ambient_fixture();raw <- .qc_get_layer(obj,"RNA","counts")
  seen <- list()
  testthat::with_mocked_bindings({
  empty <- raw[,1:3];colnames(empty)<-paste0("e",1:3)
  background <- cbind(raw,empty)
  obj$sample <- obj$library
  # ambient_batch is captured before the old run_qc sample metadata assignment.
  result <- ambient_qc(obj,decontaminate=TRUE,ambient_batch="sample",
    ambient_raw=list(L1=background,L2=background))
  expect_length(seen,2)
  expect_equal(vapply(seen,function(x)length(x$cells),integer(1)),c(30L,30L))
  expect_true(all(vapply(seen,function(x)setequal(x$background,paste0("e",1:3)),logical(1))))
  expect_equal(nrow(result$ambient),2)
  expect_identical(colnames(result$object),colnames(obj))
  expect_error(ambient_qc(obj,decontaminate=TRUE,ambient_batch="library",ambient_raw=background),"named by library")
  expect_error(ambient_qc(obj,decontaminate=TRUE,ambient_raw=raw),"No nonempty background")
  },
    .qc_ambient_backend=function(method) "Matrix",
    .qc_ambient_fit=function(counts,background,clusters,method,backend,args,seed,verbose) {
      seen[[length(seen)+1L]] <<- list(cells=colnames(counts),background=colnames(background))
      list(counts=counts*.9,contamination=rep(.1,ncol(counts)))
    },.package="scQC")
})

test_that("multiple inputs export diagnostics and preserve merged provenance", {
  testthat::with_mocked_bindings({;obj<-ambient_fixture()
  d <- tempfile();dir.create(d);on.exit(unlink(d,recursive=TRUE))
  result <- run_qc(list(A=obj,B=obj),min_feat=0,min_umi=0,max_mito=100,method="none",
    decontaminate=TRUE,merge=TRUE,outdir=d,save_object=FALSE,verbose=FALSE)
  expect_equal(nrow(result$ambient),2)
  expect_equal(nrow(result$ambient_cells),120)
  expect_true(all(result$ambient_cells$merged_cell %in% colnames(result$object)))
  expect_setequal(names(result$object@misc$scqc_ambient$inputs),c("A","B"))
  expect_true(file.exists(result$files$ambient))
  expect_true(file.exists(result$files$ambient_cells))
  expect_equal(nrow(utils::read.delim(result$files$ambient_cells)),120)
  expect_error(run_qc(list(A=obj,B=obj),decontaminate=TRUE,ambient_raw=matrix(1,2,2),outdir=NULL,verbose=FALSE),"named exactly")
  },
    .qc_ambient_backend=function(method) "Matrix",
    .qc_ambient_fit=function(counts,...) list(counts=counts*.8,contamination=rep(.2,ncol(counts))),
    .package="scQC")
})

test_that("invalid ambient arguments, matrices and backend outputs fail clearly", {
  testthat::with_mocked_bindings({;obj<-ambient_fixture()
  expect_error(ambient_qc(obj,decontaminate=TRUE,ambient_method="SoupX"),"ambient_clusters")
  expect_error(ambient_qc(obj,decontaminate=TRUE,ambient_args=list(seed=5)),"reserved")
  expect_error(ambient_qc(obj,decontaminate=TRUE,keep_raw=TRUE,ambient_assay="RNA"),"distinct")
  bad <- .qc_get_layer(obj,"RNA","counts");bad@x[1]<-.5
  expect_error(ambient_qc(bad,decontaminate=TRUE),"raw integer")
  },
    .qc_ambient_backend=function(method) "Matrix",
    .qc_ambient_fit=function(counts,...) list(counts=counts*.8,contamination=rep(.2,ncol(counts))),
    .package="scQC")
})

test_that("invalid backend totals are rejected", {
  obj <- ambient_fixture()
  testthat::with_mocked_bindings({
  expect_error(ambient_qc(obj,decontaminate=TRUE),"increased cell totals")
  },.qc_ambient_backend=function(method) "Matrix",
    .qc_ambient_fit=function(counts,...)
    list(counts=counts*2,contamination=rep(.1,ncol(counts))),.package="scQC")
})

test_that("decontX backend corrects actual counts when available", {
  if (!requireNamespace("decontX",quietly=TRUE) && !requireNamespace("celda",quietly=TRUE)) skip("decontX/celda not installed")
  obj <- ambient_fixture();raw <- .qc_get_layer(obj,"RNA","counts")
  result <- ambient_qc(obj,decontaminate=TRUE,ambient_clusters="cluster",ambient_args=list(maxIter=20))
  expect_false(any(result$ambient$backend == "Matrix"))
  corrected <- .qc_get_layer(result$object,"RNA","counts")
  expect_equal(dimnames(corrected),dimnames(raw))
  expect_true(sum(corrected)<sum(raw))
  expect_true(all(result$object$ambient_contamination>=0 & result$object$ambient_contamination<=1))
  expect_equal(sum(result$object$nCount_RNA),sum(corrected))
})

test_that("SoupX backend accepts a measured explicit fraction and raw droplets", {
  skip_if_not_installed("SoupX")
  obj <- ambient_fixture();raw <- .qc_get_layer(obj,"RNA","counts")
  set.seed(12)
  bg <- Matrix::Matrix(matrix(rpois(20*30,1),20,dimnames=list(rownames(raw),paste0("e",1:30))),sparse=TRUE)
  # SoupX 1.6.2 emits its own Matrix deprecation warning; do not silence
  # warnings from the scQC implementation itself.
  result <- suppressWarnings(ambient_qc(obj,decontaminate=TRUE,ambient_method="SoupX",
    ambient_clusters="cluster",ambient_raw=bg,ambient_args=list(contamination_fraction=.1)))
  expect_false(any(result$ambient$backend == "Matrix"))
  corrected <- .qc_get_layer(result$object,"RNA","counts")
  expect_equal(dimnames(corrected),dimnames(raw))
  expect_equal(unname(result$object$ambient_contamination),rep(.1,ncol(obj)))
  expect_lt(sum(corrected),sum(raw))
  expect_true(all(is.finite(corrected@x)))
})
