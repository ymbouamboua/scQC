make_test_matrix <- function() {
  matrix(
    c(
      10, 0, 4, 8, 2, 7,
      0, 12, 3, 6, 1, 5,
      5, 4, 9, 0, 3, 8,
      2, 6, 1, 7, 10, 4,
      8, 3, 5, 2, 6, 9
    ),
    nrow = 5,
    dimnames = list(
      c("MT-ND1", "RPL3", "GeneA", "GeneB", "GeneC"),
      paste0("Cell", 1:6)
    )
  )
}

test_that("run_qc accepts a Seurat object", {
  counts <- Matrix::Matrix(
    matrix(
      c(
        5, 0, 2,
        0, 4, 1,
        3, 1, 0
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(
        c("MT-ND1", "GeneA", "GeneB"),
        c("Cell1", "Cell2", "Cell3")
      )
    ),
    sparse = TRUE
  )

  object <- SeuratObject::CreateSeuratObject(
    counts = counts,
    project = "SeuratProject",
    min.cells = 0,
    min.features = 0
  )

  result <- run_qc(
    x = object,
    sample_id = "TestSample",
    species = "human",
    min_feat = 0,
    min_umi = 0,
    max_mito = 100,
    method = "none",
    rm_dbl = FALSE,
    save_object = FALSE,
    outdir = NULL,
    return = "result",
    verbose = FALSE
  )

  expect_s3_class(
    result,
    "scqc_result"
  )

  expect_s4_class(
    result$object,
    "Seurat"
  )

  expect_equal(
    result$summary$sample,
    "TestSample"
  )

  expect_equal(
    result$summary$input_cells,
    3
  )

  expect_equal(
    result$summary$after_qc_cells,
    3
  )

  expect_equal(
    result$summary$final_cells,
    3
  )

  expect_equal(
    ncol(result$object),
    3
  )
})


test_that("run_qc can return only the Seurat object", {
  counts <- Matrix::Matrix(
    matrix(
      c(
        5, 0, 2,
        0, 4, 1,
        3, 1, 0
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(
        c("MT-ND1", "GeneA", "GeneB"),
        c("Cell1", "Cell2", "Cell3")
      )
    ),
    sparse = TRUE
  )

  object <- SeuratObject::CreateSeuratObject(
    counts = counts,
    project = "TestSample",
    min.cells = 0,
    min.features = 0
  )

  qc_object <- run_qc(
    x = object,
    sample_id = "TestSample",
    species = "human",
    min_feat = 0,
    min_umi = 0,
    max_mito = 100,
    method = "none",
    rm_dbl = FALSE,
    save_object = FALSE,
    outdir = NULL,
    return = "object",
    verbose = FALSE
  )

  expect_s4_class(
    qc_object,
    "Seurat"
  )

  expect_equal(
    ncol(qc_object),
    3
  )
})



test_that("run_qc records the MAD multiplier", {
  test_matrix <- make_test_matrix()

  result <- run_qc(
    x = test_matrix,
    sample_id = "sample1",
    species = "human",
    method = "MAD",
    mad_n = 5,
    min_feat = 0,
    min_umi = 0,
    max_mito = 100,
    outdir = NULL,
    save_object = FALSE,
    verbose = FALSE
  )

  expect_true("mad_n" %in% names(result$thresholds))
  expect_equal(result$thresholds$mad_n, 5)
  expect_equal(result$thresholds$method, "MAD")
})


test_that("mad_n is NA when MAD thresholds are not used", {
  test_matrix <- make_test_matrix()

  result <- run_qc(
    x = test_matrix,
    sample_id = "sample1",
    species = "human",
    method = "fixed",
    fixed_thr = list(
      max_feat = 1000,
      max_umi = 10000
    ),
    min_feat = 0,
    min_umi = 0,
    max_mito = 100,
    outdir = NULL,
    save_object = FALSE,
    verbose = FALSE
  )

  expect_true("mad_n" %in% names(result$thresholds))
  expect_true(is.na(result$thresholds$mad_n))
  expect_equal(result$thresholds$method, "fixed")
})

