test_that(".infer_sample_id() identifies samples from 10x directories", {
  temp_root <- tempfile("scqc_test_")
  dir.create(temp_root)
  
  sample_dir <- file.path(temp_root, "FN_PCW6")
  outs_dir <- file.path(sample_dir, "outs")
  
  filtered_dir <- file.path(
    outs_dir,
    "filtered_feature_bc_matrix"
  )
  
  raw_dir <- file.path(
    outs_dir,
    "raw_feature_bc_matrix"
  )
  
  dir.create(
    filtered_dir,
    recursive = TRUE
  )
  
  dir.create(
    raw_dir,
    recursive = TRUE
  )
  
  expect_equal(
    .infer_sample_id(
      original_path = sample_dir,
      resolved_path = filtered_dir
    ),
    "FN_PCW6"
  )
  
  expect_equal(
    .infer_sample_id(
      original_path = outs_dir,
      resolved_path = filtered_dir
    ),
    "FN_PCW6"
  )
  
  expect_equal(
    .infer_sample_id(
      original_path = filtered_dir,
      resolved_path = filtered_dir
    ),
    "FN_PCW6"
  )
  
  expect_equal(
    .infer_sample_id(
      original_path = raw_dir,
      resolved_path = raw_dir
    ),
    "FN_PCW6"
  )
  
  unlink(
    temp_root,
    recursive = TRUE
  )
})


test_that(".infer_sample_id() identifies samples from named HDF5 files", {
  temp_root <- tempfile("scqc_test_")
  dir.create(temp_root)
  
  filtered_h5 <- file.path(
    temp_root,
    "FN_PCW6_filtered_feature_bc_matrix.h5"
  )
  
  raw_h5 <- file.path(
    temp_root,
    "FN_PCW6_raw_feature_bc_matrix.h5"
  )
  
  file.create(filtered_h5)
  file.create(raw_h5)
  
  expect_equal(
    .infer_sample_id(
      original_path = filtered_h5,
      resolved_path = filtered_h5
    ),
    "FN_PCW6"
  )
  
  expect_equal(
    .infer_sample_id(
      original_path = raw_h5,
      resolved_path = raw_h5
    ),
    "FN_PCW6"
  )
  
  unlink(
    temp_root,
    recursive = TRUE
  )
})


test_that(".infer_sample_id() identifies samples from generic HDF5 names", {
  temp_root <- tempfile("scqc_test_")
  dir.create(temp_root)
  
  sample_dir <- file.path(
    temp_root,
    "FN_PCW6"
  )
  
  dir.create(sample_dir)
  
  filtered_h5 <- file.path(
    sample_dir,
    "filtered_feature_bc_matrix.h5"
  )
  
  raw_h5 <- file.path(
    sample_dir,
    "raw_feature_bc_matrix.h5"
  )
  
  file.create(filtered_h5)
  file.create(raw_h5)
  
  expect_equal(
    .infer_sample_id(
      original_path = filtered_h5,
      resolved_path = filtered_h5
    ),
    "FN_PCW6"
  )
  
  expect_equal(
    .infer_sample_id(
      original_path = raw_h5,
      resolved_path = raw_h5
    ),
    "FN_PCW6"
  )
  
  unlink(
    temp_root,
    recursive = TRUE
  )
})


test_that(".qc_infer_input_name() identifies a Seurat sample", {
  skip_if_not_installed("SeuratObject")
  
  counts <- Matrix::Matrix(
    matrix(
      c(
        1, 0, 3,
        0, 2, 1,
        4, 1, 0
      ),
      nrow = 3,
      dimnames = list(
        c("Gene1", "Gene2", "Gene3"),
        c("Cell1", "Cell2", "Cell3")
      )
    ),
    sparse = TRUE
  )
  
  object <- SeuratObject::CreateSeuratObject(
    counts = counts,
    project = "FN_PCW6",
    min.cells = 0,
    min.features = 0
  )
  
  object$orig.ident <- "FN_PCW6"
  
  expect_equal(
    .qc_infer_input_name(object),
    "FN_PCW6"
  )
})


test_that(".qc_infer_input_name() falls back for mixed Seurat identities", {
  skip_if_not_installed("SeuratObject")
  
  counts <- Matrix::Matrix(
    matrix(
      c(
        1, 0, 3,
        0, 2, 1,
        4, 1, 0
      ),
      nrow = 3,
      dimnames = list(
        c("Gene1", "Gene2", "Gene3"),
        c("Cell1", "Cell2", "Cell3")
      )
    ),
    sparse = TRUE
  )
  
  object <- SeuratObject::CreateSeuratObject(
    counts = counts,
    min.cells = 0,
    min.features = 0
  )
  
  object$orig.ident <- c(
    "Sample_A",
    "Sample_A",
    "Sample_B"
  )
  
  expect_equal(
    .qc_infer_input_name(
      object,
      fallback = "Sample_1"
    ),
    "Sample_1"
  )
})


test_that(".qc_infer_input_name() delegates path inference", {
  temp_root <- tempfile("scqc_test_")
  
  matrix_dir <- file.path(
    temp_root,
    "FN_PCW6",
    "outs",
    "filtered_feature_bc_matrix"
  )
  
  dir.create(
    matrix_dir,
    recursive = TRUE
  )
  
  testthat::local_mocked_bindings(
    .resolve_10x_path = function(path) {
      normalizePath(
        path,
        winslash = "/",
        mustWork = TRUE
      )
    }
  )
  
  expect_equal(
    .qc_infer_input_name(matrix_dir),
    "FN_PCW6"
  )
  
  unlink(
    temp_root,
    recursive = TRUE
  )
})


test_that(".qc_infer_input_name() uses fallback for matrices", {
  counts <- Matrix::Matrix(
    matrix(
      1:9,
      nrow = 3
    ),
    sparse = TRUE
  )
  
  expect_equal(
    .qc_infer_input_name(
      counts,
      fallback = "Sample_1"
    ),
    "Sample_1"
  )
})


test_that(".qc_infer_input_name() handles missing paths", {
  missing_path <- file.path(
    tempdir(),
    "nonexistent_sample",
    "filtered_feature_bc_matrix"
  )
  
  expect_equal(
    .qc_infer_input_name(
      missing_path,
      fallback = "Sample_1"
    ),
    "Sample_1"
  )
})