printed <- function(x) {
  msg <- utils::capture.output(out <- utils::capture.output(print(x)), type = "message")
  paste(c(out, msg), collapse = "\n")
}

test_that("bio_homology_cv() never splits a cluster", {
  set.seed(1)
  folds <- bio_homology_cv(ecoli_cds, protein, v = 5, k = 2, threshold = 0.2)
  expect_s3_class(folds, c("bio_homology_cv", "rset"))
  expect_identical(nrow(folds), 5L)
  clusters <- attr(folds, "clusters")
  expect_length(clusters, nrow(ecoli_cds))
  assessed <- integer()
  for (sp in folds$splits) {
    a <- rsample::complement(sp)
    an <- sp$in_id
    expect_length(intersect(clusters[a], clusters[an]), 0L)
    assessed <- c(assessed, a)
  }
  expect_setequal(assessed, seq_len(nrow(ecoli_cds)))
  expect_match(printed(folds), "homology")
})

test_that("identical sequences share a fold", {
  d <- tibble::tibble(s = c(rep("ACGTACGTAAAC", 3), "GGGGCCCCGGGG", "TTTTATATTTTA", "CACACACACA"))
  set.seed(3)
  folds <- bio_homology_cv(d, s, v = 2, k = 2, threshold = 0.05)
  cl <- attr(folds, "clusters")
  expect_identical(length(unique(cl[1:3])), 1L)
  expect_error(bio_homology_cv(d, s, v = 10, k = 2, threshold = 0.05), "fewer than")
  expect_error(bio_homology_cv(d, s, threshold = 2), "between 0 and 1")
  expect_error(bio_homology_cv(d, s, v = 1), "at least 2")
})

test_that("bio_batch_cv() leaves whole batches out", {
  folds <- bio_batch_cv(qpcr_sim, plate)
  expect_s3_class(folds, c("bio_batch_cv", "rset"))
  expect_identical(nrow(folds), 2L)
  for (sp in folds$splits) {
    expect_length(unique(rsample::assessment(sp)$plate), 1L)
    expect_length(intersect(rsample::assessment(sp)$plate, rsample::analysis(sp)$plate), 0L)
  }
  expect_match(printed(folds), "Leave-batch-out")
  by_sample <- bio_batch_cv(qpcr_sim, sample, v = 3)
  expect_identical(nrow(by_sample), 3L)
  for (sp in by_sample$splits) {
    expect_length(intersect(rsample::assessment(sp)$sample, rsample::analysis(sp)$sample), 0L)
  }
  expect_error(bio_batch_cv(qpcr_sim, plate, v = 5), "cannot exceed")
  expect_error(bio_batch_cv(qpcr_sim[qpcr_sim$plate == "plate_1", ], plate), "at least 2")
})
