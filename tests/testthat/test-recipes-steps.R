printed <- function(x) {
  msg <- utils::capture.output(out <- utils::capture.output(print(x)), type = "message")
  paste(c(out, msg), collapse = "\n")
}

test_that("step_bio_kmers() learns a vocabulary and bakes consistently", {
  train <- ecoli_cds[1:40, c("class", "sequence")]
  test <- ecoli_cds[41:76, c("class", "sequence")]
  rec <- recipes::recipe(class ~ sequence, data = train) |>
    step_bio_kmers(sequence, k = 2)
  expect_match(printed(rec), "K-mer")
  untrained <- tidy(rec, number = 1)
  expect_true(is.na(untrained$n_kmers))
  rec <- recipes::prep(rec)
  baked <- recipes::bake(rec, new_data = test)
  expect_false("sequence" %in% names(baked))
  expect_identical(sum(startsWith(names(baked), "sequence_")), 16L)
  expect_equal(rowSums(baked[, startsWith(names(baked), "sequence_")]), rep(1, 36),
               tolerance = 1e-8)
  td <- tidy(rec, number = 1)
  expect_identical(td$n_kmers, 16L)
  expect_identical(required_pkgs(rec$steps[[1]]), "tidybiotek")
  keep <- recipes::recipe(class ~ sequence, data = train) |>
    step_bio_kmers(sequence, k = 1, normalise = FALSE, keep_original_cols = TRUE) |>
    recipes::prep()
  kb <- recipes::bake(keep, new_data = NULL)
  expect_true("sequence" %in% names(kb))
  expect_equal(kb$sequence_A + kb$sequence_C + kb$sequence_G + kb$sequence_T,
               nchar(train$sequence))
  expect_error(step_bio_kmers(recipes::recipe(class ~ sequence, data = train), sequence, k = 0))
})

test_that("step_bio_aa_props() adds protein properties", {
  rec <- recipes::recipe(class ~ protein, data = ecoli_cds) |>
    step_bio_aa_props(protein, properties = c("gravy", "isoelectric_point"))
  expect_identical(nrow(tidy(rec, number = 1)), 2L)
  rec <- recipes::prep(rec)
  baked <- recipes::bake(rec, new_data = ecoli_cds[1:3, ])
  expect_named(baked, c("class", "protein_gravy", "protein_isoelectric_point"))
  direct <- bio_aa_props(ecoli_cds[1:3, ], protein)
  expect_equal(baked$protein_gravy, direct$gravy)
  expect_match(printed(rec), "physicochemical")
  expect_identical(tidy(rec, number = 1)$property, c("gravy", "isoelectric_point"))
  expect_identical(required_pkgs(rec$steps[[1]]), "tidybiotek")
  expect_error(
    step_bio_aa_props(recipes::recipe(class ~ protein, data = ecoli_cds), protein,
      properties = "colour"
    ),
    "Unknown"
  )
})

test_that("step_bio_vst_norm() stabilises variance and uses training values", {
  wide <- bio_counts_wide(expression_sim)
  expect_identical(dim(wide), c(8L, 1002L))
  rec <- recipes::recipe(condition ~ ., data = wide) |>
    recipes::update_role(sample, new_role = "id") |>
    step_bio_vst_norm(recipes::all_numeric_predictors())
  expect_true(all(is.na(tidy(rec, number = 1)$log_geo_mean)))
  rec <- recipes::prep(rec)
  baked <- recipes::bake(rec, new_data = NULL)
  vals <- as.matrix(baked[, startsWith(names(baked), "gene_")])
  expect_true(all(is.finite(vals)))
  means <- colMeans(vals)
  sds <- apply(vals, 2, stats::sd)
  raw <- as.matrix(wide[, startsWith(names(wide), "gene_")])
  raw_sd <- apply(raw, 2, stats::sd)
  expect_lt(abs(stats::cor(means, sds, method = "spearman")),
            abs(stats::cor(colMeans(raw), raw_sd, method = "spearman")))
  one <- recipes::bake(rec, new_data = wide[1, ])
  expect_equal(unname(unlist(one[1, 2:5])), unname(unlist(baked[1, 2:5])))
  expect_identical(tidy(rec, number = 1)$method[[1]], "vst")
  expect_match(printed(rec), "Count normalisation")
  expect_identical(required_pkgs(rec$steps[[1]]), "tidybiotek")
  for (m in c("log_cpm", "anscombe")) {
    r <- recipes::recipe(condition ~ ., data = wide) |>
      recipes::update_role(sample, new_role = "id") |>
      step_bio_vst_norm(recipes::all_numeric_predictors(), method = m) |>
      recipes::prep()
    bk <- recipes::bake(r, new_data = NULL)
    expect_true(all(is.finite(as.matrix(bk[, startsWith(names(bk), "gene_")]))))
  }
  anscombe <- recipes::recipe(condition ~ ., data = wide) |>
    recipes::update_role(sample, new_role = "id") |>
    step_bio_vst_norm(gene_0001, method = "anscombe") |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)
  expect_equal(anscombe$gene_0001, 2 * sqrt(wide$gene_0001 + 3 / 8))
  neg <- wide
  neg$gene_0001[1] <- -1
  expect_error(
    recipes::recipe(condition ~ ., data = neg) |>
      recipes::update_role(sample, new_role = "id") |>
      step_bio_vst_norm(recipes::all_numeric_predictors()) |>
      recipes::prep(),
    "non-negative"
  )
})
