test_that("bio_seq_tbl() cleans sequences and records the alphabet", {
  x <- bio_seq_tbl(tibble::tibble(id = 1:2, s = c("acg t", "GG-CA")), s)
  expect_s3_class(x, "bio_seq_tbl")
  expect_true(is_bio_seq_tbl(x))
  expect_identical(x$s, c("ACGT", "GGCA"))
  expect_identical(bio_alphabet(x), "dna")
  expect_identical(bio_alphabet(bio_seq_tbl(tibble::tibble(s = "ACGU"), s)), "rna")
  expect_identical(bio_alphabet(bio_seq_tbl(tibble::tibble(s = "MKLVE"), s)), "protein")
  expect_null(bio_alphabet(tibble::tibble()))
})

test_that("bio_seq_tbl() validates its input", {
  expect_error(bio_seq_tbl(1:3, s), "data frame")
  expect_error(bio_seq_tbl(tibble::tibble(s = 1:2), s), "character")
  expect_error(
    bio_seq_tbl(tibble::tibble(s = "ACGT1"), s, alphabet = "dna"),
    "not valid"
  )
  expect_error(bio_seq_tbl(tibble::tibble(s = "A"), s, alphabet = "xyz"))
  expect_error(new_bio_seq_tbl(tibble::tibble(a = 1), seq_col = "s"), "seq_col")
})

test_that("the class and metadata survive dplyr verbs", {
  x <- bio_seq_tbl(ecoli_cds, sequence)
  f <- dplyr::filter(x, class == "ribosomal")
  expect_s3_class(f, "bio_seq_tbl")
  expect_identical(bio_alphabet(f), "dna")
  expect_s3_class(dplyr::mutate(x, n = nchar(sequence)), "bio_seq_tbl")
  expect_s3_class(dplyr::arrange(x, gene), "bio_seq_tbl")
  expect_s3_class(x[1:3, ], "bio_seq_tbl")
})

test_that("dropping the sequence column falls back to a tibble", {
  x <- bio_seq_tbl(ecoli_cds, sequence)
  s <- dplyr::select(x, gene)
  expect_false(is_bio_seq_tbl(s))
  expect_s3_class(s, "tbl_df")
  expect_false(is_bio_seq_tbl(x[, "gene"]))
})

test_that("vctrs combines like with like", {
  x <- bio_seq_tbl(ecoli_cds[1:3, ], sequence)
  y <- bio_seq_tbl(ecoli_cds[4:5, ], sequence)
  both <- vctrs::vec_rbind(x, y)
  expect_s3_class(both, "bio_seq_tbl")
  expect_identical(nrow(both), 5L)
  mixed <- vctrs::vec_rbind(x, tibble::as_tibble(ecoli_cds[4:5, ]))
  expect_false(is_bio_seq_tbl(mixed))
  expect_s3_class(mixed, "tbl_df")
  p <- bio_seq_tbl(tibble::tibble(sequence = "MKV"), sequence)
  expect_false(is_bio_seq_tbl(vctrs::vec_rbind(x[, c("gene", "sequence")], p)))
  cast <- vctrs::vec_cast(tibble::as_tibble(ecoli_cds[1:2, ]), x)
  expect_s3_class(cast, "bio_seq_tbl")
  expect_s3_class(vctrs::vec_cast(x, tibble::as_tibble(ecoli_cds[0, ])), "tbl_df")
  expect_s3_class(vctrs::vec_ptype2(x, data.frame(gene = "a")), "tbl_df")
  expect_s3_class(vctrs::vec_ptype2(tibble::tibble(gene = "a"), x), "tbl_df")
  expect_s3_class(dplyr::bind_rows(x, y), "bio_seq_tbl")
})

test_that("assay and omics classes keep metadata", {
  a <- new_bio_assay(tibble::tibble(id = 1:3, v = 3:1), "custom", required = "id")
  expect_true(is_bio_assay(a))
  expect_s3_class(dplyr::filter(a, v > 1), "bio_assay")
  expect_false(is_bio_assay(dplyr::select(a, v)))
  expect_s3_class(vctrs::vec_rbind(a, a), "bio_assay")
  o <- new_bio_omics(tibble::tibble(feature = "g", p = 0.1), "custom", required = "feature")
  expect_true(is_bio_omics(o))
  expect_s3_class(vctrs::vec_rbind(o, o), "bio_omics")
  expect_false(is_bio_omics(vctrs::vec_rbind(o, tibble::tibble(feature = "h", p = 1))))
  expect_s3_class(vctrs::vec_cast(tibble::tibble(feature = "x", p = 1), o), "bio_omics")
  expect_error(new_bio_assay(tibble::tibble(), 1), "string")
  expect_error(new_bio_omics(tibble::tibble(), c("a", "b")), "string")
})

test_that("printing shows a header", {
  x <- bio_seq_tbl(ecoli_cds[1:2, ], sequence)
  expect_output(print(x), "bio_seq_tbl")
  expect_output(print(new_bio_assay(tibble::tibble(a = 1), "custom")), "bio_assay")
  expect_output(print(new_bio_omics(tibble::tibble(a = 1), "custom")), "bio_omics")
})
