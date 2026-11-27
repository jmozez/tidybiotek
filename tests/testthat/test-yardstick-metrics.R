screen <- function() {
  tibble::tibble(
    active = factor(c("yes", "yes", "no", "no", "yes", rep("no", 15)), levels = c("yes", "no")),
    score = 20:1
  )
}

test_that("enrichment_factor() matches its definition", {
  d <- screen()
  ef <- enrichment_factor(d, active, score, fraction = 0.1)
  expect_s3_class(ef, "tbl_df")
  expect_identical(ef$.metric, "enrichment_factor")
  expect_identical(ef$.estimator, "binary")
  expect_equal(ef$.estimate, (2 / 2) / (3 / 20))
  expect_equal(enrichment_factor_vec(d$active, d$score, fraction = 0.25),
               (3 / 5) / (3 / 20))
  second <- enrichment_factor_vec(d$active, -d$score, fraction = 0.1, event_level = "second")
  expect_equal(second, (2 / 2) / (17 / 20))
  expect_identical(attr(enrichment_factor, "direction"), "maximize")
  na <- d$score
  na[1] <- NA
  expect_true(is.na(enrichment_factor_vec(d$active, na, na_rm = FALSE)))
  expect_error(enrichment_factor_vec(d$active, d$score, fraction = 2), "fraction")
  expect_error(enrichment_factor_vec(d$active, d$score, case_weights = rep(1, 20)), "case weights")
  grouped <- dplyr::group_by(dplyr::mutate(d, g = rep(1:2, 10)), g)
  expect_identical(nrow(enrichment_factor(grouped, active, score, fraction = 0.2)), 2L)
})

test_that("bedroc() is 1 for perfect and lower for poor rankings", {
  truth <- factor(rep(c("yes", "no"), c(5, 95)), levels = c("yes", "no"))
  perfect <- bedroc_vec(truth, 100:1)
  worst <- bedroc_vec(truth, 1:100)
  expect_equal(perfect, 1, tolerance = 1e-6)
  expect_lt(worst, 0.01)
  set.seed(5)
  random <- mean(replicate(200, bedroc_vec(truth, stats::runif(100))))
  expect_gt(random, worst)
  expect_lt(random, perfect)
  d <- screen()
  b <- bedroc(d, active, score, alpha = 20)
  expect_identical(b$.metric, "bedroc")
  expect_true(b$.estimate > 0 && b$.estimate <= 1)
  expect_warning(bedroc_vec(factor(rep("yes", 3), levels = c("yes", "no")), 1:3), "both")
  expect_error(bedroc_vec(truth, 100:1, alpha = -1), "positive")
})

test_that("metrics work in a metric set", {
  ms <- yardstick::metric_set(enrichment_factor, bedroc)
  out <- ms(screen(), active, score)
  expect_identical(out$.metric, c("enrichment_factor", "bedroc"))
})

test_that("assay_cv_error() is the RMS relative error in percent", {
  d <- tibble::tibble(t = c(1, 2, 4), e = c(1.1, 1.8, 4))
  out <- assay_cv_error(d, t, e)
  expect_identical(out$.metric, "assay_cv_error")
  expected <- 100 * sqrt(mean(c(0.1, -0.1, 0)^2))
  expect_equal(out$.estimate, expected)
  expect_equal(assay_cv_error_vec(d$t, d$e), expected)
  w <- assay_cv_error_vec(d$t, d$e, case_weights = c(1, 0, 0))
  expect_equal(w, 10)
  expect_equal(assay_cv_error_vec(c(1, NA), c(1.2, 3)), 20)
  expect_true(is.na(assay_cv_error_vec(c(1, NA), c(1.2, 3), na_rm = FALSE)))
  expect_error(assay_cv_error_vec(c(0, 1), c(1, 1)), "zeros")
  expect_identical(attr(assay_cv_error, "direction"), "minimize")
})
