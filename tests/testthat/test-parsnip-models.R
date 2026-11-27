one_compound <- function() {
  dose_response_sim[dose_response_sim$compound == "inhibitor_b", ]
}

test_that("bio_four_pl() creates a parsnip specification", {
  spec <- bio_four_pl()
  expect_s3_class(spec, c("bio_four_pl", "model_spec"))
  expect_identical(spec$mode, "regression")
  expect_identical(spec$engine, "nls_4pl")
  expect_output(print(spec), "Four-parameter logistic")
  expect_error(bio_four_pl(mode = "classification"), "regression")
})

test_that("bio_four_pl() fits and predicts through parsnip", {
  d <- one_compound()
  fit <- parsnip::fit(bio_four_pl(), response ~ dose, data = d)
  expect_s3_class(fit$fit, "bio_four_pl_fit")
  pred <- predict(fit, new_data = data.frame(dose = c(0.01, 2.5, 100)))
  expect_named(pred, ".pred")
  expect_true(pred$.pred[[1]] > pred$.pred[[3]])
  direct <- bio_dose_response(d, dose, response)
  td <- tidy(fit$fit)
  expect_equal(td$estimate[td$term == "ec50"], direct$ec50, tolerance = 1e-6)
  expect_identical(nrow(glance(fit$fit)), 1L)
  expect_output(print(fit$fit), "Four-parameter logistic fit")
})

test_that("engine arguments fix asymptotes", {
  d <- one_compound()
  spec <- parsnip::set_engine(bio_four_pl(), "nls_4pl", bottom = 0)
  fit <- parsnip::fit(spec, response ~ dose, data = d)
  expect_equal(fit$fit$fit$coef[["bottom"]], 0)
  expect_error(bio_four_pl_fit(response ~ dose + replicate, d), "exactly one")
  expect_error(bio_four_pl_fit(compound ~ dose, d), "numeric outcome")
})

test_that("bio_four_pl() works in a workflow-free resampling loop", {
  d <- one_compound()
  set.seed(9)
  folds <- rsample::vfold_cv(d, v = 3)
  rmse <- vapply(folds$splits, function(sp) {
    fit <- parsnip::fit(bio_four_pl(), response ~ dose, data = rsample::analysis(sp))
    a <- rsample::assessment(sp)
    yardstick::rmse_vec(a$response, predict(fit, a)$.pred)
  }, numeric(1))
  expect_true(all(rmse < 10))
})
