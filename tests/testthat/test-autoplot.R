test_that("autoplot() draws every assay type", {
  dr <- bio_dose_response(dose_response_sim, dose, response, .by = compound)
  expect_s3_class(autoplot(dr), "ggplot")
  run1 <- dnase_elisa[dnase_elisa$run == 1, ]
  run2 <- dnase_elisa[dnase_elisa$run == 2, ]
  plate <- rbind(transform(run1, known = conc), transform(run2, known = NA))
  el <- bio_elisa_quant(plate, density, known)
  expect_s3_class(autoplot(el), "ggplot")
  ek <- bio_enzyme_kinetics(puromycin_kinetics, substrate, rate, .by = state)
  expect_s3_class(autoplot(ek), "ggplot")
  mg <- bio_monod_growth(fermentation_sim[fermentation_sim$reactor == "R1", ], time_h,
    biomass_g_l,
    substrate = glucose_g_l
  )
  expect_s3_class(autoplot(mg), "ggplot")
  rq <- bio_qpcr_ddct(qpcr_sim, ct, gene, sample, treatment, c("gapA", "rpoD"), "control")
  expect_s3_class(autoplot(rq), "ggplot")
  for (p in list(autoplot(dr), autoplot(el), autoplot(ek), autoplot(mg), autoplot(rq))) {
    expect_no_error(ggplot2::ggplot_build(p))
  }
  expect_error(autoplot(new_bio_assay(tibble::tibble(a = 1), "custom")), "no plot")
})

test_that("autoplot() draws omics results", {
  de <- bio_diff_expr(expression_sim, gene, sample, count, condition,
    treatment = "stress", reference = "control", transform = "log2_cpm"
  )
  v <- autoplot(de)
  ma <- autoplot(de, type = "ma")
  expect_s3_class(v, "ggplot")
  expect_s3_class(ma, "ggplot")
  expect_identical(v$labels$x %||% "log2 fold change", "log2 fold change")
  ora <- bio_enrich(de, feature, p_adj < 0.05, gene_sets_sim)
  e <- autoplot(ora, n = 5)
  expect_s3_class(e, "ggplot")
  expect_identical(nrow(e$data), 5L)
  for (p in list(v, ma, e)) {
    expect_no_error(ggplot2::ggplot_build(p))
  }
})

test_that("broom methods return tidy tibbles", {
  dr <- bio_dose_response(dose_response_sim, dose, response, .by = compound)
  expect_identical(nrow(tidy(dr)), 12L)
  expect_identical(nrow(glance(dr)), 3L)
  aug <- augment(dr)
  expect_true(all(c(".fitted", ".resid") %in% names(aug)))
  expect_equal(aug$.resid, aug$response - aug$.fitted)
  one <- dplyr::filter(dr, compound == "inhibitor_a")
  expect_identical(nrow(tidy(one)), 4L)
  expect_identical(nrow(tidy(dr, conf.level = 0.9)), 12L)
  ek <- bio_enzyme_kinetics(puromycin_kinetics, substrate, rate, .by = state)
  expect_identical(nrow(augment(ek)), nrow(puromycin_kinetics))
  expect_identical(nrow(glance(ek)), 2L)
  run1 <- dnase_elisa[dnase_elisa$run == 1, ]
  el <- bio_elisa_quant(run1, density, conc)
  expect_true(".fitted" %in% names(augment(el)))
  expect_identical(nrow(glance(el)), 1L)
  rq <- bio_qpcr_ddct(qpcr_sim, ct, gene, sample, treatment, c("gapA", "rpoD"), "control")
  expect_identical(glance(rq)$n_targets, 4L)
  expect_identical(nrow(augment(rq)), nrow(rq))
  de <- bio_diff_expr(expression_sim, gene, sample, count, condition,
    treatment = "stress", reference = "control", transform = "log2_cpm"
  )
  expect_named(tidy(de), c("feature", "estimate", "std.error", "statistic", "p.value",
                           "p.adjusted", "ave_expr"))
  expect_true(glance(de)$n_significant > 50)
  aug_de <- augment(de, lfc = 1)
  expect_true(all(aug_de$.direction %in% c("up", "down", "none")))
  ora <- bio_enrich(de, feature, p_adj < 0.05, gene_sets_sim)
  expect_identical(nrow(tidy(ora)), nrow(ora))
  expect_identical(glance(ora)$universe_size, 1000L)
  expect_true(is.logical(augment(ora)$.significant))
  other <- new_bio_omics(tibble::tibble(a = 1), "custom")
  expect_identical(nrow(tidy(other)), 1L)
  expect_identical(glance(other)$nrow, 1L)
  custom <- new_bio_assay(tibble::tibble(a = 1:2), "custom")
  expect_identical(nrow(tidy(custom)), 2L)
  expect_identical(glance(custom)$nrow, 2L)
  expect_identical(nrow(augment(custom)), 2L)
})
