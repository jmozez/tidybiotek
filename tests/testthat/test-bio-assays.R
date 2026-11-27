test_that("bio_dose_response() recovers simulated IC50 and EC50 values", {
  fit <- bio_dose_response(dose_response_sim, dose, response, .by = compound)
  expect_s3_class(fit, "bio_assay")
  expect_identical(nrow(fit), 3L)
  truth <- c(activator_c = 1.2, inhibitor_a = 0.05, inhibitor_b = 2.5)
  est <- stats::setNames(fit$ec50, fit$compound)[names(truth)]
  expect_true(all(abs(log(est / truth)) < log(1.3)))
  expect_identical(fit$type[fit$compound == "inhibitor_a"], "IC50")
  expect_identical(fit$type[fit$compound == "activator_c"], "EC50")
  expect_true(all(fit$bottom < fit$top))
  expect_true(all(fit$ec50_low < fit$ec50 & fit$ec50 < fit$ec50_high))
  expect_true(all(fit$converged))
})

test_that("bio_dose_response() fits exact curves and fixed asymptotes", {
  x <- c(0, 0.01, 0.1, 0.3, 1, 3, 10, 100)
  y <- 2 + (90 - 2) / (1 + (0.5 / x)^1.5)
  y[1] <- 2
  fit <- bio_dose_response(data.frame(x = x, y = y), x, y)
  expect_equal(fit$ec50, 0.5, tolerance = 1e-4)
  expect_equal(fit$hill, 1.5, tolerance = 1e-4)
  expect_equal(fit$bottom, 2, tolerance = 1e-4)
  fixed <- bio_dose_response(data.frame(x = x, y = y), x, y, bottom = 2)
  expect_equal(fixed$ec50, 0.5, tolerance = 1e-4)
  expect_true(is.na(tidy(fixed)$std.error[[1]]))
  expect_error(bio_dose_response(data.frame(x = c(1, 1, 2), y = 1:3), x, y), "distinct")
  expect_error(bio_dose_response(data.frame(x = -1:5, y = 1:7), x, y), "positive")
  expect_error(bio_dose_response(data.frame(x = "a", y = 1), x, y), "numeric")
})

test_that("bio_elisa_quant() back-calculates its own standards", {
  run1 <- dnase_elisa[dnase_elisa$run == 1, ]
  run1$known <- run1$conc
  q <- bio_elisa_quant(run1, density, known)
  expect_s3_class(q, "bio_assay")
  expect_true(all(q$is_standard))
  ok <- q$in_range
  expect_gt(mean(ok), 0.8)
  rel <- abs(q$estimated_conc[ok] - q$conc[ok]) / q$conc[ok]
  expect_lt(stats::median(rel), 0.1)
  q2 <- bio_elisa_quant(dnase_elisa, density, conc, .by = run)
  expect_identical(length(attr(q2, "bio_meta")$fits), 11L)
  dil <- run1
  dil$dil <- 2
  qd <- bio_elisa_quant(dil, density, known, dilution = dil)
  expect_equal(qd$estimated_conc, 2 * q$estimated_conc)
  expect_error(bio_elisa_quant(run1[1:3, ], density, known), "at least 5")
})

test_that("four_pl_inverse() inverts the curve", {
  par <- c(bottom = 1, top = 3, log_ec50 = log(2), hill = 1.2)
  x <- c(0.5, 1, 2, 5)
  y <- tidybiotek:::four_pl_value(par, x)
  expect_equal(tidybiotek:::four_pl_inverse(par, y), x)
  expect_true(is.na(tidybiotek:::four_pl_inverse(par, 5)))
})

test_that("bio_qpcr_ddct() recovers simulated fold changes", {
  rq <- bio_qpcr_ddct(qpcr_sim, ct, gene, sample, treatment,
    reference_genes = c("gapA", "rpoD"), control = "control"
  )
  expect_s3_class(rq, "bio_assay")
  expect_identical(nrow(rq), 12L * 4L)
  expect_equal(rq$fold_change, 2^-rq$delta_delta_ct)
  s <- tidy(rq)
  expect_equal(s$fold_change[s$treatment == "control"], rep(1, 4))
  ibpa <- s$log2_fold_change[s$gene == "ibpA" & s$treatment == "heat_46C"]
  expect_equal(ibpa, 6, tolerance = 0.1)
  ompf <- s$log2_fold_change[s$gene == "ompF" & s$treatment == "heat_42C"]
  expect_lt(ompf, 0)
  expect_true(all(s$p.value[s$treatment != "control"] < 0.05))
})

test_that("bio_qpcr_ddct() supports Pfaffl efficiencies and validates input", {
  eff <- c(gapA = 1.95, rpoD = 1.9, dnaK = 1.98, groL = 2, ibpA = 1.92, ompF = 1.97)
  rq <- bio_qpcr_ddct(qpcr_sim, ct, gene, sample, treatment,
    reference_genes = "gapA", control = "control", efficiency = eff
  )
  ctrl <- rq[rq$treatment == "control" & rq$gene == "dnaK", ]
  expect_equal(exp(mean(log(ctrl$fold_change))), 1)
  expect_false(isTRUE(all.equal(rq$fold_change, 2^-rq$delta_delta_ct)))
  expect_error(
    bio_qpcr_ddct(qpcr_sim, ct, gene, sample, treatment, "nope", "control"),
    "Not found"
  )
  expect_error(
    bio_qpcr_ddct(qpcr_sim, ct, gene, sample, treatment, "gapA", "nope"),
    "control"
  )
  expect_error(
    bio_qpcr_ddct(qpcr_sim, ct, gene, sample, treatment, "gapA", "control", efficiency = 3),
    "between 1 and 2"
  )
  expect_error(
    bio_qpcr_ddct(qpcr_sim, ct, gene, sample, treatment, "gapA", "control",
      efficiency = c(gapA = 2)
    ),
    "named"
  )
  only_refs <- qpcr_sim[qpcr_sim$gene %in% c("gapA", "rpoD"), ]
  expect_error(
    bio_qpcr_ddct(only_refs, ct, gene, sample, treatment, c("gapA", "rpoD"), "control"),
    "no target"
  )
})
