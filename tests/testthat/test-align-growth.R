test_that("bio_align reproduces pwalign scores with BLOSUM62", {
  p <- ecoli_cds$protein[1:6]
  s <- ecoli_cds$protein[c(2, 3, 4, 5, 6, 1)]
  # pwalign::pairwiseAlignment(AAStringSet(p), AAStringSet(s),
  #   substitutionMatrix = "BLOSUM62", gapOpening = 10, gapExtension = 0.5)
  expect_equal(bio_align(p, s)$score, c(-0.5, -5, 5, -20, -20, 14.5))
  expect_equal(bio_align(p, s, type = "local")$score, c(62, 59, 38.5, 36.5, 42, 38))
})

test_that("bio_align reproduces pwalign scores for nucleotides", {
  d <- substr(ecoli_cds$sequence[1:5], 1, 150)
  e <- substr(ecoli_cds$sequence[c(2:5, 1)], 1, 140)
  # pwalign with nucleotideSubstitutionMatrix(match = 2, mismatch = -3),
  # gapOpening = 10 and gapExtension = 4
  expect_equal(bio_align(d, e)$score, c(-210, -220, -200, -210, -195))
  expect_equal(bio_align(d, e, type = "local")$score, c(16, 15, 14, 14, 19))
})

test_that("bio_align returns a consistent alignment", {
  res <- bio_align("HEAGAWGHEE", "PAWHEAE", type = "local")
  expect_equal(res$score, 17.5)
  expect_equal(res$aligned_pattern, "AWGHE")
  expect_equal(res$aligned_subject, "AW-HE")
  expect_equal(c(res$pattern_start, res$pattern_end), c(5L, 9L))
  expect_equal(c(res$subject_start, res$subject_end), c(2L, 5L))
  glob <- bio_align("ACGTTGCA", c("ACGTTGCA", "ACGTGCA"), match = 1, mismatch = -1,
                    gap_open = 1, gap_extend = 1)
  expect_equal(glob$identity[[1]], 1)
  expect_equal(gsub("-", "", glob$aligned_subject), c("ACGTTGCA", "ACGTGCA"))
  expect_equal(glob$score[[2]], 7 - 2)
  expect_error(bio_align("", "ACGT"))
  expect_error(bio_align(1, "ACGT"))
})

growth_sim <- function(model, p, t = 0:24, sd = 0.05, seed = 1) {
  set.seed(seed)
  y <- growth_value(model, p, t) + stats::rnorm(length(t), sd = sd)
  tibble::tibble(time = t, od = exp(y))
}

test_that("bio_growth_curve recovers known parameters", {
  truth <- c(y0 = -3, A = 4, mu_max = 0.6, lambda = 4, nu = 1)
  for (m in c("gompertz", "logistic", "baranyi")) {
    d <- growth_sim(m, truth)
    fit <- bio_growth_curve(d, time, od, model = m)
    expect_s3_class(fit, "bio_assay")
    expect_equal(fit$mu_max, 0.6, tolerance = 0.1)
    expect_equal(fit$lambda, 4, tolerance = 0.15)
    expect_equal(fit$A, 4, tolerance = 0.05)
    expect_equal(fit$doubling_time, log(2) / fit$mu_max)
  }
})

test_that("bio_growth_curve agrees with nls", {
  d <- fermentation_sim[fermentation_sim$reactor == "R1", ]
  fit <- bio_growth_curve(d, time_h, biomass_g_l, model = "gompertz")
  d$y <- log(d$biomass_g_l)
  ref <- stats::nls(y ~ y0 + A * exp(-exp(mu * exp(1) / A * (lam - time_h) + 1)), d,
                    start = list(y0 = -2.4, A = 3, mu = 0.5, lam = 0.5))
  td <- tidy(fit)
  expect_equal(td$estimate, unname(stats::coef(ref)), tolerance = 1e-4)
  expect_equal(td$std.error, unname(summary(ref)$coefficients[, 2]), tolerance = 1e-3)
  gl <- glance(fit)
  expect_equal(gl$rss, sum(stats::residuals(ref)^2), tolerance = 1e-6)
  aug <- augment(fit)
  expect_equal(aug$.resid, as.numeric(stats::residuals(ref)), tolerance = 1e-4)
})

test_that("bio_growth_curve handles groups, Richards and plotting", {
  fit <- bio_growth_curve(fermentation_sim, time_h, biomass_g_l, .by = reactor,
                          model = "richards")
  expect_equal(nrow(fit), 3L)
  expect_true("nu" %in% names(fit))
  expect_equal(nrow(tidy(fit)), 15L)
  expect_s3_class(autoplot(fit), "ggplot")
  expect_error(bio_growth_curve(fermentation_sim, time_h, glucose_g_l, model = "x"))
  bad <- fermentation_sim
  bad$biomass_g_l[1] <- 0
  expect_error(bio_growth_curve(bad, time_h, biomass_g_l, .by = reactor), "positive")
})

test_that("bio_qpcr_efficiency matches the standard curve formula", {
  std <- data.frame(
    gene = rep(c("a", "b"), each = 15),
    copies = rep(rep(10^(6:2), each = 3), 2),
    ct = c(
      15.1, 15.2, 15.0, 18.5, 18.4, 18.6, 21.8, 21.9, 21.7, 25.2, 25.1, 25.3, 28.6, 28.4, 28.7,
      16.3, 16.2, 16.4, 19.9, 19.8, 20.0, 23.4, 23.6, 23.5, 27.1, 27.0, 27.2, 30.7, 30.5, 30.8
    )
  )
  eff <- bio_qpcr_efficiency(std, ct, copies, .by = gene)
  ref <- stats::lm(ct ~ log10(copies), std[std$gene == "a", ])
  expect_equal(eff$slope[[1]], unname(stats::coef(ref)[2]))
  expect_equal(eff$amplification[[1]], 10^(-1 / unname(stats::coef(ref)[2])))
  expect_equal(eff$efficiency, eff$amplification - 1)
  expect_true(all(eff$efficiency.low < eff$efficiency & eff$efficiency < eff$efficiency.high))
  expect_equal(eff$range_log10, c(4, 4))
  expect_true(eff$acceptable[[1]])
  expect_false(eff$acceptable[[2]])
  expect_equal(nrow(tidy(eff)), 4L)
  expect_s3_class(autoplot(eff), "ggplot")
  expect_equal(nrow(augment(eff)), 30L)
  # Amplification factors feed straight into the Pfaffl model
  amp <- stats::setNames(c(1.95, 1.9), c("gapA", "rpoD"))
  expect_no_error(bio_qpcr_ddct(qpcr_sim[qpcr_sim$gene %in% names(amp), ], ct, gene, sample,
                                treatment, reference_genes = "gapA", control = "control",
                                efficiency = amp))
  expect_error(bio_qpcr_efficiency(data.frame(ct = 1:3, q = c(1, 1, 1)), ct, q))
})
