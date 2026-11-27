test_that("bio_enzyme_kinetics() matches the classic Puromycin estimates", {
  fit <- bio_enzyme_kinetics(puromycin_kinetics, substrate, rate, .by = state)
  expect_s3_class(fit, "bio_assay")
  treated <- fit[fit$state == "treated", ]
  # Bates and Watts (1988): Vm = 212.7 (SE 6.95), K = 0.0641 (SE 0.00828)
  expect_equal(treated$vmax, 212.68, tolerance = 1e-3)
  expect_equal(treated$km, 0.06412, tolerance = 1e-3)
  expect_equal(treated$vmax_se, 6.947, tolerance = 1e-2)
  untreated <- fit[fit$state == "untreated", ]
  expect_equal(untreated$vmax, 160.28, tolerance = 1e-3)
  expect_equal(untreated$km, 0.04771, tolerance = 1e-3)
})

test_that("linearisations give exact answers on noise-free data", {
  s <- c(0.5, 1, 2, 4, 8, 16)
  d <- data.frame(s = s, v = 10 * s / (3 + s))
  for (m in c("lineweaver_burk", "eadie_hofstee", "hanes_woolf", "nls")) {
    fit <- suppressWarnings(bio_enzyme_kinetics(d, s, v, method = m))
    expect_equal(fit$vmax, 10, tolerance = 1e-6)
    expect_equal(fit$km, 3, tolerance = 1e-6)
  }
})

test_that("enzyme concentration gives kcat and kcat/km", {
  fit <- bio_enzyme_kinetics(puromycin_kinetics, substrate, rate,
    .by = state, enzyme_conc = 2
  )
  expect_equal(fit$kcat, fit$vmax / 2)
  expect_equal(fit$kcat_km, fit$vmax / 2 / fit$km)
  td <- tidy(fit)
  expect_setequal(unique(td$term), c("vmax", "km", "kcat", "kcat_km"))
  expect_error(
    bio_enzyme_kinetics(puromycin_kinetics, substrate, rate, enzyme_conc = -1),
    "positive"
  )
  expect_error(bio_enzyme_kinetics(data.frame(s = c(1, 1, 2), v = 1:3), s, v), "distinct")
})

test_that("monod_solve() conserves mass and matches the ODE", {
  times <- seq(0, 10, by = 0.5)
  out <- tidybiotek:::monod_solve(times, 0.5, 0.3, 0.5, 0.1, 20)
  total <- unname(out[, "biomass"] + 0.5 * out[, "substrate"])
  expect_equal(total, rep(0.1 + 0.5 * 20, length(times)), tolerance = 1e-8)
  expect_equal(unname(out[1, "biomass"]), 0.1, tolerance = 1e-8)
  # Euler check of dX/dt at t = 2
  h <- 1e-4
  xs <- tidybiotek:::monod_solve(c(2, 2 + h), 0.5, 0.3, 0.5, 0.1, 20)
  deriv <- unname((xs[2, 1] - xs[1, 1]) / h)
  s <- unname(xs[1, 2])
  expect_equal(deriv, 0.5 * s / (0.3 + s) * unname(xs[1, 1]), tolerance = 1e-3)
})

test_that("bio_monod_growth() recovers simulated parameters", {
  fit <- bio_monod_growth(fermentation_sim, time_h, biomass_g_l,
    substrate = glucose_g_l, .by = reactor
  )
  expect_s3_class(fit, "bio_assay")
  expect_identical(nrow(fit), 3L)
  expect_true(all(abs(fit$mu_max - 0.55) < 0.06))
  expect_true(all(abs(fit$yield - 0.48) < 0.05))
  expect_equal(fit$s0, c(10, 20, 30), tolerance = 0.05)
  expect_equal(fit$doubling_time, log(2) / fit$mu_max)
  td <- tidy(fit)
  expect_identical(nrow(td), 15L)
  aug <- augment(fit)
  expect_true(all(c(".fitted", ".resid", ".fitted_substrate") %in% names(aug)))
})

test_that("bio_monod_growth() works from biomass alone with s0", {
  r2 <- fermentation_sim[fermentation_sim$reactor == "R2", ]
  fit <- bio_monod_growth(r2, time_h, biomass_g_l, s0 = 20)
  expect_equal(fit$s0, 20)
  expect_true(abs(fit$mu_max - 0.55) < 0.1)
  expect_error(bio_monod_growth(r2, time_h, biomass_g_l), "s0")
  expect_error(bio_monod_growth(r2[1:3, ], time_h, biomass_g_l, s0 = 20), "5 distinct")
  bad <- r2
  bad$biomass_g_l[1] <- -1
  expect_error(bio_monod_growth(bad, time_h, biomass_g_l, s0 = 20), "positive")
})
