# Primary growth models and qPCR standard curves.

growth_names <- c("y0", "A", "mu_max", "lambda", "nu")

# Zwietering et al. (1990) reparameterised curves and Baranyi and Roberts
# (1994), all on the log scale: y(t) = log N(t).
growth_value <- function(model, p, t) {
  y0 <- p[["y0"]]
  a <- p[["A"]]
  mu <- p[["mu_max"]]
  lam <- p[["lambda"]]
  switch(model,
    gompertz = y0 + a * exp(-exp(mu * exp(1) / a * (lam - t) + 1)),
    logistic = y0 + a / (1 + exp(4 * mu / a * (lam - t) + 2)),
    richards = {
      nu <- p[["nu"]]
      y0 + a * (1 + nu * exp(1 + nu) * exp(mu / a * (1 + nu)^(1 + 1 / nu) * (lam - t)))^(-1 / nu)
    },
    baranyi = {
      h0 <- mu * lam
      at <- t + log(exp(-mu * t) + exp(-h0) - exp(-mu * t - h0)) / mu
      y0 + mu * at - log(1 + (exp(mu * at) - 1) / exp(a))
    }
  )
}

growth_start <- function(t, y) {
  y0 <- min(y[t == min(t)])
  a <- max(max(y) - y0, 1e-3)
  ts <- sort(unique(t))
  ym <- vapply(ts, function(s) mean(y[t == s]), numeric(1))
  sl <- diff(ym) / diff(ts)
  k <- which.max(sl)
  mu <- max(sl[[k]], a / diff(range(t)) / 4, 1e-4)
  lam <- max(ts[[k]] - (ym[[k]] - y0) / mu, 0)
  c(y0 = y0, A = a, mu_max = mu, lambda = lam)
}

fit_growth <- function(time, y, model) {
  ord <- order(time)
  time <- time[ord]
  y <- y[ord]
  if (length(unique(time)) < 5L) {
    cli::cli_abort("Growth curve fits need at least 5 distinct time points.")
  }
  st <- growth_start(time, y)
  rich <- model == "richards"
  to_p <- function(q) {
    p <- c(y0 = q[[1]], A = exp(q[[2]]), mu_max = exp(q[[3]]), lambda = q[[4]],
           nu = if (rich) exp(q[[5]]) else 1)
    p
  }
  rss <- function(q) {
    r <- y - growth_value(model, to_p(q), time)
    s <- sum(r^2)
    if (is.finite(s)) s else 1e100
  }
  best <- NULL
  for (lf in c(1, 0.5, 0)) {
    q0 <- c(st[["y0"]], log(st[["A"]]), log(st[["mu_max"]]), st[["lambda"]] * lf)
    if (rich) q0 <- c(q0, 0)
    res <- stats::optim(q0, rss, method = "Nelder-Mead",
                        control = list(maxit = 5000, reltol = 1e-12))
    res <- stats::optim(res$par, rss, method = "BFGS", control = list(maxit = 1000, reltol = 1e-14))
    if (is.null(best) || res$value < best$value) best <- res
  }
  p <- to_p(best$par)
  k <- if (rich) 5L else 4L
  free_names <- growth_names[seq_len(k)]
  n <- length(y)
  df <- n - k
  if (df < 1L) cli::cli_abort("Too few observations to fit the {.val {model}} model.")
  # Standard errors from the Jacobian on the natural parameter scale.
  jac <- vapply(free_names, function(nm) {
    h <- 1e-6 * max(abs(p[[nm]]), 1e-3)
    up <- p
    dn <- p
    up[[nm]] <- p[[nm]] + h
    dn[[nm]] <- p[[nm]] - h
    (growth_value(model, up, time) - growth_value(model, dn, time)) / (2 * h)
  }, numeric(n))
  jac <- matrix(jac, nrow = n)
  sigma2 <- best$value / df
  vc <- tryCatch(sigma2 * solve(crossprod(jac)), error = function(e) NULL)
  se <- stats::setNames(rep(NA_real_, 5), growth_names)
  if (!is.null(vc)) se[free_names] <- sqrt(pmax(diag(vc), 0))
  if (!rich) p[["nu"]] <- NA_real_
  fitted <- growth_value(model, p, time)
  list(
    coef = p, se = se, free = growth_names %in% free_names, sigma = sqrt(sigma2), df = df,
    rss = best$value, nobs = n, converged = best$convergence == 0, x = time, y = y,
    fitted = fitted, order = ord, model = model
  )
}

#' Primary microbial growth curves
#'
#' @description
#' Fits the primary growth models of predictive microbiology to log
#' population size (cell counts, colony-forming units or optical density)
#' over time, and reports the lag time, maximum specific growth rate and
#' asymptotic increase for each culture:
#'
#' * `"gompertz"` and `"logistic"`: the modified forms of Zwietering et al.
#'   (1990), parameterised directly by \eqn{\mu_{max}}{mu_max},
#'   \eqn{\lambda}{lambda} and \eqn{A};
#' * `"richards"`: the Zwietering form of the Richards curve with an extra
#'   shape parameter \eqn{\nu}{nu};
#' * `"baranyi"`: the Baranyi and Roberts (1994) model, the current
#'   standard in food microbiology (for example ComBase and DMFit).
#'
#' All models are fitted by least squares to \eqn{y = \log N}{y = log N}
#' (natural logarithm), with the initial level \eqn{y_0}{y0} estimated, so
#' `mu_max` is the specific growth rate in reciprocal time units and
#' `doubling_time` is `log(2) / mu_max`. Several starting values are tried
#' and the best fit is kept. Model choice can be made with the AIC from
#' [glance()][glance.bio_assay].
#'
#' @param .data A data frame or tibble.
#' @param time,value <[`tidy-select`][dplyr::dplyr_tidy_select]> Time and
#'   population size columns.
#' @param .by <[`tidy-select`][dplyr::dplyr_tidy_select]> Optional columns
#'   identifying separate cultures or runs.
#' @param model One of `"gompertz"`, `"logistic"`, `"richards"` or
#'   `"baranyi"`.
#' @param transform `"log"` (the default) takes natural logarithms of
#'   `value`; `"none"` uses `value` as it is, for data already on a log
#'   scale (the parameters are then in those units).
#' @return A `<bio_assay>` with one row per culture: the `.by` columns,
#'   `y0`, `A` (the increase in log population), `mu_max`, `lambda` (lag
#'   time), `nu` (Richards only), `doubling_time`, `sigma`, `nobs` and
#'   `converged`. Standard errors and intervals are in
#'   [tidy()][tidy.bio_assay], fitted values in [augment()][augment.bio_assay]
#'   and the curves in [autoplot()][autoplot.bio_assay].
#' @references Zwietering, M. H., Jongenburger, I., Rombouts, F. M. and van
#'   't Riet, K. (1990). Modeling of the bacterial growth curve. Applied and
#'   Environmental Microbiology, 56, 1875-1881.
#'   \doi{10.1128/aem.56.6.1875-1881.1990}
#'
#'   Baranyi, J. and Roberts, T. A. (1994). A dynamic approach to predicting
#'   bacterial growth in food. International Journal of Food Microbiology,
#'   23, 277-294. \doi{10.1016/0168-1605(94)90157-0}
#' @export
#' @examples
#' growth <- bio_growth_curve(fermentation_sim, time_h, biomass_g_l,
#'   .by = reactor, model = "baranyi"
#' )
#' growth
#' tidy(growth)
#' glance(growth)
bio_growth_curve <- function(.data, time, value, .by = NULL,
                             model = c("gompertz", "logistic", "richards", "baranyi"),
                             transform = c("log", "none")) {
  check_data_frame(.data)
  model <- rlang::arg_match(model)
  transform <- rlang::arg_match(transform)
  t_col <- select_one(.data, {{ time }}, "time")
  v_col <- select_one(.data, {{ value }}, "value")
  by <- select_by(.data, {{ .by }})
  check_numeric_col(.data, t_col, "time")
  check_numeric_col(.data, v_col, "value")
  if (transform == "log" && any(.data[[v_col]] <= 0, na.rm = TRUE)) {
    cli::cli_abort("{.arg value} must be positive when {.code transform = \"log\"}.")
  }
  groups <- split_rows(.data, by)
  fits <- lapply(seq_along(groups$rows), function(i) {
    rows <- groups$rows[[i]]
    t <- .data[[t_col]][rows]
    v <- .data[[v_col]][rows]
    ok <- !is.na(t) & !is.na(v)
    rows <- rows[ok]
    y <- if (transform == "log") log(v[ok]) else v[ok]
    fit <- fit_growth(t[ok], y, model)
    fit$rows <- rows[fit$order]
    fit$key <- groups$key[i, ]
    fit
  })
  rows <- lapply(fits, function(fit) {
    p <- fit$coef
    tibble::tibble(
      y0 = p[["y0"]], A = p[["A"]], mu_max = p[["mu_max"]], lambda = p[["lambda"]],
      nu = p[["nu"]], doubling_time = log(2) / p[["mu_max"]], sigma = fit$sigma,
      nobs = fit$nobs, converged = fit$converged
    )
  })
  out <- bind_key(groups$key, rows)
  if (model != "richards") out$nu <- NULL
  meta <- list(fits = fits, by = by, x = t_col, y = v_col, model = model,
               transform = transform, data = tibble::as_tibble(.data))
  new_bio_assay(out, "growth_curve", meta, required = c(by, "mu_max"))
}

growth_params <- function(fit, level = 0.95) {
  tq <- t_quantile(level, fit$df)
  keep <- fit$free
  est <- fit$coef[keep]
  se <- fit$se[keep]
  tibble::tibble(
    term = growth_names[keep], estimate = unname(est), std.error = unname(se),
    statistic = estimate / std.error,
    p.value = 2 * stats::pt(-abs(statistic), df = fit$df),
    conf.low = estimate - tq * std.error, conf.high = estimate + tq * std.error
  )
}

#' qPCR standard curves and amplification efficiency
#'
#' @description
#' Fits the standard curve \eqn{C_q = a + b \log_{10}(\mathrm{quantity})}{Cq
#' = a + b log10(quantity)} to a dilution series for each assay and reports
#' the amplification factor \eqn{E = 10^{-1/b}}{E = 10^(-1/b)}, the
#' efficiency \eqn{E - 1}, the linearity (\eqn{R^2}) and the dynamic range,
#' as required by the MIQE guidelines (Bustin et al. 2009). An efficiency
#' between 90 and 110 per cent with \eqn{R^2 \ge 0.98}{R2 >= 0.98} is
#' flagged as acceptable. The amplification factors can be passed as the
#' `efficiency` of [bio_qpcr_ddct()] for the Pfaffl (2001) model.
#'
#' @param .data A data frame or tibble with one row per standard reaction.
#' @param ct,quantity <[`tidy-select`][dplyr::dplyr_tidy_select]>
#'   Quantification cycle and known template quantity (copies,
#'   concentration or relative dilution) columns.
#' @param .by <[`tidy-select`][dplyr::dplyr_tidy_select]> Optional columns
#'   identifying assays (usually the gene).
#' @param conf.level Confidence level for the efficiency interval.
#' @return A `<bio_assay>` with one row per assay: the `.by` columns,
#'   `slope`, `intercept`, `r_squared`, `amplification` (\eqn{E}),
#'   `efficiency` (\eqn{E - 1}), `efficiency.low` and `efficiency.high`
#'   (from the slope interval), `range_log10` (decades covered), `nobs` and
#'   `acceptable`.
#' @references Bustin, S. A. et al. (2009). The MIQE guidelines: minimum
#'   information for publication of quantitative real-time PCR experiments.
#'   Clinical Chemistry, 55, 611-622. \doi{10.1373/clinchem.2008.112797}
#'
#'   Pfaffl, M. W. (2001). A new mathematical model for relative
#'   quantification in real-time RT-PCR. Nucleic Acids Research, 29, e45.
#'   \doi{10.1093/nar/29.9.e45}
#' @export
#' @examples
#' # A ten-fold dilution series run in triplicate for two genes
#' std <- data.frame(
#'   gene = rep(c("gapA", "rpoD"), each = 15),
#'   copies = rep(rep(10^(6:2), each = 3), 2),
#'   ct = c(
#'     15.1, 15.2, 15.0, 18.5, 18.4, 18.6, 21.8, 21.9, 21.7,
#'     25.2, 25.1, 25.3, 28.6, 28.4, 28.7,
#'     16.3, 16.2, 16.4, 19.9, 19.8, 20.0, 23.4, 23.6, 23.5,
#'     27.1, 27.0, 27.2, 30.7, 30.5, 30.8
#'   )
#' )
#' eff <- bio_qpcr_efficiency(std, ct, copies, .by = gene)
#' eff
#' stats::setNames(eff$amplification, eff$gene)
bio_qpcr_efficiency <- function(.data, ct, quantity, .by = NULL, conf.level = 0.95) { # nolint
  check_data_frame(.data)
  ct_col <- select_one(.data, {{ ct }}, "ct")
  q_col <- select_one(.data, {{ quantity }}, "quantity")
  by <- select_by(.data, {{ .by }})
  check_numeric_col(.data, ct_col, "ct")
  check_numeric_col(.data, q_col, "quantity")
  if (any(.data[[q_col]] <= 0, na.rm = TRUE)) {
    cli::cli_abort("{.arg quantity} must be positive.")
  }
  groups <- split_rows(.data, by)
  fits <- lapply(seq_along(groups$rows), function(i) {
    rows <- groups$rows[[i]]
    q <- .data[[q_col]][rows]
    y <- .data[[ct_col]][rows]
    ok <- !is.na(q) & !is.na(y)
    rows <- rows[ok]
    q <- q[ok]
    y <- y[ok]
    if (length(unique(q)) < 3L) {
      cli::cli_abort("A standard curve needs at least 3 distinct quantities.")
    }
    lx <- log10(q)
    fit <- stats::lm(y ~ lx)
    cf <- stats::coef(fit)
    se <- sqrt(diag(stats::vcov(fit)))
    rss <- sum(stats::residuals(fit)^2)
    df <- fit$df.residual
    list(
      coef = c(intercept = cf[[1]], slope = cf[[2]]),
      se = c(intercept = se[[1]], slope = se[[2]]), free = c(TRUE, TRUE),
      sigma = sqrt(rss / df), df = df, rss = rss, nobs = length(y), converged = TRUE,
      x = q, y = y, fitted = unname(stats::fitted(fit)), rows = rows,
      r_squared = summary(fit)$r.squared, key = groups$key[i, ]
    )
  })
  rows <- lapply(fits, function(f) {
    b <- f$coef[["slope"]]
    tq <- t_quantile(conf.level, f$df)
    blim <- b + c(-1, 1) * tq * f$se[["slope"]]
    eff_lim <- sort(10^(-1 / blim) - 1)
    if (any(blim >= 0)) eff_lim <- c(NA_real_, NA_real_)
    amp <- 10^(-1 / b)
    tibble::tibble(
      slope = b, intercept = f$coef[["intercept"]], r_squared = f$r_squared,
      amplification = amp, efficiency = amp - 1,
      efficiency.low = eff_lim[[1]], efficiency.high = eff_lim[[2]],
      range_log10 = diff(range(log10(f$x))), nobs = f$nobs,
      acceptable = b < 0 && amp - 1 >= 0.9 && amp - 1 <= 1.1 && f$r_squared >= 0.98
    )
  })
  out <- bind_key(groups$key, rows)
  meta <- list(fits = fits, by = by, x = q_col, y = ct_col, level = conf.level,
               data = tibble::as_tibble(.data))
  new_bio_assay(out, "qpcr_standard", meta, required = c(by, "slope"))
}

qpcr_standard_params <- function(fit, level = 0.95) {
  tq <- t_quantile(level, fit$df)
  tibble::tibble(
    term = c("intercept", "slope"), estimate = unname(fit$coef), std.error = unname(fit$se),
    statistic = estimate / std.error,
    p.value = 2 * stats::pt(-abs(statistic), df = fit$df),
    conf.low = estimate - tq * std.error, conf.high = estimate + tq * std.error
  )
}

plot_growth <- function(object, meta) {
  fits <- active_fits(object)
  pts <- point_data(fits)
  curves <- curve_grid(fits, function(f, g) growth_value(f$model, f$coef, g), log_x = FALSE)
  ylab <- if (meta$transform == "log") paste0("log(", meta$y, ")") else meta$y
  ggplot2::ggplot(pts, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$.curve)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(data = curves, linewidth = 0.8) +
    ggplot2::labs(x = meta$x, y = ylab, colour = NULL,
                  title = paste0("Growth curves (", meta$model, " model)")) +
    ggplot2::theme_minimal()
}

plot_qpcr_standard <- function(object, meta) {
  fits <- active_fits(object)
  pts <- point_data(fits)
  curves <- curve_grid(fits, function(f, g) f$coef[["intercept"]] + f$coef[["slope"]] * log10(g))
  ggplot2::ggplot(pts, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$.curve)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(data = curves, linewidth = 0.8) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = meta$x, y = meta$y, colour = NULL, title = "qPCR standard curves") +
    ggplot2::theme_minimal()
}
