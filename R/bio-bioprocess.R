# Bioprocess kinetics: Michaelis-Menten enzyme kinetics and Monod growth in
# batch bioreactors.

mm_value <- function(par, s) par[[1]] * s / (par[[2]] + s)

#' Michaelis-Menten enzyme kinetics
#'
#' Estimates the maximum rate `vmax` and the Michaelis constant `km` from
#' initial rates measured at several substrate concentrations. The default
#' fits \eqn{v = V_{max} S / (K_m + S)}{v = Vmax S / (Km + S)} by nonlinear
#' least squares with `stats::nls()`, started from the Hanes-Woolf
#' estimates. The classic linearisations (Lineweaver-Burk, Eadie-Hofstee and
#' Hanes-Woolf) are available for teaching and comparison; their standard
#' errors come from the delta method. Give `enzyme_conc` to also get the
#' turnover number `kcat` and the specificity constant `kcat / km`.
#'
#' @param .data A data frame or tibble.
#' @param substrate,rate <[`tidy-select`][dplyr::dplyr_tidy_select]>
#'   Substrate concentration and initial rate columns.
#' @param .by <[`tidy-select`][dplyr::dplyr_tidy_select]> Optional columns
#'   defining separate enzymes or conditions.
#' @param method One of `"nls"`, `"lineweaver_burk"`, `"eadie_hofstee"` or
#'   `"hanes_woolf"`.
#' @param enzyme_conc Optional total enzyme concentration, in the same
#'   amount units as `vmax` per unit time.
#' @param level Confidence level.
#'
#' @return A `<bio_assay>` with one row per group: the `.by` columns,
#'   `vmax`, `vmax_se`, `km`, `km_se`, `kcat` and `kcat_km` (when
#'   `enzyme_conc` is given), `sigma`, `nobs` and `method`.
#' @references Michaelis, L. and Menten, M. L. (1913). Die Kinetik der
#'   Invertinwirkung. *Biochemische Zeitschrift*, 49, 333-369.
#'
#'   Bates, D. M. and Watts, D. G. (1988). *Nonlinear Regression Analysis and
#'   Its Applications*. Wiley. \doi{10.1002/9780470316757}
#' @export
#' @examples
#' bio_enzyme_kinetics(puromycin_kinetics, substrate, rate, .by = state)
#'
#' # Compare with the Lineweaver-Burk linearisation
#' bio_enzyme_kinetics(puromycin_kinetics, substrate, rate,
#'   .by = state,
#'   method = "lineweaver_burk"
#' )
bio_enzyme_kinetics <- function(.data, substrate, rate, .by = NULL,
                                method = c("nls", "lineweaver_burk", "eadie_hofstee",
                                           "hanes_woolf"),
                                enzyme_conc = NULL, level = 0.95) {
  check_data_frame(.data)
  s_col <- select_one(.data, {{ substrate }}, "substrate")
  v_col <- select_one(.data, {{ rate }}, "rate")
  by <- select_by(.data, {{ .by }})
  method <- rlang::arg_match(method)
  check_numeric_col(.data, s_col, "substrate")
  check_numeric_col(.data, v_col, "rate")
  if (!is.null(enzyme_conc)) {
    check_positive_number(enzyme_conc)
  }
  groups <- split_rows(.data, by)
  fits <- lapply(seq_along(groups$rows), function(i) {
    rows <- groups$rows[[i]]
    s <- .data[[s_col]][rows]
    v <- .data[[v_col]][rows]
    ok <- is.finite(s) & is.finite(v)
    fit <- fit_mm(s[ok], v[ok], method)
    fit$rows <- rows[ok]
    fit$key <- groups$key[i, ]
    fit$params <- mm_params(fit, level, enzyme_conc)
    fit
  })
  rows <- lapply(fits, function(fit) {
    p <- fit$params
    out <- tibble::tibble(
      vmax = p$estimate[[1]], vmax_se = p$std.error[[1]],
      km = p$estimate[[2]], km_se = p$std.error[[2]]
    )
    if (!is.null(enzyme_conc)) {
      out$kcat <- p$estimate[[3]]
      out$kcat_km <- p$estimate[[4]]
    }
    out$sigma <- fit$sigma
    out$nobs <- fit$nobs
    out$method <- method
    out
  })
  out <- bind_key(groups$key, rows)
  meta <- list(fits = fits, by = by, x = s_col, y = v_col, level = level, model = "mm",
               method = method, enzyme_conc = enzyme_conc, data = tibble::as_tibble(.data))
  new_bio_assay(out, "enzyme_kinetics", meta, required = c(by, "vmax", "km"))
}

fit_mm <- function(s, v, method) {
  if (length(unique(s)) < 3L) {
    cli::cli_abort("Michaelis-Menten fits need at least 3 distinct substrate concentrations.")
  }
  if (any(s <= 0) && method != "nls") {
    cli::cli_abort("Linearised methods need positive substrate concentrations.")
  }
  hw <- stats::lm(I(s / v) ~ s, subset = s > 0 & v > 0)
  hw_coef <- stats::coef(hw)
  start <- c(vmax = 1 / hw_coef[[2]], km = hw_coef[[1]] / hw_coef[[2]])
  if (!all(is.finite(start)) || any(start <= 0)) {
    start <- c(vmax = max(v) * 1.1, km = stats::median(s))
  }
  n <- length(s)
  if (method == "nls") {
    dat <- data.frame(s = s, v = v)
    model <- tryCatch(
      stats::nls(v ~ vmax * s / (km + s), data = dat, start = as.list(start),
                 control = stats::nls.control(maxiter = 200, scaleOffset = 1)),
      error = function(e) {
        cli::cli_abort(c("The Michaelis-Menten fit failed.", "x" = conditionMessage(e)),
                       call = NULL)
      }
    )
    coef <- stats::coef(model)
    vc <- stats::vcov(model)
    df <- n - 2
  } else {
    lin <- switch(method,
      lineweaver_burk = stats::lm(I(1 / v) ~ I(1 / s)),
      eadie_hofstee = stats::lm(v ~ I(v / s)),
      hanes_woolf = stats::lm(I(s / v) ~ s)
    )
    b <- stats::coef(lin)
    lvc <- stats::vcov(lin)
    # Map linear coefficients (a, b) to (vmax, km) and the delta-method Jacobian.
    mapped <- switch(method,
      lineweaver_burk = list(
        coef = c(1 / b[[1]], b[[2]] / b[[1]]),
        jac = rbind(c(-1 / b[[1]]^2, 0), c(-b[[2]] / b[[1]]^2, 1 / b[[1]]))
      ),
      eadie_hofstee = list(
        coef = c(b[[1]], -b[[2]]),
        jac = rbind(c(1, 0), c(0, -1))
      ),
      hanes_woolf = list(
        coef = c(1 / b[[2]], b[[1]] / b[[2]]),
        jac = rbind(c(0, -1 / b[[2]]^2), c(1 / b[[2]], -b[[1]] / b[[2]]^2))
      )
    )
    coef <- stats::setNames(mapped$coef, c("vmax", "km"))
    vc <- mapped$jac %*% lvc %*% t(mapped$jac)
    dimnames(vc) <- list(c("vmax", "km"), c("vmax", "km"))
    df <- n - 2
  }
  fitted <- mm_value(coef, s)
  rss <- sum((v - fitted)^2)
  list(
    coef = coef, vcov = vc, free = c(TRUE, TRUE), sigma = sqrt(rss / df), df = df,
    rss = rss, nobs = n, converged = TRUE, x = s, y = v, fitted = fitted,
    residuals = v - fitted
  )
}

mm_params <- function(fit, level, enzyme_conc = NULL) {
  est <- fit$coef
  se <- sqrt(pmax(diag(fit$vcov), 0))
  tq <- t_quantile(level, fit$df)
  out <- tibble::tibble(
    term = c("vmax", "km"), estimate = unname(est), std.error = unname(se),
    statistic = estimate / std.error,
    p.value = 2 * stats::pt(-abs(statistic), df = fit$df),
    conf.low = estimate - tq * std.error, conf.high = estimate + tq * std.error
  )
  if (!is.null(enzyme_conc)) {
    kcat <- est[["vmax"]] / enzyme_conc
    kcat_se <- se[["vmax"]] / enzyme_conc
    eff <- kcat / est[["km"]]
    grad <- c(1 / (enzyme_conc * est[["km"]]), -est[["vmax"]] / (enzyme_conc * est[["km"]]^2))
    eff_se <- sqrt(as.numeric(t(grad) %*% fit$vcov %*% grad))
    extra <- tibble::tibble(
      term = c("kcat", "kcat_km"), estimate = c(kcat, eff), std.error = c(kcat_se, eff_se),
      statistic = estimate / std.error,
      p.value = 2 * stats::pt(-abs(statistic), df = fit$df),
      conf.low = estimate - tq * std.error, conf.high = estimate + tq * std.error
    )
    out <- vctrs::vec_rbind(out, extra)
  }
  out
}

# ---- Monod batch growth ----------------------------------------------------

# Batch Monod model: biomass X grows at rate mu_max S / (Ks + S) and
# substrate S falls at 1 / Y of the biomass growth rate. The mass balance
# X + Y S = C (a constant) lets the equations integrate exactly: mu_max t
# equals (1 + a) log(X / X0) minus a log((C - X) / (C - X0)), where
# a = Ks Y / C. The right-hand side increases with X, so X(t) is found
# by vectorised bisection on (0, C).
monod_solve <- function(times, mu_max, ks, yield, x0, s0) {
  cc <- x0 + yield * s0
  a <- ks * yield / cc
  elapsed <- function(x) ((1 + a) * log(x / x0) - a * log((cc - x) / (cc - x0))) / mu_max
  lo <- rep(x0 * 1e-12, length(times))
  hi <- rep(cc, length(times))
  for (i in seq_len(60L)) {
    mid <- (lo + hi) / 2
    up <- elapsed(mid) < times
    lo[up] <- mid[up]
    hi[!up] <- mid[!up]
  }
  x <- (lo + hi) / 2
  cbind(biomass = x, substrate = pmax((cc - x) / yield, 0))
}

monod_names <- c("mu_max", "ks", "yield", "x0", "s0")

fit_monod <- function(time, biomass, substrate = NULL, s0 = NULL) {
  ord <- order(time)
  time <- time[ord]
  biomass <- biomass[ord]
  substrate <- if (is.null(substrate)) NULL else substrate[ord]
  if (length(unique(time)) < 5L) {
    cli::cli_abort("Monod fits need at least 5 distinct time points.")
  }
  if (any(biomass <= 0, na.rm = TRUE)) {
    cli::cli_abort("Biomass must be positive.")
  }
  sd_x <- stats::sd(biomass, na.rm = TRUE)
  sd_s <- if (is.null(substrate)) 1 else stats::sd(substrate, na.rm = TRUE)
  x0 <- biomass[[1]]
  s_start <- s0 %||% substrate[[1]]
  s_end <- if (is.null(substrate)) 0 else min(substrate, na.rm = TRUE)
  yield0 <- max((max(biomass) - x0) / max(s_start - s_end, 1e-8), 1e-4)
  lx <- log(biomass)
  slopes <- diff(lx) / diff(time)
  mu0 <- max(slopes[is.finite(slopes)], 0.05) * 1.2
  fix_s0 <- !is.null(s0)
  objective <- function(lp) {
    p <- exp(lp)
    s0_use <- if (fix_s0) s0 else p[[5]]
    pred <- monod_solve(time, p[[1]], p[[2]], p[[3]], p[[4]], s0_use)
    r <- (biomass - pred[, 1]) / sd_x
    if (!is.null(substrate)) {
      r <- c(r, (substrate - pred[, 2]) / sd_s)
    }
    sum(r^2, na.rm = TRUE)
  }
  best <- NULL
  for (ks_frac in c(0.01, 0.1, 0.5)) {
    start <- log(c(mu0, ks_frac * s_start, yield0, x0, s_start))
    if (fix_s0) {
      start <- start[1:4]
    }
    res <- stats::optim(start, objective, method = "Nelder-Mead",
                        control = list(maxit = 2000, reltol = 1e-10))
    res <- stats::optim(res$par, objective, method = "BFGS", control = list(maxit = 500))
    if (is.null(best) || res$value < best$value) {
      best <- res
    }
  }
  p <- exp(best$par)
  if (fix_s0) {
    p <- c(p, s0)
  }
  names(p) <- monod_names
  n_obs <- sum(!is.na(biomass)) + if (is.null(substrate)) 0 else sum(!is.na(substrate))
  k <- length(best$par)
  df <- n_obs - k
  sigma2 <- best$value / df
  hess <- tryCatch(stats::optimHess(best$par, objective), error = function(e) NULL)
  se <- rep(NA_real_, 5)
  if (!is.null(hess)) {
    vc_log <- tryCatch(2 * sigma2 * solve(hess), error = function(e) NULL)
    if (!is.null(vc_log)) {
      se[seq_len(k)] <- p[seq_len(k)] * sqrt(pmax(diag(vc_log), 0))
    }
  }
  names(se) <- monod_names
  pred <- monod_solve(time, p[[1]], p[[2]], p[[3]], p[[4]], p[[5]])
  list(
    coef = p, se = se, free = c(rep(TRUE, 4), !fix_s0), sigma = sqrt(sigma2), df = df,
    rss = best$value, nobs = n_obs, converged = best$convergence == 0, time = time,
    biomass = biomass, substrate = substrate, fitted = pred, order = ord
  )
}

#' Monod growth kinetics for batch bioreactors
#'
#' Fits the Monod model of substrate-limited microbial growth to batch
#' fermentation time courses:
#' \deqn{dX/dt = \mu_{max} \frac{S}{K_s + S} X, \quad dS/dt = -\frac{1}{Y}
#' \frac{dX}{dt}}{dX/dt = mu_max S / (Ks + S) X,  dS/dt = -(1 / Y) dX/dt}
#' where X is biomass, S is the limiting substrate and Y is the biomass
#' yield on substrate. The mass balance \eqn{X + Y S = X_0 + Y S_0}{X + Y S =
#' X0 + Y S0} reduces the system to its exact integrated form, which is
#' solved for X at each time by bisection. The parameters are estimated by
#' least squares on the log scale (Nelder-Mead followed by BFGS, from several starts). Biomass
#' and substrate residuals are scaled by their standard deviations so both
#' contribute. Without a substrate column the initial substrate `s0` must
#' be given.
#'
#' @param .data A data frame or tibble.
#' @param time,biomass <[`tidy-select`][dplyr::dplyr_tidy_select]> Time and
#'   biomass columns.
#' @param substrate <[`tidy-select`][dplyr::dplyr_tidy_select]> Optional
#'   substrate concentration column.
#' @param .by <[`tidy-select`][dplyr::dplyr_tidy_select]> Optional columns
#'   identifying separate runs or reactors.
#' @param s0 Initial substrate concentration, required when `substrate` is
#'   not given.
#'
#' @return A `<bio_assay>` with one row per run: the `.by` columns,
#'   `mu_max`, `ks`, `yield`, `x0`, `s0`, `doubling_time` (`log(2) /
#'   mu_max`), `sigma`, `nobs` and `converged`. Standard errors are in
#'   [tidy()][tidy.bio_assay] and fitted trajectories in
#'   [augment()][augment.bio_assay].
#' @references Monod, J. (1949). The growth of bacterial cultures. *Annual
#'   Review of Microbiology*, 3, 371-394.
#'   \doi{10.1146/annurev.mi.03.100149.002103}
#' @export
#' @examples
#' growth <- bio_monod_growth(fermentation_sim, time_h, biomass_g_l,
#'   substrate = glucose_g_l, .by = reactor
#' )
#' growth
bio_monod_growth <- function(.data, time, biomass, substrate = NULL, .by = NULL, s0 = NULL) {
  check_data_frame(.data)
  t_col <- select_one(.data, {{ time }}, "time")
  x_col <- select_one(.data, {{ biomass }}, "biomass")
  s_col <- select_optional(.data, {{ substrate }}, "substrate")
  by <- select_by(.data, {{ .by }})
  check_numeric_col(.data, t_col, "time")
  check_numeric_col(.data, x_col, "biomass")
  if (is.null(s_col)) {
    if (is.null(s0)) {
      cli::cli_abort("Give either a {.arg substrate} column or the initial substrate {.arg s0}.")
    }
    check_positive_number(s0)
  } else {
    check_numeric_col(.data, s_col, "substrate")
  }
  groups <- split_rows(.data, by)
  fits <- lapply(seq_along(groups$rows), function(i) {
    rows <- groups$rows[[i]]
    fit <- fit_monod(
      .data[[t_col]][rows], .data[[x_col]][rows],
      if (is.null(s_col)) NULL else .data[[s_col]][rows], s0
    )
    fit$rows <- rows[fit$order]
    fit$key <- groups$key[i, ]
    fit
  })
  rows <- lapply(fits, function(fit) {
    p <- fit$coef
    tibble::tibble(
      mu_max = p[["mu_max"]], ks = p[["ks"]], yield = p[["yield"]], x0 = p[["x0"]],
      s0 = p[["s0"]], doubling_time = log(2) / p[["mu_max"]], sigma = fit$sigma,
      nobs = fit$nobs, converged = fit$converged
    )
  })
  out <- bind_key(groups$key, rows)
  meta <- list(fits = fits, by = by, x = t_col, y = x_col, substrate = s_col, model = "monod",
               data = tibble::as_tibble(.data))
  new_bio_assay(out, "monod", meta, required = c(by, "mu_max"))
}

monod_params <- function(fit, level = 0.95) {
  tq <- t_quantile(level, fit$df)
  est <- fit$coef
  se <- fit$se
  tibble::tibble(
    term = monod_names, estimate = unname(est), std.error = unname(se),
    statistic = estimate / std.error,
    p.value = 2 * stats::pt(-abs(statistic), df = fit$df),
    conf.low = estimate - tq * std.error, conf.high = estimate + tq * std.error
  )
}
