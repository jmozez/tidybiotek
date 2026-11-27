# Assay quantification: qPCR, four-parameter logistic dose-response curves
# and ELISA standard curves.

# ---- Levenberg-Marquardt least squares ------------------------------------

# Minimise sum(w * (y - f(par))^2). `model(par)` returns fitted values and
# `jacobian(par)` the matrix of derivatives of the fitted values.
levenberg_marquardt <- function(par, y, model, jacobian, weights = NULL,
                                free = rep(TRUE, length(par)), maxit = 500L,
                                tol = 1e-10) {
  sw <- sqrt(weights %||% rep(1, length(y)))
  resid <- function(p) sw * (y - model(p))
  r <- resid(par)
  rss <- sum(r^2)
  lambda <- 1e-3
  converged <- FALSE
  iter <- 0L
  while (iter < maxit) {
    iter <- iter + 1L
    jac <- sw * jacobian(par)[, free, drop = FALSE]
    jac[!is.finite(jac)] <- 0
    grad <- crossprod(jac, r)
    hess <- crossprod(jac)
    improved <- FALSE
    while (lambda < 1e12) {
      damped <- hess + lambda * diag(pmax(diag(hess), 1e-12), nrow = nrow(hess))
      step <- tryCatch(solve(damped, grad), error = function(e) NULL)
      if (!is.null(step)) {
        new <- par
        new[free] <- par[free] + as.numeric(step)
        r_new <- resid(new)
        rss_new <- sum(r_new^2)
        if (is.finite(rss_new) && rss_new <= rss) {
          improved <- TRUE
          break
        }
      }
      lambda <- lambda * 10
    }
    if (!improved) {
      converged <- TRUE
      break
    }
    lambda <- max(lambda / 10, 1e-12)
    delta <- rss - rss_new
    par <- new
    r <- r_new
    rss <- rss_new
    if (delta <= tol * (rss + tol)) {
      converged <- TRUE
      break
    }
  }
  jac <- sw * jacobian(par)[, free, drop = FALSE]
  jac[!is.finite(jac)] <- 0
  list(par = par, rss = rss, converged = converged, iterations = iter, jacobian = jac)
}

# ---- four-parameter logistic ----------------------------------------------

# GraphPad parameterisation: the response is bottom plus (top - bottom)
# divided by 1 + (ec50 / x) to the power hill. The parameter vector holds
# bottom, top, log(ec50) and hill. A positive hill rises with dose and a
# negative hill falls (inhibition, IC50).
four_pl_value <- function(par, x) {
  z <- par[[4]] * (par[[3]] - log(x))
  par[[1]] + (par[[2]] - par[[1]]) * stats::plogis(-z)
}

four_pl_jacobian <- function(par, x) {
  lx <- log(x)
  z <- par[[4]] * (par[[3]] - lx)
  span <- par[[2]] - par[[1]]
  dz <- stats::dlogis(z)
  d_le <- -span * par[[4]] * dz
  d_h <- -span * (par[[3]] - lx) * dz
  d_h[!is.finite(lx)] <- 0
  cbind(stats::plogis(z), stats::plogis(-z), d_le, d_h)
}

four_pl_names <- c("bottom", "top", "log_ec50", "hill")

fit_4pl <- function(x, y, weights = NULL, bottom = NULL, top = NULL) {
  ok <- is.finite(y) & !is.na(x)
  x <- x[ok]
  y <- y[ok]
  weights <- if (is.null(weights)) NULL else weights[ok]
  if (any(x < 0)) {
    cli::cli_abort("Doses and concentrations must be zero or positive.")
  }
  pos <- x > 0
  if (length(unique(x[pos])) < 3L) {
    cli::cli_abort("A four-parameter logistic fit needs at least 3 distinct positive doses.")
  }
  free <- c(is.null(bottom), is.null(top), TRUE, TRUE)
  if (length(y) <= sum(free)) {
    cli::cli_abort("Too few observations ({length(y)}) for {sum(free)} free parameters.")
  }
  lx <- log(x[pos])
  rising <- stats::cor(lx, y[pos]) >= 0
  b0 <- bottom %||% min(y)
  t0 <- top %||% max(y)
  mid <- (b0 + t0) / 2
  le0 <- stats::median(lx[order(abs(y[pos] - mid))][seq_len(min(3L, length(lx)))])
  slopes <- c(0.5, 1, 2, 4) * if (rising) 1 else -1
  model <- function(p) four_pl_value(p, x)
  jacobian <- function(p) four_pl_jacobian(p, x)
  best <- NULL
  for (h in slopes) {
    res <- levenberg_marquardt(c(b0, t0, le0, h), y, model, jacobian, weights, free)
    if (is.null(best) || res$rss < best$rss) {
      best <- res
    }
  }
  par <- stats::setNames(best$par, four_pl_names)
  n <- length(y)
  p <- sum(free)
  df <- n - p
  sigma2 <- best$rss / df
  vc <- matrix(0, 4, 4, dimnames = list(four_pl_names, four_pl_names))
  inv <- tryCatch(solve(crossprod(best$jacobian)), error = function(e) {
    matrix(NA_real_, p, p)
  })
  vc[free, free] <- sigma2 * inv
  if (par[["top"]] < par[["bottom"]] && is.null(bottom) && is.null(top)) {
    par <- par[c(2, 1, 3, 4)]
    par[["hill"]] <- -par[["hill"]]
    names(par) <- four_pl_names
    perm <- c(2, 1, 3, 4)
    vc <- vc[perm, perm]
    sgn <- c(1, 1, 1, -1)
    vc <- vc * outer(sgn, sgn)
    dimnames(vc) <- list(four_pl_names, four_pl_names)
  }
  fitted <- four_pl_value(par, x)
  list(
    coef = par, vcov = vc, free = free, sigma = sqrt(sigma2), df = df, rss = best$rss,
    nobs = n, converged = best$converged, iterations = best$iterations, x = x, y = y,
    fitted = fitted, residuals = y - fitted
  )
}

four_pl_params <- function(fit, level = 0.95) {
  est <- fit$coef
  se <- sqrt(pmax(diag(fit$vcov), 0))
  se[!fit$free] <- NA_real_
  tq <- t_quantile(level, fit$df)
  stat <- est / se
  out <- tibble::tibble(
    term = c("bottom", "top", "ec50", "hill"),
    estimate = c(est[["bottom"]], est[["top"]], exp(est[["log_ec50"]]), est[["hill"]]),
    std.error = c(se[["bottom"]], se[["top"]], exp(est[["log_ec50"]]) * se[["log_ec50"]],
                  se[["hill"]]),
    statistic = c(stat[["bottom"]], stat[["top"]], NA_real_, stat[["hill"]]),
    p.value = 2 * stats::pt(-abs(statistic), df = fit$df),
    conf.low = c(est[1:2] - tq * se[1:2], exp(est[[3]] - tq * se[[3]]), est[[4]] - tq * se[[4]]),
    conf.high = c(est[1:2] + tq * se[1:2], exp(est[[3]] + tq * se[[3]]), est[[4]] + tq * se[[4]])
  )
  out
}

fit_stats <- function(fit) {
  n <- fit$nobs
  k <- sum(fit$free) + 1
  loglik <- -n / 2 * (log(2 * pi) + log(fit$rss / n) + 1)
  tibble::tibble(
    nobs = n, df.residual = fit$df, sigma = fit$sigma, rss = fit$rss,
    logLik = loglik, AIC = -2 * loglik + 2 * k, BIC = -2 * loglik + log(n) * k,
    converged = fit$converged
  )
}

#' Four-parameter logistic dose-response curves (IC50 and EC50)
#'
#' Fits the four-parameter logistic (4PL) model
#' \deqn{y = \mathrm{bottom} + \frac{\mathrm{top} - \mathrm{bottom}}{1 +
#' (\mathrm{EC50} / x)^{\mathrm{hill}}}}{y = bottom + (top - bottom) /
#' (1 + (EC50 / x)^hill)}
#' to each group by Levenberg-Marquardt least squares, written in base R.
#' `bottom` is always the lower asymptote and `top` the upper one. A
#' negative Hill slope means the response falls with dose, so the midpoint
#' is reported as an IC50; a positive slope gives an EC50. Zero doses are
#' allowed and are treated as true zeros. Confidence intervals for the
#' EC50 or IC50 are computed on the log scale.
#'
#' @param .data A data frame or tibble.
#' @param dose,response <[`tidy-select`][dplyr::dplyr_tidy_select]> Columns
#'   holding dose (or concentration, on the original scale) and response.
#' @param .by <[`tidy-select`][dplyr::dplyr_tidy_select]> Optional columns
#'   that define separate curves, such as compound or cell line.
#' @param weights <[`tidy-select`][dplyr::dplyr_tidy_select]> Optional
#'   column of case weights.
#' @param bottom,top Optional fixed values for the asymptotes, for example
#'   `bottom = 0` for data normalised to percent of control.
#' @param level Confidence level.
#'
#' @return A `<bio_assay>` with one row per curve: the `.by` columns,
#'   `bottom`, `top`, `ec50`, `ec50_low`, `ec50_high`, `hill`, `type`
#'   (`"IC50"` or `"EC50"`), `sigma`, `nobs` and `converged`. Use
#'   [tidy()][tidy.bio_assay] for full parameter tables and
#'   [autoplot()][autoplot.bio_assay] to draw the curves.
#' @references Ritz, C., Baty, F., Streibig, J. C. and Gerhard, D. (2015).
#'   Dose-response analysis using R. *PLOS ONE*, 10, e0146021.
#'   \doi{10.1371/journal.pone.0146021}
#' @export
#' @examples
#' fit <- bio_dose_response(dose_response_sim, dose, response, .by = compound)
#' fit
#' tidy(fit)
bio_dose_response <- function(.data, dose, response, .by = NULL, weights = NULL,
                              bottom = NULL, top = NULL, level = 0.95) {
  check_data_frame(.data)
  dose_col <- select_one(.data, {{ dose }}, "dose")
  resp_col <- select_one(.data, {{ response }}, "response")
  w_col <- select_optional(.data, {{ weights }}, "weights")
  by <- select_by(.data, {{ .by }})
  check_numeric_col(.data, dose_col, "dose")
  check_numeric_col(.data, resp_col, "response")
  groups <- split_rows(.data, by)
  fits <- lapply(seq_along(groups$rows), function(i) {
    rows <- groups$rows[[i]]
    w <- if (is.null(w_col)) NULL else .data[[w_col]][rows]
    fit <- fit_4pl(.data[[dose_col]][rows], .data[[resp_col]][rows], w, bottom, top)
    fit$rows <- rows
    fit$key <- groups$key[i, ]
    fit$params <- four_pl_params(fit, level)
    fit
  })
  rows <- lapply(fits, function(fit) {
    p <- fit$params
    if (!fit$converged) {
      cli::cli_warn("The 4PL fit did not converge for one curve; check its estimates.")
    }
    tibble::tibble(
      bottom = p$estimate[[1]], top = p$estimate[[2]], ec50 = p$estimate[[3]],
      ec50_low = p$conf.low[[3]], ec50_high = p$conf.high[[3]], hill = p$estimate[[4]],
      type = if (p$estimate[[4]] < 0) "IC50" else "EC50",
      sigma = fit$sigma, nobs = fit$nobs, converged = fit$converged
    )
  })
  out <- bind_key(groups$key, rows)
  meta <- list(fits = fits, by = by, x = dose_col, y = resp_col, level = level,
               model = "four_pl", data = tibble::as_tibble(.data))
  new_bio_assay(out, "dose_response", meta, required = c(by, "ec50"))
}

#' ELISA quantification from a standard curve
#'
#' Fits a four-parameter logistic standard curve to the wells with known
#' concentrations and back-calculates the concentration of every well
#' from its response. Wells whose response lies outside the fitted
#' asymptotes cannot be quantified and get `NA` with `in_range = FALSE`.
#' Unknowns are flagged as extrapolated when their estimate falls outside
#' the range of the standards.
#'
#' @param .data A data frame or tibble with one row per well.
#' @param response <[`tidy-select`][dplyr::dplyr_tidy_select]> Measured
#'   signal (for example optical density).
#' @param concentration <[`tidy-select`][dplyr::dplyr_tidy_select]> Known
#'   concentration for standards and `NA` for unknown samples.
#' @param .by <[`tidy-select`][dplyr::dplyr_tidy_select]> Optional columns
#'   defining separate plates or runs, each with its own standard curve.
#' @param dilution <[`tidy-select`][dplyr::dplyr_tidy_select]> Optional
#'   column of dilution factors; estimates are multiplied by it.
#' @param level Confidence level for the curve parameters.
#'
#' @return A `<bio_assay>` holding `.data` plus `estimated_conc` (after
#'   dilution correction), `in_range`, `extrapolated` and `is_standard`. The
#'   standard curves are kept for [tidy()][tidy.bio_assay],
#'   [glance()][glance.bio_assay] and [autoplot()][autoplot.bio_assay].
#' @references Findlay, J. W. A. and Dillard, R. F. (2007). Appropriate
#'   calibration curve fitting in ligand binding assays. *The AAPS Journal*,
#'   9, E260-E267. \doi{10.1208/aapsj0902029}
#' @export
#' @examples
#' # Run 1 of the DNase ELISA as the standard curve; unknowns from run 2
#' plate <- rbind(
#'   transform(subset(dnase_elisa, run == 1), known = conc),
#'   transform(subset(dnase_elisa, run == 2), known = NA)
#' )
#' quant <- bio_elisa_quant(plate, density, known)
#' quant[!quant$is_standard, c("conc", "density", "estimated_conc")]
bio_elisa_quant <- function(.data, response, concentration, .by = NULL, dilution = NULL,
                            level = 0.95) {
  check_data_frame(.data)
  resp_col <- select_one(.data, {{ response }}, "response")
  conc_col <- select_one(.data, {{ concentration }}, "concentration")
  dil_col <- select_optional(.data, {{ dilution }}, "dilution")
  by <- select_by(.data, {{ .by }})
  check_numeric_col(.data, resp_col, "response")
  check_numeric_col(.data, conc_col, "concentration")
  groups <- split_rows(.data, by)
  n <- nrow(.data)
  est <- rep(NA_real_, n)
  in_range <- rep(FALSE, n)
  extrap <- rep(FALSE, n)
  fits <- vector("list", length(groups$rows))
  for (i in seq_along(groups$rows)) {
    rows <- groups$rows[[i]]
    std <- rows[!is.na(.data[[conc_col]][rows])]
    if (length(std) < 5L) {
      cli::cli_abort("Each standard curve needs at least 5 wells with known concentrations.")
    }
    fit <- fit_4pl(.data[[conc_col]][std], .data[[resp_col]][std])
    fit$rows <- std
    fit$key <- groups$key[i, ]
    fit$params <- four_pl_params(fit, level)
    fits[[i]] <- fit
    back <- four_pl_inverse(fit$coef, .data[[resp_col]][rows])
    est[rows] <- back
    in_range[rows] <- !is.na(back)
    std_range <- range(.data[[conc_col]][std][.data[[conc_col]][std] > 0])
    extrap[rows] <- !is.na(back) & (back < std_range[[1]] | back > std_range[[2]])
  }
  if (!is.null(dil_col)) {
    est <- est * .data[[dil_col]]
  }
  out <- tibble::as_tibble(.data)
  out$estimated_conc <- est
  out$in_range <- in_range
  out$extrapolated <- extrap
  out$is_standard <- !is.na(.data[[conc_col]])
  meta <- list(fits = fits, by = by, x = conc_col, y = resp_col, level = level,
               model = "four_pl", data = tibble::as_tibble(.data))
  new_bio_assay(out, "elisa", meta, required = c(resp_col, "estimated_conc"))
}

four_pl_inverse <- function(par, y) {
  p <- (y - par[["bottom"]]) / (par[["top"]] - par[["bottom"]])
  out <- rep(NA_real_, length(y))
  ok <- !is.na(p) & p > 0 & p < 1
  u <- 1 / p[ok] - 1
  out[ok] <- exp(par[["log_ec50"]] - log(u) / par[["hill"]])
  out
}

#' Relative quantification of qPCR data (delta-delta Ct and Pfaffl)
#'
#' Calculates relative expression from quantification cycle (Ct or Cq)
#' values. Technical replicates are averaged within each sample and gene.
#' Each target gene is normalised to the reference gene(s) and then to the
#' mean of the control group. With the default amplification efficiency of
#' 2 (100 percent) this is the 2^-ddCt method of Livak and Schmittgen
#' (2001); gene-specific efficiencies give the Pfaffl (2001) ratio. Several
#' reference genes are combined by the geometric mean of their relative
#' quantities (Vandesompele et al., 2002).
#'
#' @param .data A data frame with one row per well.
#' @param ct,gene,sample,group <[`tidy-select`][dplyr::dplyr_tidy_select]>
#'   Columns holding the Ct value, gene name, biological sample and
#'   experimental group.
#' @param reference_genes Character vector of reference (housekeeping) genes.
#' @param control Level of `group` used as the calibrator.
#' @param efficiency Amplification factor per cycle, between 1 and 2: a
#'   single number, or a named vector with one value per gene.
#'
#' @return A `<bio_assay>` with one row per sample and target gene: the
#'   `sample`, `group` and `gene` columns, `ct_target`, `ct_reference`,
#'   `delta_ct`, `delta_delta_ct`, `fold_change` and `log2_fold_change`.
#'   [tidy()][tidy.bio_assay] summarises the fold change per group with
#'   confidence intervals and Welch t-tests against the control.
#' @references Livak, K. J. and Schmittgen, T. D. (2001). Analysis of
#'   relative gene expression data using real-time quantitative PCR and the
#'   2^-ddCt method. *Methods*, 25, 402-408. \doi{10.1006/meth.2001.1262}
#'
#'   Pfaffl, M. W. (2001). A new mathematical model for relative
#'   quantification in real-time RT-PCR. *Nucleic Acids Research*, 29, e45.
#'   \doi{10.1093/nar/29.9.e45}
#'
#'   Vandesompele, J. et al. (2002). Accurate normalization of real-time
#'   quantitative RT-PCR data by geometric averaging of multiple internal
#'   control genes. *Genome Biology*, 3, research0034.
#'   \doi{10.1186/gb-2002-3-7-research0034}
#' @export
#' @examples
#' rq <- bio_qpcr_ddct(qpcr_sim, ct, gene, sample, treatment,
#'   reference_genes = c("gapA", "rpoD"), control = "control"
#' )
#' rq
#' tidy(rq)
bio_qpcr_ddct <- function(.data, ct, gene, sample, group, reference_genes, control,
                          efficiency = 2) {
  check_data_frame(.data)
  ct_col <- select_one(.data, {{ ct }}, "ct")
  gene_col <- select_one(.data, {{ gene }}, "gene")
  sample_col <- select_one(.data, {{ sample }}, "sample")
  group_col <- select_one(.data, {{ group }}, "group")
  check_numeric_col(.data, ct_col, "ct")
  genes <- unique(as.character(.data[[gene_col]]))
  if (!is.character(reference_genes) || !all(reference_genes %in% genes)) {
    missing <- setdiff(reference_genes, genes)
    cli::cli_abort(c(
      "Every reference gene must be in the gene column {.field {gene_col}}.",
      "x" = "Not found: {.val {missing}}."
    ))
  }
  targets <- setdiff(genes, reference_genes)
  if (length(targets) == 0L) {
    cli::cli_abort("There are no target genes besides the reference genes.")
  }
  if (!rlang::is_string(control) || !control %in% as.character(.data[[group_col]])) {
    cli::cli_abort("{.arg control} must be a level of {.field {group_col}}.")
  }
  eff <- resolve_efficiency(efficiency, genes)

  # Mean Ct per sample and gene (technical replicates).
  key <- paste(.data[[sample_col]], .data[[gene_col]], sep = "\r")
  mean_ct <- tapply(.data[[ct_col]], key, mean, na.rm = TRUE)
  mean_ct <- stats::setNames(as.vector(mean_ct), names(mean_ct))
  sample_vec <- .data[[sample_col]]
  group_vec <- as.character(.data[[group_col]])
  samples <- unique(tibble::tibble(sample = sample_vec, group = group_vec))
  if (anyDuplicated(samples$sample)) {
    cli::cli_abort("Each sample must belong to a single group.")
  }
  get_ct <- function(s, g) {
    v <- mean_ct[paste(s, g, sep = "\r")]
    unname(ifelse(is.na(v), NA_real_, v))
  }
  is_ctrl <- samples$group == control
  ref_ct <- vapply(reference_genes, function(g) get_ct(samples$sample, g),
                   numeric(nrow(samples)))
  ref_ct <- matrix(ref_ct, nrow = nrow(samples))
  ref_ctrl_mean <- colMeans(ref_ct[is_ctrl, , drop = FALSE], na.rm = TRUE)
  ref_rq <- vapply(seq_along(reference_genes), function(j) {
    eff[[reference_genes[[j]]]]^(ref_ctrl_mean[[j]] - ref_ct[, j])
  }, numeric(nrow(samples)))
  ref_rq <- matrix(ref_rq, nrow = nrow(samples))
  ref_norm <- exp(rowMeans(log(ref_rq)))
  ref_mean_ct <- rowMeans(ref_ct)

  rows <- lapply(targets, function(g) {
    tct <- get_ct(samples$sample, g)
    dct <- tct - ref_mean_ct
    ddct <- dct - mean(dct[is_ctrl], na.rm = TRUE)
    ratio <- eff[[g]]^(mean(tct[is_ctrl], na.rm = TRUE) - tct) / ref_norm
    tibble::tibble(
      !!sample_col := samples$sample, !!group_col := samples$group, !!gene_col := g,
      ct_target = tct, ct_reference = ref_mean_ct, delta_ct = dct,
      delta_delta_ct = ddct, fold_change = ratio, log2_fold_change = log2(ratio)
    )
  })
  out <- vctrs::vec_rbind(!!!rows)
  meta <- list(sample = sample_col, group = group_col, gene = gene_col, control = control,
               reference_genes = reference_genes, efficiency = eff, data = tibble::as_tibble(.data))
  new_bio_assay(out, "qpcr", meta, required = c(group_col, gene_col, "fold_change"))
}

resolve_efficiency <- function(efficiency, genes, call = rlang::caller_env()) {
  if (!is.numeric(efficiency) || anyNA(efficiency) || any(efficiency <= 1 | efficiency > 2.2)) {
    cli::cli_abort(
      "{.arg efficiency} must hold amplification factors between 1 and 2 (for example 1.95).",
      call = call
    )
  }
  if (length(efficiency) == 1L && is.null(names(efficiency))) {
    return(stats::setNames(rep(efficiency, length(genes)), genes))
  }
  if (is.null(names(efficiency)) || !all(genes %in% names(efficiency))) {
    cli::cli_abort("A vector {.arg efficiency} must be named with every gene.", call = call)
  }
  efficiency[genes]
}

qpcr_summary <- function(x, level = 0.95) {
  meta <- bio_meta(x)
  g <- meta$group
  gene <- meta$gene
  x <- strip_bio(x)
  combos <- unique(x[c(gene, g)])
  rows <- lapply(seq_len(nrow(combos)), function(i) {
    sel <- x[[gene]] == combos[[gene]][[i]] & x[[g]] == combos[[g]][[i]]
    ctrl <- x[[gene]] == combos[[gene]][[i]] & x[[g]] == meta$control
    l2 <- x$log2_fold_change[sel]
    l2 <- l2[is.finite(l2)]
    n <- length(l2)
    m <- mean(l2)
    se <- if (n > 1) stats::sd(l2) / sqrt(n) else NA_real_
    tq <- t_quantile(level, n - 1)
    pval <- NA_real_
    l2c <- x$log2_fold_change[ctrl]
    if (combos[[g]][[i]] != meta$control && n > 1 && sum(is.finite(l2c)) > 1) {
      pval <- stats::t.test(l2, l2c[is.finite(l2c)])$p.value
    }
    tibble::tibble(
      n = n, log2_fold_change = m, std.error = se, fold_change = 2^m,
      conf.low = 2^(m - tq * se), conf.high = 2^(m + tq * se), p.value = pval
    )
  })
  tibble::as_tibble(vctrs::vec_cbind(tibble::as_tibble(combos), vctrs::vec_rbind(!!!rows)))
}
