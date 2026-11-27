# broom-style tidy(), glance() and augment() methods.

# Fits whose group key is still present in `x` (after any filtering).
active_fits <- function(x) {
  meta <- bio_meta(x)
  fits <- meta$fits
  by <- meta$by
  if (length(by) == 0L || !all(by %in% names(x))) {
    return(fits)
  }
  keys <- tibble::as_tibble(strip_bio(x))[by]
  keep <- vapply(fits, function(f) vctrs::vec_in(f$key, keys), logical(1))
  fits[keep]
}

with_key <- function(fits, f) {
  rows <- lapply(fits, f)
  keys <- vctrs::vec_rbind(!!!lapply(seq_along(fits), function(i) {
    vctrs::vec_rep(fits[[i]]$key, nrow(rows[[i]]))
  }))
  if (ncol(keys) == 0L) {
    return(tibble::as_tibble(vctrs::vec_rbind(!!!rows)))
  }
  tibble::as_tibble(vctrs::vec_cbind(keys, vctrs::vec_rbind(!!!rows)))
}

#' Tidy, glance and augment tidybiotek results
#'
#' Summarise assay, kinetic and omics results in the broom style.
#'
#' * `tidy()` returns one row per parameter (with standard errors,
#'   statistics, p-values and confidence intervals) for dose-response,
#'   ELISA, enzyme kinetic and Monod fits; per-group fold changes for qPCR;
#'   and per-feature or per-gene-set results for omics analyses.
#' * `glance()` returns one row per fitted curve or analysis with model
#'   statistics such as the residual standard deviation, log-likelihood and
#'   AIC, or the empirical Bayes prior for differential expression.
#' * `augment()` returns the original data with fitted values (`.fitted`)
#'   and residuals (`.resid`), or significance flags for omics results.
#'
#' @param x A `<bio_assay>`, `<bio_omics>` or `bio_four_pl_fit` object.
#' @param conf.level Confidence level for intervals.
#' @param alpha Adjusted p-value threshold used to flag significance.
#' @param lfc Minimum absolute log2 fold change for significance.
#' @param ... Not used.
#'
#' @return A tibble.
#' @name broom_bio
#' @examples
#' fit <- bio_dose_response(dose_response_sim, dose, response, .by = compound)
#' tidy(fit)
#' glance(fit)
#' augment(fit)
#'
#' de <- bio_diff_expr(expression_sim, gene, sample, count, condition,
#'   treatment = "stress", reference = "control", transform = "log2_cpm"
#' )
#' tidy(de)
#' glance(de)
NULL

#' @rdname broom_bio
#' @exportS3Method generics::tidy
tidy.bio_assay <- function(x, conf.level = NULL, ...) { # nolint: object_name_linter.
  meta <- bio_meta(x)
  level <- conf.level %||% meta$level %||% 0.95
  switch(meta$assay,
    dose_response = ,
    elisa = with_key(active_fits(x), function(f) four_pl_params(f, level)),
    enzyme_kinetics = with_key(active_fits(x), function(f) {
      mm_params(f, level, meta$enzyme_conc)
    }),
    monod = with_key(active_fits(x), function(f) monod_params(f, level)),
    qpcr = qpcr_summary(x, level),
    tibble::as_tibble(strip_bio(x))
  )
}

#' @rdname broom_bio
#' @exportS3Method generics::glance
glance.bio_assay <- function(x, ...) {
  meta <- bio_meta(x)
  switch(meta$assay,
    dose_response = ,
    elisa = ,
    enzyme_kinetics = ,
    monod = with_key(active_fits(x), fit_stats_any),
    qpcr = tibble::tibble(
      n_samples = length(unique(x[[meta$sample]])),
      n_targets = length(unique(x[[meta$gene]])),
      reference_genes = paste(meta$reference_genes, collapse = ", "),
      control = meta$control
    ),
    tibble::tibble(nrow = nrow(x))
  )
}

fit_stats_any <- function(f) {
  if (is.null(f$free)) {
    f$free <- rep(TRUE, length(f$coef))
  }
  fit_stats(f)
}

#' @rdname broom_bio
#' @exportS3Method generics::augment
augment.bio_assay <- function(x, ...) {
  meta <- bio_meta(x)
  if (meta$assay %in% c("dose_response", "enzyme_kinetics")) {
    data <- meta$data
    fitted <- rep(NA_real_, nrow(data))
    for (f in active_fits(x)) {
      fitted[f$rows] <- f$fitted
    }
    data$.fitted <- fitted
    data$.resid <- data[[meta$y]] - fitted
    return(data[!is.na(fitted), ])
  }
  if (meta$assay == "elisa") {
    out <- tibble::as_tibble(strip_bio(x))
    fitted <- rep(NA_real_, nrow(out))
    for (f in meta$fits) {
      fitted[f$rows] <- f$fitted
    }
    if (nrow(out) == length(fitted)) {
      out$.fitted <- fitted
      out$.resid <- out[[meta$y]] - fitted
    }
    return(out)
  }
  if (meta$assay == "monod") {
    data <- meta$data
    fb <- rep(NA_real_, nrow(data))
    fs <- rep(NA_real_, nrow(data))
    for (f in active_fits(x)) {
      fb[f$rows] <- f$fitted[, "biomass"]
      fs[f$rows] <- f$fitted[, "substrate"]
    }
    data$.fitted <- fb
    data$.resid <- data[[meta$y]] - fb
    data$.fitted_substrate <- fs
    return(data[!is.na(fb), ])
  }
  tibble::as_tibble(strip_bio(x))
}

#' @rdname broom_bio
#' @exportS3Method generics::tidy
tidy.bio_omics <- function(x, ...) {
  meta <- bio_meta(x)
  x <- tibble::as_tibble(strip_bio(x))
  if (meta$analysis == "diff_expr") {
    return(tibble::tibble(
      feature = x$feature, estimate = x$log2fc, std.error = x$std_error, statistic = x$t,
      p.value = x$p_value, p.adjusted = x$p_adj, ave_expr = x$ave_expr
    ))
  }
  if (meta$analysis == "enrichment") {
    return(tibble::tibble(
      gene_set = x$gene_set, estimate = x$fold_enrichment, overlap = x$overlap,
      set_size = x$set_size, p.value = x$p_value, p.adjusted = x$p_adj
    ))
  }
  x
}

#' @rdname broom_bio
#' @exportS3Method generics::glance
glance.bio_omics <- function(x, alpha = 0.05, ...) {
  meta <- bio_meta(x)
  if (meta$analysis == "diff_expr") {
    return(tibble::tibble(
      n_features = nrow(x), n_significant = sum(x$p_adj < alpha, na.rm = TRUE),
      df_prior = meta$df_prior, s0_sq = meta$s0_sq, df_residual = meta$df_residual,
      treatment = meta$treatment, reference = meta$reference
    ))
  }
  if (meta$analysis == "enrichment") {
    return(tibble::tibble(
      n_sets = nrow(x), n_significant = sum(x$p_adj < alpha, na.rm = TRUE),
      n_selected = meta$n_selected, universe_size = meta$universe_size
    ))
  }
  tibble::tibble(nrow = nrow(x))
}

#' @rdname broom_bio
#' @exportS3Method generics::augment
augment.bio_omics <- function(x, alpha = 0.05, lfc = 0, ...) {
  meta <- bio_meta(x)
  out <- tibble::as_tibble(strip_bio(x))
  if (meta$analysis == "diff_expr") {
    sig <- out$p_adj < alpha & abs(out$log2fc) >= lfc
    out$.significant <- sig
    out$.direction <- ifelse(sig, ifelse(out$log2fc > 0, "up", "down"), "none")
  } else if (meta$analysis == "enrichment") {
    out$.significant <- out$p_adj < alpha
  }
  out
}
