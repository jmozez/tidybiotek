# yardstick metrics for virtual screening and assay accuracy.

check_binary <- function(truth, estimator, call = rlang::caller_env()) {
  if (!identical(estimator, "binary")) {
    cli::cli_abort("This metric needs a binary {.arg truth} factor with exactly two levels.",
                   call = call)
  }
}

event_flags <- function(truth, event_level) {
  lvls <- levels(truth)
  event <- if (identical(event_level, "second")) lvls[[2]] else lvls[[1]]
  truth == event
}

#' Enrichment factor for virtual screening
#'
#' The enrichment factor (EF) compares the hit rate among the top-scoring
#' `fraction` of a ranked list with the hit rate in the whole list. An EF
#' of 10 at 1 percent means the top 1 percent holds ten times more actives
#' than random selection would give. It is widely used to assess
#' virtual screening, biomarker and protein engineering models that
#' prioritise candidates for testing.
#'
#' @param data A data frame with the `truth` and probability or score
#'   columns.
#' @param truth <[`tidy-select`][dplyr::dplyr_tidy_select]> The true class
#'   (a factor with two levels).
#' @param ... The column with the predicted probability or score for the
#'   event class.
#' @param fraction Fraction of the ranked list to examine, between 0 and 1.
#' @param na_rm Remove missing values before computing?
#' @param event_level `"first"` or `"second"`: which level of `truth` is
#'   the event (active).
#' @param case_weights Not supported; must be `NULL`.
#' @param estimate For the `_vec()` function, a numeric vector of scores.
#'
#' @return For the data frame method, a tibble with columns `.metric`,
#'   `.estimator` and `.estimate`. For `enrichment_factor_vec()`, a single
#'   number.
#' @references Truchon, J.-F. and Bayly, C. I. (2007). Evaluating virtual
#'   screening methods: good and bad metrics for the "early recognition"
#'   problem. *Journal of Chemical Information and Modeling*, 47, 488-508.
#'   \doi{10.1021/ci600426e}
#' @family tidybiotek metrics
#' @export
#' @examples
#' set.seed(2)
#' screen <- tibble::tibble(
#'   active = factor(rep(c("yes", "no"), c(20, 480)), levels = c("yes", "no")),
#'   score = c(rnorm(20, 2), rnorm(480))
#' )
#' enrichment_factor(screen, active, score, fraction = 0.05)
enrichment_factor <- function(data, ...) {
  UseMethod("enrichment_factor")
}
enrichment_factor <- yardstick::new_prob_metric(enrichment_factor, direction = "maximize")

#' @rdname enrichment_factor
#' @export
enrichment_factor.data.frame <- function(data, truth, ..., fraction = 0.01, na_rm = TRUE,
                                         event_level = "first", case_weights = NULL) {
  yardstick::prob_metric_summarizer(
    name = "enrichment_factor", fn = enrichment_factor_vec, data = data,
    truth = !!rlang::enquo(truth), ..., na_rm = na_rm, event_level = event_level,
    case_weights = !!rlang::enquo(case_weights), fn_options = list(fraction = fraction)
  )
}

#' @rdname enrichment_factor
#' @export
enrichment_factor_vec <- function(truth, estimate, fraction = 0.01, na_rm = TRUE,
                                  event_level = "first", case_weights = NULL, ...) {
  if (!is.null(case_weights)) {
    cli::cli_abort("{.fn enrichment_factor} does not support case weights.")
  }
  if (!is.numeric(fraction) || length(fraction) != 1L || fraction <= 0 || fraction > 1) {
    cli::cli_abort("{.arg fraction} must be a single number in (0, 1].")
  }
  estimator <- yardstick::finalize_estimator(truth, metric_class = "enrichment_factor")
  yardstick::check_prob_metric(truth, estimate, case_weights, estimator)
  check_binary(truth, estimator)
  if (na_rm) {
    res <- yardstick::yardstick_remove_missing(truth, estimate, case_weights)
    truth <- res$truth
    estimate <- res$estimate
  } else if (yardstick::yardstick_any_missing(truth, estimate, case_weights)) {
    return(NA_real_)
  }
  active <- event_flags(truth, event_level)
  n <- length(active)
  n_active <- sum(active)
  if (n_active == 0L) {
    cli::cli_warn("There are no actives; the enrichment factor is undefined.")
    return(NA_real_)
  }
  n_top <- max(1L, floor(fraction * n))
  top <- order(estimate, decreasing = TRUE)[seq_len(n_top)]
  (sum(active[top]) / n_top) / (n_active / n)
}

#' Boltzmann-enhanced discrimination of ROC (BEDROC)
#'
#' BEDROC (Truchon and Bayly, 2007) measures early recognition: how well a
#' score ranks actives at the very top of a list. Each active contributes
#' an exponentially decreasing weight of its rank, controlled by `alpha`
#' (the default 20 gives 80 percent of the weight to the top 8 percent of
#' the list). Values range from 0 to 1, where 1 means all actives are
#' ranked first.
#'
#' @inheritParams enrichment_factor
#' @param alpha Early-recognition parameter.
#'
#' @return For the data frame method, a tibble with columns `.metric`,
#'   `.estimator` and `.estimate`. For `bedroc_vec()`, a single number.
#' @inherit enrichment_factor references
#' @family tidybiotek metrics
#' @export
#' @examples
#' set.seed(2)
#' screen <- tibble::tibble(
#'   active = factor(rep(c("yes", "no"), c(20, 480)), levels = c("yes", "no")),
#'   score = c(rnorm(20, 2), rnorm(480))
#' )
#' bedroc(screen, active, score)
bedroc <- function(data, ...) {
  UseMethod("bedroc")
}
bedroc <- yardstick::new_prob_metric(bedroc, direction = "maximize")

#' @rdname bedroc
#' @export
bedroc.data.frame <- function(data, truth, ..., alpha = 20, na_rm = TRUE,
                              event_level = "first", case_weights = NULL) {
  yardstick::prob_metric_summarizer(
    name = "bedroc", fn = bedroc_vec, data = data, truth = !!rlang::enquo(truth), ...,
    na_rm = na_rm, event_level = event_level, case_weights = !!rlang::enquo(case_weights),
    fn_options = list(alpha = alpha)
  )
}

#' @rdname bedroc
#' @export
bedroc_vec <- function(truth, estimate, alpha = 20, na_rm = TRUE, event_level = "first",
                       case_weights = NULL, ...) {
  if (!is.null(case_weights)) {
    cli::cli_abort("{.fn bedroc} does not support case weights.")
  }
  check_positive_number(alpha)
  estimator <- yardstick::finalize_estimator(truth, metric_class = "bedroc")
  yardstick::check_prob_metric(truth, estimate, case_weights, estimator)
  check_binary(truth, estimator)
  if (na_rm) {
    res <- yardstick::yardstick_remove_missing(truth, estimate, case_weights)
    truth <- res$truth
    estimate <- res$estimate
  } else if (yardstick::yardstick_any_missing(truth, estimate, case_weights)) {
    return(NA_real_)
  }
  active <- event_flags(truth, event_level)
  big_n <- length(active)
  n <- sum(active)
  if (n == 0L || n == big_n) {
    cli::cli_warn("BEDROC needs both actives and inactives.")
    return(NA_real_)
  }
  ranks <- rank(-estimate, ties.method = "average")[active]
  ra <- n / big_n
  rie <- sum(exp(-alpha * ranks / big_n)) /
    (n / big_n * (1 - exp(-alpha)) / (exp(alpha / big_n) - 1))
  rie * ra * sinh(alpha / 2) / (cosh(alpha / 2) - cosh(alpha / 2 - alpha * ra)) +
    1 / (1 - exp(alpha * (1 - ra)))
}

#' Assay back-calculation error
#'
#' `assay_cv_error()` is the root-mean-square relative error, in percent,
#' between estimated and nominal concentrations:
#' \deqn{100 \sqrt{\frac{1}{n} \sum_i \left(\frac{\hat{y}_i - y_i}{y_i}\right)^2}}{100
#' sqrt(mean(((estimate - truth) / truth)^2))}
#' It combines bias and imprecision into one CV-like figure and is the
#' usual way to judge how well a standard curve back-calculates its
#' standards or quality controls (acceptance is often within 20 percent).
#'
#' @param data A data frame with the `truth` and `estimate` columns.
#' @param truth,estimate <[`tidy-select`][dplyr::dplyr_tidy_select]> Nominal
#'   and estimated concentrations (vectors for the `_vec()` function).
#' @param na_rm Remove missing values before computing?
#' @param case_weights <[`tidy-select`][dplyr::dplyr_tidy_select]> Optional
#'   case weights.
#' @param ... Not used.
#'
#' @return For the data frame method, a tibble with columns `.metric`,
#'   `.estimator` and `.estimate`. For `assay_cv_error_vec()`, a single
#'   number.
#' @references US Food and Drug Administration (2018). *Bioanalytical
#'   Method Validation: Guidance for Industry*.
#' @family tidybiotek metrics
#' @export
#' @examples
#' standards <- tibble::tibble(nominal = c(1, 2, 4, 8), measured = c(1.1, 1.9, 4.3, 7.6))
#' assay_cv_error(standards, nominal, measured)
assay_cv_error <- function(data, ...) {
  UseMethod("assay_cv_error")
}
assay_cv_error <- yardstick::new_numeric_metric(assay_cv_error, direction = "minimize")

#' @rdname assay_cv_error
#' @export
assay_cv_error.data.frame <- function(data, truth, estimate, na_rm = TRUE,
                                      case_weights = NULL, ...) {
  yardstick::numeric_metric_summarizer(
    name = "assay_cv_error", fn = assay_cv_error_vec, data = data,
    truth = !!rlang::enquo(truth), estimate = !!rlang::enquo(estimate), na_rm = na_rm,
    case_weights = !!rlang::enquo(case_weights)
  )
}

#' @rdname assay_cv_error
#' @export
assay_cv_error_vec <- function(truth, estimate, na_rm = TRUE, case_weights = NULL, ...) {
  yardstick::check_numeric_metric(truth, estimate, case_weights)
  if (na_rm) {
    res <- yardstick::yardstick_remove_missing(truth, estimate, case_weights)
    truth <- res$truth
    estimate <- res$estimate
    case_weights <- res$case_weights
  } else if (yardstick::yardstick_any_missing(truth, estimate, case_weights)) {
    return(NA_real_)
  }
  if (any(truth == 0)) {
    cli::cli_abort("{.arg truth} must not contain zeros for a relative error.")
  }
  rel <- ((estimate - truth) / truth)^2
  w <- if (is.null(case_weights)) rep(1, length(rel)) else as.numeric(case_weights)
  100 * sqrt(sum(w * rel) / sum(w))
}
