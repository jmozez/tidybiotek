# A parsnip model specification for four-parameter logistic regression.

#' Four-parameter logistic regression model
#'
#' `bio_four_pl()` defines a four-parameter logistic (4PL) dose-response
#' model for use with [parsnip::fit()], workflows and tuning. The only
#' engine is `"nls_4pl"`, which uses the Levenberg-Marquardt fitter behind
#' [bio_dose_response()]. The model takes a single numeric predictor (the
#' dose or concentration, on the original scale) and a numeric outcome.
#' Engine arguments `bottom` and `top` fix the asymptotes, for example
#' `set_engine("nls_4pl", bottom = 0)`.
#'
#' @param mode A single character string for the type of model. The only
#'   possible value is `"regression"`.
#' @param engine A single character string for the engine. The only
#'   possible value is `"nls_4pl"`.
#'
#' @return A model specification object with classes `bio_four_pl` and
#'   `model_spec`.
#' @export
#' @examples
#' one <- subset(dose_response_sim, compound == "inhibitor_b")
#' spec <- bio_four_pl()
#' spec
#' fitted <- parsnip::fit(spec, response ~ dose, data = one)
#' tidy(fitted$fit)
#' predict(fitted, new_data = data.frame(dose = c(0.001, 0.01, 0.1)))
bio_four_pl <- function(mode = "regression", engine = "nls_4pl") {
  if (!identical(mode, "regression")) {
    cli::cli_abort("{.arg mode} must be {.val regression}.")
  }
  parsnip::new_model_spec(
    "bio_four_pl",
    args = list(),
    eng_args = NULL,
    mode = mode,
    user_specified_mode = !missing(mode),
    method = NULL,
    engine = engine,
    user_specified_engine = !missing(engine)
  )
}

#' @export
print.bio_four_pl <- function(x, ...) {
  cat("Four-parameter logistic model specification (", x$mode, ")\n\n", sep = "")
  parsnip::model_printer(x, ...)
  if (!is.null(x$method$fit$args)) {
    cat("Model fit template:\n")
    print(parsnip::show_call(x))
  }
  invisible(x)
}

#' Fit function for the `nls_4pl` engine
#'
#' The engine fit function used by [bio_four_pl()]. It is exported so that
#' parsnip can call it, and can also be used directly with a formula of the
#' form `response ~ dose`.
#'
#' @param formula A formula with one numeric predictor.
#' @param data A data frame.
#' @param bottom,top Optional fixed asymptotes.
#' @param ... Not used.
#' @param object A `bio_four_pl_fit` object.
#' @param new_data A data frame with the predictor column.
#'
#' @return `bio_four_pl_fit()` returns an object of class `bio_four_pl_fit`.
#'   The `predict()` method returns a numeric vector.
#' @export
#' @examples
#' one <- subset(dose_response_sim, compound == "inhibitor_b")
#' fit <- bio_four_pl_fit(response ~ dose, one)
#' fit
#' predict(fit, data.frame(dose = 0.01))
bio_four_pl_fit <- function(formula, data, bottom = NULL, top = NULL, ...) {
  processed <- hardhat::mold(formula, data,
    blueprint = hardhat::default_formula_blueprint(intercept = FALSE, indicators = "none")
  )
  x <- processed$predictors
  y <- processed$outcomes
  if (ncol(x) != 1L || !is.numeric(x[[1]])) {
    cli::cli_abort("A 4PL model needs exactly one numeric predictor (the dose).")
  }
  if (ncol(y) != 1L || !is.numeric(y[[1]])) {
    cli::cli_abort("A 4PL model needs one numeric outcome.")
  }
  fit <- fit_4pl(x[[1]], y[[1]], bottom = bottom, top = top)
  fit$params <- four_pl_params(fit)
  structure(
    list(fit = fit, blueprint = processed$blueprint, predictor = names(x), outcome = names(y)),
    class = "bio_four_pl_fit"
  )
}

#' @rdname bio_four_pl_fit
#' @export
predict.bio_four_pl_fit <- function(object, new_data, ...) {
  processed <- hardhat::forge(new_data, object$blueprint)
  four_pl_value(object$fit$coef, processed$predictors[[1]])
}

#' @export
print.bio_four_pl_fit <- function(x, ...) {
  p <- x$fit$params
  cat("Four-parameter logistic fit of `", x$outcome, "` on `", x$predictor, "`\n", sep = "")
  print(p[c("term", "estimate", "std.error")], ...)
  invisible(x)
}

#' @rdname broom_bio
#' @exportS3Method generics::tidy
tidy.bio_four_pl_fit <- function(x, ...) {
  x$fit$params
}

#' @rdname broom_bio
#' @exportS3Method generics::glance
glance.bio_four_pl_fit <- function(x, ...) {
  fit_stats(x$fit)
}

make_bio_four_pl <- function() {
  current <- parsnip::get_model_env()
  if ("bio_four_pl" %in% current$models) {
    return(invisible(FALSE))
  }
  parsnip::set_new_model("bio_four_pl")
  parsnip::set_model_mode("bio_four_pl", "regression")
  parsnip::set_model_engine("bio_four_pl", mode = "regression", eng = "nls_4pl")
  parsnip::set_dependency("bio_four_pl", eng = "nls_4pl", pkg = "tidybiotek",
                          mode = "regression")
  parsnip::set_fit(
    model = "bio_four_pl", eng = "nls_4pl", mode = "regression",
    value = list(
      interface = "formula",
      protect = c("formula", "data"),
      func = c(pkg = "tidybiotek", fun = "bio_four_pl_fit"),
      defaults = list()
    )
  )
  parsnip::set_encoding(
    model = "bio_four_pl", eng = "nls_4pl", mode = "regression",
    options = list(
      predictor_indicators = "none", compute_intercept = FALSE,
      remove_intercept = FALSE, allow_sparse_x = FALSE
    )
  )
  parsnip::set_pred(
    model = "bio_four_pl", eng = "nls_4pl", mode = "regression", type = "numeric",
    value = list(
      pre = NULL, post = NULL, func = c(fun = "predict"),
      args = list(object = rlang::expr(object$fit), new_data = rlang::expr(new_data))
    )
  )
  invisible(TRUE)
}
