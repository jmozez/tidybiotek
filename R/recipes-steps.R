# recipes steps for sequence and count features.

# ---- step_bio_kmers ---------------------------------------------------------

#' K-mer features from sequences
#'
#' `step_bio_kmers()` creates a *specification* of a recipe step that turns
#' a column of DNA, RNA or protein sequences into k-mer count or frequency
#' features. The k-mer vocabulary is learnt from the training set when the
#' recipe is prepped, so new data always get the same columns.
#'
#' @param recipe A recipe object. The step is added to its sequence of
#'   operations.
#' @param ... One or more selector functions to choose sequence columns.
#'   See [recipes::selections()].
#' @param role Role for the new feature columns.
#' @param trained A logical showing whether the step has been trained.
#' @param k K-mer length.
#' @param normalise If `TRUE` (the default), use k-mer frequencies rather
#'   than counts, which removes the effect of sequence length.
#' @param vocabulary The learnt k-mers, one character vector per column.
#'   Filled in by [recipes::prep()].
#' @param keep_original_cols Keep the original sequence columns?
#' @param skip Should the step be skipped when the recipe is baked by
#'   [recipes::bake()]?
#' @param id A character string that is unique to this step.
#'
#' @return An updated version of `recipe` with the new step added. When you
#'   [tidy()][recipes::tidy.recipe()] this step, a tibble with columns
#'   `terms`, `k` and `n_kmers` is returned.
#' @family tidybiotek recipe steps
#' @export
#' @examples
#' rec <- recipes::recipe(gc ~ sequence, data = bio_gc_content(ecoli_cds, sequence)) |>
#'   step_bio_kmers(sequence, k = 2) |>
#'   recipes::prep()
#' recipes::bake(rec, new_data = NULL)
#' tidy(rec, number = 1)
step_bio_kmers <- function(recipe, ..., role = "predictor", trained = FALSE, k = 3,
                           normalise = TRUE, vocabulary = NULL, keep_original_cols = FALSE,
                           skip = FALSE, id = recipes::rand_id("bio_kmers")) {
  check_whole_number(k)
  recipes::add_step(
    recipe,
    step_bio_kmers_new(
      terms = rlang::enquos(...), role = role, trained = trained, k = k,
      normalise = normalise, vocabulary = vocabulary, columns = NULL,
      keep_original_cols = keep_original_cols, skip = skip, id = id
    )
  )
}

step_bio_kmers_new <- function(terms, role, trained, k, normalise, vocabulary, columns,
                               keep_original_cols, skip, id) {
  recipes::step(
    subclass = "bio_kmers", terms = terms, role = role, trained = trained, k = k,
    normalise = normalise, vocabulary = vocabulary, columns = columns,
    keep_original_cols = keep_original_cols, skip = skip, id = id
  )
}

#' @exportS3Method recipes::prep
prep.step_bio_kmers <- function(x, training, info = NULL, ...) {
  col_names <- recipes::recipes_eval_select(x$terms, training, info)
  for (col in col_names) {
    check_character_col(training, col, "...")
  }
  vocab <- lapply(col_names, function(col) {
    seqs <- clean_seq(training[[col]])
    kmer_vocabulary(seqs, x$k, guess_alphabet(seqs))
  })
  names(vocab) <- col_names
  step_bio_kmers_new(
    terms = x$terms, role = x$role, trained = TRUE, k = x$k, normalise = x$normalise,
    vocabulary = vocab, columns = col_names, keep_original_cols = x$keep_original_cols,
    skip = x$skip, id = x$id
  )
}

#' @exportS3Method recipes::bake
bake.step_bio_kmers <- function(object, new_data, ...) {
  col_names <- object$columns
  recipes::check_new_data(col_names, object, new_data)
  for (col in col_names) {
    vocab <- object$vocabulary[[col]]
    m <- kmer_matrix(clean_seq(new_data[[col]]), object$k, vocab, object$normalise)
    colnames(m) <- paste0(col, "_", vocab)
    new_data <- vctrs::vec_cbind(new_data, tibble::as_tibble(m))
  }
  if (!isTRUE(object$keep_original_cols)) {
    new_data <- new_data[setdiff(names(new_data), col_names)]
  }
  tibble::as_tibble(new_data)
}

#' @export
print.step_bio_kmers <- function(x, width = max(20, options()$width - 35), ...) {
  title <- paste0("K-mer (k = ", x$k, ") features from ")
  recipes::print_step(x$columns, x$terms, x$trained, title, width)
  invisible(x)
}

#' @rdname step_bio_kmers
#' @param x A `step_bio_kmers` object.
#' @exportS3Method generics::tidy
tidy.step_bio_kmers <- function(x, ...) {
  if (recipes::is_trained(x)) {
    res <- tibble::tibble(
      terms = unname(x$columns), k = x$k,
      n_kmers = unname(lengths(x$vocabulary[x$columns]))
    )
  } else {
    res <- tibble::tibble(terms = recipes::sel2char(x$terms), k = x$k, n_kmers = NA_integer_)
  }
  res$id <- x$id
  res
}

#' @rdname required_pkgs.tidybiotek
#' @exportS3Method generics::required_pkgs
required_pkgs.step_bio_kmers <- function(x, ...) {
  c("tidybiotek")
}

# ---- step_bio_aa_props ------------------------------------------------------

#' Physicochemical protein features
#'
#' `step_bio_aa_props()` creates a *specification* of a recipe step that
#' replaces a column of protein sequences with the physicochemical
#' properties computed by [bio_aa_props()]: length, molecular weight,
#' isoelectric point, charge at pH 7, GRAVY hydropathy, aromaticity and
#' the fractions of acidic, basic and hydrophobic residues.
#'
#' @inheritParams step_bio_kmers
#' @param properties Which properties to keep. Defaults to all nine.
#' @param columns The selected columns. Filled in by [recipes::prep()].
#'
#' @return An updated version of `recipe` with the new step added. When you
#'   [tidy()][recipes::tidy.recipe()] this step, a tibble with columns
#'   `terms` and `property` is returned.
#' @family tidybiotek recipe steps
#' @export
#' @examples
#' rec <- recipes::recipe(class ~ protein, data = ecoli_cds) |>
#'   step_bio_aa_props(protein) |>
#'   recipes::prep()
#' recipes::bake(rec, new_data = NULL)
step_bio_aa_props <- function(recipe, ..., role = "predictor", trained = FALSE,
                              properties = aa_prop_names, columns = NULL,
                              keep_original_cols = FALSE, skip = FALSE,
                              id = recipes::rand_id("bio_aa_props")) {
  bad <- setdiff(properties, aa_prop_names)
  if (length(bad) > 0L) {
    cli::cli_abort("Unknown propert{?y/ies}: {.val {bad}}.")
  }
  recipes::add_step(
    recipe,
    step_bio_aa_props_new(
      terms = rlang::enquos(...), role = role, trained = trained, properties = properties,
      columns = columns, keep_original_cols = keep_original_cols, skip = skip, id = id
    )
  )
}

step_bio_aa_props_new <- function(terms, role, trained, properties, columns,
                                  keep_original_cols, skip, id) {
  recipes::step(
    subclass = "bio_aa_props", terms = terms, role = role, trained = trained,
    properties = properties, columns = columns, keep_original_cols = keep_original_cols,
    skip = skip, id = id
  )
}

#' @exportS3Method recipes::prep
prep.step_bio_aa_props <- function(x, training, info = NULL, ...) {
  col_names <- recipes::recipes_eval_select(x$terms, training, info)
  for (col in col_names) {
    check_character_col(training, col, "...")
  }
  step_bio_aa_props_new(
    terms = x$terms, role = x$role, trained = TRUE, properties = x$properties,
    columns = col_names, keep_original_cols = x$keep_original_cols, skip = x$skip, id = x$id
  )
}

#' @exportS3Method recipes::bake
bake.step_bio_aa_props <- function(object, new_data, ...) {
  col_names <- object$columns
  recipes::check_new_data(col_names, object, new_data)
  for (col in col_names) {
    m <- aa_props_matrix(new_data[[col]])[, object$properties, drop = FALSE]
    colnames(m) <- paste0(col, "_", colnames(m))
    new_data <- vctrs::vec_cbind(new_data, tibble::as_tibble(m))
  }
  if (!isTRUE(object$keep_original_cols)) {
    new_data <- new_data[setdiff(names(new_data), col_names)]
  }
  tibble::as_tibble(new_data)
}

#' @export
print.step_bio_aa_props <- function(x, width = max(20, options()$width - 35), ...) {
  title <- "Protein physicochemical properties from "
  recipes::print_step(x$columns, x$terms, x$trained, title, width)
  invisible(x)
}

#' @rdname step_bio_aa_props
#' @param x A `step_bio_aa_props` object.
#' @exportS3Method generics::tidy
tidy.step_bio_aa_props <- function(x, ...) {
  terms <- if (recipes::is_trained(x)) unname(x$columns) else recipes::sel2char(x$terms)
  res <- expand_terms(terms, x$properties)
  res$id <- x$id
  res
}

expand_terms <- function(terms, properties) {
  term_col <- rep(terms, each = length(properties))
  property_col <- rep(properties, times = length(terms))
  tibble::tibble(terms = term_col, property = property_col)
}

#' @rdname required_pkgs.tidybiotek
#' @exportS3Method generics::required_pkgs
required_pkgs.step_bio_aa_props <- function(x, ...) {
  c("tidybiotek")
}

# ---- step_bio_vst_norm ------------------------------------------------------

vst_transform <- function(q, a0, a1) {
  if (a0 <= 0) {
    # Pure Poisson-like variance: the square-root transform, on a log2-like scale.
    return(2 * sqrt(q * (1 + a1)) / log(2))
  }
  log((1 + a1 + 2 * a0 * q + 2 * sqrt(a0 * q * (1 + a1 + a0 * q))) / (4 * a0)) / log(2)
}

size_factors <- function(m, log_geo_means) {
  usable <- is.finite(log_geo_means)
  if (!any(usable)) {
    cli::cli_abort("Every feature has a zero count in some training sample; cannot normalise.",
                   call = NULL)
  }
  apply(m, 1, function(row) {
    ratios <- log(row[usable]) - log_geo_means[usable]
    exp(stats::median(ratios[is.finite(ratios)]))
  })
}

#' Variance-stabilising normalisation of count features
#'
#' `step_bio_vst_norm()` creates a *specification* of a recipe step that
#' normalises count features (such as RNA-seq read counts, with samples as
#' rows and genes as columns) and stabilises their variance. Three methods
#' are offered:
#'
#' * `"vst"`: median-of-ratios size factors (Anders and Huber, 2010) and the
#'   closed-form variance-stabilising transform for a negative binomial with
#'   dispersion trend \eqn{\alpha(\mu) = a_0 + a_1 / \mu}{alpha(mu) = a0 + a1 / mu},
#'   with the trend fitted to the training data by least squares.
#' * `"log_cpm"`: log2 counts per million with a prior count.
#' * `"anscombe"`: the Anscombe transform \eqn{2 \sqrt{x + 3/8}}{2 sqrt(x + 3/8)}
#'   for Poisson counts.
#'
#' Geometric means and the dispersion trend are learnt from the training
#' set only, so new samples are transformed consistently.
#'
#' @inheritParams step_bio_kmers
#' @param ... Selector functions choosing the count columns.
#' @param method One of `"vst"`, `"log_cpm"` or `"anscombe"`.
#' @param prior_count Prior count added for `"log_cpm"`.
#' @param columns,log_geo_means,dispersion Learnt values filled in by
#'   [recipes::prep()].
#'
#' @return An updated version of `recipe` with the new step added. When you
#'   [tidy()][recipes::tidy.recipe()] this step, a tibble with columns
#'   `terms`, `method` and `log_geo_mean` is returned.
#' @references Anders, S. and Huber, W. (2010). Differential expression
#'   analysis for sequence count data. *Genome Biology*, 11, R106.
#'   \doi{10.1186/gb-2010-11-10-r106}
#'
#'   Anscombe, F. J. (1948). The transformation of Poisson, binomial and
#'   negative-binomial data. *Biometrika*, 35, 246-254.
#'   \doi{10.1093/biomet/35.3-4.246}
#' @family tidybiotek recipe steps
#' @export
#' @examples
#' wide <- bio_counts_wide(expression_sim)
#' rec <- recipes::recipe(condition ~ ., data = wide) |>
#'   recipes::update_role(sample, new_role = "id") |>
#'   step_bio_vst_norm(recipes::all_numeric_predictors()) |>
#'   recipes::prep()
#' recipes::bake(rec, new_data = NULL)[, 1:6]
step_bio_vst_norm <- function(recipe, ..., role = NA, trained = FALSE,
                              method = c("vst", "log_cpm", "anscombe"), prior_count = 0.5,
                              columns = NULL, log_geo_means = NULL, dispersion = NULL,
                              skip = FALSE, id = recipes::rand_id("bio_vst_norm")) {
  method <- rlang::arg_match(method)
  check_positive_number(prior_count)
  recipes::add_step(
    recipe,
    step_bio_vst_norm_new(
      terms = rlang::enquos(...), role = role, trained = trained, method = method,
      prior_count = prior_count, columns = columns, log_geo_means = log_geo_means,
      dispersion = dispersion, skip = skip, id = id
    )
  )
}

step_bio_vst_norm_new <- function(terms, role, trained, method, prior_count, columns,
                                  log_geo_means, dispersion, skip, id) {
  recipes::step(
    subclass = "bio_vst_norm", terms = terms, role = role, trained = trained,
    method = method, prior_count = prior_count, columns = columns,
    log_geo_means = log_geo_means, dispersion = dispersion, skip = skip, id = id
  )
}

#' @exportS3Method recipes::prep
prep.step_bio_vst_norm <- function(x, training, info = NULL, ...) {
  col_names <- recipes::recipes_eval_select(x$terms, training, info)
  recipes::check_type(training[, col_names], types = c("double", "integer"))
  m <- as.matrix(training[, col_names])
  if (any(m < 0, na.rm = TRUE)) {
    cli::cli_abort("{.fn step_bio_vst_norm} needs non-negative counts.", call = NULL)
  }
  log_geo <- colMeans(log(m))
  dispersion <- c(a0 = 0, a1 = 0)
  if (x$method == "vst") {
    sf <- size_factors(m, log_geo)
    norm <- m / sf
    mu <- colMeans(norm)
    v <- apply(norm, 2, stats::var)
    keep <- mu > 0 & is.finite(v)
    disp <- (v[keep] - mu[keep] * mean(1 / sf)) / mu[keep]^2
    fit <- stats::lm(disp ~ I(1 / mu[keep]))
    a <- stats::coef(fit)
    dispersion <- c(a0 = max(a[[1]], 1e-8), a1 = max(a[[2]], 0))
  }
  step_bio_vst_norm_new(
    terms = x$terms, role = x$role, trained = TRUE, method = x$method,
    prior_count = x$prior_count, columns = col_names, log_geo_means = log_geo,
    dispersion = dispersion, skip = x$skip, id = x$id
  )
}

#' @exportS3Method recipes::bake
bake.step_bio_vst_norm <- function(object, new_data, ...) {
  col_names <- object$columns
  recipes::check_new_data(col_names, object, new_data)
  m <- as.matrix(new_data[, col_names])
  out <- switch(object$method,
    vst = {
      sf <- size_factors(m, object$log_geo_means)
      vst_transform(m / sf, object$dispersion[["a0"]], object$dispersion[["a1"]])
    },
    log_cpm = {
      lib <- rowSums(m)
      log2((m + object$prior_count) / (lib + 2 * object$prior_count) * 1e6)
    },
    anscombe = 2 * sqrt(m + 3 / 8)
  )
  out <- matrix(out, nrow = nrow(m), dimnames = dimnames(m))
  for (j in seq_along(col_names)) {
    new_data[[col_names[[j]]]] <- out[, j]
  }
  tibble::as_tibble(new_data)
}

#' @export
print.step_bio_vst_norm <- function(x, width = max(20, options()$width - 35), ...) {
  title <- paste0("Count normalisation (", x$method, ") for ")
  recipes::print_step(x$columns, x$terms, x$trained, title, width)
  invisible(x)
}

#' @rdname step_bio_vst_norm
#' @param x A `step_bio_vst_norm` object.
#' @exportS3Method generics::tidy
tidy.step_bio_vst_norm <- function(x, ...) {
  if (recipes::is_trained(x)) {
    res <- tibble::tibble(
      terms = unname(x$columns), method = x$method, log_geo_mean = unname(x$log_geo_means)
    )
  } else {
    res <- tibble::tibble(
      terms = recipes::sel2char(x$terms), method = x$method, log_geo_mean = NA_real_
    )
  }
  res$id <- x$id
  res
}

#' @rdname required_pkgs.tidybiotek
#' @exportS3Method generics::required_pkgs
required_pkgs.step_bio_vst_norm <- function(x, ...) {
  c("tidybiotek")
}

#' S3 methods for tracking which additional packages are needed for steps
#'
#' Recipe-adjacent packages always list themselves as a required package so
#' that the steps can function properly within parallel processing schemes.
#' @param x A recipe step.
#' @param ... Not used.
#' @return A character vector of package names.
#' @rdname required_pkgs.tidybiotek
#' @keywords internal
#' @name required_pkgs.tidybiotek
NULL

#' Reshape long counts to one row per sample
#'
#' A small helper that turns long count data, such as [expression_sim],
#' into the wide layout that recipes and models expect: one row per sample
#' and one column per gene.
#'
#' @param .data A data frame in long format.
#' @param feature,sample,value <[`tidy-select`][dplyr::dplyr_tidy_select]>
#'   Feature, sample and value columns.
#' @param keep <[`tidy-select`][dplyr::dplyr_tidy_select]> Sample-level
#'   columns to carry across, such as the condition.
#'
#' @return A tibble with one row per sample.
#' @export
#' @examples
#' wide <- bio_counts_wide(expression_sim)
#' dim(wide)
bio_counts_wide <- function(.data, feature = "gene", sample = "sample", value = "count",
                            keep = "condition") {
  check_data_frame(.data)
  f_col <- select_one(.data, {{ feature }}, "feature")
  s_col <- select_one(.data, {{ sample }}, "sample")
  v_col <- select_one(.data, {{ value }}, "value")
  keep_cols <- select_by(.data, {{ keep }})
  samples <- unique(as.character(.data[[s_col]]))
  features <- unique(as.character(.data[[f_col]]))
  m <- matrix(NA_real_, length(samples), length(features),
              dimnames = list(NULL, features))
  m[cbind(match(as.character(.data[[s_col]]), samples),
          match(as.character(.data[[f_col]]), features))] <- .data[[v_col]]
  info <- .data[!duplicated(.data[[s_col]]), c(s_col, keep_cols), drop = FALSE]
  info <- info[match(samples, as.character(info[[s_col]])), , drop = FALSE]
  tibble::as_tibble(vctrs::vec_cbind(tibble::as_tibble(info), tibble::as_tibble(m)))
}
