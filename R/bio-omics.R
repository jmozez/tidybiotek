# Omics: moderated t-test differential expression and over-representation
# analysis.

trigamma_inverse <- function(x) {
  out <- x
  big <- x > 1e7
  small <- x < 1e-6
  out[big] <- 1 / sqrt(x[big])
  out[small] <- 1 / x[small]
  mid <- !big & !small
  if (any(mid)) {
    y <- 0.5 + 1 / x[mid]
    for (i in seq_len(50L)) {
      tri <- trigamma(y)
      dif <- tri * (1 - tri / x[mid]) / psigamma(y, deriv = 2)
      y <- y + dif
      if (max(-dif / y) < 1e-8) {
        break
      }
    }
    out[mid] <- y
  }
  out
}

# Fit a scaled F distribution to sample variances (Smyth, 2004).
fit_f_dist <- function(s2, df) {
  s2 <- pmax(s2, 0)
  m <- stats::median(s2)
  if (m == 0) {
    m <- 1
  }
  s2 <- pmax(s2, 1e-5 * m)
  z <- log(s2)
  e <- z - digamma(df / 2) + log(df / 2)
  n <- length(e)
  emean <- mean(e)
  evar <- sum((e - emean)^2) / (n - 1) - mean(trigamma(df / 2))
  if (n > 1 && evar > 0) {
    df0 <- 2 * trigamma_inverse(evar)
    s0_sq <- exp(emean + digamma(df0 / 2) - log(df0 / 2))
  } else {
    df0 <- Inf
    s0_sq <- exp(emean)
  }
  list(df0 = df0, s0_sq = s0_sq)
}

#' Differential expression with moderated t-statistics
#'
#' Compares two groups feature by feature (genes, transcripts, proteins or
#' metabolites) using the empirical Bayes moderated t-test of Smyth (2004),
#' the method behind limma, written in base R. A one-way linear model is
#' fitted to every feature, the residual variances are shrunk towards a
#' common prior estimated from all features, and moderated t-statistics
#' with extra degrees of freedom are formed for the chosen contrast.
#' P-values are adjusted with the Benjamini-Hochberg method.
#'
#' The input is tidy (long) data with one row per feature and sample.
#' Counts can be transformed to log2 counts per million first.
#'
#' @param .data A data frame in long format.
#' @param feature,sample,value,group <[`tidy-select`][dplyr::dplyr_tidy_select]>
#'   Columns identifying the feature, the sample, the measured value and
#'   the experimental group of each sample.
#' @param treatment,reference Levels of `group` to compare (log fold change
#'   is `treatment - reference`). By default the first level is the
#'   reference and the second the treatment.
#' @param transform `"none"` (values are already on a log scale),
#'   `"log2"` (adds 1 then takes log2) or `"log2_cpm"` (counts to log2
#'   counts per million with a prior count of 0.5).
#' @param robust_prior If `FALSE`, use ordinary t-tests (no variance
#'   moderation) instead.
#'
#' @return A `<bio_omics>` sorted by p-value, with columns `feature`,
#'   `log2fc`, `ave_expr`, `t`, `p_value`, `p_adj` and `std_error`. The
#'   prior degrees of freedom and variance are kept as metadata; see
#'   [glance()][glance.bio_omics].
#' @references Smyth, G. K. (2004). Linear models and empirical Bayes methods
#'   for assessing differential expression in microarray experiments.
#'   *Statistical Applications in Genetics and Molecular Biology*, 3, 3.
#'   \doi{10.2202/1544-6115.1027}
#'
#'   Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery
#'   rate: a practical and powerful approach to multiple testing. *Journal
#'   of the Royal Statistical Society B*, 57, 289-300.
#'   \doi{10.1111/j.2517-6161.1995.tb02031.x}
#' @export
#' @examples
#' de <- bio_diff_expr(expression_sim, gene, sample, count, condition,
#'   treatment = "stress", reference = "control", transform = "log2_cpm"
#' )
#' de
#' glance(de)
bio_diff_expr <- function(.data, feature, sample, value, group, treatment = NULL,
                          reference = NULL, transform = c("none", "log2", "log2_cpm"),
                          robust_prior = TRUE) {
  check_data_frame(.data)
  f_col <- select_one(.data, {{ feature }}, "feature")
  s_col <- select_one(.data, {{ sample }}, "sample")
  v_col <- select_one(.data, {{ value }}, "value")
  g_col <- select_one(.data, {{ group }}, "group")
  transform <- rlang::arg_match(transform)
  check_numeric_col(.data, v_col, "value")
  if (anyDuplicated(.data[c(f_col, s_col)])) {
    cli::cli_abort("Each feature must have at most one value per sample.")
  }
  sample_vec <- .data[[s_col]]
  group_vec <- .data[[g_col]]
  samples <- unique(tibble::tibble(sample = sample_vec, group = group_vec))
  if (anyDuplicated(samples$sample)) {
    cli::cli_abort("Each sample must belong to a single group.")
  }
  levels <- if (is.factor(samples$group)) levels(droplevels(samples$group)) else
    unique(as.character(samples$group))
  reference <- reference %||% levels[[1]]
  treatment <- treatment %||% setdiff(levels, reference)[[1]]
  if (!all(c(reference, treatment) %in% levels) || identical(reference, treatment)) {
    cli::cli_abort("{.arg treatment} and {.arg reference} must be two different groups.")
  }
  features <- unique(as.character(.data[[f_col]]))
  sample_ids <- as.character(samples$sample)
  mat <- matrix(NA_real_, length(features), length(sample_ids),
                dimnames = list(features, sample_ids))
  mat[cbind(match(as.character(.data[[f_col]]), features),
            match(as.character(.data[[s_col]]), sample_ids))] <- .data[[v_col]]
  if (transform == "log2_cpm") {
    if (any(mat < 0, na.rm = TRUE)) {
      cli::cli_abort("Counts must be non-negative for {.code transform = \"log2_cpm\"}.")
    }
    lib <- colSums(mat, na.rm = TRUE)
    mat <- log2(t((t(mat) + 0.5) / (lib + 1) * 1e6))
  } else if (transform == "log2") {
    mat <- log2(mat + 1)
  }
  incomplete <- rowSums(!is.finite(mat)) > 0
  if (any(incomplete)) {
    cli::cli_warn("Dropping {sum(incomplete)} feature{?s} with missing or non-finite values.")
    mat <- mat[!incomplete, , drop = FALSE]
  }
  grp <- as.character(samples$group)
  sizes <- table(grp)
  if (any(sizes[c(reference, treatment)] < 2)) {
    cli::cli_abort("Each compared group needs at least 2 samples.")
  }
  # One-way model: residual variance pooled over all groups.
  means <- vapply(names(sizes), function(g) rowMeans(mat[, grp == g, drop = FALSE]),
                  numeric(nrow(mat)))
  means <- matrix(means, nrow = nrow(mat), dimnames = list(rownames(mat), names(sizes)))
  fitted <- means[, match(grp, names(sizes)), drop = FALSE]
  df_resid <- ncol(mat) - length(sizes)
  s2 <- rowSums((mat - fitted)^2) / df_resid
  if (robust_prior) {
    prior <- fit_f_dist(s2, df_resid)
  } else {
    prior <- list(df0 = 0, s0_sq = 0)
  }
  df_total <- if (is.finite(prior$df0)) df_resid + prior$df0 else Inf
  s2_post <- if (is.finite(prior$df0)) {
    (prior$df0 * prior$s0_sq + df_resid * s2) / (prior$df0 + df_resid)
  } else {
    rep(prior$s0_sq, length(s2))
  }
  df_total <- min(df_total, nrow(mat) * df_resid)
  lfc <- means[, treatment] - means[, reference]
  se <- sqrt(s2_post * (1 / sizes[[treatment]] + 1 / sizes[[reference]]))
  tstat <- lfc / se
  p <- 2 * stats::pt(-abs(tstat), df = df_total)
  out <- tibble::tibble(
    feature = rownames(mat), log2fc = unname(lfc), ave_expr = unname(rowMeans(mat)),
    t = unname(tstat), p_value = unname(p), p_adj = stats::p.adjust(p, method = "BH"),
    std_error = unname(se)
  )
  out <- out[order(out$p_value), ]
  meta <- list(
    df_prior = prior$df0, s0_sq = prior$s0_sq, df_residual = df_resid, df_total = df_total,
    treatment = treatment, reference = reference, transform = transform,
    n_samples = ncol(mat), moderated = robust_prior
  )
  new_bio_omics(out, "diff_expr", meta, required = c("feature", "log2fc", "p_value"))
}

#' Over-representation analysis of gene sets
#'
#' Tests whether a list of selected features (for example significant
#' genes) contains more members of each gene set or pathway than expected
#' by chance, using the one-sided hypergeometric test (Fisher's exact test
#' for over-representation). Only features in the universe are counted.
#' P-values are adjusted with the Benjamini-Hochberg method.
#'
#' @param .data A data frame with one row per feature, such as the output
#'   of [bio_diff_expr()]. All of its features form the universe unless
#'   `universe` is given.
#' @param feature <[`tidy-select`][dplyr::dplyr_tidy_select]> The feature
#'   identifier column.
#' @param selected <[`data-masking`][dplyr::dplyr_data_masking]> A logical
#'   expression choosing the features of interest, for example
#'   `p_adj < 0.05 & log2fc > 1`.
#' @param gene_sets Either a named list of character vectors or a data
#'   frame with columns `gene_set` and `feature`.
#' @param universe Optional character vector of background features.
#' @param min_size,max_size Gene sets with fewer or more members in the
#'   universe are skipped.
#'
#' @return A `<bio_omics>` sorted by p-value with columns `gene_set`,
#'   `set_size`, `overlap`, `expected`, `fold_enrichment`, `p_value`,
#'   `p_adj` and `overlap_features` (a list column).
#' @references Boyle, E. I. et al. (2004). GO::TermFinder: open source
#'   software for accessing Gene Ontology information and finding
#'   significantly enriched Gene Ontology terms associated with a list of
#'   genes. *Bioinformatics*, 20, 3710-3715.
#'   \doi{10.1093/bioinformatics/bth456}
#' @export
#' @examples
#' de <- bio_diff_expr(expression_sim, gene, sample, count, condition,
#'   treatment = "stress", reference = "control", transform = "log2_cpm"
#' )
#' bio_enrich(de, feature, p_adj < 0.05 & log2fc > 0, gene_sets_sim)
bio_enrich <- function(.data, feature, selected, gene_sets, universe = NULL, min_size = 5,
                       max_size = 500) {
  check_data_frame(.data)
  f_col <- select_one(.data, {{ feature }}, "feature")
  sel <- rlang::eval_tidy(rlang::enquo(selected), data = .data)
  if (!is.logical(sel) || length(sel) != nrow(.data)) {
    cli::cli_abort("{.arg selected} must be a logical expression with one value per row.")
  }
  sel[is.na(sel)] <- FALSE
  universe <- unique(as.character(universe %||% .data[[f_col]]))
  chosen <- intersect(unique(as.character(.data[[f_col]][sel])), universe)
  sets <- normalise_gene_sets(gene_sets)
  big_n <- length(universe)
  small_n <- length(chosen)
  if (small_n == 0L) {
    cli::cli_warn("No features were selected; every gene set gets p-value 1.")
  }
  rows <- lapply(names(sets), function(nm) {
    members <- intersect(unique(sets[[nm]]), universe)
    m <- length(members)
    hits <- intersect(members, chosen)
    k <- length(hits)
    tibble::tibble(
      gene_set = nm, set_size = m, overlap = k, expected = small_n * m / big_n,
      fold_enrichment = if (small_n > 0 && m > 0) (k / small_n) / (m / big_n) else NA_real_,
      p_value = stats::phyper(k - 1, m, big_n - m, small_n, lower.tail = FALSE),
      overlap_features = list(sort(hits))
    )
  })
  out <- vctrs::vec_rbind(!!!rows)
  out <- out[out$set_size >= min_size & out$set_size <= max_size, ]
  if (nrow(out) == 0L) {
    cli::cli_warn("No gene sets had between {min_size} and {max_size} members in the universe.")
  }
  out$p_adj <- stats::p.adjust(out$p_value, method = "BH")
  out <- out[order(out$p_value), c("gene_set", "set_size", "overlap", "expected",
                                   "fold_enrichment", "p_value", "p_adj", "overlap_features")]
  meta <- list(n_selected = small_n, universe_size = big_n, n_sets = nrow(out))
  new_bio_omics(out, "enrichment", meta, required = c("gene_set", "p_value"))
}

normalise_gene_sets <- function(gene_sets, call = rlang::caller_env()) {
  if (is.data.frame(gene_sets)) {
    if (!all(c("gene_set", "feature") %in% names(gene_sets))) {
      cli::cli_abort("A gene set data frame needs columns {.field gene_set} and {.field feature}.",
                     call = call)
    }
    return(split(as.character(gene_sets$feature), as.character(gene_sets$gene_set)))
  }
  if (is.list(gene_sets) && rlang::is_named(gene_sets)) {
    return(lapply(gene_sets, as.character))
  }
  cli::cli_abort("{.arg gene_sets} must be a named list or a data frame.", call = call)
}
