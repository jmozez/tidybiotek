# Resampling schemes that respect sequence homology and batch structure.

make_group_rset <- function(.data, fold_of_row, ids, attrib, subclass) {
  n <- nrow(.data)
  folds <- sort(unique(fold_of_row))
  splits <- lapply(folds, function(f) {
    assess <- which(fold_of_row == f)
    rsample::make_splits(
      list(analysis = setdiff(seq_len(n), assess), assessment = assess),
      data = .data
    )
  })
  rsample::new_rset(
    splits = splits,
    ids = tibble::tibble(id = ids),
    attrib = attrib,
    subclass = c(subclass, "rset")
  )
}

# Assign clusters to folds: largest clusters first, each to the fold that
# currently holds the fewest rows.
balance_clusters <- function(cluster, v) {
  sizes <- sort(table(cluster), decreasing = TRUE)
  load <- numeric(v)
  fold_of_cluster <- integer(length(sizes))
  names(fold_of_cluster) <- names(sizes)
  for (cl in names(sizes)) {
    f <- which.min(load)
    fold_of_cluster[[cl]] <- f
    load[[f]] <- load[[f]] + sizes[[cl]]
  }
  unname(fold_of_cluster[as.character(cluster)])
}

#' Homology-aware cross-validation for sequences
#'
#' Similar sequences (homologues) in both the analysis and assessment sets
#' make cross-validation over-optimistic. `bio_homology_cv()` first
#' clusters the sequences by k-mer profile similarity and then assigns
#' whole clusters to folds, so no cluster is split between training and
#' testing. Sequences are clustered by average-linkage (UPGMA) hierarchical
#' clustering on the cosine distance between k-mer frequency vectors, cut
#' at `threshold`. The k-mer distance is an alignment-free stand-in for
#' sequence identity: smaller thresholds give tighter clusters.
#'
#' @param .data A data frame, tibble or [bio_seq_tbl()].
#' @param seq <[`tidy-select`][dplyr::dplyr_tidy_select]> The sequence
#'   column. Can be omitted for a `<bio_seq_tbl>`.
#' @param v Number of folds.
#' @param k K-mer length used for the similarity.
#' @param threshold Cosine distance at which the clustering tree is cut,
#'   between 0 and 1.
#'
#' @return An `rset` of class `bio_homology_cv` with columns `splits` and
#'   `id`. The cluster of every row is kept in `attr(x, "clusters")`.
#' @references Walsh, I. et al. (2021). DOME: recommendations for supervised
#'   machine learning validation in biology. *Nature Methods*, 18,
#'   1122-1127. \doi{10.1038/s41592-021-01205-4}
#' @export
#' @examples
#' set.seed(1)
#' folds <- bio_homology_cv(ecoli_cds, protein, v = 5, k = 2, threshold = 0.2)
#' folds
#' table(attr(folds, "clusters"))[1:5]
bio_homology_cv <- function(.data, seq = NULL, v = 5, k = 3, threshold = 0.3) {
  check_data_frame(.data)
  col <- resolve_seq_col(.data, {{ seq }})
  x <- clean_seq(check_character_col(.data, col, "seq"))
  v <- check_whole_number(v, min = 2)
  k <- check_whole_number(k)
  if (!is.numeric(threshold) || length(threshold) != 1L || threshold < 0 || threshold > 1) {
    cli::cli_abort("{.arg threshold} must be a single number between 0 and 1.")
  }
  if (anyNA(x)) {
    cli::cli_abort("Sequences must not be missing.")
  }
  vocab <- kmer_vocabulary(x, k, resolve_alphabet(.data, col, x))
  m <- kmer_matrix(x, k, vocab, normalise = TRUE)
  norms <- sqrt(rowSums(m^2))
  norms[norms == 0] <- 1
  unit <- m / norms
  sim <- tcrossprod(unit)
  d <- stats::as.dist(pmax(1 - sim, 0))
  clusters <- if (nrow(m) > 1) {
    stats::cutree(stats::hclust(d, method = "average"), h = threshold)
  } else {
    1L
  }
  n_clusters <- length(unique(clusters))
  if (n_clusters < v) {
    cli::cli_abort(
      c(
        "Found {n_clusters} homology cluster{?s}, fewer than the {v} folds requested.",
        "i" = "Lower {.arg threshold} or {.arg v}."
      )
    )
  }
  # Randomise the order of equal-sized clusters before balancing.
  relabel <- sample(unique(clusters))
  clusters <- match(clusters, relabel)
  fold <- balance_clusters(clusters, v)
  ids <- paste0("Fold", formatC(seq_len(v), width = nchar(v), flag = "0"))
  out <- make_group_rset(.data, fold, ids,
                         attrib = list(v = v, k = k, threshold = threshold),
                         subclass = "bio_homology_cv")
  attr(out, "clusters") <- clusters
  out
}

#' Leave-batch-out cross-validation
#'
#' Batch effects (plates, runs, sequencing lanes, fermentation campaigns)
#' are shared by samples processed together. `bio_batch_cv()` keeps every
#' batch whole: each assessment set is one batch (leave-one-batch-out) or,
#' when `v` is given, a group of batches balanced by size.
#'
#' @param .data A data frame or tibble.
#' @param batch <[`tidy-select`][dplyr::dplyr_tidy_select]> Column
#'   identifying the batch of each row.
#' @param v Number of folds. `NULL` (the default) gives one fold per batch.
#'
#' @return An `rset` of class `bio_batch_cv` with columns `splits` and `id`.
#' @references Leek, J. T. et al. (2010). Tackling the widespread and
#'   critical impact of batch effects in high-throughput data. *Nature
#'   Reviews Genetics*, 11, 733-739. \doi{10.1038/nrg2825}
#' @export
#' @examples
#' folds <- bio_batch_cv(qpcr_sim, plate)
#' folds
#' rsample::assessment(folds$splits[[1]])$plate |> unique()
bio_batch_cv <- function(.data, batch, v = NULL) {
  check_data_frame(.data)
  b_col <- select_one(.data, {{ batch }}, "batch")
  b <- as.character(.data[[b_col]])
  if (anyNA(b)) {
    cli::cli_abort("The batch column must not contain missing values.")
  }
  batches <- unique(b)
  if (length(batches) < 2L) {
    cli::cli_abort("Leave-batch-out resampling needs at least 2 batches.")
  }
  if (is.null(v)) {
    fold <- match(b, batches)
    ids <- paste0("Batch_", batches)
    v <- length(batches)
  } else {
    v <- check_whole_number(v, min = 2)
    if (v > length(batches)) {
      cli::cli_abort("{.arg v} ({v}) cannot exceed the number of batches ({length(batches)}).")
    }
    fold <- balance_clusters(b, v)
    ids <- paste0("Fold", formatC(seq_len(v), width = nchar(v), flag = "0"))
  }
  make_group_rset(.data, fold, ids, attrib = list(v = v, batch = b_col),
                  subclass = "bio_batch_cv")
}

#' @export
pretty.bio_homology_cv <- function(x, ...) {
  paste0(attr(x, "v"), "-fold homology-aware cross-validation (k = ", attr(x, "k"),
         ", threshold = ", attr(x, "threshold"), ")")
}

#' @export
pretty.bio_batch_cv <- function(x, ...) {
  paste0("Leave-batch-out cross-validation (", attr(x, "v"), " folds by `",
         attr(x, "batch"), "`)")
}

#' @export
print.bio_homology_cv <- function(x, ...) {
  cat("#", pretty(x), "\n")
  class(x) <- class(x)[!(class(x) %in% c("bio_homology_cv", "rset"))]
  print(x, ...)
  invisible(x)
}

#' @export
print.bio_batch_cv <- function(x, ...) {
  cat("#", pretty(x), "\n")
  class(x) <- class(x)[!(class(x) %in% c("bio_batch_cv", "rset"))]
  print(x, ...)
  invisible(x)
}
