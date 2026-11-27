# Tibble subclasses for sequences, assays and omics results.
#
# Each class stores its metadata in two attributes: `bio_meta` (a list) and
# `bio_required` (the columns that must survive for the class to stay
# meaningful). The vctrs and dplyr methods below keep the class and its
# metadata through subsetting, binding and dplyr verbs, and quietly fall
# back to a plain tibble once a required column is dropped.

new_bio_tbl <- function(x, meta, required, class) {
  x <- tibble::as_tibble(strip_bio(x))
  tibble::new_tibble(x, bio_meta = meta, bio_required = required, class = class)
}

strip_bio <- function(x) {
  if (!is.data.frame(x)) {
    return(x)
  }
  attr(x, "bio_meta") <- NULL
  attr(x, "bio_required") <- NULL
  class(x) <- setdiff(class(x), c("bio_seq_tbl", "bio_assay", "bio_omics"))
  x
}

bio_restore <- function(x, to, class) {
  required <- attr(to, "bio_required")
  if (!is.data.frame(x) || !all(required %in% names(x))) {
    return(tibble::as_tibble(strip_bio(x)))
  }
  new_bio_tbl(x, attr(to, "bio_meta"), required, class)
}

same_meta <- function(x, y) {
  identical(attr(x, "bio_meta"), attr(y, "bio_meta")) &&
    identical(attr(x, "bio_required"), attr(y, "bio_required"))
}

bio_ptype2_self <- function(x, y, class, ...) {
  if (same_meta(x, y)) {
    return(new_bio_tbl(vctrs::df_ptype2(x, y, ...), attr(x, "bio_meta"),
                       attr(x, "bio_required"), class))
  }
  vctrs::tib_ptype2(strip_bio(x), strip_bio(y), ...)
}

bio_cast_self <- function(x, to, class, ...) {
  new_bio_tbl(vctrs::df_cast(strip_bio(x), strip_bio(to), ...), attr(to, "bio_meta"),
              attr(to, "bio_required"), class)
}

# ---- bio_seq_tbl ------------------------------------------------------------

#' Sequence tibbles
#'
#' `bio_seq_tbl()` marks one column of a data frame as holding biological
#' sequences and records their alphabet. The sequences are cleaned (upper
#' case, whitespace and gap characters removed) and the alphabet is kept
#' through `dplyr` verbs, so later functions such as [bio_translate()] and
#' [bio_kmer_count()] do not need the column named again.
#' `new_bio_seq_tbl()` is the low-level constructor for developers.
#'
#' @param .data A data frame or tibble.
#' @param seq <[`tidy-select`][dplyr::dplyr_tidy_select]> The column holding
#'   sequences.
#' @param alphabet One of `"dna"`, `"rna"` or `"protein"`. `NULL` (the
#'   default) guesses it from the sequences.
#' @param x For `new_bio_seq_tbl()`, a data frame. For `is_bio_seq_tbl()` and
#'   `bio_alphabet()`, any object.
#' @param seq_col For `new_bio_seq_tbl()`, the name of the sequence column.
#'
#' @return `bio_seq_tbl()` and `new_bio_seq_tbl()` return a `<bio_seq_tbl>`,
#'   a tibble subclass. `is_bio_seq_tbl()` returns a single logical and
#'   `bio_alphabet()` returns the alphabet string, or `NULL`.
#' @export
#' @examples
#' seqs <- bio_seq_tbl(ecoli_cds, sequence)
#' seqs
#' bio_alphabet(seqs)
#'
#' # The class survives dplyr verbs that keep the sequence column
#' seqs |>
#'   dplyr::filter(class == "ribosomal") |>
#'   is_bio_seq_tbl()
bio_seq_tbl <- function(.data, seq, alphabet = NULL) {
  check_data_frame(.data)
  col <- select_one(.data, {{ seq }}, arg = "seq")
  x <- check_character_col(.data, col, "seq")
  x <- clean_seq(x)
  alphabet <- if (is.null(alphabet)) guess_alphabet(x) else check_alphabet(alphabet)
  allowed <- switch(alphabet,
    dna = "^[ACGTNRYSWKMBDHV]*$",
    rna = "^[ACGUNRYSWKMBDHV]*$",
    protein = "^[ACDEFGHIKLMNPQRSTVWYXBZJUO*]*$"
  )
  bad <- which(!is.na(x) & !grepl(allowed, x))
  if (length(bad) > 0L) {
    cli::cli_abort(
      c(
        "Column {.field {col}} contains characters that are not valid {alphabet}.",
        "x" = "Problem in row{?s} {utils::head(bad, 5)}."
      )
    )
  }
  .data[[col]] <- x
  new_bio_seq_tbl(.data, seq_col = col, alphabet = alphabet)
}

#' @rdname bio_seq_tbl
#' @export
new_bio_seq_tbl <- function(x = tibble::tibble(seq = character()), seq_col = "seq",
                            alphabet = "dna") {
  check_data_frame(x)
  if (!rlang::is_string(seq_col) || !seq_col %in% names(x)) {
    cli::cli_abort("{.arg seq_col} must name a column of {.arg x}.")
  }
  alphabet <- check_alphabet(alphabet)
  new_bio_tbl(x, list(seq_col = seq_col, alphabet = alphabet), seq_col, "bio_seq_tbl")
}

#' @rdname bio_seq_tbl
#' @export
is_bio_seq_tbl <- function(x) {
  inherits(x, "bio_seq_tbl")
}

#' @rdname bio_seq_tbl
#' @export
bio_alphabet <- function(x) {
  if (!is_bio_seq_tbl(x)) {
    return(NULL)
  }
  attr(x, "bio_meta")$alphabet
}

# ---- bio_assay --------------------------------------------------------------

#' Assay and omics result tibbles
#'
#' Results from the assay, bioprocess and omics functions are tibble
#' subclasses. A `<bio_assay>` carries its fitted standard curves or
#' kinetic models, and a `<bio_omics>` carries the settings and prior
#' estimates of a differential expression or enrichment analysis. The
#' metadata is used by [tidy()][tidy.bio_assay], [glance()][glance.bio_assay],
#' [augment()][augment.bio_assay] and [autoplot()][autoplot.bio_assay], and is
#' kept through `dplyr` verbs while the required columns are present.
#'
#' @param x A data frame (constructors) or any object (predicates).
#' @param assay,analysis A string naming the kind of result, for example
#'   `"dose_response"` or `"diff_expr"`.
#' @param meta A list of metadata.
#' @param required Columns that must be present for the class to be kept.
#'
#' @return The constructors return a `<bio_assay>` or `<bio_omics>` tibble.
#'   The predicates return a single logical.
#' @name bio_results
#' @examples
#' fit <- bio_dose_response(dose_response_sim, dose, response, .by = compound)
#' is_bio_assay(fit)
#' fit |> dplyr::arrange(ec50) |> is_bio_assay()
#'
#' x <- new_bio_assay(tibble::tibble(id = 1:2), "custom", required = "id")
#' x
NULL

#' @rdname bio_results
#' @export
new_bio_assay <- function(x, assay, meta = list(), required = character()) {
  check_data_frame(x)
  if (!rlang::is_string(assay)) {
    cli::cli_abort("{.arg assay} must be a single string.")
  }
  meta$assay <- assay
  new_bio_tbl(x, meta, required, "bio_assay")
}

#' @rdname bio_results
#' @export
is_bio_assay <- function(x) {
  inherits(x, "bio_assay")
}

#' @rdname bio_results
#' @export
new_bio_omics <- function(x, analysis, meta = list(), required = character()) {
  check_data_frame(x)
  if (!rlang::is_string(analysis)) {
    cli::cli_abort("{.arg analysis} must be a single string.")
  }
  meta$analysis <- analysis
  new_bio_tbl(x, meta, required, "bio_omics")
}

#' @rdname bio_results
#' @export
is_bio_omics <- function(x) {
  inherits(x, "bio_omics")
}

bio_meta <- function(x) {
  attr(x, "bio_meta")
}

# ---- printing ---------------------------------------------------------------

#' @exportS3Method pillar::tbl_sum
tbl_sum.bio_seq_tbl <- function(x, ...) {
  meta <- bio_meta(x)
  label <- c(dna = "DNA", rna = "RNA", protein = "protein")[[meta$alphabet]]
  c("A bio_seq_tbl" = paste0(nrow(x), " sequences (", label, ") in `", meta$seq_col, "`"))
}

#' @exportS3Method pillar::tbl_sum
tbl_sum.bio_assay <- function(x, ...) {
  meta <- bio_meta(x)
  c("A bio_assay" = paste0(meta$assay, ": ", nrow(x), " x ", ncol(x)))
}

#' @exportS3Method pillar::tbl_sum
tbl_sum.bio_omics <- function(x, ...) {
  meta <- bio_meta(x)
  c("A bio_omics" = paste0(meta$analysis, ": ", nrow(x), " x ", ncol(x)))
}

# ---- subsetting and dplyr ---------------------------------------------------

#' @export
`[.bio_seq_tbl` <- function(x, i, j, drop = FALSE, ...) {
  bio_restore(NextMethod(), x, "bio_seq_tbl")
}

#' @export
`[.bio_assay` <- function(x, i, j, drop = FALSE, ...) {
  bio_restore(NextMethod(), x, "bio_assay")
}

#' @export
`[.bio_omics` <- function(x, i, j, drop = FALSE, ...) {
  bio_restore(NextMethod(), x, "bio_omics")
}

#' @exportS3Method dplyr::dplyr_reconstruct
dplyr_reconstruct.bio_seq_tbl <- function(data, template) {
  bio_restore(data, template, "bio_seq_tbl")
}

#' @exportS3Method dplyr::dplyr_reconstruct
dplyr_reconstruct.bio_assay <- function(data, template) {
  bio_restore(data, template, "bio_assay")
}

#' @exportS3Method dplyr::dplyr_reconstruct
dplyr_reconstruct.bio_omics <- function(data, template) {
  bio_restore(data, template, "bio_omics")
}

# ---- vctrs methods ----------------------------------------------------------

#' vctrs methods for tidybiotek classes
#'
#' These methods let `<bio_seq_tbl>`, `<bio_assay>` and `<bio_omics>`
#' objects work with `vctrs::vec_rbind()`, `dplyr::bind_rows()` and other
#' vctrs-based tools. Two objects with the same metadata combine into the
#' same class; anything else combines into a plain tibble.
#'
#' @param x,y,to Objects to restore, combine or cast.
#' @param ... Passed on to the vctrs data frame methods.
#' @return A tibble subclass or a plain tibble.
#' @keywords internal
#' @name tidybiotek-vctrs
NULL

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_restore
vec_restore.bio_seq_tbl <- function(x, to, ...) bio_restore(x, to, "bio_seq_tbl")

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_restore
vec_restore.bio_assay <- function(x, to, ...) bio_restore(x, to, "bio_assay")

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_restore
vec_restore.bio_omics <- function(x, to, ...) bio_restore(x, to, "bio_omics")

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.bio_seq_tbl.bio_seq_tbl <- function(x, y, ...) {
  bio_ptype2_self(x, y, "bio_seq_tbl", ...)
}

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.bio_seq_tbl.tbl_df <- function(x, y, ...) vctrs::tib_ptype2(strip_bio(x), y, ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.tbl_df.bio_seq_tbl <- function(x, y, ...) vctrs::tib_ptype2(x, strip_bio(y), ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.bio_seq_tbl.data.frame <- function(x, y, ...) {
  vctrs::tib_ptype2(strip_bio(x), y, ...)
}

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.data.frame.bio_seq_tbl <- function(x, y, ...) {
  vctrs::tib_ptype2(x, strip_bio(y), ...)
}

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_cast
vec_cast.bio_seq_tbl.bio_seq_tbl <- function(x, to, ...) bio_cast_self(x, to, "bio_seq_tbl", ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_cast
vec_cast.bio_seq_tbl.tbl_df <- function(x, to, ...) bio_cast_self(x, to, "bio_seq_tbl", ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_cast
vec_cast.tbl_df.bio_seq_tbl <- function(x, to, ...) vctrs::tib_cast(strip_bio(x), to, ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_cast
vec_cast.data.frame.bio_seq_tbl <- function(x, to, ...) vctrs::df_cast(strip_bio(x), to, ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.bio_assay.bio_assay <- function(x, y, ...) bio_ptype2_self(x, y, "bio_assay", ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.bio_assay.tbl_df <- function(x, y, ...) vctrs::tib_ptype2(strip_bio(x), y, ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.tbl_df.bio_assay <- function(x, y, ...) vctrs::tib_ptype2(x, strip_bio(y), ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.bio_assay.data.frame <- function(x, y, ...) vctrs::tib_ptype2(strip_bio(x), y, ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.data.frame.bio_assay <- function(x, y, ...) vctrs::tib_ptype2(x, strip_bio(y), ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_cast
vec_cast.bio_assay.bio_assay <- function(x, to, ...) bio_cast_self(x, to, "bio_assay", ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_cast
vec_cast.bio_assay.tbl_df <- function(x, to, ...) bio_cast_self(x, to, "bio_assay", ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_cast
vec_cast.tbl_df.bio_assay <- function(x, to, ...) vctrs::tib_cast(strip_bio(x), to, ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_cast
vec_cast.data.frame.bio_assay <- function(x, to, ...) vctrs::df_cast(strip_bio(x), to, ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.bio_omics.bio_omics <- function(x, y, ...) bio_ptype2_self(x, y, "bio_omics", ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.bio_omics.tbl_df <- function(x, y, ...) vctrs::tib_ptype2(strip_bio(x), y, ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.tbl_df.bio_omics <- function(x, y, ...) vctrs::tib_ptype2(x, strip_bio(y), ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.bio_omics.data.frame <- function(x, y, ...) vctrs::tib_ptype2(strip_bio(x), y, ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.data.frame.bio_omics <- function(x, y, ...) vctrs::tib_ptype2(x, strip_bio(y), ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_cast
vec_cast.bio_omics.bio_omics <- function(x, to, ...) bio_cast_self(x, to, "bio_omics", ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_cast
vec_cast.bio_omics.tbl_df <- function(x, to, ...) bio_cast_self(x, to, "bio_omics", ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_cast
vec_cast.tbl_df.bio_omics <- function(x, to, ...) vctrs::tib_cast(strip_bio(x), to, ...)

#' @rdname tidybiotek-vctrs
#' @exportS3Method vctrs::vec_cast
vec_cast.data.frame.bio_omics <- function(x, to, ...) vctrs::df_cast(strip_bio(x), to, ...)
