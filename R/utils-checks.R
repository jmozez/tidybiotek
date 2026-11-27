# Internal validation helpers built on cli, rlang and vctrs.

check_data_frame <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is.data.frame(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a data frame or tibble.",
        "x" = "You supplied {.obj_type_friendly {x}}."
      ),
      call = call
    )
  }
  invisible(x)
}

# Resolve a tidyselect specification that must pick exactly one column.
select_one <- function(.data, col, arg, call = rlang::caller_env()) {
  quo <- rlang::enquo(col)
  if (rlang::quo_is_missing(quo) || rlang::quo_is_null(quo)) {
    cli::cli_abort("{.arg {arg}} must name a column of {.arg .data}.", call = call)
  }
  loc <- tidyselect::eval_select(quo, .data, allow_rename = FALSE, error_call = call)
  if (length(loc) != 1L) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must select exactly one column.",
        "x" = "It selected {length(loc)} column{?s}."
      ),
      call = call
    )
  }
  names(loc)
}

# Resolve an optional single column; NULL gives NULL.
select_optional <- function(.data, col, arg, call = rlang::caller_env()) {
  quo <- rlang::enquo(col)
  if (rlang::quo_is_missing(quo) || rlang::quo_is_null(quo)) {
    return(NULL)
  }
  select_one(.data, !!quo, arg = arg, call = call)
}

# Resolve zero or more grouping columns.
select_by <- function(.data, by, call = rlang::caller_env()) {
  quo <- rlang::enquo(by)
  if (rlang::quo_is_missing(quo) || rlang::quo_is_null(quo)) {
    return(character())
  }
  names(tidyselect::eval_select(quo, .data, allow_rename = FALSE, error_call = call))
}

check_numeric_col <- function(.data, col, arg, call = rlang::caller_env()) {
  val <- .data[[col]]
  if (!is.numeric(val)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a numeric column.",
        "x" = "Column {.field {col}} is {.obj_type_friendly {val}}."
      ),
      call = call
    )
  }
  invisible(val)
}

check_character_col <- function(.data, col, arg, call = rlang::caller_env()) {
  x <- .data[[col]]
  if (!is.character(x) && !is.factor(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a character column of sequences.",
        "x" = "Column {.field {col}} is {.obj_type_friendly {x}}."
      ),
      call = call
    )
  }
  invisible(as.character(x))
}

check_positive_number <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env(),
                                  allow_zero = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
    (x > 0 || (allow_zero && x == 0))
  if (!ok) {
    what <- if (allow_zero) "a single non-negative number" else "a single positive number"
    cli::cli_abort("{.arg {arg}} must be {what}.", call = call)
  }
  invisible(x)
}

check_whole_number <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env(),
                               min = 1) {
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && x == round(x) && x >= min
  if (!ok) {
    cli::cli_abort("{.arg {arg}} must be a whole number of at least {min}.", call = call)
  }
  invisible(as.integer(x))
}

check_new_name <- function(.data, name, arg = rlang::caller_arg(name),
                           call = rlang::caller_env()) {
  if (!rlang::is_string(name) || !nzchar(name)) {
    cli::cli_abort("{.arg {arg}} must be a single, non-empty string.", call = call)
  }
  invisible(name)
}

# Clean sequences: upper case, strip whitespace, gaps and stop markers.
clean_seq <- function(x) {
  x <- toupper(as.character(x))
  gsub("[[:space:]-]", "", x)
}

guess_alphabet <- function(x) {
  x <- clean_seq(x[!is.na(x)])
  if (length(x) == 0L) {
    return("dna")
  }
  chars <- unique(unlist(strsplit(paste(x, collapse = ""), "", fixed = TRUE)))
  nucleotide <- c("A", "C", "G", "T", "U", "N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V")
  if (all(chars %in% nucleotide) && mean(chars %in% c("A", "C", "G", "T", "U", "N")) > 0.5) {
    if ("U" %in% chars && !"T" %in% chars) "rna" else "dna"
  } else {
    "protein"
  }
}

check_alphabet <- function(alphabet, call = rlang::caller_env()) {
  rlang::arg_match0(alphabet, c("dna", "rna", "protein"), arg_nm = "alphabet",
                    error_call = call)
}

# Resolve the sequence column, defaulting to the one stored by bio_seq_tbl().
resolve_seq_col <- function(.data, seq, call = rlang::caller_env()) {
  quo <- rlang::enquo(seq)
  if ((rlang::quo_is_missing(quo) || rlang::quo_is_null(quo)) && is_bio_seq_tbl(.data)) {
    return(attr(.data, "bio_meta")$seq_col)
  }
  select_one(.data, !!quo, arg = "seq", call = call)
}

resolve_alphabet <- function(.data, col, x) {
  if (is_bio_seq_tbl(.data) && identical(attr(.data, "bio_meta")$seq_col, col)) {
    return(attr(.data, "bio_meta")$alphabet)
  }
  guess_alphabet(x)
}

reverse_complement <- function(x) {
  x <- chartr("ACGTUNRYKMSWBDHV", "TGCAANYRMKSWVHDB", clean_seq(x))
  out <- vapply(strsplit(x, "", fixed = TRUE), function(s) paste(rev(s), collapse = ""),
                character(1))
  out[is.na(x)] <- NA_character_
  out
}

# Split row indices by key columns using vctrs.
split_rows <- function(.data, by) {
  n <- nrow(.data)
  if (length(by) == 0L) {
    return(list(key = tibble::new_tibble(list(), nrow = 1L), rows = list(seq_len(n))))
  }
  sp <- vctrs::vec_split(seq_len(n), tibble::as_tibble(.data)[by])
  list(key = tibble::as_tibble(sp$key), rows = sp$val)
}

# Bind the key columns onto per-group result rows.
bind_key <- function(key, rows) {
  out <- vctrs::vec_rbind(!!!rows)
  if (ncol(key) == 0L) {
    return(tibble::as_tibble(out))
  }
  tibble::as_tibble(vctrs::vec_cbind(key, out))
}

t_quantile <- function(level, df) {
  stats::qt(1 - (1 - level) / 2, df = pmax(df, 1))
}
