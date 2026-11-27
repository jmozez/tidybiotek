#' Pairwise sequence alignment
#'
#' @description
#' Aligns pairs of DNA, RNA or protein sequences by dynamic programming with
#' affine gap penalties (Gotoh 1982): global alignment over the whole length
#' (Needleman and Wunsch 1970) or local alignment of the best-matching
#' regions (Smith and Waterman 1981). Proteins are scored with the BLOSUM62
#' substitution matrix (Henikoff and Henikoff 1992) by default, nucleotides
#' with a match and mismatch score. A gap of length `k` costs
#' `gap_open + k * gap_extend`, the convention of `pwalign::pairwiseAlignment()`,
#' whose scores this function reproduces.
#'
#' The alignment is written in base R, so it suits genes and proteins up to
#' a few thousand residues.
#'
#' @param pattern,subject Character vectors of sequences, recycled to a
#'   common length; each pair is aligned.
#' @param type `"global"` or `"local"`.
#' @param substitution `"blosum62"` (default for proteins) or `"nucleotide"`
#'   (match and mismatch scores; default for DNA and RNA).
#' @param match,mismatch Scores for nucleotide matches and mismatches.
#' @param gap_open,gap_extend Gap opening and extension penalties (positive
#'   numbers). Defaults are 10 and 4 for nucleotides and 10 and 0.5 for
#'   proteins.
#' @return A tibble with one row per pair: `score`, `aligned_pattern`,
#'   `aligned_subject` (with `-` for gaps), `identity` (identical positions
#'   over the alignment length), `length`, `gaps`, and `pattern_start`,
#'   `pattern_end`, `subject_start`, `subject_end` (the aligned region).
#' @references Gotoh, O. (1982). An improved algorithm for matching
#'   biological sequences. Journal of Molecular Biology, 162, 705-708.
#'   \doi{10.1016/0022-2836(82)90398-9}
#'
#'   Needleman, S. B. and Wunsch, C. D. (1970). A general method applicable
#'   to the search for similarities in the amino acid sequence of two
#'   proteins. Journal of Molecular Biology, 48, 443-453.
#'   \doi{10.1016/0022-2836(70)90057-4}
#'
#'   Smith, T. F. and Waterman, M. S. (1981). Identification of common
#'   molecular subsequences. Journal of Molecular Biology, 147, 195-197.
#'   \doi{10.1016/0022-2836(81)90087-5}
#'
#'   Henikoff, S. and Henikoff, J. G. (1992). Amino acid substitution
#'   matrices from protein blocks. Proceedings of the National Academy of
#'   Sciences, 89, 10915-10919. \doi{10.1073/pnas.89.22.10915}
#' @export
#' @examples
#' bio_align("HEAGAWGHEE", "PAWHEAE")
#' bio_align("HEAGAWGHEE", "PAWHEAE", type = "local")
#'
#' # The first 60 bases of two E. coli genes
#' cds <- substr(ecoli_cds$sequence[1:2], 1, 60)
#' bio_align(cds[1], cds[2])
bio_align <- function(pattern, subject, type = c("global", "local"),
                      substitution = NULL, match = 2, mismatch = -3,
                      gap_open = 10, gap_extend = NULL) {
  type <- rlang::arg_match(type)
  if (!is.character(pattern) || !is.character(subject) || anyNA(pattern) || anyNA(subject)) {
    cli::cli_abort(
      "{.arg pattern} and {.arg subject} must be character vectors without missing values."
    )
  }
  if (!length(pattern) || !length(subject)) {
    cli::cli_abort("{.arg pattern} and {.arg subject} must not be empty.")
  }
  pattern <- clean_seq(pattern)
  subject <- clean_seq(subject)
  alpha <- guess_alphabet(c(pattern, subject))
  substitution <- substitution %||% if (alpha == "protein") "blosum62" else "nucleotide"
  substitution <- rlang::arg_match(substitution, c("blosum62", "nucleotide"))
  gap_extend <- gap_extend %||% if (substitution == "blosum62") 0.5 else 4
  check_positive_number(gap_open, allow_zero = TRUE)
  check_positive_number(gap_extend, allow_zero = TRUE)
  n <- max(length(pattern), length(subject))
  pattern <- rep_len(pattern, n)
  subject <- rep_len(subject, n)
  rows <- lapply(seq_len(n), function(i) {
    align_pair(pattern[[i]], subject[[i]], type, substitution, match, mismatch,
               gap_open, gap_extend)
  })
  dplyr::bind_rows(rows)
}

align_score_matrix <- function(a, b, substitution, match, mismatch) {
  if (substitution == "nucleotide") {
    return(ifelse(outer(a, b, "=="), match, mismatch))
  }
  m <- blosum62()
  ia <- match(a, rownames(m))
  ib <- match(b, colnames(m))
  ia[is.na(ia)] <- match("X", rownames(m))
  ib[is.na(ib)] <- match("X", colnames(m))
  m[ia, ib, drop = FALSE]
}

align_pair <- function(x, y, type, substitution, match, mismatch, go, ge) {
  a <- strsplit(x, "")[[1]]
  b <- strsplit(y, "")[[1]]
  n <- length(a)
  m <- length(b)
  if (n == 0 || m == 0) cli::cli_abort("Sequences must not be empty.")
  s <- align_score_matrix(a, b, substitution, match, mismatch)
  local <- type == "local"
  ninf <- -Inf
  mm <- matrix(ninf, n + 1, m + 1)
  xx <- matrix(ninf, n + 1, m + 1)
  yy <- matrix(ninf, n + 1, m + 1)
  pm <- matrix(0L, n + 1, m + 1)
  px <- matrix(0L, n + 1, m + 1)
  py <- matrix(0L, n + 1, m + 1)
  mm[1, 1] <- 0
  if (local) {
    mm[, 1] <- 0
    mm[1, ] <- 0
  } else {
    xx[2:(n + 1), 1] <- -(go + seq_len(n) * ge)
    yy[1, 2:(m + 1)] <- -(go + seq_len(m) * ge)
    px[2:(n + 1), 1] <- 2L
    px[2, 1] <- 1L
    py[1, 2:(m + 1)] <- 3L
    py[1, 2] <- 1L
  }
  for (i in 2:(n + 1)) {
    for (j in 2:(m + 1)) {
      cand <- c(mm[i - 1, j - 1], xx[i - 1, j - 1], yy[i - 1, j - 1])
      k <- which.max(cand)
      v <- cand[k] + s[i - 1, j - 1]
      if (local && v <= 0) {
        mm[i, j] <- 0
        pm[i, j] <- 0L
      } else {
        mm[i, j] <- v
        pm[i, j] <- k
      }
      cx <- c(mm[i - 1, j] - go - ge, xx[i - 1, j] - ge, yy[i - 1, j] - go - ge)
      k <- which.max(cx)
      xx[i, j] <- cx[k]
      px[i, j] <- k
      cy <- c(mm[i, j - 1] - go - ge, xx[i, j - 1] - go - ge, yy[i, j - 1] - ge)
      k <- which.max(cy)
      yy[i, j] <- cy[k]
      py[i, j] <- k
    }
  }
  if (local) {
    best <- which(mm == max(mm), arr.ind = TRUE)
    best <- best[order(best[, 1], best[, 2]), , drop = FALSE][1, ]
    i <- best[[1]]
    j <- best[[2]]
    state <- 1L
    score <- mm[i, j]
  } else {
    i <- n + 1
    j <- m + 1
    fin <- c(mm[i, j], xx[i, j], yy[i, j])
    state <- which.max(fin)
    score <- fin[state]
  }
  end_i <- i - 1
  end_j <- j - 1
  ra <- character()
  rb <- character()
  repeat {
    if (local && state == 1L && (mm[i, j] == 0 || pm[i, j] == 0L)) break
    if (!local && i == 1 && j == 1) break
    if (state == 1L) {
      ra <- c(a[i - 1], ra)
      rb <- c(b[j - 1], rb)
      nxt <- pm[i, j]
      i <- i - 1
      j <- j - 1
    } else if (state == 2L) {
      ra <- c(a[i - 1], ra)
      rb <- c("-", rb)
      nxt <- px[i, j]
      i <- i - 1
    } else {
      ra <- c("-", ra)
      rb <- c(b[j - 1], rb)
      nxt <- py[i, j]
      j <- j - 1
    }
    state <- if (nxt == 0L) 1L else nxt
    if (!local && (i == 1 && j == 1)) break
  }
  len <- length(ra)
  tibble::tibble(
    score = score,
    aligned_pattern = paste(ra, collapse = ""),
    aligned_subject = paste(rb, collapse = ""),
    identity = if (len > 0) sum(ra == rb & ra != "-") / len else NA_real_,
    length = len,
    gaps = sum(ra == "-") + sum(rb == "-"),
    pattern_start = as.integer(end_i - sum(ra != "-") + 1),
    pattern_end = as.integer(end_i),
    subject_start = as.integer(end_j - sum(rb != "-") + 1),
    subject_end = as.integer(end_j)
  )
}

# BLOSUM62 (Henikoff and Henikoff 1992), as distributed by NCBI.
blosum62 <- function() {
  l <- strsplit(blosum62_letters, "")[[1]]
  matrix(blosum62_values, length(l), length(l), byrow = TRUE, dimnames = list(l, l))
}

blosum62_letters <- "ARNDCQEGHILKMFPSTWYVBZX*"
blosum62_values <- c(
  4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, -1, -4,
  -1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4,
  -2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 4, 0, -1, -4,
  -2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4,
  0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -1, -4,
  -1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 4, -1, -4,
  -1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4,
  0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4,
  -2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4,
  -1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4,
  -1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4,
  -1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4,
  -1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4,
  -2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4,
  -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -1, -4,
  1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, -1, -4,
  0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, -1, -4,
  -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -2, -1, -4,
  -2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4,
  0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4,
  -2, -1, 4, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 0, -1, -4,
  -1, 0, 0, 1, -3, 4, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -2, -2, -2, 0, 4, -1, -4,
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4,
  -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1
)
