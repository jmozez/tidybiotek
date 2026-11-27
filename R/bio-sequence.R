# Sequence manipulation: translation, GC content, codon usage, primer Tm
# and k-mer counts.

genetic_codes <- list(
  standard = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  bacterial = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  vertebrate_mito = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
  yeast_mito = "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
)

start_codons <- list(
  standard = c("TTG", "CTG", "ATG"),
  bacterial = c("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"),
  vertebrate_mito = c("ATT", "ATC", "ATA", "ATG", "GTG"),
  yeast_mito = c("ATA", "ATG", "GTG")
)

all_codons <- function() {
  b <- c("T", "C", "A", "G")
  paste0(rep(b, each = 16), rep(rep(b, each = 4), 4), rep(b, 16))
}

codon_table <- function(code = "standard") {
  aa <- strsplit(genetic_codes[[code]], "", fixed = TRUE)[[1]]
  stats::setNames(aa, all_codons())
}

split_codons <- function(s, frame = 1L) {
  n <- nchar(s)
  if (n - frame + 1L < 3L) {
    return(character())
  }
  starts <- seq.int(frame, n - 2L, by = 3L)
  substring(s, starts, starts + 2L)
}

#' Translate nucleotide sequences to protein
#'
#' Translates DNA or RNA sequences with an NCBI genetic code. Codons that
#' contain ambiguous bases become `X` and stop codons become `*`.
#'
#' @param .data A data frame, tibble or [bio_seq_tbl()].
#' @param seq <[`tidy-select`][dplyr::dplyr_tidy_select]> The sequence
#'   column. Can be omitted for a `<bio_seq_tbl>`.
#' @param frame Reading frame: `1`, `2` or `3` on the given strand, or `-1`,
#'   `-2` or `-3` on the reverse complement.
#' @param code Genetic code: `"standard"` (NCBI table 1), `"bacterial"`
#'   (table 11), `"vertebrate_mito"` (table 2) or `"yeast_mito"` (table 3).
#' @param to_stop If `TRUE`, stop translation at the first stop codon and
#'   drop it.
#' @param start_as_met If `TRUE`, translate an alternative start codon in the
#'   first position as methionine, as happens in vivo.
#' @param .into Name of the new protein column.
#'
#' @return `.data` with a new character column `.into`. A `<bio_seq_tbl>`
#'   keeps its class.
#' @references NCBI. The Genetic Codes.
#'   <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
#' @export
#' @examples
#' ecoli_cds |>
#'   bio_translate(sequence, code = "bacterial", to_stop = TRUE, start_as_met = TRUE) |>
#'   dplyr::mutate(matches = translated == protein) |>
#'   dplyr::count(matches)
bio_translate <- function(.data, seq = NULL, frame = 1, code = "standard", to_stop = FALSE,
                          start_as_met = FALSE, .into = "translated") {
  check_data_frame(.data)
  col <- resolve_seq_col(.data, {{ seq }})
  x <- clean_seq(check_character_col(.data, col, "seq"))
  code <- rlang::arg_match0(code, names(genetic_codes))
  if (!frame %in% c(1, 2, 3, -1, -2, -3)) {
    cli::cli_abort("{.arg frame} must be one of 1, 2, 3, -1, -2 or -3.")
  }
  check_new_name(.data, .into)
  x <- gsub("U", "T", x, fixed = TRUE)
  if (frame < 0) {
    x <- reverse_complement(x)
  }
  table <- codon_table(code)
  starts <- start_codons[[code]]
  .data[[.into]] <- vapply(x, function(s) {
    if (is.na(s)) {
      return(NA_character_)
    }
    codons <- split_codons(s, abs(frame))
    if (length(codons) == 0L) {
      return("")
    }
    aa <- unname(table[codons])
    aa[is.na(aa)] <- "X"
    if (start_as_met && codons[[1]] %in% starts) {
      aa[[1]] <- "M"
    }
    if (to_stop) {
      stop_at <- match("*", aa)
      if (!is.na(stop_at)) {
        aa <- aa[seq_len(stop_at - 1L)]
      }
    }
    paste(aa, collapse = "")
  }, character(1), USE.NAMES = FALSE)
  .data
}

#' GC content of nucleotide sequences
#'
#' Calculates the fraction of G and C among unambiguous bases (A, C, G and
#' T or U). With `by_position = TRUE` it also returns GC at the first,
#' second and third codon positions (GC1, GC2 and GC3), which are widely
#' used to compare codon usage between genes and genomes.
#'
#' @inheritParams bio_translate
#' @param by_position If `TRUE`, add `gc1`, `gc2` and `gc3` columns as well.
#' @param .into Name of the overall GC column.
#'
#' @return `.data` with a numeric column `.into` (and `gc1`, `gc2`, `gc3`
#'   if requested), with values between 0 and 1.
#' @export
#' @examples
#' ecoli_cds |>
#'   bio_gc_content(sequence, by_position = TRUE) |>
#'   dplyr::group_by(class) |>
#'   dplyr::summarise(dplyr::across(c(gc, gc3), mean))
bio_gc_content <- function(.data, seq = NULL, by_position = FALSE, .into = "gc") {
  check_data_frame(.data)
  col <- resolve_seq_col(.data, {{ seq }})
  x <- clean_seq(check_character_col(.data, col, "seq"))
  check_new_name(.data, .into)
  gc_frac <- function(s) {
    chars <- strsplit(s, "", fixed = TRUE)[[1]]
    valid <- chars %in% c("A", "C", "G", "T", "U")
    if (!any(valid)) {
      return(NA_real_)
    }
    sum(chars %in% c("G", "C")) / sum(valid)
  }
  .data[[.into]] <- vapply(x, function(s) if (is.na(s)) NA_real_ else gc_frac(s),
                           numeric(1), USE.NAMES = FALSE)
  if (by_position) {
    for (pos in 1:3) {
      .data[[paste0("gc", pos)]] <- vapply(x, function(s) {
        if (is.na(s) || nchar(s) < pos) {
          return(NA_real_)
        }
        idx <- seq.int(pos, nchar(s), by = 3L)
        gc_frac(paste(substring(s, idx, idx), collapse = ""))
      }, numeric(1), USE.NAMES = FALSE)
    }
  }
  .data
}

count_codons <- function(x) {
  codons <- all_codons()
  counts <- integer(length(codons))
  for (s in x[!is.na(x)]) {
    s <- gsub("U", "T", s, fixed = TRUE)
    counts <- counts + tabulate(match(split_codons(s), codons), nbins = length(codons))
  }
  stats::setNames(counts, codons)
}

#' Codon usage weights from a reference set of genes
#'
#' Counts codons in a reference set of coding sequences (usually highly
#' expressed genes) and returns relative synonymous codon usage (RSCU) and
#' the relative adaptiveness `w` of each codon, as defined by Sharp and Li
#' (1987). Codons that never occur in the reference get a count of 0.5, as
#' in the original method. These weights drive [bio_cai()] and
#' [bio_codon_optimise()].
#'
#' @inheritParams bio_translate
#' @param code Genetic code used to group synonymous codons; see
#'   [bio_translate()].
#'
#' @return A tibble with one row per codon and columns `codon`,
#'   `amino_acid`, `count`, `rscu` and `w`.
#' @references Sharp, P. M. and Li, W.-H. (1987). The codon adaptation index:
#'   a measure of directional synonymous codon usage bias, and its potential
#'   applications. *Nucleic Acids Research*, 15, 1281-1295.
#'   \doi{10.1093/nar/15.3.1281}
#' @export
#' @examples
#' ecoli_cds |>
#'   dplyr::filter(class == "ribosomal") |>
#'   bio_codon_weights(sequence) |>
#'   dplyr::filter(amino_acid == "L")
bio_codon_weights <- function(.data, seq = NULL, code = "standard") {
  check_data_frame(.data)
  col <- resolve_seq_col(.data, {{ seq }})
  x <- clean_seq(check_character_col(.data, col, "seq"))
  code <- rlang::arg_match0(code, names(genetic_codes))
  weights_from_counts(count_codons(x), code)
}

weights_from_counts <- function(counts, code = "standard") {
  table <- codon_table(code)
  counts <- as.numeric(counts[names(table)])
  adj <- ifelse(counts == 0, 0.5, counts)
  family_mean <- stats::ave(adj, table, FUN = mean)
  family_max <- stats::ave(adj, table, FUN = max)
  tibble::tibble(
    codon = names(table),
    amino_acid = unname(table),
    count = counts,
    rscu = adj / family_mean,
    w = adj / family_max
  )
}

default_weights_env <- new.env(parent = emptyenv())

default_codon_weights <- function() {
  if (is.null(default_weights_env$weights)) {
    ref <- tidybiotek::ecoli_cds
    ref <- ref$sequence[ref$class == "ribosomal"]
    default_weights_env$weights <- weights_from_counts(count_codons(ref), "bacterial")
  }
  default_weights_env$weights
}

resolve_weights <- function(reference, call = rlang::caller_env()) {
  if (is.null(reference)) {
    return(default_codon_weights())
  }
  if (is.character(reference)) {
    return(weights_from_counts(count_codons(clean_seq(reference))))
  }
  if (is.data.frame(reference) && all(c("codon", "w") %in% names(reference))) {
    return(reference)
  }
  cli::cli_abort(
    c(
      "{.arg reference} must be {.code NULL}, a character vector of coding sequences,",
      " " = "or a table from {.fn bio_codon_weights}."
    ),
    call = call
  )
}

#' Codon adaptation index
#'
#' Calculates the codon adaptation index (CAI) of Sharp and Li (1987): the
#' geometric mean of the relative adaptiveness `w` of each codon in a gene.
#' Methionine (ATG), tryptophan (TGG), stop codons and codons with
#' ambiguous bases are excluded. Values near 1 indicate codon usage like
#' that of the reference set of highly expressed genes.
#'
#' @inheritParams bio_translate
#' @param reference Codon weights. `NULL` (the default) uses weights from
#'   the 56 *Escherichia coli* K-12 ribosomal protein genes in [ecoli_cds].
#'   Can also be a character vector of reference coding sequences or a table
#'   from [bio_codon_weights()].
#' @param .into Name of the new CAI column.
#'
#' @return `.data` with a numeric column `.into`.
#' @inherit bio_codon_weights references
#' @export
#' @examples
#' ecoli_cds |>
#'   bio_cai(sequence) |>
#'   dplyr::group_by(class) |>
#'   dplyr::summarise(mean_cai = mean(cai))
bio_cai <- function(.data, seq = NULL, reference = NULL, .into = "cai") {
  check_data_frame(.data)
  col <- resolve_seq_col(.data, {{ seq }})
  x <- clean_seq(check_character_col(.data, col, "seq"))
  check_new_name(.data, .into)
  weights <- resolve_weights(reference)
  w <- stats::setNames(weights$w, weights$codon)
  excluded <- c("ATG", "TGG", "TAA", "TAG", "TGA")
  if ("amino_acid" %in% names(weights)) {
    excluded <- union(excluded, weights$codon[weights$amino_acid %in% c("M", "W", "*")])
  }
  .data[[.into]] <- vapply(x, function(s) {
    if (is.na(s)) {
      return(NA_real_)
    }
    codons <- split_codons(gsub("U", "T", s, fixed = TRUE))
    codons <- codons[codons %in% names(w) & !codons %in% excluded]
    if (length(codons) == 0L) {
      return(NA_real_)
    }
    exp(mean(log(w[codons])))
  }, numeric(1), USE.NAMES = FALSE)
  .data
}

#' Codon optimisation by back-translation
#'
#' Back-translates protein sequences into coding sequences that use the
#' most adapted codon (highest `w`) for every amino acid. This is the
#' "one amino acid, one codon" strategy, which maximises the codon
#' adaptation index for the chosen host.
#'
#' @param .data A data frame, tibble or [bio_seq_tbl()].
#' @param protein <[`tidy-select`][dplyr::dplyr_tidy_select]> The protein
#'   sequence column. Can be omitted for a protein `<bio_seq_tbl>`.
#' @param reference Codon weights for the expression host; see [bio_cai()].
#' @param add_stop If `TRUE`, append the host's preferred stop codon.
#' @param .into Name of the new coding sequence column.
#'
#' @return `.data` with a character column `.into`.
#' @export
#' @examples
#' tibble::tibble(peptide = c("MKTAYIAKQR", "MSEQNNTEMTFQIQR")) |>
#'   bio_codon_optimise(peptide) |>
#'   bio_cai(optimised_cds)
bio_codon_optimise <- function(.data, protein = NULL, reference = NULL, add_stop = TRUE,
                               .into = "optimised_cds") {
  check_data_frame(.data)
  col <- resolve_seq_col(.data, {{ protein }})
  x <- clean_seq(check_character_col(.data, col, "protein"))
  check_new_name(.data, .into)
  weights <- resolve_weights(reference)
  if (!"amino_acid" %in% names(weights)) {
    weights$amino_acid <- unname(codon_table()[weights$codon])
  }
  ord <- order(weights$amino_acid, -weights$w, -weights$count %||% 0)
  best <- weights[ord, ]
  best <- best[!duplicated(best$amino_acid), ]
  best <- stats::setNames(best$codon, best$amino_acid)
  .data[[.into]] <- vapply(x, function(s) {
    if (is.na(s)) {
      return(NA_character_)
    }
    aa <- strsplit(sub("\\*$", "", s), "", fixed = TRUE)[[1]]
    unknown <- setdiff(aa, names(best))
    if (length(unknown) > 0L) {
      cli::cli_abort("Cannot back-translate residue{?s} {.val {unknown}}.")
    }
    out <- paste(best[aa], collapse = "")
    if (add_stop) paste0(out, best[["*"]]) else out
  }, character(1), USE.NAMES = FALSE)
  .data
}

# SantaLucia (1998) unified nearest-neighbour parameters.
# dH in kcal/mol and dS in cal/(K mol).
nn_params <- data.frame(
  pair = c("AA", "AT", "TA", "CA", "GT", "CT", "GA", "CG", "GC", "GG"),
  dh = c(-7.9, -7.2, -7.2, -8.5, -8.4, -7.8, -8.2, -10.6, -9.8, -8.0),
  ds = c(-22.2, -20.4, -21.3, -22.7, -22.4, -21.0, -22.2, -27.2, -24.4, -19.9)
)

nn_lookup <- function() {
  rc <- reverse_complement(nn_params$pair)
  pairs <- c(nn_params$pair, rc)
  keep <- !duplicated(pairs)
  list(
    dh = stats::setNames(c(nn_params$dh, nn_params$dh)[keep], pairs[keep]),
    ds = stats::setNames(c(nn_params$ds, nn_params$ds)[keep], pairs[keep])
  )
}

tm_nearest_neighbour <- function(s, na_mm, primer_nm) {
  n <- nchar(s)
  lookup <- nn_lookup()
  pairs <- substring(s, seq_len(n - 1L), 2:n)
  dh <- sum(lookup$dh[pairs])
  ds <- sum(lookup$ds[pairs])
  for (end in c(substr(s, 1, 1), substr(s, n, n))) {
    if (end %in% c("G", "C")) {
      dh <- dh + 0.1
      ds <- ds - 2.8
    } else {
      dh <- dh + 2.3
      ds <- ds + 4.1
    }
  }
  self_comp <- identical(s, reverse_complement(s))
  if (self_comp) {
    ds <- ds - 1.4
  }
  ds <- ds + 0.368 * (n - 1L) * log(na_mm / 1000)
  ct <- primer_nm * 1e-9
  x <- if (self_comp) 1 else 4
  (dh * 1000) / (ds + 1.987 * log(ct / x)) - 273.15
}

#' Primer melting temperature
#'
#' Estimates the melting temperature (Tm, in degrees Celsius) of DNA
#' oligonucleotides. The default nearest-neighbour method uses the unified
#' parameters of SantaLucia (1998) with the entropy salt correction for
#' monovalent cations. The Wallace rule (2 degrees per A or T, 4 per G or C)
#' and the basic GC formula are offered for quick checks.
#'
#' @inheritParams bio_translate
#' @param method One of `"nearest_neighbour"`, `"wallace"` or `"basic"`.
#' @param na_conc Monovalent cation (Na+) concentration in mM.
#' @param primer_conc Total oligonucleotide concentration in nM.
#' @param .into Name of the new Tm column.
#'
#' @return `.data` with a numeric column `.into`. Sequences with bases
#'   other than A, C, G and T give `NA` with a warning.
#' @references SantaLucia, J. (1998). A unified view of polymer, dumbbell,
#'   and oligonucleotide DNA nearest-neighbor thermodynamics. *Proceedings
#'   of the National Academy of Sciences*, 95, 1460-1465.
#'   \doi{10.1073/pnas.95.4.1460}
#'
#'   Wallace, R. B. et al. (1979). Hybridization of synthetic
#'   oligodeoxyribonucleotides to phi chi 174 DNA: the effect of single
#'   base pair mismatch. *Nucleic Acids Research*, 6, 3543-3557.
#'   \doi{10.1093/nar/6.11.3543}
#' @export
#' @examples
#' primers <- tibble::tibble(
#'   name = c("M13F", "M13R", "T7"),
#'   primer = c("GTAAAACGACGGCCAGT", "CAGGAAACAGCTATGAC", "TAATACGACTCACTATAGGG")
#' )
#' primers |>
#'   bio_primer_tm(primer) |>
#'   bio_primer_tm(primer, method = "wallace", .into = "tm_wallace")
bio_primer_tm <- function(.data, seq = NULL, method = c("nearest_neighbour", "wallace", "basic"),
                          na_conc = 50, primer_conc = 50, .into = "tm") {
  check_data_frame(.data)
  col <- resolve_seq_col(.data, {{ seq }})
  x <- clean_seq(check_character_col(.data, col, "seq"))
  method <- rlang::arg_match(method)
  check_positive_number(na_conc)
  check_positive_number(primer_conc)
  check_new_name(.data, .into)
  x <- gsub("U", "T", x, fixed = TRUE)
  bad <- !is.na(x) & !grepl("^[ACGT]{2,}$", x)
  if (any(bad)) {
    cli::cli_warn(
      "{sum(bad)} sequence{?s} with ambiguous bases or fewer than 2 bases give NA."
    )
  }
  .data[[.into]] <- vapply(seq_along(x), function(i) {
    s <- x[[i]]
    if (is.na(s) || bad[[i]]) {
      return(NA_real_)
    }
    n <- nchar(s)
    gc <- lengths(regmatches(s, gregexpr("[GC]", s)))
    switch(method,
      nearest_neighbour = tm_nearest_neighbour(s, na_conc, primer_conc),
      wallace = 2 * (n - gc) + 4 * gc,
      basic = 64.9 + 41 * (gc - 16.4) / n
    )
  }, numeric(1))
  .data
}

kmer_vocabulary <- function(x, k, alphabet) {
  if (alphabet %in% c("dna", "rna") && k <= 6L) {
    bases <- if (alphabet == "rna") c("A", "C", "G", "U") else c("A", "C", "G", "T")
    grid <- expand.grid(rep(list(bases), k), stringsAsFactors = FALSE)
    return(sort(do.call(paste0, grid)))
  }
  observed <- unique(unlist(lapply(x[!is.na(x)], kmers_of, k = k)))
  sort(observed)
}

kmers_of <- function(s, k) {
  n <- nchar(s)
  if (n < k) {
    return(character())
  }
  substring(s, seq_len(n - k + 1L), k:n)
}

kmer_matrix <- function(x, k, vocabulary, normalise = FALSE) {
  m <- vapply(x, function(s) {
    if (is.na(s)) {
      return(rep(NA_real_, length(vocabulary)))
    }
    counts <- tabulate(match(kmers_of(s, k), vocabulary), nbins = length(vocabulary))
    if (normalise) {
      total <- max(nchar(s) - k + 1L, 0L)
      if (total > 0) counts / total else counts * 0
    } else {
      as.numeric(counts)
    }
  }, numeric(length(vocabulary)), USE.NAMES = FALSE)
  m <- matrix(m, ncol = length(x))
  t(m)
}

#' Count k-mers in sequences
#'
#' Counts overlapping k-mers (substrings of length `k`) in each sequence and
#' returns them as one column per k-mer, ready for modelling. For DNA and
#' RNA with `k` up to 6 every possible k-mer gets a column; otherwise the
#' k-mers observed in the data are used.
#'
#' @inheritParams bio_translate
#' @param k K-mer length.
#' @param normalise If `TRUE`, return frequencies (counts divided by the
#'   number of k-mers in the sequence) instead of counts.
#' @param vocabulary Optional character vector of k-mers to count.
#' @param prefix Prefix for the new column names.
#'
#' @return `.data` with one numeric column per k-mer, named
#'   `paste0(prefix, kmer)`.
#' @export
#' @examples
#' tibble::tibble(id = 1:2, seq = c("ACGTACGT", "GGGCCCAA")) |>
#'   bio_kmer_count(seq, k = 2)
bio_kmer_count <- function(.data, seq = NULL, k = 3, normalise = FALSE, vocabulary = NULL,
                           prefix = "kmer_") {
  check_data_frame(.data)
  col <- resolve_seq_col(.data, {{ seq }})
  x <- clean_seq(check_character_col(.data, col, "seq"))
  k <- check_whole_number(k)
  alphabet <- resolve_alphabet(.data, col, x)
  vocabulary <- vocabulary %||% kmer_vocabulary(x, k, alphabet)
  if (length(vocabulary) == 0L) {
    cli::cli_abort("No k-mers of length {k} were found; use a smaller {.arg k}.")
  }
  if (length(vocabulary) > 5000L) {
    cli::cli_warn("Creating {length(vocabulary)} k-mer columns; consider a smaller {.arg k}.")
  }
  m <- kmer_matrix(x, k, vocabulary, normalise)
  colnames(m) <- paste0(prefix, vocabulary)
  clash <- intersect(colnames(m), names(.data))
  if (length(clash) > 0L) {
    cli::cli_abort("Column{?s} {.field {clash}} already exist{?s/}; choose another {.arg prefix}.")
  }
  for (j in seq_len(ncol(m))) {
    .data[[colnames(m)[[j]]]] <- m[, j]
  }
  .data
}

# Amino acid property tables.
kyte_doolittle <- c(
  A = 1.8, R = -4.5, N = -3.5, D = -3.5, C = 2.5, Q = -3.5, E = -3.5, G = -0.4,
  H = -3.2, I = 4.5, L = 3.8, K = -3.9, M = 1.9, F = 2.8, P = -1.6, S = -0.8,
  T = -0.7, W = -0.9, Y = -1.3, V = 4.2
)

residue_mass <- c(
  A = 71.0788, R = 156.1875, N = 114.1038, D = 115.0886, C = 103.1388,
  E = 129.1155, Q = 128.1307, G = 57.0519, H = 137.1411, I = 113.1594,
  L = 113.1594, K = 128.1741, M = 131.1926, F = 147.1766, P = 97.1167,
  S = 87.0782, T = 101.1051, W = 186.2132, Y = 163.1760, V = 99.1326
)

# EMBOSS pKa values.
pka_values <- list(
  n_term = 8.6, c_term = 3.6,
  positive = c(K = 10.8, R = 12.5, H = 6.5),
  negative = c(D = 3.9, E = 4.1, C = 8.5, Y = 10.1)
)

net_charge <- function(counts, ph) {
  pos <- 1 / (1 + 10^(ph - pka_values$n_term)) +
    sum(counts[names(pka_values$positive)] / (1 + 10^(ph - pka_values$positive)))
  neg <- 1 / (1 + 10^(pka_values$c_term - ph)) +
    sum(counts[names(pka_values$negative)] / (1 + 10^(pka_values$negative - ph)))
  pos - neg
}

aa_props_one <- function(s) {
  aa <- strsplit(sub("\\*$", "", s), "", fixed = TRUE)[[1]]
  aa <- aa[aa %in% names(kyte_doolittle)]
  n <- length(aa)
  if (n == 0L) {
    return(rep(NA_real_, 9))
  }
  counts <- stats::setNames(tabulate(match(aa, names(kyte_doolittle)), 20),
                            names(kyte_doolittle))
  iep <- stats::uniroot(function(ph) net_charge(counts, ph), c(0, 14), tol = 1e-6)$root
  c(
    length = n,
    mol_weight = sum(residue_mass[aa]) + 18.01524,
    isoelectric_point = iep,
    charge_ph7 = net_charge(counts, 7),
    gravy = mean(kyte_doolittle[aa]),
    aromaticity = sum(counts[c("F", "W", "Y")]) / n,
    frac_acidic = sum(counts[c("D", "E")]) / n,
    frac_basic = sum(counts[c("K", "R", "H")]) / n,
    frac_hydrophobic = sum(counts[c("A", "V", "I", "L", "M", "F", "W")]) / n
  )
}

aa_prop_names <- c(
  "length", "mol_weight", "isoelectric_point", "charge_ph7", "gravy", "aromaticity",
  "frac_acidic", "frac_basic", "frac_hydrophobic"
)

aa_props_matrix <- function(x) {
  m <- vapply(clean_seq(x), function(s) {
    if (is.na(s)) rep(NA_real_, 9) else aa_props_one(s)
  }, numeric(9), USE.NAMES = FALSE)
  m <- t(matrix(m, nrow = 9))
  colnames(m) <- aa_prop_names
  m
}

#' Physicochemical properties of proteins
#'
#' Calculates sequence-based protein properties used for protein
#' engineering and expression screening: length, average molecular weight
#' (Da), isoelectric point and net charge at pH 7 (EMBOSS pKa values), GRAVY
#' hydropathy (Kyte and Doolittle, 1982), aromaticity and the fractions of
#' acidic, basic and hydrophobic residues. Non-standard residues are
#' ignored and a trailing stop (`*`) is dropped.
#'
#' @inheritParams bio_translate
#' @param prefix Prefix for the new column names.
#'
#' @return `.data` with nine numeric columns: `length`, `mol_weight`,
#'   `isoelectric_point`, `charge_ph7`, `gravy`, `aromaticity`,
#'   `frac_acidic`, `frac_basic` and `frac_hydrophobic`, each with `prefix`.
#' @references Kyte, J. and Doolittle, R. F. (1982). A simple method for
#'   displaying the hydropathic character of a protein. *Journal of
#'   Molecular Biology*, 157, 105-132. \doi{10.1016/0022-2836(82)90515-0}
#' @export
#' @examples
#' ecoli_cds |>
#'   dplyr::slice(1:5) |>
#'   bio_aa_props(protein) |>
#'   dplyr::select(gene, length, mol_weight, isoelectric_point, gravy)
bio_aa_props <- function(.data, seq = NULL, prefix = "") {
  check_data_frame(.data)
  col <- resolve_seq_col(.data, {{ seq }})
  x <- check_character_col(.data, col, "seq")
  m <- aa_props_matrix(x)
  for (j in seq_len(ncol(m))) {
    .data[[paste0(prefix, colnames(m)[[j]])]] <- m[, j]
  }
  .data
}
