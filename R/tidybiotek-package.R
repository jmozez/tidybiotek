#' tidybiotek: Tidy Tools for Applied Biotechnology Bioinformatics
#'
#' tidybiotek brings applied biotechnology bioinformatics into the
#' tidyverse and tidymodels. Every user-facing function takes a data frame
#' first, uses tidy evaluation for column arguments and returns a tibble or
#' a tibble subclass, so analyses read as pipelines.
#'
#' The package is organised into modules:
#'
#' * **Sequences**: [bio_seq_tbl()], [bio_translate()], [bio_gc_content()],
#'   [bio_codon_weights()], [bio_cai()], [bio_codon_optimise()],
#'   [bio_primer_tm()], [bio_kmer_count()] and [bio_aa_props()].
#' * **Assays**: [bio_qpcr_ddct()], [bio_dose_response()] and
#'   [bio_elisa_quant()].
#' * **Bioprocess**: [bio_enzyme_kinetics()] and [bio_monod_growth()].
#' * **Omics**: [bio_diff_expr()] and [bio_enrich()].
#' * **tidymodels**: [step_bio_kmers()], [step_bio_aa_props()],
#'   [step_bio_vst_norm()], [bio_homology_cv()], [bio_batch_cv()],
#'   [enrichment_factor()], [bedroc()], [assay_cv_error()] and
#'   [bio_four_pl()].
#' * **Reporting**: `tidy()`, `glance()`, `augment()` and `autoplot()`
#'   methods.
#'
#' All numerical methods (codon adaptation, nearest-neighbour
#' thermodynamics, Levenberg-Marquardt four-parameter logistic fitting,
#' the exact integrated solution of Monod kinetics and empirical Bayes variance
#' moderation) are written in base R and `stats`.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang .data .env := %||%
#' @importFrom generics tidy glance augment required_pkgs
#' @importFrom recipes prep bake
#' @importFrom ggplot2 autoplot

#' @importFrom stats predict
#' @importFrom dplyr dplyr_reconstruct
#' @importFrom pillar tbl_sum
#' @importFrom utils head
## usethis namespace: end
NULL

#' @export
generics::tidy

#' @export
generics::glance

#' @export
generics::augment

#' @export
generics::required_pkgs

#' @export
ggplot2::autoplot

utils::globalVariables(c(
  ".fitted", ".resid", "estimate", "std.error", "statistic", "object", "new_data"
))

.onLoad <- function(libname, pkgname) {
  make_bio_four_pl()
  invisible()
}
