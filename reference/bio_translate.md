# Translate nucleotide sequences to protein

Translates DNA or RNA sequences with an NCBI genetic code. Codons that
contain ambiguous bases become `X` and stop codons become `*`.

## Usage

``` r
bio_translate(
  .data,
  seq = NULL,
  frame = 1,
  code = "standard",
  to_stop = FALSE,
  start_as_met = FALSE,
  .into = "translated"
)
```

## Arguments

- .data:

  A data frame, tibble or
  [`bio_seq_tbl()`](https://tidybiotek.in4metrix.dev/reference/bio_seq_tbl.md).

- seq:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  The sequence column. Can be omitted for a `<bio_seq_tbl>`.

- frame:

  Reading frame: `1`, `2` or `3` on the given strand, or `-1`, `-2` or
  `-3` on the reverse complement.

- code:

  Genetic code: `"standard"` (NCBI table 1), `"bacterial"` (table 11),
  `"vertebrate_mito"` (table 2) or `"yeast_mito"` (table 3).

- to_stop:

  If `TRUE`, stop translation at the first stop codon and drop it.

- start_as_met:

  If `TRUE`, translate an alternative start codon in the first position
  as methionine, as happens in vivo.

- .into:

  Name of the new protein column.

## Value

`.data` with a new character column `.into`. A `<bio_seq_tbl>` keeps its
class.

## References

NCBI. The Genetic Codes.
<https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>

## Examples

``` r
ecoli_cds |>
  bio_translate(sequence, code = "bacterial", to_stop = TRUE, start_as_met = TRUE) |>
  dplyr::mutate(matches = translated == protein) |>
  dplyr::count(matches)
#> # A tibble: 1 × 2
#>   matches     n
#>   <lgl>   <int>
#> 1 TRUE       76
```
