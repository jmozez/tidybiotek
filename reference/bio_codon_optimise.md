# Codon optimisation by back-translation

Back-translates protein sequences into coding sequences that use the
most adapted codon (highest `w`) for every amino acid. This is the "one
amino acid, one codon" strategy, which maximises the codon adaptation
index for the chosen host.

## Usage

``` r
bio_codon_optimise(
  .data,
  protein = NULL,
  reference = NULL,
  add_stop = TRUE,
  .into = "optimised_cds"
)
```

## Arguments

- .data:

  A data frame, tibble or
  [`bio_seq_tbl()`](https://tidybiotek.in4metrix.dev/reference/bio_seq_tbl.md).

- protein:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  The protein sequence column. Can be omitted for a protein
  `<bio_seq_tbl>`.

- reference:

  Codon weights for the expression host; see
  [`bio_cai()`](https://tidybiotek.in4metrix.dev/reference/bio_cai.md).

- add_stop:

  If `TRUE`, append the host's preferred stop codon.

- .into:

  Name of the new coding sequence column.

## Value

`.data` with a character column `.into`.

## Examples

``` r
tibble::tibble(peptide = c("MKTAYIAKQR", "MSEQNNTEMTFQIQR")) |>
  bio_codon_optimise(peptide) |>
  bio_cai(optimised_cds)
#> # A tibble: 2 × 3
#>   peptide         optimised_cds                                      cai
#>   <chr>           <chr>                                            <dbl>
#> 1 MKTAYIAKQR      ATGAAAACTGCTTACATCGCTAAACAGCGTTAA                    1
#> 2 MSEQNNTEMTFQIQR ATGTCTGAACAGAACAACACTGAAATGACTTTCCAGATCCAGCGTTAA     1
```
