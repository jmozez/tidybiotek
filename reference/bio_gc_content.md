# GC content of nucleotide sequences

Calculates the fraction of G and C among unambiguous bases (A, C, G and
T or U). With `by_position = TRUE` it also returns GC at the first,
second and third codon positions (GC1, GC2 and GC3), which are widely
used to compare codon usage between genes and genomes.

## Usage

``` r
bio_gc_content(.data, seq = NULL, by_position = FALSE, .into = "gc")
```

## Arguments

- .data:

  A data frame, tibble or
  [`bio_seq_tbl()`](https://tidybiotek.in4metrix.dev/reference/bio_seq_tbl.md).

- seq:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  The sequence column. Can be omitted for a `<bio_seq_tbl>`.

- by_position:

  If `TRUE`, add `gc1`, `gc2` and `gc3` columns as well.

- .into:

  Name of the overall GC column.

## Value

`.data` with a numeric column `.into` (and `gc1`, `gc2`, `gc3` if
requested), with values between 0 and 1.

## Examples

``` r
ecoli_cds |>
  bio_gc_content(sequence, by_position = TRUE) |>
  dplyr::group_by(class) |>
  dplyr::summarise(dplyr::across(c(gc, gc3), mean))
#> # A tibble: 2 × 3
#>   class        gc   gc3
#>   <chr>     <dbl> <dbl>
#> 1 other     0.520 0.569
#> 2 ribosomal 0.505 0.490
```
