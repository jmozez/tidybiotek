# Count k-mers in sequences

Counts overlapping k-mers (substrings of length `k`) in each sequence
and returns them as one column per k-mer, ready for modelling. For DNA
and RNA with `k` up to 6 every possible k-mer gets a column; otherwise
the k-mers observed in the data are used.

## Usage

``` r
bio_kmer_count(
  .data,
  seq = NULL,
  k = 3,
  normalise = FALSE,
  vocabulary = NULL,
  prefix = "kmer_"
)
```

## Arguments

- .data:

  A data frame, tibble or
  [`bio_seq_tbl()`](https://tidybiotek.in4metrix.dev/reference/bio_seq_tbl.md).

- seq:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  The sequence column. Can be omitted for a `<bio_seq_tbl>`.

- k:

  K-mer length.

- normalise:

  If `TRUE`, return frequencies (counts divided by the number of k-mers
  in the sequence) instead of counts.

- vocabulary:

  Optional character vector of k-mers to count.

- prefix:

  Prefix for the new column names.

## Value

`.data` with one numeric column per k-mer, named `paste0(prefix, kmer)`.

## Examples

``` r
tibble::tibble(id = 1:2, seq = c("ACGTACGT", "GGGCCCAA")) |>
  bio_kmer_count(seq, k = 2)
#> # A tibble: 2 × 18
#>      id seq      kmer_AA kmer_AC kmer_AG kmer_AT kmer_CA kmer_CC kmer_CG kmer_CT
#>   <int> <chr>      <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#> 1     1 ACGTACGT       0       2       0       0       0       0       2       0
#> 2     2 GGGCCCAA       1       0       0       0       1       2       0       0
#> # ℹ 8 more variables: kmer_GA <dbl>, kmer_GC <dbl>, kmer_GG <dbl>,
#> #   kmer_GT <dbl>, kmer_TA <dbl>, kmer_TC <dbl>, kmer_TG <dbl>, kmer_TT <dbl>
```
