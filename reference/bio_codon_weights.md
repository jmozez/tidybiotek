# Codon usage weights from a reference set of genes

Counts codons in a reference set of coding sequences (usually highly
expressed genes) and returns relative synonymous codon usage (RSCU) and
the relative adaptiveness `w` of each codon, as defined by Sharp and Li
(1987). Codons that never occur in the reference get a count of 0.5, as
in the original method. These weights drive
[`bio_cai()`](https://tidybiotek.in4metrix.dev/reference/bio_cai.md) and
[`bio_codon_optimise()`](https://tidybiotek.in4metrix.dev/reference/bio_codon_optimise.md).

## Usage

``` r
bio_codon_weights(.data, seq = NULL, code = "standard")
```

## Arguments

- .data:

  A data frame, tibble or
  [`bio_seq_tbl()`](https://tidybiotek.in4metrix.dev/reference/bio_seq_tbl.md).

- seq:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  The sequence column. Can be omitted for a `<bio_seq_tbl>`.

- code:

  Genetic code used to group synonymous codons; see
  [`bio_translate()`](https://tidybiotek.in4metrix.dev/reference/bio_translate.md).

## Value

A tibble with one row per codon and columns `codon`, `amino_acid`,
`count`, `rscu` and `w`.

## References

Sharp, P. M. and Li, W.-H. (1987). The codon adaptation index: a measure
of directional synonymous codon usage bias, and its potential
applications. *Nucleic Acids Research*, 15, 1281-1295.
[doi:10.1093/nar/15.3.1281](https://doi.org/10.1093/nar/15.3.1281)

## Examples

``` r
ecoli_cds |>
  dplyr::filter(class == "ribosomal") |>
  bio_codon_weights(sequence) |>
  dplyr::filter(amino_acid == "L")
#> # A tibble: 6 × 5
#>   codon amino_acid count   rscu       w
#>   <chr> <chr>      <dbl>  <dbl>   <dbl>
#> 1 TTA   L             15 0.179  0.0356 
#> 2 TTG   L             21 0.25   0.0499 
#> 3 CTT   L             25 0.298  0.0594 
#> 4 CTC   L             20 0.238  0.0475 
#> 5 CTA   L              2 0.0238 0.00475
#> 6 CTG   L            421 5.01   1      
```
