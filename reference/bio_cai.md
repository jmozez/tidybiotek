# Codon adaptation index

Calculates the codon adaptation index (CAI) of Sharp and Li (1987): the
geometric mean of the relative adaptiveness `w` of each codon in a gene.
Methionine (ATG), tryptophan (TGG), stop codons and codons with
ambiguous bases are excluded. Values near 1 indicate codon usage like
that of the reference set of highly expressed genes.

## Usage

``` r
bio_cai(.data, seq = NULL, reference = NULL, .into = "cai")
```

## Arguments

- .data:

  A data frame, tibble or
  [`bio_seq_tbl()`](https://tidybiotek.in4metrix.dev/reference/bio_seq_tbl.md).

- seq:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  The sequence column. Can be omitted for a `<bio_seq_tbl>`.

- reference:

  Codon weights. `NULL` (the default) uses weights from the 56
  *Escherichia coli* K-12 ribosomal protein genes in
  [ecoli_cds](https://tidybiotek.in4metrix.dev/reference/ecoli_cds.md).
  Can also be a character vector of reference coding sequences or a
  table from
  [`bio_codon_weights()`](https://tidybiotek.in4metrix.dev/reference/bio_codon_weights.md).

- .into:

  Name of the new CAI column.

## Value

`.data` with a numeric column `.into`.

## References

Sharp, P. M. and Li, W.-H. (1987). The codon adaptation index: a measure
of directional synonymous codon usage bias, and its potential
applications. *Nucleic Acids Research*, 15, 1281-1295.
[doi:10.1093/nar/15.3.1281](https://doi.org/10.1093/nar/15.3.1281)

## Examples

``` r
ecoli_cds |>
  bio_cai(sequence) |>
  dplyr::group_by(class) |>
  dplyr::summarise(mean_cai = mean(cai))
#> # A tibble: 2 × 2
#>   class     mean_cai
#>   <chr>        <dbl>
#> 1 other        0.439
#> 2 ribosomal    0.695
```
