# Over-representation analysis of gene sets

Tests whether a list of selected features (for example significant
genes) contains more members of each gene set or pathway than expected
by chance, using the one-sided hypergeometric test (Fisher's exact test
for over-representation). Only features in the universe are counted.
P-values are adjusted with the Benjamini-Hochberg method.

## Usage

``` r
bio_enrich(
  .data,
  feature,
  selected,
  gene_sets,
  universe = NULL,
  min_size = 5,
  max_size = 500
)
```

## Arguments

- .data:

  A data frame with one row per feature, such as the output of
  [`bio_diff_expr()`](https://tidybiotek.in4metrix.dev/reference/bio_diff_expr.md).
  All of its features form the universe unless `universe` is given.

- feature:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  The feature identifier column.

- selected:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>
  A logical expression choosing the features of interest, for example
  `p_adj < 0.05 & log2fc > 1`.

- gene_sets:

  Either a named list of character vectors or a data frame with columns
  `gene_set` and `feature`.

- universe:

  Optional character vector of background features.

- min_size, max_size:

  Gene sets with fewer or more members in the universe are skipped.

## Value

A `<bio_omics>` sorted by p-value with columns `gene_set`, `set_size`,
`overlap`, `expected`, `fold_enrichment`, `p_value`, `p_adj` and
`overlap_features` (a list column).

## References

Boyle, E. I. et al. (2004). GO::TermFinder: open source software for
accessing Gene Ontology information and finding significantly enriched
Gene Ontology terms associated with a list of genes. *Bioinformatics*,
20, 3710-3715.
[doi:10.1093/bioinformatics/bth456](https://doi.org/10.1093/bioinformatics/bth456)

## Examples

``` r
de <- bio_diff_expr(expression_sim, gene, sample, count, condition,
  treatment = "stress", reference = "control", transform = "log2_cpm"
)
bio_enrich(de, feature, p_adj < 0.05 & log2fc > 0, gene_sets_sim)
#> # A bio_omics: enrichment: 30 x 8
#>    gene_set          set_size overlap expected fold_enrichment  p_value    p_adj
#>    <chr>                <int>   <int>    <dbl>           <dbl>    <dbl>    <dbl>
#>  1 stress_response         40      26    1.72            15.1  4.31e-30 1.29e-28
#>  2 background_set_18       52       6    2.24             2.68 2.03e- 2 3.05e- 1
#>  3 background_set_21       17       3    0.731            4.10 3.31e- 2 3.31e- 1
#>  4 background_set_12       58       5    2.49             2.00 9.74e- 2 7.31e- 1
#>  5 background_set_20       18       2    0.774            2.58 1.79e- 1 1   e+ 0
#>  6 background_set_28       21       2    0.903            2.21 2.27e- 1 1   e+ 0
#>  7 background_set_06       43       3    1.85             1.62 2.80e- 1 1   e+ 0
#>  8 background_set_25       43       3    1.85             1.62 2.80e- 1 1   e+ 0
#>  9 background_set_10       28       2    1.20             1.66 3.41e- 1 1   e+ 0
#> 10 background_set_11       30       2    1.29             1.55 3.73e- 1 1   e+ 0
#> # ℹ 20 more rows
#> # ℹ 1 more variable: overlap_features <list>
```
