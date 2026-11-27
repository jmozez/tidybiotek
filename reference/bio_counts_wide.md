# Reshape long counts to one row per sample

A small helper that turns long count data, such as
[expression_sim](https://tidybiotek.in4metrix.dev/reference/expression_sim.md),
into the wide layout that recipes and models expect: one row per sample
and one column per gene.

## Usage

``` r
bio_counts_wide(
  .data,
  feature = "gene",
  sample = "sample",
  value = "count",
  keep = "condition"
)
```

## Arguments

- .data:

  A data frame in long format.

- feature, sample, value:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Feature, sample and value columns.

- keep:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Sample-level columns to carry across, such as the condition.

## Value

A tibble with one row per sample.

## Examples

``` r
wide <- bio_counts_wide(expression_sim)
dim(wide)
#> [1]    8 1002
```
