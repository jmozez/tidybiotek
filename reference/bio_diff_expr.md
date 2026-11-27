# Differential expression with moderated t-statistics

Compares two groups feature by feature (genes, transcripts, proteins or
metabolites) using the empirical Bayes moderated t-test of Smyth (2004),
the method behind limma, written in base R. A one-way linear model is
fitted to every feature, the residual variances are shrunk towards a
common prior estimated from all features, and moderated t-statistics
with extra degrees of freedom are formed for the chosen contrast.
P-values are adjusted with the Benjamini-Hochberg method.

## Usage

``` r
bio_diff_expr(
  .data,
  feature,
  sample,
  value,
  group,
  treatment = NULL,
  reference = NULL,
  transform = c("none", "log2", "log2_cpm"),
  robust_prior = TRUE
)
```

## Arguments

- .data:

  A data frame in long format.

- feature, sample, value, group:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Columns identifying the feature, the sample, the measured value and
  the experimental group of each sample.

- treatment, reference:

  Levels of `group` to compare (log fold change is
  `treatment - reference`). By default the first level is the reference
  and the second the treatment.

- transform:

  `"none"` (values are already on a log scale), `"log2"` (adds 1 then
  takes log2) or `"log2_cpm"` (counts to log2 counts per million with a
  prior count of 0.5).

- robust_prior:

  If `FALSE`, use ordinary t-tests (no variance moderation) instead.

## Value

A `<bio_omics>` sorted by p-value, with columns `feature`, `log2fc`,
`ave_expr`, `t`, `p_value`, `p_adj` and `std_error`. The prior degrees
of freedom and variance are kept as metadata; see
[glance()](https://tidybiotek.in4metrix.dev/reference/broom_bio.md).

## Details

The input is tidy (long) data with one row per feature and sample.
Counts can be transformed to log2 counts per million first.

## References

Smyth, G. K. (2004). Linear models and empirical Bayes methods for
assessing differential expression in microarray experiments.
*Statistical Applications in Genetics and Molecular Biology*, 3, 3.
[doi:10.2202/1544-6115.1027](https://doi.org/10.2202/1544-6115.1027)

Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery
rate: a practical and powerful approach to multiple testing. *Journal of
the Royal Statistical Society B*, 57, 289-300.
[doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

## Examples

``` r
de <- bio_diff_expr(expression_sim, gene, sample, count, condition,
  treatment = "stress", reference = "control", transform = "log2_cpm"
)
de
#> # A bio_omics: diff_expr: 1000 x 7
#>    feature   log2fc ave_expr      t  p_value       p_adj std_error
#>    <chr>      <dbl>    <dbl>  <dbl>    <dbl>       <dbl>     <dbl>
#>  1 gene_0054  -3.10     6.35 -12.1  4.88e-10 0.000000245     0.256
#>  2 gene_0053  -2.75     7.42 -12.1  4.91e-10 0.000000245     0.227
#>  3 gene_0071  -2.87     7.72 -10.9  2.60e- 9 0.000000539     0.263
#>  4 gene_0051  -3.12     7.60 -10.8  2.80e- 9 0.000000539     0.288
#>  5 gene_0067  -2.70    10.9  -10.8  2.95e- 9 0.000000539     0.250
#>  6 gene_0047   2.94     9.69  10.7  3.23e- 9 0.000000539     0.273
#>  7 gene_0042   2.84    11.3   10.5  4.38e- 9 0.000000572     0.270
#>  8 gene_0062  -2.71     7.71 -10.5  4.94e- 9 0.000000572     0.260
#>  9 gene_0056  -2.91     6.20 -10.4  5.14e- 9 0.000000572     0.279
#> 10 gene_0008   2.70     8.68   9.70 1.55e- 8 0.00000155      0.279
#> # ℹ 990 more rows
glance(de)
#> # A tibble: 1 × 7
#>   n_features n_significant df_prior s0_sq df_residual treatment reference
#>        <int>         <int>    <dbl> <dbl>       <int> <chr>     <chr>    
#> 1       1000            69     11.8 0.150           6 stress    control  
```
