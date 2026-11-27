# Tidy, glance and augment tidybiotek results

Summarise assay, kinetic and omics results in the broom style.

## Usage

``` r
# S3 method for class 'bio_assay'
tidy(x, conf.level = NULL, ...)

# S3 method for class 'bio_assay'
glance(x, ...)

# S3 method for class 'bio_assay'
augment(x, ...)

# S3 method for class 'bio_omics'
tidy(x, ...)

# S3 method for class 'bio_omics'
glance(x, alpha = 0.05, ...)

# S3 method for class 'bio_omics'
augment(x, alpha = 0.05, lfc = 0, ...)

# S3 method for class 'bio_four_pl_fit'
tidy(x, ...)

# S3 method for class 'bio_four_pl_fit'
glance(x, ...)
```

## Arguments

- x:

  A `<bio_assay>`, `<bio_omics>` or `bio_four_pl_fit` object.

- conf.level:

  Confidence level for intervals.

- ...:

  Not used.

- alpha:

  Adjusted p-value threshold used to flag significance.

- lfc:

  Minimum absolute log2 fold change for significance.

## Value

A tibble.

## Details

- [`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns one
  row per parameter (with standard errors, statistics, p-values and
  confidence intervals) for dose-response, ELISA, enzyme kinetic and
  Monod fits; per-group fold changes for qPCR; and per-feature or
  per-gene-set results for omics analyses.

- [`glance()`](https://generics.r-lib.org/reference/glance.html) returns
  one row per fitted curve or analysis with model statistics such as the
  residual standard deviation, log-likelihood and AIC, or the empirical
  Bayes prior for differential expression.

- [`augment()`](https://generics.r-lib.org/reference/augment.html)
  returns the original data with fitted values (`.fitted`) and residuals
  (`.resid`), or significance flags for omics results.

## Examples

``` r
fit <- bio_dose_response(dose_response_sim, dose, response, .by = compound)
tidy(fit)
#> # A tibble: 12 × 8
#>    compound    term   estimate std.error statistic    p.value conf.low conf.high
#>    <chr>       <chr>     <dbl>     <dbl>     <dbl>      <dbl>    <dbl>     <dbl>
#>  1 activator_c bottom   7.69     0.874        8.79  1.13 e- 9   5.90      9.47  
#>  2 activator_c top     85.9      1.21        71.2   4.05 e-34  83.5      88.4   
#>  3 activator_c ec50     1.20     0.0743      NA    NA           1.06      1.36  
#>  4 activator_c hill     1.38     0.104       13.2   8.70 e-14   1.16      1.59  
#>  5 inhibitor_a bottom   4.02     1.16         3.45  1.73 e- 3   1.64      6.40  
#>  6 inhibitor_a top    102.       2.15        47.4   5.05 e-29  97.6     106.    
#>  7 inhibitor_a ec50     0.0490   0.00408     NA    NA           0.0413    0.0581
#>  8 inhibitor_a hill    -1.16     0.0930     -12.5   3.53 e-13  -1.35     -0.970 
#>  9 inhibitor_b bottom   9.31     1.56         5.96  1.78 e- 6   6.11     12.5   
#> 10 inhibitor_b top     98.6      0.958      103.    1.000e-38  96.6     101.    
#> 11 inhibitor_b ec50     2.34     0.142       NA    NA           2.07      2.65  
#> 12 inhibitor_b hill    -1.61     0.133      -12.1   7.89 e-13  -1.88     -1.34  
glance(fit)
#> # A tibble: 3 × 9
#>   compound     nobs df.residual sigma   rss logLik   AIC   BIC converged
#>   <chr>       <int>       <int> <dbl> <dbl>  <dbl> <dbl> <dbl> <lgl>    
#> 1 activator_c    33          29  3.02  264.  -81.1  172.  180. TRUE     
#> 2 inhibitor_a    33          29  4.26  526.  -92.5  195.  203. TRUE     
#> 3 inhibitor_b    33          29  3.69  395.  -87.8  186.  193. TRUE     
augment(fit)
#> # A tibble: 99 × 6
#>    compound       dose replicate response .fitted .resid
#>    <chr>         <dbl>     <int>    <dbl>   <dbl>  <dbl>
#>  1 activator_c 0               1     9.18    7.69  1.49 
#>  2 activator_c 0               2     9.62    7.69  1.93 
#>  3 activator_c 0               3    11.3     7.69  3.58 
#>  4 activator_c 0.00508         1    10.8     7.73  3.02 
#>  5 activator_c 0.00508         2     7.42    7.73 -0.308
#>  6 activator_c 0.00508         3     8.34    7.73  0.612
#>  7 activator_c 0.0152          1    10.3     7.88  2.42 
#>  8 activator_c 0.0152          2     7.5     7.88 -0.376
#>  9 activator_c 0.0152          3     3.16    7.88 -4.72 
#> 10 activator_c 0.0457          1     7.19    8.54 -1.35 
#> # ℹ 89 more rows

de <- bio_diff_expr(expression_sim, gene, sample, count, condition,
  treatment = "stress", reference = "control", transform = "log2_cpm"
)
tidy(de)
#> # A tibble: 1,000 × 7
#>    feature   estimate std.error statistic  p.value  p.adjusted ave_expr
#>    <chr>        <dbl>     <dbl>     <dbl>    <dbl>       <dbl>    <dbl>
#>  1 gene_0054    -3.10     0.256    -12.1  4.88e-10 0.000000245     6.35
#>  2 gene_0053    -2.75     0.227    -12.1  4.91e-10 0.000000245     7.42
#>  3 gene_0071    -2.87     0.263    -10.9  2.60e- 9 0.000000539     7.72
#>  4 gene_0051    -3.12     0.288    -10.8  2.80e- 9 0.000000539     7.60
#>  5 gene_0067    -2.70     0.250    -10.8  2.95e- 9 0.000000539    10.9 
#>  6 gene_0047     2.94     0.273     10.7  3.23e- 9 0.000000539     9.69
#>  7 gene_0042     2.84     0.270     10.5  4.38e- 9 0.000000572    11.3 
#>  8 gene_0062    -2.71     0.260    -10.5  4.94e- 9 0.000000572     7.71
#>  9 gene_0056    -2.91     0.279    -10.4  5.14e- 9 0.000000572     6.20
#> 10 gene_0008     2.70     0.279      9.70 1.55e- 8 0.00000155      8.68
#> # ℹ 990 more rows
glance(de)
#> # A tibble: 1 × 7
#>   n_features n_significant df_prior s0_sq df_residual treatment reference
#>        <int>         <int>    <dbl> <dbl>       <int> <chr>     <chr>    
#> 1       1000            69     11.8 0.150           6 stress    control  
```
