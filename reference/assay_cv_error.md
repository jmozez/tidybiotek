# Assay back-calculation error

`assay_cv_error()` is the root-mean-square relative error, in percent,
between estimated and nominal concentrations: \$\$100 \sqrt{\frac{1}{n}
\sum_i \left(\frac{\hat{y}\_i - y_i}{y_i}\right)^2}\$\$ It combines bias
and imprecision into one CV-like figure and is the usual way to judge
how well a standard curve back-calculates its standards or quality
controls (acceptance is often within 20 percent).

## Usage

``` r
assay_cv_error(data, ...)

# S3 method for class 'data.frame'
assay_cv_error(data, truth, estimate, na_rm = TRUE, case_weights = NULL, ...)

assay_cv_error_vec(truth, estimate, na_rm = TRUE, case_weights = NULL, ...)
```

## Arguments

- data:

  A data frame with the `truth` and `estimate` columns.

- ...:

  Not used.

- truth, estimate:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Nominal and estimated concentrations (vectors for the `_vec()`
  function).

- na_rm:

  Remove missing values before computing?

- case_weights:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Optional case weights.

## Value

For the data frame method, a tibble with columns `.metric`, `.estimator`
and `.estimate`. For `assay_cv_error_vec()`, a single number.

## References

US Food and Drug Administration (2018). *Bioanalytical Method
Validation: Guidance for Industry*.

## See also

Other tidybiotek metrics:
[`bedroc()`](https://tidybiotek.in4metrix.dev/reference/bedroc.md),
[`enrichment_factor()`](https://tidybiotek.in4metrix.dev/reference/enrichment_factor.md)

## Examples

``` r
standards <- tibble::tibble(nominal = c(1, 2, 4, 8), measured = c(1.1, 1.9, 4.3, 7.6))
assay_cv_error(standards, nominal, measured)
#> # A tibble: 1 × 3
#>   .metric        .estimator .estimate
#>   <chr>          <chr>          <dbl>
#> 1 assay_cv_error standard        7.18
```
