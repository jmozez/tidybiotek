# Enrichment factor for virtual screening

The enrichment factor (EF) compares the hit rate among the top-scoring
`fraction` of a ranked list with the hit rate in the whole list. An EF
of 10 at 1 percent means the top 1 percent holds ten times more actives
than random selection would give. It is widely used to assess virtual
screening, biomarker and protein engineering models that prioritise
candidates for testing.

## Usage

``` r
enrichment_factor(data, ...)

# S3 method for class 'data.frame'
enrichment_factor(
  data,
  truth,
  ...,
  fraction = 0.01,
  na_rm = TRUE,
  event_level = "first",
  case_weights = NULL
)

enrichment_factor_vec(
  truth,
  estimate,
  fraction = 0.01,
  na_rm = TRUE,
  event_level = "first",
  case_weights = NULL,
  ...
)
```

## Arguments

- data:

  A data frame with the `truth` and probability or score columns.

- ...:

  The column with the predicted probability or score for the event
  class.

- truth:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  The true class (a factor with two levels).

- fraction:

  Fraction of the ranked list to examine, between 0 and 1.

- na_rm:

  Remove missing values before computing?

- event_level:

  `"first"` or `"second"`: which level of `truth` is the event (active).

- case_weights:

  Not supported; must be `NULL`.

- estimate:

  For the `_vec()` function, a numeric vector of scores.

## Value

For the data frame method, a tibble with columns `.metric`, `.estimator`
and `.estimate`. For `enrichment_factor_vec()`, a single number.

## References

Truchon, J.-F. and Bayly, C. I. (2007). Evaluating virtual screening
methods: good and bad metrics for the "early recognition" problem.
*Journal of Chemical Information and Modeling*, 47, 488-508.
[doi:10.1021/ci600426e](https://doi.org/10.1021/ci600426e)

## See also

Other tidybiotek metrics:
[`assay_cv_error()`](https://tidybiotek.in4metrix.dev/reference/assay_cv_error.md),
[`bedroc()`](https://tidybiotek.in4metrix.dev/reference/bedroc.md)

## Examples

``` r
set.seed(2)
screen <- tibble::tibble(
  active = factor(rep(c("yes", "no"), c(20, 480)), levels = c("yes", "no")),
  score = c(rnorm(20, 2), rnorm(480))
)
enrichment_factor(screen, active, score, fraction = 0.05)
#> # A tibble: 1 × 3
#>   .metric           .estimator .estimate
#>   <chr>             <chr>          <dbl>
#> 1 enrichment_factor binary            12
```
