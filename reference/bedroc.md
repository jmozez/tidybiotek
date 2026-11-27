# Boltzmann-enhanced discrimination of ROC (BEDROC)

BEDROC (Truchon and Bayly, 2007) measures early recognition: how well a
score ranks actives at the very top of a list. Each active contributes
an exponentially decreasing weight of its rank, controlled by `alpha`
(the default 20 gives 80 percent of the weight to the top 8 percent of
the list). Values range from 0 to 1, where 1 means all actives are
ranked first.

## Usage

``` r
bedroc(data, ...)

# S3 method for class 'data.frame'
bedroc(
  data,
  truth,
  ...,
  alpha = 20,
  na_rm = TRUE,
  event_level = "first",
  case_weights = NULL
)

bedroc_vec(
  truth,
  estimate,
  alpha = 20,
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

- alpha:

  Early-recognition parameter.

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
and `.estimate`. For `bedroc_vec()`, a single number.

## References

Truchon, J.-F. and Bayly, C. I. (2007). Evaluating virtual screening
methods: good and bad metrics for the "early recognition" problem.
*Journal of Chemical Information and Modeling*, 47, 488-508.
[doi:10.1021/ci600426e](https://doi.org/10.1021/ci600426e)

## See also

Other tidybiotek metrics:
[`assay_cv_error()`](https://tidybiotek.in4metrix.dev/reference/assay_cv_error.md),
[`enrichment_factor()`](https://tidybiotek.in4metrix.dev/reference/enrichment_factor.md)

## Examples

``` r
set.seed(2)
screen <- tibble::tibble(
  active = factor(rep(c("yes", "no"), c(20, 480)), levels = c("yes", "no")),
  score = c(rnorm(20, 2), rnorm(480))
)
bedroc(screen, active, score)
#> # A tibble: 1 × 3
#>   .metric .estimator .estimate
#>   <chr>   <chr>          <dbl>
#> 1 bedroc  binary         0.702
```
