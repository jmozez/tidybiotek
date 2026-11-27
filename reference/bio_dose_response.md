# Four-parameter logistic dose-response curves (IC50 and EC50)

Fits the four-parameter logistic (4PL) model \$\$y = \mathrm{bottom} +
\frac{\mathrm{top} - \mathrm{bottom}}{1 + (\mathrm{EC50} /
x)^{\mathrm{hill}}}\$\$ to each group by Levenberg-Marquardt least
squares, written in base R. `bottom` is always the lower asymptote and
`top` the upper one. A negative Hill slope means the response falls with
dose, so the midpoint is reported as an IC50; a positive slope gives an
EC50. Zero doses are allowed and are treated as true zeros. Confidence
intervals for the EC50 or IC50 are computed on the log scale.

## Usage

``` r
bio_dose_response(
  .data,
  dose,
  response,
  .by = NULL,
  weights = NULL,
  bottom = NULL,
  top = NULL,
  level = 0.95
)
```

## Arguments

- .data:

  A data frame or tibble.

- dose, response:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Columns holding dose (or concentration, on the original scale) and
  response.

- .by:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Optional columns that define separate curves, such as compound or cell
  line.

- weights:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Optional column of case weights.

- bottom, top:

  Optional fixed values for the asymptotes, for example `bottom = 0` for
  data normalised to percent of control.

- level:

  Confidence level.

## Value

A `<bio_assay>` with one row per curve: the `.by` columns, `bottom`,
`top`, `ec50`, `ec50_low`, `ec50_high`, `hill`, `type` (`"IC50"` or
`"EC50"`), `sigma`, `nobs` and `converged`. Use
[tidy()](https://tidybiotek.in4metrix.dev/reference/broom_bio.md) for
full parameter tables and
[autoplot()](https://tidybiotek.in4metrix.dev/reference/autoplot_bio.md)
to draw the curves.

## References

Ritz, C., Baty, F., Streibig, J. C. and Gerhard, D. (2015).
Dose-response analysis using R. *PLOS ONE*, 10, e0146021.
[doi:10.1371/journal.pone.0146021](https://doi.org/10.1371/journal.pone.0146021)

## Examples

``` r
fit <- bio_dose_response(dose_response_sim, dose, response, .by = compound)
fit
#> # A bio_assay: dose_response: 3 x 11
#>   compound    bottom   top   ec50 ec50_low ec50_high  hill type  sigma  nobs
#>   <chr>        <dbl> <dbl>  <dbl>    <dbl>     <dbl> <dbl> <chr> <dbl> <int>
#> 1 activator_c   7.69  85.9 1.20     1.06      1.36    1.38 EC50   3.02    33
#> 2 inhibitor_a   4.02 102.  0.0490   0.0413    0.0581 -1.16 IC50   4.26    33
#> 3 inhibitor_b   9.31  98.6 2.34     2.07      2.65   -1.61 IC50   3.69    33
#> # ℹ 1 more variable: converged <lgl>
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
```
