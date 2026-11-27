# Michaelis-Menten enzyme kinetics

Estimates the maximum rate `vmax` and the Michaelis constant `km` from
initial rates measured at several substrate concentrations. The default
fits \\v = V\_{max} S / (K_m + S)\\ by nonlinear least squares with
[`stats::nls()`](https://rdrr.io/r/stats/nls.html), started from the
Hanes-Woolf estimates. The classic linearisations (Lineweaver-Burk,
Eadie-Hofstee and Hanes-Woolf) are available for teaching and
comparison; their standard errors come from the delta method. Give
`enzyme_conc` to also get the turnover number `kcat` and the specificity
constant `kcat / km`.

## Usage

``` r
bio_enzyme_kinetics(
  .data,
  substrate,
  rate,
  .by = NULL,
  method = c("nls", "lineweaver_burk", "eadie_hofstee", "hanes_woolf"),
  enzyme_conc = NULL,
  level = 0.95
)
```

## Arguments

- .data:

  A data frame or tibble.

- substrate, rate:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Substrate concentration and initial rate columns.

- .by:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Optional columns defining separate enzymes or conditions.

- method:

  One of `"nls"`, `"lineweaver_burk"`, `"eadie_hofstee"` or
  `"hanes_woolf"`.

- enzyme_conc:

  Optional total enzyme concentration, in the same amount units as
  `vmax` per unit time.

- level:

  Confidence level.

## Value

A `<bio_assay>` with one row per group: the `.by` columns, `vmax`,
`vmax_se`, `km`, `km_se`, `kcat` and `kcat_km` (when `enzyme_conc` is
given), `sigma`, `nobs` and `method`.

## References

Michaelis, L. and Menten, M. L. (1913). Die Kinetik der Invertinwirkung.
*Biochemische Zeitschrift*, 49, 333-369.

Bates, D. M. and Watts, D. G. (1988). *Nonlinear Regression Analysis and
Its Applications*. Wiley.
[doi:10.1002/9780470316757](https://doi.org/10.1002/9780470316757)

## Examples

``` r
bio_enzyme_kinetics(puromycin_kinetics, substrate, rate, .by = state)
#> # A bio_assay: enzyme_kinetics: 2 x 8
#>   state      vmax vmax_se     km   km_se sigma  nobs method
#>   <chr>     <dbl>   <dbl>  <dbl>   <dbl> <dbl> <int> <chr> 
#> 1 treated    213.    6.95 0.0641 0.00828 10.9     12 nls   
#> 2 untreated  160.    6.48 0.0477 0.00778  9.77    11 nls   

# Compare with the Lineweaver-Burk linearisation
bio_enzyme_kinetics(puromycin_kinetics, substrate, rate,
  .by = state,
  method = "lineweaver_burk"
)
#> # A bio_assay: enzyme_kinetics: 2 x 8
#>   state      vmax vmax_se     km   km_se sigma  nobs method         
#>   <chr>     <dbl>   <dbl>  <dbl>   <dbl> <dbl> <int> <chr>          
#> 1 treated    196.    27.0 0.0484 0.0117   13.9    12 lineweaver_burk
#> 2 untreated  143.    11.7 0.0308 0.00554  13.1    11 lineweaver_burk
```
