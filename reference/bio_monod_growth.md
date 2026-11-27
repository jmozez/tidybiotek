# Monod growth kinetics for batch bioreactors

Fits the Monod model of substrate-limited microbial growth to batch
fermentation time courses: \$\$dX/dt = \mu\_{max} \frac{S}{K_s + S} X,
\quad dS/dt = -\frac{1}{Y} \frac{dX}{dt}\$\$ where X is biomass, S is
the limiting substrate and Y is the biomass yield on substrate. The mass
balance \\X + Y S = X_0 + Y S_0\\ reduces the system to its exact
integrated form, which is solved for X at each time by bisection. The
parameters are estimated by least squares on the log scale (Nelder-Mead
followed by BFGS, from several starts). Biomass and substrate residuals
are scaled by their standard deviations so both contribute. Without a
substrate column the initial substrate `s0` must be given.

## Usage

``` r
bio_monod_growth(.data, time, biomass, substrate = NULL, .by = NULL, s0 = NULL)
```

## Arguments

- .data:

  A data frame or tibble.

- time, biomass:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Time and biomass columns.

- substrate:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Optional substrate concentration column.

- .by:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Optional columns identifying separate runs or reactors.

- s0:

  Initial substrate concentration, required when `substrate` is not
  given.

## Value

A `<bio_assay>` with one row per run: the `.by` columns, `mu_max`, `ks`,
`yield`, `x0`, `s0`, `doubling_time` (`log(2) / mu_max`), `sigma`,
`nobs` and `converged`. Standard errors are in
[tidy()](https://tidybiotek.in4metrix.dev/reference/broom_bio.md) and
fitted trajectories in
[augment()](https://tidybiotek.in4metrix.dev/reference/broom_bio.md).

## References

Monod, J. (1949). The growth of bacterial cultures. *Annual Review of
Microbiology*, 3, 371-394.
[doi:10.1146/annurev.mi.03.100149.002103](https://doi.org/10.1146/annurev.mi.03.100149.002103)

## Examples

``` r
growth <- bio_monod_growth(fermentation_sim, time_h, biomass_g_l,
  substrate = glucose_g_l, .by = reactor
)
growth
#> # A bio_assay: monod: 3 x 10
#>   reactor mu_max    ks yield     x0    s0 doubling_time  sigma  nobs converged
#>   <chr>    <dbl> <dbl> <dbl>  <dbl> <dbl>         <dbl>  <dbl> <int> <lgl>    
#> 1 R1       0.556 0.335 0.512 0.0749  9.97          1.25 0.0536    30 TRUE     
#> 2 R2       0.529 0.145 0.473 0.0866 20.1           1.31 0.0372    30 TRUE     
#> 3 R3       0.561 0.839 0.466 0.0801 30.1           1.24 0.0333    30 TRUE     
```
