# Puromycin enzyme kinetics

Initial reaction velocity against substrate concentration for an enzyme
in cells treated or not treated with puromycin. A tidy copy of
[datasets::Puromycin](https://rdrr.io/r/datasets/Puromycin.html), a
standard example of Michaelis-Menten kinetics.

## Usage

``` r
puromycin_kinetics
```

## Format

A tibble with 23 rows and 3 columns:

- state:

  `"treated"` or `"untreated"`.

- substrate:

  Substrate concentration (ppm).

- rate:

  Initial reaction rate (counts/min/min).

## Source

R Core Team,
[datasets::Puromycin](https://rdrr.io/r/datasets/Puromycin.html).
Original data from Treloar, M. A. (1974). *Effects of Puromycin on
Galactosyltransferase in Golgi Membranes*. MSc thesis, University of
Toronto; analysed in Bates, D. M. and Watts, D. G. (1988). *Nonlinear
Regression Analysis and Its Applications*. Wiley.
[doi:10.1002/9780470316757](https://doi.org/10.1002/9780470316757)

## Examples

``` r
puromycin_kinetics
#> # A tibble: 23 × 3
#>    state   substrate  rate
#>    <chr>       <dbl> <dbl>
#>  1 treated      0.02    76
#>  2 treated      0.02    47
#>  3 treated      0.06    97
#>  4 treated      0.06   107
#>  5 treated      0.11   123
#>  6 treated      0.11   139
#>  7 treated      0.22   159
#>  8 treated      0.22   152
#>  9 treated      0.56   191
#> 10 treated      0.56   201
#> # ℹ 13 more rows
```
