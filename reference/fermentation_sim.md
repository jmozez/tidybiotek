# Simulated batch fermentations

**Simulated data.** Biomass and glucose time courses from three batch
bioreactors with different starting glucose (10, 20 and 30 g/l),
generated from the exact solution of the Monod model with \\\mu\_{max}\\
= 0.55 per hour, Ks = 0.4 g/l, yield 0.48 g/g and starting biomass 0.08
g/l, with 4 percent multiplicative noise on biomass and 0.25 g/l
additive noise on glucose.

## Usage

``` r
fermentation_sim
```

## Format

A tibble with 45 rows and 4 columns:

- reactor:

  Reactor (`"R1"`, `"R2"`, `"R3"`).

- time_h:

  Time (hours).

- biomass_g_l:

  Dry cell weight (g/l).

- glucose_g_l:

  Residual glucose (g/l).

## Source

Simulated with a fixed seed in `data-raw/make_data.R`.

## Examples

``` r
fermentation_sim
#> # A tibble: 45 × 4
#>    reactor time_h biomass_g_l glucose_g_l
#>    <chr>    <dbl>       <dbl>       <dbl>
#>  1 R1           0      0.0877      10.1  
#>  2 R1           1      0.129        9.66 
#>  3 R1           2      0.224        9.61 
#>  4 R1           3      0.384        9.35 
#>  5 R1           4      0.636        9.04 
#>  6 R1           5      1.08         8.05 
#>  7 R1           6      1.94         6.42 
#>  8 R1           7      3.11         3.98 
#>  9 R1           8      4.80         0    
#> 10 R1           9      5.33         0.318
#> # ℹ 35 more rows
```
