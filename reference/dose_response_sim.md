# Simulated dose-response screen

**Simulated data.** Responses (percent of control) of three compounds
over a ten-point, three-fold dilution series from 100 uM plus a vehicle
control, in triplicate. Generated from four-parameter logistic curves
with known values (inhibitor_a IC50 0.05 uM, Hill -1.1; inhibitor_b IC50
2.5 uM, Hill -1.6; activator_c EC50 1.2 uM, Hill 1.3) plus Gaussian
noise (SD 3.5). The compound names are placeholders, not real drugs.

## Usage

``` r
dose_response_sim
```

## Format

A tibble with 99 rows and 4 columns:

- compound:

  Compound name.

- dose:

  Dose (uM); 0 is the vehicle control.

- replicate:

  Replicate (1 to 3).

- response:

  Response (percent of control).

## Source

Simulated with a fixed seed in `data-raw/make_data.R`.

## Examples

``` r
dose_response_sim
#> # A tibble: 99 × 4
#>    compound       dose replicate response
#>    <chr>         <dbl>     <int>    <dbl>
#>  1 activator_c 0               1     9.18
#>  2 activator_c 0               2     9.62
#>  3 activator_c 0               3    11.3 
#>  4 activator_c 0.00508         1    10.8 
#>  5 activator_c 0.00508         2     7.42
#>  6 activator_c 0.00508         3     8.34
#>  7 activator_c 0.0152          1    10.3 
#>  8 activator_c 0.0152          2     7.5 
#>  9 activator_c 0.0152          3     3.16
#> 10 activator_c 0.0457          1     7.19
#> # ℹ 89 more rows
```
