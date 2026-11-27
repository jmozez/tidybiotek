# DNase ELISA standard curves

Optical density from an ELISA for recombinant DNase in rat serum, with
eleven runs of a serial dilution standard curve. A tidy copy of
[datasets::DNase](https://rdrr.io/r/datasets/DNase.html).

## Usage

``` r
dnase_elisa
```

## Format

A tibble with 176 rows and 3 columns:

- run:

  Assay run (1 to 11).

- conc:

  Known protein concentration (ng/ml).

- density:

  Measured optical density.

## Source

R Core Team, [datasets::DNase](https://rdrr.io/r/datasets/DNase.html).
Original data from Davidian, M. and Giltinan, D. M. (1995). *Nonlinear
Models for Repeated Measurement Data*. Chapman and Hall.
[doi:10.1201/9780203745502](https://doi.org/10.1201/9780203745502)

## Examples

``` r
dnase_elisa
#> # A tibble: 176 × 3
#>      run   conc density
#>    <int>  <dbl>   <dbl>
#>  1     1 0.0488   0.017
#>  2     1 0.0488   0.018
#>  3     1 0.195    0.121
#>  4     1 0.195    0.124
#>  5     1 0.391    0.206
#>  6     1 0.391    0.215
#>  7     1 0.781    0.377
#>  8     1 0.781    0.374
#>  9     1 1.56     0.614
#> 10     1 1.56     0.609
#> # ℹ 166 more rows
```
