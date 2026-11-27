# Simulated qPCR experiment

**Simulated data.** Ct values for an *E. coli* heat-shock experiment
with two reference genes (`gapA`, `rpoD`) and four targets, three
treatments, four biological samples per treatment and three technical
replicates per sample, run on two plates. Ct values were generated in
`data-raw/make_data.R` from known log2 fold changes (control 0; heat 42
C: dnaK 3, groL 2, ibpA 4, ompF -1; heat 46 C: dnaK 4.5, groL 3, ibpA 6,
ompF -2) with sample-level input variation and technical noise.
Amplification efficiency is 2 for every gene.

## Usage

``` r
qpcr_sim
```

## Format

A tibble with 432 rows and 6 columns:

- plate:

  Plate (`"plate_1"` or `"plate_2"`).

- sample:

  Biological sample.

- treatment:

  `"control"`, `"heat_42C"` or `"heat_46C"`.

- gene:

  Gene name.

- replicate:

  Technical replicate (1 to 3).

- ct:

  Quantification cycle.

## Source

Simulated with a fixed seed in `data-raw/make_data.R`.

## Examples

``` r
qpcr_sim
#> # A tibble: 216 × 6
#>    plate   sample    treatment gene  replicate    ct
#>    <chr>   <chr>     <chr>     <chr>     <int> <dbl>
#>  1 plate_1 control_1 control   gapA          1  17.8
#>  2 plate_1 control_1 control   gapA          2  17.6
#>  3 plate_1 control_1 control   gapA          3  17.8
#>  4 plate_1 control_1 control   rpoD          1  23.1
#>  5 plate_1 control_1 control   rpoD          2  23.1
#>  6 plate_1 control_1 control   rpoD          3  23.0
#>  7 plate_1 control_1 control   dnaK          1  25.0
#>  8 plate_1 control_1 control   dnaK          2  25.2
#>  9 plate_1 control_1 control   dnaK          3  25.3
#> 10 plate_1 control_1 control   groL          1  22.0
#> # ℹ 206 more rows
```
