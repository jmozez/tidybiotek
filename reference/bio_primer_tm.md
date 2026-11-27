# Primer melting temperature

Estimates the melting temperature (Tm, in degrees Celsius) of DNA
oligonucleotides. The default nearest-neighbour method uses the unified
parameters of SantaLucia (1998) with the entropy salt correction for
monovalent cations. The Wallace rule (2 degrees per A or T, 4 per G or
C) and the basic GC formula are offered for quick checks.

## Usage

``` r
bio_primer_tm(
  .data,
  seq = NULL,
  method = c("nearest_neighbour", "wallace", "basic"),
  na_conc = 50,
  primer_conc = 50,
  .into = "tm"
)
```

## Arguments

- .data:

  A data frame, tibble or
  [`bio_seq_tbl()`](https://tidybiotek.in4metrix.dev/reference/bio_seq_tbl.md).

- seq:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  The sequence column. Can be omitted for a `<bio_seq_tbl>`.

- method:

  One of `"nearest_neighbour"`, `"wallace"` or `"basic"`.

- na_conc:

  Monovalent cation (Na+) concentration in mM.

- primer_conc:

  Total oligonucleotide concentration in nM.

- .into:

  Name of the new Tm column.

## Value

`.data` with a numeric column `.into`. Sequences with bases other than
A, C, G and T give `NA` with a warning.

## References

SantaLucia, J. (1998). A unified view of polymer, dumbbell, and
oligonucleotide DNA nearest-neighbor thermodynamics. *Proceedings of the
National Academy of Sciences*, 95, 1460-1465.
[doi:10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460)

Wallace, R. B. et al. (1979). Hybridization of synthetic
oligodeoxyribonucleotides to phi chi 174 DNA: the effect of single base
pair mismatch. *Nucleic Acids Research*, 6, 3543-3557.
[doi:10.1093/nar/6.11.3543](https://doi.org/10.1093/nar/6.11.3543)

## Examples

``` r
primers <- tibble::tibble(
  name = c("M13F", "M13R", "T7"),
  primer = c("GTAAAACGACGGCCAGT", "CAGGAAACAGCTATGAC", "TAATACGACTCACTATAGGG")
)
primers |>
  bio_primer_tm(primer) |>
  bio_primer_tm(primer, method = "wallace", .into = "tm_wallace")
#> # A tibble: 3 × 4
#>   name  primer                  tm tm_wallace
#>   <chr> <chr>                <dbl>      <dbl>
#> 1 M13F  GTAAAACGACGGCCAGT     49.2         52
#> 2 M13R  CAGGAAACAGCTATGAC     43.6         50
#> 3 T7    TAATACGACTCACTATAGGG  44.6         56
```
