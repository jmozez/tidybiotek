# Physicochemical properties of proteins

Calculates sequence-based protein properties used for protein
engineering and expression screening: length, average molecular weight
(Da), isoelectric point and net charge at pH 7 (EMBOSS pKa values),
GRAVY hydropathy (Kyte and Doolittle, 1982), aromaticity and the
fractions of acidic, basic and hydrophobic residues. Non-standard
residues are ignored and a trailing stop (`*`) is dropped.

## Usage

``` r
bio_aa_props(.data, seq = NULL, prefix = "")
```

## Arguments

- .data:

  A data frame, tibble or
  [`bio_seq_tbl()`](https://tidybiotek.in4metrix.dev/reference/bio_seq_tbl.md).

- seq:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  The sequence column. Can be omitted for a `<bio_seq_tbl>`.

- prefix:

  Prefix for the new column names.

## Value

`.data` with nine numeric columns: `length`, `mol_weight`,
`isoelectric_point`, `charge_ph7`, `gravy`, `aromaticity`,
`frac_acidic`, `frac_basic` and `frac_hydrophobic`, each with `prefix`.

## References

Kyte, J. and Doolittle, R. F. (1982). A simple method for displaying the
hydropathic character of a protein. *Journal of Molecular Biology*, 157,
105-132.
[doi:10.1016/0022-2836(82)90515-0](https://doi.org/10.1016/0022-2836%2882%2990515-0)

## Examples

``` r
ecoli_cds |>
  dplyr::slice(1:5) |>
  bio_aa_props(protein) |>
  dplyr::select(gene, length, mol_weight, isoelectric_point, gravy)
#> # A tibble: 5 × 5
#>   gene  length mol_weight isoelectric_point  gravy
#>   <chr>  <dbl>      <dbl>             <dbl>  <dbl>
#> 1 rplA     234     24730.              10.3 -0.109
#> 2 rplB     273     29860.              11.6 -0.699
#> 3 rplC     209     22244.              10.7 -0.235
#> 4 rplD     201     22087.              10.4 -0.235
#> 5 rplE     179     20302.              10.0 -0.282
```
