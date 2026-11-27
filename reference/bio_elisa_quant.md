# ELISA quantification from a standard curve

Fits a four-parameter logistic standard curve to the wells with known
concentrations and back-calculates the concentration of every well from
its response. Wells whose response lies outside the fitted asymptotes
cannot be quantified and get `NA` with `in_range = FALSE`. Unknowns are
flagged as extrapolated when their estimate falls outside the range of
the standards.

## Usage

``` r
bio_elisa_quant(
  .data,
  response,
  concentration,
  .by = NULL,
  dilution = NULL,
  level = 0.95
)
```

## Arguments

- .data:

  A data frame or tibble with one row per well.

- response:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Measured signal (for example optical density).

- concentration:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Known concentration for standards and `NA` for unknown samples.

- .by:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Optional columns defining separate plates or runs, each with its own
  standard curve.

- dilution:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Optional column of dilution factors; estimates are multiplied by it.

- level:

  Confidence level for the curve parameters.

## Value

A `<bio_assay>` holding `.data` plus `estimated_conc` (after dilution
correction), `in_range`, `extrapolated` and `is_standard`. The standard
curves are kept for
[tidy()](https://tidybiotek.in4metrix.dev/reference/broom_bio.md),
[glance()](https://tidybiotek.in4metrix.dev/reference/broom_bio.md) and
[autoplot()](https://tidybiotek.in4metrix.dev/reference/autoplot_bio.md).

## References

Findlay, J. W. A. and Dillard, R. F. (2007). Appropriate calibration
curve fitting in ligand binding assays. *The AAPS Journal*, 9,
E260-E267.
[doi:10.1208/aapsj0902029](https://doi.org/10.1208/aapsj0902029)

## Examples

``` r
# Run 1 of the DNase ELISA as the standard curve; unknowns from run 2
plate <- rbind(
  transform(subset(dnase_elisa, run == 1), known = conc),
  transform(subset(dnase_elisa, run == 2), known = NA)
)
quant <- bio_elisa_quant(plate, density, known)
quant[!quant$is_standard, c("conc", "density", "estimated_conc")]
#> # A bio_assay: elisa: 16 x 3
#>       conc density estimated_conc
#>      <dbl>   <dbl>          <dbl>
#>  1  0.0488   0.045         0.0808
#>  2  0.0488   0.05          0.0891
#>  3  0.195    0.137         0.246 
#>  4  0.195    0.123         0.219 
#>  5  0.391    0.225         0.425 
#>  6  0.391    0.207         0.387 
#>  7  0.781    0.401         0.846 
#>  8  0.781    0.383         0.799 
#>  9  1.56     0.672         1.70  
#> 10  1.56     0.681         1.73  
#> 11  3.12     1.12          3.99  
#> 12  3.12     1.08          3.73  
#> 13  6.25     1.55          8.92  
#> 14  6.25     1.53          8.44  
#> 15 12.5      1.93         21.6   
#> 16 12.5      1.91         20.5   
```
