# Assay and omics result tibbles

Results from the assay, bioprocess and omics functions are tibble
subclasses. A `<bio_assay>` carries its fitted standard curves or
kinetic models, and a `<bio_omics>` carries the settings and prior
estimates of a differential expression or enrichment analysis. The
metadata is used by
[tidy()](https://tidybiotek.in4metrix.dev/reference/broom_bio.md),
[glance()](https://tidybiotek.in4metrix.dev/reference/broom_bio.md),
[augment()](https://tidybiotek.in4metrix.dev/reference/broom_bio.md) and
[autoplot()](https://tidybiotek.in4metrix.dev/reference/autoplot_bio.md),
and is kept through `dplyr` verbs while the required columns are
present.

## Usage

``` r
new_bio_assay(x, assay, meta = list(), required = character())

is_bio_assay(x)

new_bio_omics(x, analysis, meta = list(), required = character())

is_bio_omics(x)
```

## Arguments

- x:

  A data frame (constructors) or any object (predicates).

- assay, analysis:

  A string naming the kind of result, for example `"dose_response"` or
  `"diff_expr"`.

- meta:

  A list of metadata.

- required:

  Columns that must be present for the class to be kept.

## Value

The constructors return a `<bio_assay>` or `<bio_omics>` tibble. The
predicates return a single logical.

## Examples

``` r
fit <- bio_dose_response(dose_response_sim, dose, response, .by = compound)
is_bio_assay(fit)
#> [1] TRUE
fit |> dplyr::arrange(ec50) |> is_bio_assay()
#> [1] TRUE

x <- new_bio_assay(tibble::tibble(id = 1:2), "custom", required = "id")
x
#> # A bio_assay: custom: 2 x 1
#>      id
#>   <int>
#> 1     1
#> 2     2
```
