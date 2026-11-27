# Physicochemical protein features

`step_bio_aa_props()` creates a *specification* of a recipe step that
replaces a column of protein sequences with the physicochemical
properties computed by
[`bio_aa_props()`](https://tidybiotek.in4metrix.dev/reference/bio_aa_props.md):
length, molecular weight, isoelectric point, charge at pH 7, GRAVY
hydropathy, aromaticity and the fractions of acidic, basic and
hydrophobic residues.

## Usage

``` r
step_bio_aa_props(
  recipe,
  ...,
  role = "predictor",
  trained = FALSE,
  properties = aa_prop_names,
  columns = NULL,
  keep_original_cols = FALSE,
  skip = FALSE,
  id = recipes::rand_id("bio_aa_props")
)

# S3 method for class 'step_bio_aa_props'
tidy(x, ...)
```

## Arguments

- recipe:

  A recipe object. The step is added to its sequence of operations.

- ...:

  One or more selector functions to choose sequence columns. See
  [`recipes::selections()`](https://recipes.tidymodels.org/reference/selections.html).

- role:

  Role for the new feature columns.

- trained:

  A logical showing whether the step has been trained.

- properties:

  Which properties to keep. Defaults to all nine.

- columns:

  The selected columns. Filled in by
  [`recipes::prep()`](https://recipes.tidymodels.org/reference/prep.html).

- keep_original_cols:

  Keep the original sequence columns?

- skip:

  Should the step be skipped when the recipe is baked by
  [`recipes::bake()`](https://recipes.tidymodels.org/reference/bake.html)?

- id:

  A character string that is unique to this step.

- x:

  A `step_bio_aa_props` object.

## Value

An updated version of `recipe` with the new step added. When you
[tidy()](https://recipes.tidymodels.org/reference/tidy.recipe.html) this
step, a tibble with columns `terms` and `property` is returned.

## See also

Other tidybiotek recipe steps:
[`step_bio_kmers()`](https://tidybiotek.in4metrix.dev/reference/step_bio_kmers.md),
[`step_bio_vst_norm()`](https://tidybiotek.in4metrix.dev/reference/step_bio_vst_norm.md)

## Examples

``` r
rec <- recipes::recipe(class ~ protein, data = ecoli_cds) |>
  step_bio_aa_props(protein) |>
  recipes::prep()
recipes::bake(rec, new_data = NULL)
#> # A tibble: 76 × 10
#>    class     protein_length protein_mol_weight protein_isoelectric_point
#>    <fct>              <dbl>              <dbl>                     <dbl>
#>  1 ribosomal            234             24730.                     10.3 
#>  2 ribosomal            273             29860.                     11.6 
#>  3 ribosomal            209             22244.                     10.7 
#>  4 ribosomal            201             22087.                     10.4 
#>  5 ribosomal            179             20302.                     10.0 
#>  6 ribosomal            177             18904.                     10.3 
#>  7 ribosomal            149             15769.                      6.55
#>  8 ribosomal            165             17712.                      9.51
#>  9 ribosomal            142             14875.                     10.3 
#> 10 ribosomal            121             12295.                      4.31
#> # ℹ 66 more rows
#> # ℹ 6 more variables: protein_charge_ph7 <dbl>, protein_gravy <dbl>,
#> #   protein_aromaticity <dbl>, protein_frac_acidic <dbl>,
#> #   protein_frac_basic <dbl>, protein_frac_hydrophobic <dbl>
```
