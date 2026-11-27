# K-mer features from sequences

`step_bio_kmers()` creates a *specification* of a recipe step that turns
a column of DNA, RNA or protein sequences into k-mer count or frequency
features. The k-mer vocabulary is learnt from the training set when the
recipe is prepped, so new data always get the same columns.

## Usage

``` r
step_bio_kmers(
  recipe,
  ...,
  role = "predictor",
  trained = FALSE,
  k = 3,
  normalise = TRUE,
  vocabulary = NULL,
  keep_original_cols = FALSE,
  skip = FALSE,
  id = recipes::rand_id("bio_kmers")
)

# S3 method for class 'step_bio_kmers'
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

- k:

  K-mer length.

- normalise:

  If `TRUE` (the default), use k-mer frequencies rather than counts,
  which removes the effect of sequence length.

- vocabulary:

  The learnt k-mers, one character vector per column. Filled in by
  [`recipes::prep()`](https://recipes.tidymodels.org/reference/prep.html).

- keep_original_cols:

  Keep the original sequence columns?

- skip:

  Should the step be skipped when the recipe is baked by
  [`recipes::bake()`](https://recipes.tidymodels.org/reference/bake.html)?

- id:

  A character string that is unique to this step.

- x:

  A `step_bio_kmers` object.

## Value

An updated version of `recipe` with the new step added. When you
[tidy()](https://recipes.tidymodels.org/reference/tidy.recipe.html) this
step, a tibble with columns `terms`, `k` and `n_kmers` is returned.

## See also

Other tidybiotek recipe steps:
[`step_bio_aa_props()`](https://tidybiotek.in4metrix.dev/reference/step_bio_aa_props.md),
[`step_bio_vst_norm()`](https://tidybiotek.in4metrix.dev/reference/step_bio_vst_norm.md)

## Examples

``` r
rec <- recipes::recipe(gc ~ sequence, data = bio_gc_content(ecoli_cds, sequence)) |>
  step_bio_kmers(sequence, k = 2) |>
  recipes::prep()
recipes::bake(rec, new_data = NULL)
#> # A tibble: 76 × 17
#>       gc sequence_AA sequence_AC sequence_AG sequence_AT sequence_CA sequence_CC
#>    <dbl>       <dbl>       <dbl>       <dbl>       <dbl>       <dbl>       <dbl>
#>  1 0.516      0.114       0.0710      0.0582      0.0341      0.0611      0.0469
#>  2 0.527      0.0914      0.0658      0.0487      0.0402      0.0658      0.0609
#>  3 0.514      0.100       0.0636      0.0700      0.0223      0.0445      0.0556
#>  4 0.530      0.0843      0.0545      0.0628      0.0364      0.0479      0.0413
#>  5 0.498      0.0928      0.0631      0.0501      0.0538      0.0594      0.0501
#>  6 0.521      0.0807      0.0525      0.0582      0.0413      0.0469      0.0450
#>  7 0.498      0.111       0.0535      0.0668      0.0356      0.0423      0.0245
#>  8 0.524      0.0966      0.0463      0.0604      0.0382      0.0423      0.0463
#>  9 0.522      0.103       0.0491      0.0701      0.0421      0.0701      0.0607
#> 10 0.462      0.148       0.0356      0.0795      0.0411      0.0384      0.0164
#> # ℹ 66 more rows
#> # ℹ 10 more variables: sequence_CG <dbl>, sequence_CT <dbl>, sequence_GA <dbl>,
#> #   sequence_GC <dbl>, sequence_GG <dbl>, sequence_GT <dbl>, sequence_TA <dbl>,
#> #   sequence_TC <dbl>, sequence_TG <dbl>, sequence_TT <dbl>
tidy(rec, number = 1)
#> # A tibble: 1 × 4
#>   terms        k n_kmers id             
#>   <chr>    <dbl>   <int> <chr>          
#> 1 sequence     2      16 bio_kmers_dF7Sp
```
