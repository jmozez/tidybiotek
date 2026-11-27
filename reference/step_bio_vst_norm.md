# Variance-stabilising normalisation of count features

`step_bio_vst_norm()` creates a *specification* of a recipe step that
normalises count features (such as RNA-seq read counts, with samples as
rows and genes as columns) and stabilises their variance. Three methods
are offered:

## Usage

``` r
step_bio_vst_norm(
  recipe,
  ...,
  role = NA,
  trained = FALSE,
  method = c("vst", "log_cpm", "anscombe"),
  prior_count = 0.5,
  columns = NULL,
  log_geo_means = NULL,
  dispersion = NULL,
  skip = FALSE,
  id = recipes::rand_id("bio_vst_norm")
)

# S3 method for class 'step_bio_vst_norm'
tidy(x, ...)
```

## Arguments

- recipe:

  A recipe object. The step is added to its sequence of operations.

- ...:

  Selector functions choosing the count columns.

- role:

  Role for the new feature columns.

- trained:

  A logical showing whether the step has been trained.

- method:

  One of `"vst"`, `"log_cpm"` or `"anscombe"`.

- prior_count:

  Prior count added for `"log_cpm"`.

- columns, log_geo_means, dispersion:

  Learnt values filled in by
  [`recipes::prep()`](https://recipes.tidymodels.org/reference/prep.html).

- skip:

  Should the step be skipped when the recipe is baked by
  [`recipes::bake()`](https://recipes.tidymodels.org/reference/bake.html)?

- id:

  A character string that is unique to this step.

- x:

  A `step_bio_vst_norm` object.

## Value

An updated version of `recipe` with the new step added. When you
[tidy()](https://recipes.tidymodels.org/reference/tidy.recipe.html) this
step, a tibble with columns `terms`, `method` and `log_geo_mean` is
returned.

## Details

- `"vst"`: median-of-ratios size factors (Anders and Huber, 2010) and
  the closed-form variance-stabilising transform for a negative binomial
  with dispersion trend \\\alpha(\mu) = a_0 + a_1 / \mu\\, with the
  trend fitted to the training data by least squares.

- `"log_cpm"`: log2 counts per million with a prior count.

- `"anscombe"`: the Anscombe transform \\2 \sqrt{x + 3/8}\\ for Poisson
  counts.

Geometric means and the dispersion trend are learnt from the training
set only, so new samples are transformed consistently.

## References

Anders, S. and Huber, W. (2010). Differential expression analysis for
sequence count data. *Genome Biology*, 11, R106.
[doi:10.1186/gb-2010-11-10-r106](https://doi.org/10.1186/gb-2010-11-10-r106)

Anscombe, F. J. (1948). The transformation of Poisson, binomial and
negative-binomial data. *Biometrika*, 35, 246-254.
[doi:10.1093/biomet/35.3-4.246](https://doi.org/10.1093/biomet/35.3-4.246)

## See also

Other tidybiotek recipe steps:
[`step_bio_aa_props()`](https://tidybiotek.in4metrix.dev/reference/step_bio_aa_props.md),
[`step_bio_kmers()`](https://tidybiotek.in4metrix.dev/reference/step_bio_kmers.md)

## Examples

``` r
wide <- bio_counts_wide(expression_sim)
rec <- recipes::recipe(condition ~ ., data = wide) |>
  recipes::update_role(sample, new_role = "id") |>
  step_bio_vst_norm(recipes::all_numeric_predictors()) |>
  recipes::prep()
recipes::bake(rec, new_data = NULL)[, 1:6]
#> # A tibble: 8 × 6
#>   sample gene_0001 gene_0002 gene_0003 gene_0004 gene_0005
#>   <chr>      <dbl>     <dbl>     <dbl>     <dbl>     <dbl>
#> 1 S01         8.10      7.62     10.3       10.4      10.7
#> 2 S02         8.30      8.09     10.00      10.8      10.5
#> 3 S03         8.87      8.26      9.88      11.0      10.9
#> 4 S04         8.60      7.42     10.6       10.8      10.8
#> 5 S05        10.9       9.90     11.3       12.7      13.0
#> 6 S06        11.1       9.74     10.5       11.9      12.5
#> 7 S07        11.0      10.3      11.6       11.8      13.1
#> 8 S08        10.4       9.53     11.6       11.6      12.9
```
