# Four-parameter logistic regression model

`bio_four_pl()` defines a four-parameter logistic (4PL) dose-response
model for use with
[`parsnip::fit()`](https://generics.r-lib.org/reference/fit.html),
workflows and tuning. The only engine is `"nls_4pl"`, which uses the
Levenberg-Marquardt fitter behind
[`bio_dose_response()`](https://tidybiotek.in4metrix.dev/reference/bio_dose_response.md).
The model takes a single numeric predictor (the dose or concentration,
on the original scale) and a numeric outcome. Engine arguments `bottom`
and `top` fix the asymptotes, for example
`set_engine("nls_4pl", bottom = 0)`.

## Usage

``` r
bio_four_pl(mode = "regression", engine = "nls_4pl")
```

## Arguments

- mode:

  A single character string for the type of model. The only possible
  value is `"regression"`.

- engine:

  A single character string for the engine. The only possible value is
  `"nls_4pl"`.

## Value

A model specification object with classes `bio_four_pl` and
`model_spec`.

## Examples

``` r
one <- subset(dose_response_sim, compound == "inhibitor_b")
spec <- bio_four_pl()
spec
#> Four-parameter logistic model specification (regression)
#> 
#> Computational engine: nls_4pl 
#> 
fitted <- parsnip::fit(spec, response ~ dose, data = one)
tidy(fitted$fit)
#> # A tibble: 4 × 7
#>   term   estimate std.error statistic    p.value conf.low conf.high
#>   <chr>     <dbl>     <dbl>     <dbl>      <dbl>    <dbl>     <dbl>
#> 1 bottom     9.31     1.56       5.96  1.78 e- 6     6.11     12.5 
#> 2 top       98.6      0.958    103.    1.000e-38    96.6     101.  
#> 3 ec50       2.34     0.142     NA    NA             2.07      2.65
#> 4 hill      -1.61     0.133    -12.1   7.89 e-13    -1.88     -1.34
predict(fitted, new_data = data.frame(dose = c(0.001, 0.01, 0.1)))
#> # A tibble: 3 × 1
#>   .pred
#>   <dbl>
#> 1  98.6
#> 2  98.5
#> 3  98.0
```
