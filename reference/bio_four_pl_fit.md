# Fit function for the `nls_4pl` engine

The engine fit function used by
[`bio_four_pl()`](https://tidybiotek.in4metrix.dev/reference/bio_four_pl.md).
It is exported so that parsnip can call it, and can also be used
directly with a formula of the form `response ~ dose`.

## Usage

``` r
bio_four_pl_fit(formula, data, bottom = NULL, top = NULL, ...)

# S3 method for class 'bio_four_pl_fit'
predict(object, new_data, ...)
```

## Arguments

- formula:

  A formula with one numeric predictor.

- data:

  A data frame.

- bottom, top:

  Optional fixed asymptotes.

- ...:

  Not used.

- object:

  A `bio_four_pl_fit` object.

- new_data:

  A data frame with the predictor column.

## Value

`bio_four_pl_fit()` returns an object of class `bio_four_pl_fit`. The
[`predict()`](https://rdrr.io/r/stats/predict.html) method returns a
numeric vector.

## Examples

``` r
one <- subset(dose_response_sim, compound == "inhibitor_b")
fit <- bio_four_pl_fit(response ~ dose, one)
fit
#> Four-parameter logistic fit of `response` on `dose`
#> # A tibble: 4 × 3
#>   term   estimate std.error
#>   <chr>     <dbl>     <dbl>
#> 1 bottom     9.31     1.56 
#> 2 top       98.6      0.958
#> 3 ec50       2.34     0.142
#> 4 hill      -1.61     0.133
predict(fit, data.frame(dose = 0.01))
#> [1] 98.54976
```
