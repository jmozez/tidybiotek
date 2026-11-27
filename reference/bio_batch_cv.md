# Leave-batch-out cross-validation

Batch effects (plates, runs, sequencing lanes, fermentation campaigns)
are shared by samples processed together. `bio_batch_cv()` keeps every
batch whole: each assessment set is one batch (leave-one-batch-out) or,
when `v` is given, a group of batches balanced by size.

## Usage

``` r
bio_batch_cv(.data, batch, v = NULL)
```

## Arguments

- .data:

  A data frame or tibble.

- batch:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Column identifying the batch of each row.

- v:

  Number of folds. `NULL` (the default) gives one fold per batch.

## Value

An `rset` of class `bio_batch_cv` with columns `splits` and `id`.

## References

Leek, J. T. et al. (2010). Tackling the widespread and critical impact
of batch effects in high-throughput data. *Nature Reviews Genetics*, 11,
733-739. [doi:10.1038/nrg2825](https://doi.org/10.1038/nrg2825)

## Examples

``` r
folds <- bio_batch_cv(qpcr_sim, plate)
folds
#> # Leave-batch-out cross-validation (2 folds by `plate`) 
#> # A tibble: 2 × 2
#>   splits            id           
#>   <list>            <chr>        
#> 1 <split [108/108]> Batch_plate_1
#> 2 <split [108/108]> Batch_plate_2
rsample::assessment(folds$splits[[1]])$plate |> unique()
#> [1] "plate_1"
```
