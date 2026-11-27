# Homology-aware cross-validation for sequences

Similar sequences (homologues) in both the analysis and assessment sets
make cross-validation over-optimistic. `bio_homology_cv()` first
clusters the sequences by k-mer profile similarity and then assigns
whole clusters to folds, so no cluster is split between training and
testing. Sequences are clustered by average-linkage (UPGMA) hierarchical
clustering on the cosine distance between k-mer frequency vectors, cut
at `threshold`. The k-mer distance is an alignment-free stand-in for
sequence identity: smaller thresholds give tighter clusters.

## Usage

``` r
bio_homology_cv(.data, seq = NULL, v = 5, k = 3, threshold = 0.3)
```

## Arguments

- .data:

  A data frame, tibble or
  [`bio_seq_tbl()`](https://tidybiotek.in4metrix.dev/reference/bio_seq_tbl.md).

- seq:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  The sequence column. Can be omitted for a `<bio_seq_tbl>`.

- v:

  Number of folds.

- k:

  K-mer length used for the similarity.

- threshold:

  Cosine distance at which the clustering tree is cut, between 0 and 1.

## Value

An `rset` of class `bio_homology_cv` with columns `splits` and `id`. The
cluster of every row is kept in `attr(x, "clusters")`.

## References

Walsh, I. et al. (2021). DOME: recommendations for supervised machine
learning validation in biology. *Nature Methods*, 18, 1122-1127.
[doi:10.1038/s41592-021-01205-4](https://doi.org/10.1038/s41592-021-01205-4)

## Examples

``` r
set.seed(1)
folds <- bio_homology_cv(ecoli_cds, protein, v = 5, k = 2, threshold = 0.2)
folds
#> # 5-fold homology-aware cross-validation (k = 2, threshold = 0.2) 
#> # A tibble: 5 × 2
#>   splits          id   
#>   <list>          <chr>
#> 1 <split [60/16]> Fold1
#> 2 <split [61/15]> Fold2
#> 3 <split [61/15]> Fold3
#> 4 <split [61/15]> Fold4
#> 5 <split [61/15]> Fold5
table(attr(folds, "clusters"))[1:5]
#> 
#> 1 2 3 4 5 
#> 1 1 1 1 1 
```
