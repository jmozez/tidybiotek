# Simulated RNA-seq counts and gene sets

**Simulated data.** `expression_sim` holds negative binomial read counts
for 1000 genes in 8 samples (4 control, 4 stress) with random library
sizes. Genes `gene_0001` to `gene_0050` are truly up-regulated under
stress (log2 fold change 1 to 3), `gene_0051` to `gene_0075` truly
down-regulated and the rest unchanged. `gene_sets_sim` holds 30 gene
sets in long format: `stress_response` (30 up-regulated genes and 10
unchanged ones), `repressed_metabolism` (15 down-regulated and 15
unchanged) and 28 random background sets.

## Usage

``` r
expression_sim

gene_sets_sim
```

## Format

`expression_sim` is a tibble with 8000 rows and 4 columns:

- gene:

  Gene identifier.

- sample:

  Sample identifier.

- condition:

  `"control"` or `"stress"`.

- count:

  Read count.

`gene_sets_sim` is a tibble with columns `gene_set` and `feature`.

An object of class `tbl_df` (inherits from `tbl`, `data.frame`) with
1081 rows and 2 columns.

## Source

Simulated with a fixed seed in `data-raw/make_data.R`.

## Examples

``` r
expression_sim
#> # A tibble: 8,000 × 4
#>    gene      sample condition count
#>    <chr>     <chr>  <chr>     <dbl>
#>  1 gene_0001 S01    control     316
#>  2 gene_0002 S01    control     221
#>  3 gene_0003 S01    control    1584
#>  4 gene_0004 S01    control    1635
#>  5 gene_0005 S01    control    2086
#>  6 gene_0006 S01    control    2513
#>  7 gene_0007 S01    control    2105
#>  8 gene_0008 S01    control     111
#>  9 gene_0009 S01    control     491
#> 10 gene_0010 S01    control     120
#> # ℹ 7,990 more rows
gene_sets_sim
#> # A tibble: 1,081 × 2
#>    gene_set        feature  
#>    <chr>           <chr>    
#>  1 stress_response gene_0002
#>  2 stress_response gene_0015
#>  3 stress_response gene_0008
#>  4 stress_response gene_0036
#>  5 stress_response gene_0047
#>  6 stress_response gene_0025
#>  7 stress_response gene_0040
#>  8 stress_response gene_0034
#>  9 stress_response gene_0032
#> 10 stress_response gene_0022
#> # ℹ 1,071 more rows
```
