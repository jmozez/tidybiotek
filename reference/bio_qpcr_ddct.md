# Relative quantification of qPCR data (delta-delta Ct and Pfaffl)

Calculates relative expression from quantification cycle (Ct or Cq)
values. Technical replicates are averaged within each sample and gene.
Each target gene is normalised to the reference gene(s) and then to the
mean of the control group. With the default amplification efficiency of
2 (100 percent) this is the 2^-ddCt method of Livak and Schmittgen
(2001); gene-specific efficiencies give the Pfaffl (2001) ratio. Several
reference genes are combined by the geometric mean of their relative
quantities (Vandesompele et al., 2002).

## Usage

``` r
bio_qpcr_ddct(
  .data,
  ct,
  gene,
  sample,
  group,
  reference_genes,
  control,
  efficiency = 2
)
```

## Arguments

- .data:

  A data frame with one row per well.

- ct, gene, sample, group:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Columns holding the Ct value, gene name, biological sample and
  experimental group.

- reference_genes:

  Character vector of reference (housekeeping) genes.

- control:

  Level of `group` used as the calibrator.

- efficiency:

  Amplification factor per cycle, between 1 and 2: a single number, or a
  named vector with one value per gene.

## Value

A `<bio_assay>` with one row per sample and target gene: the `sample`,
`group` and `gene` columns, `ct_target`, `ct_reference`, `delta_ct`,
`delta_delta_ct`, `fold_change` and `log2_fold_change`.
[tidy()](https://tidybiotek.in4metrix.dev/reference/broom_bio.md)
summarises the fold change per group with confidence intervals and Welch
t-tests against the control.

## References

Livak, K. J. and Schmittgen, T. D. (2001). Analysis of relative gene
expression data using real-time quantitative PCR and the 2^-ddCt method.
*Methods*, 25, 402-408.
[doi:10.1006/meth.2001.1262](https://doi.org/10.1006/meth.2001.1262)

Pfaffl, M. W. (2001). A new mathematical model for relative
quantification in real-time RT-PCR. *Nucleic Acids Research*, 29, e45.
[doi:10.1093/nar/29.9.e45](https://doi.org/10.1093/nar/29.9.e45)

Vandesompele, J. et al. (2002). Accurate normalization of real-time
quantitative RT-PCR data by geometric averaging of multiple internal
control genes. *Genome Biology*, 3, research0034.
[doi:10.1186/gb-2002-3-7-research0034](https://doi.org/10.1186/gb-2002-3-7-research0034)

## Examples

``` r
rq <- bio_qpcr_ddct(qpcr_sim, ct, gene, sample, treatment,
  reference_genes = c("gapA", "rpoD"), control = "control"
)
rq
#> # A bio_assay: qpcr: 48 x 9
#>    sample     treatment gene  ct_target ct_reference delta_ct delta_delta_ct
#>    <chr>      <chr>     <chr>     <dbl>        <dbl>    <dbl>          <dbl>
#>  1 control_1  control   dnaK       25.2         20.4  4.78            0.402 
#>  2 control_2  control   dnaK       23.6         19.7  3.92           -0.457 
#>  3 control_3  control   dnaK       24.6         20.3  4.30           -0.0717
#>  4 control_4  control   dnaK       23.8         19.3  4.50            0.127 
#>  5 heat_42C_1 heat_42C  dnaK       20.6         19.6  1.03           -3.34  
#>  6 heat_42C_2 heat_42C  dnaK       21.9         20.0  1.95           -2.42  
#>  7 heat_42C_3 heat_42C  dnaK       19.9         18.4  1.43           -2.94  
#>  8 heat_42C_4 heat_42C  dnaK       20.5         19.0  1.54           -2.83  
#>  9 heat_46C_1 heat_46C  dnaK       19.5         19.9 -0.420          -4.80  
#> 10 heat_46C_2 heat_46C  dnaK       19.3         19.3 -0.00333        -4.38  
#> # ℹ 38 more rows
#> # ℹ 2 more variables: fold_change <dbl>, log2_fold_change <dbl>
tidy(rq)
#> # A tibble: 12 × 9
#>    gene  treatment     n log2_fold_change std.error fold_change conf.low
#>    <chr> <chr>     <int>            <dbl>     <dbl>       <dbl>    <dbl>
#>  1 dnaK  control       4        -5.13e-16    0.181        1        0.672
#>  2 dnaK  heat_42C      4         2.89e+ 0    0.190        7.39     4.86 
#>  3 dnaK  heat_46C      4         4.60e+ 0    0.182       24.2     16.2  
#>  4 groL  control       4        -4.86e-16    0.140        1        0.735
#>  5 groL  heat_42C      4         2.08e+ 0    0.183        4.22     2.82 
#>  6 groL  heat_46C      4         2.84e+ 0    0.1000       7.15     5.73 
#>  7 ibpA  control       4         2.18e-15    0.268        1.00     0.554
#>  8 ibpA  heat_42C      4         3.82e+ 0    0.184       14.1      9.41 
#>  9 ibpA  heat_46C      4         5.90e+ 0    0.276       59.8     32.5  
#> 10 ompF  control       4         1.25e-15    0.106        1.00     0.791
#> 11 ompF  heat_42C      4        -7.60e- 1    0.248        0.590    0.342
#> 12 ompF  heat_46C      4        -1.76e+ 0    0.287        0.295    0.157
#> # ℹ 2 more variables: conf.high <dbl>, p.value <dbl>
```
