# Escherichia coli K-12 coding sequences

Coding sequences and protein translations of the 56 ribosomal protein
genes and 20 regulatory or metabolic genes of *Escherichia coli* K-12
substr. MG1655. Ribosomal protein genes are among the most highly
expressed in the cell, so they are the classic reference set for codon
adaptation (Sharp and Li, 1987) and form the default reference of
[`bio_cai()`](https://tidybiotek.in4metrix.dev/reference/bio_cai.md).

## Usage

``` r
ecoli_cds
```

## Format

A tibble with 76 rows and 7 columns:

- gene:

  Gene name.

- locus_tag:

  Blattner locus tag (b number).

- protein_id:

  RefSeq protein accession.

- product:

  Protein product name.

- class:

  `"ribosomal"` for ribosomal protein genes, otherwise `"other"`.

- sequence:

  Coding sequence (DNA, start to stop codon).

- protein:

  RefSeq protein sequence.

## Source

NCBI RefSeq NC_000913.3, *Escherichia coli* str. K-12 substr. MG1655
complete genome, downloaded with the NCBI E-utilities as CDS nucleotide
and protein FASTA (`data-raw/make_data.R`).
<https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3>

## References

Blattner, F. R. et al. (1997). The complete genome sequence of
*Escherichia coli* K-12. *Science*, 277, 1453-1462.
[doi:10.1126/science.277.5331.1453](https://doi.org/10.1126/science.277.5331.1453)

## Examples

``` r
ecoli_cds
#> # A tibble: 76 × 7
#>    gene  locus_tag protein_id  product                    class sequence protein
#>    <chr> <chr>     <chr>       <chr>                      <chr> <chr>    <chr>  
#>  1 rplA  b3984     NP_418411.1 50S ribosomal subunit pro… ribo… ATGGCTA… MAKLTK…
#>  2 rplB  b3317     NP_417776.1 50S ribosomal subunit pro… ribo… ATGGCAG… MAVVKC…
#>  3 rplC  b3320     NP_417779.1 50S ribosomal subunit pro… ribo… ATGATTG… MIGLVG…
#>  4 rplD  b3319     NP_417778.1 50S ribosomal subunit pro… ribo… ATGGAAT… MELVLK…
#>  5 rplE  b3308     NP_417767.1 50S ribosomal subunit pro… ribo… ATGGCGA… MAKLHD…
#>  6 rplF  b3305     NP_417764.1 50S ribosomal subunit pro… ribo… ATGTCTC… MSRVAK…
#>  7 rplI  b4203     NP_418624.1 50S ribosomal subunit pro… ribo… ATGCAAG… MQVILL…
#>  8 rplJ  b3985     NP_418412.1 50S ribosomal subunit pro… ribo… ATGGCTT… MALNLQ…
#>  9 rplK  b3983     NP_418410.1 50S ribosomal subunit pro… ribo… ATGGCTA… MAKKVQ…
#> 10 rplL  b3986     NP_418413.1 50S ribosomal subunit pro… ribo… ATGTCTA… MSITKD…
#> # ℹ 66 more rows
```
