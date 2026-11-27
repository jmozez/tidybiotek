# Sequence tibbles

`bio_seq_tbl()` marks one column of a data frame as holding biological
sequences and records their alphabet. The sequences are cleaned (upper
case, whitespace and gap characters removed) and the alphabet is kept
through `dplyr` verbs, so later functions such as
[`bio_translate()`](https://tidybiotek.in4metrix.dev/reference/bio_translate.md)
and
[`bio_kmer_count()`](https://tidybiotek.in4metrix.dev/reference/bio_kmer_count.md)
do not need the column named again. `new_bio_seq_tbl()` is the low-level
constructor for developers.

## Usage

``` r
bio_seq_tbl(.data, seq, alphabet = NULL)

new_bio_seq_tbl(
  x = tibble::tibble(seq = character()),
  seq_col = "seq",
  alphabet = "dna"
)

is_bio_seq_tbl(x)

bio_alphabet(x)
```

## Arguments

- .data:

  A data frame or tibble.

- seq:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  The column holding sequences.

- alphabet:

  One of `"dna"`, `"rna"` or `"protein"`. `NULL` (the default) guesses
  it from the sequences.

- x:

  For `new_bio_seq_tbl()`, a data frame. For `is_bio_seq_tbl()` and
  `bio_alphabet()`, any object.

- seq_col:

  For `new_bio_seq_tbl()`, the name of the sequence column.

## Value

`bio_seq_tbl()` and `new_bio_seq_tbl()` return a `<bio_seq_tbl>`, a
tibble subclass. `is_bio_seq_tbl()` returns a single logical and
`bio_alphabet()` returns the alphabet string, or `NULL`.

## Examples

``` r
seqs <- bio_seq_tbl(ecoli_cds, sequence)
seqs
#> # A bio_seq_tbl: 76 sequences (DNA) in `sequence`
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
bio_alphabet(seqs)
#> [1] "dna"

# The class survives dplyr verbs that keep the sequence column
seqs |>
  dplyr::filter(class == "ribosomal") |>
  is_bio_seq_tbl()
#> [1] TRUE
```
