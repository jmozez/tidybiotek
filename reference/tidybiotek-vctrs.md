# vctrs methods for tidybiotek classes

These methods let `<bio_seq_tbl>`, `<bio_assay>` and `<bio_omics>`
objects work with
[`vctrs::vec_rbind()`](https://vctrs.r-lib.org/reference/vec_bind.html),
[`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html)
and other vctrs-based tools. Two objects with the same metadata combine
into the same class; anything else combines into a plain tibble.

## Usage

``` r
# S3 method for class 'bio_seq_tbl'
vec_restore(x, to, ...)

# S3 method for class 'bio_assay'
vec_restore(x, to, ...)

# S3 method for class 'bio_omics'
vec_restore(x, to, ...)

# S3 method for class 'bio_seq_tbl.bio_seq_tbl'
vec_ptype2(x, y, ...)

# S3 method for class 'bio_seq_tbl.tbl_df'
vec_ptype2(x, y, ...)

# S3 method for class 'tbl_df.bio_seq_tbl'
vec_ptype2(x, y, ...)

# S3 method for class 'bio_seq_tbl.data.frame'
vec_ptype2(x, y, ...)

# S3 method for class 'data.frame.bio_seq_tbl'
vec_ptype2(x, y, ...)

# S3 method for class 'bio_seq_tbl.bio_seq_tbl'
vec_cast(x, to, ...)

# S3 method for class 'bio_seq_tbl.tbl_df'
vec_cast(x, to, ...)

# S3 method for class 'tbl_df.bio_seq_tbl'
vec_cast(x, to, ...)

# S3 method for class 'data.frame.bio_seq_tbl'
vec_cast(x, to, ...)

# S3 method for class 'bio_assay.bio_assay'
vec_ptype2(x, y, ...)

# S3 method for class 'bio_assay.tbl_df'
vec_ptype2(x, y, ...)

# S3 method for class 'tbl_df.bio_assay'
vec_ptype2(x, y, ...)

# S3 method for class 'bio_assay.data.frame'
vec_ptype2(x, y, ...)

# S3 method for class 'data.frame.bio_assay'
vec_ptype2(x, y, ...)

# S3 method for class 'bio_assay.bio_assay'
vec_cast(x, to, ...)

# S3 method for class 'bio_assay.tbl_df'
vec_cast(x, to, ...)

# S3 method for class 'tbl_df.bio_assay'
vec_cast(x, to, ...)

# S3 method for class 'data.frame.bio_assay'
vec_cast(x, to, ...)

# S3 method for class 'bio_omics.bio_omics'
vec_ptype2(x, y, ...)

# S3 method for class 'bio_omics.tbl_df'
vec_ptype2(x, y, ...)

# S3 method for class 'tbl_df.bio_omics'
vec_ptype2(x, y, ...)

# S3 method for class 'bio_omics.data.frame'
vec_ptype2(x, y, ...)

# S3 method for class 'data.frame.bio_omics'
vec_ptype2(x, y, ...)

# S3 method for class 'bio_omics.bio_omics'
vec_cast(x, to, ...)

# S3 method for class 'bio_omics.tbl_df'
vec_cast(x, to, ...)

# S3 method for class 'tbl_df.bio_omics'
vec_cast(x, to, ...)

# S3 method for class 'data.frame.bio_omics'
vec_cast(x, to, ...)
```

## Arguments

- x, y, to:

  Objects to restore, combine or cast.

- ...:

  Passed on to the vctrs data frame methods.

## Value

A tibble subclass or a plain tibble.
