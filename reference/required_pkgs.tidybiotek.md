# S3 methods for tracking which additional packages are needed for steps

Recipe-adjacent packages always list themselves as a required package so
that the steps can function properly within parallel processing schemes.

## Usage

``` r
# S3 method for class 'step_bio_kmers'
required_pkgs(x, ...)

# S3 method for class 'step_bio_aa_props'
required_pkgs(x, ...)

# S3 method for class 'step_bio_vst_norm'
required_pkgs(x, ...)
```

## Arguments

- x:

  A recipe step.

- ...:

  Not used.

## Value

A character vector of package names.
