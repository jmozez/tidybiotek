# tidybiotek: Tidy Tools for Applied Biotechnology Bioinformatics

tidybiotek brings applied biotechnology bioinformatics into the
tidyverse and tidymodels. Every user-facing function takes a data frame
first, uses tidy evaluation for column arguments and returns a tibble or
a tibble subclass, so analyses read as pipelines.

## Details

The package is organised into modules:

- **Sequences**:
  [`bio_seq_tbl()`](https://tidybiotek.in4metrix.dev/reference/bio_seq_tbl.md),
  [`bio_translate()`](https://tidybiotek.in4metrix.dev/reference/bio_translate.md),
  [`bio_gc_content()`](https://tidybiotek.in4metrix.dev/reference/bio_gc_content.md),
  [`bio_codon_weights()`](https://tidybiotek.in4metrix.dev/reference/bio_codon_weights.md),
  [`bio_cai()`](https://tidybiotek.in4metrix.dev/reference/bio_cai.md),
  [`bio_codon_optimise()`](https://tidybiotek.in4metrix.dev/reference/bio_codon_optimise.md),
  [`bio_primer_tm()`](https://tidybiotek.in4metrix.dev/reference/bio_primer_tm.md),
  [`bio_kmer_count()`](https://tidybiotek.in4metrix.dev/reference/bio_kmer_count.md)
  and
  [`bio_aa_props()`](https://tidybiotek.in4metrix.dev/reference/bio_aa_props.md).

- **Assays**:
  [`bio_qpcr_ddct()`](https://tidybiotek.in4metrix.dev/reference/bio_qpcr_ddct.md),
  [`bio_dose_response()`](https://tidybiotek.in4metrix.dev/reference/bio_dose_response.md)
  and
  [`bio_elisa_quant()`](https://tidybiotek.in4metrix.dev/reference/bio_elisa_quant.md).

- **Bioprocess**:
  [`bio_enzyme_kinetics()`](https://tidybiotek.in4metrix.dev/reference/bio_enzyme_kinetics.md)
  and
  [`bio_monod_growth()`](https://tidybiotek.in4metrix.dev/reference/bio_monod_growth.md).

- **Omics**:
  [`bio_diff_expr()`](https://tidybiotek.in4metrix.dev/reference/bio_diff_expr.md)
  and
  [`bio_enrich()`](https://tidybiotek.in4metrix.dev/reference/bio_enrich.md).

- **tidymodels**:
  [`step_bio_kmers()`](https://tidybiotek.in4metrix.dev/reference/step_bio_kmers.md),
  [`step_bio_aa_props()`](https://tidybiotek.in4metrix.dev/reference/step_bio_aa_props.md),
  [`step_bio_vst_norm()`](https://tidybiotek.in4metrix.dev/reference/step_bio_vst_norm.md),
  [`bio_homology_cv()`](https://tidybiotek.in4metrix.dev/reference/bio_homology_cv.md),
  [`bio_batch_cv()`](https://tidybiotek.in4metrix.dev/reference/bio_batch_cv.md),
  [`enrichment_factor()`](https://tidybiotek.in4metrix.dev/reference/enrichment_factor.md),
  [`bedroc()`](https://tidybiotek.in4metrix.dev/reference/bedroc.md),
  [`assay_cv_error()`](https://tidybiotek.in4metrix.dev/reference/assay_cv_error.md)
  and
  [`bio_four_pl()`](https://tidybiotek.in4metrix.dev/reference/bio_four_pl.md).

- **Reporting**:
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`glance()`](https://generics.r-lib.org/reference/glance.html),
  [`augment()`](https://generics.r-lib.org/reference/augment.html) and
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  methods.

All numerical methods (codon adaptation, nearest-neighbour
thermodynamics, Levenberg-Marquardt four-parameter logistic fitting, the
exact integrated solution of Monod kinetics and empirical Bayes variance
moderation) are written in base R and `stats`.

## See also

Useful links:

- <https://github.com/jmozez/tidybiotek>

- Report bugs at <https://github.com/jmozez/tidybiotek/issues>

## Author

**Maintainer**: Jimmy Moses <jimmy.moses@pnguot.ac.pg>
([ORCID](https://orcid.org/0000-0003-3342-2384)) \[copyright holder\]

Authors:

- Jimmy Moses <jimmy.moses@pnguot.ac.pg>
  ([ORCID](https://orcid.org/0000-0003-3342-2384)) \[copyright holder\]
