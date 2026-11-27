# tidybiotek (development version)

* New `bio_align()` for pairwise global (Needleman-Wunsch) and local
  (Smith-Waterman) alignment of DNA, RNA or protein sequences with affine
  gap penalties (Gotoh 1982) and the BLOSUM62 matrix. Scores match
  `pwalign::pairwiseAlignment()`.
* New `bio_growth_curve()` fits the modified Gompertz, logistic and
  Richards growth models (Zwietering et al. 1990) and the Baranyi and
  Roberts (1994) model to log population size, reporting the lag time,
  maximum specific growth rate and doubling time, with `tidy()`,
  `glance()`, `augment()` and `autoplot()` methods.
* New `bio_qpcr_efficiency()` fits qPCR standard curves and reports the
  amplification efficiency with a confidence interval, the linearity and
  the dynamic range, flagged against the MIQE guidelines. The
  amplification factors feed into `bio_qpcr_ddct()`.

# tidybiotek 0.1.0

* First version, with modules for sequences, assays, bioprocess kinetics
  and omics, tibble subclasses built with vctrs, recipes steps, rsample
  resampling schemes, yardstick metrics, a parsnip 4PL model, and
  `tidy()`, `glance()`, `augment()` and `autoplot()` methods.
