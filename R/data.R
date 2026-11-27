#' Escherichia coli K-12 coding sequences
#'
#' Coding sequences and protein translations of the 56 ribosomal protein
#' genes and 20 regulatory or metabolic genes of *Escherichia coli* K-12
#' substr. MG1655. Ribosomal protein genes are among the most highly
#' expressed in the cell, so they are the classic reference set for codon
#' adaptation (Sharp and Li, 1987) and form the default reference of
#' [bio_cai()].
#'
#' @format A tibble with 76 rows and 7 columns:
#' \describe{
#'   \item{gene}{Gene name.}
#'   \item{locus_tag}{Blattner locus tag (b number).}
#'   \item{protein_id}{RefSeq protein accession.}
#'   \item{product}{Protein product name.}
#'   \item{class}{`"ribosomal"` for ribosomal protein genes, otherwise `"other"`.}
#'   \item{sequence}{Coding sequence (DNA, start to stop codon).}
#'   \item{protein}{RefSeq protein sequence.}
#' }
#' @source NCBI RefSeq NC_000913.3, *Escherichia coli* str. K-12 substr.
#'   MG1655 complete genome, downloaded with the NCBI E-utilities as CDS
#'   nucleotide and protein FASTA (`data-raw/make_data.R`).
#'   <https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3>
#' @references Blattner, F. R. et al. (1997). The complete genome sequence
#'   of *Escherichia coli* K-12. *Science*, 277, 1453-1462.
#'   \doi{10.1126/science.277.5331.1453}
#' @examples
#' ecoli_cds
"ecoli_cds"

#' DNase ELISA standard curves
#'
#' Optical density from an ELISA for recombinant DNase in rat serum, with
#' eleven runs of a serial dilution standard curve. A tidy copy of
#' [datasets::DNase].
#'
#' @format A tibble with 176 rows and 3 columns:
#' \describe{
#'   \item{run}{Assay run (1 to 11).}
#'   \item{conc}{Known protein concentration (ng/ml).}
#'   \item{density}{Measured optical density.}
#' }
#' @source R Core Team, [datasets::DNase]. Original data from Davidian, M.
#'   and Giltinan, D. M. (1995). *Nonlinear Models for Repeated Measurement
#'   Data*. Chapman and Hall. \doi{10.1201/9780203745502}
#' @examples
#' dnase_elisa
"dnase_elisa"

#' Puromycin enzyme kinetics
#'
#' Initial reaction velocity against substrate concentration for an enzyme
#' in cells treated or not treated with puromycin. A tidy copy of
#' [datasets::Puromycin], a standard example of Michaelis-Menten kinetics.
#'
#' @format A tibble with 23 rows and 3 columns:
#' \describe{
#'   \item{state}{`"treated"` or `"untreated"`.}
#'   \item{substrate}{Substrate concentration (ppm).}
#'   \item{rate}{Initial reaction rate (counts/min/min).}
#' }
#' @source R Core Team, [datasets::Puromycin]. Original data from Treloar,
#'   M. A. (1974). *Effects of Puromycin on Galactosyltransferase in Golgi
#'   Membranes*. MSc thesis, University of Toronto; analysed in Bates, D. M.
#'   and Watts, D. G. (1988). *Nonlinear Regression Analysis and Its
#'   Applications*. Wiley. \doi{10.1002/9780470316757}
#' @examples
#' puromycin_kinetics
"puromycin_kinetics"

#' Simulated qPCR experiment
#'
#' **Simulated data.** Ct values for an *E. coli* heat-shock experiment
#' with two reference genes (`gapA`, `rpoD`) and four targets, three
#' treatments, four biological samples per treatment and three technical
#' replicates per sample, run on two plates. Ct values were generated in
#' `data-raw/make_data.R` from known log2 fold changes (control 0; heat
#' 42 C: dnaK 3, groL 2, ibpA 4, ompF -1; heat 46 C: dnaK 4.5, groL 3,
#' ibpA 6, ompF -2) with sample-level input variation and technical noise.
#' Amplification efficiency is 2 for every gene.
#'
#' @format A tibble with 432 rows and 6 columns:
#' \describe{
#'   \item{plate}{Plate (`"plate_1"` or `"plate_2"`).}
#'   \item{sample}{Biological sample.}
#'   \item{treatment}{`"control"`, `"heat_42C"` or `"heat_46C"`.}
#'   \item{gene}{Gene name.}
#'   \item{replicate}{Technical replicate (1 to 3).}
#'   \item{ct}{Quantification cycle.}
#' }
#' @source Simulated with a fixed seed in `data-raw/make_data.R`.
#' @examples
#' qpcr_sim
"qpcr_sim"

#' Simulated dose-response screen
#'
#' **Simulated data.** Responses (percent of control) of three compounds
#' over a ten-point, three-fold dilution series from 100 uM plus a vehicle
#' control, in triplicate. Generated from four-parameter logistic curves
#' with known values (inhibitor_a IC50 0.05 uM, Hill -1.1; inhibitor_b IC50
#' 2.5 uM, Hill -1.6; activator_c EC50 1.2 uM, Hill 1.3) plus Gaussian
#' noise (SD 3.5). The compound names are placeholders, not real drugs.
#'
#' @format A tibble with 99 rows and 4 columns:
#' \describe{
#'   \item{compound}{Compound name.}
#'   \item{dose}{Dose (uM); 0 is the vehicle control.}
#'   \item{replicate}{Replicate (1 to 3).}
#'   \item{response}{Response (percent of control).}
#' }
#' @source Simulated with a fixed seed in `data-raw/make_data.R`.
#' @examples
#' dose_response_sim
"dose_response_sim"

#' Simulated batch fermentations
#'
#' **Simulated data.** Biomass and glucose time courses from three batch
#' bioreactors with different starting glucose (10, 20 and 30 g/l),
#' generated from the exact solution of the Monod model with
#' \eqn{\mu_{max}}{mu_max} = 0.55 per hour, Ks = 0.4 g/l, yield 0.48 g/g and
#' starting biomass 0.08 g/l, with 4 percent multiplicative noise on
#' biomass and 0.25 g/l additive noise on glucose.
#'
#' @format A tibble with 45 rows and 4 columns:
#' \describe{
#'   \item{reactor}{Reactor (`"R1"`, `"R2"`, `"R3"`).}
#'   \item{time_h}{Time (hours).}
#'   \item{biomass_g_l}{Dry cell weight (g/l).}
#'   \item{glucose_g_l}{Residual glucose (g/l).}
#' }
#' @source Simulated with a fixed seed in `data-raw/make_data.R`.
#' @examples
#' fermentation_sim
"fermentation_sim"

#' Simulated RNA-seq counts and gene sets
#'
#' **Simulated data.** `expression_sim` holds negative binomial read counts
#' for 1000 genes in 8 samples (4 control, 4 stress) with random library
#' sizes. Genes `gene_0001` to `gene_0050` are truly up-regulated under
#' stress (log2 fold change 1 to 3), `gene_0051` to `gene_0075` truly
#' down-regulated and the rest unchanged. `gene_sets_sim` holds 30 gene
#' sets in long format: `stress_response` (30 up-regulated genes and 10
#' unchanged ones), `repressed_metabolism` (15 down-regulated and 15
#' unchanged) and 28 random background sets.
#'
#' @format `expression_sim` is a tibble with 8000 rows and 4 columns:
#' \describe{
#'   \item{gene}{Gene identifier.}
#'   \item{sample}{Sample identifier.}
#'   \item{condition}{`"control"` or `"stress"`.}
#'   \item{count}{Read count.}
#' }
#' `gene_sets_sim` is a tibble with columns `gene_set` and `feature`.
#' @source Simulated with a fixed seed in `data-raw/make_data.R`.
#' @examples
#' expression_sim
#' gene_sets_sim
"expression_sim"

#' @rdname expression_sim
"gene_sets_sim"
