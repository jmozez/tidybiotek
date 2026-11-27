# Builds the datasets shipped with tidybiotek.
#
# Run from the package root with: Rscript data-raw/make_data.R
#
# Public data:
#   ecoli_cds           NCBI RefSeq NC_000913.3 (Escherichia coli K-12 MG1655)
#   dnase_elisa         R datasets::DNase (Davidian and Giltinan, 1995)
#   puromycin_kinetics  R datasets::Puromycin (Treloar, 1974; Bates and Watts, 1988)
#
# Simulated data (clearly labelled, generated below with fixed seeds):
#   qpcr_sim, dose_response_sim, fermentation_sim, expression_sim, gene_sets_sim

stopifnot(file.exists("DESCRIPTION"))
dir.create("data", showWarnings = FALSE)

save_data <- function(obj, name) {
  assign(name, obj)
  save(list = name, file = file.path("data", paste0(name, ".rda")), compress = "xz",
       version = 2)
}

# ---- ecoli_cds: E. coli K-12 MG1655 coding sequences from NCBI -------------

eutils <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
fetch_fasta <- function(rettype) {
  url <- paste0(eutils, "?db=nuccore&id=NC_000913.3&rettype=", rettype, "&retmode=text")
  dest <- tempfile(fileext = ".fa")
  utils::download.file(url, dest, quiet = TRUE, mode = "wb")
  lines <- readLines(dest)
  head_idx <- grep("^>", lines)
  ends <- c(head_idx[-1] - 1L, length(lines))
  data.frame(
    header = lines[head_idx],
    seq = vapply(seq_along(head_idx), function(i) {
      paste(lines[(head_idx[i] + 1L):ends[i]], collapse = "")
    }, character(1)),
    stringsAsFactors = FALSE
  )
}
tag <- function(header, key) {
  m <- regmatches(header, regexpr(paste0("\\[", key, "=[^]]*\\]"), header))
  out <- rep(NA_character_, length(header))
  has <- grepl(paste0("\\[", key, "="), header)
  out[has] <- sub(paste0("^\\[", key, "=(.*)\\]$"), "\\1", m)
  out
}

cds <- fetch_fasta("fasta_cds_na")
aa <- fetch_fasta("fasta_cds_aa")
cds$gene <- tag(cds$header, "gene")
cds$locus_tag <- tag(cds$header, "locus_tag")
cds$protein_id <- tag(cds$header, "protein_id")
cds$product <- tag(cds$header, "protein")
aa$protein_id <- tag(aa$header, "protein_id")

ribosomal <- grepl("^(30S|50S) ribosomal subunit protein", cds$product)
other_genes <- c(
  "lacI", "lacZ", "lacY", "lacA", "araC", "trpR", "galR", "lexA", "crp", "fnr",
  "cyaA", "rpoS", "rpoH", "dnaA", "metJ", "purR", "argR", "malT", "gcvA", "ompR"
)
keep <- ribosomal | cds$gene %in% other_genes
ecoli_cds <- cds[keep, ]
ecoli_cds <- ecoli_cds[!duplicated(ecoli_cds$gene), ]
ecoli_cds$class <- ifelse(grepl("ribosomal subunit protein", ecoli_cds$product),
                          "ribosomal", "other")
ecoli_cds$protein <- aa$seq[match(ecoli_cds$protein_id, aa$protein_id)]
ecoli_cds <- ecoli_cds[order(ecoli_cds$class != "ribosomal", ecoli_cds$gene), ]
ecoli_cds <- tibble::tibble(
  gene = ecoli_cds$gene,
  locus_tag = ecoli_cds$locus_tag,
  protein_id = ecoli_cds$protein_id,
  product = ecoli_cds$product,
  class = ecoli_cds$class,
  sequence = toupper(ecoli_cds$seq),
  protein = ecoli_cds$protein
)
stopifnot(sum(ecoli_cds$class == "ribosomal") >= 50, !anyNA(ecoli_cds$protein))
save_data(ecoli_cds, "ecoli_cds")

# ---- dnase_elisa: R datasets::DNase -----------------------------------------

dnase <- datasets::DNase
dnase_elisa <- tibble::tibble(
  run = as.integer(as.character(dnase$Run)),
  conc = dnase$conc,
  density = dnase$density
)
dnase_elisa <- dnase_elisa[order(dnase_elisa$run, dnase_elisa$conc), ]
save_data(dnase_elisa, "dnase_elisa")

# ---- puromycin_kinetics: R datasets::Puromycin ------------------------------

puro <- datasets::Puromycin
puromycin_kinetics <- tibble::tibble(
  state = as.character(puro$state),
  substrate = puro$conc,
  rate = puro$rate
)
save_data(puromycin_kinetics, "puromycin_kinetics")

# ---- qpcr_sim: SIMULATED qPCR Ct values -------------------------------------

set.seed(20261127)
true_log2fc <- list(
  control = c(dnaK = 0, groL = 0, ibpA = 0, ompF = 0),
  heat_42C = c(dnaK = 3, groL = 2, ibpA = 4, ompF = -1),
  heat_46C = c(dnaK = 4.5, groL = 3, ibpA = 6, ompF = -2)
)
baseline_ct <- c(gapA = 17, rpoD = 22, dnaK = 24, groL = 21, ibpA = 29, ompF = 19)
qpcr_sim <- do.call(rbind, lapply(names(true_log2fc), function(trt) {
  do.call(rbind, lapply(1:4, function(rep_id) {
    smp <- paste0(trt, "_", rep_id)
    input_shift <- stats::rnorm(1, 0, 0.6)
    plate <- if (rep_id <= 2) "plate_1" else "plate_2"
    do.call(rbind, lapply(names(baseline_ct), function(g) {
      lfc <- if (g %in% names(true_log2fc[[trt]])) true_log2fc[[trt]][[g]] else 0
      bio_ct <- baseline_ct[[g]] - lfc + input_shift + stats::rnorm(1, 0, 0.25)
      data.frame(
        plate = plate, sample = smp, treatment = trt, gene = g,
        replicate = 1:3, ct = round(bio_ct + stats::rnorm(3, 0, 0.12), 2),
        stringsAsFactors = FALSE
      )
    }))
  }))
}))
qpcr_sim <- tibble::as_tibble(qpcr_sim)
save_data(qpcr_sim, "qpcr_sim")

# ---- dose_response_sim: SIMULATED 4PL dose-response data --------------------

set.seed(42)
four_pl <- function(x, bottom, top, ec50, hill) bottom + (top - bottom) / (1 + (ec50 / x)^hill)
compounds <- data.frame(
  compound = c("inhibitor_a", "inhibitor_b", "activator_c"),
  bottom = c(4, 10, 8),
  top = c(100, 98, 85),
  ec50 = c(0.05, 2.5, 1.2),
  hill = c(-1.1, -1.6, 1.3),
  stringsAsFactors = FALSE
)
doses <- c(0, 100 / 3^(0:9))
dose_response_sim <- do.call(rbind, lapply(seq_len(nrow(compounds)), function(i) {
  p <- compounds[i, ]
  grid <- expand.grid(dose = doses, replicate = 1:3)
  mu <- ifelse(grid$dose == 0,
               if (p$hill < 0) p$top else p$bottom,
               four_pl(grid$dose, p$bottom, p$top, p$ec50, p$hill))
  data.frame(
    compound = p$compound, dose = grid$dose, replicate = grid$replicate,
    response = round(mu + stats::rnorm(nrow(grid), 0, 3.5), 2),
    stringsAsFactors = FALSE
  )
}))
dose_response_sim <- tibble::as_tibble(dose_response_sim[order(dose_response_sim$compound,
                                                               dose_response_sim$dose), ])
save_data(dose_response_sim, "dose_response_sim")

# ---- fermentation_sim: SIMULATED Monod batch fermentations ------------------

set.seed(7)
monod_exact <- function(times, mu_max, ks, yield, x0, s0) {
  cc <- x0 + yield * s0
  a <- ks * yield / cc
  elapsed <- function(x) ((1 + a) * log(x / x0) - a * log((cc - x) / (cc - x0))) / mu_max
  vapply(times, function(t) {
    hi <- cc * (1 - 1e-12)
    if (t <= 0) return(x0)
    if (elapsed(hi) <= t) return(hi)
    stats::uniroot(function(x) elapsed(x) - t, c(x0, hi), tol = 1e-12)$root
  }, numeric(1))
}
reactors <- data.frame(reactor = c("R1", "R2", "R3"), s0 = c(10, 20, 30))
fermentation_sim <- do.call(rbind, lapply(seq_len(nrow(reactors)), function(i) {
  times <- seq(0, 14, by = 1)
  x <- monod_exact(times, mu_max = 0.55, ks = 0.4, yield = 0.48, x0 = 0.08,
                   s0 = reactors$s0[i])
  s <- pmax(reactors$s0[i] - (x - 0.08) / 0.48, 0)
  data.frame(
    reactor = reactors$reactor[i], time_h = times,
    biomass_g_l = round(x * exp(stats::rnorm(length(x), 0, 0.04)), 4),
    glucose_g_l = round(pmax(s + stats::rnorm(length(s), 0, 0.25), 0), 3),
    stringsAsFactors = FALSE
  )
}))
fermentation_sim <- tibble::as_tibble(fermentation_sim)
save_data(fermentation_sim, "fermentation_sim")

# ---- expression_sim and gene_sets_sim: SIMULATED RNA-seq counts --------------

set.seed(2026)
n_genes <- 1000
genes <- sprintf("gene_%04d", seq_len(n_genes))
log2fc <- c(stats::runif(50, 1, 3), -stats::runif(25, 1, 3), rep(0, n_genes - 75))
base_mean <- exp(stats::rnorm(n_genes, log(200), 1.3))
dispersion <- 0.05 + 2 / base_mean
samples <- data.frame(
  sample = sprintf("S%02d", 1:8),
  condition = rep(c("control", "stress"), each = 4),
  size_factor = stats::runif(8, 0.7, 1.4),
  stringsAsFactors = FALSE
)
expression_sim <- do.call(rbind, lapply(seq_len(nrow(samples)), function(j) {
  mu <- base_mean * samples$size_factor[j] *
    if (samples$condition[j] == "stress") 2^log2fc else 1
  data.frame(
    gene = genes, sample = samples$sample[j], condition = samples$condition[j],
    count = stats::rnbinom(n_genes, mu = mu, size = 1 / dispersion),
    stringsAsFactors = FALSE
  )
}))
expression_sim <- tibble::as_tibble(expression_sim)
save_data(expression_sim, "expression_sim")

up <- genes[1:50]
down <- genes[51:75]
null_genes <- genes[76:n_genes]
set_list <- list(
  stress_response = c(sample(up, 30), sample(null_genes, 10)),
  repressed_metabolism = c(sample(down, 15), sample(null_genes, 15))
)
for (i in seq_len(28)) {
  set_list[[sprintf("background_set_%02d", i)]] <- sample(genes, sample(15:60, 1))
}
gene_sets_sim <- tibble::tibble(
  gene_set = rep(names(set_list), lengths(set_list)),
  feature = unlist(set_list, use.names = FALSE)
)
save_data(gene_sets_sim, "gene_sets_sim")
