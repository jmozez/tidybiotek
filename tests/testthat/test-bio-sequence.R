test_that("bio_translate() reproduces NCBI protein sequences", {
  x <- bio_translate(ecoli_cds, sequence,
    code = "bacterial", to_stop = TRUE,
    start_as_met = TRUE
  )
  expect_identical(x$translated, x$protein)
})

test_that("bio_translate() handles frames, codes and ambiguity", {
  d <- tibble::tibble(s = "ATGGCCTGATAA")
  expect_identical(bio_translate(d, s)$translated, "MA**")
  expect_identical(bio_translate(d, s, to_stop = TRUE)$translated, "MA")
  expect_identical(bio_translate(d, s, code = "vertebrate_mito")$translated, "MAW*")
  expect_identical(bio_translate(tibble::tibble(s = "ATGNNN"), s)$translated, "MX")
  expect_identical(bio_translate(tibble::tibble(s = "AUGGCC"), s)$translated, "MA")
  rc <- bio_translate(tibble::tibble(s = "TTATCAGGCCAT"), s, frame = -1)
  expect_identical(rc$translated, "MA**")
  expect_identical(bio_translate(tibble::tibble(s = "AATGGCC"), s, frame = 2)$translated, "MA")
  expect_identical(bio_translate(tibble::tibble(s = "AT"), s)$translated, "")
  expect_error(bio_translate(d, s, frame = 4), "frame")
  expect_error(bio_translate(d, s, code = "nope"))
  seqs <- bio_seq_tbl(ecoli_cds[1:2, ], sequence)
  expect_s3_class(bio_translate(seqs), "bio_seq_tbl")
})

test_that("bio_gc_content() is exact on simple sequences", {
  d <- tibble::tibble(s = c("GGCC", "ATAT", "GCAT", "ACGNNN", NA))
  out <- bio_gc_content(d, s)
  expect_equal(out$gc, c(1, 0, 0.5, 2 / 3, NA))
  pos <- bio_gc_content(tibble::tibble(s = "GAAGCCGTT"), s, by_position = TRUE)
  expect_equal(pos$gc1, 1)
  expect_equal(pos$gc2, 1 / 3)
  expect_equal(pos$gc3, 1 / 3)
})

test_that("codon weights follow Sharp and Li", {
  w <- bio_codon_weights(tibble::tibble(s = c("TTTTTTTTC", "AAAAAAAAG")), s)
  expect_named(w, c("codon", "amino_acid", "count", "rscu", "w"))
  expect_identical(nrow(w), 64L)
  phe <- w[w$amino_acid == "F", ]
  expect_equal(phe$w[phe$codon == "TTT"], 1)
  expect_equal(phe$w[phe$codon == "TTC"], 0.5)
  expect_equal(sum(phe$rscu), 2)
  lys <- w[w$amino_acid == "K", ]
  expect_equal(lys$w[lys$codon == "AAG"], 0.5)
})

test_that("bio_cai() is 1 for optimal genes and higher for ribosomal genes", {
  w <- bio_codon_weights(dplyr::filter(ecoli_cds, class == "ribosomal"), sequence)
  best <- tibble::tibble(s = "ATGCTGAAAGGT")
  ref <- bio_cai(best, s, reference = w)$cai
  expect_true(ref > 0 && ref <= 1)
  out <- bio_cai(ecoli_cds, sequence)
  expect_true(all(out$cai > 0 & out$cai <= 1))
  expect_gt(mean(out$cai[out$class == "ribosomal"]), mean(out$cai[out$class == "other"]))
  from_seqs <- bio_cai(best, s, reference = ecoli_cds$sequence[1:10])
  expect_true(is.numeric(from_seqs$cai))
  expect_true(is.na(bio_cai(tibble::tibble(s = "ATGTGG"), s)$cai))
  expect_error(bio_cai(best, s, reference = 1:3), "reference")
})

test_that("bio_codon_optimise() round-trips and maximises CAI", {
  d <- tibble::tibble(p = c("MKTAYIAKQR", "MSEQNNTEMTFQIQR"))
  out <- bio_codon_optimise(d, p)
  back <- bio_translate(out, optimised_cds, to_stop = TRUE)
  expect_identical(back$translated, d$p)
  cai <- bio_cai(out, optimised_cds)$cai
  expect_equal(cai, c(1, 1))
  expect_false(grepl("TAA$|TAG$|TGA$", bio_codon_optimise(d, p, add_stop = FALSE)$optimised_cds[1]))
  expect_error(bio_codon_optimise(tibble::tibble(p = "MK1"), p), "back-translate")
})

test_that("bio_primer_tm() gives sensible temperatures", {
  d <- tibble::tibble(p = c("GTAAAACGACGGCCAGT", "AAAAAAAAAATTTTTTTTTT", "GCGCGCGCGCGCGCGCGCGC"))
  nn <- bio_primer_tm(d, p)$tm
  expect_true(all(nn > 20 & nn < 100))
  expect_gt(nn[[3]], nn[[2]])
  w <- bio_primer_tm(d, p, method = "wallace")$tm
  expect_equal(w, c(52, 40, 80))
  b <- bio_primer_tm(d, p, method = "basic")$tm
  expect_equal(b[[3]], 64.9 + 41 * (20 - 16.4) / 20)
  more_salt <- bio_primer_tm(d, p, na_conc = 200)$tm
  expect_true(all(more_salt > nn))
  expect_warning(out <- bio_primer_tm(tibble::tibble(p = "ACGTN"), p), "ambiguous")
  expect_true(is.na(out$tm))
  expect_error(bio_primer_tm(d, p, na_conc = -1), "positive")
})

test_that("bio_primer_tm() matches a hand calculation", {
  s <- "ACGTACGTAC"
  lookup <- tidybiotek:::nn_lookup()
  pairs <- substring(s, 1:9, 2:10)
  dh <- sum(lookup$dh[pairs]) + 2.3 + 0.1
  ds <- sum(lookup$ds[pairs]) + 4.1 - 2.8 + 0.368 * 9 * log(0.05)
  expected <- dh * 1000 / (ds + 1.987 * log(50e-9 / 4)) - 273.15
  expect_equal(bio_primer_tm(tibble::tibble(p = s), p)$tm, expected)
})

test_that("bio_kmer_count() counts overlapping k-mers", {
  d <- tibble::tibble(id = 1:2, s = c("AAAA", "ACGT"))
  out <- bio_kmer_count(d, s, k = 2)
  expect_identical(sum(grepl("^kmer_", names(out))), 16L)
  expect_equal(out$kmer_AA, c(3, 0))
  expect_equal(out$kmer_CG, c(0, 1))
  freq <- bio_kmer_count(d, s, k = 2, normalise = TRUE)
  expect_equal(freq$kmer_AA[[1]], 1)
  voc <- bio_kmer_count(d, s, k = 1, vocabulary = c("A", "T"), prefix = "n_")
  expect_named(voc, c("id", "s", "n_A", "n_T"))
  prot <- bio_kmer_count(tibble::tibble(s = c("MKV", "MKK")), s, k = 2)
  expect_named(prot, c("s", "kmer_KK", "kmer_KV", "kmer_MK"))
  expect_error(bio_kmer_count(d, s, k = 10, vocabulary = character()), "No k-mers")
  expect_error(bio_kmer_count(out, s, k = 2), "already exist")
})

test_that("bio_aa_props() computes known properties", {
  out <- bio_aa_props(tibble::tibble(p = c("GGG", "KKKK", "DDDD*")), p)
  expect_equal(out$length, c(3, 4, 4))
  expect_equal(out$mol_weight[[1]], 3 * 57.0519 + 18.01524)
  expect_equal(out$gravy[[1]], -0.4)
  expect_gt(out$isoelectric_point[[2]], 9)
  expect_lt(out$isoelectric_point[[3]], 4)
  expect_gt(out$charge_ph7[[2]], 0)
  expect_equal(out$frac_basic[[2]], 1)
  expect_equal(out$frac_acidic[[3]], 1)
  na <- bio_aa_props(tibble::tibble(p = c(NA, "XX")), p, prefix = "aa_")
  expect_true(all(is.na(na$aa_length)))
})
