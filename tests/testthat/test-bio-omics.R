de_fit <- function(...) {
  bio_diff_expr(expression_sim, gene, sample, count, condition,
    treatment = "stress", reference = "control", transform = "log2_cpm", ...
  )
}

test_that("bio_diff_expr() finds the simulated differentially expressed genes", {
  de <- de_fit()
  expect_s3_class(de, "bio_omics")
  expect_identical(nrow(de), 1000L)
  sig <- de$feature[de$p_adj < 0.05]
  truth <- sprintf("gene_%04d", 1:75)
  expect_gt(mean(truth %in% sig), 0.7)
  expect_lt(mean(!sig %in% truth), 0.1)
  up <- de$log2fc[de$feature %in% sprintf("gene_%04d", 1:50)]
  expect_gt(mean(up > 0), 0.95)
  expect_true(is.finite(attr(de, "bio_meta")$df_prior))
  expect_false(is.unsorted(de$p_value))
})

test_that("moderation reduces to the pooled t-test without a prior", {
  plain <- de_fit(robust_prior = FALSE)
  g <- "gene_0001"
  lc <- expression_sim
  lib <- tapply(lc$count, lc$sample, sum)
  lc$y <- log2((lc$count + 0.5) / (lib[lc$sample] + 1) * 1e6)
  sub <- lc[lc$gene == g, ]
  tt <- stats::t.test(y ~ factor(condition, levels = c("stress", "control")),
    data = sub, var.equal = TRUE
  )
  row <- plain[plain$feature == g, ]
  expect_equal(row$t, unname(tt$statistic), tolerance = 1e-8)
  expect_equal(row$p_value, tt$p.value, tolerance = 1e-8)
})

test_that("trigamma_inverse() inverts trigamma", {
  x <- c(0.01, 0.5, 2, 10)
  expect_equal(trigamma(tidybiotek:::trigamma_inverse(x)), x, tolerance = 1e-6)
})

test_that("bio_diff_expr() validates its input", {
  dup <- rbind(expression_sim[1, ], expression_sim)
  expect_error(bio_diff_expr(dup, gene, sample, count, condition), "one value")
  expect_error(
    bio_diff_expr(expression_sim, gene, sample, count, condition, treatment = "x"),
    "two different"
  )
  small <- expression_sim[expression_sim$sample %in% c("S01", "S05", "S06"), ]
  expect_error(bio_diff_expr(small, gene, sample, count, condition), "at least 2")
  neg <- expression_sim
  neg$count[1] <- -1
  expect_error(bio_diff_expr(neg, gene, sample, count, condition, transform = "log2_cpm"),
               "non-negative")
  miss <- expression_sim
  miss$count[1] <- NA
  expect_warning(bio_diff_expr(miss, gene, sample, count, condition, transform = "log2"),
                 "Dropping")
})

test_that("bio_enrich() detects the planted gene set", {
  de <- de_fit()
  ora <- bio_enrich(de, feature, p_adj < 0.05 & log2fc > 0, gene_sets_sim)
  expect_s3_class(ora, "bio_omics")
  expect_identical(ora$gene_set[[1]], "stress_response")
  expect_lt(ora$p_adj[[1]], 1e-10)
  expect_gt(ora$fold_enrichment[[1]], 5)
  expect_true(all(ora$overlap <= ora$set_size))
  expect_type(ora$overlap_features, "list")
})

test_that("bio_enrich() matches a hand hypergeometric test", {
  d <- tibble::tibble(g = paste0("g", 1:100), hit = rep(c(TRUE, FALSE), c(10, 90)))
  sets <- list(a = paste0("g", 1:20), b = paste0("g", 50:80))
  ora <- bio_enrich(d, g, hit, sets)
  a <- ora[ora$gene_set == "a", ]
  expect_identical(a$overlap, 10L)
  expect_equal(a$p_value, stats::phyper(9, 20, 80, 10, lower.tail = FALSE))
  expect_equal(a$expected, 2)
  df_sets <- tibble::tibble(gene_set = "a", feature = paste0("g", 1:20))
  expect_equal(bio_enrich(d, g, hit, df_sets)$p_value, a$p_value)
  expect_warning(bio_enrich(d, g, hit, sets, min_size = 50), "No gene sets")
  expect_warning(bio_enrich(d, g, rep(FALSE, 100), sets), "No features")
  expect_error(bio_enrich(d, g, 1, sets), "logical")
  expect_error(bio_enrich(d, g, hit, list(1:3)), "named list")
  expect_error(bio_enrich(d, g, hit, data.frame(x = 1)), "gene_set")
})
