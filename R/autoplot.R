# ggplot2 autoplot() methods.

curve_label <- function(key) {
  if (ncol(key) == 0L) {
    return("all")
  }
  do.call(paste, c(lapply(key, as.character), sep = " / "))
}

curve_grid <- function(fits, predict_fn, log_x = TRUE, n = 200L) {
  rows <- lapply(fits, function(f) {
    x <- f$x %||% f$time
    pos <- x[x > 0]
    grid <- if (log_x) {
      exp(seq(log(min(pos)), log(max(pos)), length.out = n))
    } else {
      seq(min(x), max(x), length.out = n)
    }
    tibble::tibble(.curve = curve_label(f$key), x = grid, y = predict_fn(f, grid))
  })
  vctrs::vec_rbind(!!!rows)
}

point_data <- function(fits) {
  vctrs::vec_rbind(!!!lapply(fits, function(f) {
    tibble::tibble(.curve = curve_label(f$key), x = f$x, y = f$y)
  }))
}

#' Plot tidybiotek results
#'
#' `autoplot()` methods draw the standard figure for each kind of result:
#'
#' * dose-response fits: data and fitted 4PL curves on a log dose axis;
#' * ELISA: the standard curve with back-calculated unknowns;
#' * enzyme kinetics: rates and fitted Michaelis-Menten curves;
#' * Monod growth: biomass and substrate time courses with fitted growth
#'   curves;
#' * primary growth curves: log population size with the fitted Gompertz,
#'   logistic, Richards or Baranyi curves;
#' * qPCR standard curves: quantification cycle against log quantity;
#' * qPCR: fold changes with confidence intervals on a log2 scale;
#' * differential expression: a volcano plot (`type = "volcano"`) or an MA
#'   plot (`type = "ma"`);
#' * enrichment: the most significant gene sets.
#'
#' @param object A `<bio_assay>` or `<bio_omics>` result.
#' @param type For differential expression, `"volcano"` or `"ma"`.
#' @param alpha Adjusted p-value threshold used to highlight features or
#'   gene sets.
#' @param lfc Absolute log2 fold change threshold for highlighting.
#' @param n For enrichment results, the number of gene sets to show.
#' @param ... Not used.
#'
#' @return A [ggplot2::ggplot()] object.
#' @name autoplot_bio
#' @examples
#' fit <- bio_dose_response(dose_response_sim, dose, response, .by = compound)
#' autoplot(fit)
#'
#' growth <- bio_monod_growth(fermentation_sim, time_h, biomass_g_l,
#'   substrate = glucose_g_l, .by = reactor
#' )
#' autoplot(growth)
#'
#' de <- bio_diff_expr(expression_sim, gene, sample, count, condition,
#'   treatment = "stress", reference = "control", transform = "log2_cpm"
#' )
#' autoplot(de)
#' autoplot(de, type = "ma")
NULL

#' @rdname autoplot_bio
#' @exportS3Method ggplot2::autoplot
autoplot.bio_assay <- function(object, ...) {
  meta <- bio_meta(object)
  switch(meta$assay,
    dose_response = plot_dose_response(object, meta),
    elisa = plot_elisa(object, meta),
    enzyme_kinetics = plot_enzyme(object, meta),
    monod = plot_monod(object, meta),
    growth_curve = plot_growth(object, meta),
    qpcr_standard = plot_qpcr_standard(object, meta),
    qpcr = plot_qpcr(object, meta),
    cli::cli_abort("There is no plot for assay type {.val {meta$assay}}.")
  )
}

plot_dose_response <- function(object, meta) {
  fits <- active_fits(object)
  pts <- point_data(fits)
  pts <- pts[pts$x > 0, ]
  curves <- curve_grid(fits, function(f, g) four_pl_value(f$coef, g))
  ggplot2::ggplot(pts, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$.curve)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_line(data = curves, linewidth = 0.8) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = meta$x, y = meta$y, colour = NULL,
                  title = "Four-parameter logistic dose-response") +
    ggplot2::theme_minimal()
}

plot_elisa <- function(object, meta) {
  fits <- meta$fits
  pts <- point_data(fits)
  pts <- pts[pts$x > 0, ]
  curves <- curve_grid(fits, function(f, g) four_pl_value(f$coef, g))
  unknown <- tibble::as_tibble(strip_bio(object))
  unknown <- unknown[!unknown$is_standard & unknown$in_range, ]
  unknown <- tibble::tibble(x = unknown$estimated_conc, y = unknown[[meta$y]])
  p <- ggplot2::ggplot(pts, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_line(data = curves, ggplot2::aes(group = .data$.curve), colour = "grey40") +
    ggplot2::geom_point(ggplot2::aes(shape = "Standards"), size = 2)
  if (nrow(unknown) > 0L) {
    p <- p + ggplot2::geom_point(data = unknown, ggplot2::aes(shape = "Unknowns"),
                                 colour = "firebrick", size = 2)
  }
  p +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = meta$x, y = meta$y, shape = NULL, title = "ELISA standard curve") +
    ggplot2::theme_minimal()
}

plot_enzyme <- function(object, meta) {
  fits <- active_fits(object)
  pts <- point_data(fits)
  curves <- curve_grid(fits, function(f, g) mm_value(f$coef, g), log_x = FALSE)
  ggplot2::ggplot(pts, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$.curve)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(data = curves, linewidth = 0.8) +
    ggplot2::labs(x = meta$x, y = meta$y, colour = NULL,
                  title = "Michaelis-Menten kinetics") +
    ggplot2::theme_minimal()
}

plot_monod <- function(object, meta) {
  fits <- active_fits(object)
  obs <- vctrs::vec_rbind(!!!lapply(fits, function(f) {
    lab <- curve_label(f$key)
    out <- tibble::tibble(.curve = lab, time = f$time, value = f$biomass, variable = "Biomass")
    if (!is.null(f$substrate)) {
      out <- vctrs::vec_rbind(out, tibble::tibble(.curve = lab, time = f$time,
                                                  value = f$substrate, variable = "Substrate"))
    }
    out
  }))
  fitted <- vctrs::vec_rbind(!!!lapply(fits, function(f) {
    grid <- seq(min(f$time), max(f$time), length.out = 200L)
    p <- f$coef
    pred <- monod_solve(grid, p[["mu_max"]], p[["ks"]], p[["yield"]], p[["x0"]], p[["s0"]])
    lab <- curve_label(f$key)
    vctrs::vec_rbind(
      tibble::tibble(.curve = lab, time = grid, value = pred[, "biomass"], variable = "Biomass"),
      tibble::tibble(.curve = lab, time = grid, value = pred[, "substrate"],
                     variable = "Substrate")
    )
  }))
  ggplot2::ggplot(obs, ggplot2::aes(x = .data$time, y = .data$value, colour = .data$.curve)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(data = fitted, linewidth = 0.8) +
    ggplot2::facet_wrap(ggplot2::vars(.data$variable), scales = "free_y") +
    ggplot2::labs(x = meta$x, y = NULL, colour = NULL, title = "Monod batch growth") +
    ggplot2::theme_minimal()
}

plot_qpcr <- function(object, meta) {
  summ <- qpcr_summary(object)
  ggplot2::ggplot(summ, ggplot2::aes(x = .data[[meta$group]], y = .data$fold_change)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high),
                           width = 0.2, na.rm = TRUE) +
    ggplot2::facet_wrap(ggplot2::vars(.data[[meta$gene]])) +
    ggplot2::scale_y_continuous(trans = "log2") +
    ggplot2::labs(x = NULL, y = "Fold change (log2 scale)",
                  title = "Relative expression by qPCR") +
    ggplot2::theme_minimal()
}

#' @rdname autoplot_bio
#' @exportS3Method ggplot2::autoplot
autoplot.bio_omics <- function(object, type = c("volcano", "ma"), alpha = 0.05, lfc = 1,
                               n = 15, ...) {
  meta <- bio_meta(object)
  if (meta$analysis == "enrichment") {
    return(plot_enrichment(object, alpha, n))
  }
  type <- rlang::arg_match(type)
  d <- augment.bio_omics(object, alpha = alpha, lfc = lfc)
  colours <- c(up = "firebrick", down = "steelblue", none = "grey70")
  if (type == "volcano") {
    ggplot2::ggplot(d, ggplot2::aes(x = .data$log2fc, y = -log10(.data$p_value),
                                    colour = .data$.direction)) +
      ggplot2::geom_point(alpha = 0.7, size = 1.2) +
      ggplot2::geom_vline(xintercept = c(-lfc, lfc), linetype = 2, colour = "grey50") +
      ggplot2::scale_colour_manual(values = colours) +
      ggplot2::labs(x = "log2 fold change", y = "-log10 p-value", colour = NULL,
                    title = paste(meta$treatment, "vs", meta$reference)) +
      ggplot2::theme_minimal()
  } else {
    ggplot2::ggplot(d, ggplot2::aes(x = .data$ave_expr, y = .data$log2fc,
                                    colour = .data$.direction)) +
      ggplot2::geom_point(alpha = 0.7, size = 1.2) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
      ggplot2::scale_colour_manual(values = colours) +
      ggplot2::labs(x = "Average expression", y = "log2 fold change", colour = NULL,
                    title = paste(meta$treatment, "vs", meta$reference)) +
      ggplot2::theme_minimal()
  }
}

plot_enrichment <- function(object, alpha, n) {
  d <- tibble::as_tibble(strip_bio(object))
  d <- d[order(d$p_value), ][seq_len(min(n, nrow(d))), ]
  d$gene_set <- factor(d$gene_set, levels = rev(d$gene_set))
  d$significant <- d$p_adj < alpha
  ggplot2::ggplot(d, ggplot2::aes(x = -log10(.data$p_adj), y = .data$gene_set,
                                  fill = .data$significant)) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept = -log10(alpha), linetype = 2, colour = "grey50") +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey70")) +
    ggplot2::labs(x = "-log10 adjusted p-value", y = NULL, fill = "Significant",
                  title = "Over-representation analysis") +
    ggplot2::theme_minimal()
}
