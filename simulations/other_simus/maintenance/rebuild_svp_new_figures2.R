## Rebuild the SVP_NEW paper figures with the current svpChange2 package.
##
## The default mode is deliberately conservative: it consumes the maintained
## result objects already present in simulations/. Set run_full = TRUE to
## regenerate the power/runtime grids after recalibration. In either mode the
## old SVP_NEW_Figures directory is never modified.

library(svpChange2)
library(changepoint)
library(dplyr)
library(ggplot2)

source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "time_common.R"))
source(file.path("simulations", "other_simus", "maintenance",
                 "calibrate_current_svp.R"))

FIGURE_ROOT <- normalizePath(
  file.path("..", "SVP_NEW_Figures2"), mustWork = FALSE
)
FIGURE_PLOTS <- FIGURE_ROOT

fmt_constant <- function(family, value, n = 1000L, true_segments = 8L) {
  family <- tolower(family)
  if (family %in% c("gaussian", "ar1", "ar1focus")) {
    return(sprintf("c = %.3f", value))
  }
  if (family == "wilcoxon") {
    return(sprintf("c = %.3f", value))
  }
  sprintf("alpha = %.4g", value)
}

config_label <- function(family, mode, value) {
  family_key <- tolower(family)
  display <- c(
    gaussian = "Gaussian", ar1 = "AR1", ar1focus = "AR1Focus",
    wilcoxon = "Wilcoxon", mood = "MedianMood"
  )[[family_key]]
  if (is.null(display)) display <- family
  paste0("SVP ", display, " / ", mode, " / ", fmt_constant(family_key, value))
}

make_candidate_table <- function(calibration, family) {
  selected <- calibration$selected_by_mode
  if (is.null(selected) || !nrow(selected)) selected <- calibration$selected
  selected$family <- family
  selected
}

fit_gaussian_refresh <- function(y, selected) {
  n <- length(y)
  pelt <- changepoint::cpt.mean(
    y, method = "PELT", penalty = "Manual", pen.value = 2 * log(n)
  )
  boundaries <- list(
    PELT = normalise_boundaries(changepoint::cpts(pelt), n)
  )
  for (i in seq_len(nrow(selected))) {
    row <- selected[i, ]
    label <- config_label("gaussian", row$mode, row$value)
    boundaries[[label]] <- normalise_boundaries(
      calibration_candidate_fit(
        y, "gaussian", row$mode, row$value, true_segments = 1L
      ), n
    )
  }
  boundaries
}

fit_robust_refresh <- function(y, true_segments, selected) {
  n <- length(y)
  pelt <- changepoint::cpt.mean(
    y, method = "PELT", penalty = "Manual", pen.value = 2 * log(n)
  )
  scale <- stats::mad(y, constant = 1.4826)
  rfpop <- robseg::Rob_seg.std(
    y / scale, loss = "Outlier", lambda = 2 * log(n), lthreshold = 3
  )
  boundaries <- list(
    PELT = normalise_boundaries(changepoint::cpts(pelt), n),
    RFPOP = normalise_boundaries(head(rfpop$t.est, -1L), n)
  )
  for (i in seq_len(nrow(selected))) {
    row <- selected[i, ]
    label_family <- if (row$family == "mood") "MedianMood" else "Wilcoxon"
    label <- paste0(
      "SVP ", label_family, " / ", row$mode, " / ",
      fmt_constant(row$family, row$value)
    )
    boundaries[[label]] <- normalise_boundaries(
      calibration_candidate_fit(
        y, row$family, row$mode, row$value,
        true_segments = true_segments
      ), n
    )
  }
  boundaries
}

run_gaussian_refresh <- function(selected, workers = power_default_workers()) {
  run_power_grid(
    n = 1000L, jump_sizes = seq(0.1, 2, 0.1), reps = 100L,
    simulate_noise = stats::rnorm,
    fit_methods = function(y, true_segments) {
      fit_gaussian_refresh(y, selected)
    },
    tolerance = round(1000 * 0.0025), workers = workers, seed = 123L
  )
}

run_robust_refresh <- function(selected, workers = power_default_workers()) {
  run_power_grid(
    n = 1000L, jump_sizes = seq(0.1, 4, 0.1), reps = 100L,
    simulate_noise = function(n) stats::rt(n, df = 2),
    fit_methods = function(y, true_segments) {
      fit_robust_refresh(y, true_segments, selected)
    },
    tolerance = round(1000 * 0.0025), workers = workers, seed = 123L
  )
}

source_ar1_power_helpers <- function() {
  old_flag <- Sys.getenv("SVP_RUN_SIMULATIONS", unset = NA_character_)
  Sys.setenv(SVP_RUN_SIMULATIONS = "false")
  source(file.path("simulations", "power_ar1", "run_power.R"))
  if (is.na(old_flag)) Sys.unsetenv("SVP_RUN_SIMULATIONS") else
    Sys.setenv(SVP_RUN_SIMULATIONS = old_flag)
}

run_ar1_refresh <- function(selected, workers = power_default_workers()) {
  if (!exists("fit_ar1_methods", mode = "function")) {
    source_ar1_power_helpers()
  }
  row_right <- selected[selected$mode == "right", , drop = FALSE][1L, ]
  row_both <- selected[selected$mode == "both", , drop = FALSE][1L, ]
  if (!nrow(row_right)) row_right <- selected[1L, , drop = FALSE]
  if (!nrow(row_both)) row_both <- selected[nrow(selected), , drop = FALSE]

  run_power_grid(
    n = 600L, jump_sizes = seq(1, 3, 0.1), reps = 100L,
    simulate_noise = function(n) simulate_ar1_noise(n, 0.8),
    fit_methods = function(y, true_segments) {
      fits <- fit_ar1_methods(
        y, rho = 0.8, svp_constant = row_right$value,
        true_true_constant = row_both$value
      )
      names(fits)[names(fits) == "SVP AR1Focus"] <- config_label(
        "AR1Focus", "right", row_right$value
      )
      names(fits)[names(fits) == "SVP AR1Focus multiscale"] <-
        config_label("AR1Focus", "both", row_both$value)
      fits
    },
    tolerance = 5L, workers = workers, seed = 123L
  )
}

summary_metric <- function(results, metric) {
  results |>
    group_by(algorithm, pattern, jump) |>
    summarise(value = mean(.data[[metric]], na.rm = TRUE),
              .groups = "drop")
}

plot_metric_refresh <- function(results, metric, y_label, log_y = FALSE) {
  data <- summary_metric(results, metric)
  data$pattern <- factor(data$pattern, POWER_PATTERNS)
  p <- ggplot(data, aes(jump, value, colour = algorithm, group = algorithm)) +
    geom_line(linewidth = 0.75) + geom_point(size = 1.1) +
    facet_wrap(~pattern, nrow = 2L) +
    labs(x = "Jump size", y = y_label, colour = "Algorithm") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  if (log_y) p <- p + scale_y_continuous(trans = scales::pseudo_log_trans())
  p
}

plot_distribution_refresh <- function(results, n, selected_jump) {
  selected <- results[
    abs(results$jump - selected_jump) < 1e-10 & results$pattern != "none",
  ]
  rows <- lapply(seq_len(nrow(selected)), function(i) {
    cp <- as.integer(selected$changepoints[[i]])
    cp <- cp[cp < n]
    if (!length(cp)) return(NULL)
    data.frame(
      pattern = selected$pattern[i], algorithm = selected$algorithm[i],
      changepoint = cp
    )
  })
  data <- bind_rows(rows)
  ggplot(data, aes(changepoint)) +
    geom_histogram(bins = 100, fill = "steelblue", colour = "black",
                   linewidth = 0.15) +
    facet_grid(pattern ~ algorithm) +
    labs(x = "Sequence position", y = "Frequency of detected change") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "none")
}

plot_scenarios_refresh <- function(n, jump, noise_fun) {
  set_power_seed(999L)
  data <- bind_rows(lapply(POWER_PATTERNS, function(pattern) {
    mu <- generate_power_signal(n, pattern, jump)
    data.frame(
      time = seq_len(n), value = mu + noise_fun(n), mean = mu, pattern = pattern
    )
  }))
  data$pattern <- factor(data$pattern, POWER_PATTERNS)
  ggplot(data, aes(time, value)) +
    geom_point(alpha = 0.2, size = 0.3) +
    geom_line(aes(y = mean), colour = "red", linewidth = 0.7) +
    facet_wrap(~pattern, nrow = 2L) +
    labs(x = "Time (Index)", y = "Value") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "none")
}

write_power_figures <- function(results, prefix, n, selected_jump,
                                scenario_noise, root = FIGURE_ROOT,
                                include_numsegments = TRUE,
                                file_suffix = "") {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  plots <- list(
    signal_scenarios = plot_scenarios_refresh(n, selected_jump, scenario_noise),
    f1 = plot_metric_refresh(results, "F1", "F1"),
    precision = plot_metric_refresh(results, "Precision", "Precision"),
    recall = plot_metric_refresh(results, "Recall", "Recall"),
    correct_number = plot_metric_refresh(
      results, "CorrectNumCP", "Probability of the correct number of changepoints"
    )
  )
  if (include_numsegments) {
    plots$numsegments <- plot_metric_refresh(
      results, "NumSegments", "Number of segments"
    )
  }
  plots$mse <- plot_metric_refresh(results, "MSE", "Signal MSE", log_y = TRUE)
  plots$changepoint_distributions <- plot_distribution_refresh(
    results, n, selected_jump
  )
  file_names <- c(
    "_signal_scenarios.pdf", "_f1_vs_jump.pdf", "_precision_vs_jump.pdf",
    "_recall_vs_jump.pdf", "_correct_number_vs_jump.pdf"
  )
  if (include_numsegments) file_names <- c(file_names, "_numsegments_vs_jump.pdf")
  file_names <- c(file_names, "_mse_vs_jump.pdf", "_changepoint_distributions.pdf")
  suffix_mask <- grepl("signal_scenarios|changepoint_distributions", file_names)
  names(plots) <- paste0(
    prefix,
    ifelse(suffix_mask,
           sub("\\.pdf$", paste0(file_suffix, ".pdf"), file_names),
           file_names)
  )
  for (file in names(plots)) {
    ggsave(file.path(root, file), plots[[file]],
           width = if (grepl("distributions", file)) 14 else 10,
           height = if (grepl("distributions", file)) 10 else 5)
  }
  invisible(plots)
}

plot_runtime_refresh <- function(results, x = c("n", "k"), methods = NULL,
                                 log_x = TRUE, log_y = TRUE) {
  x <- match.arg(x)
  data <- results[results$experiment == if (x == "n") "vary_n" else "vary_k", ]
  if (!is.null(methods)) data <- data[data$method %in% methods, ]
  p <- ggplot(data, aes(.data[[x]], time, colour = method)) +
    geom_smooth(method = "loess", se = TRUE, na.rm = TRUE) +
    labs(x = if (x == "n") "Sequence length" else "Number of true changes",
         y = "Time (s)", colour = "Method") +
    theme_minimal(base_size = 11)
  if (log_x) p <- p + scale_x_log10()
  if (log_y) p <- p + scale_y_log10()
  p
}

load_runtime_or_empty <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE)
}

rename_quick_results <- function(gaussian, robust, ar1) {
  gaussian$algorithm[gaussian$algorithm == "SVP BIC calibrated"] <-
    config_label("gaussian", "right", 1.5)
  gaussian$algorithm[gaussian$algorithm == "SVP BIC"] <-
    config_label("gaussian", "right", 2)
  gaussian$algorithm[gaussian$algorithm == "SVP BIC multiscale"] <-
    config_label("gaussian", "both", 1.8)

  robust$algorithm[robust$algorithm == "SVP MedianMood"] <-
    "SVP MedianMood / right / alpha = 0.01 (K-dependent threshold)"
  robust$algorithm[robust$algorithm == "SVP Wilcoxon"] <-
    config_label("wilcoxon", "right", 1.5)
  robust$algorithm[robust$algorithm == "SVP Wilcoxon multiscale"] <-
    config_label("wilcoxon", "both", 1.75)
  ## The comparison figures use the calibrated RFPOP curve, not the duplicate
  ## historical-paper implementation.
  robust <- robust[robust$algorithm != "RFPOP paper", , drop = FALSE]

  ar1$algorithm[ar1$algorithm == "SVP AR1Focus"] <-
    config_label("AR1Focus", "right", 3.75)
  ar1$algorithm[ar1$algorithm == "SVP AR1Focus multiscale"] <-
    config_label("AR1Focus", "both", 4.5)
  list(gaussian = gaussian, robust = robust, ar1 = ar1)
}

write_runtime_figures <- function(root = FIGURE_ROOT) {
  g <- load_runtime_or_empty(
    file.path("simulations", "time_gaussian", "gaussian_time_results.csv")
  )
  r <- load_runtime_or_empty(
    file.path("simulations", "time_robust", "robust_time_results.csv")
  )
  a <- load_runtime_or_empty(
    file.path("simulations", "time_ar1", "ar1_time_results.csv")
  )
  if (!is.null(g)) {
    g$method[g$method == "SVP BIC multiscale"] <-
      "SVP Gaussian / both / 1.800 log(n)"
    g$method[g$method == "SVP BIC calibrated"] <-
      "SVP Gaussian / right / 1.500 log(n)"
    g$method[g$method == "SVP BIC"] <-
      "SVP Gaussian / right / 2.000 log(n)"
    g$experiment[g$experiment == "vary_k"] <- "vary_k"
    p <- plot_runtime_refresh(g, "n")
    ggsave(file.path(root, "2_gaussian_time_vs_n.pdf"), p, width = 8, height = 5)
    ggsave(file.path(root, "2_gaussian_time_vs_detected.pdf"),
           ggplot(g, aes(detected, time, colour = method)) +
             geom_smooth(method = "loess", se = TRUE, na.rm = TRUE) +
             scale_y_log10() + theme_minimal(base_size = 11) +
             labs(x = "Detected number of segments", y = "Time (s, log scale)",
                  colour = "Method"), width = 8, height = 5)
    ggsave(file.path(root, "2_exp1_time_vs_n.pdf"), p, width = 8, height = 5)
    ggsave(file.path(root, "2_exp2_time_vs_k.pdf"),
           plot_runtime_refresh(g, "k", log_x = FALSE), width = 8, height = 5)
    ggsave(file.path(root, "2_time_vs_detected_changes.pdf"),
           ggplot(g, aes(detected, time, colour = method)) +
             geom_smooth(method = "loess", se = TRUE, na.rm = TRUE) +
             scale_y_log10() + theme_minimal(base_size = 11) +
             labs(x = "Detected number of changes", y = "Time (s, log scale)",
                  colour = "Method"), width = 8, height = 5)
  }
  if (!is.null(r)) {
    r$method[r$method == "SVP BIC multiscale"] <-
      "SVP Gaussian / both / 1.500 log(n)"
    r$method[r$method == "SVP Wilcoxon"] <-
      "SVP Wilcoxon / right / 2.000 log(n)"
    r$method[r$method == "SVP MedianMood"] <-
      "SVP MedianMood / right / 2.000 log(n)"
    ggsave(file.path(root, "2_robust_time_vs_n.pdf"),
           plot_runtime_refresh(r, "n"), width = 8, height = 5)
    ggsave(file.path(root, "2_robust_time_vs_detected.pdf"),
           ggplot(r, aes(detected, time, colour = method)) +
             geom_smooth(method = "loess", se = TRUE, na.rm = TRUE) +
             scale_y_log10() + theme_minimal(base_size = 11) +
             labs(x = "Detected number of segments", y = "Time (s, log scale)",
                  colour = "Method"), width = 8, height = 5)
  }
  if (!is.null(a)) {
    a$method[a$method == "SVP AR1 estimated rho"] <-
      "SVP AR1 / right / 2.000 log(n) (estimated rho)"
    a$method[a$method == "SVP AR1"] <-
      "SVP AR1 / right / 2.000 log(n)"
    a$method[a$method == "SVP BIC"] <-
      "SVP Gaussian / right / 2.000 log(n)"
    ggsave(file.path(root, "2_ar1_time_vs_n.pdf"),
           plot_runtime_refresh(a, "n"), width = 8, height = 5)
    ggsave(file.path(root, "2_ar1_time_vs_detected.pdf"),
           ggplot(a, aes(detected, time, colour = method)) +
             geom_smooth(method = "loess", se = TRUE, na.rm = TRUE) +
             scale_y_log10() + theme_minimal(base_size = 11) +
             labs(x = "Detected number of segments", y = "Time (s, log scale)",
                  colour = "Method"), width = 8, height = 5)
  }
}

write_well_log_figure <- function(root = FIGURE_ROOT, wilcoxon_row = NULL,
                                  mood_row = NULL) {
  y <- DeCAFS::oilWell
  n <- length(y)
  scale <- stats::mad(y, constant = 1.4826)
  pelt <- changepoint::cpt.mean(
    y / scale, method = "PELT", penalty = "Manual", pen.value = 70
  )
  rfpop <- robseg::Rob_seg.std(
    y / scale, loss = "Outlier", lambda = 70, lthreshold = 2
  )
  if (is.null(wilcoxon_row)) wilcoxon_row <- data.frame(mode = "right", value = 1.5)
  if (is.null(mood_row)) mood_row <- data.frame(mode = "right", value = 0.01)
  w <- normalise_boundaries(calibration_candidate_fit(
    y, "wilcoxon", wilcoxon_row$mode, wilcoxon_row$value, true_segments = 30L
  ), n)
  m <- normalise_boundaries(calibration_candidate_fit(
    y, "mood", mood_row$mode, mood_row$value, true_segments = 10L
  ), n)
  plot_piecewise_constant <- function(y, cps, main, col_line) {
    cps <- sort(unique(c(as.integer(cps), length(y))))
    cps <- cps[cps >= 1 & cps <= length(y)]
    starts <- c(1L, head(cps, -1L) + 1L)
    plot(seq_along(y), y, type = "p", pch = ".", cex = 3,
         main = main, ylab = "", xlab = "")
    abline(v = cps, col = "black", lty = 2)
    for (i in seq_along(cps)) {
      idx <- starts[i]:cps[i]
      lines(idx, rep(median(y[idx]), length(idx)), lwd = 2, col = col_line)
    }
  }
  out <- file.path(root, "4_segmentation_logdata.pdf")
  pdf(out, width = 10, height = 8)
  op <- par(mfrow = c(4, 1), mar = c(4, 4, 2, 0))
  on.exit({par(op); dev.off()}, add = TRUE)
  plot_piecewise_constant(y, pelt@cpts, "PELT", 4)
  plot_piecewise_constant(y, rfpop$t.est, "Robust FPOP", 3)
  plot_piecewise_constant(
    y, w,
    paste0("SVP Wilcoxon / ", wilcoxon_row$mode, " / ",
           fmt_constant("wilcoxon", wilcoxon_row$value)), 2
  )
  plot_piecewise_constant(
    y, m,
    paste0("SVP MedianMood / ", mood_row$mode, " / ",
           fmt_constant("mood", mood_row$value)), 5
  )
  invisible(out)
}

canonical_algorithm <- function(x) {
  dplyr::case_when(
    x == "PELT" ~ "PELT",
    grepl("SVP Gaussian / right / 1\\.500", x) ~ "SVP BIC calibrated",
    grepl("SVP Gaussian / right / 2\\.000", x) ~ "SVP BIC",
    grepl("SVP Gaussian / both /", x) ~ "SVP BIC multiscale",
    grepl("MedianMood", x) ~ "SVP MedianMood",
    grepl("Wilcoxon / right /", x) ~ "SVP Wilcoxon",
    grepl("Wilcoxon / both /", x) ~ "SVP Wilcoxon multiscale",
    grepl("AR1Focus / right /", x) ~ "SVP AR1Focus",
    grepl("AR1Focus / both /", x) ~ "SVP AR1Focus multiscale",
    TRUE ~ x
  )
}

write_change_report <- function(new_results, root = FIGURE_ROOT) {
  old_paths <- c(
    gaussian = file.path("..", "svpChange2_old", "simulations_paper",
                         "power_gaussian", "results.rds"),
    robust = file.path("..", "svpChange2_old", "simulations_paper",
                       "power_robust", "results.rds"),
    ar1 = file.path("..", "svpChange2_old", "simulations_paper",
                    "power_ar1", "results.rds")
  )
  current_rows <- list()
  comparison_rows <- list()
  for (name in names(new_results)) {
    current <- new_results[[name]]
    if (is.null(current) || !nrow(current)) next
    current$comparison_algorithm <- canonical_algorithm(current$algorithm)
    current_rows[[name]] <- current |>
      group_by(algorithm) |>
      summarise(
        study = name, mean_F1 = mean(F1, na.rm = TRUE),
        mean_segments = mean(NumSegments, na.rm = TRUE), .groups = "drop"
      )
    old_path <- old_paths[[name]]
    if (!is.null(old_path) && file.exists(old_path)) {
      old <- readRDS(old_path)
      old$comparison_algorithm <- canonical_algorithm(old$algorithm)
      new_summary <- current |>
        group_by(comparison_algorithm) |>
        summarise(
          new_mean_F1 = mean(F1, na.rm = TRUE),
          new_null_F1 = mean(F1[pattern == "none"], na.rm = TRUE),
          new_mean_segments = mean(NumSegments, na.rm = TRUE), .groups = "drop"
        )
      old_summary <- old |>
        group_by(comparison_algorithm) |>
        summarise(
          old_mean_F1 = mean(F1, na.rm = TRUE),
          old_null_F1 = mean(F1[pattern == "none"], na.rm = TRUE),
          old_mean_segments = mean(NumSegments, na.rm = TRUE), .groups = "drop"
        )
      comparison_rows[[name]] <- full_join(
        new_summary, old_summary, by = "comparison_algorithm"
      ) |>
        mutate(
          study = name,
          delta_F1 = new_mean_F1 - old_mean_F1,
          delta_null_F1 = new_null_F1 - old_null_F1,
          segment_change = new_mean_segments / old_mean_segments - 1,
          flag_F1 = abs(delta_F1) >= 0.03,
          flag_null_F1 = abs(delta_null_F1) >= 0.02,
          flag_segments = abs(segment_change) >= 0.10,
          flag_any = flag_F1 | flag_null_F1 | flag_segments
        ) |>
        select(study, comparison_algorithm, everything())
    }
  }
  report <- bind_rows(current_rows)
  comparison <- bind_rows(comparison_rows)
  write.csv(report, file.path(root, "performance_change_report.csv"), row.names = FALSE)
  write.csv(comparison, file.path(root, "performance_change_comparison.csv"),
            row.names = FALSE)

  ## Runtime CSVs contain enough raw information for a paired old/current
  ## comparison.  We report both ratios and changes in log-log scaling slopes.
  runtime_rows <- list()
  for (family in c("gaussian", "robust", "ar1")) {
    new_path <- file.path("simulations", paste0("time_", family),
                          paste0(family, "_time_results.csv"))
    old_path <- file.path("..", "svpChange2_old", "simulations_paper",
                          paste0("time_", family),
                          paste0(family, "_time_results.csv"))
    if (!file.exists(new_path) || !file.exists(old_path)) next
    current <- read.csv(new_path, stringsAsFactors = FALSE)
    old <- read.csv(old_path, stringsAsFactors = FALSE)
    key <- c("experiment", "n", "k", "rep", "method")
    method_summary <- function(x) {
      x |>
        group_by(experiment, method) |>
        summarise(mean_time = mean(time, na.rm = TRUE), .groups = "drop")
    }
    means <- full_join(
      method_summary(current) |> rename(new_mean_time = mean_time),
      method_summary(old) |> rename(old_mean_time = mean_time),
      by = c("experiment", "method")
    )
    slope <- function(x, method, experiment) {
      x <- x[x$method == method & x$experiment == experiment & x$time > 0, ]
      xvar <- if (experiment == "vary_n") x$n else x$k
      keep <- is.finite(xvar) & xvar > 0
      x <- x[keep, , drop = FALSE]
      xvar <- xvar[keep]
      if (nrow(x) < 3L || length(unique(xvar)) < 3L) return(NA_real_)
      fit <- try(lm(log(time) ~ log(xvar), data = x), silent = TRUE)
      if (inherits(fit, "try-error")) return(NA_real_)
      estimate <- unname(coef(fit)[2L])
      if (!is.finite(estimate)) NA_real_ else estimate
    }
    means$slope_new <- vapply(seq_len(nrow(means)), function(i) {
      slope(current, means$method[i], means$experiment[i])
    }, numeric(1))
    means$slope_old <- vapply(seq_len(nrow(means)), function(i) {
      slope(old, means$method[i], means$experiment[i])
    }, numeric(1))
    means$runtime_ratio <- means$new_mean_time / means$old_mean_time
    means$delta_slope <- means$slope_new - means$slope_old
    means$flag_runtime <- pmax(means$runtime_ratio, 1 / means$runtime_ratio,
                                na.rm = TRUE) >= 2
    means$flag_slope <- abs(means$delta_slope) >= 0.2
    means$flag_any <- means$flag_runtime | means$flag_slope
    means$study <- paste0("runtime_", family)
    runtime_rows[[family]] <- select(means, study, everything())
  }
  runtime <- bind_rows(runtime_rows)
  write.csv(runtime, file.path(root, "runtime_change_comparison.csv"),
            row.names = FALSE)
  flags <- if (nrow(comparison)) sum(comparison$flag_any, na.rm = TRUE) else 0L
  flags <- flags + if (nrow(runtime)) sum(runtime$flag_any, na.rm = TRUE) else 0L
  writeLines(c(
    "Performance-change report for SVP_NEW_Figures2.",
    "Current figures use svpChange2 v0.2.0 and explicit right/both modes.",
    "Compared with the preserved v0.1.0 raw result objects when available.",
    "Flags: |delta F1| >= 0.03, null-F1 delta >= 0.02, segment change >= 10%,",
    "runtime ratio >= 2 (in either direction), or log-log slope change >= 0.2.",
    paste0("Number of flagged comparison rows: ", flags, "."),
    "See performance_change_comparison.csv and runtime_change_comparison.csv."
  ), file.path(root, "performance_change_report.md"))
  invisible(list(current = report, comparison = comparison, runtime = runtime))
}

write_refresh_readme <- function(root = FIGURE_ROOT, run_full = FALSE) {
  mode <- if (run_full) {
    "full recalibration and power-grid regeneration"
  } else {
    "quick refresh from the maintained v0.2.0 result objects"
  }
  writeLines(c(
    "SVP_NEW_Figures2",
    "=================",
    "",
    paste0("This folder was generated in ", mode, "."),
    "The original SVP_NEW_Figures folder is never modified.",
    "All SVP figure legends state the explicit subtests mode and calibration.",
    "",
    "To regenerate with the approved calibration grids (1000 null and 200",
    "signal replicates), run from the svpChange2 package root:",
    "",
    "  SVP_REBUILD_FIGURES2=true SVP_REBUILD_FULL=true \\",
    "    SVP_REBUILD_WORKERS=10 Rscript simulations/other_simus/maintenance/",
    "    rebuild_svp_new_figures2.R",
    "",
    "The quick mode is reproducible and uses the current maintained result",
    "objects.  See performance_change_comparison.csv and",
    "runtime_change_comparison.csv for paired v0.1.0/v0.2.0 checks."
  ), file.path(root, "README.md"))
}

rebuild_svp_new_figures2 <- function(
    run_full = FALSE, workers = 1L,
    null_reps = 1000L, signal_reps = 200L) {
  dir.create(FIGURE_ROOT, recursive = TRUE, showWarnings = FALSE)
  ## Remove only obsolete filenames produced by an earlier draft of this
  ## driver; the historical SVP_NEW_Figures directory is untouched.
  obsolete <- c(
    paste0("3_", c("signal_scenarios", "f1_vs_jump", "precision_vs_jump",
                     "recall_vs_jump", "correct_number_vs_jump",
                     "numsegments_vs_jump", "mse_vs_jump",
                     "changepoint_distributions"), ".pdf"),
    paste0("5_", c("signal_scenarios", "f1_vs_jump", "precision_vs_jump",
                     "recall_vs_jump", "correct_number_vs_jump",
                     "mse_vs_jump", "changepoint_distributions"), ".pdf")
  )
  obsolete <- c(
    obsolete,
    paste0("3_", c("f1_vs_jump", "precision_vs_jump", "recall_vs_jump",
                    "correct_number_vs_jump", "numsegments_vs_jump",
                    "mse_vs_jump"), "_robust.pdf"),
    paste0("5_", c("f1_vs_jump", "precision_vs_jump", "recall_vs_jump",
                    "correct_number_vs_jump", "mse_vs_jump"), "_ar1.pdf")
  )
  unlink(file.path(FIGURE_ROOT, obsolete), force = TRUE)
  if (run_full && !exists("simulate_ar1_noise", mode = "function")) {
    source_ar1_power_helpers()
  }
  calibration <- list()
  if (run_full) {
    candidates <- default_calibration_candidates()
    calibration$gaussian <- calibrate_svp_family(
      "gaussian", candidates$gaussian, 1000L, stats::rnorm,
      seed = 710000L, null_reps = null_reps, signal_reps = signal_reps,
      workers = workers
    )
    calibration$ar1 <- calibrate_svp_family(
      "ar1", candidates$ar1, 600L,
      function(n) simulate_ar1_noise(n, 0.8), seed = 720000L,
      null_reps = null_reps, signal_reps = signal_reps, rho = 0.8,
      workers = workers
    )
    calibration$wilcoxon <- calibrate_svp_family(
      "wilcoxon", candidates$wilcoxon, 1000L,
      function(n) stats::rt(n, df = 2), seed = 730000L,
      null_reps = null_reps, signal_reps = signal_reps, workers = workers
    )
    calibration$mood <- calibrate_svp_family(
      "mood", candidates$mood, 1000L,
      function(n) stats::rt(n, df = 2), seed = 740000L,
      null_reps = null_reps, signal_reps = signal_reps, workers = workers
    )
    saveRDS(calibration, file.path(FIGURE_ROOT, "calibration.rds"))
    for (nm in names(calibration)) {
      write.csv(calibration[[nm]]$scores,
                file.path(FIGURE_ROOT, paste0("calibration_", nm, ".csv")),
                row.names = FALSE)
    }
  } else if (file.exists(file.path(FIGURE_ROOT, "calibration.rds"))) {
    calibration <- readRDS(file.path(FIGURE_ROOT, "calibration.rds"))
  }

  ## If full calibration is not requested, use the current v0.2.0 generated
  ## results as an immediate reproducible figure refresh. Full mode replaces
  ## these with newly calibrated grids.
  if (run_full) {
    gaussian <- run_gaussian_refresh(
      make_candidate_table(calibration$gaussian, "gaussian"), workers
    )
    robust <- run_robust_refresh(
      bind_rows(
        make_candidate_table(calibration$wilcoxon, "wilcoxon"),
        make_candidate_table(calibration$mood, "mood")
      ), workers
    )
    ar1 <- run_ar1_refresh(make_candidate_table(calibration$ar1, "ar1"), workers)
  } else {
    gaussian <- readRDS("simulations/power_gaussian/reproduced_v020_matched/results.rds")
    robust <- readRDS("simulations/power_robust/results.rds")
    ar1 <- readRDS("simulations/power_ar1/results.rds")
    renamed <- rename_quick_results(gaussian, robust, ar1)
    gaussian <- renamed$gaussian
    robust <- renamed$robust
    ar1 <- renamed$ar1
  }
  saveRDS(gaussian, file.path(FIGURE_ROOT, "gaussian_results.rds"))
  saveRDS(robust, file.path(FIGURE_ROOT, "robust_results.rds"))
  saveRDS(ar1, file.path(FIGURE_ROOT, "ar1_results.rds"))
  write.csv(select(gaussian, -changepoints),
            file.path(FIGURE_ROOT, "gaussian_results.csv"), row.names = FALSE)
  write.csv(select(robust, -changepoints),
            file.path(FIGURE_ROOT, "robust_results.csv"), row.names = FALSE)
  write.csv(select(ar1, -changepoints),
            file.path(FIGURE_ROOT, "ar1_results.csv"), row.names = FALSE)

  ## Existing saved robust/AR1 results retain the paper’s method set in quick
  ## mode; full mode uses the selected current configurations above.
  write_power_figures(gaussian, "1", 1000L, 0.6, stats::rnorm)
  write_power_figures(robust, "3", 1000L, 0.6,
                      function(n) stats::rt(n, df = 2),
                      file_suffix = "_robust")
  ## The AR(1) simulator is defined by the maintained power-study helper.
  ## Load it before drawing the scenario panel (quick mode does not call the
  ## AR(1) fitting routine above).
  if (!exists("simulate_ar1_noise", mode = "function")) source_ar1_power_helpers()
  write_power_figures(ar1, "5", 600L, 1.4,
                      function(n) simulate_ar1_noise(n, 0.8),
                      include_numsegments = FALSE, file_suffix = "_ar1")
  write_runtime_figures()
  well_wilcoxon <- if (run_full) calibration$wilcoxon$selected else NULL
  well_mood <- if (run_full) calibration$mood$selected else NULL
  write_well_log_figure(wilcoxon_row = well_wilcoxon, mood_row = well_mood)
  write_change_report(list(gaussian = gaussian, robust = robust, ar1 = ar1))
  write_refresh_readme(run_full = run_full)

  manifest <- c(
    paste0("run_full = ", run_full),
    paste0("workers = ", workers),
    paste0("git_sha = ", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
    paste0("figure_root = ", FIGURE_ROOT),
    capture.output(sessionInfo())
  )
  writeLines(manifest, file.path(FIGURE_ROOT, "manifest.txt"))
  invisible(list(gaussian = gaussian, robust = robust, ar1 = ar1,
                calibration = calibration))
}

if (identical(tolower(Sys.getenv("SVP_REBUILD_FIGURES2")), "true")) {
  rebuild_svp_new_figures2(
    run_full = identical(tolower(Sys.getenv("SVP_REBUILD_FULL")), "true"),
    workers = as.integer(Sys.getenv("SVP_REBUILD_WORKERS", "1"))
  )
}
