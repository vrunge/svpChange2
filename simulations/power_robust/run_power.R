## Heavy-tail power study: paper Figures 3 and 7 plus companion figures.

library(svpChange2)
library(changepoint)
library(robseg)
library(dplyr)
library(ggplot2)
library(future)
library(future.apply)

source(file.path("simulations", "power_common.R"))

ROBUST_ROOT <- file.path("simulations", "power_robust")

mood_threshold <- function(n, segments, alpha = 0.01) {
  splits <- n / segments - 1
  single_alpha <- 1 - (1 - alpha)^(1 / splits)
  stats::qchisq(1 - single_alpha, df = segments)
}

fit_rfpop <- function(y, penalty_constant) {
  scale <- stats::mad(y, constant = 1.4826)
  fit <- robseg::Rob_seg.std(
    y / scale, loss = "Outlier",
    lambda = penalty_constant * log(length(y)), lthreshold = 3
  )
  list(
    boundaries = normalise_boundaries(head(fit$t.est, -1L), length(y)),
    fitted = as.numeric(fit$smt) * scale
  )
}

calibrate_rfpop_penalty <- function(
    n = 1000L, grid = seq(2, 6, 0.25), reps = 200L,
    target_no_change = 0.98, seed = 90210L) {
  scores <- lapply(grid, function(constant) {
    correct <- logical(reps)
    for (rep in seq_len(reps)) {
      set.seed(seed + rep)
      fit <- fit_rfpop(stats::rt(n, df = 2), constant)
      correct[rep] <- length(fit$boundaries) == 1L
    }
    data.frame(constant = constant, no_change_probability = mean(correct))
  }) |>
    dplyr::bind_rows()
  admissible <- dplyr::filter(scores, no_change_probability >= target_no_change)
  if (!nrow(admissible)) stop("RFPOP calibration grid does not attain null target")
  list(constant = min(admissible$constant), table = scores)
}

fit_robust_methods <- function(y, true_segments, rfpop_constant) {
  n <- length(y)
  pelt <- changepoint::cpt.mean(
    y, method = "PELT", penalty = "Manual", pen.value = 2 * log(n)
  )
  wilcoxon_gamma <- sqrt((n / true_segments)^3 / 12)
  paper_rfpop <- fit_rfpop(y, 2)
  calibrated_rfpop <- fit_rfpop(y, rfpop_constant)
  calibrated_label <- sprintf("RFPOP (null-matched c=%g)", rfpop_constant)
  boundaries <- list(
    "PELT" = normalise_boundaries(changepoint::cpts(pelt), n),
    "RFPOP (paper)" = paper_rfpop$boundaries,
    "SVP (MedianMood)" = normalise_boundaries(
      SVP(y, mood_threshold(n, true_segments), "MedianMoodCost")$changepoints, n
    ),
    "SVP (Wilcoxon)" = normalise_boundaries(
      SVP(y, 1.5 * wilcoxon_gamma, "WilcoxonCost")$changepoints, n
    ),
    "SVP Wilcoxon TRUE/TRUE c=1.8" = normalise_boundaries(
      SVP(y, 1.8 * wilcoxon_gamma, "WilcoxonCost",
          prune_after_if_unvalid = TRUE,
          prune_before_if_invalid = TRUE)$changepoints,
      n
    )
  )
  boundaries[[calibrated_label]] <- calibrated_rfpop$boundaries
  fitted <- lapply(boundaries, function(x) fitted_piecewise_mean(y, x))
  fitted[["RFPOP (paper)"]] <- paper_rfpop$fitted
  fitted[[calibrated_label]] <- calibrated_rfpop$fitted
  list(boundaries = boundaries, fitted = fitted)
}

run_robust_power <- function(
    n = 1000L, jump_sizes = seq(0.1, 4, 0.1), reps = 100L,
    workers = 1L, seed = 123L,
    rfpop_constant = NULL) {
  if (is.null(rfpop_constant)) {
    calibration <- calibrate_rfpop_penalty(n = n)
    rfpop_constant <- calibration$constant
    utils::write.csv(calibration$table,
                     file.path(ROBUST_ROOT, "rfpop_null_calibration.csv"),
                     row.names = FALSE)
  }
  design <- expand.grid(pattern = POWER_PATTERNS, jump = jump_sizes,
                        rep = seq_len(reps), stringsAsFactors = FALSE)
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  rows <- future.apply::future_lapply(seq_len(nrow(design)), function(i) {
    d <- design[i, ]
    set.seed(seed + i)
    mu <- generate_power_signal(n, d$pattern, d$jump)
    y <- mu + stats::rt(n, df = 2)
    truth <- normalise_boundaries(which(diff(mu) != 0), n)
    fits <- fit_robust_methods(y, length(truth), rfpop_constant)
    output <- dplyr::bind_rows(lapply(names(fits$boundaries), function(method) {
      result_row(d$pattern, d$jump, d$rep, method,
                 fits$boundaries[[method]], truth, mu, fits$fitted[[method]],
                 tolerance = round(n * 0.0025))
    }))
    output$changepoints <- unname(fits$boundaries)
    output
  }, future.seed = TRUE)
  dplyr::bind_rows(rows)
}

robust_scenario_plot <- function(n = 1000L, jump = 0.6, seed = 999L) {
  set.seed(seed)
  data <- dplyr::bind_rows(lapply(POWER_PATTERNS, function(pattern) {
    mu <- generate_power_signal(n, pattern, jump)
    data.frame(time = seq_len(n), value = mu + stats::rt(n, df = 2),
               mean = mu, pattern = pattern)
  }))
  data$pattern <- factor(data$pattern, POWER_PATTERNS)
  ggplot2::ggplot(data, ggplot2::aes(time, value)) +
    ggplot2::geom_point(alpha = 0.2, size = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = mean), colour = "red", linewidth = 0.8) +
    ggplot2::facet_wrap(~pattern, nrow = 2L) +
    ggplot2::labs(x = "Time (Index)", y = "Value") +
    ggplot2::coord_cartesian(ylim = c(-10, 10)) +
    paper_theme() + ggplot2::theme(legend.position = "none")
}

run_and_save_robust <- function(workers = 1L) {
  results <- run_robust_power(workers = workers)
  save_power_outputs(results, ROBUST_ROOT, robust_scenario_plot())
  invisible(results)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_and_save_robust()
}
