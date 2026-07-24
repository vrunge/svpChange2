## Heavy-tail power study: paper Figures 3 and 7 plus companion figures.

library(svpChange2)
library(changepoint)
library(robseg)
library(dplyr)
library(ggplot2)

source(file.path("simulations", "power_common.R"))

ROBUST_ROOT <- file.path("simulations", "power_robust")
ROBUST_TRUE_TRUE_CONSTANT <- 1.75

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
    target_no_change = 0.98, seed = 90210L,
    workers = power_default_workers()) {
  scores <- power_lapply(grid, function(constant) {
    correct <- logical(reps)
    for (rep in seq_len(reps)) {
      set_power_seed(seed + rep)
      fit <- fit_rfpop(stats::rt(n, df = 2), constant)
      correct[rep] <- length(fit$boundaries) == 1L
    }
    data.frame(constant = constant, no_change_probability = mean(correct))
  }, workers = workers) |>
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
  boundaries <- list(
    "PELT" = normalise_boundaries(changepoint::cpts(pelt), n),
    "RFPOP paper" = paper_rfpop$boundaries,
    "SVP MedianMood" = normalise_boundaries(
      SVP(y, mood_threshold(n, true_segments), "MedianMoodCost")$changepoints, n
    ),
    "SVP Wilcoxon" = normalise_boundaries(
      SVP(y, 1.5 * wilcoxon_gamma, "WilcoxonCost")$changepoints, n
    ),
    "SVP Wilcoxon multiscale" = normalise_boundaries(
      SVP(y, ROBUST_TRUE_TRUE_CONSTANT * wilcoxon_gamma, "WilcoxonCost",
          prune_after_if_unvalid = TRUE,
          prune_before_if_invalid = TRUE)$changepoints,
      n
    )
  )
  boundaries[["RFPOP"]] <- calibrated_rfpop$boundaries
  fitted <- lapply(boundaries, function(x) fitted_piecewise_mean(y, x))
  fitted[["RFPOP paper"]] <- paper_rfpop$fitted
  fitted[["RFPOP"]] <- calibrated_rfpop$fitted
  list(boundaries = boundaries, fitted = fitted)
}

run_robust_power <- function(
    n = 1000L, jump_sizes = seq(0.1, 4, 0.1), reps = 100L,
    workers = power_default_workers(), seed = 123L,
    rfpop_constant = NULL) {
  if (is.null(rfpop_constant)) {
    calibration <- calibrate_rfpop_penalty(n = n, workers = workers)
    rfpop_constant <- calibration$constant
    utils::write.csv(calibration$table,
                     file.path(ROBUST_ROOT, "rfpop_null_calibration.csv"),
                     row.names = FALSE)
  }
  run_power_grid(
    n, jump_sizes, reps,
    simulate_noise = function(size) stats::rt(size, df = 2),
    fit_methods = function(y, true_segments) {
      fit_robust_methods(y, true_segments, rfpop_constant)
    },
    tolerance = round(n * 0.0025),
    workers = workers,
    seed = seed
  )
}

run_and_save_robust <- function(workers = power_default_workers()) {
  results <- run_robust_power(workers = workers)
  save_power_outputs(
    results, ROBUST_ROOT,
    power_scenario_plot(
      1000L, 0.6, function(size) stats::rt(size, df = 2),
      y_limits = c(-10, 10)
    ),
    plot_results = set_algorithm_order(
      results,
      c("PELT", "RFPOP", "SVP MedianMood", "SVP Wilcoxon",
        "SVP Wilcoxon multiscale"),
      drop = "RFPOP paper"
    )
  )
  invisible(results)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_and_save_robust()
}
