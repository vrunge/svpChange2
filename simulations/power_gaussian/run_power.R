## Gaussian power study: SVP paper Figures 1, 2, 6, and 8, plus two metrics.

library(svpChange2)
library(changepoint)
library(dplyr)
library(ggplot2)

source(file.path("simulations", "power_common.R"))

GAUSSIAN_ROOT <- file.path("simulations", "power_gaussian")
GAUSSIAN_TRUE_TRUE_CONSTANT <- 1.8

fit_gaussian_methods <- function(y) {
  n <- length(y)
  pelt <- changepoint::cpt.mean(
    y, method = "PELT", penalty = "Manual", pen.value = 2 * log(n)
  )
  list(
    "PELT" = normalise_boundaries(changepoint::cpts(pelt), n),
    "SVP BIC calibrated" = normalise_boundaries(
      SVP(y, 1.5 * log(n), "gaussian_mean")$changepoints, n
    ),
    "SVP BIC" = normalise_boundaries(
      SVP(y, 2 * log(n), "gaussian_mean")$changepoints, n
    ),
    "SVP BIC multiscale" = normalise_boundaries(
      SVP(y, GAUSSIAN_TRUE_TRUE_CONSTANT * log(n), "gaussian_mean",
          prune_after_if_unvalid = TRUE,
          prune_before_if_invalid = TRUE)$changepoints,
      n
    )
  )
}

run_gaussian_power <- function(
    n = 1000L, jump_sizes = seq(0.1, 2, 0.1), reps = 100L,
    workers = power_default_workers(), seed = 123L) {
  run_power_grid(
    n, jump_sizes, reps,
    simulate_noise = stats::rnorm,
    fit_methods = function(y, true_segments) fit_gaussian_methods(y),
    tolerance = round(n * 0.0025),
    workers = workers,
    seed = seed
  )
}

run_and_save_gaussian <- function(workers = power_default_workers()) {
  results <- run_gaussian_power(workers = workers)
  save_power_outputs(
    results, GAUSSIAN_ROOT,
    power_scenario_plot(1000L, 0.6, stats::rnorm),
    plot_results = set_algorithm_order(
      results,
      c("PELT", "SVP BIC", "SVP BIC calibrated",
        "SVP BIC multiscale")
    )
  )
  invisible(results)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_and_save_gaussian()
}
