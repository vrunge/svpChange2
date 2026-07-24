## Runtime scaling under independent Gaussian noise.

library(svpChange2)
library(changepoint)
library(dplyr)
library(ggplot2)

source(file.path("simulations", "time_common.R"))

GAUSSIAN_TIME_ROOT <- file.path("simulations", "time_gaussian")
GAUSSIAN_TIME_METHODS <- c(
  "PELT", "SVP BIC", "SVP BIC calibrated", "SVP BIC multiscale"
)

simulate_gaussian_time_data <- function(n, changes = 0L, jump = 10) {
  stats::rnorm(n, alternating_mean(n, changes, jump), 1)
}

fit_gaussian_time_method <- function(y, method) {
  n <- length(y)
  fit <- switch(
    method,
    "PELT" = changepoint::cpt.mean(
      y, method = "PELT", penalty = "Manual", pen.value = 2 * log(n)
    ),
    "SVP BIC" = SVP(
      y, 2 * log(n), "gaussian_mean",
      prune_after_if_unvalid = TRUE,
      prune_before_if_invalid = FALSE
    ),
    "SVP BIC calibrated" = SVP(
      y, 1.5 * log(n), "gaussian_mean",
      prune_after_if_unvalid = TRUE,
      prune_before_if_invalid = FALSE
    ),
    "SVP BIC multiscale" = SVP(
      y, 1.5 * log(n), "gaussian_mean",
      prune_after_if_unvalid = TRUE,
      prune_before_if_invalid = TRUE
    )
  )
  if (method == "PELT") {
    return(sort(unique(c(changepoint::cpts(fit), n))))
  }
  extract_svp_boundaries(fit, n)
}

run_time_gaussian <- function(
    n_values = round(2^(seq(8, 14, length.out = 20))),
    k_values = 0:100,
    fixed_n = 10000L,
    reps_n = 20L,
    reps_k = 10L,
    seed = 456L) {
  benchmark_time_study(
    n_values, k_values, fixed_n, reps_n, reps_k,
    methods = GAUSSIAN_TIME_METHODS,
    simulate_data = simulate_gaussian_time_data,
    fit_method = fit_gaussian_time_method,
    seed = seed
  )
}

run_and_save_time_gaussian <- function() {
  results <- run_time_gaussian()
  save_time_outputs(results, GAUSSIAN_TIME_ROOT, "gaussian_time")
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_and_save_time_gaussian()
}
