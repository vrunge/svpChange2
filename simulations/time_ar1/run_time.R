## Runtime scaling under stationary AR(1) noise.

library(svpChange2)
library(dplyr)
library(ggplot2)

source(file.path("simulations", "time_common.R"))

AR1_TIME_ROOT <- file.path("simulations", "time_ar1")
AR1_TIME_METHODS <- c(
  "SVP BIC", "SVP AR1", "SVP AR1 estimated rho"
)
AR1_TIME_RHO <- 0.6

simulate_ar1_time_data <- function(n, changes = 0L, rho = AR1_TIME_RHO,
                                   jump = 4) {
  mean <- alternating_mean(n, changes, jump)
  y <- numeric(n)
  y[1L] <- mean[1L] + stats::rnorm(1L) / sqrt(1 - rho^2)
  for (i in 2:n) {
    y[i] <- mean[i] + rho * (y[i - 1L] - mean[i - 1L]) + stats::rnorm(1L)
  }
  y
}

fit_ar1_time_method <- function(y, method, rho = AR1_TIME_RHO) {
  n <- length(y)
  fit <- switch(
    method,
    "SVP BIC" = SVP(y, 2 * log(n), "gaussian_mean"),
    "SVP AR1" =
      SVP(y, 2 * log(n), "AR1", rho = rho, sigma2 = 1),
    "SVP AR1 estimated rho" = SVP(
      y, 2 * log(n), "AR1", rho = AR1_rho(y), sigma2 = 1
    )
  )
  extract_svp_boundaries(fit, n)
}

run_time_ar1 <- function(
    n_values = round(2^(seq(8, 13, length.out = 12))),
    k_values = 0:40,
    fixed_n = 2000L,
    reps_n = 10L,
    reps_k = 5L,
    seed = 810L) {
  benchmark_time_study(
    n_values, k_values, fixed_n, reps_n, reps_k,
    methods = AR1_TIME_METHODS,
    simulate_data = simulate_ar1_time_data,
    fit_method = fit_ar1_time_method,
    seed = seed
  )
}

run_and_save_time_ar1 <- function() {
  results <- run_time_ar1()
  save_time_outputs(results, AR1_TIME_ROOT, "ar1_time")
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_and_save_time_ar1()
}
