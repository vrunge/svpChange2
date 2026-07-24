## Runtime scaling under Student-t(2) noise.

library(svpChange2)
library(changepoint)
library(dplyr)
library(ggplot2)

source(file.path("simulations", "time_common.R"))

ROBUST_TIME_ROOT <- file.path("simulations", "time_robust")
ROBUST_TIME_METHODS <- c(
  "PELT", "SVP Wilcoxon", "SVP MedianMood",
  "SVP BIC multiscale"
)
ROBUST_TIME_MAX_N <- 1000L

simulate_robust_time_data <- function(n, changes = 0L, jump = 4) {
  alternating_mean(n, changes, jump) + stats::rt(n, df = 2)
}

fit_robust_time_method <- function(y, method) {
  n <- length(y)
  fit <- switch(
    method,
    "PELT" = changepoint::cpt.mean(
      y, method = "PELT", penalty = "Manual", pen.value = 2 * log(n)
    ),
    "SVP Wilcoxon" = SVP(y, 2 * log(n), "WilcoxonCost"),
    "SVP MedianMood" = SVP(y, 2 * log(n), "MedianMoodCost"),
    "SVP BIC multiscale" =
      SVP(y, 1.5 * log(n), "gaussian_mean", TRUE, TRUE)
  )
  if (method == "PELT") {
    return(sort(unique(c(changepoint::cpts(fit), n))))
  }
  extract_svp_boundaries(fit, n)
}

run_time_robust <- function(
    n_values = round(2^(seq(8, log2(ROBUST_TIME_MAX_N), length.out = 12))),
    k_values = 0:40,
    fixed_n = ROBUST_TIME_MAX_N,
    reps_n = 10L,
    reps_k = 5L,
    seed = 910L) {
  if (any(n_values > ROBUST_TIME_MAX_N) || fixed_n > ROBUST_TIME_MAX_N) {
    stop("the robust time study is limited to n <= ", ROBUST_TIME_MAX_N)
  }
  benchmark_time_study(
    n_values, k_values, fixed_n, reps_n, reps_k,
    methods = ROBUST_TIME_METHODS,
    simulate_data = simulate_robust_time_data,
    fit_method = fit_robust_time_method,
    seed = seed
  )
}

run_and_save_time_robust <- function() {
  results <- run_time_robust()
  save_time_outputs(results, ROBUST_TIME_ROOT, "robust_time")
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_and_save_time_robust()
}
