## AR(1) power study parallel to the Gaussian paper simulation.

library(svpChange2)
library(changepoint)
library(dplyr)
library(DeCAFS)
library(ggplot2)

source(file.path("simulations", "power_common.R"))

AR1_ROOT <- file.path("simulations", "power_ar1")
AR1_TRUE_TRUE_CONSTANT <- 4.5
AR1_DECAFS_CONSTANT <- 2

simulate_ar1_noise <- function(n, rho = 0.8, marginal_sd = 1) {
  dataRWAR(n = n, phi = rho, type = "none")$y
}

ar1_cost_matrix <- function(y, rho, innovation_variance) {
  n <- length(y)
  transition <- y[-1L] - rho * y[-n]
  sum1 <- c(0, cumsum(transition))
  sum2 <- c(0, cumsum(transition^2))
  first_weight <- 1 - rho^2
  mean_weight <- 1 - rho
  cost <- matrix(0, n, n)
  for (s in seq_len(n)) {
    t <- s:n
    count <- t - s
    transition_sum <- sum1[t] - sum1[s]
    transition_sum2 <- sum2[t] - sum2[s]
    denominator <- first_weight + count * mean_weight^2
    numerator <- first_weight * y[s] + mean_weight * transition_sum
    constant <- first_weight * y[s]^2 + transition_sum2
    cost[s, t] <- pmax(0, constant - numerator^2 / denominator) /
      innovation_variance
  }
  cost
}

ar1_approximate_partition <- function(y, rho, innovation_variance, penalty) {
  n <- length(y)
  cost <- ar1_cost_matrix(y, rho, innovation_variance)
  best <- rep(Inf, n + 1L)
  previous <- integer(n)
  best[1L] <- -penalty
  for (t in seq_len(n)) {
    candidates <- 0:(t - 1L)
    values <- best[candidates + 1L] +
      cost[cbind(candidates + 1L, rep(t, t))] + penalty
    winner <- which.min(values)
    best[t + 1L] <- values[winner]
    previous[t] <- candidates[winner]
  }
  boundaries <- integer(0)
  endpoint <- n
  while (endpoint > 0L) {
    boundaries <- c(endpoint, boundaries)
    endpoint <- previous[endpoint]
  }
  boundaries
}

decafs_boundaries <- function(y, rho, innovation_variance,
                              penalty = AR1_DECAFS_CONSTANT * log(length(y))) {
  n <- length(y)
  fit <- DeCAFS::DeCAFS(
    y,
    beta = penalty,
    modelParam = list(
      sdEta = 0,
      #sdNu = sqrt(innovation_variance),
      sdNu = 1,
      phi = rho
    )
  )
  normalise_boundaries(fit$changepoints, n)
}

fit_ar1_methods <- function(y, rho = 0.8, svp_constant = 3.75,
                            true_true_constant = AR1_TRUE_TRUE_CONSTANT) {
  if (!requireNamespace("DeCAFS", quietly = TRUE)) {
    stop("the optional package 'DeCAFS' is required for the AR(1) power study")
  }
  n <- length(y)
  innovation_variance <- 1 - rho^2
  inflated <- changepoint::cpt.mean(
    y, method = "PELT", penalty = "Manual",
    pen.value = 3 * log(n) * (1 + rho) / (1 - rho)
  )
  list(
    "PELT inflated" =
      normalise_boundaries(changepoint::cpts(inflated), n),
    "PELT AR1 approximate" =
      ar1_approximate_partition(y, rho, 1, 2.2 * log(n)),
    "DeCAFS AR1" =
      decafs_boundaries(y, rho, innovation_variance),
    "SVP AR1Focus" = normalise_boundaries(
      SVP(y, svp_constant * log(n), "AR1Focus",
          prune_after_if_unvalid = TRUE,
          prune_before_if_invalid = FALSE,
          rho = rho, sigma2 = innovation_variance)$changepoints,
      n
    ),
    "SVP AR1Focus multiscale" = normalise_boundaries(
      SVP(y, true_true_constant * log(n), "AR1Focus",
          prune_after_if_unvalid = TRUE,
          prune_before_if_invalid = TRUE,
          rho = rho, sigma2 = innovation_variance)$changepoints,
      n
    )
  )
}

run_ar1_power <- function(
    n = 600L, rho = 0.8, jump_sizes = seq(1, 3, 0.1), reps = 100L,
    workers = power_default_workers(), seed = 123L) {
  run_power_grid(
    n, jump_sizes, reps,
    simulate_noise = function(size) simulate_ar1_noise(size, rho),
    fit_methods = function(y, true_segments) fit_ar1_methods(y, rho),
    tolerance = max(5L, round(n * 0.0025)),
    workers = workers,
    seed = seed
  )
}

run_and_save_ar1 <- function(workers = power_default_workers()) {
  results <- run_ar1_power(workers = workers)
  save_power_outputs(
    results, AR1_ROOT,
    power_scenario_plot(
      600L, 0.8, function(size) simulate_ar1_noise(size, 0.8)
    ),
    selected_jump = 1.4,
    plot_results = set_algorithm_order(
      results,
      c("PELT AR1 approximate", "PELT inflated", "DeCAFS AR1",
        "SVP AR1Focus", "SVP AR1Focus multiscale")
    )
  )
  invisible(results)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_and_save_ar1()
}
