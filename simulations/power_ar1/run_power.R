## AR(1) power study parallel to the Gaussian paper simulation.

library(svpChange2)
library(changepoint)
library(dplyr)
library(ggplot2)

source(file.path("simulations", "power_common.R"))

AR1_ROOT <- file.path("simulations", "power_ar1")
AR1_TRUE_TRUE_CONSTANT <- 1.75

simulate_ar1_noise <- function(n, rho = 0.8, marginal_sd = 1) {
  innovation_sd <- marginal_sd * sqrt(1 - rho^2)
  noise <- numeric(n)
  noise[1L] <- stats::rnorm(1L, sd = marginal_sd)
  for (i in 2:n) {
    noise[i] <- rho * noise[i - 1L] + stats::rnorm(1L, sd = innovation_sd)
  }
  noise
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

refine_ar1_boundaries <- function(y, boundaries, rho, innovation_variance,
                                  minimum_segment = 5L) {
  boundaries <- normalise_boundaries(boundaries, length(y))
  if (length(boundaries) <= 1L) return(boundaries)
  original <- boundaries
  previous <- c(0L, head(original, -1L))
  refined <- original
  for (j in seq_len(length(original) - 1L)) {
    left <- max(previous[j] + minimum_segment,
                floor((previous[j] + original[j]) / 2))
    right <- min(original[j + 1L] - minimum_segment,
                 ceiling((original[j] + original[j + 1L]) / 2))
    if (left > right) next
    first <- previous[j] + 1L
    last <- original[j + 1L]
    values <- y[first:last]
    costs <- ar1_cost_matrix(values, rho, innovation_variance)
    candidates <- left:right
    relative <- candidates - first + 1L
    full_cost <- costs[1L, length(values)]
    gain <- vapply(relative, function(split) {
      full_cost - costs[1L, split] - costs[split + 1L, length(values)]
    }, numeric(1))
    refined[j] <- candidates[which.max(gain)]
  }
  normalise_boundaries(refined, length(y))
}

fit_ar1_methods <- function(y, rho = 0.8, svp_constant = 3.75,
                            true_true_constant = AR1_TRUE_TRUE_CONSTANT) {
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
      ar1_approximate_partition(y, rho, innovation_variance, 3 * log(n)),
    "SVP AR1Focus" = normalise_boundaries(
      SVP(y, svp_constant * log(n), "AR1Focus",
          prune_after_if_unvalid = TRUE,
          prune_before_if_invalid = FALSE,
          rho = rho, sigma2 = innovation_variance)$changepoints,
      n
    ),
    "SVP AR1Focus Multiscale" = refine_ar1_boundaries(
      y,
      SVP(y, true_true_constant * log(n), "AR1Focus",
          prune_after_if_unvalid = TRUE,
          prune_before_if_invalid = TRUE,
          rho = rho, sigma2 = innovation_variance)$changepoints,
      rho, innovation_variance
    )
  )
}

run_ar1_power <- function(
    n = 600L, rho = 0.8, jump_sizes = seq(1, 3, 0.1), reps = 100L,
    workers = power_default_workers(), seed = 123L) {
  design <- expand.grid(pattern = POWER_PATTERNS, jump = jump_sizes,
                        rep = seq_len(reps), stringsAsFactors = FALSE)
  rows <- power_lapply(seq_len(nrow(design)), function(i) {
    d <- design[i, ]
    set_power_seed(seed + i)
    mu <- generate_power_signal(n, d$pattern, d$jump)
    y <- mu + simulate_ar1_noise(n, rho)
    truth <- normalise_boundaries(which(diff(mu) != 0), n)
    fits <- fit_ar1_methods(y, rho)
    output <- dplyr::bind_rows(lapply(names(fits), function(method) {
      boundaries <- fits[[method]]
      result_row(d$pattern, d$jump, d$rep, method, boundaries, truth,
                 mu, fitted_piecewise_mean(y, boundaries),
                 tolerance = max(5L, round(n * 0.0025)))
    }))
    output$changepoints <- unname(fits)
    output
  }, workers = workers)
  dplyr::bind_rows(rows)
}

ar1_scenario_plot <- function(n = 600L, rho = 0.8, jump = 0.8,
                              seed = 999L) {
  set.seed(seed)
  data <- dplyr::bind_rows(lapply(POWER_PATTERNS, function(pattern) {
    mu <- generate_power_signal(n, pattern, jump)
    data.frame(time = seq_len(n), value = mu + simulate_ar1_noise(n, rho),
               mean = mu, pattern = pattern)
  }))
  data$pattern <- factor(data$pattern, POWER_PATTERNS)
  ggplot2::ggplot(data, ggplot2::aes(time, value)) +
    ggplot2::geom_point(alpha = 0.2, size = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = mean), colour = "red", linewidth = 0.8) +
    ggplot2::facet_wrap(~pattern, nrow = 2L) +
    ggplot2::labs(x = "Time (Index)", y = "Value") +
    paper_theme() + ggplot2::theme(legend.position = "none")
}

run_and_save_ar1 <- function(workers = power_default_workers()) {
  results <- run_ar1_power(workers = workers)
  save_power_outputs(results, AR1_ROOT, ar1_scenario_plot(),
                     selected_jump = 1.4)
  invisible(results)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_and_save_ar1()
}
