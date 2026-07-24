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
    "SVP" = normalise_boundaries(
      SVP(y, 1.5 * log(n), "gaussian_mean")$changepoints, n
    ),
    "SVP (BIC)" = normalise_boundaries(
      SVP(y, 2 * log(n), "gaussian_mean")$changepoints, n
    ),
    "SVP FOCUS TRUE/TRUE c=1.8 refined" = refine_svp_boundaries(
      y,
      SVP(y, GAUSSIAN_TRUE_TRUE_CONSTANT * log(n), "gaussian_mean",
          prune_after_if_unvalid = TRUE,
          prune_before_if_invalid = TRUE)$changepoints,
      robust = FALSE
    )
  )
}

run_gaussian_power <- function(
    n = 1000L, jump_sizes = seq(0.1, 2, 0.1), reps = 100L,
    workers = power_default_workers(), seed = 123L) {
  design <- expand.grid(pattern = POWER_PATTERNS, jump = jump_sizes,
                        rep = seq_len(reps), stringsAsFactors = FALSE)
  rows <- power_lapply(seq_len(nrow(design)), function(i) {
    d <- design[i, ]
    set_power_seed(seed + i)
    mu <- generate_power_signal(n, d$pattern, d$jump)
    y <- mu + stats::rnorm(n)
    truth <- normalise_boundaries(which(diff(mu) != 0), n)
    fits <- fit_gaussian_methods(y)
    output <- dplyr::bind_rows(lapply(names(fits), function(method) {
      boundaries <- fits[[method]]
      result_row(d$pattern, d$jump, d$rep, method, boundaries, truth,
                 mu, fitted_piecewise_mean(y, boundaries),
                 tolerance = round(n * 0.0025))
    }))
    output$changepoints <- unname(fits)
    output
  }, workers = workers)
  dplyr::bind_rows(rows)
}

gaussian_scenario_plot <- function(n = 1000L, jump = 0.6, seed = 999L) {
  set.seed(seed)
  data <- dplyr::bind_rows(lapply(POWER_PATTERNS, function(pattern) {
    mu <- generate_power_signal(n, pattern, jump)
    data.frame(time = seq_len(n), value = mu + stats::rnorm(n),
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

gaussian_plot_results <- function(results) {
  set_algorithm_display(
    results,
    labels = c(
      "PELT" = "PELT",
      "SVP (BIC)" = "SVP (BIC)",
      "SVP" = "SVP (BIC calibrated)",
      "SVP FOCUS TRUE/TRUE c=1.8 refined" = "SVP multiscale"
    ),
    order = c("PELT", "SVP (BIC)", "SVP (BIC calibrated)",
              "SVP multiscale")
  )
}

run_and_save_gaussian <- function(workers = power_default_workers()) {
  results <- run_gaussian_power(workers = workers)
  save_power_outputs(
    results, GAUSSIAN_ROOT, gaussian_scenario_plot(),
    plot_results = gaussian_plot_results(results)
  )
  invisible(results)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_and_save_gaussian()
}
