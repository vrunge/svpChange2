## Test the former minimum-segment-length-five boundary refinement on every
## method in the three power studies. Run from the svpChange2 package root:
##   Rscript simulations/other_simus/all_method_refinement/run_refinement_test.R

.libPaths(c(normalizePath(".r-library"), .libPaths()))

library(dplyr)
library(ggplot2)

source(file.path("simulations", "power_common.R"))

TEST_ROOT <- file.path(
  "simulations", "other_simus", "all_method_refinement"
)
MIN_SEGMENT_LENGTH <- 5L

segment_candidates <- function(boundaries, index, minimum = MIN_SEGMENT_LENGTH) {
  previous <- if (index == 1L) 0L else boundaries[index - 1L]
  following <- boundaries[index + 1L]
  lower <- previous + minimum
  upper <- following - minimum
  if (lower > upper) return(integer())
  seq.int(lower, upper)
}

refine_mean_boundaries <- function(y, boundaries) {
  boundaries <- normalise_boundaries(boundaries, length(y))
  if (length(boundaries) == 1L) return(boundaries)
  refined <- boundaries
  for (index in seq_len(length(boundaries) - 1L)) {
    candidates <- segment_candidates(boundaries, index)
    if (!length(candidates)) next
    first <- if (index == 1L) 1L else boundaries[index - 1L] + 1L
    last <- boundaries[index + 1L]
    cumulative <- c(0, cumsum(y[first:last]))
    left_n <- candidates - first + 1L
    right_n <- last - candidates
    left_sum <- cumulative[left_n + 1L]
    right_sum <- cumulative[length(cumulative)] - left_sum
    score <- left_n * right_n / (left_n + right_n) *
      (left_sum / left_n - right_sum / right_n)^2
    refined[index] <- candidates[which.max(score)]
  }
  normalise_boundaries(refined, length(y))
}

refine_rank_boundaries <- function(y, boundaries) {
  boundaries <- normalise_boundaries(boundaries, length(y))
  if (length(boundaries) == 1L) return(boundaries)
  refined <- boundaries
  for (index in seq_len(length(boundaries) - 1L)) {
    candidates <- segment_candidates(boundaries, index)
    if (!length(candidates)) next
    first <- if (index == 1L) 1L else boundaries[index - 1L] + 1L
    last <- boundaries[index + 1L]
    ranks <- rank(y[first:last], ties.method = "average")
    centered <- cumsum(ranks - mean(ranks))
    offsets <- candidates - first + 1L
    scale <- sqrt(offsets * (length(ranks) - offsets))
    refined[index] <- candidates[which.max(abs(centered[offsets]) / scale)]
  }
  normalise_boundaries(refined, length(y))
}

ar1_segment_cost <- function(y, first, last, rho, innovation_variance) {
  if (first == last) return(0)
  transition <- y[(first + 1L):last] - rho * y[first:(last - 1L)]
  first_weight <- 1 - rho^2
  mean_weight <- 1 - rho
  denominator <- first_weight + length(transition) * mean_weight^2
  numerator <- first_weight * y[first] + mean_weight * sum(transition)
  constant <- first_weight * y[first]^2 + sum(transition^2)
  max(0, constant - numerator^2 / denominator) / innovation_variance
}

refine_ar1_boundaries <- function(y, boundaries, rho) {
  boundaries <- normalise_boundaries(boundaries, length(y))
  if (length(boundaries) == 1L) return(boundaries)
  innovation_variance <- 1 - rho^2
  refined <- boundaries
  for (index in seq_len(length(boundaries) - 1L)) {
    candidates <- segment_candidates(boundaries, index)
    if (!length(candidates)) next
    first <- if (index == 1L) 1L else boundaries[index - 1L] + 1L
    last <- boundaries[index + 1L]
    cost <- vapply(candidates, function(split) {
      ar1_segment_cost(y, first, split, rho, innovation_variance) +
        ar1_segment_cost(y, split + 1L, last, rho, innovation_variance)
    }, numeric(1))
    refined[index] <- candidates[which.min(cost)]
  }
  normalise_boundaries(refined, length(y))
}

recreate_noise <- function(framework, n, rho) {
  switch(
    framework,
    gaussian = stats::rnorm(n),
    robust = stats::rt(n, df = 2),
    ar1 = simulate_ar1_noise(n, rho),
    stop("Unknown framework: ", framework)
  )
}

refine_saved_study <- function(framework, source_root, n, tolerance,
                               rho = NA_real_, workers = 1L, seed = 123L) {
  input <- readRDS(file.path(source_root, "results.rds"))
  design <- expand.grid(
    pattern = POWER_PATTERNS,
    jump = sort(unique(input$jump)),
    rep = sort(unique(input$rep)),
    stringsAsFactors = FALSE
  )
  rows <- power_lapply(seq_len(nrow(design)), function(i) {
    task <- design[i, ]
    selected <- input[
      input$pattern == task$pattern & input$jump == task$jump &
        input$rep == task$rep,
    ]
    set_power_seed(seed + i)
    mu <- generate_power_signal(n, task$pattern, task$jump)
    y <- mu + recreate_noise(framework, n, rho)
    truth <- normalise_boundaries(which(diff(mu) != 0), n)
    refined <- lapply(selected$changepoints, function(boundaries) {
      switch(
        framework,
        gaussian = refine_mean_boundaries(y, boundaries),
        robust = refine_rank_boundaries(y, boundaries),
        ar1 = refine_ar1_boundaries(y, boundaries, rho)
      )
    })
    output <- dplyr::bind_rows(Map(function(method, boundaries) {
      result_row(
        task$pattern, task$jump, task$rep, as.character(method),
        boundaries, truth, mu, fitted_piecewise_mean(y, boundaries), tolerance
      )
    }, selected$algorithm, refined))
    output$changepoints <- unname(refined)
    output
  }, workers = workers)
  dplyr::bind_rows(rows)
}

save_refinement_study <- function(name, results, n, scenario_plot,
                                  algorithm_order, selected_jump) {
  root <- file.path(TEST_ROOT, name)
  save_power_outputs(
    results, root, scenario_plot, selected_jump = selected_jump,
    plot_results = set_algorithm_order(results, algorithm_order)
  )
}

run_refinement_test <- function(workers = power_default_workers()) {
  source(file.path("simulations", "power_ar1", "run_power.R"))

  cat("Refining Gaussian methods...\n")
  gaussian <- refine_saved_study(
    "gaussian", file.path("simulations", "power_gaussian"),
    n = 1000L, tolerance = round(1000L * 0.0025), workers = workers
  )
  cat("Refining robust methods...\n")
  robust <- refine_saved_study(
    "robust", file.path("simulations", "power_robust"),
    n = 1000L, tolerance = round(1000L * 0.0025), workers = workers
  )
  cat("Refining AR1 methods...\n")
  ar1 <- refine_saved_study(
    "ar1", file.path("simulations", "power_ar1"),
    n = 600L, tolerance = max(5L, round(600L * 0.0025)),
    rho = 0.8, workers = workers
  )
  studies <- list(gaussian = gaussian, robust = robust, ar1 = ar1)

  cat("Saving test results and figures...\n")
  save_refinement_study(
    "power_gaussian", studies$gaussian, 1000L,
    power_scenario_plot(1000L, 0.6, stats::rnorm),
    c("PELT", "SVP BIC", "SVP BIC calibrated", "SVP BIC multiscale"),
    0.6
  )
  save_refinement_study(
    "power_robust", studies$robust, 1000L,
    power_scenario_plot(
      1000L, 0.6, function(size) stats::rt(size, df = 2),
      y_limits = c(-10, 10)
    ),
    c("PELT", "RFPOP paper", "RFPOP", "SVP MedianMood", "SVP Wilcoxon",
      "SVP Wilcoxon multiscale"),
    0.6
  )
  save_refinement_study(
    "power_ar1", studies$ar1, 600L,
    power_scenario_plot(
      600L, 0.8, function(size) simulate_ar1_noise(size, 0.8)
    ),
    c("PELT AR1 approximate", "PELT inflated", "DeCAFS AR1",
      "SVP AR1Focus", "SVP AR1Focus multiscale"),
    1.4
  )
  cat("Done: ", TEST_ROOT, "\n", sep = "")
  invisible(studies)
}

run_refinement_test()
