## Runtime complexity of the four SVP pruning configurations versus PELT.

library(svpChange2)
library(changepoint)
library(ggplot2)
library(dplyr)

options(
  svpChange2.simulation_root =
    file.path("simulations", "other_simus", "time_pruning")
)
source(file.path("simulations", "time_common.R"))

PRUNING_STRATEGIES <- data.frame(
  method = c("SVP after=FALSE before=FALSE",
             "SVP after=TRUE before=FALSE",
             "SVP after=FALSE before=TRUE",
             "SVP after=TRUE before=TRUE"),
  prune_after_if_unvalid = c(FALSE, TRUE, FALSE, TRUE),
  prune_before_if_invalid = c(FALSE, FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

generate_pruning_signal <- function(n, k = 0, jump = 10) {
  alternating_mean(n, k, jump) + rnorm(n)
}

time_pruning_strategy <- function(data, strategy) {
  n <- length(data)
  elapsed <- system.time({
    fit <- if (strategy$method == "PELT") {
      cpt.mean(data, method = "PELT", penalty = "Manual",
               pen.value = 2 * log(n))
    } else {
      SVP(data, gamma = 2 * log(n), test = "gaussian_mean",
          prune_after_if_unvalid = strategy$prune_after_if_unvalid,
          prune_before_if_invalid = strategy$prune_before_if_invalid)
    }
  })[["elapsed"]]

  detected <- if (strategy$method == "PELT") {
    length(cpts(fit)) + 1L
  } else {
    length(extract_svp_boundaries(fit, n))
  }
  data.frame(method = strategy$method, time = elapsed, detected = detected)
}

all_strategies <- function() {
  rbind(
    data.frame(method = "PELT",
               prune_after_if_unvalid = NA,
               prune_before_if_invalid = NA,
               stringsAsFactors = FALSE),
    PRUNING_STRATEGIES
  )
}

run_time_pruning <- function(
    n_values = round(2^(seq(8, 14, length.out = 16))),
    k_values = 0:100,
    fixed_n = 10000,
    reps_n = 20,
    reps_k = 10,
    jump = 10,
    seed = 1010) {
  strategies <- all_strategies()
  rows <- list()
  row_index <- 0L

  for (n in n_values) for (rep in seq_len(reps_n)) {
    set.seed(seed + n + rep)
    data <- generate_pruning_signal(n, k = 0, jump = jump)
    for (i in seq_len(nrow(strategies))) {
      row_index <- row_index + 1L
      rows[[row_index]] <- cbind(
        experiment = "vary_n", n = n, k = 0, rep = rep,
        time_pruning_strategy(data, strategies[i, , drop = FALSE])
      )
    }
  }

  for (k in k_values) for (rep in seq_len(reps_k)) {
    set.seed(seed + 100000 + k * 100 + rep)
    data <- generate_pruning_signal(fixed_n, k = k, jump = jump)
    for (i in seq_len(nrow(strategies))) {
      row_index <- row_index + 1L
      rows[[row_index]] <- cbind(
        experiment = "vary_k", n = fixed_n, k = k, rep = rep,
        time_pruning_strategy(data, strategies[i, , drop = FALSE])
      )
    }
  }

  dplyr::bind_rows(rows)
}

plot_time_pruning <- function(results) {
  plot_time_study(
    results,
    getOption("svpChange2.simulation_root"),
    "pruning_time"
  )
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  results <- run_time_pruning()
  save_time_outputs(
    results,
    getOption("svpChange2.simulation_root"),
    "pruning_time"
  )
}


