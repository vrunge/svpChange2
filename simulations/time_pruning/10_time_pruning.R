## Runtime complexity of the four SVP pruning configurations versus PELT.

library(svpChange2)
library(changepoint)
library(ggplot2)
library(dplyr)

options(svpChange2.simulation_root = file.path("simulations", "time_pruning"))
source(file.path("simulations", "0_simulation_helpers.R"))

PRUNING_STRATEGIES <- data.frame(
  method = c("SVP after=FALSE before=FALSE",
             "SVP after=TRUE before=FALSE",
             "SVP after=FALSE before=TRUE",
             "SVP after=TRUE before=TRUE"),
  prune_after_if_unvalid = c(FALSE, TRUE, FALSE, TRUE),
  prune_before_if_invalid = c(FALSE, FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

generate_signal <- function(n, k = 0, jump = 10) {
  segment_count <- k + 1L
  sizes <- rep(n %/% segment_count, segment_count)
  remainder <- n - sum(sizes)
  if (remainder > 0) sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1L
  means <- rep(c(0, jump), length.out = segment_count)
  rep(means, sizes) + rnorm(n)
}

time_one <- function(data, strategy) {
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
    length(extract_boundaries(fit, n))
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
    data <- generate_signal(n, k = 0, jump = jump)
    for (i in seq_len(nrow(strategies))) {
      row_index <- row_index + 1L
      rows[[row_index]] <- cbind(
        experiment = "vary_n", n = n, k = 0, rep = rep,
        time_one(data, strategies[i, , drop = FALSE])
      )
    }
  }

  for (k in k_values) for (rep in seq_len(reps_k)) {
    set.seed(seed + 100000 + k * 100 + rep)
    data <- generate_signal(fixed_n, k = k, jump = jump)
    for (i in seq_len(nrow(strategies))) {
      row_index <- row_index + 1L
      rows[[row_index]] <- cbind(
        experiment = "vary_k", n = fixed_n, k = k, rep = rep,
        time_one(data, strategies[i, , drop = FALSE])
      )
    }
  }

  dplyr::bind_rows(rows)
}

plot_time_pruning <- function(results) {
  ensure_simulation_dirs()
  vary_n <- dplyr::filter(results, experiment == "vary_n")
  vary_k <- dplyr::filter(results, experiment == "vary_k")

  plot_n <- ggplot(vary_n, aes(n, time, colour = method)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "Sequence length (log scale)", y = "Time (s, log scale)") +
    theme_minimal()

  plot_k <- ggplot(vary_k, aes(k, time, colour = method)) +
    geom_smooth(method = "loess", se = TRUE) +
    labs(x = "Number of changes", y = "Time (s)") +
    theme_minimal()

  plot_detected <- ggplot(results, aes(detected, time, colour = method)) +
    geom_smooth(method = "loess", se = TRUE) +
    scale_y_log10() +
    labs(x = "Detected number of segments", y = "Time (s, log scale)") +
    theme_minimal()

  ggsave(simulation_path("plots", "10_pruning_time_vs_n.pdf"),
         plot_n, width = 9, height = 5)
  ggsave(simulation_path("plots", "10_pruning_time_vs_k.pdf"),
         plot_k, width = 9, height = 5)
  ggsave(simulation_path("plots", "10_pruning_time_vs_detected_changes.pdf"),
         plot_detected, width = 9, height = 5)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  results <- run_time_pruning()
  saveRDS(results, simulation_path("10_pruning_time_results.rds"))
  write.csv(results, simulation_path("10_pruning_time_results.csv"),
            row.names = FALSE)
  plot_time_pruning(results)
}




