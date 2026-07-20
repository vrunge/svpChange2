## Complexity plots for the baseline and two new SVP configurations.

library(svpChange2)
options(svpChange2.simulation_root = file.path("simulations", "time_gaussian"))
library(ggplot2)
library(dplyr)
library(tidyr)

source(file.path("simulations", "0_simulation_helpers.R"))

time_svp_method <- function(data, method, q_alpha_n = 3) {
  fit <- NULL
  elapsed <- system.time({
    fit <- if (method == "PELT") {
      changepoint::cpt.mean(data, method = "PELT", penalty = "Manual",
                            pen.value = 2 * log(length(data)))
    } else {
      run_svp_method(data, method, q_alpha_n, stats::var(data))
    }
  })[["elapsed"]]
  boundaries <- if (method == "PELT") {
    c(changepoint::cpts(fit), length(data))
  } else {
    extract_boundaries(fit, length(data))
  }
  data.frame(method = method, time = elapsed,
             detected = length(unique(boundaries)))
}

run_new_methods_complexity <- function(
    n_values = round(2^(seq(8, 14, length.out = 20))),
    k_values = 0:100,
    fixed_n = 10000,
    reps_n = 20,
    reps_k = 10,
    gap = 10,
    q_alpha_n = 3,
    seed = 456) {
  if (!requireNamespace("changepoint", quietly = TRUE)) {
    stop("the optional package 'changepoint' is required")
  }
  methods <- c("PELT", names(svp_method_specs(q_alpha_n)))
  rows_n <- list()
  index <- 0L
  for (n in n_values) for (rep in seq_len(reps_n)) {
    set.seed(seed + n + rep)
    data <- stats::rnorm(n)
    for (method in methods) {
      index <- index + 1L
      rows_n[[index]] <- cbind(experiment = "vary_n", n = n, k = 0,
                               rep = rep,
                               time_svp_method(data, method, q_alpha_n))
    }
  }

  rows_k <- list()
  index <- 0L
  for (k in k_values) for (rep in seq_len(reps_k)) {
    set.seed(seed + 100000 + k * 100 + rep)
    segment_count <- k + 1L
    sizes <- rep(fixed_n %/% segment_count, segment_count)
    sizes[seq_len(fixed_n - sum(sizes))] <-
      sizes[seq_len(fixed_n - sum(sizes))] + 1L
    means <- rep(c(0, gap), length.out = segment_count)
    data <- stats::rnorm(fixed_n, rep(means, sizes), 1)
    for (method in methods) {
      index <- index + 1L
      rows_k[[index]] <- cbind(experiment = "vary_k", n = fixed_n, k = k,
                               rep = rep,
                               time_svp_method(data, method, q_alpha_n))
    }
  }
  dplyr::bind_rows(rows_n, rows_k)
}

plot_new_methods_complexity <- function(results) {
  ensure_simulation_dirs()
  vary_n <- dplyr::filter(results, experiment == "vary_n")
  plot_n <- ggplot2::ggplot(vary_n, ggplot2::aes(n, time, colour = method)) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Sequence length (log scale)", y = "Time (s, log scale)") +
    ggplot2::theme_minimal()

  vary_k <- dplyr::filter(results, experiment == "vary_k")
  plot_k <- ggplot2::ggplot(vary_k, ggplot2::aes(k, time, colour = method)) +
    ggplot2::geom_smooth(method = "loess", se = TRUE) +
    ggplot2::labs(x = "True number of changes", y = "Time (s)") +
    ggplot2::theme_minimal()

  plot_detected <- ggplot2::ggplot(
    results,
    ggplot2::aes(detected, time, colour = method)
  ) + ggplot2::geom_smooth(method = "loess", se = TRUE) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Detected number of segments", y = "Time (s, log scale)") +
    ggplot2::theme_minimal()

  ggplot2::ggsave(simulation_path("plots", "7_new_methods_time_vs_n.pdf"),
                  plot_n, width = 8, height = 5)
  ggplot2::ggsave(simulation_path("plots", "7_new_methods_time_vs_k.pdf"),
                  plot_k, width = 8, height = 5)
  ggplot2::ggsave(simulation_path("plots",
                                  "7_new_methods_time_vs_detected_changes.pdf"),
                  plot_detected, width = 8, height = 5)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  complexity_results <- run_new_methods_complexity()
  saveRDS(complexity_results,
          simulation_path("7_new_methods_complexity_results.rds"))
  utils::write.csv(complexity_results,
                   simulation_path("7_new_methods_complexity_results.csv"),
                   row.names = FALSE)
  plot_new_methods_complexity(complexity_results)
}
