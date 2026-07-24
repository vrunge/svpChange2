## Shared infrastructure for the three runtime studies.

balanced_segment_sizes <- function(n, changes) {
  segments <- changes + 1L
  sizes <- rep(n %/% segments, segments)
  remainder <- n - sum(sizes)
  if (remainder > 0L) {
    sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1L
  }
  sizes
}

alternating_mean <- function(n, changes = 0L, jump = 4) {
  sizes <- balanced_segment_sizes(n, changes)
  rep(rep(c(0, jump), length.out = changes + 1L), sizes)
}

extract_svp_boundaries <- function(fit, n) {
  if (is.null(fit) || is.null(fit$changepoints)) return(n)
  sort(unique(c(as.integer(unlist(fit$changepoints)), n)))
}

benchmark_time_study <- function(
    n_values, k_values, fixed_n, reps_n, reps_k, methods,
    simulate_data, fit_method, seed) {
  total <- length(methods) *
    (length(n_values) * reps_n + length(k_values) * reps_k)
  if (total == 0L) return(data.frame())

  progress <- utils::txtProgressBar(min = 0, max = total, style = 3)
  on.exit(close(progress), add = TRUE)
  rows <- vector("list", total)
  completed <- 0L

  add_tasks <- function(experiment, values, reps, n_for, k_for, seed_for) {
    for (value in values) {
      for (rep in seq_len(reps)) {
        n <- n_for(value)
        k <- k_for(value)
        set.seed(seed_for(value, rep))
        data <- simulate_data(n, k)
        for (method in methods) {
          elapsed <- system.time({
            boundaries <- fit_method(data, method)
          })[["elapsed"]]
          completed <<- completed + 1L
          rows[[completed]] <<- data.frame(
            experiment = experiment,
            n = n,
            k = k,
            rep = rep,
            method = method,
            time = elapsed,
            detected = length(boundaries)
          )
          utils::setTxtProgressBar(progress, completed)
        }
      }
    }
  }

  add_tasks(
    "vary_n", n_values, reps_n,
    n_for = function(value) value,
    k_for = function(value) 0L,
    seed_for = function(value, rep) seed + value + rep
  )
  add_tasks(
    "vary_k", k_values, reps_k,
    n_for = function(value) fixed_n,
    k_for = function(value) value,
    seed_for = function(value, rep) seed + 100000L + value * 100L + rep
  )
  dplyr::bind_rows(rows)
}

plot_time_study <- function(results, root, prefix) {
  dir.create(file.path(root, "plots"), recursive = TRUE, showWarnings = FALSE)
  vary_n <- dplyr::filter(results, experiment == "vary_n")
  vary_k <- dplyr::filter(results, experiment == "vary_k")
  positive_vary_n <- dplyr::filter(vary_n, is.finite(time), time > 0)
  positive_results <- dplyr::filter(results, is.finite(time), time > 0)

  time_vs_n <- ggplot2::ggplot(
    positive_vary_n, ggplot2::aes(n, time, colour = method)
  ) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Sequence length (log scale)", y = "Time (s, log scale)"
    ) +
    ggplot2::theme_minimal()

  time_vs_k <- ggplot2::ggplot(
    vary_k, ggplot2::aes(k, time, colour = method)
  ) +
    ggplot2::geom_smooth(method = "loess", se = TRUE) +
    ggplot2::labs(x = "True number of changes", y = "Time (s)") +
    ggplot2::theme_minimal()

  time_vs_detected <- ggplot2::ggplot(
    positive_results, ggplot2::aes(detected, time, colour = method)
  ) +
    ggplot2::geom_smooth(method = "loess", se = TRUE) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Detected number of segments", y = "Time (s, log scale)"
    ) +
    ggplot2::theme_minimal()

  plots <- list(
    time_vs_n = time_vs_n,
    time_vs_k = time_vs_k,
    time_vs_detected = time_vs_detected
  )
  for (name in names(plots)) {
    ggplot2::ggsave(
      file.path(root, "plots", paste0(prefix, "_", name, ".pdf")),
      plots[[name]], width = 8, height = 5
    )
  }
  invisible(plots)
}

save_time_outputs <- function(results, root, prefix) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, file.path(root, paste0(prefix, "_results.rds")))
  utils::write.csv(
    results, file.path(root, paste0(prefix, "_results.csv")),
    row.names = FALSE
  )
  plot_time_study(results, root, prefix)
  invisible(results)
}
