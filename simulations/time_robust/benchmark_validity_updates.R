## Empirical scaling of the Section 4.2 robust validity updates.
##
## The complete update used by SVP is update(y_t) followed by statistic().
## Linear per-update work implies quadratic cumulative stream time and a
## log-log slope of one for total_time / n versus n.

library(Rcpp)
library(dplyr)
library(ggplot2)

ROOT <- file.path("simulations", "time_robust")
Rcpp::sourceCpp(file.path(ROOT, "benchmark_validity_updates.cpp"))

run_validity_update_benchmark <- function(
    lengths = as.integer(2^(7:14)), trials = 7L, repetitions = 3L,
    seed = 20260723L) {
  design <- expand.grid(
    n = lengths,
    method = c("Wilcoxon", "MedianMood"),
    complete_update = c(TRUE, FALSE),
    trial = seq_len(trials),
    stringsAsFactors = FALSE
  )
  rows <- lapply(seq_len(nrow(design)), function(i) {
    d <- design[i, ]
    set.seed(seed + d$n + d$trial)
    data <- rt(d$n, df = 2)
    timing <- benchmark_robust_validity_stream(
      data, d$method, repetitions, d$complete_update
    )
    data.frame(
      n = d$n,
      method = d$method,
      operation = if (d$complete_update) "update + statistic" else "update only",
      trial = d$trial,
      total_seconds = timing$seconds,
      seconds_per_update = timing$seconds / d$n
    )
  })
  bind_rows(rows)
}

summarise_validity_update_benchmark <- function(results) {
  timing <- results |>
    group_by(n, method, operation) |>
    summarise(
      total_seconds = median(total_seconds),
      seconds_per_update = median(seconds_per_update),
      .groups = "drop"
    )
  slopes <- timing |>
    group_by(method, operation) |>
    group_modify(~{
      total_fit <- lm(log(total_seconds) ~ log(n), data = .x)
      update_fit <- lm(log(seconds_per_update) ~ log(n), data = .x)
      largest <- tail(arrange(.x, n), 4L)
      asymptotic_fit <- lm(log(seconds_per_update) ~ log(n), data = largest)
      data.frame(
        total_time_slope = unname(coef(total_fit)[2]),
        per_update_slope = unname(coef(update_fit)[2]),
        largest_four_per_update_slope = unname(coef(asymptotic_fit)[2]),
        total_r_squared = summary(total_fit)$r.squared,
        update_r_squared = summary(update_fit)$r.squared,
        largest_four_r_squared = summary(asymptotic_fit)$r.squared
      )
    }) |>
    ungroup()
  list(timing = timing, slopes = slopes)
}

plot_validity_update_benchmark <- function(timing) {
  ggplot(timing,
         aes(n, seconds_per_update, colour = method, shape = operation)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Segment length", y = "Median seconds per streamed update",
         colour = "Validity test", shape = "Measured operation") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  results <- run_validity_update_benchmark()
  summary <- summarise_validity_update_benchmark(results)
  write.csv(results, file.path(ROOT, "validity_update_benchmark_raw.csv"),
            row.names = FALSE)
  write.csv(summary$timing,
            file.path(ROOT, "validity_update_benchmark_summary.csv"),
            row.names = FALSE)
  write.csv(summary$slopes,
            file.path(ROOT, "validity_update_benchmark_slopes.csv"),
            row.names = FALSE)
  ggsave(
    file.path(ROOT, "validity_update_benchmark.pdf"),
    plot_validity_update_benchmark(summary$timing),
    width = 8, height = 5
  )
}
