## Runtime scaling for robust and heavy-tail SVP tests.
library(svpChange2)
library(changepoint)
library(ggplot2)
library(dplyr)
source(file.path("simulations", "0_simulation_helpers.R"))
options(svpChange2.simulation_root = file.path("simulations", "time_robust"))

simulate_heavy_tail <- function(n, k = 0, jump = 4) {
  sizes <- rep(n %/% (k + 1), k + 1)
  sizes[seq_len(n - sum(sizes))] <- sizes[seq_len(n - sum(sizes))] + 1L
  rep(rep(c(0, jump), length.out = k + 1), sizes) + rt(n, df = 2)
}

time_one <- function(y, method) {
  n <- length(y)
  elapsed <- system.time({
    fit <- switch(method,
      PELT = changepoint::cpt.mean(y, method = "PELT", penalty = "Manual", pen.value = 2 * log(n)),
      Wilcoxon = SVP(y, 2 * log(n), "WilcoxonCost"),
      MedianMood = SVP(y, 2 * log(n), "MedianMoodCost"),
      Focus = SVP(y, 1.5 * log(n), "gaussian_mean", TRUE, TRUE),
      Multiscale = SVP(y, 0, "gaussian_mean", TRUE, TRUE, TRUE, mad(y)^2, 3)
    )
  })[["elapsed"]]
  data.frame(method = method, time = elapsed,
             detected = if (method == "PELT") length(cpts(fit)) + 1 else length(extract_boundaries(fit, n)))
}

run_time_robust <- function(n_values = round(2^(seq(8, 13, length.out = 12))),
                              k_values = 0:40, fixed_n = 2000,
                              reps_n = 10, reps_k = 5, seed = 910) {
  methods <- c("PELT", "Wilcoxon", "MedianMood", "Focus", "Multiscale")
  rows <- list(); z <- 0L
  for (n in n_values) for (rep in seq_len(reps_n)) for (method in methods) {
    set.seed(seed + n + rep)
    z <- z + 1L; rows[[z]] <- cbind(experiment = "vary_n", n = n, k = 0,
                                    rep = rep, time_one(simulate_heavy_tail(n), method))
  }
  for (k in k_values) for (rep in seq_len(reps_k)) for (method in methods) {
    set.seed(seed + 100000 + k * 100 + rep)
    z <- z + 1L; rows[[z]] <- cbind(experiment = "vary_k", n = fixed_n, k = k,
                                    rep = rep, time_one(simulate_heavy_tail(fixed_n, k), method))
  }
  bind_rows(rows)
}

plot_time_robust <- function(results) {
  dir.create(simulation_path("plots"), recursive = TRUE, showWarnings = FALSE)
  p1 <- ggplot(filter(results, experiment == "vary_n"), aes(n, time, colour = method)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE) + scale_x_log10() + scale_y_log10() + theme_minimal()
  p2 <- ggplot(filter(results, experiment == "vary_k"), aes(k, time, colour = method)) +
    geom_smooth(method = "loess", se = TRUE) + theme_minimal()
  ggsave(simulation_path("plots", "9_robust_time_vs_n.pdf"), p1, width = 8, height = 5)
  ggsave(simulation_path("plots", "9_robust_time_vs_k.pdf"), p2, width = 8, height = 5)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  results <- run_time_robust()
  saveRDS(results, simulation_path("9_robust_time_results.rds"))
  write.csv(results, simulation_path("9_robust_time_results.csv"), row.names = FALSE)
  plot_time_robust(results)
}
