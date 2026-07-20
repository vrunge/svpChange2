## Runtime scaling for SVP under AR(1) noise.
library(svpChange2)
library(ggplot2)
library(dplyr)
source(file.path("simulations", "0_simulation_helpers.R"))
options(svpChange2.simulation_root = file.path("simulations", "time_ar1"))

simulate_ar1 <- function(n, k = 0, rho = 0.6, jump = 4) {
  sizes <- rep(n %/% (k + 1), k + 1)
  sizes[seq_len(n - sum(sizes))] <- sizes[seq_len(n - sum(sizes))] + 1L
  mu <- rep(rep(c(0, jump), length.out = k + 1), sizes)
  y <- numeric(n); y[1] <- mu[1] + rnorm(1) / sqrt(1 - rho^2)
  for (i in 2:n) y[i] <- mu[i] + rho * (y[i - 1] - mu[i - 1]) + rnorm(1)
  y
}

time_one <- function(y, method, rho = 0.6) {
  elapsed <- system.time({
    fit <- switch(method,
      independent = SVP(y, 2 * log(length(y)), "gaussian_mean"),
      known = SVP(y, 2 * log(length(y)), "AR1", rho = rho, sigma2 = 1),
      estimated = SVP(y, 2 * log(length(y)), "AR1", rho = AR1_rho(y), sigma2 = 1)
    )
  })[["elapsed"]]
  data.frame(method = method, time = elapsed,
             detected = length(extract_boundaries(fit, length(y))))
}

run_time_ar1 <- function(n_values = round(2^(seq(8, 13, length.out = 12))),
                           k_values = 0:40, fixed_n = 2000,
                           reps_n = 10, reps_k = 5, rho = 0.6, seed = 810) {
  methods <- c("independent", "known", "estimated")
  rows <- list(); z <- 0L
  for (n in n_values) for (rep in seq_len(reps_n)) for (method in methods) {
    set.seed(seed + n + rep)
    z <- z + 1L; rows[[z]] <- cbind(experiment = "vary_n", n = n, k = 0,
                                    rep = rep, time_one(simulate_ar1(n, rho = rho), method, rho))
  }
  for (k in k_values) for (rep in seq_len(reps_k)) for (method in methods) {
    set.seed(seed + 100000 + k * 100 + rep)
    z <- z + 1L; rows[[z]] <- cbind(experiment = "vary_k", n = fixed_n, k = k,
                                    rep = rep, time_one(simulate_ar1(fixed_n, k, rho), method, rho))
  }
  bind_rows(rows)
}

plot_time_ar1 <- function(results) {
  dir.create(simulation_path("plots"), recursive = TRUE, showWarnings = FALSE)
  p1 <- ggplot(filter(results, experiment == "vary_n"), aes(n, time, colour = method)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE) + scale_x_log10() + scale_y_log10() + theme_minimal()
  p2 <- ggplot(filter(results, experiment == "vary_k"), aes(k, time, colour = method)) +
    geom_smooth(method = "loess", se = TRUE) + theme_minimal()
  ggsave(simulation_path("plots", "8_ar1_time_vs_n.pdf"), p1, width = 8, height = 5)
  ggsave(simulation_path("plots", "8_ar1_time_vs_k.pdf"), p2, width = 8, height = 5)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  results <- run_time_ar1()
  saveRDS(results, simulation_path("8_ar1_time_results.rds"))
  write.csv(results, simulation_path("8_ar1_time_results.csv"), row.names = FALSE)
  plot_time_ar1(results)
}
