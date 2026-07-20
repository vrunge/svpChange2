## Power study for piecewise-constant marginal means under AR(1) noise.

library(svpChange2)
options(svpChange2.simulation_root = file.path("simulations", "power_ar1"))
library(ggplot2)
library(dplyr)
library(tidyr)

source(file.path("simulations", "0_simulation_helpers.R"))

simulate_piecewise_ar1 <- function(mu, rho, innovation_sd = 1) {
  if (!is.numeric(rho) || length(rho) != 1L || abs(rho) >= 1) {
    stop("rho must be a single number strictly between -1 and 1")
  }
  n <- length(mu)
  innovations <- stats::rnorm(n, sd = innovation_sd)
  data <- numeric(n)
  data[1] <- mu[1] + innovations[1] / sqrt(1 - rho^2)
  for (time in 2:n) {
    data[time] <- mu[time] +
      rho * (data[time - 1L] - mu[time - 1L]) + innovations[time]
  }
  data
}

run_ar1_study <- function(n = 1000,
                          patterns = c("none", "up", "updown", "rand1"),
                          jump_sizes = seq(0.1, 2, 0.1),
                          rho_values = c(0.3, 0.6, 0.9),
                          reps = 100,
                          nb_segments = 8,
                          innovation_sd = 1,
                          gamma_multiplier = 1.5,
                          seed = 321) {
  methods <- c("SVP independent", "SVP AR1 (known rho)",
               "SVP AR1 (estimated rho)")
  design <- expand.grid(pattern = patterns, jump = jump_sizes,
                        rho = rho_values, rep = seq_len(reps),
                        stringsAsFactors = FALSE)
  rows <- vector("list", nrow(design) * length(methods))
  row_index <- 0L

  for (i in seq_len(nrow(design))) {
    set.seed(seed + i)
    mu <- generate_signal(n, design$pattern[i], nb_segments, design$jump[i])
    data <- simulate_piecewise_ar1(mu, design$rho[i], innovation_sd)
    truth <- c(which(diff(mu) != 0), n)
    rho_estimate <- AR1_rho(data)

    fits <- list(
      "SVP independent" = SVP(data, gamma_multiplier * log(n),
                               "gaussian_mean"),
      "SVP AR1 (known rho)" = SVP(
        data, gamma_multiplier * log(n), "AR1",
        prune_after_if_unvalid = TRUE,
        prune_before_if_invalid = TRUE,
        rho = design$rho[i], sigma2 = innovation_sd^2
      ),
      "SVP AR1 (estimated rho)" = SVP(
        data, gamma_multiplier * log(n), "AR1",
        prune_after_if_unvalid = TRUE,
        prune_before_if_invalid = TRUE,
        rho = rho_estimate, sigma2 = innovation_sd^2
      )
    )

    for (method in methods) {
      row_index <- row_index + 1L
      boundaries <- extract_boundaries(fits[[method]], n)
      metrics <- cp_metrics(truth, boundaries, round(n * 0.0025))
      estimate <- piecewise_mean(data, boundaries)
      rows[[row_index]] <- data.frame(
        pattern = design$pattern[i], jump = design$jump[i],
        rho = design$rho[i], rho_estimate = rho_estimate,
        rep = design$rep[i], algorithm = method,
        Precision = metrics[["Precision"]],
        Recall = metrics[["Recall"]], F1 = metrics[["F1"]],
        NumSegments = length(boundaries), MSE = mean((mu - estimate)^2)
      )
      rows[[row_index]]$changepoints <- list(boundaries)
    }
  }
  dplyr::bind_rows(rows)
}

plot_ar1_metrics <- function(results, stem = "6_ar1") {
  ## Reuse the paper metric plots, with one plot set for each autocorrelation.
  for (rho_value in sort(unique(results$rho))) {
    subset <- dplyr::filter(results, rho == rho_value)
    suffix <- gsub("\\.", "_", format(rho_value, trim = TRUE))
    summarise_and_plot_metrics(subset, paste0(stem, "_rho_", suffix))
    plot_cp_distributions(subset, paste0(stem, "_rho_", suffix),
                          selected_jump = 0.6,
                          n = max(unlist(subset$changepoints)))
  }
}

plot_ar1_scenarios <- function(n = 1000, rho_values = c(0.3, 0.6, 0.9),
                               jump_size = 0.6, nb_segments = 8) {
  scenarios <- lapply(seq_along(rho_values), function(i) {
    set.seed(900 + i)
    mu <- generate_signal(n, "updown", nb_segments, jump_size)
    data.frame(time = seq_len(n), data = simulate_piecewise_ar1(mu, rho_values[i]),
               mean = mu, rho = factor(rho_values[i]))
  })
  scenarios <- dplyr::bind_rows(scenarios)
  plot <- ggplot2::ggplot(scenarios, ggplot2::aes(time, data)) +
    ggplot2::geom_point(alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = mean), colour = "red") +
    ggplot2::facet_wrap(~rho) + ggplot2::theme_minimal()
  ggplot2::ggsave(simulation_path("plots", "6_ar1_signal_scenarios.pdf"),
                  plot, width = 11, height = 5)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  ar1_results <- run_ar1_study()
  saveRDS(ar1_results, simulation_path("6_ar1_results.rds"))
  utils::write.csv(dplyr::select(ar1_results, -changepoints),
                   simulation_path("6_ar1_results.csv"), row.names = FALSE)
  plot_ar1_metrics(ar1_results)
  plot_ar1_scenarios()
}
