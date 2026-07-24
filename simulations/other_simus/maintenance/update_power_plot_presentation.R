## Regenerate power-study plots from saved results without rerunning simulations.

source(file.path("simulations", "power_gaussian", "run_power.R"))
source(file.path("simulations", "power_robust", "run_power.R"))
source(file.path("simulations", "power_ar1", "run_power.R"))

regenerate_power_plots <- function(root, n, scenario, order,
                                   selected_jump = 0.6,
                                   drop = character()) {
  results <- readRDS(file.path(root, "results.rds"))
  results <- add_localization_error(results, n)
  display <- set_algorithm_order(results, order, drop)
  save_power_outputs(
    results, root, scenario, selected_jump = selected_jump,
    plot_results = display, save_results = FALSE
  )
}

regenerate_power_plots(
  GAUSSIAN_ROOT,
  1000L,
  power_scenario_plot(1000L, 0.6, stats::rnorm),
  c("PELT", "SVP BIC", "SVP BIC calibrated", "SVP BIC multiscale")
)
regenerate_power_plots(
  ROBUST_ROOT,
  1000L,
  power_scenario_plot(
    1000L, 0.6, function(size) stats::rt(size, df = 2),
    y_limits = c(-10, 10)
  ),
  c("PELT", "RFPOP", "SVP MedianMood", "SVP Wilcoxon",
    "SVP Wilcoxon multiscale"),
  drop = "RFPOP paper"
)
regenerate_power_plots(
  AR1_ROOT,
  600L,
  power_scenario_plot(
    600L, 0.8, function(size) simulate_ar1_noise(size, 0.8)
  ),
  c("PELT AR1 approximate", "PELT inflated", "DeCAFS AR1",
    "SVP AR1Focus", "SVP AR1Focus multiscale"),
  selected_jump = 1.4
)
