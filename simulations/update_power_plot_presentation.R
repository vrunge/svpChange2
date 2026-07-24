## Regenerate power-study plots from saved results without running simulations.

source(file.path("simulations", "power_gaussian", "run_power.R"))
source(file.path("simulations", "power_robust", "run_power.R"))
source(file.path("simulations", "power_ar1", "run_power.R"))

regenerate_power_plots <- function(root, n, scenario_plot, display,
                                   selected_jump = 0.6) {
  results <- readRDS(file.path(root, "results.rds"))
  results <- add_localization_error(results, n)
  save_power_outputs(
    results, root, scenario_plot, selected_jump = selected_jump,
    plot_results = display(results), save_results = FALSE
  )
}

regenerate_power_plots(
  GAUSSIAN_ROOT, 1000L, gaussian_scenario_plot(), gaussian_plot_results
)
regenerate_power_plots(
  ROBUST_ROOT, 1000L, robust_scenario_plot(), robust_plot_results
)
regenerate_power_plots(
  AR1_ROOT, 600L, ar1_scenario_plot(), ar1_plot_results,
  selected_jump = 1.4
)
