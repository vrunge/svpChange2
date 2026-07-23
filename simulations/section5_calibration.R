## Calibration utilities for the extended Section 5 simulations.
##
## Tuning is deliberately separated from the final Monte Carlo samples.  The
## selected multiplier maximizes calibration-set power subject to the requested
## null control.  The final curves use separate random seeds.

select_null_controlled_multiplier <- function(calibration_results,
                                              target_null_f1 = 0.97) {
  required <- c("multiplier", "pattern", "F1")
  if (!all(required %in% names(calibration_results))) {
    stop("calibration_results must contain: ", paste(required, collapse = ", "))
  }
  score <- calibration_results |>
    dplyr::group_by(multiplier) |>
    dplyr::summarise(
      null_f1 = mean(F1[pattern == "none"], na.rm = TRUE),
      signal_f1 = mean(F1[pattern != "none"], na.rm = TRUE),
      .groups = "drop"
    )
  admissible <- dplyr::filter(score, null_f1 >= target_null_f1)
  if (!nrow(admissible)) {
    warning("No multiplier attained target null F1; selecting best null control")
    admissible <- dplyr::filter(score, null_f1 == max(null_f1))
  }
  admissible |>
    dplyr::arrange(dplyr::desc(signal_f1), multiplier) |>
    dplyr::slice(1L)
}

paper_plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "plain")
    )
}

summarise_paper_metric <- function(results, metric = "F1") {
  results |>
    dplyr::group_by(algorithm, pattern, jump) |>
    dplyr::summarise(
      mean = mean(.data[[metric]], na.rm = TRUE),
      se = stats::sd(.data[[metric]], na.rm = TRUE) /
        sqrt(sum(is.finite(.data[[metric]]))),
      lower = pmax(0, mean - 1.96 * se),
      upper = pmin(1, mean + 1.96 * se),
      .groups = "drop"
    )
}

plot_paper_power <- function(results, metric = "F1") {
  summary <- summarise_paper_metric(results, metric)
  summary$pattern <- factor(summary$pattern,
                            levels = c("none", "up", "updown", "rand1"))
  ggplot2::ggplot(
    summary,
    ggplot2::aes(jump, mean, colour = algorithm, fill = algorithm)
  ) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         alpha = 0.10, colour = NA) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~pattern, nrow = 2L) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Jump size", y = metric, colour = "algorithm",
                  fill = "algorithm") +
    paper_plot_theme()
}
