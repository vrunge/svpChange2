## Shared, deliberately small helpers for the svpChange2 simulation studies.

simulation_path <- function(...) {
  file.path(getOption("svpChange2.simulation_root", "simulations"), ...)
}

ensure_simulation_dirs <- function() {
  dir.create(simulation_path("plots"), recursive = TRUE, showWarnings = FALSE)
}

generate_signal <- function(n,
                            pattern = c("none", "up", "updown", "rand1"),
                            nb_segments = 8,
                            jump_size = 1) {
  pattern <- match.arg(pattern)
  if (n %% nb_segments != 0 && pattern != "rand1") {
    stop("n must be divisible by nb_segments for regular scenarios")
  }

  if (pattern == "rand1") {
    ## Recreate the historical random segment lengths and means deterministically.
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(42)
    raw_sizes <- rpois(nb_segments, 10)
    sizes <- pmax(round(raw_sizes * n / sum(raw_sizes)), 1)
    sizes[nb_segments] <- sizes[nb_segments] + n - sum(sizes)
    set.seed(43)
    means <- runif(nb_segments, 0.5, 1) *
      sample(c(-1, 1), nb_segments, replace = TRUE)
    return(rep(means * jump_size, sizes))
  }

  if (pattern == "none") return(rep(0, n))

  means <- switch(
    pattern,
    up = seq(0, nb_segments - 1) * jump_size,
    updown = (seq(0, nb_segments - 1) %% 2) * jump_size
  )
  rep(means, each = n / nb_segments)
}

interior_changepoints <- function(boundaries, n) {
  setdiff(as.integer(boundaries), n)
}

cp_metrics <- function(true_boundaries, estimated_boundaries, tolerance) {
  truth <- interior_changepoints(true_boundaries, max(true_boundaries))
  estimate <- interior_changepoints(estimated_boundaries, max(true_boundaries))
  matched <- rep(FALSE, length(estimate))
  true_positive <- 0L

  for (cp in truth) {
    available <- which(!matched & abs(estimate - cp) <= tolerance)
    if (length(available)) {
      best <- available[which.min(abs(estimate[available] - cp))]
      matched[best] <- TRUE
      true_positive <- true_positive + 1L
    }
  }

  precision <- if (length(estimate)) true_positive / length(estimate) else {
    if (length(truth)) 0 else 1
  }
  recall <- if (length(truth)) true_positive / length(truth) else 1
  f1 <- if (precision + recall > 0) {
    2 * precision * recall / (precision + recall)
  } else 0

  c(Precision = precision, Recall = recall, F1 = f1)
}

piecewise_mean <- function(data, boundaries) {
  starts <- c(1L, head(boundaries, -1L) + 1L)
  unlist(Map(function(first, last) {
    rep(mean(data[first:last]), last - first + 1L)
  }, starts, boundaries), use.names = FALSE)
}

extract_boundaries <- function(result, n) {
  if (is.null(result) || is.null(result$changepoints)) return(n)
  sort(unique(c(as.integer(unlist(result$changepoints)), n)))
}

svp_method_specs <- function(q_alpha_n = 3, sigma2 = 1) {
  list(
    "SVP" = list(
      gamma_multiplier = 1.5,
      args = list(prune_after_if_unvalid = TRUE,
                  prune_before_if_invalid = FALSE)
    ),
    "SVP FOCUS (after/before)" = list(
      gamma_multiplier = 1.5,
      args = list(prune_after_if_unvalid = TRUE,
                  prune_before_if_invalid = TRUE)
    ),
    "SVP multiscale" = list(
      gamma_multiplier = 1.5,
      args = list(prune_after_if_unvalid = TRUE,
                  prune_before_if_invalid = TRUE,
                  use_multiscale_gamma = TRUE,
                  sigma2 = sigma2,
                  q_alpha_n = q_alpha_n)
    )
  )
}

run_svp_method <- function(data, method, q_alpha_n = 3, sigma2 = 1) {
  spec <- svp_method_specs(q_alpha_n, sigma2)[[method]]
  do.call(
    SVP,
    c(list(data = data,
           gamma = spec$gamma_multiplier * log(length(data)),
           test = "gaussian_mean"),
      spec$args)
  )
}

summarise_and_plot_metrics <- function(results, prefix) {
  ensure_simulation_dirs()
  long <- tidyr::pivot_longer(
    results,
    cols = c("Precision", "Recall", "F1", "NumSegments", "MSE"),
    names_to = "metric",
    values_to = "value"
  )
  summary <- dplyr::summarise(
    dplyr::group_by(long, algorithm, pattern, jump, metric),
    mean = mean(value, na.rm = TRUE),
    se = stats::sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se,
    .groups = "drop"
  )

  for (metric_name in unique(summary$metric)) {
    plot_data <- dplyr::filter(summary, metric == metric_name)
    plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = jump, y = mean, colour = algorithm, fill = algorithm)
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper),
        alpha = 0.12,
        colour = NA
      ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::facet_wrap(~pattern) +
      ggplot2::labs(x = "Jump size", y = metric_name) +
      ggplot2::theme_minimal()
    filename <- paste0(prefix, "_", tolower(metric_name), "_vs_jump.pdf")
    ggplot2::ggsave(simulation_path("plots", filename), plot,
                    width = 11, height = 6)
  }
  invisible(summary)
}

plot_cp_distributions <- function(results, prefix, selected_jump, n) {
  rows <- lapply(seq_len(nrow(results)), function(i) {
    cps <- interior_changepoints(results$changepoints[[i]], n)
    if (!length(cps)) return(NULL)
    data.frame(pattern = results$pattern[i],
               algorithm = results$algorithm[i],
               jump = results$jump[i],
               changepoint = cps)
  })
  cp_data <- dplyr::bind_rows(rows)
  cp_data <- dplyr::filter(cp_data, jump == selected_jump, pattern != "none")
  if (!nrow(cp_data)) return(invisible(NULL))

  plot <- ggplot2::ggplot(cp_data, ggplot2::aes(changepoint)) +
    ggplot2::geom_histogram(bins = 100, fill = "steelblue") +
    ggplot2::facet_grid(pattern ~ algorithm) +
    ggplot2::labs(x = "Sequence position", y = "Detected changes") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(simulation_path("plots", paste0(prefix,
                                                   "_changepoint_distributions.pdf")),
                  plot, width = 15, height = 10)
  invisible(plot)
}
