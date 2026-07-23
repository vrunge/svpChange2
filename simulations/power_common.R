## Shared infrastructure for the three paper power studies.

POWER_PATTERNS <- c("none", "up", "updown", "rand1")

power_default_workers <- function() {
  detected <- parallel::detectCores(logical = FALSE)
  if (is.na(detected)) return(1L)
  max(1L, as.integer(detected) - 1L)
}

power_lapply <- function(X, FUN, workers = power_default_workers(), ...) {
  workers <- as.integer(workers)
  if (length(workers) != 1L || is.na(workers) || workers < 1L) {
    stop("workers must be one positive integer")
  }
  if (workers == 1L || .Platform$OS.type == "windows") {
    return(lapply(X, FUN, ...))
  }
  parallel::mclapply(
    X, FUN, ..., mc.cores = workers, mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
}

set_power_seed <- function(seed) {
  set.seed(seed, kind = "Mersenne-Twister",
           normal.kind = "Inversion", sample.kind = "Rejection")
}

generate_power_signal <- function(n, pattern, jump, segments = 8L) {
  pattern <- match.arg(pattern, POWER_PATTERNS)
  if (pattern == "none") return(rep(0, n))
  if (pattern == "rand1") {
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
    raw_sizes <- stats::rpois(segments, 10)
    sizes <- pmax(round(raw_sizes * n / sum(raw_sizes)), 1)
    sizes[segments] <- sizes[segments] + n - sum(sizes)
    set.seed(43)
    levels <- stats::runif(segments, 0.5, 1) *
      sample(c(-1, 1), segments, replace = TRUE)
    return(rep(levels * jump, sizes))
  }
  if (n %% segments != 0) stop("n must be divisible by segments")
  levels <- if (pattern == "up") {
    seq(0, segments - 1) * jump
  } else {
    (seq(0, segments - 1) %% 2) * jump
  }
  rep(levels, each = n / segments)
}

normalise_boundaries <- function(boundaries, n) {
  sort(unique(c(as.integer(unlist(boundaries)), n)))
}

refine_svp_boundaries <- function(y, boundaries, robust = FALSE,
                                  minimum_segment = 5L) {
  boundaries <- normalise_boundaries(boundaries, length(y))
  if (length(boundaries) <= 1L) return(boundaries)
  original <- boundaries
  previous <- c(0L, head(original, -1L))
  refined <- original

  for (j in seq_len(length(original) - 1L)) {
    left <- max(previous[j] + minimum_segment,
                floor((previous[j] + original[j]) / 2))
    right <- min(original[j + 1L] - minimum_segment,
                 ceiling((original[j] + original[j + 1L]) / 2))
    if (left > right) next
    first <- previous[j] + 1L
    last <- original[j + 1L]
    values <- y[first:last]
    if (robust) {
      values <- rank(values, ties.method = "average")
    }
    sums <- c(0, cumsum(values))
    candidates <- left:right
    relative <- candidates - first + 1L
    left_n <- relative
    right_n <- length(values) - relative
    left_sum <- sums[relative + 1L]
    right_sum <- sums[length(sums)] - left_sum
    gain <- left_sum^2 / left_n + right_sum^2 / right_n
    refined[j] <- candidates[which.max(gain)]
  }
  normalise_boundaries(refined, length(y))
}

paper_metrics <- function(true_boundaries, estimated_boundaries, tolerance) {
  truth <- as.integer(true_boundaries)
  estimate <- as.integer(estimated_boundaries)
  distances <- abs(outer(truth, estimate, "-"))
  used <- rep(FALSE, length(estimate))
  matched <- 0L
  for (i in seq_along(truth)) {
    available <- which(!used & distances[i, ] <= tolerance)
    if (length(available)) {
      best <- available[which.min(distances[i, available])]
      used[best] <- TRUE
      matched <- matched + 1L
    }
  }
  precision <- matched / length(estimate)
  recall <- matched / length(truth)
  f1 <- if (precision + recall) 2 * precision * recall / (precision + recall) else 0
  c(Precision = precision, Recall = recall, F1 = f1,
    CorrectNumCP = as.numeric(length(estimate) == length(truth)))
}

fitted_piecewise_mean <- function(y, boundaries) {
  starts <- c(1L, head(boundaries, -1L) + 1L)
  unlist(Map(function(first, last) {
    rep(mean(y[first:last]), last - first + 1L)
  }, starts, boundaries), use.names = FALSE)
}

result_row <- function(pattern, jump, rep, algorithm, boundaries,
                       truth, mu, fitted, tolerance) {
  metric <- paper_metrics(truth, boundaries, tolerance)
  data.frame(
    pattern = pattern, jump = jump, rep = rep, algorithm = algorithm,
    Precision = metric[["Precision"]], Recall = metric[["Recall"]],
    F1 = metric[["F1"]], CorrectNumCP = metric[["CorrectNumCP"]],
    MSE = mean((mu - fitted)^2), NumSegments = length(boundaries),
    stringsAsFactors = FALSE
  )
}

paper_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "plain"),
      legend.position = "right"
    )
}

summarise_metric <- function(results, metric) {
  summary <- results |>
    dplyr::group_by(algorithm, pattern, jump) |>
    dplyr::summarise(
      mean = mean(.data[[metric]], na.rm = TRUE),
      se = stats::sd(.data[[metric]], na.rm = TRUE) /
        sqrt(sum(is.finite(.data[[metric]]))),
      lower = mean - 1.96 * se, upper = mean + 1.96 * se,
      .groups = "drop"
    )
  if (metric %in% c("Precision", "Recall", "F1", "CorrectNumCP")) {
    summary$lower <- pmax(0, summary$lower)
    summary$upper <- pmin(1, summary$upper)
  } else if (metric == "MSE") {
    summary$lower <- pmax(0, summary$lower)
  }
  summary
}

metric_plot <- function(results, metric, y_label = metric,
                        probability_scale = FALSE) {
  data <- summarise_metric(results, metric)
  data$pattern <- factor(data$pattern, POWER_PATTERNS)
  plot <- ggplot2::ggplot(
    data, ggplot2::aes(jump, mean, colour = algorithm, fill = algorithm)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.1) +
    ggplot2::facet_wrap(~pattern, nrow = 2L) +
    ggplot2::labs(x = "Jump size", y = y_label,
                  colour = "algorithm", fill = "algorithm") +
    paper_theme()
  if (probability_scale) plot <- plot + ggplot2::coord_cartesian(ylim = c(0, 1))
  plot
}

precision_recall_plot <- function(results) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork is required for the combined precision/recall figure")
  }
  precision <- metric_plot(results, "Precision", "Precision", TRUE) +
    ggplot2::theme(legend.position = "none")
  recall <- metric_plot(results, "Recall", "Recall", TRUE)
  patchwork::wrap_plots(precision, recall, ncol = 1L,
                        guides = "collect") +
    patchwork::plot_annotation(tag_levels = "a") &
    ggplot2::theme(legend.position = "right")
}

changepoint_distribution_plot <- function(results, n, selected_jump = 0.6,
                                          bins = 100L) {
  selected <- results[
    abs(results$jump - selected_jump) < 1e-10 & results$pattern != "none",
  ]
  rows <- lapply(seq_len(nrow(selected)), function(i) {
    cps <- selected$changepoints[[i]]
    cps <- cps[cps < n]
    if (!length(cps)) return(NULL)
    data.frame(pattern = selected$pattern[i],
               algorithm = selected$algorithm[i], changepoint = cps)
  })
  data <- dplyr::bind_rows(rows)
  data$pattern <- factor(data$pattern, c("rand1", "up", "updown"))
  ggplot2::ggplot(data, ggplot2::aes(changepoint)) +
    ggplot2::geom_histogram(bins = bins, fill = "steelblue",
                            alpha = 0.75, colour = "black", linewidth = 0.15) +
    ggplot2::facet_grid(pattern ~ algorithm) +
    ggplot2::labs(x = "Sequence Position (Time)",
                  y = "Frequency of Detected Change") +
    paper_theme() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

save_power_outputs <- function(results, root, scenario_plot,
                               selected_jump = 0.6) {
  dir.create(file.path(root, "plots"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, file.path(root, "results.rds"))
  utils::write.csv(dplyr::select(results, -changepoints),
                   file.path(root, "results.csv"), row.names = FALSE)
  mse_plot <- metric_plot(results, "MSE", "Signal MSE")
  mse_summary <- summarise_metric(results, "MSE")
  positive <- mse_summary$mean[is.finite(mse_summary$mean) & mse_summary$mean > 0]
  if (length(positive) && max(positive) > 20 * stats::median(positive)) {
    mse_plot <- mse_plot +
      ggplot2::scale_y_continuous(
        trans = scales::pseudo_log_trans(sigma = 0.01, base = 10),
        breaks = c(0, 0.1, 1, 10, 100, 1000, 10000),
        labels = c("0", "0.1", "1", "10", "100", "1,000", "10,000")
      ) +
      ggplot2::labs(y = "Signal MSE (pseudo-log scale)")
  }
  plots <- list(
    "01_scenarios.pdf" = scenario_plot,
    "02_f1.pdf" = metric_plot(results, "F1", "F1", TRUE),
    "03_precision_recall.pdf" = precision_recall_plot(results),
    "04_changepoint_distributions.pdf" =
      changepoint_distribution_plot(results, max(unlist(results$changepoints)),
                                    selected_jump),
    "05_correct_number_probability.pdf" =
      metric_plot(results, "CorrectNumCP",
                  "Probability of the correct number of changepoints", TRUE),
    "06_signal_mse.pdf" = mse_plot
  )
  sizes <- list(c(10, 5), c(10, 5), c(10, 10), c(15, 10), c(10, 5), c(10, 5))
  for (i in seq_along(plots)) {
    ggplot2::ggsave(file.path(root, "plots", names(plots)[i]), plots[[i]],
                    width = sizes[[i]][1], height = sizes[[i]][2])
  }
  invisible(plots)
}
