## Reproduce the historical Experiment 1 timing study with the current SVP.
## The data-generating process and plotting convention match
## time_gaussian_legacy/2_complexity_study.R exactly.

library(svpChange2)
library(changepoint)
library(ggplot2)

legacy_time_pelt <- function(data, penalty) {
  elapsed <- system.time({
    fit <- tryCatch(
      changepoint::cpt.mean(
        data, method = "PELT", penalty = "Manual", pen.value = penalty
      ),
      error = function(e) NULL
    )
  })[["elapsed"]]
  detected <- if (is.null(fit)) NA_integer_ else {
    length(as.integer(unlist(changepoint::cpts(fit)))) + 1L
  }
  list(time = if (is.finite(elapsed)) elapsed else NA_real_,
       nchanges = detected)
}

legacy_time_svp <- function(data, penalty, subtests = "both") {
  elapsed <- system.time({
    fit <- tryCatch(
      SVP(data, penalty, "gaussian_mean", subtests = subtests),
      error = function(e) NULL
    )
  })[["elapsed"]]
  detected <- if (is.null(fit)) NA_integer_ else {
    length(as.integer(unlist(fit$changepoints)))
  }
  list(time = if (is.finite(elapsed)) elapsed else NA_real_,
       nchanges = detected)
}

run_current_exp1 <- function(
    n_values = round(2^(seq(10, 14, length.out = 20))),
    reps = 20L,
    penalty_large = 1000,
    seed = 123L) {
  set.seed(seed)
  rows <- vector("list", length(n_values) * reps)
  row_id <- 0L
  for (n in n_values) {
    message("[legacy exp1] n = ", n)
    for (replicate_id in seq_len(reps)) {
      data <- rnorm(n, mean = 0, sd = 1)
      pelt <- legacy_time_pelt(data, penalty_large)
      svp <- legacy_time_svp(data, penalty_large, "both")
      svp_right <- legacy_time_svp(data, penalty_large, "right")
      row_id <- row_id + 1L
      rows[[row_id]] <- data.frame(
        n = as.integer(n), replicate = replicate_id,
        PELT = pelt$time, PELT_nchanges = pelt$nchanges,
        SVP = svp$time, SVP_nchanges = svp$nchanges,
        SVP_right = svp_right$time,
        SVP_right_nchanges = svp_right$nchanges
      )
    }
  }
  do.call(rbind, rows)
}

results <- run_current_exp1()
legacy_root <- file.path("simulations", "other_simus", "time_gaussian_legacy")
write.csv(
  results,
  file.path(legacy_root, "complexity_exp1_nochange_current.csv"),
  row.names = FALSE
)

valid <- is.finite(results$PELT) & is.finite(results$SVP) &
  results$PELT > 0 & results$SVP > 0
if (sum(valid) > 2L) {
  slope_pelt <- unname(coef(lm(log(PELT) ~ log(n), data = results[valid, ]))[2])
  slope_svp <- unname(coef(lm(log(SVP) ~ log(n), data = results[valid, ]))[2])
} else {
  slope_pelt <- NA_real_
  slope_svp <- NA_real_
}

plot_data <- results[valid, , drop = FALSE]
p <- ggplot(plot_data, aes(x = n, y = PELT, colour = "PELT")) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.2) +
  geom_smooth(
    data = plot_data, aes(x = n, y = SVP, colour = "SVP"),
    method = "lm", formula = y ~ x, se = TRUE, alpha = 0.2
  ) +
  geom_smooth(
    data = plot_data[is.finite(plot_data$SVP_right) & plot_data$SVP_right > 0, ],
    aes(x = n, y = SVP_right, colour = "SVP (right)"),
    method = "lm", formula = y ~ x, se = TRUE, alpha = 0.2
  ) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "n (log)", y = "Time (s, log)") +
  theme_minimal()
if (is.finite(slope_pelt)) {
  p <- p + annotate(
    "text", x = Inf, y = Inf,
    label = paste("PELT slope:", round(slope_pelt, 3)),
    hjust = 1.1, vjust = 20, colour = "#F8766D", size = 4
  )
}
if (is.finite(slope_svp)) {
  p <- p + annotate(
    "text", x = Inf, y = Inf,
    label = paste("SVP slope:", round(slope_svp, 3)),
    hjust = 1.1, vjust = 22, colour = "#00BFC4", size = 4
  )
}

figure_root <- file.path("..", "SVP_NEW_Figures2")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
ggsave(
  file.path(figure_root, "2_exp1_time_vs_n.pdf"), p,
  width = 8, height = 5
)

summary <- data.frame(
  n = sort(unique(results$n)),
  PELT_mean = as.numeric(tapply(results$PELT, results$n, mean, na.rm = TRUE)),
  SVP_both_mean = as.numeric(tapply(results$SVP, results$n, mean, na.rm = TRUE)),
  SVP_right_mean = as.numeric(
    tapply(results$SVP_right, results$n, mean, na.rm = TRUE)
  )
)
write.csv(
  summary,
  file.path(figure_root, "2_exp1_time_vs_n_current_summary.csv"),
  row.names = FALSE
)
message(sprintf("Current slopes: PELT = %.3f, SVP = %.3f", slope_pelt, slope_svp))
