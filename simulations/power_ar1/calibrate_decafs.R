## Null calibration for the pure-AR(1) DeCAFS penalty.
## Run from the svpChange2 package root.

saved_run_flag <- Sys.getenv("SVP_RUN_SIMULATIONS", unset = NA_character_)
Sys.setenv(SVP_RUN_SIMULATIONS = "false")
source(file.path("simulations", "power_ar1", "run_power.R"))
if (is.na(saved_run_flag)) {
  Sys.unsetenv("SVP_RUN_SIMULATIONS")
} else {
  Sys.setenv(SVP_RUN_SIMULATIONS = saved_run_flag)
}
rm(saved_run_flag)

calibrate_decafs_penalty <- function(
    n = 600L, rho = 0.8, grid = seq(2, 10, 0.25), reps = 1000L,
    target_no_change = 0.95, seed = 73129L) {
  set_power_seed(seed)
  null_series <- replicate(
    reps, simulate_ar1_noise(n, rho), simplify = FALSE
  )
  scores <- lapply(grid, function(constant) {
    segments <- integer(reps)
    for (rep in seq_len(reps)) {
      segments[rep] <- length(decafs_boundaries(
        null_series[[rep]], rho, 1 - rho^2, constant * log(n)
      ))
    }
    message("DeCAFS calibration: completed multiplier ", constant)
    data.frame(
      constant = constant,
      no_change_rate = mean(segments == 1L),
      mean_segments = mean(segments),
      false_positive_rate = mean(segments > 1L)
    )
  })
  table <- dplyr::bind_rows(scores)
  admissible <- dplyr::filter(table, no_change_rate >= target_no_change)
  if (!nrow(admissible)) {
    stop("DeCAFS calibration grid does not attain the null target")
  }
  list(constant = min(admissible$constant), table = table)
}

evaluate_decafs_power_tradeoff <- function(
    n = 600L, rho = 0.8,
    grid = c(3.5, 4, 4.25, 4.5, 4.75, 5, 5.5, 6, 6.25),
    jump_sizes = seq(1, 3, 0.1), reps = 100L, seed = 123L) {
  design <- expand.grid(
    pattern = POWER_PATTERNS,
    jump = jump_sizes,
    rep = seq_len(reps),
    stringsAsFactors = FALSE
  )
  series <- lapply(seq_len(nrow(design)), function(i) {
    set_power_seed(seed + i)
    mu <- generate_power_signal(n, design$pattern[i], design$jump[i])
    list(
      y = mu + simulate_ar1_noise(n, rho),
      truth = normalise_boundaries(which(diff(mu) != 0), n)
    )
  })
  scores <- lapply(grid, function(constant) {
    rows <- lapply(seq_along(series), function(i) {
      boundaries <- decafs_boundaries(
        series[[i]]$y, rho, 1 - rho^2, constant * log(n)
      )
      metric <- paper_metrics(
        series[[i]]$truth, boundaries, max(5L, round(n * 0.0025))
      )
      data.frame(
        pattern = design$pattern[i],
        F1 = metric[["F1"]],
        CorrectNumCP = metric[["CorrectNumCP"]],
        NumSegments = length(boundaries)
      )
    })
    message("DeCAFS power tradeoff: completed multiplier ", constant)
    dplyr::bind_rows(rows) |>
      dplyr::group_by(pattern) |>
      dplyr::summarise(
        mean_F1 = mean(F1),
        correct_number_rate = mean(CorrectNumCP),
        mean_segments = mean(NumSegments),
        .groups = "drop"
      ) |>
      dplyr::mutate(constant = constant, .before = 1L)
  })
  dplyr::bind_rows(scores)
}

calibration <- calibrate_decafs_penalty(target_no_change = 0.95)
utils::write.csv(
  calibration$table,
  file.path(AR1_ROOT, "decafs_null_calibration.csv"),
  row.names = FALSE
)
tradeoff <- evaluate_decafs_power_tradeoff()
utils::write.csv(
  tradeoff,
  file.path(AR1_ROOT, "decafs_power_tradeoff.csv"),
  row.names = FALSE
)
message("Selected DeCAFS multiplier: ", calibration$constant)
