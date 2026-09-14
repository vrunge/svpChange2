## Calibration utilities for the v0.2.0 figure refresh.
##
## The calibration target is the null false-positive probability of PELT.
## Signal F1 is used only as a deterministic tie-breaker.  All data are
## generated before candidate evaluation, so every candidate sees identical
## observations.

if (!exists("POWER_PATTERNS", inherits = TRUE)) {
  source(file.path("simulations", "power_common.R"))
}

calibration_seed <- function(seed) {
  set.seed(seed, kind = "Mersenne-Twister",
           normal.kind = "Inversion", sample.kind = "Rejection")
}

calibration_pelt <- function(y) {
  fit <- changepoint::cpt.mean(
    y, method = "PELT", penalty = "Manual",
    pen.value = 2 * log(length(y))
  )
  length(unique(c(changepoint::cpts(fit), length(y)))) > 1L
}

calibration_candidate_fit <- function(y, family, mode, value,
                                      true_segments = 1L,
                                      rho = 0.8) {
  n <- length(y)
  if (family == "gaussian") {
    gamma <- value * log(n)
    return(SVP(y, gamma, "gaussian_mean", subtests = mode)$changepoints)
  }
  if (family == "ar1") {
    gamma <- value * log(n)
    return(SVP(
      y, gamma, "AR1Focus", subtests = mode,
      rho = rho, sigma2 = 1 - rho^2
    )$changepoints)
  }
  if (family == "wilcoxon") {
    gamma <- value * sqrt((n / true_segments)^3 / 12)
    return(SVP(y, gamma, "WilcoxonCost", subtests = mode)$changepoints)
  }
  if (family == "mood") {
    splits <- n / true_segments - 1
    single_alpha <- 1 - (1 - value)^(1 / splits)
    gamma <- stats::qchisq(1 - single_alpha, df = true_segments)
    return(SVP(y, gamma, "MedianMoodCost", subtests = mode)$changepoints)
  }
  stop("Unknown calibration family: ", family)
}

make_calibration_data <- function(n, reps, noise, seed,
                                  signal = FALSE, jump_values = NULL) {
  if (!signal) {
    return(lapply(seq_len(reps), function(i) {
      calibration_seed(seed + i)
      list(y = noise(n), pattern = "none", jump = 0,
           truth = n, true_segments = 1L)
    }))
  }

  if (is.null(jump_values)) {
    jump_values <- c(0.2, 0.6, 1.0, 1.4, 1.8)
  }
  design <- expand.grid(
    pattern = POWER_PATTERNS[POWER_PATTERNS != "none"],
    jump = jump_values, rep = seq_len(max(1L, ceiling(reps / 15L))),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  design <- design[seq_len(min(reps, nrow(design))), , drop = FALSE]
  lapply(seq_len(nrow(design)), function(i) {
    task <- design[i, ]
    calibration_seed(seed + i)
    mu <- generate_power_signal(n, task$pattern, task$jump)
    list(
      y = mu + noise(n), pattern = task$pattern, jump = task$jump,
      truth = normalise_boundaries(which(diff(mu) != 0), n),
      true_segments = sum(diff(mu) != 0) + 1L
    )
  })
}

select_calibration <- function(scores) {
  ordered <- scores[order(
    scores$absolute_fpr_gap,
    -scores$signal_f1,
    scores$elapsed_median,
    scores$value,
    scores$mode
  ), , drop = FALSE]
  ordered[1L, , drop = FALSE]
}

select_calibration_by_mode <- function(scores) {
  modes <- unique(scores$mode)
  dplyr::bind_rows(lapply(modes, function(mode) {
    select_calibration(scores[scores$mode == mode, , drop = FALSE])
  }))
}

calibrate_svp_family <- function(
    family, candidates, n, noise, seed,
    null_reps = 1000L, signal_reps = 200L,
    rho = 0.8, workers = power_default_workers()) {
  null_data <- make_calibration_data(n, null_reps, noise, seed)
  pelt_fpr <- mean(vapply(null_data, function(x) calibration_pelt(x$y), logical(1)))
  signal_data <- make_calibration_data(
    n, signal_reps, noise, seed + 100000L, signal = TRUE
  )

  rows <- power_lapply(seq_len(nrow(candidates)), function(i) {
    candidate <- candidates[i, , drop = FALSE]
    null_fit <- lapply(null_data, function(x) {
      t0 <- proc.time()[["elapsed"]]
      boundaries <- calibration_candidate_fit(
        x$y, family, candidate$mode, candidate$value,
        true_segments = x$true_segments, rho = rho
      )
      data.frame(
        false_positive = as.numeric(length(boundaries) > 1L),
        elapsed = proc.time()[["elapsed"]] - t0
      )
    })
    signal_f1 <- vapply(signal_data, function(x) {
      boundaries <- calibration_candidate_fit(
        x$y, family, candidate$mode, candidate$value,
        true_segments = x$true_segments, rho = rho
      )
      paper_metrics(x$truth, boundaries, round(0.0025 * n))[["F1"]]
    }, numeric(1))
    null_frame <- dplyr::bind_rows(null_fit)
    data.frame(
      family = family, mode = candidate$mode, value = candidate$value,
      pelt_fpr = pelt_fpr,
      svp_fpr = mean(null_frame$false_positive),
      absolute_fpr_gap = abs(mean(null_frame$false_positive) - pelt_fpr),
      signal_f1 = mean(signal_f1, na.rm = TRUE),
      elapsed_median = stats::median(null_frame$elapsed),
      null_reps = null_reps, signal_reps = length(signal_data),
      seed = seed
    )
  }, workers = workers)
  scores <- dplyr::bind_rows(rows)
  list(
    scores = scores,
    selected = select_calibration(scores),
    ## Keep one calibrated value for each admissible mode: the refreshed
    ## paper figures explicitly compare right and both.
    selected_by_mode = select_calibration_by_mode(scores)
  )
}

default_calibration_candidates <- function() {
  list(
    gaussian = expand.grid(
      mode = c("right", "both"), value = seq(1, 5, 0.1),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    ),
    ar1 = expand.grid(
      mode = c("right", "both"), value = seq(1, 5, 0.1),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    ),
    wilcoxon = expand.grid(
      mode = c("right", "both"), value = seq(0.5, 3, 0.1),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    ),
    mood = expand.grid(
      mode = c("right", "both"),
      value = c(0.001, 0.0025, 0.005, 0.01, 0.02, 0.05),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
  )
}
