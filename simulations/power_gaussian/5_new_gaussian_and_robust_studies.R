## Gaussian and heavy-tail power studies with the two new svpChange2 methods.

library(svpChange2)
options(svpChange2.simulation_root = file.path("simulations", "power_gaussian"))
library(changepoint)
library(ggplot2)
library(dplyr)
library(tidyr)

source(file.path("simulations", "0_simulation_helpers.R"))

run_method <- function(data,
                       method,
                       noise = c("gaussian", "student"),
                       nb_segments = 8,
                       q_alpha_n = 3) {
  noise <- match.arg(noise)
  n <- length(data)

  if (method == "PELT") {
    fit <- changepoint::cpt.mean(
      data,
      method = "PELT",
      penalty = "Manual",
      pen.value = 2 * log(n)
    )
    return(sort(unique(c(as.integer(changepoint::cpts(fit)), n))))
  }

  if (method %in% names(svp_method_specs())) {
    ## For heavy tails, a robust scale prevents isolated extremes from making
    ## the multiscale threshold depend on the non-existent variance of t_2.
    sigma2 <- if (noise == "student") stats::mad(data)^2 else 1
    fit <- run_svp_method(data, method, q_alpha_n, sigma2)
    return(extract_boundaries(fit, n))
  }

  if (method == "SVP (BIC)") {
    fit <- SVP(data, gamma = 2 * log(n), test = "gaussian_mean")
    return(extract_boundaries(fit, n))
  }

  if (method == "SVP (Wilcoxon)") {
    fit <- SVP(
      data,
      gamma = 1.5 * sqrt((n / nb_segments)^3 / 12),
      test = "WilcoxonCost",
      prune_after_if_unvalid = TRUE
    )
    return(extract_boundaries(fit, n))
  }

  if (method == "SVP (MedianMood)") {
    intervals <- max(n / nb_segments - 1, 1)
    single_alpha <- 1 - (1 - 0.01)^(1 / intervals)
    fit <- SVP(
      data,
      gamma = stats::qchisq(1 - single_alpha, df = nb_segments),
      test = "MedianMoodCost",
      prune_after_if_unvalid = TRUE
    )
    return(extract_boundaries(fit, n))
  }

  if (method == "rfpop") {
    if (!requireNamespace("robseg", quietly = TRUE)) {
      stop("the optional package 'robseg' is required for rfpop")
    }
    scale <- stats::mad(data)
    fit <- robseg::Rob_seg.std(data / scale, loss = "Outlier",
                              lambda = 2 * log(n), lthreshold = 3)
    return(sort(unique(c(as.integer(head(fit$t.est, -1L)), n))))
  }

  stop("Unknown method: ", method)
}

run_power_study_v2 <- function(n = 1000,
                               patterns = c("none", "up", "updown", "rand1"),
                               jump_sizes = seq(0.1, 2, 0.1),
                               reps = 100,
                               noise = c("gaussian", "student"),
                               nb_segments = 8,
                               q_alpha_n = 3,
                               seed = 123) {
  noise <- match.arg(noise)
  methods <- if (noise == "gaussian") {
    c("PELT", names(svp_method_specs(q_alpha_n)), "SVP (BIC)")
  } else {
    robust_methods <- c("PELT", "SVP (Wilcoxon)", "SVP (MedianMood)",
                        "SVP FOCUS (after/before)", "SVP multiscale")
    if (requireNamespace("robseg", quietly = TRUE)) {
      append(robust_methods, "rfpop", after = 3L)
    } else {
      warning("Optional package 'robseg' is unavailable; omitting rfpop")
      robust_methods
    }
  }
  design <- expand.grid(pattern = patterns,
                        jump = jump_sizes,
                        rep = seq_len(reps),
                        stringsAsFactors = FALSE)

  rows <- vector("list", nrow(design) * length(methods))
  row_index <- 0L
  for (i in seq_len(nrow(design))) {
    set.seed(seed + i)
    mu <- generate_signal(n, design$pattern[i], nb_segments, design$jump[i])
    errors <- if (noise == "gaussian") stats::rnorm(n) else stats::rt(n, df = 2)
    data <- mu + errors
    truth <- c(which(diff(mu) != 0), n)

    for (method in methods) {
      row_index <- row_index + 1L
      boundaries <- tryCatch(
        run_method(data, method, noise, length(truth), q_alpha_n),
        error = function(error) {
          warning(method, " failed: ", conditionMessage(error))
          n
        }
      )
      metrics <- cp_metrics(truth, boundaries, round(n * 0.0025))
      estimate <- piecewise_mean(data, boundaries)
      rows[[row_index]] <- data.frame(
        pattern = design$pattern[i],
        jump = design$jump[i],
        rep = design$rep[i],
        algorithm = method,
        Precision = metrics[["Precision"]],
        Recall = metrics[["Recall"]],
        F1 = metrics[["F1"]],
        NumSegments = length(boundaries),
        MSE = mean((mu - estimate)^2)
      )
      rows[[row_index]]$changepoints <- list(boundaries)
    }
  }
  dplyr::bind_rows(rows)
}

save_power_study_v2 <- function(results, stem) {
  old_root <- getOption("svpChange2.simulation_root", "simulations")
  target_root <- if (grepl("robust", stem, ignore.case = TRUE)) {
    file.path("simulations", "power_robust")
  } else {
    file.path("simulations", "power_gaussian")
  }
  options(svpChange2.simulation_root = target_root)
  on.exit(options(svpChange2.simulation_root = old_root), add = TRUE)
  ensure_simulation_dirs()
  saveRDS(results, simulation_path(paste0(stem, ".rds")))
  utils::write.csv(dplyr::select(results, -changepoints),
                   simulation_path(paste0(stem, ".csv")), row.names = FALSE)
  summarise_and_plot_metrics(results, stem)
  plot_cp_distributions(results, stem, selected_jump = 0.6,
                        n = max(unlist(results$changepoints)))
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  gaussian_results <- run_power_study_v2(noise = "gaussian")
  save_power_study_v2(gaussian_results, "5_gaussian_new_methods")

  robust_results <- run_power_study_v2(
    jump_sizes = seq(0.1, 4, 0.1),
    noise = "student"
  )
  save_power_study_v2(robust_results, "5_robust_new_methods")
}
