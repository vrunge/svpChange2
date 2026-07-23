## Compare AR1, AR1Profile, and AR1Focus partitions on identical simulations.
## No files or figures are produced.

library(svpChange2)

simulate_ar1 <- function(n = 1000, tau = c(300, 650), means = c(0, 2, -1),
                         rho = 0.7, sigma = 1) {
  stopifnot(length(tau) + 1L == length(means), all(diff(tau) > 0))
  stopifnot(all(tau > 0), all(tau < n))
  mu <- rep(means, diff(c(0, tau, n)))
  y <- numeric(n)
  y[1] <- mu[1] + rnorm(1, sd = sigma / sqrt(1 - rho^2))
  for (i in 2:n)
    y[i] <- mu[i] + rho * (y[i - 1] - mu[i - 1]) + rnorm(1, sd = sigma)
  y
}

fit_three_methods <- function(y, rho = 0.7, sigma2 = 1, gamma = 8) {
  methods <- c("AR1Profile", "AR1", "AR1Focus")
  fits <- setNames(lapply(methods, function(method) {
    SVP(y, gamma = gamma, test = method, rho = rho, sigma2 = sigma2,
        profile_sigma = method == "AR1Profile",
        prune_after_if_unvalid = TRUE,
        prune_before_if_invalid = FALSE)
  }), methods)
  lapply(fits, function(fit) as.integer(fit$changepoints))
}

partition_diagnostics <- function(partition, true_changepoints) {
  estimated <- partition[partition < max(partition)]
  list(
    changepoints = estimated,
    number = length(estimated),
    correct_number = length(estimated) == length(true_changepoints),
    localization_error = if (length(estimated) == length(true_changepoints))
      max(abs(sort(estimated) - sort(true_changepoints))) else NA_real_
  )
}

compare_replicates <- function(B = 10, n = 1000, rho = 0.7,
                               sigma = 1, gamma = 8) {
  out <- vector("list", B)
  for (b in seq_len(B)) {
    y <- simulate_ar1(n = n, rho = rho, sigma = sigma)
    p <- fit_three_methods(y, rho = rho, sigma2 = sigma^2, gamma = gamma)
    out[[b]] <- list(replicate = b, partitions = p)
    cat("replicate", b, "\n")
    for (method in names(p))
      cat(" ", method, ": ", paste(p[[method]], collapse = ", "), "\n", sep = "")
    for (method in names(p)) {
      d <- partition_diagnostics(p[[method]], c(300L, 650L))
      cat("  ", method, " localization=", d$localization_error,
          " correct_number=", d$correct_number, "\n", sep = "")
    }
  }
  invisible(out)
}

## Run with an unfixed seed so every call simulates new data:
compare_replicates(B = 10)
