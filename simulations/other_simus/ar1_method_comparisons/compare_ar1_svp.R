## Compare exact and Gaetano/FOCuS AR(1) SVP tests.
## No files are written; the optional panel is drawn in the active device.

if (!requireNamespace("svpChange2", quietly = TRUE))
  stop("Install svpChange2 first")
library(svpChange2)

simulate_ar1_mean_change <- function(n = 2000, tau = n %/% 2,
                                     mu1 = 0, mu2 = 2, rho = .7,
                                     sigma = 1) {
  mu <- c(rep(mu1, tau), rep(mu2, n - tau))
  y <- numeric(n); y[1] <- mu[1] + rnorm(1, sd = sigma / sqrt(1-rho^2))
  for (i in 2:n) y[i] <- mu[i] + rho * (y[i-1] - mu[i-1]) + rnorm(1, sd=sigma)
  y
}

run_ar1_svp_comparison <- function(y, gamma = 8, rho = .7, sigma2 = 1) {
  methods <- c("AR1", "AR1Profile", "AR1Focus")
  fits <- lapply(methods, function(m) {
    system.time(SVP(y, gamma = gamma, test = m, rho = rho,
                    sigma2 = sigma2, profile_sigma = (m == "AR1Profile"),
                    prune_after_if_unvalid = TRUE,
                    prune_before_if_invalid = FALSE))
  })
  ## Re-run to retain the actual results and measure elapsed time separately.
  results <- lapply(methods, function(m) SVP(
    y, gamma = gamma, test = m, rho = rho, sigma2 = sigma2,
    profile_sigma = (m == "AR1Profile"), prune_after_if_unvalid = TRUE,
    prune_before_if_invalid = FALSE))
  times <- sapply(seq_along(methods), function(i) unname(fits[[i]][["elapsed"]]))
  list(partitions = setNames(lapply(results, `[[`, "changepoints"), methods),
       elapsed = setNames(times, methods), results = setNames(results, methods))
}

benchmark_ar1_svp <- function(n_values = c(200, 500, 1000, 2000),
                              repetitions = 3, rho = .7, sigma = 1,
                              gamma = 8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  rows <- list(); z <- 1L
  for (n in n_values) for (b in seq_len(repetitions)) {
    y <- simulate_ar1_mean_change(n=n, tau=n %/% 2, rho=rho, sigma=sigma)
    for (m in c("AR1", "AR1Profile", "AR1Focus")) {
      tm <- system.time(SVP(y, gamma=gamma, test=m, rho=rho,
                             sigma2=sigma^2, profile_sigma=(m == "AR1Profile"),
                             prune_after_if_unvalid=TRUE,
                             prune_before_if_invalid=FALSE))
      rows[[z]] <- data.frame(n=n, replicate=b, method=m,
                              elapsed=unname(tm[["elapsed"]])); z <- z+1L
    }
  }
  do.call(rbind, rows)
}

## Example interactive run:
## y <- simulate_ar1_mean_change(); cmp <- run_ar1_svp_comparison(y)
## cmp$partitions; cmp$elapsed
## benchmark_ar1_svp(seed = NULL)
