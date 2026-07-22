## Self-contained comparison of the exact svpChange2 AR1 statistic and the
## unified_focus ARP/FOCuS statistic. No figures or output files are created.

if (!requireNamespace("svpChange2", quietly = TRUE))
  stop("Install svpChange2 first")
if (!requireNamespace("focus", quietly = TRUE))
  stop("Install with remotes::install_github('gtromano/unified_focus', subdir='focus')")

simulate_ar1 <- function(n = 1000, tau = NULL, mu1 = 0, mu2 = 0,
                         rho = 0.7, sigma = 1) {
  mu <- if (is.null(tau)) rep(mu1, n) else c(rep(mu1, tau), rep(mu2, n-tau))
  y <- numeric(n)
  y[1] <- mu[1] + stats::rnorm(1, sd = sigma / sqrt(1-rho^2))
  for (i in 2:n)
    y[i] <- mu[i] + rho * (y[i-1] - mu[i-1]) + stats::rnorm(1, sd=sigma)
  y
}

compare_statistics <- function(y, gamma = 8, rho = 0.7, sigma2 = 1,
                               mu0_arp = NULL) {
  exact <- svpChange2::AR1_single_change(
    y, gamma = gamma, rho = rho, sigma2 = sigma2,
    profile_sigma = FALSE
  )
  focus_fit <- if (is.null(mu0_arp)) {
    focus::focus_offline(y, threshold=gamma, type="arp", family="arp", rho=rho)
  } else {
    focus::focus_offline(y, threshold=gamma, type="arp", family="arp",
                         rho=rho, mu0_arp=mu0_arp)
  }
  focus_stat <- as.numeric(focus_fit$stat)
  keep <- is.finite(focus_stat)
  focus_stat <- focus_stat[keep]
  focus_time <- which(keep)[which.max(focus_stat)]
  focus_cp <- as.integer(focus_fit$changepoint[focus_time])
  data.frame(
    exact_statistic = exact$statistic,
    exact_changepoint = exact$changepoint,
    focus_statistic = max(focus_stat, 0),
    focus_changepoint = focus_cp,
    statistic_ratio = max(focus_stat, 0) / exact$statistic,
    same_changepoint = identical(as.integer(exact$changepoint), focus_cp)
  )
}

run_comparison <- function(seed = 10, n = 1000, rho = 0.7, sigma = 1,
                           gamma = 8) {
  set.seed(seed)
  scenarios <- list(
    no_change = simulate_ar1(n=n, rho=rho, sigma=sigma),
    weak_change = simulate_ar1(n=n, tau=n/2, mu2=1, rho=rho, sigma=sigma),
    strong_change = simulate_ar1(n=n, tau=n/2, mu2=3, rho=rho, sigma=sigma)
  )
  result <- do.call(rbind, lapply(names(scenarios), function(s) {
    cbind(scenario=s, compare_statistics(scenarios[[s]], gamma, rho, sigma^2))
  }))
  rownames(result) <- NULL
  result
}

## Run explicitly:
## print(run_comparison())
