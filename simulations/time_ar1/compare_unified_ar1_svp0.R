## Compare exact AR1 SVP, AR1Focus SVP, and svp0 with unified_focus ARP
## validity. Install with remotes::install_github("gtromano/unified_focus",
## subdir = "focus") before running.

library(svpChange2)
if (!requireNamespace("focus", quietly = TRUE))
  stop("Install unified_focus with remotes::install_github('gtromano/unified_focus', subdir='focus')")

unified_ar1_validity <- function(rho) {
  force(rho)
  function(segment, gamma) {
    ans <- focus::focus_offline(
      segment, threshold = gamma, type = "arp", family = "arp",
      rho = rho
    )
    stat <- as.numeric(ans$stat)
    stat <- stat[is.finite(stat)]
    !length(stat) || max(stat) < gamma
  }
}

simulate_ar1 <- function(n = 1000, tau = c(300, 650), means = c(0, 2, -1),
                         rho = .7, sigma = 1) {
  mu <- rep(means, diff(c(0, tau, n))); y <- numeric(n)
  y[1] <- mu[1] + rnorm(1, sd = sigma / sqrt(1-rho^2))
  for (i in 2:n)
    y[i] <- mu[i] + rho * (y[i-1] - mu[i-1]) + rnorm(1, sd=sigma)
  y
}

compare_one <- function(y, rho = .7, sigma2 = 1, gamma = 8) {
  tm <- system.time(exact <- SVP(y, gamma, "AR1", rho=rho, sigma2=sigma2,
                                  prune_after_if_unvalid=TRUE,
                                  prune_before_if_invalid=FALSE))
  tf <- system.time(focus <- SVP(y, gamma, "AR1Focus", rho=rho,
                                 sigma2=sigma2,
                                 prune_after_if_unvalid=TRUE,
                                 prune_before_if_invalid=FALSE))
  test <- unified_ar1_validity(rho)
  tv <- system.time(unified <- svp0(y, gamma, test,
                                    prune_after_if_unvalid=TRUE,
                                    prune_if_PELT=FALSE))
  list(partitions = list(AR1=exact$changepoints,
                         AR1Focus=focus$changepoints,
                         unified_svp0=unified$changepoints),
       elapsed = c(AR1=tm[["elapsed"]], AR1Focus=tf[["elapsed"]],
                   unified_svp0=tv[["elapsed"]]))
}

## Example:
## set.seed(1); z <- compare_one(simulate_ar1()); z$partitions; z$elapsed
