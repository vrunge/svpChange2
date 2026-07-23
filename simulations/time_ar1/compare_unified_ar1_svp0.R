## Compare exact AR1 SVP, AR1Focus SVP, and svp0 with the exact AR1 validity.

library(svpChange2)

## Exact AR1 validity adapter. This is the fair svp0 comparison: svp0 still
## uses the same ordinary SSE segment cost as SVP, while the validity decision
## is exactly the package AR1 likelihood test.
exact_ar1_validity <- function(rho, sigma2 = 1) {
  force(rho); force(sigma2)
  function(segment, gamma) {
    AR1_single_change(segment, gamma = gamma, rho = rho,
                      sigma2 = sigma2, profile_sigma = FALSE)$valid
  }
}

simulate_ar1 <- function(n = 1000, tau = c(300, 650), means = c(0, 2, -1),
                         rho = .7, sigma = 3) {
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
  exact_test <- exact_ar1_validity(rho, sigma2)
  te <- system.time(exact_svp0 <- svp0(
    y, gamma, exact_test, prune_after_if_unvalid=TRUE,
    prune_if_PELT=FALSE
  ))
  list(partitions = list(AR1=exact$changepoints,
                         AR1Focus=focus$changepoints,
                         exact_svp0=exact_svp0$changepoints),
       exact_svp0 = exact_svp0$changepoints,
       elapsed = c(AR1=tm[["elapsed"]], AR1Focus=tf[["elapsed"]],
                   exact_svp0=te[["elapsed"]]))
}

## Example:
z <- compare_one(simulate_ar1()); z$partitions; z$elapsed
