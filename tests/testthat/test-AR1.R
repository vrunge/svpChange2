simulate_ar1_change <- function(n, changepoint, means, rho, innovation_sd = 1) {
  innovations <- rnorm(n, sd = innovation_sd)
  mu <- rep(means, c(changepoint, n - changepoint))
  data <- numeric(n)
  data[1] <- mu[1] + innovations[1] / sqrt(1 - rho^2)
  for (i in 2:n) {
    data[i] <- mu[i] + rho * (data[i - 1] - mu[i - 1]) + innovations[i]
  }
  data
}

exact_ar1_reference <- function(data, rho, sigma2 = 1, profile_sigma = FALSE) {
  n <- length(data)
  z <- data[-1] - rho * data[-n]
  rss0 <- sum((z - mean(z))^2)
  taus <- 2:(n - 2)
  rss1 <- vapply(taus, function(tau) {
    design <- matrix(0, nrow = n - 1, ncol = 2)
    design[seq_len(tau - 1), 1] <- 1 - rho
    design[tau, ] <- c(-rho, 1)
    if (tau + 1 <= n - 1) {
      design[(tau + 1):(n - 1), 2] <- 1 - rho
    }
    sum(lm.fit(design, z)$residuals^2)
  }, numeric(1))
  best <- which.min(rss1)
  statistic <- if (profile_sigma) {
    (n - 1) * log(rss0 / rss1[best])
  } else {
    (rss0 - rss1[best]) / sigma2
  }
  list(
    statistic = statistic,
    changepoint = taus[best],
    rss0 = rss0,
    rss1 = rss1[best]
  )
}

test_that("robust AR1 estimate matches equation 3", {
  set.seed(1)
  data <- simulate_ar1_change(301, 150, c(0, 1.5), 0.6)
  lag1_median <- median(abs(diff(data)))
  lag2_median <- median(abs(data[3:length(data)] - data[1:(length(data) - 2)]))
  reference <- lag2_median^2 / lag1_median^2 - 1

  expect_equal(AR1_rho(data), max(-0.999, min(0.999, reference)), tolerance = 1e-12)
})

test_that("AR1 single-change test locates a simulated mean change", {
  set.seed(42)
  data <- simulate_ar1_change(400, 200, c(0, 2), 0.7)
  result <- AR1_single_change(data, gamma = 10, rho = 0.7)

  expect_false(result$valid)
  expect_lte(abs(result$changepoint - 200), 15)
})

test_that("AR1 statistic equals the exact fixed-rho reference", {
  set.seed(7)
  data <- simulate_ar1_change(120, 53, c(-0.5, 1.25), 0.65)
  reference <- exact_ar1_reference(data, rho = 0.65, sigma2 = 1)
  result <- AR1_single_change(data, gamma = Inf, rho = 0.65, sigma2 = 1)

  expect_equal(result$statistic, reference$statistic, tolerance = 1e-9)
  expect_equal(result$rss0, reference$rss0, tolerance = 1e-9)
  expect_equal(result$rss1, reference$rss1, tolerance = 1e-9)
  expect_equal(result$changepoint, reference$changepoint)
})

test_that("AR1 profiled statistic equals the exact reference", {
  set.seed(8)
  data <- simulate_ar1_change(100, 47, c(0, 1), 0.4, innovation_sd = 1.3)
  reference <- exact_ar1_reference(data, rho = 0.4, profile_sigma = TRUE)
  result <- AR1_single_change(
    data, gamma = Inf, rho = 0.4, profile_sigma = TRUE
  )

  expect_equal(result$statistic, reference$statistic, tolerance = 1e-9)
  expect_equal(result$changepoint, reference$changepoint)
})

test_that("AR1 validity test works inside SVP", {
  set.seed(42)
  data <- simulate_ar1_change(400, 200, c(0, 2), 0.7)
  result <- SVP(data, gamma = 10, test = "AR1", rho = 0.7)

  expect_equal(result$rho, 0.7)
  expect_equal(tail(result$changepoints, 1), length(data))
  expect_true(any(abs(head(result$changepoints, -1) - 200) <= 10))
})

test_that("AR1 validity test works inside main SVP with pruning options", {
  set.seed(43)
  data <- simulate_ar1_change(400, 200, c(0, 2), 0.7)
  result <- SVP(
    data,
    gamma = 10,
    test = "AR1",
    prune_after_if_unvalid = TRUE,
    prune_before_if_invalid = TRUE,
    rho = 0.7,
    sigma2 = 1
  )

  expect_equal(tail(result$changepoints, 1), length(data))
  expect_true(any(abs(head(result$changepoints, -1) - 200) <= 10))
})
