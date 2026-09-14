library(testthat)
library(svpChange2)

########## test svp0 result = SVP result with right pruning ##########
########## test svp0 result = SVP result with right pruning ##########
########## test svp0 result = SVP result with right pruning ##########

test_that("svp0 and SVP agree with right pruning", {
  n <- 500
  gap <- 1
  chpts <- c(0.1, 0.3, 0.4, 0.45, 0.55, 0.7, 0.75, 0.95, 1) * n
  data <- ts_generator(
    chpts = chpts,
    parameters = c(0, gap, 0, gap, 0, gap, 0, gap, 0),
    sd_noise = 1
  )
  gamma <- 5

  res_svp0 <- svp0(data,
    gamma,
    test = valid_FOCUS_last,
    subtests = "right",
    PELT_pruning = FALSE
  )

  res_svp <- SVP(
    data = data,
    gamma = gamma,
    test = "gaussian_mean",
    subtests = "right"
  )

  expect_equal(res_svp0$changepoints, res_svp$changepoints)
  expect_equal(res_svp0$R[, 1], res_svp$R[, 1])
})

test_that("svp0 and SVP return the same Gaussian FOCUS result", {
  n <- 100
  gap <- 5
  chpts <- c(0.1, 0.3, 0.4, 0.45, 0.55, 0.7, 0.75, 0.95, 1) * n
  data <- ts_generator(
    chpts = chpts,
    parameters = c(0, gap, 0, gap, 0, gap, 0, gap, 0),
    sd_noise = 1
  )
  gamma <- 7

  res_svp0 <- svp0(data,
    gamma,
    test = valid_FOCUS_last,
    subtests = "right",
    PELT_pruning = FALSE
  )
  res_svp <- SVP(
    data = data,
    gamma = gamma,
    test = "gaussian_mean",
    subtests = "right"
  )
  expect_equal(res_svp$changepoints, res_svp0$changepoints)
  expect_equal(res_svp$R, res_svp0$R)
})

test_that("SVP accepts all before/after pruning combinations", {
  n <- 80
  data <- rep(c(0, 1, -0.5, 0.8), each = n / 4) + rnorm(n)
  gamma <- 2 * log(n)

  opts <- data.frame(subtests = c("none", "right", "left", "both"))

  for (i in seq_len(nrow(opts))) {
    res <- SVP(
      data = data,
      gamma = gamma,
      test = "gaussian_mean",
      subtests = opts$subtests[i]
    )

    expect_equal(tail(res$changepoints, 1), n)
  }
})

test_that("SVP evaluates Gaussian FOCUS at the current endpoint", {
  data <- c(2.002482730292829, 0.066700870930183,
            1.866851844706863)
  gamma <- 0.4873983399942518

  expect_false(valid_FOCUS(data, gamma))
  expect_true(valid_FOCUS_last(data, gamma))

  result <- SVP(
    data,
    gamma = gamma,
    test = "gaussian_mean",
    subtests = "none"
  )
  reference <- svp0(
    data,
    gamma = gamma,
    test = valid_FOCUS_last,
    subtests = "none",
    PELT_pruning = FALSE
  )

  expect_equal(result$changepoints, reference$changepoints)
  expect_equal(result$R, reference$R)
  expect_equal(result$changepoints, length(data))
})

test_that("SVP matches svp0 for every pruning mode", {
  set.seed(27)
  data <- c(rnorm(20), rnorm(20, 2), rnorm(20, -1))
  gamma <- 1.5

  for (mode in c("none", "right", "left", "both")) {
    result <- SVP(
      data,
      gamma = gamma,
      test = "gaussian_mean",
      subtests = mode
    )
    reference <- svp0(
      data,
      gamma = gamma,
      test = valid_FOCUS_last,
      subtests = mode,
      PELT_pruning = FALSE
    )

    expect_equal(result$changepoints, reference$changepoints, info = mode)
    expect_equal(result$R, reference$R, info = mode)
    expect_equal(result$nb, reference$nb, info = mode)
    expect_equal(result$lastIndexSet, reference$lastIndexSet, info = mode)
  }
})

test_that("SVP computes stable squared-error costs for large offsets", {
  data <- 1e12 + c(0, 1, -1, 0, 2, -2)
  result <- SVP(
    data,
    gamma = 1e12,
    test = "varCost",
    subtests = "none"
  )
  expected <- sum((data - mean(data))^2)

  expect_equal(result$changepoints, length(data))
  expect_equal(result$R[length(data), 1], expected, tolerance = 1e-8)
  expect_true(all(result$R[, 1] >= 0))
})

test_that("SVP defaults to Gaussian mean with both pruning rules", {
  set.seed(17)
  data <- rnorm(80)
  gamma <- 1.5 * log(length(data))

  default <- SVP(data, gamma)
  explicit <- SVP(data, gamma, test = "gaussian_mean")
  explicit_both <- SVP(
    data, gamma, test = "gaussian_mean", subtests = "both"
  )

  expect_equal(default$changepoints, explicit$changepoints)
  expect_equal(default$R, explicit$R)
  expect_equal(default$changepoints, explicit_both$changepoints)
  expect_equal(default$R, explicit_both$R)
})

test_that("both Gaussian pruning returns a valid final endpoint", {
  set.seed(2)
  n <- 100
  data <- rnorm(n)
  gamma <- 0.3

  forward <- SVP(
    data,
    gamma,
    test = "gaussian_mean",
    sigma2 = 1,
    subtests = "both"
  )
  expect_equal(tail(forward$changepoints, 1), n)
})

test_that("both keeps the one-segment result for an unrestrictive threshold", {
  set.seed(3)
  data <- rnorm(1000)

  forward <- SVP(
    data,
    gamma = 1e12,
    test = "gaussian_mean",
    sigma2 = 1,
    subtests = "both"
  )
  expect_equal(forward$changepoints, 1000)
})

test_that("SVP validates its public arguments", {
  expect_error(SVP(numeric(), 1), "at least one observation")
  expect_error(SVP(c(1, NA_real_), 1), "finite")
  expect_error(SVP(1:3, 0), "gamma")
  expect_error(SVP(1:3, 1, subtests = "invalid"), "subtests")
  expect_error(SVP(1:3, 1, test = "unknown"), "Unknown test")
  expect_error(SVP(c(1, 0, 2), 1, test = "gamma_rate"),
               "strictly positive")
  expect_error(SVP(1:3, 1, test = "quantile", quantile = 0.9), "0.5")
  expect_error(SVP(1:3, 1, test = "AR1", rho = 0, sigma2 = 0), "sigma2")
  expect_silent(SVP(1:3, 1, test = "gaussian_mean", sigma2 = 0))
})
