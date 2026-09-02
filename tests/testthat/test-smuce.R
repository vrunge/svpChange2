.smuce_test_partitions <- function(n) {
  lapply(seq_len(2^(n - 1L)) - 1L, function(mask) {
    cuts <- which(as.logical(intToBits(mask)[seq_len(n - 1L)]))
    as.integer(c(cuts, n))
  })
}

.smuce_test_partition_cost <- function(y, ends, gamma, sigma2, n) {
  starts <- c(1L, head(ends, -1L) + 1L)
  sum(vapply(seq_along(ends), function(i) {
    smuce_cost(y[starts[i]:ends[i]], gamma, sigma2, n)
  }, numeric(1)))
}

test_that("SMUCE validity uses the admissible constant-mean interval", {
  y <- c(-0.2, 0.1, 0.05, -0.1)
  interval <- smuce_theta_interval(y, gamma = 2, sigma2 = 1, n = 10)

  expect_identical(
    valid_SMUCE(y, gamma = 2, sigma2 = 1, n = 10),
    !anyNA(interval)
  )
  expect_true(valid_SMUCE(y, 2, 1, 10, theta = mean(y)))
  expect_false(valid_SMUCE(y, 2, 1, 10, theta = interval[2] + 1))
  expect_false(valid_SMUCE(y, 2, 1, 10, theta = c(mean(y), 0)))
})

test_that("R SMUCE dynamic programming is globally optimal", {
  set.seed(202)
  y <- c(rnorm(4), 3 + rnorm(4))
  gamma <- 1.5
  sigma2 <- 1
  n <- length(y)

  parts <- .smuce_test_partitions(n)
  valid <- vapply(parts, function(ends) {
    starts <- c(1L, head(ends, -1L) + 1L)
    all(vapply(seq_along(ends), function(i) {
      is.finite(smuce_cost(y[starts[i]:ends[i]], gamma, sigma2, n))
    }, logical(1)))
  }, logical(1))
  segment_counts <- vapply(parts, length, integer(1))
  min_segments <- min(segment_counts[valid])
  candidates <- parts[valid & segment_counts == min_segments]
  costs <- vapply(candidates, .smuce_test_partition_cost, numeric(1),
                  y = y, gamma = gamma, sigma2 = sigma2, n = n)

  actual <- svp_smuce(y, gamma, sigma2)
  actual_cost <- .smuce_test_partition_cost(
    y, actual, gamma, sigma2, n
  )

  expect_equal(length(actual), min_segments)
  expect_equal(actual_cost, min(costs), tolerance = 1e-10)
  expect_equal(tail(actual, 1), n)
})

test_that("R and C++ SMUCE implementations agree", {
  set.seed(203)
  y <- c(rnorm(5), 2 + rnorm(5))

  for (sigma2 in c(0.5, 1, 2)) {
    for (gamma in c(0, 0.75, 1.5, 3)) {
      expect_equal(
        svp_smuce(y, gamma, sigma2),
        svp_smuce_cpp(y, gamma, sigma2)
      )
    }
  }
})
