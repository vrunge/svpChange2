test_that("SVP dispatches every former cost test", {
  set.seed(123)
  y <- rep(c(0, 1, -0.5), each = 20) + rnorm(60)
  tests <- c("quantileExact", "quantile", "varCost",
             "WilcoxonCost", "MedianMoodCost")
  for (test_name in tests) {
    result <- SVP(y, gamma = 10, test = test_name,
                  prune_after_if_unvalid = TRUE,
                  prune_before_if_invalid = FALSE, quantile = 0.1)
    expect_true(is.list(result), info = test_name)
    expect_true(all(c("changepoints", "nb", "R") %in% names(result)),
                info = test_name)
    expect_equal(tail(result$changepoints, 1), length(y), info = test_name)
  }
})

test_that("quantile is only accepted on its valid parameter range", {
  y <- rnorm(30)
  expect_error(SVP(y, 5, "quantile", quantile = 0))
  expect_error(SVP(y, 5, "quantileExact", quantile = 1))
  expect_error(SVP(y, 5, "quantile", quantile = NA_real_))
  expect_silent(SVP(y, 5, "quantile", quantile = 0.1))
})

test_that("rank cost tests retain exact tie handling", {
  wilcoxon_valid <- function(x, gamma) {
    if (length(x) < 2) return(TRUE)
    best <- 0
    for (u in seq_len(length(x) - 1)) {
      score <- outer(
        x[seq_len(u)], x[(u + 1):length(x)],
        function(a, b) ifelse(a < b, 0.5, ifelse(a > b, -0.5, 0))
      )
      best <- max(best, abs(sum(score)))
    }
    best < gamma
  }

  mood_valid <- function(x, gamma) {
    n <- length(x)
    if (n < 2) return(TRUE)
    med <- sort(x, partial = n %/% 2 + 1)[n %/% 2 + 1]
    below <- x < med
    above <- x > med
    total_below <- sum(below)
    total_above <- sum(above)
    best <- 0
    for (u in seq_len(n - 1)) {
      a11 <- sum(below[seq_len(u)])
      a12 <- sum(above[seq_len(u)])
      a21 <- total_below - a11
      a22 <- total_above - a12
      nA <- a11 + a12
      nB <- a21 + a22
      if (nA > 0 && nB > 0 && total_below > 0 && total_above > 0) {
        determinant <- a11 * a22 - a12 * a21
        best <- max(
          best,
          (nA + nB) * determinant^2 /
            (nA * nB * total_below * total_above)
        )
      }
    }
    best < gamma
  }

  set.seed(1)
  data <- round(c(rnorm(20), rnorm(20, 1)), 1)
  for (gamma in c(3, 7, 15)) {
    expect_equal(
      SVP(data, gamma, "WilcoxonCost")$changepoints,
      svp0(data, gamma, wilcoxon_valid)$changepoints
    )
    expect_equal(
      SVP(data, gamma, "MedianMoodCost")$changepoints,
      svp0(data, gamma, mood_valid)$changepoints
    )
  }
})

test_that("integrated AR1 cost tests retain rho metadata", {
  set.seed(321)
  y <- numeric(80)
  mu <- rep(c(0, 1), each = 40)
  y[1] <- mu[1] + rnorm(1) / sqrt(1 - 0.6^2)
  for (i in 2:length(y))
    y[i] <- mu[i] + 0.6 * (y[i - 1] - mu[i - 1]) + rnorm(1)

  for (test_name in c("AR1", "AR1Profile", "AR1Focus")) {
    result <- SVP(y, 10, test_name, rho = 0.6, sigma2 = 1)
    expect_equal(result$rho, 0.6, info = test_name)
    expect_equal(result$sigma2, 1, info = test_name)
    expect_equal(tail(result$changepoints, 1), length(y), info = test_name)
  }
})
