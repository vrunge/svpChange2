source(testthat::test_path("..", "..", "simulations", "power_common.R"))

test_that("paper metrics exclude the terminal endpoint", {
  perfect_null <- paper_metrics(10L, 10L, tolerance = 2L, terminal = 10L)
  false_positive <- paper_metrics(
    10L, c(4L, 10L), tolerance = 2L, terminal = 10L
  )

  expect_equal(unname(perfect_null[c("Precision", "Recall", "F1")]), c(1, 1, 1))
  expect_equal(unname(false_positive[c("Precision", "Recall", "F1")]), c(0, 0, 0))
  expect_equal(false_positive[["CorrectNumCP"]], 0)
})

test_that("paper metrics use one-to-one tolerance matching", {
  metric <- paper_metrics(
    c(10L, 20L, 30L), c(9L, 11L, 30L),
    tolerance = 2L, terminal = 30L
  )

  expect_equal(metric[["Precision"]], 0.5)
  expect_equal(metric[["Recall"]], 0.5)
  expect_equal(metric[["F1"]], 0.5)
})
