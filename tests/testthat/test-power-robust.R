library(testthat)
library(svpChange2)

package_root <- normalizePath(
  file.path(testthat::test_path(), "..", ".."), mustWork = TRUE
)
old_working_directory <- setwd(package_root)
sys.source(
  file.path(package_root, "simulations", "power_robust", "run_power.R"),
  envir = .GlobalEnv
)
setwd(old_working_directory)

test_that("Median-Mood uses the Pearson chi-square one df threshold", {
  n <- 1000L
  oracle_segments <- 8L
  alpha <- 0.01
  alpha_s <- 1 - (1 - alpha)^(1 / (n / oracle_segments - 1))
  expected <- qchisq(1 - alpha_s, df = 1)
  expect_equal(mood_threshold(n, oracle_segments, alpha), expected)
  expect_false(isTRUE(all.equal(
    mood_threshold(n, oracle_segments, alpha),
    qchisq(1 - alpha_s, df = oracle_segments)
  )))
})

test_that("Student-t generation is deterministic under the study seed", {
  first <- student_null(100L, 20260915L)
  second <- student_null(100L, 20260915L)
  expect_identical(first, second)
  expect_length(first, 100L)
})

test_that("robust fitting exposes only the four canonical methods", {
  calibration <- list(selected = list(
    rfpop_constant = 3,
    wilcoxon_constant = 1.5,
    mood_alpha = 0.01
  ))
  set_power_seed(20260915L)
  y <- ts_generator(
    chpts = c(25L, 50L, 75L, 100L),
    parameters = c(0, 0.5, 0, 0.5),
    type = "student", df = 2, scale = 1
  )
  fit <- fit_robust_methods(y, oracle_segments = 4L, calibration)
  expect_setequal(names(fit$boundaries), robust_method_labels(calibration))
  expect_false(any(grepl("both|multiscale", names(fit$boundaries), ignore.case = TRUE)))
  expect_true(grepl("SVP MedianMood / right / alpha = 0.01",
                    robust_method_labels(calibration)[[3L]], fixed = TRUE))
})

test_that("the heavy-tail study contains right-mode SVP calls only", {
  source_lines <- readLines(file.path(
    package_root, "simulations", "power_robust", "run_power.R"
  ))
  expect_true(any(grepl('subtests = "right"', source_lines, fixed = TRUE)))
  expect_false(any(grepl('subtests = "both"', source_lines, fixed = TRUE)))
  expect_false(any(grepl("multiscale", source_lines, ignore.case = TRUE)))
})

test_that("calibration smoke run is reproducible", {
  args <- list(
    n = 80L, calibration_reps = 10L, validation_reps = 8L,
    rfpop_target = 0.5, rfpop_grid = c(2, 4),
    wilcoxon_grid = c(0.75, 1.5), mood_grid = c(0.01, 0.1),
    workers = 1L, write_outputs = FALSE
  )
  first <- do.call(run_robust_calibration, args)
  second <- do.call(run_robust_calibration, args)
  expect_identical(first$selected, second$selected)
  expect_identical(first$table, second$table)
})
