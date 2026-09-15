project_root <- normalizePath(testthat::test_path("..", ".."))
setwd(project_root)
source(file.path(project_root, "simulations", "power_ar1", "run_power.R"))

test_that("AR(1) generator is deterministic and has the documented scale", {
  first <- simulate_ar1_null(600L, rho = 0.8, seed = 930001L)
  second <- simulate_ar1_null(600L, rho = 0.8, seed = 930001L)
  expect_identical(first, second)
  expect_lt(abs(mean(vapply(
    seq_len(100L),
    function(i) var(simulate_ar1_null(600L, 0.8, 930100L + i)),
    numeric(1)
  )) - 1), 0.05)
})

test_that("AR(1) methods contain only the four requested algorithms", {
  skip_if_not_installed("DeCAFS")
  calibration <- run_ar1_calibration(
    n = 96L,
    calibration_reps = 20L,
    validation_reps = 20L,
    calibration_seed = 930100L,
    validation_seed = 940100L,
    decafs_grid = seq(4, 10, 2),
    approximate_grid = seq(1, 3, 1),
    svp_grid = seq(2, 6, 2),
    workers = 1L,
    write_outputs = FALSE
  )
  y <- simulate_ar1_null(96L, 0.8, 930200L)
  fits <- fit_ar1_methods(y, rho = 0.8, calibration = calibration)
  expect_length(fits, 4L)
  expect_identical(
    names(fits),
    c(
      "PELT AR1 approximate", "PELT inflated", "DeCAFS AR1",
      ar1_method_labels(calibration)[[4L]]
    )
  )
  expect_false(grepl("subtests\\s*=\\s*['\"]both['\"]",
                     paste(readLines(file.path(
                       project_root, "simulations", "power_ar1", "run_power.R"
                     )), collapse = "\n")))
})

test_that("AR(1) calibration selection is reproducible", {
  skip_if_not_installed("DeCAFS")
  args <- list(
    n = 96L,
    calibration_reps = 20L,
    validation_reps = 20L,
    calibration_seed = 930300L,
    validation_seed = 940300L,
    decafs_grid = seq(4, 10, 2),
    approximate_grid = seq(1, 3, 1),
    svp_grid = seq(2, 6, 2),
    workers = 1L,
    write_outputs = FALSE
  )
  first <- do.call(run_ar1_calibration, args)
  second <- do.call(run_ar1_calibration, args)
  expect_identical(first$selected, second$selected)
  expect_identical(first$table, second$table)
})

test_that("small AR(1) power run has four methods and corrected metrics", {
  skip_if_not_installed("DeCAFS")
  calibration <- run_ar1_calibration(
    n = 96L,
    calibration_reps = 20L,
    validation_reps = 20L,
    calibration_seed = 930400L,
    validation_seed = 940400L,
    decafs_grid = seq(4, 10, 2),
    approximate_grid = seq(1, 3, 1),
    svp_grid = seq(2, 6, 2),
    workers = 1L,
    write_outputs = FALSE
  )
  results <- run_ar1_power(
    n = 96L, jump_sizes = c(1, 2), reps = 2L, workers = 1L,
    calibration = calibration
  )
  expect_equal(nrow(results), 64L)
  expect_setequal(
    unique(as.character(results$algorithm)),
    ar1_method_labels(calibration)
  )
  expect_true(all(results$NumSegments >= 1L))
  expect_equal(
    paper_metrics(96L, c(5L, 96L), tolerance = 5L, terminal = 96L)[["F1"]],
    0
  )
})
