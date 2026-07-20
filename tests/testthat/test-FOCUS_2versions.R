

library(testthat)
library(svpChange2)

########## test svp0 result = SVP result with prune_after_if_unvalid == TRUE ##########
########## test svp0 result = SVP result with prune_after_if_unvalid == TRUE ##########
########## test svp0 result = SVP result with prune_after_if_unvalid == TRUE ##########

test_that("test svp0 result = SVP result with prune_after_if_unvalid == TRUE",
          {
            n <- 500
            gap <- 1
            chpts = c(0.1,0.3,0.4,0.45,0.55,0.7,0.75,0.95,1)*n
            data <- tsGenerator(chpts = chpts,
                                parameters = c(0,gap,0,gap,0,gap,0,gap,0),
                                sdNoise = 1)
            gamma <- 5
            bool <- TRUE

            res_svp0 <- svp0(data,
                             gamma,
                             test = valid_FOCUS, #valid_FOCUS_last
                             prune_after_if_unvalid = bool)

            res_svp <- SVP(data = data,
                           gamma = gamma,
                           test = "gaussian_mean",
                           prune_after_if_unvalid = bool)

            expect_equal(res_svp0$changepoints, res_svp$changepoints)
            expect_equal(res_svp0$R[,1], res_svp$R[,1])
          })

test_that("SVP accepts multiscale gamma",
          {
            n <- 80
            data <- rep(c(0, 1, -0.5, 0.8), each = n / 4) + rnorm(n)
            gamma <- 2 * log(n)

            res <- SVP(
              data = data,
              gamma = gamma,
              test = "gaussian_mean",
              prune_after_if_unvalid = TRUE,
              prune_before_if_invalid = TRUE,
              use_multiscale_gamma = TRUE,
              sigma2 = 1,
              q_alpha_n = 3
            )

            expect_equal(tail(res$changepoints, 1), n)
          })


test_that("svp0 and SVP return the same Gaussian FOCUS result",
          {
            n <- 50
            gap <- 5
            chpts = c(0.1,0.3,0.4,0.45,0.55,0.7,0.75,0.95,1)*n
            data <- tsGenerator(chpts = chpts,
                                parameters = c(0,gap,0,gap,0,gap,0,gap,0),
                                sdNoise = 1)
            bool <- TRUE

            gamma <- 7

            res_svp0 <- svp0(data,
                             gamma,
                             test = valid_FOCUS, #valid_FOCUS_last
                             prune_after_if_unvalid = bool)
            res_svp <- SVP(data = data,
                           gamma = gamma,
                           test = "gaussian_mean",
                           prune_after_if_unvalid = bool)
            expect_equal(res_svp$changepoints, res_svp0$changepoints)
            expect_equal(res_svp$R, res_svp0$R)

          })

test_that("SVP accepts all before/after pruning combinations",
          {
            n <- 80
            data <- rep(c(0, 1, -0.5, 0.8), each = n / 4) + rnorm(n)
            gamma <- 2 * log(n)

            opts <- expand.grid(
              prune_after_if_unvalid = c(FALSE, TRUE),
              prune_before_if_invalid = c(FALSE, TRUE)
            )

            for (i in seq_len(nrow(opts))) {
              res <- SVP(
                data = data,
                gamma = gamma,
                test = "gaussian_mean",
                prune_after_if_unvalid = opts$prune_after_if_unvalid[i],
                prune_before_if_invalid = opts$prune_before_if_invalid[i]
              )

              expect_equal(tail(res$changepoints, 1), n)
            }
          })
