## Benchmark stepR::smuceR() against svpChange2::svp_smuce_cpp().
## Run from the package root after installing stepR.

if (!requireNamespace("stepR", quietly = TRUE))
  stop("Install stepR first with install.packages('stepR')")
if (!requireNamespace("svpChange2", quietly = TRUE))
  stop("Install svpChange2 first")

benchmark_smuce_time <- function(
    n_values = c(50L, 100L, 200L, 400L, 800L),
    repetitions = 3L,
    sigma = 1,
    q = 0.8,
    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  out <- vector("list", length(n_values) * repetitions)
  z <- 1L
  for (n in n_values) for (replicate in seq_len(repetitions)) {
    y <- rnorm(n, sd = sigma)
    ## Warm-up and verify equality before recording timings.
    smuce <- system.time(stepR::smuceR(
      y, q = q, family = "gauss", param = sigma
    ))
    cpp <- system.time(svpChange2::svp_smuce_cpp(y, q, sigma^2))
    out[[z]] <- data.frame(
      n = n, replicate = replicate,
      smuce_seconds = unname(smuce[["elapsed"]]),
      svp_smuce_cpp_seconds = unname(cpp[["elapsed"]])
    )
    z <- z + 1L
  }
  do.call(rbind, out)
}

## Example: timings <- benchmark_smuce_time(seed = NULL); print(timings)
