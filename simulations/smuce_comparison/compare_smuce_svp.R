## Compare SMUCE from stepR with the multiscale Gaussian SVP.
##
## This script compares the two methods on a white-noise null with an
## deliberately under-calibrated q. Both methods therefore produce many
## spurious segments; the returned partitions can be compared directly.

if (!requireNamespace("stepR", quietly = TRUE)) {
  stop("Install stepR first with install.packages('stepR')")
}
if (!requireNamespace("svpChange2", quietly = TRUE)) {
  stop("Install svpChange2 first")
}

library(svpChange2)

extract_step_boundaries <- function(fit, n) {
  ## stepR versions expose the fitted step function either through fitted()
  ## or through a stepblock object with rightEnd.  Supporting both keeps this
  ## comparison script usable across current CRAN releases.
  fitted_values <- tryCatch(stats::fitted(fit), error = function(e) NULL)
  if (!is.null(fitted_values)) {
    fitted_values <- as.numeric(fitted_values)
    if (length(fitted_values) == n) {
      return(sort(unique(c(as.integer(which(diff(fitted_values) != 0)),
                           as.integer(n)))))
    }
  }

  if (!is.null(fit$rightEnd)) {
    return(sort(unique(c(as.integer(fit$rightEnd), as.integer(n)))))
  }
  if (!is.null(fit$step$rightEnd)) {
    return(sort(unique(c(as.integer(fit$step$rightEnd), as.integer(n)))))
  }

  stop("Could not extract SMUCE segment ends; inspect str(fit) for this stepR version")
}

compare_smuce_svp <- function(
                              n = 100,
                              sigma = 0.5,
                              q_scale = 0.8,
                              seed = NULL,
                              n_segments = 10) {

  stopifnot(n >= 10, sigma > 0, q_scale > 0)
  if (!is.null(seed)) set.seed(seed)
  ## Random, heterogeneous piecewise-constant Gaussian signal.  No seed is
  ## used by default, so every invocation generates a different rich
  ## partition.
  cuts <- sort(sample(2:(n - 1), n_segments - 1L))
  ends <- c(cuts, n)
  lengths <- diff(c(0, ends))
  means <- cumsum(stats::rnorm(n_segments, mean = 0, sd = 1.5))
  signal <- rep(means, lengths)
  y <- signal + stats::rnorm(n, sd = sigma)

  ## SMUCE reports sqrt(2 T_i - 2 log(en/l)) <= q.  SVP stores the
  ## corresponding squared threshold, hence q_alpha_n = q^2.
  q <- q_scale * stepR::thresh.smuceR(n)

  smuce_fit <- stepR::smuceR(
    y,
    x = seq_len(n),
    family = "gauss",
    q = q,
    param = sigma
  )

  constrained <- svpChange2::svp_smuce(y, q, sigma^2)

  smuce_cp <- extract_step_boundaries(smuce_fit, n)

  list(
    smuce = smuce_cp,
    svp_smuce = constrained
  )
}

#' Run repeated random SMUCE/SVP comparisons.
#' @return A data.frame with one row per simulation and an equality indicator.
run_smuce_replicates <- function(B = 20, n = 120, ...) {
  out <- vector("list", B)
  for (b in seq_len(B)) {
    z <- compare_smuce_svp(n = n, seed = NULL, ...)
    out[[b]] <- data.frame(
      replicate = b,
      equal = identical(z$smuce, z$svp_smuce),
      smuce = paste(z$smuce, collapse = ","),
      svp_smuce = paste(z$svp_smuce, collapse = ","),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

## Run explicitly with: source("simulations/smuce_comparison/compare_smuce_svp.R")
answer <- compare_smuce_svp()
answer
