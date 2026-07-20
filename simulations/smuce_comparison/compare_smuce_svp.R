## Compare SMUCE from stepR with the multiscale Gaussian SVP.
##
## This script follows the equivalence discussed in
## smuce_comparison/SVP_SMUCE_report.tex.  It deliberately uses independent
## Gaussian observations, known sigma^2, the same multiscale threshold q, and
## no directional candidate pruning.  In this setting the ordinary Gaussian
## segment MLE is admissible, so the data-only SVP and SMUCE partitions agree.

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
      return(sort(unique(c(which(diff(fitted_values) != 0), n))))
    }
  }

  if (!is.null(fit$rightEnd)) {
    return(sort(unique(c(as.integer(fit$rightEnd), n))))
  }
  if (!is.null(fit$step$rightEnd)) {
    return(sort(unique(c(as.integer(fit$step$rightEnd), n))))
  }

  stop("Could not extract SMUCE segment ends; inspect str(fit) for this stepR version")
}

compare_smuce_svp <- function(seed = 20260720,
                              n_per_segment = 100,
                              means = c(0, 2, -1, 1),
                              sigma = 1,
                              alpha = 0.1,
                              plot = TRUE) {
  stopifnot(length(means) >= 2, n_per_segment >= 10, sigma > 0,
            alpha > 0, alpha < 1)
  set.seed(seed)
  n <- length(means) * n_per_segment
  truth <- rep(means, each = n_per_segment)
  y <- truth + stats::rnorm(n, sd = sigma)

  ## stepR calibrates q for the global multiscale SMUCE constraint.  Passing
  ## exactly this q to SVP's multiscale threshold is the important step.
  q <- stepR::thresh.smuceR(n, alpha = alpha)

  smuce_fit <- stepR::smuceR(
    y,
    x = seq_len(n),
    family = "gauss",
    q = q,
    sd = sigma
  )

  ## No directional pruning: this compares the global constrained objective,
  ## rather than an online pruning approximation.
  svp_fit <- svpChange2::SVP(
    y,
    gamma = 0,
    test = "gaussian_mean",
    prune_after_if_unvalid = FALSE,
    prune_before_if_invalid = FALSE,
    use_multiscale_gamma = TRUE,
    sigma2 = sigma^2,
    q_alpha_n = q
  )

  smuce_cp <- extract_step_boundaries(smuce_fit, n)
  svp_cp <- sort(unique(as.integer(svp_fit$changepoints)))
  equal_partition <- identical(smuce_cp, svp_cp)

  result <- list(
    n = n,
    sigma = sigma,
    alpha = alpha,
    q = q,
    smuce_changepoints = smuce_cp,
    svp_changepoints = svp_cp,
    equal_partition = equal_partition,
    smuce_fit = smuce_fit,
    svp_fit = svp_fit
  )

  if (plot) {
    pdf("smuce_svp_comparison.pdf", width = 10, height = 5)
    plot(y, type = "p", pch = 16, cex = 0.45,
         main = sprintf("SMUCE versus multiscale SVP (q = %.3f)", q),
         xlab = "Index", ylab = "Observation")
    lines(truth, col = "grey40", lwd = 2)
    abline(v = smuce_cp[-length(smuce_cp)], col = "firebrick", lwd = 2)
    abline(v = svp_cp[-length(svp_cp)], col = "navy", lwd = 1, lty = 2)
    legend("topright", legend = c("truth", "SMUCE", "SVP"),
           col = c("grey40", "firebrick", "navy"),
           lwd = c(2, 2, 1), lty = c(1, 1, 2), bty = "n")
    dev.off()
  }

  if (!equal_partition) {
    warning(
      "Partitions differ. This is a useful diagnostic: the report's exact " ,
      "equivalence requires the same constrained segment-parameter fit, " ,
      "not only the same multiscale threshold."
    )
  }
  result
}

## Run explicitly with: source("simulations/smuce_comparison/compare_smuce_svp.R")
## answer <- compare_smuce_svp()
## stopifnot(answer$equal_partition)

