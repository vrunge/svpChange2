#' Internal helper: admissible interval of constant means
#' @param y Numeric observations.
#' @param gamma SMUCE threshold.
#' @param sigma2 Known Gaussian variance.
#' @param n Total series length.
#' @return Numeric lower and upper admissible bounds.
#' @keywords internal
smuce_theta_interval <- function(y, gamma, sigma2 = 1, n = length(y)) {
  lo <- -Inf
  hi <- Inf
  cs <- c(0, cumsum(y))
  for (u in seq_along(y)) {
    for (v in u:length(y)) {
      m <- v - u + 1
      centre <- (cs[v + 1] - cs[u]) / m
      r <- sqrt(sigma2 / m) *
        (gamma + sqrt(2 * log(exp(1) * n / m)))
      lo <- max(lo, centre - r)
      hi <- min(hi, centre + r)
      if (lo > hi) {
        return(c(NA_real_, NA_real_))
      }
    }
  }
  c(lo, hi)
}

#' Constrained Gaussian cost for a SMUCE-valid segment
#' @param y Numeric observations in the candidate segment.
#' @param gamma SMUCE threshold.
#' @param sigma2 Known Gaussian variance.
#' @param n Total series length.
#' @return Numeric constrained Gaussian residual cost.
#' @export
smuce_cost <- function(y, gamma, sigma2 = 1, n = length(y)) {
  interval <- smuce_theta_interval(y, gamma, sigma2, n)
  if (anyNA(interval)) {
    return(Inf)
  }
  theta <- min(max(mean(y), interval[1]), interval[2])
  sum((y - theta)^2) / sigma2
}

#' SVP with SMUCE validity and constrained Gaussian cost
#' @param y Numeric observations.
#' @param gamma SMUCE threshold.
#' @param sigma2 Known Gaussian variance.
#' @return Integer segment-end indices.
#' @export
svp_smuce <- function(y, gamma, sigma2 = 1) {
  n <- length(y)
  K <- rep(Inf, n + 1)
  C <- rep(Inf, n + 1)
  prev <- integer(n + 1)
  K[1] <- 0
  C[1] <- 0
  for (t in seq_len(n)) {
    for (s in 0:(t - 1L)) {
      cc <- smuce_cost(y[(s + 1L):t], gamma, sigma2, n)
      if (is.finite(cc)) {
        kk <- K[s + 1L] + 1L
        val <- C[s + 1L] + cc
        if (kk < K[t + 1L] || (kk == K[t + 1L] && val < C[t + 1L])) {
          K[t + 1L] <- kk
          C[t + 1L] <- val
          prev[t + 1L] <- s
        }
      }
    }
  }
  ends <- integer(0)
  t <- n
  while (t > 0L) {
    ends <- c(t, ends)
    t <- prev[t + 1L]
  }
  as.integer(ends)
}

### here we used 2*gamma to be coherent with FOCUS
