##  GPL-3 License
## Copyright (c) 2025 Alexandre Combeau

#' Validity Test Based on FOCuS
#'
#' @title Validity Test Based on FOCuS
#' @description Checks whether a given segment is valid (i.e., has no changepoint).
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric value used as a threshold in the validation function and as a penalty for each segment.
#' @return A logical value indicating whether the segment is considered valid (TRUE if no changepoint detected).
#' @export
valid_FOCUS <- function(y, gamma)
{
  .focus_valid_cpp(y, gamma, check_all_prefixes = TRUE)
}

################################################
#' Validity Test Based on FOCuS last statistic
#'
#' @title Validity Test Based on FOCuS last statistic
#' @description Checks whether a given segment is valid (i.e., has no changepoint).
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric value used as a threshold in the validation function and as a penalty for each segment.
#' @return A logical value indicating whether the segment is considered valid (TRUE if no changepoint detected).
#' @export
valid_FOCUS_last <- function(y, gamma)
{
  .focus_valid_cpp(y, gamma, check_all_prefixes = FALSE)
}

################################################

#' Validity Test Based on Variance
#'
#' @title Validity Test Based on Variance
#' @description Checks if the variance of a segment is lower than or equal to a threshold.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the maximum allowed variance.
#' @return TRUE if the variance is less than or equal to gamma, or if the segment is very short.
#' @export
valid_VAR <- function(y, gamma)
{
  if(length(y) < 6)
  {
    return(TRUE)
  }
  sum((y - mean(y))^2) <= gamma
}

################################################

#' Validity Test Based on Range
#'
#' @title Validity Test Based on Range
#' @description Checks whether the range (max - min) of the segment is below a threshold.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the maximum allowed range.
#' @return TRUE if the range is less than or equal to gamma.
#' @export
valid_RANGE <- function(y, gamma)
{
  max(y) - min(y) <= gamma
}

################################################

#' Validity Test Based on Trimmed Range
#'
#' @title Validity Test Based on Trimmed Range
#' @description Applies a slack range test using trimmed minimum and maximum (ignores smallest and largest 3 values).
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the trimmed range.
#' @return TRUE if the trimmed range is less than or equal to gamma, or if the segment is too short.
#' @export
valid_RANGE_SLACK <- function(y, gamma)
{
  if(length(y) < 7)
  {
    return(TRUE)
  }
  sortedY <- sort(y)
  sortedY[length(y) - 3] - sortedY[3] <= gamma
}

################################################

#' Validity Test Based on Interquantile Range
#'
#' @title Validity Test Based on Interquantile Range
#' @description Checks whether the interquantile range (95% - 5%) is below a threshold.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the interquantile range.
#' @return TRUE if the interquantile range is less than or equal to gamma.
#' @export
valid_QUANTILE <- function(y, gamma)
{
  unname(diff(quantile(y, probs = c(0.05, 0.95)))) <= gamma
}

################################################
#' Optimal Partitioning Cost Test
#'
#' @title Optimal Partitioning Cost Test
#' @description Tests whether the total cost of the segment is smaller than the best two-part penalized segmentation cost.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for penalizing the introduction of a new segment
#' @return TRUE if the segment cannot be split into two parts with lower (penalized) cost.
#' @export
valid_OP <- function(y, gamma)
{
  len <- length(y)
  if(len == 1){return(TRUE)}
  total <- sum((y - mean(y))^2)
  val <- Inf
  for(i in 2:len)
  {
    segment1 <- y[1:(i-1)]
    segment2 <- y[i:len]
    mean_seg1 <- mean(segment1)
    mean_seg2 <- mean(segment2)
    temp <- (sum((segment1 - mean_seg1)^2) + sum((segment2 - mean_seg2)^2))
    if(temp < val){val <- temp}
  }
  return(test = (total < (val + 2*gamma)))
}

################################################
#' SMUCE multiscale validity test for a constant Gaussian segment
#'
#' @description Tests a segment against the SMUCE multiscale constraint. The
#' candidate mean defaults to the segment mean; supply `theta` to test another
#' constant mean. `gamma` is interpreted as the SMUCE q threshold.
#' @param y Numeric observations in the candidate segment.
#' @param gamma SMUCE threshold q (not q squared).
#' @param sigma2 Known Gaussian variance.
#' @param n Total series length used in the multiscale penalty.
#' @param theta Optional candidate constant mean. If omitted, the function
#' searches for any admissible mean by intersecting all subinterval constraints.
#' @return Logical validity indicator.
#' @export
valid_SMUCE <- function(y, gamma, sigma2 = 1, n = length(y), theta = NULL)
{
  stopifnot(length(y) > 0, sigma2 > 0, gamma >= 0, n >= length(y))
  cs <- c(0, cumsum(y))
  lower <- -Inf
  upper <- Inf
  for (u in seq_along(y)) {
    for (v in u:length(y)) {
      m <- v - u + 1
      centre <- (cs[v + 1] - cs[u]) / m
      radius <- sqrt(sigma2 / m) *
        (gamma + sqrt(2 * log(exp(1) * n / m)))
      lower <- max(lower, centre - radius)
      upper <- min(upper, centre + radius)
      if (lower > upper) return(FALSE)
    }
  }
  if (is.null(theta)) return(TRUE)
  isTRUE(theta >= lower && theta <= upper)
}

#' Internal helper: admissible interval of constant means
#' @param y Numeric observations.
#' @param gamma SMUCE threshold.
#' @param sigma2 Known Gaussian variance.
#' @param n Total series length.
#' @return Numeric lower and upper admissible bounds.
#' @keywords internal
smuce_theta_interval <- function(y, gamma, sigma2 = 1, n = length(y)) {
  lo <- -Inf; hi <- Inf; cs <- c(0, cumsum(y))
  for (u in seq_along(y)) for (v in u:length(y)) {
    m <- v-u+1; centre <- (cs[v+1]-cs[u])/m
    r <- sqrt(sigma2/m) *
      (gamma + sqrt(2*log(exp(1)*n/m)))
    lo <- max(lo, centre-r); hi <- min(hi, centre+r)
    if (lo > hi) return(c(NA_real_, NA_real_))
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
  if (anyNA(interval)) return(Inf)
  theta <- min(max(mean(y), interval[1]), interval[2])
  sum((y-theta)^2) / sigma2
}

#' SVP with SMUCE validity and constrained Gaussian cost
#' @param y Numeric observations.
#' @param gamma SMUCE threshold.
#' @param sigma2 Known Gaussian variance.
#' @return Integer segment-end indices.
#' @export
svp_smuce <- function(y, gamma, sigma2 = 1) {
  n <- length(y); K <- rep(Inf, n + 1); C <- rep(Inf, n + 1)
  prev <- integer(n + 1); K[1] <- 0; C[1] <- 0
  for (t in seq_len(n)) for (s in 0:(t - 1L)) {
    cc <- smuce_cost(y[(s + 1L):t], gamma, sigma2, n)
    if (is.finite(cc)) {
      kk <- K[s + 1L] + 1L; val <- C[s + 1L] + cc
      if (kk < K[t + 1L] || (kk == K[t + 1L] && val < C[t + 1L])) {
        K[t + 1L] <- kk; C[t + 1L] <- val; prev[t + 1L] <- s
      }
    }
  }
  ends <- integer(0); t <- n
  while (t > 0L) { ends <- c(t, ends); t <- prev[t + 1L] }
  as.integer(ends)
}

### here we used 2*gamma to be coherent with FOCUS
