## Common change-point metrics used in the SVP/PELT power studies.
## LocalizationError is the maximum matched distance when every true change
## point is recovered; otherwise it is NA. CorrectNumCP records exact model
## size recovery and is therefore directly estimable as a proportion.

localization_error <- function(cp_true, cp_est) {
  cp_true <- sort(as.integer(cp_true))
  cp_est <- sort(as.integer(cp_est))
  if (!length(cp_true) || length(cp_true) != length(cp_est)) return(NA_real_)
  max(abs(cp_est - cp_true))
}

change_point_metrics <- function(cp_true, cp_est, tol = 5L) {
  cp_true <- sort(as.integer(cp_true)); cp_est <- sort(as.integer(cp_est))
  if (!length(cp_true) && !length(cp_est))
    return(list(Precision = NA_real_, Recall = NA_real_, F1 = NA_real_,
                LocalizationError = NA_real_, CorrectNumCP = TRUE))
  if (!length(cp_true) || !length(cp_est))
    return(list(Precision = if (length(cp_est)) 0 else NA_real_,
                Recall = if (length(cp_true)) 0 else NA_real_, F1 = NA_real_,
                LocalizationError = NA_real_,
                CorrectNumCP = length(cp_true) == length(cp_est)))
  D <- abs(outer(cp_true, cp_est, "-")); used <- rep(FALSE, length(cp_est))
  distances <- numeric(0)
  for (i in seq_along(cp_true)) {
    ord <- order(D[i, ])
    j <- ord[which(!used[ord])[1L]]
    if (!is.na(j) && D[i, j] <= tol) { used[j] <- TRUE; distances <- c(distances, D[i,j]) }
  }
  tp <- length(distances); fp <- sum(!used); fn <- length(cp_true) - tp
  precision <- if (tp + fp) tp/(tp+fp) else NA_real_
  recall <- if (tp + fn) tp/(tp+fn) else NA_real_
  f1 <- if (is.finite(precision) && is.finite(recall) && precision + recall > 0)
    2*precision*recall/(precision+recall) else NA_real_
  list(Precision = precision, Recall = recall, F1 = f1,
       LocalizationError = localization_error(cp_true, cp_est),
       CorrectNumCP = length(cp_true) == length(cp_est))
}
