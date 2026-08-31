###############################################
#############    tsGenerator   ################
###############################################

#' tsGenerator
#'
#' @description Generating univariate time series for multiple change-point
#' detection based on uni-parametric models of the exponential family. The
#' \code{"gaussAR1"} model generates pure AR(1) noise around the
#' piecewise-constant signal, matching the DeCAFS model with \code{sdEta = 0}.
#' @param chpts a vector of increasing change-point indices (the last value is
#' data length)
#' @param parameters vector of successive segment parameters (as many parameters
#' as values in \code{chpts} vector); for \code{"exp"} these are rates and for
#' \code{"variance"} these are standard deviations
#' @param sdNoise (types \code{"gauss"} and \code{"gaussAR1"}) standard
#' deviation of each Gaussian innovation
#' @param df (type \code{"student"}) positive degrees of freedom for the
#' Student-t noise
#' @param scale (type \code{"student"}) non-negative scale of the Student-t
#' noise
#' @param rho (type \code{"gaussAR1"}) AR(1) coefficient for the Gaussian noise
#' @param nbTrials (type \code{"binom"}) number of trials
#' @param nbSuccess (type \code{"negbin"}) number of successes
#' @param type the model: \code{"gauss"}, \code{"gaussAR1"},
#' \code{"student"}, \code{"exp"}, \code{"poisson"}, \code{"geom"},
#' \code{"bern"}, \code{"binom"}, \code{"negbin"}, \code{"variance"}
#' @return a univariate time series following the chosen model type and
#' parameters
#' @details For \code{type = "gaussAR1"}, the generated series is
#' \code{y[t] = parameters[t] + e[t]}, where
#' \code{e[t] = rho * e[t - 1] + eta[t]} and
#' \code{eta[t] ~ N(0, sdNoise^2)}. The first residual is drawn from the
#' stationary distribution. Consequently, the marginal noise standard
#' deviation is \code{sdNoise / sqrt(1 - rho^2)}, not \code{sdNoise}. To
#' obtain marginal standard deviation \code{target_sd}, use
#' \code{sdNoise = target_sd * sqrt(1 - rho^2)}. The \code{"student"} model
#' uses \code{parameters + scale * T}, where \code{T} has a Student-t
#' distribution with \code{df} degrees of freedom; for \code{df <= 2}, its
#' variance is not finite.
#' @examples
#' set.seed(1)
#' # Independent Gaussian noise with a few mean changes.
#' tsGenerator(
#'   chpts = c(50, 100, 150, 200),
#'   parameters = c(0, 2, -1, 1),
#'   sdNoise = 1,
#'   type = "gauss"
#' )
#'
#' # Pure AR(1) Gaussian noise around the same type of signal.
#' tsGenerator(
#'   chpts = c(50, 100, 150, 200),
#'   parameters = c(0, 2, -1, 1),
#'   sdNoise = 1,
#'   rho = 0.8,
#'   type = "gaussAR1"
#' )
#'
#' # Heavy-tailed Student-t noise.
#' tsGenerator(
#'   chpts = c(50, 100),
#'   parameters = c(0, 2),
#'   df = 2,
#'   scale = 1,
#'   type = "student"
#' )
#'
#' # Other supported distributions.
#' tsGenerator(chpts = c(50, 100), parameters = c(2, 7), type = "exp")
#' tsGenerator(chpts = c(50, 100), parameters = c(3, 5), type = "poisson")
#' tsGenerator(chpts = c(50, 100), parameters = c(0.6, 0.3), type = "geom")
#' tsGenerator(chpts = c(50, 100), parameters = c(0.7, 0.2), type = "bern")
#' tsGenerator(
#'   chpts = c(50, 100), parameters = c(0.7, 0.3), nbTrials = 5,
#'   type = "binom"
#' )
#' tsGenerator(
#'   chpts = c(50, 100), parameters = c(0.4, 0.7), nbSuccess = 10,
#'   type = "negbin"
#' )
#' tsGenerator(
#'   chpts = c(50, 70, 120, 200), parameters = c(0, 3, -1, 1),
#'   type = "gauss"
#' )
#' tsGenerator(
#'   chpts = c(50, 100, 180), parameters = c(3, 1, 6), type = "variance"
#' )
tsGenerator <- function(chpts = 100,
                        parameters = 0.5,
                        sdNoise = 1,
                        rho = 0,
                        nbTrials = 10,
                        nbSuccess = 10,
                        type = "gauss",
                        df = 2,
                        scale = 1) {
  ############
  ### STOP ###
  ############

  if (!is.numeric(chpts) || length(chpts) == 0 || any(!is.finite(chpts))) {
    stop("chpts must be a non-empty finite numeric vector")
  }
  if (any(
    abs(chpts - round(chpts)) >
      sqrt(.Machine$double.eps) * pmax(1, abs(chpts))
  )) {
    stop("chpts must contain integer positions")
  }
  chpts <- round(chpts)
  if (!all(chpts > 0)) {
    stop("chpts values are not all positives")
  }
  if (is.unsorted(chpts, strictly = TRUE)) {
    stop(paste(
      "chpts should be a strictly increasing vector of",
      "change-point positions (indices)"
    ))
  }

  if (!is.numeric(parameters) || any(!is.finite(parameters))) {
    stop("parameters must be finite numeric values")
  }
  if (length(chpts) != length(parameters)) {
    stop("chpts and parameters vectors are of different size")
  }

  allowed.types <- c(
    "gauss", "gaussAR1", "student", "exp", "poisson", "geom", "bern",
    "binom", "negbin", "variance"
  )
  if (
    !is.character(type) || length(type) != 1 || is.na(type) ||
      !type %in% allowed.types
  ) {
    stop("type must be one of: ", paste(allowed.types, collapse = ", "))
  }

  ###################################
  ### Distribution specific stops ###
  ###################################

  if (type %in% c("gauss", "gaussAR1")) {
    if (length(sdNoise) != 1 || !is.numeric(sdNoise) || !is.finite(sdNoise)) {
      stop("sdNoise must be one finite numeric value")
    }
    if (sdNoise < 0) {
      stop("sdNoise cannot be negative")
    }
    if (
      type == "gaussAR1" &&
        (length(rho) != 1 || !is.numeric(rho) || !is.finite(rho) ||
          abs(rho) >= 1)
    ) {
      stop("rho must be one finite numeric value strictly between -1 and 1")
    }
    if (
      type == "gauss" &&
        (!is.numeric(rho) || length(rho) != 1 || !is.finite(rho) || rho != 0)
    ) {
      stop("rho is only supported for the Gaussian AR(1) model")
    }
  } else if (
    length(rho) != 1 || !is.numeric(rho) || !is.finite(rho) || rho != 0
  ) {
    stop("rho is only supported for the Gaussian AR(1) model")
  }

  if (type == "student") {
    if (length(df) != 1 || !is.numeric(df) || !is.finite(df) || df <= 0) {
      stop("df must be one finite positive numeric value")
    }
    if (
      length(scale) != 1 || !is.numeric(scale) || !is.finite(scale) || scale < 0
    ) {
      stop("scale must be one finite non-negative numeric value")
    }
  }

  if (type == "binom") {
    if (
      length(nbTrials) != 1 || !is.numeric(nbTrials) || !is.finite(nbTrials)
    ) {
      stop("nbTrials must be one finite numeric value")
    }
    if ((nbTrials %% 1 != 0) || (nbTrials <= 0)) {
      stop("nbTrials cannot be non-positive or non-integer")
    }
  }

  if (type == "negbin") {
    if (
      length(nbSuccess) != 1 || !is.numeric(nbSuccess) || !is.finite(nbSuccess)
    ) {
      stop("nbSuccess must be one finite numeric value")
    }
    if ((nbSuccess %% 1 != 0) || (nbSuccess <= 0)) {
      stop("nbSuccess cannot be non-positive or non-integer")
    }
  }

  #############################
  ### parameter constraints ###
  #############################

  if (type == "variance") {
    if (min(parameters) < 0) {
      stop("no negative parameters (variances) possible for Variance model")
    }
  }
  if (type == "exp") {
    if (min(parameters) <= 0) {
      stop(
        "parameters must be strictly positive rates for the exponential",
        " model"
      )
    }
  }
  if (type == "poisson") {
    if (min(parameters) <= 0) {
      stop("parameters must be strictly positive means for the Poisson model")
    }
  }

  if (type == "bern" || type == "binom") {
    if (any(parameters > 1) || any(parameters < 0)) {
      stop("parameters should be probabilities between 0 and 1 (included)")
    }
  }

  if (type == "geom" || type == "negbin") {
    if (any(parameters > 1) || any(parameters <= 0)) {
      stop("parameters should be probabilities between 0 and 1 (0 excluded)")
    }
  }

  ############  data generation   ############

  n <- chpts[length(chpts)]
  repetition <- c(chpts[1], diff(chpts))
  mu <- rep(parameters, repetition)

  y <- switch(type,
    gauss = rnorm(n, mean = mu, sd = sdNoise),
    student = mu + scale * stats::rt(n, df = df),
    gaussAR1 = {
      innovations <- rnorm(n, sd = sdNoise)
      y <- numeric(n)
      y[1] <- mu[1] + innovations[1] / sqrt(1 - rho^2)
      if (n > 1) {
        for (index in 2:n) {
          y[index] <- mu[index] + rho * (y[index - 1] - mu[index - 1]) +
            innovations[index]
        }
      }
      y
    },
    variance = rnorm(n, mean = 0, sd = mu),
    exp = rexp(n = n, rate = mu),
    poisson = rpois(n = n, lambda = mu),
    geom = rgeom(n = n, prob = mu) + 1,
    bern = rbinom(n = n, size = 1, prob = mu),
    binom = rbinom(n = n, size = nbTrials, prob = mu),
    negbin = rnbinom(n = n, size = nbSuccess, prob = mu)
  )
  y
}
