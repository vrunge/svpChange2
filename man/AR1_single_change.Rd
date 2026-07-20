\name{AR1_single_change}
\alias{AR1_single_change}
\title{Single Mean Change Test for an AR(1) Process}
\usage{
AR1_single_change(
  data,
  gamma,
  rho = NA_real_,
  sigma2 = 1,
  profile_sigma = FALSE
)
}
\arguments{
  \item{data}{A numeric AR(1) series with at least three observations.}
  \item{gamma}{The FOCuS validity threshold.}
  \item{rho}{An optional known AR(1) coefficient. If missing, it is estimated
  robustly.}
  \item{sigma2}{A positive innovation variance used for the known-variance
  likelihood-ratio statistic.}
  \item{profile_sigma}{If \code{TRUE}, profile the innovation variance and use
  the log RSS-ratio statistic.}
}
\value{
A list with the coefficient and variance used, test statistic, null and
alternative residual sums of squares, estimated changepoint, and a logical
validity indicator.
}
\description{
Computes the exact fixed-\code{rho} likelihood-ratio scan for a change in the
marginal mean of an AR(1) process. The residual crossing the candidate change
uses both segment means.
}
\references{
Chakar, S., Lebarbier, E., Levy-Leduc, C. and Robin, S. (2017).
A robust approach for estimating change-points in the mean of an AR(1)
process. \emph{Bernoulli}, 23(2), 1408--1447.
}
