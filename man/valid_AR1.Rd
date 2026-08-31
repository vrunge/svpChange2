\name{valid_AR1}
\alias{valid_AR1}
\title{AR(1) Mean-Change Validity Test}
\usage{
valid_AR1(
  y,
  gamma,
  rho = NA_real_,
  sigma2 = 1,
  profile_sigma = FALSE
)
}
\arguments{
  \item{y}{A numeric AR(1) segment. At least four observations are needed for
  a non-trivial single-change scan.}
  \item{gamma}{A threshold for the AR(1) likelihood-ratio statistic.}
  \item{rho}{An optional known AR(1) coefficient. If \code{NA}, it is
  estimated robustly from \code{y}.}
  \item{sigma2}{A positive innovation variance used by the known-variance
  statistic.}
  \item{profile_sigma}{If \code{TRUE}, profile out the innovation variance.}
}
\value{
\code{TRUE} if the AR(1) statistic is strictly below \code{gamma}; otherwise,
\code{FALSE}.
}
\description{
Tests whether a segment is valid under an exact Gaussian AR(1)
single-mean-change statistic implemented in R.
}
\details{
The function transforms the observations into AR(1) innovations and scans
change locations that leave at least two observations on each side. The
statistic is one half of the Gaussian likelihood-ratio statistic, using the
same convention as \code{\link{AR1_single_change}}. If \code{rho} is
\code{NA}, the robust estimator requires at least three observations. The
function can be passed to \code{\link{svp0}} directly when \code{rho} is fixed,
or through a closure when other optional arguments are needed. The caller
should enforce any minimum segment length required by the chosen SVP procedure.
}
\seealso{
\code{\link{AR1_single_change}}, \code{\link{svp0}}
}
