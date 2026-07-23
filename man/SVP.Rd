% Generated from the documentation in src/svp.cpp
\name{SVP}
\alias{SVP}
\title{Smallest Valid Partitioning with Incremental Validity Tests}
\usage{
SVP(
  data,
  gamma,
  test = "gaussian_mean",
  prune_after_if_unvalid = TRUE,
  prune_before_if_invalid = FALSE,
  use_multiscale_gamma = FALSE,
  sigma2 = 1,
  q_alpha_n = 0,
  rho = NA_real_,
  profile_sigma = FALSE,
  quantile = 0.01
)
}
\arguments{
\item{data}{Numeric vector containing the univariate series. Missing and
non-finite values are not supported.}

\item{gamma}{Scalar validity threshold when \code{use_multiscale_gamma =
FALSE}. Larger values accept less homogeneous segments and generally produce
fewer changes.}

\item{test}{Character scalar selecting a validity test; see Details. Defaults
to \code{"gaussian_mean"}.}

\item{prune_after_if_unvalid}{If \code{TRUE}, discard a candidate boundary
after its current segment fails the test.}

\item{prune_before_if_invalid}{If \code{TRUE}, when a candidate segment fails,
also discard candidate boundaries older than its start.}

\item{use_multiscale_gamma}{If \code{TRUE}, use the length-dependent threshold
described in Details instead of \code{gamma}.}

\item{sigma2}{Positive finite variance. It scales the multiscale threshold and
is the innovation variance for \code{test = "AR1"}.}

\item{q_alpha_n}{Global multiscale calibration constant. This must be
calibrated externally for the desired error level and sample size.}

\item{rho}{AR(1) coefficient for the AR1 tests. It must be finite and strictly
  between -1 and 1. If omitted, it is estimated robustly from the full series.}

\item{profile_sigma}{If \code{TRUE}, profile the innovation variance for
\code{test = "AR1"}; \code{test = "AR1Profile"} has the same effect.}

\item{quantile}{Quantile level used by \code{"quantile"} and
\code{"quantileExact"}; ignored by other tests.}
}
\value{
A list containing:
\describe{
  \item{changepoints}{Inclusive end of every segment, including
  \code{length(data)}.}
  \item{lastIndexSet}{Candidate boundaries remaining at termination.}
  \item{nb}{Number of candidate boundaries considered at each time.}
  \item{costQ}{Currently \code{NULL}.}
  \item{R}{Dynamic-programming matrix. Row \code{t + 1} contains the best
  cumulative squared-error cost, number of segments, and previous boundary at
  time \code{t}.}
}
}
\description{
Segments a univariate series into the smallest number of segments that pass a
selected validity test. Among partitions with the same number of segments, it
chooses the one with the smallest within-segment squared-error cost. Validity
statistics are maintained incrementally in C++.
}
\details{
The available validity tests are:
\itemize{
  \item \code{"gaussian_mean"}: Gaussian FOCUS likelihood-ratio test for a
  change in mean. Use for independent Gaussian data with constant variance.
  \item \code{"gamma_rate"}: likelihood-ratio test for a rate change in
  positive Gamma data (unit shape in the current implementation).
  \item \code{"gaussian_variance"}: the Gamma-rate test applied to squared
  observations, for variance changes around a known zero mean.
  \item \code{"AR1"}: exact fixed-\code{rho} Gaussian likelihood-ratio scan
  for a change in the marginal AR(1) mean. \code{sigma2} is the innovation
  variance.
  \item \code{"AR1Profile"}: the same exact AR(1) scan, profiling out the
  innovation variance.
  \item \code{"AR1Focus"}: faster approximation applying Gaussian FOCUS to
  innovations \code{x[t] - rho * x[t - 1]}. Prefer \code{"AR1"} when exact
  boundary treatment matters.
  \item \code{"quantile"} and \code{"quantileExact"}: quantile-based robust
  tests, using \code{quantile} as their quantile level.
  \item \code{"varCost"}: robust variance cost test.
  \item \code{"WilcoxonCost"}: Wilcoxon rank-based cost test.
  \item \code{"MedianMoodCost"}: Median-Mood rank-based cost test.
}

The usual pruning configuration is \code{TRUE/FALSE}. The \code{TRUE/TRUE}
configuration enables both after- and before-invalid pruning and is the most
aggressive FOCUS variant. Disabling pruning retains more candidates and can
substantially increase run time.

With multiscale thresholds, a candidate segment of length \code{l} in a series
of length \code{n} uses
\deqn{sigma2 \{q_alpha_n + 2 \log(e n/l)\}.}
The function does not calibrate \code{q_alpha_n}.

The available nonparametric and robust cost tests are \code{"quantile"},
\code{"quantileExact"}, \code{"varCost"}, \code{"WilcoxonCost"}, and
\code{"MedianMoodCost"}. They are selected through the same \code{test}
argument and use the same pruning options as the Gaussian tests.
}
\examples{
set.seed(1)
x <- rep(c(0, 2, -1), each = 40) + rnorm(120)
SVP(x, gamma = 1.5 * log(length(x)))$changepoints

SVP(x, gamma = 0, test = "gaussian_mean",
    prune_after_if_unvalid = TRUE,
    prune_before_if_invalid = TRUE,
    use_multiscale_gamma = TRUE, sigma2 = 1, q_alpha_n = 3)
}
\seealso{
\code{\link{svp0}}, \code{\link{AR1_rho}}, \code{\link{AR1_single_change}}
}
