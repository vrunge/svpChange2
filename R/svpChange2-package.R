

#' svpChange2: Optimized Smallest Valid Partitioning for multiple
#' change-point detection
#'
#' @description
#' Tools for multiple change-point detection and segmentation of
#' **univariate** sequences using the *smallest valid partitioning* (SVP)
#' algorithm.
#'
#' SVP builds a segmentation by enforcing a **local notion of segment
#' validity**: a candidate segment is retained only if it passes a
#' user-specified validity test (for instance, a single change-point test).
#' If a change is detected, the segment is deemed **invalid** and must be
#' split, preventing invalid segments from appearing in the final partition.
#'
#' The computational core is implemented in C++ and exposed to R via **Rcpp**
#' for performance.
#'
#' @details
#' SVP is designed for situations where you want **explicit, test-based
#' control** of within-segment homogeneity rather than relying solely on a
#' global penalized objective.
#'
#'
#' Validity-test helpers include the Gaussian FOCuS likelihood-ratio test, the
#' exact AR(1) mean-change test, robust scale and range tests, and the SMUCE
#' multiscale constraint. The original methods are described in the references
#' below.
#'
#' @references
#' Romano, G., Eckley, I. A., Fearnhead, P., and Rigaill, G. (2023).
#' Fast Online Changepoint Detection via Functional Pruning CUSUM Statistics.
#' \emph{Journal of Machine Learning Research}, 24(81), 1--36.
#' \url{https://www.jmlr.org/papers/v24/21-1230.html}
#'
#' Frick, K., Munk, A., and Sieling, H. (2014). Multiscale Change Point
#' Inference. \emph{Journal of the Royal Statistical Society: Series B},
#' 76(3), 495--580. \url{https://doi.org/10.1111/rssb.12047}
#'
#' Chakar, S., Lebarbier, E., Levy-Leduc, C. and Robin, S. (2017).
#' A robust approach for estimating change-points in the mean of an AR(1)
#' process. \emph{Bernoulli}, 23(2), 1408--1447.
#'
#' @section Functions:
#' - \code{\link{svp0}}:
#' - More functions here.
#'
#' @docType package
#' @name svpChange2
#' @keywords internal
"_PACKAGE"
