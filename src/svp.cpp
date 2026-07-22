#include "tests.h"
#include <Rcpp.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>
#include <utility>

using namespace Rcpp;

namespace {

inline double segment_cost(const std::vector<double>& S1,
                           const std::vector<double>& S2,
                           size_t s,
                           size_t t)
{
  const double sum = S1[t] - S1[s];
  return (S2[t] - S2[s]) - sum * sum / static_cast<double>(t - s);
}

double median_in_place(std::vector<double>& values)
{
  const size_t middle = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + middle, values.end());
  const double upper = values[middle];
  if (values.size() % 2 == 1) return upper;
  const double lower = *std::max_element(values.begin(), values.begin() + middle);
  return 0.5 * (lower + upper);
}

double robust_ar1_rho(const std::vector<double>& data)
{
  if (data.size() < 3) stop("AR(1) estimation requires at least 3 observations");
  std::vector<double> lag1(data.size() - 1), lag2(data.size() - 2);
  for (size_t i = 0; i < lag1.size(); ++i) lag1[i] = std::fabs(data[i + 1] - data[i]);
  for (size_t i = 0; i < lag2.size(); ++i) lag2[i] = std::fabs(data[i + 2] - data[i]);
  const double med1 = median_in_place(lag1);
  const double med2 = median_in_place(lag2);
  if (!std::isfinite(med1) || !std::isfinite(med2) || med1 <= 0.0)
    stop("Cannot estimate AR(1) correlation from constant or non-finite data");
  const double estimate = (med2 * med2) / (med1 * med1) - 1.0;
  return std::max(-0.999, std::min(0.999, estimate));
}

NumericMatrix build_R_matrix(const std::vector<double>& Q,
                             const std::vector<size_t>& K,
                             const std::vector<size_t>& previous)
{
  const size_t n = Q.size() - 1;
  NumericMatrix R(n, 3);
  for (size_t t = 1; t <= n; ++t) {
    const size_t row = t - 1;
    R(row, 0) = Q[t];
    R(row, 1) = static_cast<double>(K[t]);
    R(row, 2) = static_cast<double>(previous[t]);
  }
  return R;
}

inline double validity_threshold(size_t n,
                                 size_t segment_length,
                                 double gamma,
                                 bool use_multiscale_gamma,
                                 double sigma2,
                                 double q_alpha_n)
{
  if (!use_multiscale_gamma) return gamma;

  return sigma2 * (
    q_alpha_n +
    2.0 * std::log(std::exp(1.0) * static_cast<double>(n) /
                   static_cast<double>(segment_length))
  );
}

template <typename Test, typename... Args>
List svp_impl(const std::vector<double>& data,
              double gamma,
              bool prune_after_if_unvalid,
              bool prune_before_if_invalid,
              bool use_multiscale_gamma,
              double sigma2,
              double q_alpha_n,
              Args&&... args)
{
  const size_t n = data.size();

  std::vector<double> Q(n + 1, std::numeric_limits<double>::infinity());
  std::vector<size_t> K(n + 1, std::numeric_limits<size_t>::max());
  std::vector<size_t> previous(n + 1, 0);
  Q[0] = 0.0;
  K[0] = 0;

  std::vector<size_t> nb(n);

  std::vector<double> S1(n + 1, 0.0);
  std::vector<double> S2(n + 1, 0.0);
  for (size_t i = 0; i < n; ++i) {
    S1[i + 1] = S1[i] + data[i];
    S2[i + 1] = S2[i] + data[i] * data[i];
  }

  std::vector<size_t> index;
  std::vector<Test> tests;
  std::vector<size_t> last_updates;
  index.reserve(n + 1);
  tests.reserve(n + 1);
  last_updates.reserve(n + 1);
  index.push_back(0);
  tests.emplace_back(std::forward<Args>(args)...);
  last_updates.push_back(0);

  for (size_t t = 1; t <= n; ++t) {
    const size_t m = index.size();
    nb[t - 1] = m;

    double best_Q = std::numeric_limits<double>::infinity();
    size_t best_K = std::numeric_limits<size_t>::max();
    size_t best_s = 0;

    if (!prune_before_if_invalid && prune_after_if_unvalid) {
      size_t write = 0;

      for (size_t k = 0; k < m; ++k) {
        const size_t s = index[k];
        const size_t candidate_K = K[s] + 1;

        if (candidate_K > best_K) {
          for (size_t j = k; j < m; ++j) {
            size_t last_up = last_updates[j];
            for (size_t u = last_up + 1; u <= t; ++u) {
              tests[j].update(data[u - 1]);
            }
            last_up = t;

            const double threshold = validity_threshold(
              n,
              t - index[j],
              gamma,
              use_multiscale_gamma,
              sigma2,
              q_alpha_n
            );

            if (tests[j].statistic() < threshold) {
              if (write != j) {
                index[write] = index[j];
                tests[write] = std::move(tests[j]);
              }
              last_updates[write] = last_up;
              ++write;
            }
          }
          break;
        }

        size_t last_up = last_updates[k];
        for (size_t u = last_up + 1; u <= t; ++u) {
          tests[k].update(data[u - 1]);
        }
        last_up = t;

        const double threshold = validity_threshold(
          n,
          t - s,
          gamma,
          use_multiscale_gamma,
          sigma2,
          q_alpha_n
        );

        if (tests[k].statistic() < threshold) {
          const double candidate_Q = Q[s] + segment_cost(S1, S2, s, t);

          if (candidate_K < best_K ||
              (candidate_K == best_K && candidate_Q < best_Q)) {
            best_Q = candidate_Q;
            best_K = candidate_K;
            best_s = s;
          }

          if (write != k) {
            index[write] = s;
            tests[write] = std::move(tests[k]);
          }
          last_updates[write] = last_up;
          ++write;
        }
      }

      index.resize(write);
      tests.erase(tests.begin() + write, tests.end());
      last_updates.resize(write);
    } else if (!prune_before_if_invalid && !prune_after_if_unvalid) {
      for (size_t k = 0; k < m; ++k) {
        const size_t s = index[k];

        for (size_t u = last_updates[k] + 1; u <= t; ++u) {
          tests[k].update(data[u - 1]);
        }
        last_updates[k] = t;

        const double threshold = validity_threshold(
          n,
          t - s,
          gamma,
          use_multiscale_gamma,
          sigma2,
          q_alpha_n
        );

        if (tests[k].statistic() < threshold) {
          const size_t candidate_K = K[s] + 1;
          const double candidate_Q = Q[s] + segment_cost(S1, S2, s, t);

          if (candidate_K < best_K ||
              (candidate_K == best_K && candidate_Q < best_Q)) {
            best_Q = candidate_Q;
            best_K = candidate_K;
            best_s = s;
          }
        }
      }
    } else {
      size_t write = 0;

      for (size_t k = 0; k < m; ++k) {
        const size_t s = index[k];

        size_t last_up = last_updates[k];
        for (size_t u = last_up + 1; u <= t; ++u) {
          tests[k].update(data[u - 1]);
        }
        last_up = t;

        const double threshold = validity_threshold(
          n,
          t - s,
          gamma,
          use_multiscale_gamma,
          sigma2,
          q_alpha_n
        );
        const bool valid = tests[k].statistic() < threshold;

        if (!valid) {
          write = 0;
        }

        if (valid) {
          const size_t candidate_K = K[s] + 1;
          // Equation (4) optimization: once a valid candidate requires
          // strictly more segments than the current lexicographic optimum,
          // its cost cannot affect the result. We still retain/update it
          // below when pruning requires the candidate state.
          if (candidate_K <= best_K) {
            const double candidate_Q = Q[s] + segment_cost(S1, S2, s, t);

            if (candidate_K < best_K ||
                (candidate_K == best_K && candidate_Q < best_Q)) {
              best_Q = candidate_Q;
              best_K = candidate_K;
              best_s = s;
            }
          }
        }

        if (valid || !prune_after_if_unvalid) {
          if (write != k) {
            index[write] = s;
            tests[write] = std::move(tests[k]);
          }
          last_updates[write] = last_up;
          ++write;
        }
      }

      index.resize(write);
      tests.erase(tests.begin() + write, tests.end());
      last_updates.resize(write);
    }

    Q[t] = best_Q;
    K[t] = best_K;
    previous[t] = best_s;

    index.push_back(t);
    tests.emplace_back(std::forward<Args>(args)...);
    last_updates.push_back(t);
  }

  std::vector<size_t> changepoints;
  changepoints.reserve(n / 2 + 1);
  for (size_t i = n; i > 0; i = previous[i]) {
    changepoints.push_back(i);
  }
  std::reverse(changepoints.begin(), changepoints.end());
  std::reverse(index.begin(), index.end());

  return List::create(
    _["changepoints"] = changepoints,
    _["lastIndexSet"] = index,
    _["nb"] = nb,
    _["costQ"] = NULL,
    _["R"] = build_R_matrix(Q, K, previous)
  );
}

template <typename Test, typename... Args>
List svp_cost_impl(const std::vector<double>& data,
                   double gamma,
                   bool prune_after_if_unvalid,
                   bool prune_before_if_invalid,
                   bool use_multiscale_gamma,
                   double sigma2,
                   double q_alpha_n,
                   Args&&... args)
{
  return svp_impl<Test>(data, gamma, prune_after_if_unvalid,
                        prune_before_if_invalid, use_multiscale_gamma,
                        sigma2, q_alpha_n, std::forward<Args>(args)...);
}

} // namespace

// Internal backend for the R validity wrappers.  Checking every prefix matches
// the online pruning rule used by SVP; checking only the final statistic gives
// the legacy valid_FOCUS_last behaviour.
// [[Rcpp::export(name = ".focus_valid_cpp")]]
bool focus_valid_cpp(std::vector<double> data,
                     double gamma,
                     bool check_all_prefixes = true)
{
  GaussianMean test;
  for (double value : data) {
    test.update(value);
    if (check_all_prefixes && test.statistic() >= gamma) return false;
  }
  return test.statistic() < gamma;
}

//' Smallest Valid Partitioning with Incremental Validity Tests
//'
//' Segments a univariate series into the smallest number of segments that pass
//' a selected validity test. Among partitions with the same number of segments,
//' the function chooses the one with the smallest within-segment squared-error
//' cost. The validity statistics are maintained incrementally in C++, so this
//' interface is considerably faster than supplying an R validity function to
//' [svp0()].
//'
//' A candidate segment is valid while its test statistic is strictly below its
//' threshold. The available values of `test` are:
//'
//' * `"gaussian_mean"`: Gaussian FOCUS likelihood-ratio test for a change in
//'   mean. Use this for independent Gaussian observations with constant
//'   variance.
//' * `"gamma_rate"`: likelihood-ratio test for a change in the rate of
//'   positive Gamma observations (unit shape in the current implementation).
//' * `"gaussian_variance"`: Gamma-rate test applied to squared observations,
//'   for changes in Gaussian variance around a known zero mean.
//' * `"AR1"`: exact fixed-`rho` Gaussian likelihood-ratio scan for a change in
//'   the marginal mean of an AR(1) series. `sigma2` is the innovation variance.
//' * `"AR1Profile"`: the same exact AR(1) scan, profiling out the innovation
//'   variance. This is useful when its scale is unknown.
//' * `"AR1Focus"`: faster approximate AR(1) test that applies Gaussian FOCUS
//'   to the innovations `x[t] - rho * x[t - 1]`. The exact `"AR1"` test is the
//'   preferred choice when boundary accuracy matters.
//'
//' The main pruning setting is `prune_after_if_unvalid = TRUE`: once a segment
//' beginning at a candidate boundary becomes invalid, that candidate is not
//' extended further. Setting `prune_before_if_invalid = TRUE` additionally
//' removes all older candidate boundaries when a later candidate fails. Thus
//' `TRUE/TRUE` is the most aggressive FOCUS pruning configuration. Set either
//' option to `FALSE` only when comparing pruning rules; doing so can retain more
//' candidates and increase run time.
//'
//' With `use_multiscale_gamma = FALSE`, every segment uses the scalar `gamma`.
//' With it set to `TRUE`, a candidate of length `l` in a series of length `n`
//' instead uses
//' `sigma2 * (q_alpha_n + 2 * log(exp(1) * n / l))`.
//' The user must supply a calibration `q_alpha_n` appropriate to the desired
//' global error level and sample size; it is not estimated by this function.
//'
//' @param data Numeric vector containing the univariate series. Missing or
//'   non-finite values are not supported.
//' @param gamma Positive scalar validity threshold used when
//'   `use_multiscale_gamma = FALSE`. A larger value accepts longer or less
//'   homogeneous segments and therefore generally produces fewer changes.
//' @param test Character scalar selecting one of the validity tests listed in
//'   Details.
//' @param prune_after_if_unvalid Logical; discard a candidate boundary after
//'   its current segment fails the test.
//' @param prune_before_if_invalid Logical; when a candidate segment fails,
//'   also discard candidate boundaries older than its start.
//' @param use_multiscale_gamma Logical; use a segment-length-dependent threshold
//'   instead of `gamma`.
//' @param sigma2 Positive finite variance. It scales the multiscale threshold
//'   and is the innovation variance for `test = "AR1"`.
//' @param q_alpha_n Numeric global calibration constant for the multiscale
//'   threshold. It must be calibrated externally under the relevant null model.
//' @param rho AR(1) coefficient for the three AR1 tests. It must be finite and
//'   strictly between -1 and 1. Use [AR1_rho()] to obtain a robust estimate if
//'   `rho` is unknown, it is estimated robustly from the full series.
//' @param profile_sigma Logical; profile the innovation variance when
//'   `test = "AR1"`. Using `test = "AR1Profile"` has the same effect.
//' @param quantile Quantile level used by `"quantile"` and
//'   `"quantileExact"`. It is ignored by other tests.
//'
//' @return A list with `changepoints` (the inclusive end of every segment,
//'   including `length(data)`), `lastIndexSet` (candidate boundaries remaining
//'   at termination), `nb` (candidate count at each time), `costQ` (currently
//'   `NULL`), and `R`. Row `t + 1` of matrix `R` stores the best cumulative
//'   squared-error cost, number of segments, and previous boundary at time `t`.
//'
//' @examples
//' set.seed(1)
//' x <- rep(c(0, 2, -1), each = 40) + rnorm(120)
//' SVP(x, gamma = 1.5 * log(length(x)), test = "gaussian_mean")$changepoints
//'
//' # Aggressive pruning and a multiscale threshold:
//' SVP(x, gamma = 0, test = "gaussian_mean",
//'     prune_after_if_unvalid = TRUE,
//'     prune_before_if_invalid = TRUE,
//'     use_multiscale_gamma = TRUE, sigma2 = 1, q_alpha_n = 3)
//'
//' @seealso [svp0()] for arbitrary R validity functions, [AR1_rho()], and
//'   [AR1_single_change()].
//' @export
// [[Rcpp::export]]
List SVP(std::vector<double> data,
         double gamma,
         std::string test,
         bool prune_after_if_unvalid = true,
         bool prune_before_if_invalid = false,
         bool use_multiscale_gamma = false,
         double sigma2 = 1.0,
         double q_alpha_n = 0.0,
         double rho = NA_REAL,
         bool profile_sigma = false,
         double quantile = 0.01)
{
  if (!std::isfinite(sigma2) || sigma2 <= 0.0) {
    stop("sigma2 must be finite and positive");
  }
  if ((test == "quantile" || test == "quantileExact") &&
      (!std::isfinite(quantile) || quantile <= 0.0 || quantile >= 1.0)) {
    stop("quantile must be finite and strictly between 0 and 1");
  }

  if (test == "gaussian_mean") {
    return svp_impl<GaussianMean>(
      data,
      gamma,
      prune_after_if_unvalid,
      prune_before_if_invalid,
      use_multiscale_gamma,
      sigma2,
      q_alpha_n
    );
  }
  if (test == "gamma_rate") {
    return svp_impl<GammaRate>(
      data,
      gamma,
      prune_after_if_unvalid,
      prune_before_if_invalid,
      use_multiscale_gamma,
      sigma2,
      q_alpha_n
    );
  }
  if (test == "gaussian_variance") {
    return svp_impl<GaussianVariance>(
      data,
      gamma,
      prune_after_if_unvalid,
      prune_before_if_invalid,
      use_multiscale_gamma,
      sigma2,
      q_alpha_n
    );
  }
  if (test == "quantileExact") {
    return svp_cost_impl<QuantileCostExact>(
      data, gamma, prune_after_if_unvalid, prune_before_if_invalid,
      use_multiscale_gamma, sigma2, q_alpha_n, quantile
    );
  }
  if (test == "quantile") {
    return svp_cost_impl<QuantileCost>(
      data, gamma, prune_after_if_unvalid, prune_before_if_invalid,
      use_multiscale_gamma, sigma2, q_alpha_n, quantile
    );
  }
  if (test == "varCost") {
    return svp_cost_impl<varCost>(
      data, gamma, prune_after_if_unvalid, prune_before_if_invalid,
      use_multiscale_gamma, sigma2, q_alpha_n
    );
  }
  if (test == "WilcoxonCost") {
    return svp_cost_impl<WilcoxonCost>(
      data, gamma, prune_after_if_unvalid, prune_before_if_invalid,
      use_multiscale_gamma, sigma2, q_alpha_n
    );
  }
  if (test == "MedianMoodCost") {
    return svp_cost_impl<MedianMoodCost>(
      data, gamma, prune_after_if_unvalid, prune_before_if_invalid,
      use_multiscale_gamma, sigma2, q_alpha_n
    );
  }
  if (test == "AR1" || test == "AR1Profile" || test == "AR1Focus") {
    const double rho_used = NumericVector::is_na(rho) ? robust_ar1_rho(data) : rho;
    if (!std::isfinite(rho_used) || std::fabs(rho_used) >= 1.0)
      stop("rho must be finite and strictly between -1 and 1");
    if (test == "AR1Focus") {
      List result = svp_impl<AR1MeanChange>(
        data,
        gamma,
        prune_after_if_unvalid,
        prune_before_if_invalid,
        use_multiscale_gamma,
        sigma2,
        q_alpha_n,
        rho_used
      );
      result["rho"] = rho_used;
      result["sigma2"] = sigma2;
      return result;
    }
    List result = svp_impl<AR1ExactMeanChange>(
      data,
      gamma,
      prune_after_if_unvalid,
      prune_before_if_invalid,
      use_multiscale_gamma,
      sigma2,
      q_alpha_n,
        rho_used,
      sigma2,
      test == "AR1Profile" || profile_sigma
    );
    result["rho"] = rho_used;
    result["sigma2"] = sigma2;
    return result;
  }

  stop("Unknown test type");
}
