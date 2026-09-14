#ifndef SVPCHANGE2_SVP_INTERNAL_H
#define SVPCHANGE2_SVP_INTERNAL_H

#include "validity_tests.h"
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

// Private implementation of the SVP dynamic program. The templates are kept
// in this header so that svp.cpp can instantiate svp_impl() for each validity
// test used by the public dispatcher.
namespace svp_detail {

// Return the within-segment sum of squared errors for data[s + 1:t]. The
// boundaries use the SVP convention: s is zero-based and t is one-based.
inline double segment_cost(const std::vector<double>& S1,
                           const std::vector<double>& S2,
                           std::size_t s,
                           std::size_t t)
{
  const double sum = S1[t] - S1[s];
  const double length = static_cast<double>(t - s);
  const double cost = (S2[t] - S2[s]) - sum * sum / length;

  if (cost <= 0.0) return 0.0;
  return cost;
}

// State associated with one possible previous boundary in the dynamic
// program. Each candidate owns an incremental validity-test object for its
// segment.
template <typename Test>
struct CandidateState
{
  std::size_t boundary;    // zero-based start boundary of the candidate segment
  Test test;               // incremental validity-test state
  std::size_t last_update; // last endpoint already processed by test
  bool valid;               // validity at the current endpoint

  template <typename... Args>
  CandidateState(std::size_t boundary_, Args&&... args)
    : boundary(boundary_),
      test(std::forward<Args>(args)...),
      last_update(boundary_),
      valid(true)
  {}

  CandidateState(const CandidateState&) = delete;
  CandidateState& operator=(const CandidateState&) = delete;
  CandidateState(CandidateState&&) = default;
  CandidateState& operator=(CandidateState&&) = default;
};

// Extend a candidate to endpoint t by feeding it only observations that have
// not been processed yet, then evaluate its validity statistic at t. Singleton
// segments are always valid, as in svp0().
template <typename Test>
void update_to_endpoint(CandidateState<Test>& candidate,
                        const std::vector<double>& data,
                        std::size_t t,
                        double gamma)
{
  // All built-in SVP tests are evaluated at the current endpoint. A previous
  // failure does not permanently invalidate a candidate, which is required
  // when subtests == "none".
  candidate.valid = true;

  for (std::size_t u = candidate.last_update + 1; u <= t; ++u) {
    candidate.test.update(data[u - 1]);

    // As in svp0(), a one-observation segment is always valid.
    const bool check_statistic =
      u - candidate.boundary > 1 &&
      u == t;
    if (check_statistic) {
      const double statistic = candidate.test.statistic();
      if (!std::isfinite(statistic) || statistic >= gamma) {
        candidate.valid = false;
        candidate.last_update = t;
        return;
      }
    }
  }
  candidate.last_update = t;
}

// Compare two valid partitions lexicographically: fewer segments are
// preferred first, and SSE is used only when the segment counts are equal.
inline bool is_better_candidate(std::size_t candidate_K,
                                double candidate_Q,
                                std::size_t best_K,
                                double best_Q)
{
  return candidate_K < best_K ||
    (candidate_K == best_K && candidate_Q < best_Q);
}

// Convert the internal dynamic-programming vectors into the three-column R
// matrix: cumulative cost Q, segment count K, and previous boundary.
inline Rcpp::NumericMatrix build_R_matrix(
    const std::vector<double>& Q,
    const std::vector<std::size_t>& K,
    const std::vector<std::size_t>& previous)
{
  const std::size_t n = Q.size() - 1;
  Rcpp::NumericMatrix R(n, 3);
  auto output = R.begin();

  // R matrices are stored column-major. Fill each complete column so that
  // writes are contiguous in memory; index zero is the DP initialization row.
  std::copy(Q.begin() + 1, Q.end(), output);
  std::transform(
    K.begin() + 1, K.end(), output + n,
    [](std::size_t value) {
      return static_cast<double>(value);
    }
  );
  std::transform(
    previous.begin() + 1, previous.end(), output + 2 * n,
    [](std::size_t value) {
      return static_cast<double>(value);
    }
  );

  return R;
}

// Run SVP for one concrete validity-test type. Test-specific constructor
// arguments are forwarded through Args (for example quantile, rho, or sigma2).
// The public SVP() function validates arguments and selects the Test type before
// calling this implementation.
template <typename Test, typename... Args>
Rcpp::List svp_impl(const std::vector<double>& data,
                    double gamma,
                    const std::string& subtests,
                    Args&&... args)
{
  const std::size_t n = data.size();

  std::vector<double> Q(n + 1, std::numeric_limits<double>::infinity());
  std::vector<std::size_t> K(
    n + 1, std::numeric_limits<std::size_t>::max()
  );
  std::vector<std::size_t> previous(n + 1, 0);
  Q[0] = 0.0;
  K[0] = 0;

  std::vector<std::size_t> nb(n);

  // Center only the cost calculation. This does not alter the validity test
  // or the data passed to it, but avoids cancellation in sum(y^2) - sum(y)^2/n.
  double center = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    center += (data[i] - center) / static_cast<double>(i + 1);
  }
  if (!std::isfinite(center)) {
    Rcpp::stop("Unable to compute a stable center for the segment costs.");
  }

  std::vector<double> S1(n + 1, 0.0);
  std::vector<double> S2(n + 1, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double centered = data[i] - center;
    if (!std::isfinite(centered)) {
      Rcpp::stop("The data range is too large for stable segment costs.");
    }
    S1[i + 1] = S1[i] + centered;
    S2[i + 1] = S2[i] + centered * centered;
    if (!std::isfinite(S1[i + 1]) || !std::isfinite(S2[i + 1])) {
      Rcpp::stop("The data range is too large for stable segment costs.");
    }
  }

  using Candidate = CandidateState<Test>;
  std::vector<Candidate> candidates;
  std::vector<Candidate> retained;
  candidates.reserve(n + 1);
  retained.reserve(n + 1);
  candidates.emplace_back(0, std::forward<Args>(args)...);

  // Dynamic-programming loop over endpoints t = 1, ..., n.
  for (std::size_t t = 1; t <= n; ++t) {
    nb[t - 1] = candidates.size();

    double best_Q = std::numeric_limits<double>::infinity();
    std::size_t best_K = std::numeric_limits<std::size_t>::max();
    std::size_t best_s = 0;
    std::size_t max_invalid_boundary = 0;

    for (Candidate& candidate : candidates) {
      update_to_endpoint(candidate, data, t, gamma);

      if (!candidate.valid) {
        max_invalid_boundary = std::max(
          max_invalid_boundary, candidate.boundary
        );
        continue;
      }

      if (K[candidate.boundary] == std::numeric_limits<std::size_t>::max()) {
        continue;
      }

      const std::size_t candidate_K = K[candidate.boundary] + 1;
      const double candidate_Q = Q[candidate.boundary] +
        segment_cost(S1, S2, candidate.boundary, t);

      if (std::isfinite(candidate_Q) &&
          is_better_candidate(candidate_K, candidate_Q,
                              best_K, best_Q)) {
        best_Q = candidate_Q;
        best_K = candidate_K;
        best_s = candidate.boundary;
      }
    }

    if (!std::isfinite(best_Q)) {
      Rcpp::stop("No valid finite-cost partition was found at the current endpoint.");
    }

    Q[t] = best_Q;
    K[t] = best_K;
    previous[t] = best_s;

    if (subtests != "none") {
      retained.clear();
      for (Candidate& candidate : candidates) {
        bool keep = true;

        if ((subtests == "right" || subtests == "both") &&
            !candidate.valid) {
          keep = false;
        }
        if ((subtests == "left" || subtests == "both") &&
            candidate.boundary < max_invalid_boundary) {
          keep = false;
        }

        if (keep) retained.push_back(std::move(candidate));
      }
      candidates.swap(retained);
    }

    candidates.emplace_back(t, std::forward<Args>(args)...);
  }

  std::vector<std::size_t> changepoints;
  changepoints.reserve(n / 2 + 1);
  for (std::size_t i = n; i > 0;) {
    if (previous[i] >= i) {
      Rcpp::stop("Internal error while reconstructing the SVP partition.");
    }
    changepoints.push_back(i);
    i = previous[i];
  }
  std::reverse(changepoints.begin(), changepoints.end());

  std::vector<std::size_t> last_index_set;
  last_index_set.reserve(candidates.size());
  for (auto it = candidates.rbegin(); it != candidates.rend(); ++it) {
    last_index_set.push_back(it->boundary);
  }

  return Rcpp::List::create(
    Rcpp::Named("changepoints") = changepoints,
    Rcpp::Named("lastIndexSet") = last_index_set,
    Rcpp::Named("nb") = nb,
    Rcpp::Named("costQ") = R_NilValue,
    Rcpp::Named("R") = build_R_matrix(Q, K, previous)
  );
}

} // namespace svp_detail

#endif
