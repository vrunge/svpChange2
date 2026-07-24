// [[Rcpp::plugins(cpp20)]]

#include <Rcpp.h>
#include <chrono>
#include <memory>
#include <string>
#include "../../src/tests.h"

using namespace Rcpp;

// Benchmark a complete streaming validity update. When evaluate_statistic is
// true, every new observation is followed by statistic(), matching SVP.
// [[Rcpp::export]]
List benchmark_robust_validity_stream(const NumericVector& data,
                                      const std::string& method,
                                      int repetitions = 5,
                                      bool evaluate_statistic = true) {
  if (repetitions < 1) stop("repetitions must be positive");
  double elapsed = 0.0;
  volatile double sink = 0.0;

  for (int rep = 0; rep < repetitions; ++rep) {
    std::unique_ptr<TestBase> test;
    if (method == "Wilcoxon") {
      test = std::make_unique<WilcoxonCost>();
    } else if (method == "MedianMood") {
      test = std::make_unique<MedianMoodCost>();
    } else {
      stop("method must be 'Wilcoxon' or 'MedianMood'");
    }

    const auto start = std::chrono::steady_clock::now();
    for (double value : data) {
      test->update(value);
      if (evaluate_statistic) sink = sink + test->statistic();
    }
    const auto end = std::chrono::steady_clock::now();
    elapsed += std::chrono::duration<double>(end - start).count();
  }

  return List::create(
    _["seconds"] = elapsed / static_cast<double>(repetitions),
    _["sink"] = static_cast<double>(sink)
  );
}
