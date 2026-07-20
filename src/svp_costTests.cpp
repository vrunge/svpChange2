#include "tests.h"
#include <Rcpp.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <utility>

using namespace Rcpp;

namespace {

double median_in_place(std::vector<double>& values)
{
  const size_t n = values.size();
  const size_t middle = n / 2;
  std::nth_element(values.begin(), values.begin() + middle, values.end());
  const double upper = values[middle];
  if (n % 2 == 1) return upper;
  const double lower = *std::max_element(values.begin(), values.begin() + middle);
  return 0.5 * (lower + upper);
}

double robust_ar1_rho(const std::vector<double>& data)
{
  if (data.size() < 3) stop("AR(1) estimation requires at least 3 observations");

  std::vector<double> lag1(data.size() - 1);
  std::vector<double> lag2(data.size() - 2);
  for (size_t i = 0; i < lag1.size(); ++i) {
    lag1[i] = std::fabs(data[i + 1] - data[i]);
  }
  for (size_t i = 0; i < lag2.size(); ++i) {
    lag2[i] = std::fabs(data[i + 2] - data[i]);
  }

  const double med1 = median_in_place(lag1);
  const double med2 = median_in_place(lag2);
  if (!std::isfinite(med1) || !std::isfinite(med2) || med1 <= 0.0) {
    stop("Cannot estimate AR(1) correlation from constant or non-finite data");
  }

  const double estimate = (med2 * med2) / (med1 * med1) - 1.0;
  // The model in the paper assumes stationarity. Clamping also prevents
  // numerical instability for short or highly contaminated series.
  return std::max(-0.999, std::min(0.999, estimate));
}

} // namespace

// [[Rcpp::export]]
double AR1_rho(std::vector<double> data)
{
  return robust_ar1_rho(data);
}

// [[Rcpp::export]]
List AR1_single_change(std::vector<double> data,
                       double gamma,
                       double rho = NA_REAL,
                       double sigma2 = 1.0,
                       bool profile_sigma = false)
{
  const double rho_used = NumericVector::is_na(rho) ? robust_ar1_rho(data) : rho;
  if (!std::isfinite(rho_used) || std::fabs(rho_used) >= 1.0) {
    stop("rho must be finite and strictly between -1 and 1");
  }
  if (!std::isfinite(sigma2) || sigma2 <= 0.0) {
    stop("sigma2 must be finite and positive");
  }

  AR1ExactMeanChange test(rho_used, sigma2, profile_sigma);
  for (double value : data) test.update(value);
  const double statistic = test.statistic();
  const int changepoint = test.changepoint();

  return List::create(
    _["rho"] = rho_used,
    _["sigma2"] = sigma2,
    _["profile_sigma"] = profile_sigma,
    _["statistic"] = statistic,
    _["rss0"] = test.rss0(),
    _["rss1"] = test.rss1(),
    _["changepoint"] = changepoint > 0 ? changepoint : NA_INTEGER,
    _["valid"] = statistic < gamma
  );
}
