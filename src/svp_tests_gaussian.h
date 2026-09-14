#ifndef SVPCHANGE2_SVP_TESTS_GAUSSIAN_H
#define SVPCHANGE2_SVP_TESTS_GAUSSIAN_H

#include "gaussian_focus_state.h"
#include "svp_test_base.h"
#include "svp_tests_gamma.h"
class GaussianMean : public TestBase
{
  public:
    changepoint::GaussianFocusState info;

    GaussianMean() = default;

    void update(double y) override
    {
      info.update(y);
    }

    double statistic() const override
    {
      // The unified Gaussian cost is twice the half-log-likelihood scale used
      // by SVP and by the former FOCuS implementation.
      return 0.5 * info.max_cost();
    }

    bool passes(double gamma) const override
    {
      // GaussianMean::statistic() is half the unified Gaussian cost. The
      // threshold can therefore be checked directly while scanning costs.
      return !info.exceeds_cost(2.0 * gamma);
    }

};

// Gaussian mean validity is invariant under reversal of the observations, so
// the bidirectional whole-series preflight can certify the one-segment null
// case without changing the dynamic-programming result.
template <>
struct svp_supports_both_preflight<GaussianMean> : std::true_type {};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class GaussianVariance : public TestBase
{
  public:
    GammaRate gamma_rate;

    explicit GaussianVariance(double shape = 1.0)
      : gamma_rate(shape) {}

    void update(double y) override
    {
      gamma_rate.update(y * y); // Use squared data
    }

    double statistic() const override
    {
      return gamma_rate.statistic();
    }
};

////////////////////////////////////////////////////////////////////////////////

#endif
