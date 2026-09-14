#ifndef SVPCHANGE2_SVP_TEST_BASE_H
#define SVPCHANGE2_SVP_TEST_BASE_H

#include <cmath>
#include <type_traits>

// Tests opt in only when a bidirectional whole-series preflight is both
// mathematically valid and useful. The generic default is deliberately
// conservative: a test can be stateful, directional, or too expensive for
// the preflight to improve the computation.
template <typename Test>
struct svp_supports_both_preflight : std::false_type {};

class TestBase
{
  public:
    virtual ~TestBase() = default;
    virtual void update(double y) = 0;
    virtual double statistic() const = 0;

    // Validity-only interface. Tests with an incremental threshold check can
    // override this to stop before computing a full statistic.
    virtual bool passes(double gamma) const
    {
      const double value = statistic();
      return std::isfinite(value) && value < gamma;
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#endif
