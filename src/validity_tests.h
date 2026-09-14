#pragma once
#include "focus.h"
#include "unified_Costs.h"
#include "unified_CostsArp.h"
#include <memory>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdint>


class TestBase
{
  public:
    virtual ~TestBase() = default;
    virtual void update(double y) = 0;
    virtual double statistic() const = 0;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class GaussianMean : public TestBase
{
  public:
    std::unique_ptr<changepoint::UnivariateInfo> info;

    GaussianMean()
      : info(std::make_unique<changepoint::UnivariateInfo>()) {}

    void update(double y) override
    {
      info->update(y);
    }

    double statistic() const override
    {
      // The unified Gaussian cost is twice the half-log-likelihood scale used
      // by SVP and by the former FOCuS implementation.
      return 0.5 * info->gaussian_max_cost();
    }

    int changepoint() const
    {
      const auto result = changepoint::compute_costs_gaussian(*info, {});
      return result.changepoint.has_value()
        ? static_cast<int>(*result.changepoint) : -1;
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Gamma change-in-rate test
class GammaRate : public TestBase
{
  public:
    std::unique_ptr<Info> info;

    GammaRate(double shape = 1.0)
    {
      auto newP = [shape](double St, int tau, double m0){
        auto p = std::make_unique<PieceGam>();
        p->St = St;
        p->tau = tau;
        p->m0 = m0;
        p->set_shape(shape);
        return p;
      };
      info = std::make_unique<Info>(newP, NAN);
    }

    void update(double y) override
    {
      info->update(y);
    }

    double statistic() const override
    {
      return info->statistic();
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class GaussianVariance : public TestBase
{
  public:
    std::unique_ptr<GammaRate> gamma_rate;

    GaussianVariance(double shape = 1.0)
    {
      gamma_rate = std::make_unique<GammaRate>(shape);
    }

    void update(double y) override
    {
      gamma_rate->update(y * y); // Use squared data
    }

    double statistic() const override
    {
      return gamma_rate->statistic();
    }
};

////////////////////////////////////////////////////////////////////////////////
// Exact fixed-rho AR(1) FOCuS likelihood recurrence.
////////////////////////////////////////////////////////////////////////////////

class AR1MeanChange : public TestBase
{
public:
  explicit AR1MeanChange(double rho, double sigma2 = 1.0)
    : sigma2_(sigma2), observation_count_(0),
      info_(std::make_unique<changepoint::ARpInfo>(
        std::vector<double>{rho}, false, 0.0)) {}

  void update(double y) override
  {
    ++observation_count_;
    info_->update(y);
  }

  double statistic() const override
  {
    if (observation_count_ < 4) return 0.0;
    const auto result = changepoint::compute_costs_arp_typed(*info_);
    if (!result.stat.has_value()) return 0.0;
    // Unified AR FOCuS returns the conventional 2 log-likelihood ratio;
    // SVP's gamma and exact AR1 test use the half-log-likelihood scale.
    return std::max(0.0, 0.5 * std::get<double>(*result.stat) / sigma2_);
  }

  int changepoint() const
  {
    const auto result = changepoint::compute_costs_arp_typed(*info_);
    return result.changepoint.has_value()
      ? static_cast<int>(*result.changepoint) : -1;
  }

private:
  double sigma2_;
  std::size_t observation_count_;
  std::unique_ptr<changepoint::ARpInfo> info_;
};

////////////////////////////////////////////////////////////////////////////////
// Exact fixed-rho likelihood scan for a change in the marginal AR(1) mean.
//
// For a change after observation tau, the transition tau -> tau + 1 has residual
//   x[tau + 1] - rho*x[tau] - mu2 + rho*mu1,
// so it cannot be represented by an ordinary two-mean split of the innovations.
// Prefix sums make each candidate tau O(1); scanning all tau is O(n).
////////////////////////////////////////////////////////////////////////////////

class AR1ExactMeanChange : public TestBase
{
public:
  AR1ExactMeanChange(double rho, double sigma2 = 1.0, bool profile_sigma = false)
    : rho_(rho),
      sigma2_(sigma2),
      profile_sigma_(profile_sigma),
      previous_(0.0),
      has_previous_(false),
      observation_count_(0),
      cached_statistic_(0.0),
      cached_rss0_(0.0),
      cached_rss1_(0.0),
      cached_changepoint_(-1),
      dirty_(false)
  {
    prefix_sum_.push_back(0.0);
    prefix_sq_.push_back(0.0);
  }

  void update(double y) override
  {
    ++observation_count_;
    if (has_previous_) {
      const double innovation = y - rho_ * previous_;
      innovations_.push_back(innovation);
      prefix_sum_.push_back(prefix_sum_.back() + innovation);
      prefix_sq_.push_back(prefix_sq_.back() + innovation * innovation);
    }
    previous_ = y;
    has_previous_ = true;
    dirty_ = true;
  }

  double statistic() const override
  {
    compute();
    return cached_statistic_;
  }

  int changepoint() const
  {
    compute();
    return cached_changepoint_;
  }

  double rss0() const
  {
    compute();
    return cached_rss0_;
  }

  double rss1() const
  {
    compute();
    return cached_rss1_;
  }

private:
  void compute() const
  {
    if (!dirty_) return;

    cached_statistic_ = 0.0;
    cached_rss0_ = 0.0;
    cached_rss1_ = 0.0;
    cached_changepoint_ = -1;

    // Candidate changes leave at least two observations on each side.
    if (observation_count_ < 4) {
      dirty_ = false;
      return;
    }

    const std::size_t n = observation_count_;
    const std::size_t transition_count = n - 1;
    const double total_sum = prefix_sum_[transition_count];
    const double total_sq = prefix_sq_[transition_count];
    cached_rss0_ = std::max(
      0.0,
      total_sq - total_sum * total_sum /
        static_cast<double>(transition_count)
    );

    const double one_minus_rho = 1.0 - rho_;
    double best_rss = std::numeric_limits<double>::infinity();
    int best_tau = -1;

    for (std::size_t tau = 2; tau <= n - 2; ++tau) {
      const std::size_t pre_count = tau - 1;
      const std::size_t post_count = n - tau - 1;
      const double pre_sum = prefix_sum_[pre_count];
      const double transition = innovations_[tau - 1];
      const double post_sum = total_sum - prefix_sum_[tau];

      const double a =
        static_cast<double>(pre_count) * one_minus_rho * one_minus_rho +
        rho_ * rho_;
      const double b = -rho_;
      const double d =
        1.0 +
        static_cast<double>(post_count) * one_minus_rho * one_minus_rho;
      const double v1 = one_minus_rho * pre_sum - rho_ * transition;
      const double v2 = transition + one_minus_rho * post_sum;
      const double determinant = a * d - b * b;

      if (determinant <= std::numeric_limits<double>::epsilon()) continue;

      const double fitted =
        (d * v1 * v1 - 2.0 * b * v1 * v2 + a * v2 * v2) /
        determinant;
      const double rss = std::max(0.0, total_sq - fitted);
      if (rss < best_rss) {
        best_rss = rss;
        best_tau = static_cast<int>(tau);
      }
    }

    if (best_tau < 0) {
      dirty_ = false;
      return;
    }

    cached_rss1_ = best_rss;
    cached_changepoint_ = best_tau;
    if (profile_sigma_) {
      if (cached_rss0_ > 0.0 && best_rss <=
          std::numeric_limits<double>::epsilon() * cached_rss0_) {
        cached_statistic_ = std::numeric_limits<double>::infinity();
      } else if (best_rss > 0.0 && cached_rss0_ > best_rss) {
        // Use the log-likelihood-ratio scale (FOCuS uses this convention).
        cached_statistic_ =
          0.5 * static_cast<double>(transition_count) *
          std::log(cached_rss0_ / best_rss);
      }
    } else {
      // Half the usual 2 log-likelihood ratio, matching Gaussian FOCuS.
      cached_statistic_ =
        std::max(0.0, (cached_rss0_ - best_rss) / (2.0 * sigma2_));
    }
    dirty_ = false;
  }

  double rho_;
  double sigma2_;
  bool profile_sigma_;
  double previous_;
  bool has_previous_;
  std::size_t observation_count_;
  std::vector<double> innovations_;
  std::vector<double> prefix_sum_;
  std::vector<double> prefix_sq_;
  mutable double cached_statistic_;
  mutable double cached_rss0_;
  mutable double cached_rss1_;
  mutable int cached_changepoint_;
  mutable bool dirty_;
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// ----------------------------------------------------------
// Order-statistics tree (treap) for doubles, with multiplicities
// ----------------------------------------------------------

class OrderStatisticTree
{
private:
  struct Node {
    double key;
    int priority;   // heap key
    int count;      // multiplicity of this key
    int size;       // total size of subtree (with multiplicities)
    Node* left;
    Node* right;

    Node(double k, int pr)
      : key(k), priority(pr), count(1), size(1),
        left(nullptr), right(nullptr) {}
  };

  Node* root_;
  std::uint64_t rng_state_;

  static int get_size(Node* n) {
    return n ? n->size : 0;
  }

  static void recalc(Node* n) {
    if (!n) return;
    n->size = n->count + get_size(n->left) + get_size(n->right);
  }

  Node* rotate_right(Node* y) {
    Node* x = y->left;
    Node* T2 = x->right;
    x->right = y;
    y->left = T2;
    recalc(y);
    recalc(x);
    return x;
  }

  Node* rotate_left(Node* x) {
    Node* y = x->right;
    Node* T2 = y->left;
    y->left = x;
    x->right = T2;
    recalc(x);
    recalc(y);
    return y;
  }

  int random_priority() {
    // xorshift64*: fast deterministic priorities without constructing a
    // hardware-seeded random engine for every SVP candidate.
    rng_state_ ^= rng_state_ >> 12;
    rng_state_ ^= rng_state_ << 25;
    rng_state_ ^= rng_state_ >> 27;
    return static_cast<int>((rng_state_ * 2685821657736338717ULL) >> 33);
  }

  Node* insert(Node* node, double key) {
    if (!node)
      return new Node(key, random_priority());

    if (key == node->key) {
      node->count += 1;
    } else if (key < node->key) {
      node->left = insert(node->left, key);
      if (node->left->priority > node->priority)
        node = rotate_right(node);
    } else {
      node->right = insert(node->right, key);
      if (node->right->priority > node->priority)
        node = rotate_left(node);
    }
    recalc(node);
    return node;
  }

  double kth(Node* node, int k) const {
    // k is 0-based, assume 0 <= k < size(root)
    int left_size = get_size(node->left);
    if (k < left_size) {
      return kth(node->left, k);
    }
    if (k < left_size + node->count) {
      return node->key;
    }
    return kth(node->right, k - left_size - node->count);
  }

  void clear(Node* node) {
    if (!node) return;
    clear(node->left);
    clear(node->right);
    delete node;
  }


public:
  OrderStatisticTree()
    : root_(nullptr),
      rng_state_(0x9e3779b97f4a7c15ULL) {}

  OrderStatisticTree(const OrderStatisticTree&) = delete;
  OrderStatisticTree& operator=(const OrderStatisticTree&) = delete;

  OrderStatisticTree(OrderStatisticTree&& other) noexcept
    : root_(other.root_), rng_state_(other.rng_state_) {
    other.root_ = nullptr;
  }

  OrderStatisticTree& operator=(OrderStatisticTree&& other) noexcept {
    if (this != &other) {
      clear(root_);
      root_ = other.root_;
      rng_state_ = other.rng_state_;
      other.root_ = nullptr;
    }
    return *this;
  }

  ~OrderStatisticTree() {
    clear(root_);
  }

  void insert(double key) {
    root_ = insert(root_, key);
  }

  int size() const {
    return get_size(root_);
  }

  int count_less(double key) const {
    int result = 0;
    Node* node = root_;
    while (node) {
      if (key <= node->key) {
        node = node->left;
      } else {
        result += get_size(node->left) + node->count;
        node = node->right;
      }
    }
    return result;
  }

  int count_equal(double key) const {
    Node* node = root_;
    while (node) {
      if (key < node->key) node = node->left;
      else if (key > node->key) node = node->right;
      else return node->count;
    }
    return 0;
  }

  double kth(int k) const {
    // No range check here; caller must ensure 0 <= k < size()
    return kth(root_, k);
  }
};

////////////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------
// Exact quantile-based cost: C^{quant}(y_{a..b}; x) = q_{1-x} - q_x
// with O(log n) update via order-statistics tree
// ----------------------------------------------------------

class QuantileCostExact : public TestBase
{
public:
  // x in (0, 0.5] typically
  QuantileCostExact(double x)
    : x_(x),
      tree_(),
      n_(0),
      q_low_(std::numeric_limits<double>::quiet_NaN()),
      q_high_(std::numeric_limits<double>::quiet_NaN())
  {}

  void update(double y) override
  {
    tree_.insert(y);
    ++n_;

    if (n_ == 1) {
      // Single point: both quantiles are that point
      q_low_ = q_high_ = y;
      return;
    }

    // 0-based indices in [0, n_-1]
    int idx_low  = static_cast<int>(std::floor(x_ * (n_ - 1)));
    int idx_high = static_cast<int>(std::floor((1.0 - x_) * (n_ - 1)));

    if (idx_low < 0) idx_low = 0;
    if (idx_high < 0) idx_high = 0;
    if (idx_low >= n_) idx_low = n_ - 1;
    if (idx_high >= n_) idx_high = n_ - 1;

    q_low_  = tree_.kth(idx_low);
    q_high_ = tree_.kth(idx_high);
  }

  double statistic() const override
  {
    if (n_ == 0)
      return std::numeric_limits<double>::quiet_NaN();
    if (n_ <= 20)
      return 0.0; // q_{1-x} = q_x

    return q_high_ - q_low_;
  }

private:
  double x_;
  OrderStatisticTree tree_;
  int n_;
  double q_low_;
  double q_high_;
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////





// ----------------------------------------------------------
// Streaming quantile estimator (P² algorithm) for one p
// ----------------------------------------------------------

class P2Quantile
{
public:
  explicit P2Quantile(double p) : p_(p)
  {
    reset();
  }

  void reset()
  {
    count_ = 0;
  }

  void add(double x)
  {
    // Collect first 5 points
    if (count_ < 5) {
      q_[count_] = x;
      ++count_;
      if (count_ == 5) {
        std::sort(q_, q_ + 5);

        for (int i = 0; i < 5; ++i)
          n_[i] = i;

        ns_[0]  = 0.0;
        ns_[1]  = 2.0 * p_;
        ns_[2]  = 4.0 * p_;
        ns_[3]  = 2.0 + 2.0 * p_;
        ns_[4]  = 4.0;
        dns_[0] = 0.0;
        dns_[1] = p_ / 2.0;
        dns_[2] = p_;
        dns_[3] = (1.0 + p_) / 2.0;
        dns_[4] = 1.0;
      }
      return;
    }

    // General case: update markers
    int k;
    if (x < q_[0]) {
      q_[0] = x;
      k = 0;
    } else if (x < q_[1]) {
      k = 0;
    } else if (x < q_[2]) {
      k = 1;
    } else if (x < q_[3]) {
      k = 2;
    } else if (x < q_[4]) {
      k = 3;
    } else {
      q_[4] = x;
      k = 3;
    }

    // Update positions
    for (int i = k + 1; i < 5; ++i)
      ++n_[i];
    for (int i = 0; i < 5; ++i)
      ns_[i] += dns_[i];

    // Adjust interior markers
    for (int i = 1; i <= 3; ++i) {
      double d = ns_[i] - n_[i];
      if ((d >= 1.0 && n_[i + 1] - n_[i] > 1) ||
          (d <= -1.0 && n_[i - 1] - n_[i] < -1)) {
        int dInt = (d > 0.0) ? 1 : -1;
        double qs = parabolic(i, dInt);
        if (q_[i - 1] < qs && qs < q_[i + 1])
          q_[i] = qs;
        else
          q_[i] = linear(i, dInt);
        n_[i] += dInt;
      }
    }

    ++count_;
  }

  double quantile() const
  {
    if (count_ == 0)
      return std::numeric_limits<double>::quiet_NaN();

    // For the first few points, just use the empirical quantile
    if (count_ <= 5) {
      double tmp[5];
      for (int i = 0; i < count_; ++i)
        tmp[i] = q_[i];
      std::sort(tmp, tmp + count_);

      double pos  = (count_ - 1) * p_;
      int index   = static_cast<int>(std::round(pos));
      if (index < 0) index = 0;
      if (index >= count_) index = count_ - 1;
      return tmp[index];
    }

    // Steady state: p-quantile is q_[2]
    return q_[2];
  }

private:
  double parabolic(int i, int d) const
  {
    return q_[i] + d / static_cast<double>(n_[i + 1] - n_[i - 1]) *
      ( (n_[i] - n_[i - 1] + d) * (q_[i + 1] - q_[i]) / static_cast<double>(n_[i + 1] - n_[i]) +
      (n_[i + 1] - n_[i] - d) * (q_[i] - q_[i - 1]) / static_cast<double>(n_[i] - n_[i - 1]) );
  }

  double linear(int i, int d) const
  {
    return q_[i] + d * (q_[i + d] - q_[i]) / static_cast<double>(n_[i + d] - n_[i]);
  }

  double p_;       // target quantile in (0,1)
  int    n_[5];    // marker positions
  double ns_[5];   // desired positions
  double dns_[5];  // position increments
  double q_[5];    // marker heights
  int    count_{}; // number of observations seen
};

// ----------------------------------------------------------
// Quantile-based cost: C^{quant}(y_{a..b}; x) = q_{1-x} - q_x
// ----------------------------------------------------------

class QuantileCost : public TestBase
{
public:
  // x is the lower quantile level in (0, 0.5], typically
  explicit QuantileCost(double x)
    : lower_(x),
      upper_(1.0 - x),
      q_low_(std::numeric_limits<double>::quiet_NaN()),
      q_high_(std::numeric_limits<double>::quiet_NaN()),
      n_(0)
  {}

  void update(double y) override
  {
    lower_.add(y);
    upper_.add(y);
    ++n_;

    // Cache current quantiles so statistic() is trivial
    q_low_  = lower_.quantile();
    q_high_ = upper_.quantile();
  }

  double statistic() const override
  {
    if (n_ == 0)
      return std::numeric_limits<double>::quiet_NaN();
    if (n_ <= 20)
      return 0.0; // q_{1-x} = q_x = y_1

    // O(1): just use the cached values
    return q_high_ - q_low_;
  }

private:
  P2Quantile lower_; // estimator for q_x
  P2Quantile upper_; // estimator for q_{1-x}
  double q_low_;
  double q_high_;
  int    n_;        // number of points seen
};




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// GvarCost cost:
/// C_n = 1/n sum (y_t - mean)^2, statistic = C_n / sigma^2 ~ chi^2_{n-1} under H0.
/// n <= 10 cannot be rejected
class varCost : public TestBase
{
public:
  varCost()
    : n_(0),
      mean_(0.0),
      M2_(0.0)   // will store sum (y_t - mean)^2
  {}

  void update(double y) override
  {
    ++n_;
    // Welford's online update
    double delta  = y - mean_;
    mean_        += delta / static_cast<double>(n_);
    double delta2 = y - mean_;
    M2_          += delta * delta2; // accumulated SSE around the running mean
  }

  // Returns chi-square statistic: C_n / sigma^2
  double statistic() const override
  {
    if (n_ == 0)
      return std::numeric_limits<double>::quiet_NaN();
    if (n_ <= 10)
      // With one point, df = 0, SSE = 0, so statistic is 0.
      return 0.0;

    return M2_/n_;
  }

  // Optional helper: degrees of freedom of the chi^2
  int dof() const
  {
    return (n_ > 0) ? (n_ - 1) : 0;
  }

private:
  int    n_;      // number of points
  double mean_;   // running mean
  double M2_;     // running sum of squared deviations from mean (SSE)
};





/// WilcoxonCost:
/// - Stores the data y_1, ..., y_n via update(y).
/// - update() maintains, for each split t = 1..n-1,
///     A = {y_1, ..., y_t},  B = {y_{t+1}, ..., y_n}
///   the Wilcoxon rank-sum statistic W_t (centered),
///   and returns max_t |W_t|.
///
/// New pairwise rank contributions update all existing splits in O(n),
/// including exact half-weight handling of ties. statistic() is O(1).

class WilcoxonCost : public TestBase
{
public:
  WilcoxonCost() : doubled_statistic_(0) {}

  void update(double y) override
  {
    const std::size_t old_n = values_.size();
    std::int64_t contribution = 0;
    for (std::size_t i = 0; i < old_n; ++i) {
      contribution += static_cast<std::int64_t>(values_[i] > y) -
        static_cast<std::int64_t>(values_[i] < y);

      // Existing split t = i + 1 gains the contribution of y to its
      // right-hand group. The final contribution is the new last split.
      if (i < split_statistics_.size()) {
        split_statistics_[i] += contribution;
      }
    }

    values_.push_back(y);
    if (old_n > 0) {
      split_statistics_.push_back(contribution);
    }

    doubled_statistic_ = 0;
    for (std::int64_t value : split_statistics_) {
      const std::int64_t magnitude = value < 0 ? -value : value;
      doubled_statistic_ = std::max(doubled_statistic_, magnitude);
    }
  }

  double statistic() const override
  {
    if (values_.empty()) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return 0.5 * static_cast<double>(doubled_statistic_);
  }

private:
  std::vector<double> values_;
  std::vector<std::int64_t> split_statistics_;
  std::int64_t doubled_statistic_;
};


/**
 * MedianMoodCost:
 *
 * Implements a scan version of Mood's median test for a single changepoint.
 *
 * - Data are added via update(y).
 * - statistic() computes:
 *      * global median m of all y's,
 *      * for each split t, counts below/above m in the two groups,
 *      * builds a 2x2 table and computes chi-square,
 *      * returns the maximum chi-square over t.
 *
 * This is robust (median + only above/below information) and
 * detects changes in location (median).
 */
class MedianMoodCost : public TestBase {
public:
  MedianMoodCost() = default;

  void update(double y) override {
    values_.push_back(y);
    tree_.insert(y);
  }

  double statistic() const override {
    const std::size_t n = values_.size();
    if (n < 2) {
      return 0.0;  // not enough data to split
    }

    // The maintained tree returns the same upper median as nth_element(n / 2)
    // without copying and partitioning the complete segment on every call.
    const double med = tree_.kth(static_cast<int>(n / 2));

    const int total_below = tree_.count_less(med);
    const int total_above = static_cast<int>(n) - total_below -
      tree_.count_equal(med);

    const int N_effective   = total_below + total_above;

    // If everything is exactly equal to the median, no information.
    if (N_effective == 0) {
      return 0.0;
    }

    // --- 3) Scan all splits and compute chi-square ---

    double best_chisq = 0.0;
    int prefix_below = 0;
    int prefix_above = 0;

    for (std::size_t t = 1; t < n; ++t) {
      const double value = values_[t - 1];
      if (value < med) ++prefix_below;
      else if (value > med) ++prefix_above;

      // Group A: indices [0, t-1]
      // Group B: indices [t, n-1]
      int a11 = prefix_below;               // A, below
      int a12 = prefix_above;               // A, above
      int a21 = total_below - a11;         // B, below
      int a22 = total_above - a12;         // B, above

      int nA = a11 + a12;
      int nB = a21 + a22;

      // no "signal" in one side -> skip
      if (nA == 0 || nB == 0) {
        continue;
      }

      if (total_below == 0 || total_above == 0) continue;

      // Algebraically identical Pearson chi-square formula for a 2x2 table.
      const double determinant =
        static_cast<double>(a11) * a22 - static_cast<double>(a12) * a21;
      const double N = static_cast<double>(nA + nB);
      const double denominator = static_cast<double>(nA) * nB *
        total_below * total_above;
      const double chisq = N * determinant * determinant / denominator;

      if (chisq > best_chisq) {
        best_chisq = chisq;
      }
    }

    return best_chisq;
  }

private:
  std::vector<double> values_;
  OrderStatisticTree tree_;
};
