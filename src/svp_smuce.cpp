#include <Rcpp.h>
using namespace Rcpp;

//' C++ SVP with SMUCE validity and constrained Gaussian cost
//' @param y Numeric observations.
//' @param q SMUCE threshold.
//' @param sigma2 Known Gaussian variance.
//' @return Integer segment-end indices.
// [[Rcpp::export]]
IntegerVector svp_smuce_cpp(NumericVector y, double q, double sigma2 = 1.0) {
  int n = y.size();
  if (n < 1 || !R_finite(q) || q < 0 || !R_finite(sigma2) || sigma2 <= 0)
    stop("invalid arguments");
  std::vector<double> cs(n + 1, 0.0), cs2(n + 1, 0.0);
  for (int i = 0; i < n; ++i) {
    cs[i+1] = cs[i] + y[i];
    cs2[i+1] = cs2[i] + y[i]*y[i];
  }
  std::vector<int> K(n+1, n+1), prev(n+1, -1);
  std::vector<double> C(n+1, R_PosInf);
  K[0] = 0; C[0] = 0.0;
  // For each start s, extend the segment and update constraints ending at t.
  for (int s = 0; s < n; ++s) {
    // BEFORE-invalid pruning: a start with no finite optimal prefix can
    // never contribute to a later partition.
    if (K[s] > n || !R_finite(C[s])) continue;
    double lo = R_NegInf, hi = R_PosInf;
    for (int t = s; t < n; ++t) {
      int m = t - s + 1;
      // Add every new subinterval [u,t].
      for (int u = s; u <= t; ++u) {
        int len = t-u+1;
        double mean = (cs[t+1]-cs[u])/len;
        double rad = std::sqrt(sigma2/len) *
          (q + std::sqrt(2.0*std::log(std::exp(1.0)*n/len)));
        lo = std::max(lo, mean-rad);
        hi = std::min(hi, mean+rad);
      }
      // AFTER-invalid pruning: SMUCE feasibility is hereditary under
      // extension. Once the admissible theta interval is empty, every
      // longer segment with this same start is invalid as well.
      if (lo > hi) break;
      double mean_seg = (cs[t+1]-cs[s])/m;
      double theta = std::min(std::max(mean_seg, lo), hi);
      double rss = (cs2[t+1]-cs2[s]) - 2.0*theta*(cs[t+1]-cs[s]) + m*theta*theta;
      double cost = C[s] + rss/sigma2;
      int nk = K[s] + 1;
      if (nk < K[t+1] || (nk == K[t+1] && cost < C[t+1])) {
        K[t+1] = nk; C[t+1] = cost; prev[t+1] = s;
      }
    }
  }
  std::vector<int> ends; int t = n;
  while (t > 0) { if (prev[t] < 0) stop("no valid partition"); ends.push_back(t); t = prev[t]; }
  std::reverse(ends.begin(), ends.end());
  return IntegerVector(ends.begin(), ends.end());
}
