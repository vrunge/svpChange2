SMUCE/SVP COMPARISON

compare_smuce_svp.R compares the SMUCE estimator from the CRAN stepR package
with `svp_smuce()`, the constrained Gaussian SVP implementation. It uses
random heterogeneous Gaussian piecewise-constant signals and compares the
returned partitions directly.

Install the dependency first:

  install.packages("stepR")

Then run from the svpChange2 package root:

  source("simulations/smuce_comparison/compare_smuce_svp.R")
  answer <- compare_smuce_svp()
  answer$smuce
  answer$svp_smuce

The random seed is unset by default, so every run generates a new signal.
`run_smuce_replicates()` repeats the comparison and reports exact equality.

`time_smuce_svp.R` benchmarks `stepR::smuceR()` and
`svp_smuce_cpp()` for increasing sample sizes. The C++ routine updates
SMUCE envelopes incrementally and prunes invalid extensions; its worst-case
complexity is O(n^3), although pruning can substantially reduce runtime.
The stepR implementation uses specialized multiscale algorithms and is
typically approximately quadratic for n <= 1000 (larger n uses dyadic
interval lengths). Timings depend on the signal and hardware.

The function returns only:

  $smuce  the SMUCE segment-end indices
  $svp_smuce  the constrained SVP segment-end indices

The report in smuce_comparison/SVP_SMUCE_report.tex explains why the basic
data-only SVP is not universally identical to SMUCE: exact equivalence also
requires the same constrained segment-parameter optimization.
