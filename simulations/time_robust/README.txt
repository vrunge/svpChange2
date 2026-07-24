ROBUST TIME-COMPLEXITY STUDY

This folder is reserved for paper-style runtime experiments under Student-t
noise (two degrees of freedom).  The two experiments mirror the Gaussian timing
study: (1) vary sequence length with no changes, and (2) vary the number of
changes at fixed sequence length.  The planned methods are PELT, Wilcoxon,
Median-Mood, FOCUS TRUE/TRUE, and multiscale SVP.

Run 9_time_robust.R from the svpChange2 package root.  It is opt-in through
SVP_RUN_SIMULATIONS=true and writes CSV results and PDFs to this folder.

SECTION 4.2 VALIDITY-UPDATE MICROBENCHMARK

benchmark_validity_updates.R compiles benchmark_validity_updates.cpp and times
the actual WilcoxonCost and MedianMoodCost C++ classes.  The complete operation
is update(new_observation) followed by statistic(), as used by SVP.  Run:

SVP_RUN_SIMULATIONS=true Rscript \
  simulations/time_robust/benchmark_validity_updates.R

The benchmark writes raw timings, median summaries, fitted log-log slopes, and
validity_update_benchmark.pdf.  It also times update() alone to show that Mood's
tree insertion is sublinear, while its complete validity update is linear
because statistic() scans all split positions.
