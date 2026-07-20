SMUCE/SVP COMPARISON

compare_smuce_svp.R compares the SMUCE estimator from the CRAN stepR package
with the multiscale Gaussian SVP implementation. It uses independent Gaussian
data, known variance, the same global multiscale calibration q, and disables
directional pruning.

Install the dependency first:

  install.packages("stepR")

Then run from the svpChange2 package root:

  source("simulations/smuce_comparison/compare_smuce_svp.R")
  answer <- compare_smuce_svp()
  answer$equal_partition

The function returns only:

  $smuce  the SMUCE segment-end indices
  $svp    the SVP segment-end indices

The report in smuce_comparison/SVP_SMUCE_report.tex explains why the basic
data-only SVP is not universally identical to SMUCE: exact equivalence also
requires the same constrained segment-parameter optimization.
