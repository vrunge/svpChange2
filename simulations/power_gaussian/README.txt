GAUSSIAN POWER STUDY

This folder contains independent-Gaussian piecewise-constant signal studies.
Scenarios are none, up, updown, and rand1, with varying jump size and repeated
Monte Carlo samples.  Metrics are precision, recall, F1, number of segments,
MSE, and changepoint distributions.

Scripts:
  1_power_study.R                          original SVP/PELT study
  5_new_gaussian_and_robust_studies.R      new FOCUS and multiscale methods

The second script also contains the heavy-tail extension used by the robust
study.  Existing CSV/RDS files and 1_*.pdf files are paper outputs.  Run from
the svpChange2 package root; large runs require SVP_RUN_SIMULATIONS=true.
The shared `simulations/metrics.R` utility adds localization error and exact
number-of-change-points recovery metrics, alongside precision, recall, and F1.
