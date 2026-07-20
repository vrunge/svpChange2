TIME COMPLEXITY SIMULATIONS

This folder contains the paper-style runtime studies.  The experiments vary
the sequence length and the number of true changes, and compare PELT with the
SVP Gaussian, FOCUS TRUE/TRUE, and multiscale configurations.  The all-cost
script also compares robust/nonparametric tests.

Scripts:
  2_complexity_study.R                         original PELT/SVP benchmark
  2_complexity_study_allCost.R                 cost-test benchmark
  2_complexity_study_focus_pelt_figures.R      paper figure reproduction
  7_new_methods_complexity.R                   new svpChange2 methods

CSV files contain previously generated benchmark results; PDF files contain
the corresponding paper figures.  Run from the svpChange2 package root.  The
new script is opt-in through SVP_RUN_SIMULATIONS=true.
