SVP PRUNING TIME-COMPLEXITY STUDY

This folder isolates the computational effect of the four SVP pruning modes:

  subtests
  --------
  none
  right
  left
  both

Each strategy is compared with PELT, following the two runtime experiments
used in the SVP paper:

  1. no-change Gaussian signals with increasing sequence length;
  2. fixed sequence length with an increasing number of changes.

Run 10_time_pruning.R from the svpChange2 package root. The full benchmark is
opt-in through SVP_RUN_SIMULATIONS=true. Results are written as CSV/RDS files
and runtime figures are written to plots/.
