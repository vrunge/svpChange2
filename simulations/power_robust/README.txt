ROBUST POWER STUDY

This folder contains heavy-tailed (Student-t with two degrees of freedom)
power studies.  It compares PELT, Wilcoxon, Median-Mood, rfpop when the
optional robseg package is available, and the new FOCUS TRUE/TRUE and
multiscale SVP methods.

Scripts:
  3_robust_study.R  original heavy-tail paper analysis

The new-method implementation is in the sibling Gaussian-power script
5_new_gaussian_and_robust_studies.R; its robust run writes results into this
folder. Existing CSV/RDS files and 3_*.pdf files are retained paper outputs.
Run from the svpChange2 package root.
