ROBUST TIME-COMPLEXITY STUDY

This folder is reserved for paper-style runtime experiments under Student-t
noise (two degrees of freedom).  The two experiments mirror the Gaussian timing
study: (1) vary sequence length with no changes, and (2) vary the number of
changes at fixed sequence length.  The planned methods are PELT, Wilcoxon,
Median-Mood, FOCUS TRUE/TRUE, and multiscale SVP.

Run 9_time_robust.R from the svpChange2 package root.  It is opt-in through
SVP_RUN_SIMULATIONS=true and writes CSV results and PDFs to this folder.
