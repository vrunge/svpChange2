AR1 TIME-COMPLEXITY STUDY

This folder is reserved for paper-style runtime experiments under stationary
AR(1) noise.  The two experiments mirror the Gaussian timing study: (1) vary
the sequence length with no changes, and (2) fix the sequence length while
varying the number of changes.  The planned comparison is independent SVP,
AR1 SVP with known rho, and AR1 SVP with robustly estimated rho.

Run 8_time_ar1.R from the svpChange2 package root.  It is opt-in through
SVP_RUN_SIMULATIONS=true and writes CSV results and PDFs to this folder.
