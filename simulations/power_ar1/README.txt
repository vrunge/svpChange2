AR1 POWER STUDY

This folder contains power studies for piecewise-constant marginal means with
stationary AR(1) noise.  It compares independent-Gaussian SVP, exact AR1 SVP
with known rho, and exact AR1 SVP with robustly estimated rho.

The study varies signal scenario, jump size, autocorrelation, and replication.
Outputs include accuracy/MSE/segment-count plots, changepoint distributions,
and AR1 signal scenarios.  Run 6_ar1_study.R from the svpChange2 package root.
Large Monte Carlo runs are opt-in through SVP_RUN_SIMULATIONS=true.
