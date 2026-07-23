# AR(1) power study

Run from the `svpChange2` package root:

```r
Sys.setenv(SVP_RUN_SIMULATIONS = "true")
source("simulations/power_ar1/run_power.R")
```

This clean study follows `svp_ar1_power_study.pdf`: `n=600`, `rho=0.8`,
marginal noise variance one, four signal scenarios, jumps from 1 to 3 by 0.2,
and 30 replications. It compares inflated Gaussian PELT, approximate
AR(1)-likelihood optimal partitioning, AR1Focus SVP, and AR1Focus SVP with
TRUE/TRUE pruning. Both SVP methods use `3.75*log(n)`.

The six PDFs parallel the Gaussian and heavy-tail outputs. Changepoint
distributions use jump 1.4, which lies on the AR(1) simulation grid.
