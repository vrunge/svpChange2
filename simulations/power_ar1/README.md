# AR(1) power study

Run from the `svpChange2` package root:

```r
Sys.setenv(SVP_RUN_SIMULATIONS = "true")
source("simulations/power_ar1/run_power.R")
```

On macOS/Linux the default uses all detected physical cores except one through
`parallel::mclapply()`. Set `run_and_save_ar1(workers = 1L)` for a sequential
run or pass another explicit worker count.

The separately calibrated TRUE/TRUE AR1Focus constant is 1.75; see
`true_true_calibration.csv`. Its boundaries are refined with the local exact
AR1 likelihood while preserving their number. In the final study this method
has rand1 F1 0.789 (AR1 approximate PELT: 0.778) and no-change F1 0.994.
The final design uses 100 replications for each of the four patterns and 21
jump sizes from 1 to 3 in steps of 0.1, giving 8,400 simulated series and
33,600 method fits.

This clean study follows `svp_ar1_power_study.pdf`: `n=600`, `rho=0.8`,
marginal noise variance one, four signal scenarios, jumps from 1 to 3 by 0.2,
and 30 replications. It compares inflated Gaussian PELT, approximate
AR(1)-likelihood optimal partitioning, AR1Focus SVP, and AR1Focus SVP with
TRUE/TRUE pruning. Both SVP methods use `3.75*log(n)`.

The six PDFs parallel the Gaussian and heavy-tail outputs. Changepoint
distributions use jump 1.4, which lies on the AR(1) simulation grid.
