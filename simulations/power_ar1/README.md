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
`true_true_calibration.csv`. The reported boundaries are the direct output of
`SVP()`, without post-processing. In the final study this method has rand1 F1
0.751 (AR1 approximate PELT: 0.778) and no-change F1 0.994.
The final design uses 100 replications for each of the four patterns and 21
jump sizes from 1 to 3 in steps of 0.1, giving 8,400 simulated series and
42,000 method fits.

This clean study follows `svp_ar1_power_study.pdf`: `n=600`, `rho=0.8`,
marginal noise variance one, four signal scenarios, jumps from 1 to 3 by 0.1,
and 100 replications. It compares inflated Gaussian PELT, approximate
AR(1)-likelihood optimal partitioning, DeCAFS with its random-walk component
disabled (`sdEta=0`), AR1Focus SVP, and AR1Focus SVP with TRUE/TRUE pruning.
The TRUE/TRUE method is labeled `SVP AR1Focus multiscale`.
DeCAFS uses `4.75*log(n)`, the smallest multiplier attaining at least 95%
no-change recovery in an independent 1,000-replicate calibration. The
calibration is reproduced by `calibrate_decafs.R`, and its complete grid is
saved in `decafs_null_calibration.csv`. The null-versus-power comparison on
the complete study design is saved in `decafs_power_tradeoff.csv`. Relative
to the more conservative multiplier 6.25, multiplier 4.75 improves mean F1
from 0.701 to 0.740 for `up`, 0.335 to 0.454 for `updown`, and 0.639 to 0.697
for `rand1`, while retaining 95.4% correct no-change recovery on the study
samples. The pure AR(1) specification remains
`sdEta=0`, `sdNu=sqrt(1-rho^2)`, and `phi=rho`.
The standard AR1Focus method uses `3.75*log(n)` and its multiscale variant uses
the separately calibrated `1.75*log(n)`.

The seven PDFs parallel the Gaussian and heavy-tail outputs. Changepoint
distributions use jump 1.4, which lies on the AR(1) simulation grid.
