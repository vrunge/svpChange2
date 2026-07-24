# Heavy-tail power study

Run from the `svpChange2` package root:

```r
Sys.setenv(SVP_RUN_SIMULATIONS = "true")
source("simulations/power_robust/run_power.R")
```

On macOS/Linux the default uses all detected physical cores except one through
`parallel::mclapply()`. Set `run_and_save_robust(workers = 1L)` for a
sequential run or pass another explicit worker count.

The separately calibrated TRUE/TRUE Wilcoxon constant is 1.75; see
`true_true_calibration.csv`. The reported boundaries are the direct output of
`SVP()`, without post-processing. In the final study this method has rand1 F1
0.595 (RFPOP paper: 0.688) and
no-change F1 0.998.

The design reproduces the Student-t(2) study in `SVP_Paper.pdf`: `n=1000`,
four scenarios, jumps from 0.1 to 4, 100 replications, and the paper's PELT,
RFPOP, Median-Mood SVP, and Wilcoxon SVP configurations.

`SVP Wilcoxon multiscale` is the Wilcoxon method with TRUE/TRUE pruning. The original
`RFPOP paper` result is
retained in the result files but excluded from the primary legend. `RFPOP`
uses the smallest penalty multiplier
attaining at least 98% no-change recovery in an independent 200-replicate null
calibration; its calibration table is saved as `rfpop_null_calibration.csv`.
This makes the false-positive trade-off explicit rather than altering the paper
curve without disclosure.

The seven PDFs and result files have the same meanings as in the Gaussian
folder.
