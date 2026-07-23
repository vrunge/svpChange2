# Heavy-tail power study

Run from the `svpChange2` package root:

```r
Sys.setenv(SVP_RUN_SIMULATIONS = "true")
source("simulations/power_robust/run_power.R")
```

The design reproduces the Student-t(2) study in `SVP_Paper.pdf`: `n=1000`,
four scenarios, jumps from 0.1 to 4, 100 replications, and the paper's PELT,
RFPOP, Median-Mood SVP, and Wilcoxon SVP configurations.

`SVP Wilcoxon TRUE/TRUE c=1.8` is added. The original `RFPOP (paper)` curve is
retained exactly. A second RFPOP curve uses the smallest penalty multiplier
attaining at least 98% no-change recovery in an independent 200-replicate null
calibration; its calibration table is saved as `rfpop_null_calibration.csv`.
This makes the false-positive trade-off explicit rather than altering the paper
curve without disclosure.

The six PDFs and result files have the same meanings as in the Gaussian folder.
