# Simulations for `svpChange2`

The three power-study folders are self-contained reproductions and extensions
of the studies in `SVP_Paper.pdf`:

- `power_gaussian/`: Gaussian Figures 1, 2, 6, and 8, with the fixed
  `SVP FOCUS TRUE/TRUE c=1.8` curve.
- `power_robust/`: Student-t(2) Figures 3 and 7 plus companion figures, with
  fixed `SVP Wilcoxon TRUE/TRUE c=1.8` and both paper/null-matched RFPOP.
- `power_ar1/`: the analogous AR(1) study with inflated PELT, approximate
  AR(1) optimal partitioning, AR1Focus SVP, and AR1Focus TRUE/TRUE SVP.

Each folder contains only `run_power.R`, `README.md`, `results.csv`,
`results.rds`, and six PDFs under `plots/` (plus the RFPOP calibration table in
the robust folder). Shared signal generation, metrics, and plotting code lives
in `power_common.R`.

Run a study from the package root:

```r
Sys.setenv(SVP_RUN_SIMULATIONS = "true")
source("simulations/power_gaussian/run_power.R")
source("simulations/power_robust/run_power.R")
source("simulations/power_ar1/run_power.R")
```

The remaining folders contain runtime, pruning, SMUCE, and real-data studies.
