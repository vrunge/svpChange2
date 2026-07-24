# Gaussian power study

Run from the `svpChange2` package root:

```r
Sys.setenv(SVP_RUN_SIMULATIONS = "true")
source("simulations/power_gaussian/run_power.R")
```

On macOS/Linux the default uses all detected physical cores except one through
`parallel::mclapply()`. Set `run_and_save_gaussian(workers = 1L)` for a
sequential run or pass another explicit worker count.

The design reproduces the Gaussian power study in `SVP_Paper.pdf`: `n=1000`,
four signal scenarios, jumps from 0.1 to 2, 100 replications, Gaussian variance
one, and tolerance `round(0.0025*n)`. The canonical method names are `PELT`,
`SVP BIC`, `SVP BIC calibrated`, and `SVP BIC multiscale`.
Every design row uses an explicit Mersenne-Twister seed, so sequential and
fork-parallel runs generate the same Gaussian samples and all methods are
compared on exactly the same observations.

The TRUE/TRUE constant was selected on separate calibration seeds subject to
high no-change F1; the calibration table is `true_true_calibration.csv`.
The reported boundaries are the direct output of `SVP()`, without
post-processing. In the final common-sample study its rand1 F1 is 0.647 versus
0.658 for PELT, with no-change F1 0.997.

`results.rds` preserves detected changepoints; `results.csv` contains scalar
metrics. The seven PDFs in `plots/` are scenarios, F1, precision/recall,
changepoint distributions at jump 0.6, probability of recovering the correct
number of changepoints, signal MSE, and conditional localization error.
