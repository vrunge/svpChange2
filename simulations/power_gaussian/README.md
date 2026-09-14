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
The two Section 5.1 FOCuS configurations explicitly use `subtests = "right"`,
while the multiscale extension uses `subtests = "both"`.
Every design row uses an explicit Mersenne-Twister seed, so sequential and
fork-parallel runs generate the same Gaussian samples and all methods are
compared on exactly the same observations.

To reproduce into a separate directory without replacing the historical files,
source `simulations/power_gaussian/run_reproduce_v020.R` and call
`run_and_save_v020_matched(workers = 1L)`; the output is written to
`reproduced_v020_matched/`.

The TRUE/TRUE constant was selected on separate calibration seeds subject to
high no-change F1; the calibration table is `true_true_calibration.csv`.
The reported boundaries are the direct output of `SVP()`, without
post-processing. In the current-source matched reproduction, rand1 F1 is
0.658 for PELT, 0.633 for SVP BIC, 0.636 for the calibrated right-only fit,
and 0.647 for the multiscale both-sided fit; the corresponding no-change F1
values are 0.995, 1.000, 0.988, and 0.984.

`results.rds` preserves detected changepoints; `results.csv` contains scalar
metrics. The seven PDFs in `plots/` are scenarios, F1, precision/recall,
changepoint distributions at jump 0.6, probability of recovering the correct
number of changepoints, signal MSE, and conditional localization error.

`run_both_calibrated.R` and `reproduced_v020_both_calibrations/` contain an
independent null calibration and full-grid comparison for the current
`subtests = "both"` implementation. It selects `1.9 log(n)` as the smallest
tested multiplier attaining null F1 0.99; `1.8 log(n)` remains the historical
reference. `run_and_save_both_comparison(workers = 1L)` regenerates the
side-by-side `1.9` and `2.0` grids in that directory.
