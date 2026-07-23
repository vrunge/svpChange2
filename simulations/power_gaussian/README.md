# Gaussian power study

Run from the `svpChange2` package root:

```r
Sys.setenv(SVP_RUN_SIMULATIONS = "true")
source("simulations/power_gaussian/run_power.R")
```

The design reproduces the Gaussian power study in `SVP_Paper.pdf`: `n=1000`,
four signal scenarios, jumps from 0.1 to 2, 100 replications, Gaussian variance
one, and tolerance `round(0.0025*n)`. The original PELT, SVP, and SVP (BIC)
methods are retained. `SVP FOCUS TRUE/TRUE c=1.8` is added to every figure.

`results.rds` preserves detected changepoints; `results.csv` contains scalar
metrics. The six PDFs in `plots/` are scenarios, F1, precision/recall,
changepoint distributions at jump 0.6, probability of recovering the correct
number of changepoints, and signal MSE.
