# Gaussian runtime study

`run_time.R` benchmarks the same canonical method names used in the Gaussian
power-study legend:

- `PELT`
- `SVP BIC`
- `SVP BIC calibrated`
- `SVP BIC multiscale` (TRUE/TRUE pruning)

It varies sequence length under no change and the number of changes at fixed
length. Run it from the package root:

```r
Sys.setenv(SVP_RUN_SIMULATIONS = "true")
source("simulations/time_gaussian/run_time.R")
```

The progress bar counts completed method fits. Results and three runtime plots
are saved in this folder. Superseded timing implementations and their outputs
are retained under `simulations/other_simus/time_gaussian_legacy/`.
