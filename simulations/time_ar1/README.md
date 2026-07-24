# AR(1) runtime study

`run_time.R` compares `SVP BIC`, `SVP AR1`, and
`SVP AR1 estimated rho` under stationary AR(1) noise. It varies sequence
length under no change and the number of changes at fixed length.

Run from the package root:

```r
Sys.setenv(SVP_RUN_SIMULATIONS = "true")
source("simulations/time_ar1/run_time.R")
```

The progress bar counts completed method fits. Development-time AR(1)
diagnostics are under `simulations/other_simus/ar1_method_comparisons/`.
