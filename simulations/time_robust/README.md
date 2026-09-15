# Robust runtime study

`run_time.R` benchmarks PELT and SVP validity tests under Student-t(2) noise.
Method labels follow the power-study legend style:

- `PELT`
- `SVP Wilcoxon / right / c = 2.000`
- `SVP MedianMood / right / c = 2.000`
- `SVP Gaussian / both / c = 1.500`

The Wilcoxon and Median-Mood calls explicitly use `subtests = "right"`; the
Gaussian call uses `subtests = "both"`.

Both timing experiments enforce a maximum sequence length of 1,000.

Run from the package root:

```r
Sys.setenv(SVP_RUN_SIMULATIONS = "true")
source("simulations/time_robust/run_time.R")
```

The separate C++ validity-update benchmark is under
`simulations/other_simus/robust_validity_benchmark/`.
