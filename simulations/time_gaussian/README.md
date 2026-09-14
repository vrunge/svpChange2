# Gaussian runtime study

`run_time.R` benchmarks the same canonical method names used in the Gaussian
power-study legend:

- `PELT`
- `SVP BIC`
- `SVP BIC calibrated`
- `SVP BIC multiscale` (`subtests = "both"`)

It varies sequence length under no change and the number of changes at fixed
length. Run it from the package root:

```r
Sys.setenv(SVP_RUN_SIMULATIONS = "true")
source("simulations/time_gaussian/run_time.R")
```

The progress bar counts completed method fits. Results and three runtime plots
are saved in this folder. The checked-in `gaussian_time_*` files are the
pre-shortcut baseline; `run_shortcut_comparison.R` writes the version 0.2.0
benchmark under a separate `gaussian_time_shortcuts_*` prefix and checks that
PELT and the two right-mode methods retain their detected segment counts. The
multiscale comparison is reported separately because the baseline predates the
current `subtests` API and later multiscale candidate corrections. Superseded
timing implementations and their outputs are retained under
`simulations/other_simus/time_gaussian_legacy/`.
