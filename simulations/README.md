# Reproducible simulation studies

Run every script from the `svpChange2` package root. The six primary studies
are deliberately arranged as three power studies and three runtime studies:

- `power_gaussian/`
- `power_robust/`
- `power_ar1/`
- `time_gaussian/`
- `time_robust/`
- `time_ar1/`

Each primary folder has one entry point (`run_power.R` or `run_time.R`), one
README, and its generated results. Supporting experiments and historical
scripts live under `other_simus/`.

Shared code is defined once:

- `power_common.R` contains power-study simulation, metrics, plotting, and
  saving infrastructure.
- `time_common.R` contains timing grids, progress reporting, signal helpers,
  plotting, and saving infrastructure.

Method names stored in result files are the same names used in plot legends.
Every stochastic task receives a deterministic seed based on the study seed
and its design-row index. Consequently, sequential and fork-parallel power
runs generate the same simulated datasets.

## Running a study

Set the opt-in environment variable and source one entry point:

```r
Sys.setenv(SVP_RUN_SIMULATIONS = "true")
source("simulations/power_gaussian/run_power.R")
source("simulations/time_gaussian/run_time.R")
```

For a quick or sequential run, source the script with simulations disabled and
call its `run_*()` function with smaller arguments. Power studies accept
`workers = 1L`; timing studies display a progress bar.

Each `run_power.R` writes `results.rds`, `results.csv`, and all seven PDFs in
its own `plots/` directory. Each `run_time.R` writes its RDS/CSV result pair
and all three PDFs in its own `plots/` directory. Timing values equal to zero
at the operating-system timer resolution remain in the result files but are
excluded from the two logarithmic-time plots.

| Study folders | Entry point | PDFs produced |
|---|---|---|
| `power_gaussian`, `power_robust`, `power_ar1` | `run_power.R` | `01_scenarios.pdf` through `07_loc_error.pdf` |
| `time_gaussian` | `run_time.R` | `gaussian_time_time_vs_n.pdf`, `gaussian_time_time_vs_k.pdf`, `gaussian_time_time_vs_detected.pdf` |
| `time_robust` | `run_time.R` | `robust_time_time_vs_n.pdf`, `robust_time_time_vs_k.pdf`, `robust_time_time_vs_detected.pdf` |
| `time_ar1` | `run_time.R` | `ar1_time_time_vs_n.pdf`, `ar1_time_time_vs_k.pdf`, `ar1_time_time_vs_detected.pdf` |

The AR(1) power study requires the optional `DeCAFS` package. Robust power
comparisons additionally require `robseg`.
