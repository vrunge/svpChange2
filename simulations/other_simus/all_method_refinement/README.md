# All-method boundary-refinement test

This isolated experiment applies the former boundary-refinement idea to every
method in each power framework. It does not change or write into any primary
`power_*` or `time_*` study.

Each internal boundary is independently relocated between its two neighboring
detected boundaries, with candidate splits at least five observations from
those neighbors. Gaussian data use the mean-CUSUM objective, robust data use
the rank-CUSUM objective, and AR(1) data use the pure-AR(1) likelihood cost with
the study value `rho = 0.8`.

Because neighboring boundaries are refined independently, two of them can
occasionally move to the same location. Boundary normalization then collapses
the duplicate, reducing the number of segments. This is rare in the Gaussian
and AR(1) studies but common for PELT under the heavy-tailed robust study,
where many initial detections are tightly spaced.

The test starts from the saved primary-study detections and reconstructs the
exact simulated observations from their documented deterministic seeds. This
isolates boundary refinement from detection.

From the package root, run:

```sh
Rscript simulations/other_simus/all_method_refinement/run_refinement_test.R
```

Outputs are written to this directory under `power_gaussian/`,
`power_robust/`, and `power_ar1/`. Each contains `results.rds`, `results.csv`,
and a `plots/` directory with the seven standard power-study PDF figures.
