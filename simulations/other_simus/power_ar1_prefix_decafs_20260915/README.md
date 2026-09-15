# AR(1) power study outputs from before the DeCAFS fix

Snapshot of `simulations/power_ar1/` as it stood on 2026-09-15, immediately
before two changes:

1. DeCAFS 3.3.6 corrected the precision it puts on the first observation in its
   dynamic programme, from `1/(sdNu^2 (1-phi^2))` to the stationary AR(1) value
   `(1-phi^2)/sdNu^2` required by equation (4) of the DeCAFS paper. Against the
   older build, DeCAFS calibrated to `beta = 7.75*log(n)` at the study's null
   target and trailed PELT; with the fix it calibrates to `2.75*log(n)` and
   leads the study.
2. `SVP()` gained a `cost` argument, and the study switched its SVP row from the
   Gaussian squared-error cost to the AR(1) innovation cost.

Files:

* `results.rds` — the 33,600-row power grid produced before both changes.
* `ar1_calibration.rds` — the matching calibration, with `decafs_constant` 7.75,
  `approximate_constant` 0.8, `svp_constant` 1.85, null target F1 0.99.
* `decafs_null_calibration.csv`, `decafs_power_tradeoff.csv` — standalone DeCAFS
  tuning sweeps against the older build. No current script regenerates them and
  their constants (3.5 to 12) only make sense for the uncorrected cost, so they
  are kept here rather than in the active study folder.
