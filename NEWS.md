# svpChange2 0.2.0

- Removed the redundant public `subtests = "left"` option from `SVP()` and
  `svp0()`. The supported values are now `"none"`, `"right"`, and `"both"`.
- Added exact lexicographic shortcuts for the three surviving `SVP()` modes.
  These reduce validity-test or SSE work without changing the returned
  partitions, dynamic-programming table, or candidate diagnostics.
- Updated the pruning documentation and runtime simulations.
