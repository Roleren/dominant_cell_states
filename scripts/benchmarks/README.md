# One-off benchmarks and migrations

Not part of the pipeline: the runner never sources these files. They are kept to reproduce the dated reports in
`md/` (ORFik accessor/mapping benchmarks 2026-10-09/10, qs2 migration and serialization benchmark 2026-10-09,
incremental-measurement benchmark 2026-10-09). Run them from the project root, e.g.
`env -u LC_ALL R --vanilla -q -f scripts/benchmarks/benchmark_dominant_incremental_measurements.R`.
