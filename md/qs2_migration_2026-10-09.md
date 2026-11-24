# qs2 migration and verification

Date: 2026-10-09.

All 2,389 discovered local RDS cache/checkpoint objects were read, exported using qs2, and read back with attributes checked before their originals were moved into an ignored recovery directory. Two identical signature aliases share an output. Converted source-file sizes total 398.55 MiB; the corresponding qs files total 317.21 MiB, about 20.4% smaller. These totals count the two source aliases separately.

Production transcript coverage, per-gene derived measurements, manual activity/ranges, pipeline manifests, calibration checkpoints and relevant read-length diagnostics now use `.qs`. Reads and atomic writes use a common serialization helper, with four threads by default and `DOMINANT_QS_THREADS` configurable. External Bio_data files were not converted: their loaders prefer an existing current `.qs` sibling, otherwise retain compatible external reading. The MegaFST index and used count-table update markers remain part of cache invalidation.

A tiny in-memory legacy coverage-key calculation still writes a temporary RDS solely to preserve the existing key convention. Benchmark and migration utilities intentionally retain RDS support for comparison/recovery. There are no active project RDS cache outputs outside the recovery artifacts.

The prior FMR1/BOD1L1/FN1 serialization benchmark found summed median read/write times of 0.716/4.607 seconds for gzip RDS and 0.106/0.180 seconds for qs2 with four threads: approximately 6.75× faster reading and 25.6× faster writing. All 54 benchmark round trips preserved values and attributes. This measures serialization on warm local storage, not total biological-analysis runtime; extraction, coverage mapping and model fitting remain separate costs.

Verification: 162 package assertions pass with no failures, warnings or skips. Script/package parsing and whitespace checks pass. The independent output audit verifies unchanged core scientific CSVs and retained feature values, including FMR1 and BOD1L1. Pipeline cache metadata and derived object signatures were migrated only after conversion/read-back and existing scientific/external-input checks.

During producer validation, an external ORFik update changed the reported version from 1.29.12 to 1.29.14 and interrupted a long-lived session. Its version guard then correctly attempted a cold cache refresh. The user confirmed: “Yes, just pushed a fix, it changes nothing for you.” This specific compatibility acknowledgement is recorded in the local migration artifacts and used to rekey the affected metadata without changing cached measurements. Future package versions still invalidate normally. The interrupted cold refresh was stopped; it is not a benchmark of qs2 performance.

Local artifacts: `results/qs2_migration_2026-10-09/`, including conversion manifest, source snapshots, test logs, output equivalence log, runtime compatibility acknowledgement and producer/warm-run verification logs. Original objects remain recoverable under `legacy_rds/`. No generated cache files were staged or committed.

Final producer verification completed successfully: manual translons 3.3 s, FST/transcript selection 31.3 s, clean CDS 35.9 s, uORF inference 73.7 s, and RDG tables 84.4 s. A subsequent normal run skipped all 50 core/RDG steps in 6.05 s. The first normal run in a fresh R session also skipped all 50 steps, taking 20.95 s including runtime/package setup.

The remaining producer costs are gene iteration, library/transcript bookkeeping, inference and table/model construction. qs2 accelerates object serialization, but cannot remove those costs; normal warm runs avoid them through the step manifests. All five forced producers and the final output audit passed.
