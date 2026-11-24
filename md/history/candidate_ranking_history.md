# Candidate ranking history

The current algorithm combines geometry/novelty, buffering, matched-control
support, design/study replication, specificity, bootstrap certainty, and
shrunk effects in one additive candidate priority score. Its formulas and
weights have not changed in the October 2026 cleanup.

Earlier development exported geometry-only rankings, then buffering (v2),
matched-control (v3), and replication-adjusted (v4) rankings. These were
successive implementations of one priority ranking, not independent validation
models. The intermediate copies, rank columns, score aliases, and exports were
removed on 2026-10-09. The final table is now
`dominant_next_model_candidate_rankings.csv`; `candidate_rank` and
`candidate_score` replace the final v4 names. Initial geometry rank remains
`preliminary_candidate_rank` for interpretation. Intermediate score sorts are
retained to preserve tie ordering in the final ranking.

Other statistical layers remain separate because branch allocation, pooling,
held-out prediction, and cross-study robustness answer different questions.
Old output files are archived with the benchmark artifacts rather than used
as live inputs. Historical scientific notes may cite their original names.
