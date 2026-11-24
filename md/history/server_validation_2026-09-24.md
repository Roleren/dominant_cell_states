# Historical server validation record

This is a dated record, not the current pipeline specification.

## 2026-09-24 First Full-Server Run and Iteration Review

This repository was developed and validated against a subset of `.fst` files on
a local machine. This date marks the first attempt to run the full pipeline
against the complete backing dataset on the analysis server. Getting there
required fixing several real environment bugs that were invisible on the local
subset (a non-existent `ORFik::txdbFile()` call, `.qs`-versus-`.rds` count
table extensions, a stale ORFik config cached under isolated XDG/BiocFileCache
state, and a project-specific `collection_tables_indexed` path convention that
seven scripts had silently missed) plus one real circular bootstrap dependency
between `manual_translons` and `clean_cds` on a fresh checkout. The EIF4G2 raw
tRNA/UV reprocessing was also independently reproduced on the server by reusing
massiveNGSpipe's already-trimmed reads instead of re-downloading, reproducing
the manuscript's qualitative result.

A full forward-looking review -- what is good, what is not optimal, what is
missing, and prioritized next steps for both the EIF4G2 and post-viral-fatigue
stories -- is in `dominant_cell_state_iteration_review_2026-09-24.md`. The
short version: EIF4G2 is ready for its first lab experiment today; the
fatigue/Long-COVID story has comparably strong computational candidates but
still has no lab go/no-go plan, and writing one is the single highest-priority
next step.
