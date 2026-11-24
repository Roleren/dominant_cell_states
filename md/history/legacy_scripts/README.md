# Retired count-table workflow

These two scripts implemented the original count-table-only state scorer and
its enrichment scan. They were not registered in the current pipeline and had
no active callers. They are retained here as historical code, outside the
active script directory. Use `scripts/dominant_cell_states_fst_clean_cds.R`
and the adjusted/control-aware enrichment steps for current analyses.

They are not supported entry points. Historical relative paths and outputs
have deliberately not been modernized. The HTML-widget project-root detector
now identifies the active clean-CDS producer instead of these legacy files.
