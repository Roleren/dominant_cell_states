#!/usr/bin/env Rscript
#
# Loads and validates the curated transcription-vs-translation gene sidecar
# table (gap #6 from md/dominant_state_pipeline_lessons_2026-09-28.md), and
# reports how much of the current 129-gene marker panel it actually
# annotates. The table itself is hand-curated (curated_inputs/
# gene_translation_regulation_notes.csv, at the repo root -- NOT under
# results/curated_inputs/, since that directory is gitignored and this
# table is meant to be a durable, version-controlled reference, not a
# regenerable output); this script's job is only to validate it stays in
# sync with the current gene panel and surface a coverage number, not to
# regenerate its content.

suppressPackageStartupMessages(library(data.table))

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (basename(here) == "dominant_cell_states" &&
        dir.exists(file.path(here, "scripts"))) {
      return(here)
    }
    candidate <- file.path(here, "dominant_cell_states")
    if (dir.exists(file.path(candidate, "scripts"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find dominant_cell_states analysis directory from: ",
       start, call. = FALSE)
}

analysis_dir <- find_analysis_dir()
results_dir <- Sys.getenv("DOMINANT_CELL_STATES_RESULTS_DIR",
                          unset = file.path(analysis_dir, "results"))
result_file <- function(...) file.path(results_dir, ...)
source(file.path(analysis_dir, "scripts", "dominant_column_dictionary.R"))

sidecar_file <- file.path(analysis_dir, "curated_inputs", "gene_translation_regulation_notes.csv")
if (!file.exists(sidecar_file)) {
  stop("Missing curated sidecar table: ", sidecar_file)
}
sidecar <- fread(sidecar_file)

required_cols <- c("gene_symbol", "regulatory_level", "mechanism_note", "source_finding", "confidence")
missing_cols <- setdiff(required_cols, names(sidecar))
if (length(missing_cols) > 0) {
  stop("Sidecar table is missing required columns: ", paste(missing_cols, collapse = ", "))
}

valid_levels <- c(
  "translational_efficiency_switch", "transcriptional_with_translational_component",
  "translational_repression_atypical", "transcript_translation_decoupled",
  "translational_gatekeeping_under_infection", "epitranscriptomic_writer_translation_readout",
  "unstable_reference_gene_not_translation_finding", "unknown"
)
bad_levels <- setdiff(sidecar$regulatory_level, valid_levels)
if (length(bad_levels) > 0) {
  message("NOTE: sidecar table has regulatory_level values not in the documented",
          " vocabulary (not necessarily an error, but worth a second look): ",
          paste(unique(bad_levels), collapse = ", "))
}

source(file.path(analysis_dir, "scripts", "dominant_state_module_definitions.R"))
modules <- load_dominant_state_modules()
panel_genes <- dominant_state_target_genes(modules)

sidecar[, in_current_panel := gene_symbol %chin% panel_genes]
n_annotated_in_panel <- sum(sidecar$in_current_panel)
n_panel_total <- length(panel_genes)

output_file <- result_file("dominant_state_gene_translation_regulation_audit.csv")
fwrite(sidecar, output_file)
message("Saved: ", output_file)

message(sprintf(
  "Gene translation-regulation sidecar: %d rows, %d/%d in the current marker panel (%.1f%% panel coverage).",
  nrow(sidecar), n_annotated_in_panel, n_panel_total,
  100 * n_annotated_in_panel / n_panel_total
))
if (any(!sidecar$in_current_panel)) {
  message("Rows for genes no longer in the current marker panel (check if still relevant): ",
          paste(sidecar$gene_symbol[!sidecar$in_current_panel], collapse = ", "))
}

write_column_dictionary(
  output_file,
  c(
    gene_symbol = "Gene this row documents.",
    regulatory_level = paste(
      "One of:", paste(valid_levels, collapse = ", "),
      "-- see curated_inputs/gene_translation_regulation_notes.csv for the",
      "hand-curated source; extend the vocabulary in this script's",
      "valid_levels if a new category is needed."
    ),
    mechanism_note = "Free-text explanation of why this gene's translation-level behavior does or doesn't match transcript-level literature expectations.",
    source_finding = "Which document/session finding this row is based on.",
    confidence = "Curator's own confidence (high/medium/low) in this row, not a statistical measure.",
    in_current_panel = "TRUE if gene_symbol is currently in the atlas's scored marker panel (dominant_state_target_genes()); FALSE for perturbation-target or removed-reference genes documented here for context but not themselves scored."
  ),
  title = "Curated transcription-vs-translation regulatory notes, cross-checked against the current marker panel."
)
