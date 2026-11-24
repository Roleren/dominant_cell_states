#!/usr/bin/env Rscript
#
# Flags disagreements between the RDG consensus table's headline direction
# and the context-interaction model's properly-pooled global effect, for
# the same gene x design_family x clean-CDS comparison.
#
# This exists because of a real, caught error: VEGFA's consensus-table
# headline (`atlas_clean_cds_pct_change = +46.95%`, framed as induction
# under viral infection) had the OPPOSITE sign from the context-interaction
# table's own pooled `global_axis_delta` for the clean_CDS branch (-0.027)
# once checked directly -- see
# papers/to_be_explored_more_later/VEGFA_uorf_viral_infection/README.md and
# md/dominant_state_pipeline_lessons_2026-09-28.md's gap #1. Nothing in
# either source file flagged the disagreement; it was only found by manual
# cross-checking. This script automates that cross-check so future
# consensus-table rows carry a visible warning instead of requiring the
# same manual dig every time.

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

consensus_file <- result_file("dominant_rdg_consensus", "dominant_rdg_consensus_candidates.csv")
context_file <- result_file("dominant_rdg_context_interactions",
                            "dominant_rdg_context_interaction_branch_context_effects.csv")
output_file <- result_file("dominant_rdg_consensus_context_consistency_audit.csv")
summary_output <- result_file("dominant_rdg_consensus_context_consistency_summary.txt")

if (!file.exists(consensus_file)) stop("Missing ", consensus_file)
if (!file.exists(context_file)) stop("Missing ", context_file)

consensus <- fread(consensus_file)
context <- fread(context_file)

# The context-interaction table's global_axis_delta is constant across all
# context rows for a given gene/design_family/branch (it is the pooled,
# non-context-specific estimate) -- take one row per combination.
context_clean_cds <- unique(
  context[branch_class == "clean_CDS",
          .(gene_symbol, design_family, global_axis_delta, global_axis_delta_p = global_axis_delta_p)]
)
context_clean_cds <- unique(context_clean_cds, by = c("gene_symbol", "design_family"))

audit <- merge(
  consensus[, .(gene_symbol, design_family, consensus_score, consensus_class,
               atlas_clean_cds_pct_change, atlas_direction, n_studies,
               strict_supported_studies, evidence_tags)],
  context_clean_cds,
  by = c("gene_symbol", "design_family"),
  all.x = TRUE
)

negligible_threshold <- 0.02  # |global_axis_delta| below this: treat as "no clear pooled effect" rather than force a sign call
audit[, consensus_sign := sign(atlas_clean_cds_pct_change)]
audit[, pooled_sign := sign(global_axis_delta)]
audit[, agreement_status := fifelse(
  is.na(global_axis_delta), "no_context_interaction_row",
  fifelse(abs(global_axis_delta) < negligible_threshold, "pooled_effect_negligible",
          fifelse(consensus_sign == pooled_sign, "agree", "DISAGREE_SIGN"))
)]

fwrite(audit, output_file)
message("Saved: ", output_file)

n_total <- nrow(audit)
n_disagree <- sum(audit$agreement_status == "DISAGREE_SIGN", na.rm = TRUE)
n_negligible <- sum(audit$agreement_status == "pooled_effect_negligible", na.rm = TRUE)
n_no_row <- sum(audit$agreement_status == "no_context_interaction_row", na.rm = TRUE)
n_agree <- sum(audit$agreement_status == "agree", na.rm = TRUE)

disagreements <- audit[agreement_status == "DISAGREE_SIGN"][order(-consensus_score)]

context_design_families <- sort(unique(context_clean_cds$design_family))
n_no_row_covered_family <- sum(audit$agreement_status == "no_context_interaction_row" &
                               audit$design_family %chin% context_design_families, na.rm = TRUE)

summary_lines <- c(
  sprintf("Consensus vs. context-interaction consistency audit (%s)", Sys.Date()),
  sprintf("Total consensus rows checked: %d", n_total),
  sprintf("  agree: %d", n_agree),
  sprintf("  DISAGREE_SIGN: %d", n_disagree),
  sprintf("  pooled_effect_negligible (|delta| < %.2f): %d", negligible_threshold, n_negligible),
  sprintf("  no_context_interaction_row (not evaluable): %d", n_no_row),
  sprintf(
    "    of which %d are in a design_family the context-interaction model doesn't cover at all yet (%s);",
    n_no_row - n_no_row_covered_family,
    paste(context_design_families, collapse = ", ")
  ),
  sprintf(
    "    only %d are genuinely missing despite being in a covered design_family.",
    n_no_row_covered_family
  ),
  "  NOTE: as of this run, the context-interaction model only covers",
  sprintf("  %s -- 'not evaluable' mostly reflects that coverage gap,", paste(context_design_families, collapse = ", ")),
  "  not a join failure. Extending context-interaction coverage to more",
  "  design families would directly increase how many consensus rows this",
  "  audit can actually check.",
  "",
  "Top disagreements by consensus_score (these are the ones most likely to",
  "be quoted, so most worth checking first):"
)
writeLines(summary_lines, summary_output)
if (nrow(disagreements) > 0) {
  disagree_lines <- capture.output(
    print(disagreements[1:min(20, .N),
                        .(gene_symbol, design_family, consensus_score,
                          atlas_clean_cds_pct_change, global_axis_delta)])
  )
  write(disagree_lines, summary_output, append = TRUE)
}
message("Saved: ", summary_output)

message(sprintf(
  "Consistency audit: %d/%d rows agree, %d DISAGREE, %d negligible pooled effect, %d not evaluable.",
  n_agree, n_total, n_disagree, n_negligible, n_no_row
))
if (nrow(disagreements) > 0) {
  message("Top disagreements (by consensus_score):")
  print(disagreements[1:min(10, .N), .(gene_symbol, design_family, consensus_score,
                                       atlas_clean_cds_pct_change, global_axis_delta)])
}

write_column_dictionary(
  output_file,
  c(
    gene_symbol = "Gene this row audits.",
    design_family = "Design family (e.g. viral_infection) this consensus row is scored under.",
    consensus_score = "The consensus table's own headline score for this row (higher = stronger consensus).",
    consensus_class = "The consensus table's own tier label for this row.",
    atlas_clean_cds_pct_change = "The consensus table's headline percent-change number -- this is the number that disagreed with global_axis_delta for VEGFA; treat with caution when agreement_status != 'agree'.",
    atlas_direction = "The consensus table's own direction label for atlas_clean_cds_pct_change.",
    n_studies = "Number of studies supporting this consensus row (from the consensus table).",
    strict_supported_studies = "Number of studies meeting the consensus table's own stricter support bar.",
    evidence_tags = "The consensus table's own semicolon-separated evidence/caution tags for this row.",
    global_axis_delta = "The context-interaction table's pooled, non-context-specific estimate for the clean_CDS branch of the same gene/design_family -- the more carefully-pooled number to trust when it disagrees with atlas_clean_cds_pct_change.",
    global_axis_delta_p = "P-value for global_axis_delta.",
    consensus_sign = "sign(atlas_clean_cds_pct_change).",
    pooled_sign = "sign(global_axis_delta).",
    agreement_status = "agree / DISAGREE_SIGN / pooled_effect_negligible (|global_axis_delta| below the negligible_threshold in this script) / no_context_interaction_row (the context-interaction model doesn't cover this gene/design_family -- currently only viral_infection is covered at all, so this is the majority case, not a bug)."
  ),
  title = "Consensus-table vs. context-interaction-table sign agreement, per gene x design_family."
)
