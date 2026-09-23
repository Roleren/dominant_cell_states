#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(here, "DESCRIPTION")) &&
        dir.exists(file.path(here, "results"))) return(here)
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find dominant_cell_states analysis root.", call. = FALSE)
}

analysis_dir <- find_analysis_dir()
input_file <- file.path(
  analysis_dir, "results",
  "dominant_rdg_eif4g2_trna_glu_lane_collapsed",
  "union_run_metrics.csv"
)
output_dir <- file.path(
  analysis_dir, "results",
  "dominant_rdg_eif4g2_trna_lane_subset_audit"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) {
  stop("Run the exact tRNA-Glu validator first: ", input_file,
       call. = FALSE)
}

runs <- fread(input_file)
required <- c(
  "group", "replicate_label", "fraction", "uorf_counts",
  "clean_cds_counts", "clean_cds_allocation", "uorf_phase0_counts",
  "uorf_phase1_counts", "uorf_phase2_counts", "top_readlength"
)
missing <- setdiff(required, names(runs))
if (length(missing)) {
  stop("Run metrics are missing: ", paste(missing, collapse = ", "),
       call. = FALSE)
}

subsets <- list(
  list(
    subset = "all_lanes",
    read_length_design = "partially_linked_in_replicate2_lanes1_2",
    keep = rep(TRUE, nrow(runs))
  ),
  list(
    subset = "replicate1_all_lanes",
    read_length_design = "condition_matched_within_each_lane",
    keep = runs$replicate_label == 1L
  ),
  list(
    subset = "lane3_both_replicates",
    read_length_design = "30nt_dominant_in_every_arm",
    keep = runs$fraction == "lane3"
  ),
  list(
    subset = "replicate1_lanes1_2",
    read_length_design = "28nt_dominant_in_every_arm",
    keep = runs$replicate_label == 1L & runs$fraction %chin% c("lane1", "lane2")
  ),
  list(
    subset = "replicate2_lanes1_2",
    read_length_design = "case_30nt_control_28nt_partially_confounded",
    keep = runs$replicate_label == 2L & runs$fraction %chin% c("lane1", "lane2")
  )
)

summarise_subset <- function(specification) {
  x <- runs[specification$keep]
  grouped <- x[, .(
    n_runs = .N,
    uorf_counts = sum(uorf_counts),
    clean_cds_counts = sum(clean_cds_counts),
    allocation_min = min(clean_cds_allocation),
    allocation_max = max(clean_cds_allocation),
    phase0_counts = sum(uorf_phase0_counts),
    phase1_counts = sum(uorf_phase1_counts),
    phase2_counts = sum(uorf_phase2_counts),
    dominant_read_lengths = paste(sort(unique(top_readlength)), collapse = ";")
  ), by = group]
  if (!all(c("case", "control") %chin% grouped$group)) {
    stop("Subset does not contain both conditions: ", specification$subset,
         call. = FALSE)
  }
  grouped[, clean_cds_allocation := clean_cds_counts /
            pmax(clean_cds_counts + uorf_counts, 1L)]
  grouped[, uorf_phase0_fraction := phase0_counts /
            pmax(phase0_counts + phase1_counts + phase2_counts, 1L)]
  case <- grouped[group == "case"]
  control <- grouped[group == "control"]
  data.table(
    subset = specification$subset,
    read_length_design = specification$read_length_design,
    n_case_runs = case$n_runs,
    n_control_runs = control$n_runs,
    case_union_uorf_counts = case$uorf_counts,
    control_union_uorf_counts = control$uorf_counts,
    case_clean_cds_counts = case$clean_cds_counts,
    control_clean_cds_counts = control$clean_cds_counts,
    case_clean_cds_allocation = case$clean_cds_allocation,
    control_clean_cds_allocation = control$clean_cds_allocation,
    clean_cds_allocation_delta =
      case$clean_cds_allocation - control$clean_cds_allocation,
    run_level_range_separation =
      control$allocation_min - case$allocation_max,
    case_uorf_phase0_fraction = case$uorf_phase0_fraction,
    control_uorf_phase0_fraction = control$uorf_phase0_fraction,
    case_dominant_read_lengths = case$dominant_read_lengths,
    control_dominant_read_lengths = control$dominant_read_lengths
  )
}

summary <- rbindlist(lapply(subsets, summarise_subset), fill = TRUE)
summary[, direction_matches_full :=
  sign(clean_cds_allocation_delta) ==
    sign(clean_cds_allocation_delta[subset == "all_lanes"])]
summary[, matched_length_subset := !grepl(
  "confound|linked", read_length_design, ignore.case = TRUE
)]

fwrite(
  summary,
  file.path(output_dir, "eif4g2_trna_glu_lane_subset_audit.csv")
)
print(summary)
message("Saved tRNA-Glu lane-subset audit: ", output_dir)
