find_analysis_dir_local <- function(start = getwd()) {
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
    if (identical(parent, here)) {
      stop("Could not find dominant_cell_states from: ", start, call. = FALSE)
    }
    here <- parent
  }
}

expect <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

expect_file <- function(path) {
  expect(file.exists(path), paste0("Missing file: ", path))
}

analysis_dir <- find_analysis_dir_local()
repo_root <- analysis_dir
results_dir <- file.path(analysis_dir, "results")
generated_root_names <- c(
  "^dominant_",
  "^human_dominant_",
  "^required_fst_",
  "^postviral_fatigue_outputs$",
  "^htmlwidget_libs$",
  "^readlength_inputs$"
)
generated_manual_files <- c(
  "manual_translon_candidates.csv",
  "manual_translon_requested_gene_uoorf_candidates.csv",
  "manual_translons_manifest.csv",
  "manual_translons_ranges.rds"
)
is_generated_root_entry <- function(x) {
  is.character(x) &&
    length(x) >= 1L &&
    any(vapply(generated_root_names, grepl, logical(1), x = x[[1]]))
}
is_analysis_root_arg <- function(x) {
  is.character(x) &&
    length(x) == 1L &&
    identical(normalizePath(x, mustWork = FALSE), analysis_dir)
}
is_generated_manual_entry <- function(parts) {
  length(parts) >= 3L &&
    identical(parts[[2]], "manual_translons") &&
    is.character(parts[[3]]) &&
    parts[[3]][[1]] %in% generated_manual_files
}
file.path <- function(..., fsep = .Platform$file.sep) {
  parts <- list(...)
  if (length(parts) >= 2L &&
      is_analysis_root_arg(parts[[1]]) &&
      (is_generated_root_entry(parts[[2]]) ||
       is_generated_manual_entry(parts))) {
    parts[[1]] <- results_dir
  }
  do.call(base::file.path, c(parts, list(fsep = fsep)))
}
script_dir <- file.path(analysis_dir, "scripts")

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("data.table is required for this invariant test.", call. = FALSE)
}

message("Checking script parsing")
scripts <- list.files(script_dir, pattern = "\\.R$", full.names = TRUE)
bad <- lapply(scripts, function(path) {
  tryCatch({
    parse(path)
    NULL
  }, error = function(e) {
    data.frame(file = path, error = conditionMessage(e))
  })
})
bad <- do.call(rbind, bad)
expect(is.null(bad), paste(c("Script parse failures:", capture.output(print(bad))),
                           collapse = "\n"))

message("Checking Shiny app parsing")
rdg_atlas_app <- file.path(
  analysis_dir, "inst", "shiny", "rdg_atlas_browser", "app.R"
)
expect_file(rdg_atlas_app)
tryCatch(
  parse(rdg_atlas_app),
  error = function(e) {
    stop("RDG atlas browser app parse failure: ", conditionMessage(e),
         call. = FALSE)
  }
)
message("Checking Shiny app review-workbench startup")
app_env <- new.env(parent = globalenv())
source(rdg_atlas_app, local = app_env)
expect(exists("review_workbench", envir = app_env, inherits = FALSE),
       "RDG atlas browser should build a review_workbench table at startup.")
app_workbench <- get("review_workbench", envir = app_env)
expect(exists("sdrive_transcript_readiness", envir = app_env,
              inherits = FALSE),
       "RDG atlas browser should load S-drive transcript readiness at startup.")
expect(exists("sdrive_readiness_summary_metrics", envir = app_env,
              inherits = FALSE),
       "RDG atlas browser should load S-drive readiness summary metrics at startup.")
sdrive_app_readiness <- get("sdrive_transcript_readiness", envir = app_env)
sdrive_app_metrics <- get("sdrive_readiness_summary_metrics", envir = app_env)
expect(is.data.frame(sdrive_app_readiness),
       "S-drive transcript readiness should be a table-like object.")
expect(is.data.frame(sdrive_app_metrics),
       "S-drive readiness metrics should be a table-like object.")
workbench_required <- c(
  "queue_source", "gene_symbol", "tx_id", "design_family",
  "target_label_class", "label_gap_class", "label_acquisition_score",
  "protocol_pooling_need",
  "pooling_proof_tier", "pooling_hazard_class", "pooling_hazard_score"
)
expect(all(workbench_required %in% names(app_workbench)),
       "Review workbench is missing label-target or pooling/protocol columns.")
expect(any(app_workbench$queue_source == "Review batch" &
             app_workbench$batch_role == "protocol_pooling_audit" &
             suppressWarnings(as.numeric(app_workbench$protocol_pooling_need)) >=
               0.55,
           na.rm = TRUE),
       "Review workbench should expose B03 protocol/pooling audit rows.")
expect(any(app_workbench$queue_source == "Label gap" &
             nzchar(app_workbench$target_label_class),
           na.rm = TRUE),
       "Review workbench should expose label-gap target classes.")
expect(exists("review_label_draft_cols", envir = app_env, inherits = FALSE),
       "RDG atlas browser should expose review_label_draft_cols at startup.")
app_draft_cols <- get("review_label_draft_cols", envir = app_env)
expect(exists("review_workbench_gene_choices", envir = app_env,
              inherits = FALSE),
       "RDG atlas browser should expose Workbench gene choices at startup.")
app_workbench_gene_choices <- get("review_workbench_gene_choices",
                                  envir = app_env)
expect(length(app_workbench_gene_choices) > 0,
       "Workbench gene selector should have at least one gene choice.")
expect(all(app_workbench_gene_choices %in% app_workbench$gene_symbol),
       "Workbench gene choices should come from review_workbench genes.")
expect(exists("review_workbench_display_columns", envir = app_env,
              inherits = FALSE),
       "RDG atlas browser should expose Workbench display columns at startup.")
app_workbench_display_columns <- get("review_workbench_display_columns",
                                     envir = app_env)
expect(identical(
  app_workbench_display_columns[seq_len(5L)],
  c("Links", "review_destination", "gene_symbol", "tx_id", "design_family")
), "Workbench table should display destination, gene, transcript, and design columns first.")
draft_context_required <- c(
  "target_label_class", "label_gap_class", "label_acquisition_score",
  "protocol_pooling_need", "pooling_proof_tier", "pooling_hazard_class",
  "pooling_hazard_score", "pooling_hazard_sources"
)
expect(all(draft_context_required %in% app_draft_cols),
       "Review label drafts are missing target-label or pooling context columns.")
expect(exists("review_label_open_rows", envir = app_env, inherits = FALSE),
       "RDG atlas browser should expose review_label_open_rows helper.")
expect(exists("choose_next_review_label_row", envir = app_env,
              inherits = FALSE),
       "RDG atlas browser should expose choose_next_review_label_row helper.")
expect(exists("review_workbench_selected_destination", envir = app_env,
              inherits = FALSE),
       "RDG atlas browser should expose selected Workbench row routing helper.")
expect(exists("review_workbench_row_destination", envir = app_env,
              inherits = FALSE),
       "RDG atlas browser should expose vector Workbench row destinations.")
expect(exists("review_workbench_destination_counts", envir = app_env,
              inherits = FALSE),
       "RDG atlas browser should expose Workbench destination counts.")
expect(exists("review_workbench_protocol_audit_mask", envir = app_env,
              inherits = FALSE),
       "RDG atlas browser should expose Workbench protocol-audit mask.")
helper_rows <- data.table::data.table(
  queue_source = "Review batch",
  queue_rank = c(1, 2, 3, 4),
  queue_score = c(4, 3, 2, 1),
  active_review_order = c(1L, 2L, 3L, 4L),
  draft_label_class = c("", "positive", "", ""),
  feedback_label_class = c("", "", "negative", ""),
  is_open = c(TRUE, TRUE, TRUE, FALSE)
)
open_helper_rows <- app_env$review_label_open_rows(helper_rows)
expect(identical(open_helper_rows$active_review_order, 1L),
       "Open review-label helper must skip drafted, ingested, and closed rows.")
label_gap_shape_rows <- data.table::data.table(
  label_gap_rank = c(1, 2, 3),
  label_acquisition_score = c(3, 2, 1),
  active_review_order = c(10L, 11L, 12L),
  draft_label_class = c("", "negative", ""),
  feedback_label_class = c("", "", "confounded")
)
label_gap_open_rows <- app_env$review_label_open_rows(label_gap_shape_rows)
expect(identical(label_gap_open_rows$active_review_order, 10L),
       "Open review-label helper must support Label Dashboard row schema.")
next_helper_rows <- data.table::data.table(
  queue_source = "Review batch",
  queue_rank = c(1, 2, 3),
  queue_score = c(3, 2, 1),
  active_review_order = c(1L, 2L, 3L),
  draft_label_class = "",
  feedback_label_class = "",
  is_open = TRUE
)
expect(identical(
  app_env$choose_next_review_label_row(
    next_helper_rows, current_order = 1L
  )$active_review_order,
  2L
), "Save-and-next helper should advance to the next open task.")
expect(identical(
  app_env$choose_next_review_label_row(
    next_helper_rows, current_order = 3L
  )$active_review_order,
  1L
), "Save-and-next helper should wrap within the current open-task set.")
expect(identical(
  app_env$review_workbench_selected_destination(data.table::data.table(
    active_review_order = 5L,
    inspection_source_row = 12L
  )),
  "review_labels"
), "Selected Workbench active tasks should route to Review Labels.")
expect(identical(
  app_env$review_workbench_selected_destination(data.table::data.table(
    active_review_order = NA_integer_,
    inspection_source_row = 12L
  )),
  "inspection_assets"
), "Selected Workbench inspection-only rows should route to Inspection Assets.")
expect(identical(
  app_env$review_workbench_selected_destination(data.table::data.table(
    active_review_order = NA_integer_,
    inspection_source_row = NA_integer_
  )),
  "selected_card"
), "Selected Workbench rows without active tasks or inspection packets should route to Selected Card.")
route_rows <- data.table::data.table(
  active_review_order = c(5L, NA_integer_, NA_integer_),
  inspection_source_row = c(12L, 12L, NA_integer_),
  protocol_pooling_need = c(0.1, 0.6, NA_real_),
  batch_role = c("", "", "protocol_pooling_audit"),
  recommended_action = c("", "", ""),
  review_reason = c("", "", "")
)
expect(identical(
  app_env$review_workbench_row_destination(route_rows),
  c("review_labels", "inspection_assets", "selected_card")
), "Workbench row destination vector should prioritize labels over inspection over selected card.")
route_counts <- app_env$review_workbench_destination_counts(route_rows)
expect(identical(route_counts$rows, c(1L, 1L, 1L)),
       "Workbench destination counts should count all route classes.")
expect(identical(
  app_env$review_workbench_protocol_audit_mask(route_rows),
  c(FALSE, TRUE, TRUE)
), "Protocol-audit mask should include high protocol need and explicit protocol-pooling audit rows.")
local({
  old_draft_file <- get("review_label_drafts_file", envir = app_env)
  temp_draft_file <- tempfile(fileext = ".csv")
  assign("review_label_drafts_file", temp_draft_file, envir = app_env)
  on.exit({
    assign("review_label_drafts_file", old_draft_file, envir = app_env)
    unlink(temp_draft_file)
  }, add = TRUE)
  make_draft_row <- function(order_id, label) {
    out <- data.table::as.data.table(as.list(stats::setNames(
      rep("", length(app_draft_cols)),
      app_draft_cols
    )))
    out[, global_review_order := as.integer(order_id)]
    out[, review_label_class_preview := label]
    out
  }
  data.table::fwrite(data.table::rbindlist(list(
    make_draft_row(1L, "positive"),
    make_draft_row(2L, "negative")
  ), fill = TRUE), temp_draft_file)
  remaining <- app_env$delete_review_label_draft(1L)
  expect(nrow(remaining) == 1L &&
           identical(remaining$global_review_order, 2L),
         "Deleting one review-label draft must preserve other draft rows.")
})

message("Checking pipeline registry")
source(file.path(script_dir, "run_dominant_cell_state_pipeline.R"))
steps <- pipeline_steps(analysis_dir)
expect(!anyDuplicated(steps$step), "Pipeline step names must be unique.")
expect(all(file.exists(steps$script)),
       paste("Missing pipeline scripts:",
             paste(steps$script[!file.exists(steps$script)], collapse = ", ")))
outputs <- pipeline_step_outputs(analysis_dir)
expect(all(steps$step %in% names(outputs)),
       "Every pipeline step needs declared outputs.")
expect("fst_pages" %in% steps$step,
       "fst_pages step is missing from pipeline registry.")
expect(!any(c("fatigue_fst_ready", "expanded_oxphos_fst_ready",
              "literature_gap_fst_ready") %in% steps$step),
       "Obsolete subset FST-readiness steps should not be registered.")
expect_file(file.path(analysis_dir, "required_fst_pages.csv"))
expect_file(file.path(analysis_dir, "required_fst_pages.txt"))
expect_file(file.path(analysis_dir, "required_fst_page_basenames.txt"))
message("Checking S-drive readiness audit outputs")
sdrive_readiness_dir <- file.path(analysis_dir, "dominant_rdg_sdrive_readiness")
sdrive_status_file <- file.path(
  sdrive_readiness_dir, "sdrive_experiment_status.csv"
)
sdrive_page_file <- file.path(
  sdrive_readiness_dir, "sdrive_fst_page_readiness.csv"
)
sdrive_transcript_file <- file.path(
  sdrive_readiness_dir, "sdrive_transcript_readiness.csv"
)
sdrive_metrics_file <- file.path(
  sdrive_readiness_dir, "sdrive_readiness_summary_metrics.csv"
)
sdrive_class_file <- file.path(
  sdrive_readiness_dir, "sdrive_readiness_class_counts.csv"
)
sdrive_summary_file <- file.path(
  sdrive_readiness_dir, "sdrive_readiness_summary.md"
)
expect_file(sdrive_status_file)
expect_file(sdrive_page_file)
expect_file(sdrive_transcript_file)
expect_file(sdrive_metrics_file)
expect_file(sdrive_class_file)
expect_file(sdrive_summary_file)
sdrive_status <- data.table::fread(sdrive_status_file, nThread = 1)
sdrive_pages <- data.table::fread(sdrive_page_file, nThread = 1)
sdrive_transcripts <- data.table::fread(sdrive_transcript_file, nThread = 1)
sdrive_metrics <- data.table::fread(sdrive_metrics_file, nThread = 1)
expect(all(c("check", "status", "ok", "path", "detail") %in%
             names(sdrive_status)),
       "S-drive experiment status table is missing required columns.")
expect(all(c("fst_basename", "local_file_exists_now",
             "readiness_status", "file_size_mb") %in%
             names(sdrive_pages)),
       "S-drive FST page readiness table is missing required columns.")
expect(nrow(sdrive_pages) > 0,
       "S-drive FST page readiness table should not be empty.")
expect(all(sdrive_pages$local_file_exists_now),
       "All required S-drive FST pages should be present in the current audit.")
expect(all(c(
  "source_row", "gene_symbol", "tx_id", "design_family",
  "sdrive_readiness_class", "review_support_class", "fst_pages_ready",
  "required_fst_pages_missing", "mrna_sequence_available",
  "display_sequence_available", "reference_ready", "all_merged_core_ready"
) %in% names(sdrive_transcripts)),
"S-drive transcript readiness table is missing required columns.")
expect(nrow(sdrive_transcripts) > 0,
       "S-drive transcript readiness table should not be empty.")
expect(all(sdrive_transcripts$sdrive_readiness_class == "sdrive_ready"),
       "All qualified inspection source rows should currently be backed by S-drive data.")
expect(any(sdrive_transcripts$review_support_class == "coverage_reviewable"),
       "S-drive readiness should preserve coverage-reviewable rows.")
expect(any(sdrive_transcripts$review_support_class != "coverage_reviewable"),
       "S-drive readiness should distinguish backing-data availability from coverage support.")
expect(all(sdrive_transcripts$fst_pages_ready),
       "All qualified transcript rows should have required FST pages ready.")
expect(all(c("metric", "value") %in% names(sdrive_metrics)),
       "S-drive readiness metric table is missing required columns.")
required_sdrive_metrics <- c(
  "required_fst_pages_present", "required_fst_pages_missing",
  "qualified_source_rows", "qualified_transcripts",
  "sdrive_ready_source_rows", "coverage_reviewable_source_rows"
)
expect(all(required_sdrive_metrics %in% sdrive_metrics$metric),
       "S-drive readiness metrics are missing required metric rows.")
expect("branch_allocation" %in% steps$step,
       "branch_allocation step is missing from pipeline registry.")
expect("hierarchical_branch_covariates" %in% steps$step,
       "hierarchical_branch_covariates step is missing from pipeline registry.")
expect("graph_geometry_model" %in% steps$step,
       "graph_geometry_model step is missing from pipeline registry.")
expect("mechanistic_motif_model" %in% steps$step,
       "mechanistic_motif_model step is missing from pipeline registry.")
expect("motif_prior_calibration" %in% steps$step,
       "motif_prior_calibration step is missing from pipeline registry.")
expect("joint_branch_allocation" %in% steps$step,
       "joint_branch_allocation step is missing from pipeline registry.")
expect("hierarchical_joint_allocation" %in% steps$step,
       "hierarchical_joint_allocation step is missing from pipeline registry.")
expect("dirichlet_multinomial_allocation" %in% steps$step,
       "dirichlet_multinomial_allocation step is missing from pipeline registry.")
expect("model_agreement" %in% steps$step,
       "model_agreement step is missing from pipeline registry.")
expect("model_agreement_loo" %in% steps$step,
       "model_agreement_loo step is missing from pipeline registry.")
expect("model_agreement_transport" %in% steps$step,
       "model_agreement_transport step is missing from pipeline registry.")
expect("context_interactions" %in% steps$step,
       "context_interactions step is missing from pipeline registry.")
expect("context_partial_pooling" %in% steps$step,
       "context_partial_pooling step is missing from pipeline registry.")
expect("pooling_bias_audit" %in% steps$step,
       "pooling_bias_audit step is missing from pipeline registry.")
expect("therapeutic_hypotheses" %in% steps$step,
       "therapeutic_hypotheses step is missing from pipeline registry.")
expect("therapeutic_perturbation_predictions" %in% steps$step,
       "therapeutic_perturbation_predictions step is missing from pipeline registry.")
expect("allocation_output_disagreements" %in% steps$step,
       "allocation_output_disagreements step is missing from pipeline registry.")
expect("auxiliary_state_screen" %in% steps$step,
       "auxiliary_state_screen step is missing from pipeline registry.")
expect("auxiliary_state_rdg_audit" %in% steps$step,
       "auxiliary_state_rdg_audit step is missing from pipeline registry.")
expect("rdg_consensus" %in% steps$step,
       "rdg_consensus step is missing from pipeline registry.")
expect("rdg_validation_dossiers" %in% steps$step,
       "rdg_validation_dossiers step is missing from pipeline registry.")
expect("rdg_evidence_gaps" %in% steps$step,
       "rdg_evidence_gaps step is missing from pipeline registry.")
expect("merged_replicate_translon_shift" %in% steps$step,
       "merged_replicate_translon_shift step is missing from pipeline registry.")
expect("rdg_review_batches" %in% steps$step,
       "rdg_review_batches step is missing from pipeline registry.")
expect("rdg_review_links" %in% steps$step,
       "rdg_review_links step is missing from pipeline registry.")
expect("rdg_review_feedback" %in% steps$step,
       "rdg_review_feedback step is missing from pipeline registry.")
expect("rdg_validation_labels" %in% steps$step,
       "rdg_validation_labels step is missing from pipeline registry.")
cfg_skip_plots <- parse_pipeline_args(c("--skip-rdg-plots"))
expect(identical(cfg_skip_plots$rdg_plot_mode, "tables"),
       "--skip-rdg-plots should select table-only RDG mode.")
cfg_full_plots <- parse_pipeline_args(c("--rdg-plot-mode", "all"))
expect(identical(cfg_full_plots$rdg_plot_mode, "all"),
       "--rdg-plot-mode all should select full RDG plotting.")
rdg_dependencies <- pipeline_step_dependencies("rdg_graphs", analysis_dir)
rdg_mode_file <- pipeline_runtime_option_file(analysis_dir, "rdg_plot_mode")
expect(normalizePath(rdg_mode_file, mustWork = FALSE) %in%
         normalizePath(rdg_dependencies, mustWork = FALSE),
       "rdg_graphs cache dependencies must include the RDG plot-mode option.")

message("Checking manual translon manifest")
manual_txt <- file.path(analysis_dir, "manual_translons",
                        "translons_to_be_added_to_manual.txt")
manual_manifest <- file.path(analysis_dir, "manual_translons",
                             "manual_translons_manifest.csv")
expect_file(manual_txt)
expect_file(manual_manifest)
manual_dt <- data.table::fread(manual_manifest, nThread = 1)
expect(nrow(manual_dt) > 0, "Manual translon manifest should not be empty.")
if ("manual_selected" %in% names(manual_dt)) {
  expect(any(manual_dt$manual_selected == TRUE, na.rm = TRUE),
         "Manual manifest has no selected manual translons.")
}

message("Checking candidate rank table")
candidate_file <- file.path(analysis_dir, "dominant_next_model",
                            "dominant_next_model_candidate_rankings_v4.csv")
expect_file(candidate_file)
candidate <- data.table::fread(candidate_file, nThread = 1)
expect("candidate_rank_v4" %in% names(candidate),
       "candidate_rank_v4 is missing.")
expect(all(candidate$candidate_rank_v4 == seq_len(nrow(candidate))),
       "candidate_rank_v4 must match row order.")
if ("candidate_rank" %in% names(candidate)) {
  expect(all(candidate$candidate_rank == candidate$candidate_rank_v4),
         "candidate_rank should be the current rank alias.")
}
expect(all(diff(candidate$next_iteration_candidate_score_v4) <= 1e-12),
       "candidate table must be sorted by descending v4 score.")

message("Checking auxiliary dominant-state screen")
aux_dir <- file.path(analysis_dir, "dominant_state_auxiliary_modules")
aux_summary_file <- file.path(aux_dir,
                              "dominant_state_auxiliary_module_summary.csv")
aux_genes_file <- file.path(aux_dir,
                            "dominant_state_auxiliary_module_genes.csv")
aux_scores_file <- file.path(aux_dir,
                             "dominant_state_auxiliary_module_scores.csv")
aux_metrics_file <- file.path(aux_dir,
                              "dominant_state_auxiliary_summary_metrics.csv")
expect_file(aux_summary_file)
expect_file(aux_genes_file)
expect_file(aux_scores_file)
expect_file(aux_metrics_file)
aux_summary <- data.table::fread(aux_summary_file, nThread = 1)
aux_genes <- data.table::fread(aux_genes_file, nThread = 1)
aux_scores <- data.table::fread(aux_scores_file, nThread = 1)
expect(nrow(aux_summary) >= 8,
       "Auxiliary module summary should contain the curated screen modules.")
expect(nrow(aux_scores) >= 8 * 3000,
       "Auxiliary module score table should contain per-run scores.")
expect(all(c("candidate_state", "novelty_score", "recommendation") %in%
             names(aux_summary)),
       "Auxiliary module summary is missing required columns.")
expect(all(is.finite(aux_summary$novelty_score)),
       "Auxiliary novelty scores must be finite.")
expect(all(aux_summary$novelty_score >= 0 & aux_summary$novelty_score <= 1),
       "Auxiliary novelty scores must be bounded in [0, 1].")
expect("available_in_clean_cds" %in% names(aux_genes),
       "Auxiliary marker table must expose current clean-CDS availability.")
pipeline_candidates <- aux_summary[
  recommendation == "candidate_dominant_state_for_pipeline",
  candidate_state
]
expect("Antigen presentation / immune composition" %in% pipeline_candidates,
       "Antigen-presentation auxiliary state should remain a pipeline candidate.")
expect("DNA damage / p53 checkpoint" %in% pipeline_candidates,
       "DNA-damage auxiliary state should remain a pipeline candidate.")

message("Checking auxiliary-state RDG audit")
aux_rdg_dir <- file.path(analysis_dir, "dominant_rdg_auxiliary_state_audit")
aux_rdg_metrics_file <- file.path(
  aux_rdg_dir,
  "dominant_rdg_auxiliary_state_summary_metrics.csv"
)
aux_rdg_gene_file <- file.path(
  aux_rdg_dir,
  "dominant_rdg_auxiliary_state_gene_design_summary.csv"
)
aux_rdg_atlas_file <- file.path(
  aux_rdg_dir,
  "dominant_rdg_auxiliary_state_atlas_annotations.csv"
)
aux_rdg_stratum_file <- file.path(
  aux_rdg_dir,
  "dominant_rdg_auxiliary_state_stratum_scores.csv"
)
aux_rdg_branch_file <- file.path(
  aux_rdg_dir,
  "dominant_rdg_auxiliary_state_branch_rows.csv"
)
expect_file(aux_rdg_metrics_file)
expect_file(aux_rdg_gene_file)
expect_file(aux_rdg_atlas_file)
expect_file(aux_rdg_stratum_file)
expect_file(aux_rdg_branch_file)
aux_rdg_metrics <- data.table::fread(aux_rdg_metrics_file, nThread = 1)
metric_value <- function(name) aux_rdg_metrics[metric == name, value][1]
strata_with_aux <- as.numeric(metric_value("strata_with_aux_case_control"))
branch_aux_rows <- as.numeric(metric_value("branch_aux_rows"))
branch_aux_rows_written <- as.numeric(metric_value("branch_aux_rows_written"))
expect(is.finite(strata_with_aux) && strata_with_aux > 0,
       "Auxiliary RDG audit must reconstruct at least one matched stratum.")
expect(is.finite(branch_aux_rows) && branch_aux_rows > 1000,
       "Auxiliary RDG audit full branch denominator is unexpectedly small.")
expect(is.finite(branch_aux_rows_written) &&
         branch_aux_rows_written > 0 &&
         branch_aux_rows_written <= branch_aux_rows,
       "Auxiliary RDG branch output should be a nonempty compact subset.")
aux_rdg_gene <- data.table::fread(aux_rdg_gene_file, nThread = 1)
expect(nrow(aux_rdg_gene) > 0,
       "Auxiliary RDG gene/design summary is empty.")
gene_required <- c(
  "candidate_state", "gene_symbol", "tx_id", "design_family",
  "max_abs_aux_delta_z", "max_aux_branch_overlap_score",
  "max_weighted_aux_branch_overlap_score", "aux_overlap_class",
  "auxiliary_recommendation_label"
)
expect(all(gene_required %in% names(aux_rdg_gene)),
       "Auxiliary RDG gene/design summary is missing required columns.")
expect(any(aux_rdg_gene$aux_overlap_class ==
             "supported_rdg_in_strong_aux_context", na.rm = TRUE),
       "Auxiliary RDG audit should contain strong context-overlap rows.")
aux_rdg_atlas <- data.table::fread(aux_rdg_atlas_file, nThread = 1)
atlas_required <- c(
  "gene_symbol", "tx_id", "design_family", "aux_top_state",
  "aux_top_abs_delta_z", "aux_top_weighted_overlap_score",
  "aux_audit_priority", "aux_audit_warning"
)
expect(all(atlas_required %in% names(aux_rdg_atlas)),
       "Auxiliary RDG atlas annotations are missing required columns.")
priority <- suppressWarnings(as.numeric(aux_rdg_atlas$aux_audit_priority))
expect(all(is.finite(priority)),
       "Auxiliary RDG atlas audit priority must be finite.")
expect(all(priority >= -1e-12 & priority <= 1 + 1e-12),
       "Auxiliary RDG atlas audit priority must be bounded in [0, 1].")
expect(any(aux_rdg_atlas$aux_top_state %in% c(
  "Antigen presentation / immune composition",
  "Mitochondrial UPR / mitonuclear stress",
  "DNA damage / p53 checkpoint"
), na.rm = TRUE),
       "Auxiliary RDG atlas annotations should include screen-ready states.")

message("Checking atlas browser-priority queue")
atlas_queue_file <- file.path(analysis_dir, "dominant_rdg_atlas",
                              "dominant_rdg_atlas_browser_validation_queue.csv")
expect_file(atlas_queue_file)
atlas_queue <- data.table::fread(atlas_queue_file, nThread = 1)
expect("atlas_browser_review_priority_score" %in% names(atlas_queue),
       "atlas browser priority column is missing.")
expect(!any(is.na(atlas_queue$atlas_browser_review_priority_score)),
       "atlas browser priority score must not contain NA.")
expect(all(is.finite(atlas_queue$atlas_browser_review_priority_score)),
       "atlas browser priority score must be finite.")
expect(all(diff(atlas_queue$atlas_browser_review_priority_score) <= 1e-12),
       "atlas browser queue must be sorted by descending priority.")

atlas_cards_file <- file.path(analysis_dir, "dominant_rdg_atlas",
                              "dominant_rdg_atlas_cards.csv")
expect_file(atlas_cards_file)
atlas_cards <- data.table::fread(atlas_cards_file, nThread = 1)
expect("count_aware_branch_allocation_component" %in% names(atlas_cards),
       "atlas cards must include count-aware branch-allocation component.")
expect("count_aware_top_support_class" %in% names(atlas_cards),
       "atlas cards must include count-aware support class.")
expect("joint_branch_allocation_component" %in% names(atlas_cards),
       "atlas cards must include joint branch-allocation component.")
expect("hierarchical_joint_allocation_component" %in% names(atlas_cards),
       "atlas cards must include hierarchical joint-allocation component.")
expect("dirichlet_multinomial_component" %in% names(atlas_cards),
       "atlas cards must include Dirichlet-multinomial component.")
expect("model_agreement_component" %in% names(atlas_cards),
       "atlas cards must include model-agreement component.")
expect("model_disagreement_component" %in% names(atlas_cards),
       "atlas cards must include model-disagreement component.")
expect(!any(is.na(atlas_cards$count_aware_branch_allocation_component)),
       "count-aware atlas component must not contain NA.")
expect(all(is.finite(atlas_cards$count_aware_branch_allocation_component)),
       "count-aware atlas component must be finite.")
expect(!any(is.na(atlas_cards$joint_branch_allocation_component)),
       "joint branch-allocation atlas component must not contain NA.")
expect(all(is.finite(atlas_cards$joint_branch_allocation_component)),
       "joint branch-allocation atlas component must be finite.")
expect(!any(is.na(atlas_cards$hierarchical_joint_allocation_component)),
       "hierarchical joint-allocation atlas component must not contain NA.")
expect(all(is.finite(atlas_cards$hierarchical_joint_allocation_component)),
       "hierarchical joint-allocation atlas component must be finite.")
expect(!any(is.na(atlas_cards$dirichlet_multinomial_component)),
       "Dirichlet-multinomial atlas component must not contain NA.")
expect(all(is.finite(atlas_cards$dirichlet_multinomial_component)),
       "Dirichlet-multinomial atlas component must be finite.")
expect(!any(is.na(atlas_cards$model_agreement_component)),
       "Model-agreement atlas component must not contain NA.")
expect(all(is.finite(atlas_cards$model_agreement_component)),
       "Model-agreement atlas component must be finite.")
expect(!any(is.na(atlas_cards$model_disagreement_component)),
       "Model-disagreement atlas component must not contain NA.")
expect(all(is.finite(atlas_cards$model_disagreement_component)),
       "Model-disagreement atlas component must be finite.")

message("Checking atlas pooling/protocol proof integration")
pooling_atlas_required <- c(
  "pooling_proof_tier", "pooling_claim_language",
  "pooling_proof_component", "pooling_hazard_class",
  "pooling_hazard_score", "pooling_hazard_sources",
  "pooling_exact_context_fraction", "pooling_joint_strict_fraction",
  "pooling_top_protocol_classes"
)
expect(all(pooling_atlas_required %in% names(atlas_cards)),
       "Atlas cards are missing pooling/protocol proof columns.")
for (column in c("pooling_proof_component", "pooling_hazard_score",
                 "pooling_exact_context_fraction",
                 "pooling_joint_strict_fraction")) {
  value <- suppressWarnings(as.numeric(atlas_cards[[column]]))
  expect(all(is.finite(value)),
         paste0("Atlas ", column, " must be finite."))
  expect(all(value >= -1e-12 & value <= 1 + 1e-12),
         paste0("Atlas ", column, " must be bounded in [0, 1]."))
}
expect(any(atlas_cards$pooling_proof_tier ==
             "A_exact_strict_project_matched_control"),
       "Atlas cards should include exact strict project-matched pooling proof.")
expect(any(atlas_cards$pooling_hazard_class %in%
             c("high_pooling_confounding_hazard",
               "moderate_pooling_confounding_hazard")),
       "Atlas cards should include moderate/high pooling hazard rows.")
expect(any(grepl("pooling_", atlas_cards$warnings, fixed = TRUE)),
       "Atlas warnings should include pooling proof/hazard tags.")

atlas_branch_contexts_file <- file.path(
  analysis_dir,
  "dominant_rdg_atlas",
  "dominant_rdg_atlas_branch_contexts.csv"
)
expect_file(atlas_branch_contexts_file)
atlas_branch_contexts <- data.table::fread(atlas_branch_contexts_file,
                                           nThread = 1)
expect(nrow(atlas_branch_contexts) > 0,
       "Atlas branch-context table should not be empty.")
branch_context_required <- c(
  "gene_symbol", "tx_id", "design_family", "branch_context_rank",
  "branch_context_label", "joint_review_class", "joint_review_priority",
  "joint_top_branch_class", "joint_l1_allocation_shift",
  "pooling_proof_tier", "pooling_hazard_class", "pooling_hazard_score"
)
expect(all(branch_context_required %in% names(atlas_branch_contexts)),
       "Atlas branch-context table is missing required columns.")
expect(!any(is.na(atlas_branch_contexts$branch_context_rank)),
       "Atlas branch-context rank must not contain NA.")
expect(all(atlas_branch_contexts$branch_context_rank >= 1),
       "Atlas branch-context rank must be positive.")
atlas_card_keys <- unique(atlas_cards[, .(gene_symbol, tx_id, design_family)])
missing_context_keys <- atlas_branch_contexts[
  !atlas_card_keys,
  on = .(gene_symbol, tx_id, design_family)
]
expect(nrow(missing_context_keys) == 0,
       "Atlas branch-context rows must all map to atlas cards.")
rank_min <- atlas_branch_contexts[, .(min_rank = min(branch_context_rank)),
                                  by = .(gene_symbol, tx_id, design_family)]
expect(all(rank_min$min_rank == 1),
       "Each atlas branch-context group must start at rank 1.")

message("Checking RDG consensus atlas")
consensus_dir <- file.path(analysis_dir, "dominant_rdg_consensus")
consensus_file <- file.path(
  consensus_dir,
  "dominant_rdg_consensus_candidates.csv"
)
consensus_review_file <- file.path(
  consensus_dir,
  "dominant_rdg_consensus_review_queue.csv"
)
consensus_evidence_file <- file.path(
  consensus_dir,
  "dominant_rdg_consensus_evidence_matrix.csv"
)
consensus_gene_file <- file.path(
  consensus_dir,
  "dominant_rdg_consensus_gene_summary.csv"
)
consensus_metrics_file <- file.path(
  consensus_dir,
  "dominant_rdg_consensus_summary_metrics.csv"
)
expect_file(consensus_file)
expect_file(consensus_review_file)
expect_file(consensus_evidence_file)
expect_file(consensus_gene_file)
expect_file(consensus_metrics_file)
consensus <- data.table::fread(consensus_file, nThread = 1)
consensus_review <- data.table::fread(consensus_review_file, nThread = 1)
consensus_evidence <- data.table::fread(consensus_evidence_file, nThread = 1)
consensus_metrics <- data.table::fread(consensus_metrics_file, nThread = 1)
expect(nrow(consensus) == nrow(atlas_cards),
       "RDG consensus candidate table must have one row per atlas card.")
consensus_required <- c(
  "consensus_rank", "gene_symbol", "tx_id", "design_family",
  "consensus_class", "consensus_score", "evidence_tags",
  "branch_component", "model_component", "robustness_component",
  "pooling_proof_tier", "pooling_hazard_class",
  "pooling_component", "pooling_hazard_component",
  "fragility_component", "consensus_review_question"
)
expect(all(consensus_required %in% names(consensus)),
       "RDG consensus candidate table is missing required columns.")
expect(all(consensus$consensus_rank == seq_len(nrow(consensus))),
       "RDG consensus rank must match sorted row order.")
expect(!any(is.na(consensus$consensus_score)),
       "RDG consensus score must not contain NA.")
expect(all(is.finite(consensus$consensus_score)),
       "RDG consensus score must be finite.")
expect(all(consensus$consensus_score >= 0 &
             consensus$consensus_score <= 1),
       "RDG consensus score must be bounded in [0, 1].")
expect(all(diff(consensus$consensus_score) <= 1e-12),
       "RDG consensus table must be sorted by descending score.")
expect(!any(is.na(consensus$consensus_class) |
              !nzchar(consensus$consensus_class)),
       "RDG consensus class must not be blank.")
for (column in c("pooling_component", "pooling_hazard_component")) {
  value <- suppressWarnings(as.numeric(consensus[[column]]))
  expect(all(is.finite(value)),
         paste0("RDG consensus ", column, " must be finite."))
  expect(all(value >= -1e-12 & value <= 1 + 1e-12),
         paste0("RDG consensus ", column, " must be bounded in [0, 1]."))
}
expect(any(grepl("exact_project_matched", consensus$evidence_tags,
                 fixed = TRUE)),
       "RDG consensus tags should expose exact project-matched rows.")
expect(any(grepl("pooling_hazard", consensus$evidence_tags, fixed = TRUE)),
       "RDG consensus tags should expose pooling-hazard rows.")
expect(any(consensus$consensus_class ==
             "A_browser_validated_consensus"),
       "RDG consensus should preserve at least one browser-validated anchor.")
expect(nrow(consensus_review) > 0,
       "RDG consensus review queue should not be empty.")
expect("review_rank" %in% names(consensus_review),
       "RDG consensus review queue must include review_rank.")
expect(all(consensus_review$review_rank == seq_len(nrow(consensus_review))),
       "RDG consensus review rank must match row order.")
evidence_required <- c(
  "gene_symbol", "tx_id", "design_family", "component",
  "component_label", "component_score"
)
expect(all(evidence_required %in% names(consensus_evidence)),
       "RDG consensus evidence matrix is missing required columns.")
expect(all(consensus_evidence$component_score >= 0 &
             consensus_evidence$component_score <= 1),
       "RDG consensus evidence component scores must be bounded in [0, 1].")
expect(any(consensus_evidence$component == "pooling_component"),
       "RDG consensus evidence matrix must include pooling proof.")
expect(any(consensus_evidence$component == "pooling_hazard_component"),
       "RDG consensus evidence matrix must include pooling hazard.")
metric_value <- function(name) consensus_metrics[metric == name, value][1]
expect(as.numeric(metric_value("atlas_rows")) == nrow(atlas_cards),
       "RDG consensus metrics atlas_rows must match atlas cards.")
expect(as.numeric(metric_value("review_queue_rows")) == nrow(consensus_review),
       "RDG consensus metrics review_queue_rows must match review queue.")
expect(as.numeric(metric_value("pooling_exact_project_matched_rows")) > 0,
       "RDG consensus metrics should count exact project-matched pooling rows.")
expect(as.numeric(metric_value("pooling_hazard_ge_0_40_rows")) > 0,
       "RDG consensus metrics should count pooling-hazard rows.")

message("Checking RDG validation dossiers")
validation_dir <- file.path(analysis_dir, "dominant_rdg_validation_dossiers")
validation_candidate_file <- file.path(
  validation_dir,
  "dominant_rdg_validation_dossier_candidates.csv"
)
validation_context_file <- file.path(
  validation_dir,
  "dominant_rdg_validation_dossier_branch_contexts.csv"
)
validation_template_file <- file.path(
  validation_dir,
  "dominant_rdg_validation_dossier_manual_review_template.csv"
)
validation_gene_file <- file.path(
  validation_dir,
  "dominant_rdg_validation_dossier_gene_summary.csv"
)
validation_metrics_file <- file.path(
  validation_dir,
  "dominant_rdg_validation_dossier_summary_metrics.csv"
)
validation_report_file <- file.path(
  validation_dir,
  "dominant_rdg_validation_dossier_report.md"
)
expect_file(validation_candidate_file)
expect_file(validation_context_file)
expect_file(validation_template_file)
expect_file(validation_gene_file)
expect_file(validation_metrics_file)
expect_file(validation_report_file)
validation_candidate <- data.table::fread(validation_candidate_file,
                                          nThread = 1)
validation_context <- data.table::fread(validation_context_file, nThread = 1)
validation_template <- data.table::fread(validation_template_file, nThread = 1)
validation_metrics <- data.table::fread(validation_metrics_file, nThread = 1)
expect(nrow(validation_candidate) == nrow(consensus),
       "Validation dossier candidate table must have one row per consensus candidate.")
validation_required <- c(
  "review_packet_id", "selected_for_primary_review",
  "review_priority_tier", "validation_priority_score",
  "consensus_rank", "gene_symbol", "tx_id", "design_family",
  "expected_browser_pattern", "consensus_review_question",
  "pooling_proof_tier", "pooling_claim_language",
  "pooling_hazard_class", "pooling_hazard_score",
  "pooling_hazard_sources"
)
expect(all(validation_required %in% names(validation_candidate)),
       "Validation dossier candidate table is missing required columns.")
validation_hazard <- suppressWarnings(
  as.numeric(validation_candidate$pooling_hazard_score)
)
expect(all(is.finite(validation_hazard)),
       "Validation dossier pooling hazard scores must be finite.")
expect(all(validation_hazard >= -1e-12 & validation_hazard <= 1 + 1e-12),
       "Validation dossier pooling hazard scores must be bounded in [0, 1].")
expect(!anyDuplicated(validation_candidate$review_packet_id),
       "Validation dossier packet IDs must be unique.")
expect(all(is.finite(validation_candidate$validation_priority_score)),
       "Validation priority scores must be finite.")
expect(all(validation_candidate$validation_priority_score >= 0 &
             validation_candidate$validation_priority_score <= 1),
       "Validation priority scores must be bounded in [0, 1].")
primary_rows <- validation_candidate[
  selected_for_primary_review == TRUE
]
expect(nrow(primary_rows) > 0,
       "Validation dossier must select at least one primary review row.")
expect(any(primary_rows$review_priority_tier ==
             "anchor_positive_control"),
       "Validation dossier should preserve an anchor positive-control row.")
expect(nrow(validation_template) == nrow(primary_rows),
       "Manual review template must have one row per primary review row.")
expect(all(validation_template$primary_review_rank == seq_len(
  nrow(validation_template)
)),
       "Manual review template primary ranks must be continuous.")
expect(all(validation_template$review_packet_id %in%
             primary_rows$review_packet_id),
       "Manual review template rows must map to selected primary packets.")
expect(nrow(validation_context) > 0,
       "Validation branch-context packet table should not be empty.")
context_required <- c(
  "review_packet_id", "branch_context_rank", "branch_context_label",
  "context_count_gate", "replicate_gate", "expected_context_pattern"
)
expect(all(context_required %in% names(validation_context)),
       "Validation branch-context packet table is missing required columns.")
expect(all(validation_context$review_packet_id %in%
             primary_rows$review_packet_id),
       "Validation branch-context rows must map to selected primary packets.")
expect(all(validation_context$branch_context_rank >= 1 &
             validation_context$branch_context_rank <= 3),
       "Validation branch-context packets should keep only top contexts.")
metric_value_validation <- function(name) {
  validation_metrics[metric == name, value][1]
}
expect(as.numeric(metric_value_validation("candidate_rows")) ==
         nrow(validation_candidate),
       "Validation metrics candidate_rows must match candidate table.")
expect(as.numeric(metric_value_validation("manual_template_rows")) ==
         nrow(validation_template),
       "Validation metrics manual_template_rows must match template table.")
expect(as.numeric(metric_value_validation("branch_context_rows")) ==
         nrow(validation_context),
       "Validation metrics branch_context_rows must match context table.")
expect(as.numeric(metric_value_validation("pooling_exact_primary_rows")) > 0,
       "Validation metrics should count exact pooling-proof primary rows.")
expect(as.numeric(metric_value_validation("pooling_protocol_hazard_primary_rows")) > 0,
       "Validation metrics should count pooling/protocol-hazard primary rows.")

message("Checking RDG evidence gaps and active-review batches")
evidence_gap_dir <- file.path(analysis_dir, "dominant_rdg_evidence_gaps")
evidence_priority_file <- file.path(
  evidence_gap_dir,
  "dominant_rdg_evidence_gap_priorities.csv"
)
evidence_metrics_file <- file.path(
  evidence_gap_dir,
  "dominant_rdg_evidence_gap_summary_metrics.csv"
)
evidence_action_file <- file.path(
  evidence_gap_dir,
  "dominant_rdg_evidence_gap_action_summary.csv"
)
expect_file(evidence_priority_file)
expect_file(evidence_metrics_file)
expect_file(evidence_action_file)
evidence_priority <- data.table::fread(evidence_priority_file, nThread = 1)
evidence_metrics <- data.table::fread(evidence_metrics_file, nThread = 1)
evidence_action <- data.table::fread(evidence_action_file, nThread = 1)
evidence_required <- c(
  "acquisition_rank", "action_label", "claim_stage",
  "protocol_pooling_need", "review_reason",
  "pooling_proof_tier", "pooling_claim_language",
  "pooling_hazard_class", "pooling_hazard_score",
  "pooling_hazard_sources", "pooling_top_protocol_classes"
)
expect(all(evidence_required %in% names(evidence_priority)),
       "Evidence-gap priority table is missing pooling/protocol columns.")
expect(all(evidence_priority$acquisition_rank == seq_len(nrow(evidence_priority))),
       "Evidence-gap acquisition rank must match row order.")
for (column in c("protocol_pooling_need", "pooling_hazard_score")) {
  values <- suppressWarnings(as.numeric(evidence_priority[[column]]))
  expect(all(is.finite(values)),
         paste("Evidence-gap column must be finite:", column))
  expect(all(values >= -1e-12 & values <= 1 + 1e-12),
         paste("Evidence-gap column must be bounded in [0, 1]:", column))
}
expect(any(evidence_priority$protocol_pooling_need >= 0.55, na.rm = TRUE),
       "Evidence-gap priorities should expose rows needing pooling/protocol audit.")
expect("protocol_pooling_audit" %in% evidence_action$action_label,
       "Evidence-gap action summary should include protocol_pooling_audit.")
metric_value_evidence <- function(name) {
  evidence_metrics[metric == name, value][1]
}
expect(as.numeric(metric_value_evidence("candidate_rows")) ==
         nrow(evidence_priority),
       "Evidence-gap metrics candidate_rows must match priority table.")
expect(as.numeric(metric_value_evidence("protocol_pooling_need_rows")) > 0,
       "Evidence-gap metrics should count pooling/protocol-need rows.")
expect(as.numeric(metric_value_evidence("pooling_hazard_rows")) > 0,
       "Evidence-gap metrics should count pooling-hazard rows.")

message("Checking merged-replicate translon-shift algorithm")
merged_shift_dir <- file.path(
  analysis_dir,
  "dominant_rdg_merged_replicate_translon_shift"
)
merged_shift_rows_file <- file.path(
  merged_shift_dir,
  "merged_replicate_translon_shift_rows.csv"
)
merged_shift_queue_file <- file.path(
  merged_shift_dir,
  "merged_replicate_translon_shift_review_queue.csv"
)
merged_shift_gene_file <- file.path(
  merged_shift_dir,
  "merged_replicate_translon_shift_gene_design_summary.csv"
)
merged_shift_candidate_file <- file.path(
  merged_shift_dir,
  "merged_replicate_translon_shift_candidate_queue.csv"
)
merged_shift_pooling_file <- file.path(
  merged_shift_dir,
  "merged_replicate_translon_shift_pooling_summary.csv"
)
merged_shift_metrics_file <- file.path(
  merged_shift_dir,
  "merged_replicate_translon_shift_summary_metrics.csv"
)
expect_file(merged_shift_rows_file)
expect_file(merged_shift_queue_file)
expect_file(merged_shift_gene_file)
expect_file(merged_shift_candidate_file)
expect_file(merged_shift_pooling_file)
expect_file(merged_shift_metrics_file)
merged_shift_rows <- data.table::fread(merged_shift_rows_file, nThread = 1)
merged_shift_queue <- data.table::fread(merged_shift_queue_file, nThread = 1)
merged_shift_gene <- data.table::fread(merged_shift_gene_file, nThread = 1)
merged_shift_candidate <- data.table::fread(merged_shift_candidate_file,
                                            nThread = 1)
merged_shift_pooling <- data.table::fread(merged_shift_pooling_file,
                                          nThread = 1)
merged_shift_metrics <- data.table::fread(merged_shift_metrics_file,
                                          nThread = 1)
expect(nrow(merged_shift_rows) > 0,
       "Merged-replicate translon-shift rows should not be empty.")
expect(nrow(merged_shift_queue) > 0,
       "Merged-replicate translon-shift review queue should not be empty.")
expect(nrow(merged_shift_gene) > 0,
       "Merged-replicate translon-shift gene summary should not be empty.")
expect(nrow(merged_shift_candidate) > 0,
       "Merged-replicate translon-shift candidate queue should not be empty.")
expect(nrow(merged_shift_pooling) > 0,
       "Merged-replicate translon-shift pooling summary should not be empty.")
merged_shift_required <- c(
  "merged_translon_shift_rank", "gene_symbol", "tx_id", "design_family",
  "dominant_uorf_branch_class", "dominant_uorf_delta",
  "clean_uorf_opposition", "clean_uorf_shift_class",
  "merged_replicate_support_class", "pooling_scope_class",
  "inspection_possible", "discovery_priority_score",
  "claim_readiness_class", "missing_evidence_flags"
)
expect(all(merged_shift_required %in% names(merged_shift_rows)),
       "Merged-replicate translon-shift rows are missing required columns.")
expect(all(merged_shift_required %in% names(merged_shift_queue)),
       "Merged-replicate translon-shift queue is missing required columns.")
candidate_required <- c(
  "candidate_queue_rank", "gene_symbol", "tx_id", "design_family",
  "clean_uorf_shift_class", "candidate_rows", "candidate_studies",
  "exact_context_rows", "opposite_direction_rows",
  "candidate_priority_score", "recurrence_class",
  "pooling_benefit_class", "pooling_loss_class",
  "candidate_interpretation_class", "next_review_action",
  "clean_cds_output_allocation_alignment",
  "allocation_output_disagreement_rows"
)
expect(all(candidate_required %in% names(merged_shift_candidate)),
       "Merged-replicate translon-shift candidate queue is missing required columns.")
pooling_required <- c(
  "design_family", "pooling_scope_class", "clean_uorf_shift_class",
  "claim_readiness_class", "rows", "review_queue_rows",
  "clean_uorf_opposition_rows", "exact_merged_rows",
  "protocol_or_pooling_audit_rows"
)
expect(all(pooling_required %in% names(merged_shift_pooling)),
       "Merged-replicate translon-shift pooling summary is missing required columns.")
expect(all(merged_shift_queue$review_queue_rank ==
             seq_len(nrow(merged_shift_queue))),
       "Merged-replicate translon-shift queue rank must match row order.")
expect(all(merged_shift_candidate$candidate_queue_rank ==
             seq_len(nrow(merged_shift_candidate))),
       "Merged-replicate translon-shift candidate rank must match row order.")
score_values <- suppressWarnings(
  as.numeric(merged_shift_rows$discovery_priority_score)
)
expect(all(is.finite(score_values)),
       "Merged-replicate translon-shift scores must be finite.")
expect(all(score_values >= -1e-12 & score_values <= 1 + 1e-12),
       "Merged-replicate translon-shift scores must be bounded in [0, 1].")
expect(all(merged_shift_queue$inspection_possible == TRUE),
       "Merged-replicate translon-shift review queue should be inspectable.")
expect(any(merged_shift_queue$clean_uorf_shift_class ==
             "clean_CDS_up_uORF_down", na.rm = TRUE),
       "Merged-replicate translon-shift queue should include CDS-up/uORF-down rows.")
expect(any(merged_shift_queue$clean_uorf_shift_class ==
             "clean_CDS_down_uORF_up", na.rm = TRUE),
       "Merged-replicate translon-shift queue should include CDS-down/uORF-up rows.")
expect(any(merged_shift_queue$design_family == "viral_infection" &
             merged_shift_queue$gene_symbol %in%
               c("ATF4", "DDIT3", "PPP1R15A"), na.rm = TRUE),
       "Merged-replicate translon-shift queue should retain viral ISR anchors.")
candidate_scores <- suppressWarnings(
  as.numeric(merged_shift_candidate$candidate_priority_score)
)
expect(all(is.finite(candidate_scores)),
       "Merged-replicate translon-shift candidate scores must be finite.")
expect(all(candidate_scores >= -1e-12 & candidate_scores <= 1 + 1e-12),
       "Merged-replicate translon-shift candidate scores must be bounded in [0, 1].")
expect(any(merged_shift_candidate$candidate_rows > 1, na.rm = TRUE),
       "Merged-replicate translon-shift candidate queue should preserve recurrent rows.")
expect(any(merged_shift_candidate$design_family == "viral_infection" &
             merged_shift_candidate$gene_symbol %in%
               c("ATF4", "DDIT3", "PPP1R15A"), na.rm = TRUE),
       "Merged-replicate translon-shift candidate queue should retain viral ISR anchors.")
bad_claim_rows <- merged_shift_queue[
  claim_readiness_class == "claim_candidate_exact_merged_review_ready" &
    grepl(
      "missing_browser_label|missing_static_coverage_support|protocol_pooling_audit_needed|missing_sdrive_backing_data",
      missing_evidence_flags
    )
]
expect(nrow(bad_claim_rows) == 0,
       "Claim-ready merged-translon rows must not carry missing label/coverage/protocol flags.")
metric_value_shift <- function(name) {
  merged_shift_metrics[metric == name, value][1]
}
expect(as.numeric(metric_value_shift("review_queue_rows")) ==
         nrow(merged_shift_queue),
       "Merged-replicate translon-shift metrics review_queue_rows must match queue table.")
expect(as.numeric(metric_value_shift("candidate_queue_rows")) ==
         nrow(merged_shift_candidate),
       "Merged-replicate translon-shift metrics candidate_queue_rows must match candidate table.")
expect(as.numeric(metric_value_shift("clean_uorf_opposition_rows")) > 0,
       "Merged-replicate translon-shift metrics should count clean/uORF opposition rows.")

review_batch_dir <- file.path(analysis_dir, "dominant_rdg_review_batches")
review_batch_file <- file.path(
  review_batch_dir,
  "dominant_rdg_review_batches.csv"
)
review_batch_template_file <- file.path(
  review_batch_dir,
  "dominant_rdg_review_batch_template.csv"
)
review_batch_metrics_file <- file.path(
  review_batch_dir,
  "dominant_rdg_review_batch_metrics.csv"
)
expect_file(review_batch_file)
expect_file(review_batch_template_file)
expect_file(review_batch_metrics_file)
review_batch <- data.table::fread(review_batch_file, nThread = 1)
review_batch_template <- data.table::fread(review_batch_template_file,
                                           nThread = 1)
review_batch_metrics <- data.table::fread(review_batch_metrics_file,
                                          nThread = 1)
review_required <- c(
  "global_review_order", "batch_id", "batch_role",
  "protocol_pooling_need", "pooling_proof_tier",
  "pooling_hazard_class", "pooling_hazard_score",
  "pooling_hazard_sources", "review_instruction"
)
expect(all(review_required %in% names(review_batch)),
       "Review-batch table is missing pooling/protocol review columns.")
expect(all(review_required %in% names(review_batch_template)),
       "Review-batch template is missing pooling/protocol review columns.")
expect(any(review_batch$batch_role == "protocol_pooling_audit"),
       "Active-review batches should include a protocol/pooling audit batch.")
expect(all(review_batch$global_review_order == seq_len(nrow(review_batch))),
       "Review-batch global order must match row order.")
review_protocol <- review_batch[batch_role == "protocol_pooling_audit"]
expect(all(review_protocol$protocol_pooling_need >= 0.55 |
             suppressWarnings(as.numeric(review_protocol$pooling_hazard_score)) >= 0.40,
           na.rm = TRUE),
       "Protocol/pooling audit rows should have explicit pooling need or hazard.")
metric_value_review <- function(name) {
  review_batch_metrics[metric == name, value][1]
}
expect(as.numeric(metric_value_review("protocol_pooling_audit_rows")) > 0,
       "Review-batch metrics should count protocol/pooling audit rows.")

message("Checking review-feedback context columns")
review_feedback_dir <- file.path(analysis_dir, "dominant_rdg_review_feedback")
review_feedback_label_file <- file.path(
  review_feedback_dir,
  "dominant_rdg_review_feedback_labels.csv"
)
review_feedback_training_file <- file.path(
  review_feedback_dir,
  "dominant_rdg_review_feedback_training_set.csv"
)
expect_file(review_feedback_label_file)
expect_file(review_feedback_training_file)
review_feedback_labels <- data.table::fread(review_feedback_label_file,
                                            nThread = 1)
review_feedback_training <- data.table::fread(review_feedback_training_file,
                                              nThread = 1)
review_feedback_context_required <- c(
  "target_label_class", "label_gap_class", "label_acquisition_score",
  "protocol_pooling_need", "pooling_proof_tier", "pooling_hazard_class",
  "pooling_hazard_score", "pooling_hazard_sources"
)
expect(all(review_feedback_context_required %in%
             names(review_feedback_labels)),
       "Feedback labels should preserve review target-label/pooling context.")
expect(all(review_feedback_context_required %in%
             names(review_feedback_training)),
       "Feedback training set should preserve review target-label/pooling context.")

message("Checking grouped branch denominator invariants")
grouped_contrast_file <- file.path(
  analysis_dir,
  "dominant_rdg_grouped_branch_usage",
  "dominant_rdg_grouped_branch_usage_contrast_table.csv"
)
expect_file(grouped_contrast_file)
grouped <- data.table::fread(grouped_contrast_file, nThread = 1)
expect(nrow(grouped) > 0, "Grouped branch contrast table is empty.")
expect(all(grouped$case_sum_raw_counts + grouped$case_other_branch_counts ==
             grouped$case_branch_total_counts),
       "Case feature + other counts must equal case branch total.")
expect(all(grouped$control_sum_raw_counts +
             grouped$control_other_branch_counts ==
             grouped$control_branch_total_counts),
       "Control feature + other counts must equal control branch total.")

message("Checking count-aware branch-allocation outputs")
allocation_file <- file.path(
  analysis_dir,
  "dominant_rdg_branch_allocation",
  "dominant_rdg_branch_allocation_gene_design_summary.csv"
)
allocation_metrics_file <- file.path(
  analysis_dir,
  "dominant_rdg_branch_allocation",
  "dominant_rdg_branch_allocation_summary_metrics.csv"
)
expect_file(allocation_file)
expect_file(allocation_metrics_file)
allocation <- data.table::fread(allocation_file, nThread = 1)
expect(nrow(allocation) > 0,
       "Count-aware branch-allocation gene-design summary is empty.")
expect(!any(is.na(allocation$count_aware_top_score)),
       "Count-aware branch-allocation scores must not be NA.")
expect(all(is.finite(allocation$count_aware_top_score)),
       "Count-aware branch-allocation scores must be finite.")
expect(all(diff(allocation$count_aware_top_score) <= 1e-12),
       "Count-aware branch-allocation summary must be score sorted.")

message("Checking hierarchical count-aware integration")
hierarchical_file <- file.path(
  analysis_dir,
  "dominant_rdg_hierarchical",
  "dominant_rdg_hierarchical_gene_design_effects.csv"
)
expect_file(hierarchical_file)
hierarchical <- data.table::fread(hierarchical_file, nThread = 1)
expect("hierarchical_count_aware_branch_component" %in% names(hierarchical),
       "hierarchical output must include count-aware branch component.")
expect("hierarchical_mechanistic_support_class" %in% names(hierarchical),
       "hierarchical output must include mechanistic support class.")
expect(!any(is.na(hierarchical$hierarchical_count_aware_branch_component)),
       "hierarchical count-aware branch component must not contain NA.")
expect(all(is.finite(hierarchical$hierarchical_count_aware_branch_component)),
       "hierarchical count-aware branch component must be finite.")
expect(!any(is.na(hierarchical$hierarchical_mechanistic_support_class) |
              !nzchar(hierarchical$hierarchical_mechanistic_support_class)),
       "hierarchical mechanistic support class must not be blank.")

message("Checking hierarchical branch-covariate residual model")
branch_cov_dir <- file.path(
  analysis_dir,
  "dominant_rdg_hierarchical_branch_covariates"
)
branch_cov_summary_file <- file.path(
  branch_cov_dir,
  "dominant_rdg_hierarchical_branch_covariate_model_summary.csv"
)
branch_cov_gene_file <- file.path(
  branch_cov_dir,
  "dominant_rdg_hierarchical_branch_covariate_model_gene_design.csv"
)
branch_cov_training_file <- file.path(
  branch_cov_dir,
  "dominant_rdg_hierarchical_branch_covariates_training_covariates.csv"
)
expect_file(branch_cov_summary_file)
expect_file(branch_cov_gene_file)
expect_file(branch_cov_training_file)
branch_cov_summary <- data.table::fread(branch_cov_summary_file, nThread = 1)
expect("all" %in% branch_cov_summary$scope,
       "Branch-covariate summary must include all-scope metrics.")
branch_required_cols <- c(
  "branch_model_rmse_log2fc",
  "metadata_model_rmse_log2fc",
  "branch_vs_metadata_rmse_delta",
  "branch_vs_metadata_mae_delta"
)
expect(all(branch_required_cols %in% names(branch_cov_summary)),
       "Branch-covariate summary is missing required metric columns.")
branch_all <- branch_cov_summary[scope == "all"][1]
expect(all(is.finite(unlist(branch_all[, ..branch_required_cols]))),
       "All-scope branch-covariate metrics must be finite.")
branch_cov_training <- data.table::fread(branch_cov_training_file, nThread = 1)
training_required_cols <- c(
  "heldout_study",
  "train_eval_rows_clean_CDS",
  "train_eval_rows_leader_uORF",
  "train_eval_rows_overlapping_uORF",
  "branch_covariate_total_evaluable_rows"
)
expect(all(training_required_cols %in% names(branch_cov_training)),
       "Branch training covariates are missing leakage-control columns.")

message("Checking graph-geometry generalization model")
graph_geometry_dir <- file.path(
  analysis_dir,
  "dominant_rdg_graph_geometry_model"
)
graph_geometry_summary_file <- file.path(
  graph_geometry_dir,
  "dominant_rdg_graph_geometry_model_summary.csv"
)
graph_geometry_gene_file <- file.path(
  graph_geometry_dir,
  "dominant_rdg_graph_geometry_model_gene_design.csv"
)
graph_geometry_covariate_file <- file.path(
  graph_geometry_dir,
  "dominant_rdg_graph_geometry_model_covariates.csv"
)
expect_file(graph_geometry_summary_file)
expect_file(graph_geometry_gene_file)
expect_file(graph_geometry_covariate_file)
graph_geometry_summary <- data.table::fread(graph_geometry_summary_file,
                                            nThread = 1)
expect("all" %in% graph_geometry_summary$scope,
       "Graph-geometry model summary must include all-scope metrics.")
graph_required_cols <- c(
  "context_nogene_rmse_log2fc",
  "graph_geometry_rmse_log2fc",
  "graph_vs_context_rmse_delta",
  "graph_vs_context_mae_delta"
)
expect(all(graph_required_cols %in% names(graph_geometry_summary)),
       "Graph-geometry summary is missing required metric columns.")
graph_all <- graph_geometry_summary[scope == "all"][1]
expect(all(is.finite(unlist(graph_all[, ..graph_required_cols]))),
       "All-scope graph-geometry metrics must be finite.")
graph_geometry_covariates <- data.table::fread(graph_geometry_covariate_file,
                                               nThread = 1)
graph_cov_required_cols <- c(
  "gene_symbol", "tx_id", "heldout_study", "n_uorfs",
  "n_overlapping_uorfs", "overlap_class"
)
expect(all(graph_cov_required_cols %in% names(graph_geometry_covariates)),
       "Graph-geometry covariates are missing core architecture columns.")

message("Checking mechanistic RDG motif model")
motif_dir <- file.path(
  analysis_dir,
  "dominant_rdg_mechanistic_motif_model"
)
motif_summary_file <- file.path(
  motif_dir,
  "dominant_rdg_mechanistic_motif_model_summary.csv"
)
motif_gene_file <- file.path(
  motif_dir,
  "dominant_rdg_mechanistic_motif_model_gene_design.csv"
)
motif_scores_file <- file.path(
  motif_dir,
  "dominant_rdg_mechanistic_motif_scores.csv"
)
expect_file(motif_summary_file)
expect_file(motif_gene_file)
expect_file(motif_scores_file)
motif_summary <- data.table::fread(motif_summary_file, nThread = 1)
expect("all" %in% motif_summary$scope,
       "Mechanistic motif model summary must include all-scope metrics.")
motif_required_cols <- c(
  "context_nogene_rmse_log2fc",
  "motif_model_rmse_log2fc",
  "motif_design_model_rmse_log2fc",
  "motif_model_vs_context_rmse_delta",
  "best_motif_model_prefix"
)
expect(all(motif_required_cols %in% names(motif_summary)),
       "Mechanistic motif summary is missing required metric columns.")
motif_all <- motif_summary[scope == "all"][1]
finite_motif_cols <- setdiff(motif_required_cols, "best_motif_model_prefix")
expect(all(is.finite(unlist(motif_all[, ..finite_motif_cols]))),
       "All-scope mechanistic motif metrics must be finite.")
motif_scores <- data.table::fread(motif_scores_file, nThread = 1)
motif_score_required_cols <- c(
  "gene_symbol", "tx_id", "inhibitory_overlap_burden",
  "reinitiation_distance_opportunity", "dense_uorf_collision_burden",
  "strong_start_leakage_burden", "downstream_rescue_opportunity",
  "motif_complexity_index"
)
expect(all(motif_score_required_cols %in% names(motif_scores)),
       "Mechanistic motif scores are missing core score columns.")
bounded_cols <- setdiff(motif_score_required_cols, c("gene_symbol", "tx_id"))
for (column in bounded_cols) {
  values <- suppressWarnings(as.numeric(motif_scores[[column]]))
  expect(all(is.finite(values)),
         paste("Motif score must be finite:", column))
  expect(all(values >= -1e-12 & values <= 1 + 1e-12),
         paste("Motif score must be bounded in [0, 1]:", column))
}

message("Checking motif-prior branch calibration")
motif_prior_dir <- file.path(
  analysis_dir,
  "dominant_rdg_motif_prior_calibration"
)
motif_prior_summary_file <- file.path(
  motif_prior_dir,
  "dominant_rdg_motif_prior_calibration_summary.csv"
)
motif_prior_gene_file <- file.path(
  motif_prior_dir,
  "dominant_rdg_motif_prior_gene_design_calibration.csv"
)
motif_prior_review_file <- file.path(
  motif_prior_dir,
  "dominant_rdg_motif_prior_review_queue.csv"
)
motif_prior_viral_file <- file.path(
  motif_prior_dir,
  "dominant_rdg_motif_prior_viral_review.csv"
)
expect_file(motif_prior_summary_file)
expect_file(motif_prior_gene_file)
expect_file(motif_prior_review_file)
expect_file(motif_prior_viral_file)
motif_prior_summary <- data.table::fread(motif_prior_summary_file,
                                         nThread = 1)
expect("all" %in% motif_prior_summary$scope,
       "Motif-prior calibration summary must include all-scope metrics.")
motif_prior_required_cols <- c(
  "spearman_prior_vs_abs_delta",
  "spearman_prior_vs_evidence_score",
  "concordant_rows",
  "discordant_rows",
  "high_prior_weak_evidence_rows"
)
expect(all(motif_prior_required_cols %in% names(motif_prior_summary)),
       "Motif-prior summary is missing required calibration columns.")
motif_prior_gene <- data.table::fread(motif_prior_gene_file, nThread = 1)
gene_required_cols <- c(
  "motif_prior_review_class",
  "motif_prior_review_priority",
  "motif_concordance_score",
  "motif_discordance_score",
  "motif_unseen_prior_score"
)
expect(all(gene_required_cols %in% names(motif_prior_gene)),
       "Motif-prior gene-design table is missing review columns.")
priority <- suppressWarnings(
  as.numeric(motif_prior_gene$motif_prior_review_priority)
)
expect(all(is.finite(priority)),
       "Motif-prior review priority must be finite.")
expect(all(priority >= -1e-12 & priority <= 1 + 1e-12),
       "Motif-prior review priority must be bounded in [0, 1].")

message("Checking joint branch-allocation model")
joint_dir <- file.path(
  analysis_dir,
  "dominant_rdg_joint_branch_allocation"
)
joint_rows_file <- file.path(
  joint_dir,
  "dominant_rdg_joint_branch_allocation_rows.csv"
)
joint_gene_file <- file.path(
  joint_dir,
  "dominant_rdg_joint_branch_allocation_gene_design_summary.csv"
)
joint_summary_file <- file.path(
  joint_dir,
  "dominant_rdg_joint_branch_allocation_summary_metrics.csv"
)
expect_file(joint_rows_file)
expect_file(joint_gene_file)
expect_file(joint_summary_file)
joint_rows <- data.table::fread(joint_rows_file, nThread = 1)
expect(nrow(joint_rows) > 0, "Joint branch-allocation rows are empty.")
case_prop_cols <- c(
  "joint_case_prop_clean_CDS",
  "joint_case_prop_leader_uORF",
  "joint_case_prop_overlapping_uORF",
  "joint_case_prop_other"
)
control_prop_cols <- c(
  "joint_control_prop_clean_CDS",
  "joint_control_prop_leader_uORF",
  "joint_control_prop_overlapping_uORF",
  "joint_control_prop_other"
)
expect(all(case_prop_cols %in% names(joint_rows)),
       "Joint rows are missing case posterior proportion columns.")
expect(all(control_prop_cols %in% names(joint_rows)),
       "Joint rows are missing control posterior proportion columns.")
case_sums <- rowSums(joint_rows[, ..case_prop_cols], na.rm = FALSE)
control_sums <- rowSums(joint_rows[, ..control_prop_cols], na.rm = FALSE)
expect(all(is.finite(case_sums)) && max(abs(case_sums - 1)) < 1e-8,
       "Joint case posterior proportions must sum to one.")
expect(all(is.finite(control_sums)) && max(abs(control_sums - 1)) < 1e-8,
       "Joint control posterior proportions must sum to one.")
expect("joint_review_priority" %in% names(joint_rows),
       "Joint rows must include review priority.")
joint_priority <- suppressWarnings(as.numeric(joint_rows$joint_review_priority))
expect(all(is.finite(joint_priority)),
       "Joint row review priority must be finite.")
expect(all(joint_priority >= -1e-12 & joint_priority <= 1 + 1e-12),
       "Joint row review priority must be bounded in [0, 1].")
joint_gene <- data.table::fread(joint_gene_file, nThread = 1)
expect("joint_top_review_priority" %in% names(joint_gene),
       "Joint gene-design summary must include top review priority.")
expect(all(is.finite(as.numeric(joint_gene$joint_top_review_priority))),
       "Joint gene-design review priority must be finite.")
joint_summary <- data.table::fread(joint_summary_file, nThread = 1)
expect("all" %in% joint_summary$scope,
       "Joint branch-allocation summary must include all-scope metrics.")

message("Checking hierarchical joint branch-allocation model")
hier_joint_dir <- file.path(
  analysis_dir,
  "dominant_rdg_hierarchical_joint_allocation"
)
hier_joint_branch_file <- file.path(
  hier_joint_dir,
  "dominant_rdg_hierarchical_joint_allocation_branch_effects.csv"
)
hier_joint_gene_file <- file.path(
  hier_joint_dir,
  "dominant_rdg_hierarchical_joint_allocation_gene_design_summary.csv"
)
hier_joint_summary_file <- file.path(
  hier_joint_dir,
  "dominant_rdg_hierarchical_joint_allocation_summary_metrics.csv"
)
hier_joint_study_file <- file.path(
  hier_joint_dir,
  "dominant_rdg_hierarchical_joint_allocation_study_effects.csv"
)
expect_file(hier_joint_branch_file)
expect_file(hier_joint_gene_file)
expect_file(hier_joint_summary_file)
expect_file(hier_joint_study_file)
hier_joint_branch <- data.table::fread(hier_joint_branch_file, nThread = 1)
expect(nrow(hier_joint_branch) > 0,
       "Hierarchical joint branch-effect table is empty.")
hier_joint_required <- c(
  "hierarchical_joint_delta",
  "hierarchical_joint_delta_lower",
  "hierarchical_joint_delta_upper",
  "hierarchical_joint_delta_q",
  "hierarchical_joint_support_class",
  "hierarchical_joint_priority"
)
expect(all(hier_joint_required %in% names(hier_joint_branch)),
       "Hierarchical joint branch-effect table is missing required columns.")
hier_joint_priority <- suppressWarnings(
  as.numeric(hier_joint_branch$hierarchical_joint_priority)
)
expect(all(is.finite(hier_joint_priority)),
       "Hierarchical joint priority must be finite.")
expect(all(hier_joint_priority >= -1e-12 & hier_joint_priority <= 1 + 1e-12),
       "Hierarchical joint priority must be bounded in [0, 1].")
hier_joint_gene <- data.table::fread(hier_joint_gene_file, nThread = 1)
expect("hierarchical_joint_top_priority" %in% names(hier_joint_gene),
       "Hierarchical joint gene-design table must include top priority.")
expect(all(is.finite(as.numeric(
  hier_joint_gene$hierarchical_joint_top_priority
))),
"Hierarchical joint gene-design top priority must be finite.")
hier_joint_summary <- data.table::fread(hier_joint_summary_file, nThread = 1)
expect("all" %in% hier_joint_summary$scope,
       "Hierarchical joint summary must include all-scope metrics.")

message("Checking empirical-Bayes Dirichlet-multinomial branch model")
dm_dir <- file.path(
  analysis_dir,
  "dominant_rdg_dirichlet_multinomial"
)
dm_rows_file <- file.path(
  dm_dir,
  "dominant_rdg_dirichlet_multinomial_branch_rows.csv"
)
dm_branch_file <- file.path(
  dm_dir,
  "dominant_rdg_dirichlet_multinomial_branch_effects.csv"
)
dm_gene_file <- file.path(
  dm_dir,
  "dominant_rdg_dirichlet_multinomial_gene_design_summary.csv"
)
dm_summary_file <- file.path(
  dm_dir,
  "dominant_rdg_dirichlet_multinomial_summary_metrics.csv"
)
dm_study_file <- file.path(
  dm_dir,
  "dominant_rdg_dirichlet_multinomial_study_effects.csv"
)
expect_file(dm_rows_file)
expect_file(dm_branch_file)
expect_file(dm_gene_file)
expect_file(dm_summary_file)
expect_file(dm_study_file)
dm_rows <- data.table::fread(dm_rows_file, nThread = 1)
expect(nrow(dm_rows) > 0,
       "Dirichlet-multinomial branch-row table is empty.")
dm_row_required <- c(
  "dm_prior_source",
  "dm_posterior_baseline_prop",
  "dm_case_prop",
  "dm_delta",
  "dm_delta_p",
  "dm_delta_q",
  "dm_row_priority",
  "dm_row_support_class"
)
expect(all(dm_row_required %in% names(dm_rows)),
       "Dirichlet-multinomial branch-row table is missing required columns.")
dm_row_priority <- suppressWarnings(as.numeric(dm_rows$dm_row_priority))
expect(all(is.finite(dm_row_priority)),
       "Dirichlet-multinomial row priority must be finite.")
expect(all(dm_row_priority >= -1e-12 & dm_row_priority <= 1 + 1e-12),
       "Dirichlet-multinomial row priority must be bounded in [0, 1].")
baseline <- suppressWarnings(as.numeric(dm_rows$dm_posterior_baseline_prop))
expect(all(!is.na(baseline)),
       "Dirichlet-multinomial posterior baselines must not be NA.")
expect(all(baseline >= -1e-12 & baseline <= 1 + 1e-12),
       "Dirichlet-multinomial posterior baselines must be probabilities.")
dm_branch <- data.table::fread(dm_branch_file, nThread = 1)
expect(nrow(dm_branch) > 0,
       "Dirichlet-multinomial branch-effect table is empty.")
dm_branch_required <- c(
  "dm_delta",
  "dm_delta_lower",
  "dm_delta_upper",
  "dm_delta_q",
  "dm_support_class",
  "dm_priority"
)
expect(all(dm_branch_required %in% names(dm_branch)),
       "Dirichlet-multinomial branch-effect table is missing required columns.")
dm_priority <- suppressWarnings(as.numeric(dm_branch$dm_priority))
expect(all(is.finite(dm_priority)),
       "Dirichlet-multinomial branch priority must be finite.")
expect(all(dm_priority >= -1e-12 & dm_priority <= 1 + 1e-12),
       "Dirichlet-multinomial branch priority must be bounded in [0, 1].")
dm_gene <- data.table::fread(dm_gene_file, nThread = 1)
expect("dm_top_priority" %in% names(dm_gene),
       "Dirichlet-multinomial gene-design table must include top priority.")
expect(all(is.finite(as.numeric(dm_gene$dm_top_priority))),
       "Dirichlet-multinomial gene-design top priority must be finite.")
dm_summary <- data.table::fread(dm_summary_file, nThread = 1)
expect("all" %in% dm_summary$scope,
       "Dirichlet-multinomial summary must include all-scope metrics.")

message("Checking RDG model-agreement layer")
agreement_dir <- file.path(analysis_dir, "dominant_rdg_model_agreement")
agreement_branch_file <- file.path(
  agreement_dir,
  "dominant_rdg_model_agreement_branch_comparison.csv"
)
agreement_branch_summary_file <- file.path(
  agreement_dir,
  "dominant_rdg_model_agreement_branch_summary.csv"
)
agreement_gene_file <- file.path(
  agreement_dir,
  "dominant_rdg_model_agreement_gene_design.csv"
)
agreement_summary_file <- file.path(
  agreement_dir,
  "dominant_rdg_model_agreement_summary_metrics.csv"
)
expect_file(agreement_branch_file)
expect_file(agreement_branch_summary_file)
expect_file(agreement_gene_file)
expect_file(agreement_summary_file)
agreement_branch <- data.table::fread(agreement_branch_file, nThread = 1)
expect(nrow(agreement_branch) > 0,
       "Model-agreement branch-comparison table is empty.")
expect(all(c(
  "branch_agreement_class",
  "hierarchical_joint_supported",
  "dm_supported",
  "direction_agree",
  "direction_opposite"
) %in% names(agreement_branch)),
"Model-agreement branch-comparison table is missing required columns.")
expect(!any(agreement_branch$direction_agree &
              agreement_branch$direction_opposite, na.rm = TRUE),
       "A branch cannot agree and disagree in direction simultaneously.")
agreement_gene <- data.table::fread(agreement_gene_file, nThread = 1)
expect(nrow(agreement_gene) > 0,
       "Model-agreement gene-design table is empty.")
expect(all(c(
  "model_agreement_tier",
  "model_agreement_score",
  "model_agreement_review_priority",
  "model_support_count",
  "model_support_signature",
  "model_disagreement_flags"
) %in% names(agreement_gene)),
"Model-agreement gene-design table is missing required columns.")
agreement_score <- suppressWarnings(as.numeric(
  agreement_gene$model_agreement_score
))
expect(all(is.finite(agreement_score)),
       "Model-agreement scores must be finite.")
expect(all(agreement_score >= -1e-12 & agreement_score <= 1 + 1e-12),
       "Model-agreement scores must be bounded in [0, 1].")
expect(all(agreement_gene$model_support_count >= 0 &
             agreement_gene$model_support_count <= 4),
       "Model support count must be bounded in [0, 4].")
agreement_summary <- data.table::fread(agreement_summary_file, nThread = 1)
expect("all" %in% agreement_summary$scope,
       "Model-agreement summary must include all-scope metrics.")

message("Checking leave-one-study-out model-agreement layer")
loo_dir <- file.path(analysis_dir, "dominant_rdg_model_agreement_loo")
loo_study_file <- file.path(
  loo_dir,
  "dominant_rdg_model_agreement_study_evidence.csv"
)
loo_rows_file <- file.path(
  loo_dir,
  "dominant_rdg_model_agreement_loo_rows.csv"
)
loo_branch_file <- file.path(
  loo_dir,
  "dominant_rdg_model_agreement_loo_branch_summary.csv"
)
loo_gene_file <- file.path(
  loo_dir,
  "dominant_rdg_model_agreement_loo_gene_design.csv"
)
loo_summary_file <- file.path(
  loo_dir,
  "dominant_rdg_model_agreement_loo_summary_metrics.csv"
)
expect_file(loo_study_file)
expect_file(loo_rows_file)
expect_file(loo_branch_file)
expect_file(loo_gene_file)
expect_file(loo_summary_file)
loo_study <- data.table::fread(loo_study_file, nThread = 1)
expect(nrow(loo_study) > 0, "LOO study-evidence table is empty.")
expect(all(c(
  "study", "study_top_branch_class", "study_top_evidence_class",
  "study_top_hierarchical_delta", "study_top_dm_delta",
  "study_branch_summary"
) %in% names(loo_study)),
"LOO study-evidence table is missing required columns.")
loo_rows <- data.table::fread(loo_rows_file, nThread = 1)
expect(nrow(loo_rows) > 0, "LOO omission-row table is empty.")
expect(all(c(
  "omitted_study", "model_studies", "remaining_studies",
  "loo_same_direction", "loo_support_retained", "loo_support_class"
) %in% names(loo_rows)),
"LOO omission-row table is missing required columns.")
loo_gene <- data.table::fread(loo_gene_file, nThread = 1)
expect(nrow(loo_gene) == nrow(agreement_gene),
       "LOO gene-design table must preserve all model-agreement rows.")
expect(all(c(
  "loo_agreement_class", "loo_consensus_class", "loo_stability_score",
  "loo_min_direction_fraction", "loo_min_support_retention",
  "loo_min_model_studies"
) %in% names(loo_gene)),
"LOO gene-design table is missing required columns.")
loo_score <- suppressWarnings(as.numeric(loo_gene$loo_stability_score))
expect(all(is.finite(loo_score)), "LOO stability scores must be finite.")
expect(all(loo_score >= -1e-12 & loo_score <= 1 + 1e-12),
       "LOO stability scores must be bounded in [0, 1].")
loo_summary <- data.table::fread(loo_summary_file, nThread = 1)
expect("all" %in% loo_summary$scope,
       "LOO summary must include all-scope metrics.")
expect(any(loo_gene$loo_consensus_class == "consensus_loo_support_robust"),
       "LOO output should contain at least one support-robust consensus row.")

message("Checking grouped-omission model-agreement transportability layer")
transport_dir <- file.path(
  analysis_dir, "dominant_rdg_model_agreement_transport"
)
transport_rows_file <- file.path(
  transport_dir,
  "dominant_rdg_model_agreement_transport_omission_rows.csv"
)
transport_axis_file <- file.path(
  transport_dir,
  "dominant_rdg_model_agreement_transport_gene_axis_summary.csv"
)
transport_gene_file <- file.path(
  transport_dir,
  "dominant_rdg_model_agreement_transport_gene_design.csv"
)
transport_summary_file <- file.path(
  transport_dir,
  "dominant_rdg_model_agreement_transport_summary_metrics.csv"
)
expect_file(transport_rows_file)
expect_file(transport_axis_file)
expect_file(transport_gene_file)
expect_file(transport_summary_file)
transport_rows <- data.table::fread(transport_rows_file, nThread = 1)
expect(nrow(transport_rows) > 0,
       "Grouped-omission transport row table is empty.")
expect(all(c(
  "transport_axis", "omitted_group", "omitted_studies",
  "remaining_studies", "transport_same_direction",
  "transport_support_retained"
) %in% names(transport_rows)),
"Grouped-omission transport rows are missing required columns.")
expect(all(transport_rows$omitted_studies >= 2),
       "Grouped omission must remove at least two studies.")
expect(all(transport_rows$remaining_studies >= 2),
       "Grouped omission must leave at least two studies.")
transport_gene <- data.table::fread(transport_gene_file, nThread = 1)
expect(nrow(transport_gene) == nrow(agreement_gene),
       "Transport gene-design table must preserve all agreement rows.")
expect(all(c(
  "transport_agreement_class", "transport_consensus_class",
  "transport_stability_score", "transport_min_direction_fraction",
  "transport_min_support_retention", "transport_axis_coverage_fraction"
) %in% names(transport_gene)),
"Transport gene-design table is missing required columns.")
transport_score <- suppressWarnings(as.numeric(
  transport_gene$transport_stability_score
))
expect(all(is.finite(transport_score)),
       "Transport stability scores must be finite.")
expect(all(transport_score >= -1e-12 & transport_score <= 1 + 1e-12),
       "Transport stability scores must be bounded in [0, 1].")
transport_axis_coverage <- suppressWarnings(as.numeric(
  transport_gene$transport_axis_coverage_fraction
))
expect(all(is.finite(transport_axis_coverage)),
       "Transport axis-coverage fractions must be finite.")
expect(all(transport_axis_coverage >= -1e-12 &
             transport_axis_coverage <= 1 + 1e-12),
       "Transport axis-coverage fractions must be bounded in [0, 1].")
transport_axis <- data.table::fread(transport_axis_file, nThread = 1)
expect(any(
  transport_axis$gene_symbol == "IFIH1" &
    transport_axis$design_family == "viral_infection" &
    transport_axis$transport_axis == "tissue_context" &
    transport_axis$transport_axis_class == "transport_direction_fragile"
), "IFIH1 viral tissue-family omission regression result is missing.")

message("Checking viral context-interaction layer")
context_interaction_dir <- file.path(
  analysis_dir, "dominant_rdg_context_interactions"
)
context_interaction_branch_file <- file.path(
  context_interaction_dir,
  "dominant_rdg_context_interaction_branch_context_effects.csv"
)
context_interaction_gene_file <- file.path(
  context_interaction_dir,
  "dominant_rdg_context_interaction_gene_design_summary.csv"
)
context_interaction_review_file <- file.path(
  context_interaction_dir,
  "dominant_rdg_context_interaction_review.csv"
)
context_interaction_summary_file <- file.path(
  context_interaction_dir,
  "dominant_rdg_context_interaction_summary_metrics.csv"
)
expect_file(context_interaction_branch_file)
expect_file(context_interaction_gene_file)
expect_file(context_interaction_review_file)
expect_file(context_interaction_summary_file)
context_interaction_branch <- data.table::fread(
  context_interaction_branch_file, nThread = 1
)
expect(nrow(context_interaction_branch) > 0,
       "Context-interaction branch table is empty.")
expect(all(c(
  "context_axis", "context_value", "branch_class",
  "context_studies", "complement_studies", "context_delta",
  "complement_delta", "interaction_delta", "interaction_delta_q",
  "context_interaction_class", "context_interaction_score"
) %in% names(context_interaction_branch)),
"Context-interaction branch table is missing required columns.")
context_score <- suppressWarnings(as.numeric(
  context_interaction_branch$context_interaction_score
))
expect(all(is.finite(context_score)),
       "Context-interaction scores must be finite.")
expect(all(context_score >= -1e-12 & context_score <= 1 + 1e-12),
       "Context-interaction scores must be bounded in [0, 1].")
ifih1_lung <- context_interaction_branch[
  gene_symbol == "IFIH1" &
    design_family == "viral_infection" &
    context_axis == "tissue_context" &
    context_value == "lung" &
    branch_class %in% c("clean_CDS", "leader_uORF")
]
expect(nrow(ifih1_lung) == 2,
       "IFIH1 viral lung clean-CDS/leader-uORF context rows are missing.")
expect(all(ifih1_lung$context_studies >= 2 &
             ifih1_lung$complement_studies >= 2),
       "IFIH1 lung context test must include context and complement studies.")
expect(all(ifih1_lung$context_interaction_class ==
             "context_direction_reversal"),
       "IFIH1 lung context rows should be direction-reversal calls.")
expect(any(ifih1_lung$branch_class == "clean_CDS" &
             ifih1_lung$context_delta > 0 &
             ifih1_lung$complement_delta < 0),
       "IFIH1 clean-CDS lung effect should oppose the non-lung complement.")
context_interaction_gene <- data.table::fread(
  context_interaction_gene_file, nThread = 1
)
expect(nrow(context_interaction_gene) > 0,
       "Context-interaction gene-design summary is empty.")
expect(any(context_interaction_gene$context_interaction_gene_class ==
             "context_direction_reversal"),
       "Context-interaction summary should contain direction reversals.")

message("Checking partial-pooled context-interaction layer")
partial_pool_dir <- file.path(
  analysis_dir, "dominant_rdg_context_partial_pooling"
)
partial_pool_rows_file <- file.path(
  partial_pool_dir,
  "dominant_rdg_context_partial_pooling_rows.csv"
)
partial_pool_gene_file <- file.path(
  partial_pool_dir,
  "dominant_rdg_context_partial_pooling_gene_design_summary.csv"
)
partial_pool_viral_file <- file.path(
  partial_pool_dir,
  "dominant_rdg_context_partial_pooling_viral_review.csv"
)
partial_pool_summary_file <- file.path(
  partial_pool_dir,
  "dominant_rdg_context_partial_pooling_summary_metrics.csv"
)
expect_file(partial_pool_rows_file)
expect_file(partial_pool_gene_file)
expect_file(partial_pool_viral_file)
expect_file(partial_pool_summary_file)
partial_pool_rows <- data.table::fread(partial_pool_rows_file, nThread = 1)
expect(nrow(partial_pool_rows) > 0,
       "Partial-pooled context rows are empty.")
partial_pool_required <- c(
  "pooled_interaction_delta",
  "pooled_interaction_q",
  "partial_pooled_context_class",
  "partial_pooling_survival_class",
  "partial_pooling_review_priority",
  "partial_pooling_weight",
  "pooling_group_contexts"
)
expect(all(partial_pool_required %in% names(partial_pool_rows)),
       "Partial-pooled context rows are missing required columns.")
partial_pool_q <- suppressWarnings(as.numeric(
  partial_pool_rows$pooled_interaction_q
))
expect(all(is.finite(partial_pool_q)),
       "Partial-pooled context q-values must be finite.")
expect(all(partial_pool_q >= -1e-12 & partial_pool_q <= 1 + 1e-12),
       "Partial-pooled context q-values must be bounded in [0, 1].")
for (column in c("partial_pooling_review_priority",
                 "partial_pooling_weight")) {
  value <- suppressWarnings(as.numeric(partial_pool_rows[[column]]))
  expect(all(is.finite(value)),
         paste0("Partial-pooled context ", column, " must be finite."))
  expect(all(value >= -1e-12 & value <= 1 + 1e-12),
         paste0("Partial-pooled context ", column,
                " must be bounded in [0, 1]."))
}
expect(any(partial_pool_rows$partial_pooled_context_class ==
             "partial_pooled_supported_context_shift"),
       "Partial-pooling should retain at least one supported context shift.")
ifih1_pool <- partial_pool_rows[
  gene_symbol == "IFIH1" &
    design_family == "viral_infection" &
    context_axis == "tissue_context" &
    context_value == "lung" &
    branch_class %in% c("clean_CDS", "leader_uORF")
]
expect(nrow(ifih1_pool) == 2,
       "IFIH1 viral lung partial-pooling rows are missing.")
expect(all(ifih1_pool$partial_pooled_context_class ==
             "single_context_supported_not_pooled"),
       "IFIH1 lung rows should be labeled single-context, not pooled support.")
expect(all(ifih1_pool$partial_pooling_survival_class ==
             "single_context_raw_signal_not_pooled"),
       "IFIH1 lung rows should preserve raw single-context evidence.")
expect(max(abs(ifih1_pool$pooled_interaction_delta -
                 ifih1_pool$interaction_delta)) < 1e-10,
       "Single-context IFIH1 estimates should equal raw interaction deltas.")
partial_pool_gene <- data.table::fread(partial_pool_gene_file, nThread = 1)
expect(nrow(partial_pool_gene) > 0,
       "Partial-pooled gene-design summary is empty.")
expect(any(partial_pool_gene$partial_pooling_gene_class ==
             "partial_pooled_context_supported"),
       "Partial-pooled gene-design summary should contain supported rows.")

message("Checking leave-one-study-out atlas integration")
for (column in c("loo_stability_component", "loo_fragility_component")) {
  expect(column %in% names(atlas_cards),
         paste0("Atlas cards must include ", column, "."))
  value <- suppressWarnings(as.numeric(atlas_cards[[column]]))
  expect(all(is.finite(value)),
         paste0("Atlas ", column, " must be finite."))
  expect(all(value >= -1e-12 & value <= 1 + 1e-12),
         paste0("Atlas ", column, " must be bounded in [0, 1]."))
}
expect("loo_consensus_class" %in% names(atlas_cards),
       "Atlas cards must include LOO consensus class.")

message("Checking grouped-omission transport atlas integration")
for (column in c(
  "transport_stability_component", "transport_fragility_component"
)) {
  expect(column %in% names(atlas_cards),
         paste0("Atlas cards must include ", column, "."))
  value <- suppressWarnings(as.numeric(atlas_cards[[column]]))
  expect(all(is.finite(value)),
         paste0("Atlas ", column, " must be finite."))
  expect(all(value >= -1e-12 & value <= 1 + 1e-12),
         paste0("Atlas ", column, " must be bounded in [0, 1]."))
}
expect("transport_consensus_class" %in% names(atlas_cards),
       "Atlas cards must include transport consensus class.")

message("Checking context-interaction atlas integration")
expect("context_interaction_component" %in% names(atlas_cards),
       "Atlas cards must include context-interaction component.")
context_component <- suppressWarnings(as.numeric(
  atlas_cards$context_interaction_component
))
expect(all(is.finite(context_component)),
       "Atlas context-interaction component must be finite.")
expect(all(context_component >= -1e-12 & context_component <= 1 + 1e-12),
       "Atlas context-interaction component must be bounded in [0, 1].")
expect("context_interaction_gene_class" %in% names(atlas_cards),
       "Atlas cards must include context-interaction gene class.")
expect(any(
  atlas_cards$gene_symbol == "IFIH1" &
    atlas_cards$design_family == "viral_infection" &
    atlas_cards$context_interaction_gene_class ==
	    "context_direction_reversal"
), "IFIH1 viral context-interaction class is missing from atlas cards.")

message("Checking partial-pooling atlas integration")
partial_pool_atlas_required <- c(
  "partial_pooling_gene_class",
  "partial_pooling_top_class",
  "partial_pooling_review_priority",
  "partial_pooling_component",
  "partial_pooling_top_axis",
  "partial_pooling_top_value",
  "partial_pooling_top_branch_class",
  "partial_pooling_top_delta",
  "partial_pooling_top_q"
)
expect(all(partial_pool_atlas_required %in% names(atlas_cards)),
       "Atlas cards are missing partial-pooling columns.")
partial_component <- suppressWarnings(as.numeric(
  atlas_cards$partial_pooling_component
))
expect(all(is.finite(partial_component)),
       "Atlas partial-pooling component must be finite.")
expect(all(partial_component >= -1e-12 & partial_component <= 1 + 1e-12),
       "Atlas partial-pooling component must be bounded in [0, 1].")
expect(any(
  atlas_cards$gene_symbol == "ATF4" &
    atlas_cards$design_family == "viral_infection" &
    atlas_cards$partial_pooling_gene_class ==
    "partial_pooled_context_supported"
), "ATF4 viral partial-pooled support is missing from atlas cards.")
expect(any(
  atlas_cards$gene_symbol == "TRIB3" &
    atlas_cards$design_family == "viral_infection" &
    atlas_cards$partial_pooling_gene_class ==
    "partial_pooled_context_supported"
), "TRIB3 viral partial-pooled support is missing from atlas cards.")
expect(any(
  atlas_cards$gene_symbol == "IFIH1" &
    atlas_cards$design_family == "viral_infection" &
    atlas_cards$partial_pooling_gene_class ==
    "single_context_raw_signal_not_pooled"
), "IFIH1 viral single-context partial-pooling label is missing from atlas cards.")
expect(any(grepl("partial_pooled_context_supported",
                 atlas_cards$warnings, fixed = TRUE)),
       "Atlas warnings should include partial-pooled context support.")
expect(any(grepl("single_context_context_signal",
                 atlas_cards$warnings, fixed = TRUE)),
       "Atlas warnings should include single-context context signals.")

message("Checking pooling and protocol-bias audit")
pooling_bias_dir <- file.path(
  analysis_dir, "dominant_rdg_pooling_bias_audit"
)
pooling_claim_file <- file.path(pooling_bias_dir,
                                "pooling_claim_language.csv")
normalization_limits_file <- file.path(pooling_bias_dir,
                                       "normalization_limits.csv")
pooling_gain_loss_file <- file.path(pooling_bias_dir,
                                    "pooling_gain_loss_summary.csv")
pooling_tier_file <- file.path(pooling_bias_dir,
                               "pooling_contrast_tier_summary.csv")
design_pool_file <- file.path(pooling_bias_dir,
                              "design_family_pooling_bias_summary.csv")
protocol_run_file <- file.path(pooling_bias_dir,
                               "protocol_bias_run_summary.csv")
partial_tradeoff_file <- file.path(pooling_bias_dir,
                                   "partial_pooling_tradeoff_summary.csv")
pooling_metrics_file <- file.path(pooling_bias_dir,
                                  "pooling_bias_summary_metrics.csv")
pooling_report_file <- file.path(pooling_bias_dir, "pooling_bias_report.md")
expect_file(pooling_claim_file)
expect_file(normalization_limits_file)
expect_file(pooling_gain_loss_file)
expect_file(pooling_tier_file)
expect_file(design_pool_file)
expect_file(protocol_run_file)
expect_file(partial_tradeoff_file)
expect_file(pooling_metrics_file)
expect_file(pooling_report_file)

pooling_claim <- data.table::fread(pooling_claim_file, nThread = 1)
expect("A_project_specific_matched_contrast" %in%
         pooling_claim$claim_tier,
       "Pooling claim language must include project-specific matched contrasts.")
expect("F_all_merged_visual_background" %in% pooling_claim$claim_tier,
       "Pooling claim language must constrain all-merged profiles.")

normalization_limits <- data.table::fread(normalization_limits_file,
                                          nThread = 1)
expect("branch_allocation_delta" %in% normalization_limits$quantity,
       "Normalization limits must describe branch allocation deltas.")
expect(any(grepl("does not control|protocol|inhibitor",
                 paste(normalization_limits$does_not_control, collapse = " "),
                 ignore.case = TRUE)),
       "Normalization limits must mention remaining protocol/inhibitor bias.")

pooling_gain_loss <- data.table::fread(pooling_gain_loss_file, nThread = 1)
single_ge30 <- pooling_gain_loss[
  metric == "single_run_fraction_ge_30", value_numeric
][1]
merged_ge30 <- pooling_gain_loss[
  metric == "merged_n2_fraction_ge_30", value_numeric
][1]
expect(is.finite(single_ge30) && is.finite(merged_ge30) &&
         merged_ge30 > single_ge30,
       "Merged replicate groups should improve small-uORF >=30 coverage.")

pooling_tier <- data.table::fread(pooling_tier_file, nThread = 1)
expect(nrow(pooling_tier) > 0,
       "Pooling contrast tier summary is empty.")
expect(any(pooling_tier$proof_tier ==
             "A_exact_strict_project_matched_control"),
       "Pooling tier summary must include exact strict project-matched rows.")
expect(any(pooling_tier$proof_tier ==
             "C_relaxed_study_cell_gene_pool"),
       "Pooling tier summary must include relaxed pooled rows.")

design_pool <- data.table::fread(design_pool_file, nThread = 1)
expect(nrow(design_pool) > 0,
       "Design-family pooling-bias summary is empty.")
expect(all(c("pooling_hazard_score", "pooling_hazard_class",
             "top_protocol_classes") %in% names(design_pool)),
       "Design-family pooling-bias summary is missing required columns.")
hazard_score <- suppressWarnings(as.numeric(design_pool$pooling_hazard_score))
expect(all(is.finite(hazard_score)),
       "Pooling hazard scores must be finite.")
expect(all(hazard_score >= -1e-12 & hazard_score <= 1 + 1e-12),
       "Pooling hazard scores must be bounded in [0, 1].")

protocol_run <- data.table::fread(protocol_run_file, nThread = 1)
expect(nrow(protocol_run) > 0,
       "Protocol-bias run summary is empty.")
expect(all(c("protocol_bias_class", "condition_bias_family",
             "median_tis_peak_strength",
             "median_global_canonical_frame_fraction") %in%
             names(protocol_run)),
       "Protocol-bias run summary is missing required columns.")

partial_tradeoff <- data.table::fread(partial_tradeoff_file, nThread = 1)
expect(nrow(partial_tradeoff) > 0,
       "Partial-pooling tradeoff summary is empty.")
expect(all(c("raw_supported_rows", "pooled_supported_rows",
             "single_context_rows", "raw_shrunk_rows") %in%
             names(partial_tradeoff)),
       "Partial-pooling tradeoff summary is missing required columns.")

pooling_metrics <- data.table::fread(pooling_metrics_file, nThread = 1)
expect(all(c("metric", "value") %in% names(pooling_metrics)),
       "Pooling-bias metrics must have metric/value columns.")
pooling_metric_value <- function(metric_name) {
  suppressWarnings(as.numeric(pooling_metrics[metric == metric_name, value][1]))
}
expect(pooling_metric_value("grouped_contrast_rows") > 0,
       "Pooling-bias metrics should report grouped contrast rows.")
expect(pooling_metric_value("partial_pooling_rows") > 0,
       "Pooling-bias metrics should report partial-pooling rows.")

message("Checking allocation-versus-output disagreement audit")
allocation_output_dir <- file.path(
  analysis_dir, "dominant_rdg_allocation_output_disagreements"
)
allocation_output_audit_file <- file.path(
  allocation_output_dir,
  "dominant_rdg_allocation_output_disagreement_audit.csv"
)
allocation_output_review_file <- file.path(
  allocation_output_dir,
  "dominant_rdg_allocation_output_disagreement_review.csv"
)
allocation_output_context_file <- file.path(
  allocation_output_dir,
  "dominant_rdg_allocation_output_disagreement_contexts.csv"
)
allocation_output_summary_file <- file.path(
  allocation_output_dir,
  "dominant_rdg_allocation_output_disagreement_summary_metrics.csv"
)
expect_file(allocation_output_audit_file)
expect_file(allocation_output_review_file)
expect_file(allocation_output_context_file)
expect_file(allocation_output_summary_file)
allocation_output_review <- data.table::fread(
  allocation_output_review_file, nThread = 1
)
expect(nrow(allocation_output_review) > 0,
       "Allocation-output disagreement review is empty.")
expect(all(c(
  "gene_symbol", "design_family", "disagreement_audit_class",
  "disagreement_review_priority", "compensation_evidence_score",
  "heterogeneity_risk_score", "clean_cds_output_pct",
  "clean_cds_allocation_delta"
) %in% names(allocation_output_review)),
"Allocation-output disagreement review is missing required columns.")
for (column in c(
  "disagreement_review_priority",
  "compensation_evidence_score",
  "heterogeneity_risk_score"
)) {
  value <- suppressWarnings(as.numeric(allocation_output_review[[column]]))
  expect(all(is.finite(value)),
         paste0("Allocation-output ", column, " must be finite."))
  expect(all(value >= -1e-12 & value <= 1 + 1e-12),
         paste0("Allocation-output ", column, " must be bounded in [0, 1]."))
}
expect(any(allocation_output_review$gene_symbol == "SESN2" &
             allocation_output_review$design_family ==
             "nutrient_starvation"),
       "SESN2 nutrient-starvation disagreement row is missing.")
expect(any(allocation_output_review$gene_symbol == "ATF3" &
             allocation_output_review$design_family ==
             "nutrient_starvation"),
       "ATF3 nutrient-starvation disagreement row is missing.")
expect(any(allocation_output_review$clean_cds_output_pct > 0 &
             allocation_output_review$clean_cds_allocation_delta < 0),
       "Audit should contain output-up/allocation-down rows.")

message("Checking therapeutic-hypothesis layer")
therapeutic_dir <- file.path(
  analysis_dir, "dominant_rdg_therapeutic_hypotheses"
)
therapeutic_scores_file <- file.path(
  therapeutic_dir,
  "dominant_rdg_therapeutic_hypothesis_scores.csv"
)
therapeutic_review_file <- file.path(
  therapeutic_dir,
  "dominant_rdg_therapeutic_hypothesis_review.csv"
)
therapeutic_gene_file <- file.path(
  therapeutic_dir,
  "dominant_rdg_therapeutic_hypothesis_gene_evidence.csv"
)
therapeutic_summary_file <- file.path(
  therapeutic_dir,
  "dominant_rdg_therapeutic_hypothesis_summary_metrics.csv"
)
expect_file(therapeutic_scores_file)
expect_file(therapeutic_review_file)
expect_file(therapeutic_gene_file)
expect_file(therapeutic_summary_file)
therapeutic_scores <- data.table::fread(
  therapeutic_scores_file, nThread = 1
)
expect(nrow(therapeutic_scores) > 0,
       "Therapeutic-hypothesis score table is empty.")
expect(all(c(
  "therapeutic_hypothesis", "hypothesis_label",
  "therapeutic_hypothesis_score", "trial_readiness_class",
  "intervention_stage", "example_interventions", "primary_genes"
) %in% names(therapeutic_scores)),
"Therapeutic-hypothesis scores are missing required columns.")
therapeutic_score <- suppressWarnings(as.numeric(
  therapeutic_scores$therapeutic_hypothesis_score
))
expect(all(is.finite(therapeutic_score)),
       "Therapeutic-hypothesis scores must be finite.")
expect(all(therapeutic_score >= -1e-12 & therapeutic_score <= 1 + 1e-12),
       "Therapeutic-hypothesis scores must be bounded in [0, 1].")
expect(any(therapeutic_scores$therapeutic_hypothesis ==
             "JAK_STAT_IFN_attenuation"),
       "JAK/STAT IFN hypothesis is missing.")
expect(any(therapeutic_scores$therapeutic_hypothesis ==
             "ISR_eIF2B_ATF4_attenuation"),
       "ISR/eIF2B-ATF4 hypothesis is missing.")
main_scores <- therapeutic_scores[
  therapeutic_hypothesis %in% c(
    "JAK_STAT_IFN_attenuation",
    "ISR_eIF2B_ATF4_attenuation"
  )
]
expect(all(main_scores$therapeutic_hypothesis_score > 0.3),
       "Primary therapeutic hypotheses should have nontrivial evidence.")

message("Checking therapeutic perturbation-prediction layer")
perturbation_dir <- file.path(
  analysis_dir, "dominant_rdg_therapeutic_perturbation_predictions"
)
perturbation_gene_file <- file.path(
  perturbation_dir,
  "dominant_rdg_therapeutic_perturbation_gene_predictions.csv"
)
perturbation_summary_file <- file.path(
  perturbation_dir,
  "dominant_rdg_therapeutic_perturbation_hypothesis_summary.csv"
)
perturbation_review_file <- file.path(
  perturbation_dir,
  "dominant_rdg_therapeutic_perturbation_review.csv"
)
perturbation_metrics_file <- file.path(
  perturbation_dir,
  "dominant_rdg_therapeutic_perturbation_summary_metrics.csv"
)
expect_file(perturbation_gene_file)
expect_file(perturbation_summary_file)
expect_file(perturbation_review_file)
expect_file(perturbation_metrics_file)
perturbation_gene <- data.table::fread(
  perturbation_gene_file, nThread = 1
)
expect(nrow(perturbation_gene) > 0,
       "Therapeutic perturbation gene-prediction table is empty.")
expect(all(c(
  "therapeutic_hypothesis", "gene_symbol",
  "quantitative_prediction_confidence",
  "prediction_support_class",
  "disease_like_cds_pct_change",
  "half_normalization_cds_pct_change"
) %in% names(perturbation_gene)),
"Perturbation gene-prediction table is missing required columns.")
perturbation_confidence <- suppressWarnings(as.numeric(
  perturbation_gene$quantitative_prediction_confidence
))
expect(all(is.finite(perturbation_confidence)),
       "Perturbation prediction confidence must be finite.")
expect(all(perturbation_confidence >= -1e-12 &
             perturbation_confidence <= 1 + 1e-12),
       "Perturbation prediction confidence must be bounded in [0, 1].")
isr_targets <- perturbation_gene[
  therapeutic_hypothesis == "ISR_eIF2B_ATF4_attenuation" &
    gene_symbol %in% c("DDIT3", "ATF4", "PPP1R15A")
]
expect(nrow(isr_targets) == 3,
       "Core ISR perturbation targets DDIT3/ATF4/PPP1R15A are missing.")
expect(all(isr_targets$prediction_support_class %in% c(
  "quantitative_review_ready",
  "directional_prediction",
  "exploratory_prediction"
)),
"Core ISR perturbation targets should have reviewable prediction support.")
expect(all(isr_targets$half_normalization_cds_pct_change < 0),
       "ISR attenuation should predict reduced disease-like CDS output for core targets.")
perturbation_summary <- data.table::fread(
  perturbation_summary_file, nThread = 1
)
expect(nrow(perturbation_summary) > 0,
       "Perturbation hypothesis summary is empty.")
expect(any(perturbation_summary$hypothesis_label ==
             "ISR/eIF2B-ATF4 attenuation"),
       "Perturbation summary is missing ISR/eIF2B-ATF4 attenuation.")

message("Dominant-state pipeline invariant tests passed.")
