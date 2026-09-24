#!/usr/bin/env Rscript

# Single entry point for the dominant-cell-state analysis.
# Run from the historical RiboCrypt repo root, or from a standalone
# dominant_cell_states checkout with RIBOCRYPT_REPO pointing to RiboCrypt.

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

source_dominant_package_helpers <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    helper_dir <- file.path(here, "R")
    helper_files <- file.path(
      helper_dir,
      c("paths.R", "results-links.R", "pipeline.R")
    )
    if (all(file.exists(helper_files))) {
      for (helper in helper_files) sys.source(helper, envir = globalenv())
      return(invisible(TRUE))
    }
    candidate <- file.path(here, "dominant_cell_states", "R")
    helper_files <- file.path(
      candidate,
      c("paths.R", "results-links.R", "pipeline.R")
    )
    if (all(file.exists(helper_files))) {
      for (helper in helper_files) sys.source(helper, envir = globalenv())
      return(invisible(TRUE))
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  invisible(FALSE)
}

source_dominant_package_helpers()

analysis_results_dir <- function(analysis_dir, create = FALSE) {
  if (exists("dcs_results_dir", mode = "function")) {
    return(dcs_results_dir(analysis_dir, create = create))
  }
  results_dir <- file.path(analysis_dir, "results")
  if (isTRUE(create)) {
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(results_dir, mustWork = FALSE)
}

find_repo_root <- function(start = getwd()) {
  env_repo <- Sys.getenv("RIBOCRYPT_REPO", unset = "")
  if (nzchar(env_repo)) {
    env_repo <- normalizePath(env_repo, mustWork = FALSE)
    if (file.exists(file.path(env_repo, "DESCRIPTION"))) {
      return(env_repo)
    }
    stop("RIBOCRYPT_REPO does not point to a R package root: ", env_repo,
         call. = FALSE)
  }
  if (exists("dcs_ribocrypt_repo", mode = "function")) {
    helper_repo <- dcs_ribocrypt_repo(start = start, required = FALSE)
    if (nzchar(helper_repo)) {
      return(helper_repo)
    }
  }
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(here, "DESCRIPTION")) &&
        dir.exists(file.path(here, "dominant_cell_states"))) {
      return(here)
    }
    parent <- dirname(here)
    if (identical(parent, here)) {
      stop("Could not find RiboCrypt repo root from: ", start,
           ". Set RIBOCRYPT_REPO=/path/to/RiboCrypt when running from a ",
           "standalone dominant_cell_states checkout.",
           call. = FALSE)
    }
    here <- parent
  }
}

find_analysis_dir <- function(start = getwd(), repo_root = NULL) {
  if (exists("dcs_project_root", mode = "function")) {
    helper_root <- tryCatch(dcs_project_root(start), error = function(e) NULL)
    if (!is.null(helper_root)) return(helper_root)
  }
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (basename(here) == "dominant_cell_states" &&
        dir.exists(file.path(here, "scripts")) &&
        file.exists(file.path(here, "md",
                              "dominant_cell_state_scientific_note.md"))) {
      return(here)
    }
    candidate <- file.path(here, "dominant_cell_states")
    if (dir.exists(file.path(candidate, "scripts")) &&
        file.exists(file.path(candidate, "md",
                              "dominant_cell_state_scientific_note.md"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  if (!is.null(repo_root)) {
    candidate <- file.path(repo_root, "dominant_cell_states")
    if (dir.exists(file.path(candidate, "scripts"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  stop("Could not find dominant_cell_states analysis directory from: ",
       start, call. = FALSE)
}

load_ribocrypt <- function(repo_root) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Package 'devtools' is required to load RiboCrypt with devtools::load_all().",
         call. = FALSE)
  }
  devtools::load_all(repo_root, quiet = TRUE)
}

setup_runtime_env <- function(analysis_dir) {
  # Capture ORFik's already-correct config (resolved via the default,
  # un-redirected BiocFileCache) BEFORE isolating XDG/BFC state below.
  # Once isolated, ORFik::config() can no longer see that cached entry and
  # regenerates one from its own hardcoded, environment-naive default
  # ("~/Bio_data/..."), which is wrong on any server where Bio_data lives
  # elsewhere (e.g. "~/livemount/Bio_data" here). Verified directly: this
  # reproduces with nothing more than pointing XDG_CACHE_HOME/
  # XDG_CONFIG_HOME at a fresh empty directory.
  orfik_config <- if (requireNamespace("ORFik", quietly = TRUE)) {
    tryCatch(ORFik::config(), error = function(e) NULL)
  } else NULL

  runtime_root <- file.path(analysis_dir, ".runtime")
  cache_dir <- file.path(runtime_root, "xdg-cache")
  config_dir <- file.path(runtime_root, "xdg-config")
  bfc_dir <- file.path(runtime_root, "biocfilecache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(bfc_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(
    XDG_CACHE_HOME = cache_dir,
    XDG_CONFIG_HOME = config_dir,
    BFC_CACHE = bfc_dir
  )

  # Re-seed the now-isolated BiocFileCache with the real config captured
  # above, so ORFik::config() inside this isolated environment still
  # resolves to the correct project paths instead of a wrong default.
  if (!is.null(orfik_config) && length(orfik_config) == 4) {
    conf_dt <- data.frame(
      type = names(orfik_config),
      directory = unname(orfik_config),
      stringsAsFactors = FALSE
    )
    tryCatch(ORFik:::config.save(conf = conf_dt), error = function(e) NULL)
  }
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    BiocParallel::register(
      BiocParallel::SerialParam(progressbar = FALSE),
      default = TRUE
    )
  }
}

source_pipeline_cache <- function(analysis_dir) {
  cache_file <- file.path(analysis_dir, "scripts", "dominant_pipeline_cache.R")
  if (!file.exists(cache_file)) {
    stop("Missing pipeline cache helper: ", cache_file, call. = FALSE)
  }
  source(cache_file)
}

normalize_rdg_plot_mode <- function(mode) {
  mode <- tolower(trimws(as.character(mode)[1]))
  if (!nzchar(mode) || is.na(mode)) mode <- "all"
  if (mode %in% c("all", "full", "plots", "plot")) return("all")
  if (mode %in% c("tables", "table", "none", "skip", "skip_plots",
                  "no_plots", "no-plot", "no-plots")) {
    return("tables")
  }
  stop("Unknown RDG plot mode: ", mode,
       ". Use 'all' or 'tables'.", call. = FALSE)
}

pipeline_runtime_option_file <- function(analysis_dir, name) {
  file.path(analysis_dir, ".runtime", paste0(name, ".txt"))
}

write_pipeline_runtime_options <- function(analysis_dir, rdg_plot_mode) {
  runtime_dir <- file.path(analysis_dir, ".runtime")
  dir.create(runtime_dir, recursive = TRUE, showWarnings = FALSE)
  option_file <- pipeline_runtime_option_file(analysis_dir, "rdg_plot_mode")
  writeLines(paste0("rdg_plot_mode=", rdg_plot_mode), option_file)
  option_file
}

pipeline_steps <- function(analysis_dir) {
  data.frame(
    step = c(
      "manual_translons",
      "fst_pages",
      "clean_cds",
      "score_audit",
      "state_activity_metadata",
      "auxiliary_state_screen",
      "adjusted_enrichment",
      "control_aware",
      "broad_metadata",
      "glmnet",
      "latent",
      "known_uorf",
      "uorf_regulation",
      "score_semantics",
      "score_semantics_design",
      "rdg_graphs",
      "uorf_structure",
      "next_model",
      "branch_usage",
      "branch_coverage_qc",
      "grouped_branch_usage",
      "branch_allocation",
      "ouorf_cds_coupling",
      "downstream_tis_rescue",
      "fatigue_summary",
      "hierarchical_model",
      "hierarchical_residual_audit",
      "hierarchical_context_calibration",
      "hierarchical_context_model",
      "hierarchical_branch_covariates",
      "graph_geometry_model",
      "mechanistic_motif_model",
      "motif_prior_calibration",
      "joint_branch_allocation",
      "hierarchical_joint_allocation",
      "dirichlet_multinomial_allocation",
      "model_agreement",
      "model_agreement_loo",
      "model_agreement_transport",
      "context_interactions",
      "context_partial_pooling",
      "rdg_atlas_cards",
      "pooling_bias_audit",
      "auxiliary_state_rdg_audit",
      "allocation_output_disagreements",
      "therapeutic_hypotheses",
      "therapeutic_perturbation_predictions",
      "rdg_consensus",
      "rdg_validation_dossiers",
      "rdg_evidence_gaps",
      "merged_replicate_translon_shift",
      "rdg_review_batches",
      "rdg_review_links",
      "rdg_review_feedback",
      "rdg_validation_labels",
      "cds_frame_bumpiness_qc",
      "infected_off_frame_scan",
      "infected_off_frame_frame_risk",
      "infected_off_frame_hotspots",
      "infected_off_frame_review",
      "infected_off_frame_frame_diagnostics",
      "infected_off_frame_readlength_diagnostics"
    ),
    script = file.path(
      analysis_dir,
      "scripts",
      c(
        "dominant_manual_translons.R",
        "dominant_cell_state_required_fst_pages.R",
        "dominant_cell_states_fst_clean_cds.R",
        "dominant_state_score_audit.R",
        "dominant_state_activity_metadata_fst.R",
        "dominant_state_auxiliary_module_screen.R",
        "dominant_cell_state_adjusted_enrichment.R",
        "dominant_cell_state_control_aware_enrichment.R",
        "dominant_state_broad_metadata_model.R",
        "dominant_state_glmnet_metadata_model.R",
        "dominant_state_latent_program_summary.R",
        "known_uorf_effects_with_citations.R",
        "dominant_uorf_regulation_inference.R",
        "dominant_state_score_semantics_comparison.R",
        "dominant_state_score_semantics_design_aware.R",
        "dominant_ribosome_decision_graphs.R",
        "dominant_uorf_structure_rules.R",
        "dominant_next_model_candidate_discovery.R",
        "dominant_rdg_branch_usage_model.R",
        "dominant_rdg_branch_coverage_qc.R",
        "dominant_rdg_grouped_branch_usage_model.R",
        "dominant_rdg_count_aware_branch_allocation.R",
        "dominant_ouorf_cds_coupling_analysis.R",
        "dominant_downstream_tis_rescue_audit.R",
        "postviral_fatigue_pipeline_summary.R",
        "dominant_rdg_hierarchical_effect_model.R",
        "dominant_rdg_hierarchical_residual_audit.R",
        "dominant_rdg_hierarchical_context_calibration.R",
        "dominant_rdg_hierarchical_context_model.R",
        "dominant_rdg_hierarchical_branch_covariate_model.R",
        "dominant_rdg_graph_geometry_generalization_model.R",
        "dominant_rdg_mechanistic_motif_model.R",
        "dominant_rdg_motif_prior_branch_calibration.R",
        "dominant_rdg_joint_branch_allocation_model.R",
        "dominant_rdg_hierarchical_joint_allocation_model.R",
        "dominant_rdg_dirichlet_multinomial_branch_model.R",
        "dominant_rdg_model_agreement.R",
        "dominant_rdg_model_agreement_loo.R",
        "dominant_rdg_model_agreement_transportability.R",
        "dominant_rdg_context_interaction_model.R",
        "dominant_rdg_context_partial_pooling_model.R",
        "dominant_rdg_atlas_cards.R",
        "dominant_rdg_pooling_bias_audit.R",
        "dominant_rdg_auxiliary_state_audit.R",
        "dominant_rdg_allocation_output_disagreement_audit.R",
        "dominant_rdg_therapeutic_hypotheses.R",
        "dominant_rdg_therapeutic_perturbation_predictions.R",
        "dominant_rdg_consensus_atlas.R",
        "dominant_rdg_validation_dossiers.R",
        "dominant_rdg_evidence_gap_prioritization.R",
        "dominant_rdg_merged_replicate_translon_shift.R",
        "dominant_rdg_active_review_batches.R",
        "dominant_rdg_review_batch_links.R",
        "dominant_rdg_review_feedback_ingest.R",
        "dominant_rdg_validation_label_bridge.R",
        "dominant_cds_frame_bumpiness_qc.R",
        "dominant_infected_off_frame_cds_scan.R",
        "dominant_infected_off_frame_frame_risk_hotspots.R",
        "dominant_infected_off_frame_hotspots.R",
        "dominant_infected_off_frame_review_queue.R",
        "dominant_infected_off_frame_frame_diagnostics.R",
        "dominant_b2m_readlength_frame_diagnostics.R"
      )
    ),
    stringsAsFactors = FALSE
  )
}

pipeline_step_outputs <- function(analysis_dir) {
  results_dir <- analysis_results_dir(analysis_dir, create = TRUE)
  p <- function(...) file.path(results_dir, ...)
  list(
    # dominant_manual_translons.R deliberately writes under
    # <analysis_dir>/manual_translons/ (not results/) because these outputs
    # are curated/regenerable-from-a-tracked-queue, not purely disposable --
    # translons_to_be_added_to_manual.txt is git-tracked from that same
    # directory. Using results_dir here (like every other step) made this
    # step a permanent cache miss, since its real outputs were never where
    # the cache checked.
    manual_translons = c(
      file.path(analysis_dir, "manual_translons", "manual_translons_manifest.csv"),
      file.path(analysis_dir, "manual_translons", "manual_translons_ranges.rds"),
      file.path(analysis_dir, "manual_translons", "manual_translon_candidates.csv")
    ),
    fst_pages = c(
      p("required_fst_pages.csv"),
      p("required_fst_pages_by_gene.csv"),
      p("required_fst_selected_transcripts.csv"),
      p("required_fst_pages.txt"),
      p("required_fst_page_basenames.txt"),
      p("required_fst_pages_summary.txt")
    ),
    clean_cds = c(
      p("human_dominant_cell_states_clean_cds.csv"),
      p("human_dominant_cell_states_clean_cds_gene_expression.csv"),
      p("human_dominant_cell_states_clean_cds_gene_diagnostics.csv"),
      p("human_dominant_cell_states_cds_exon_coverage_qc.csv")
    ),
    score_audit = c(
      p("dominant_state_score_audit",
        "dominant_state_score_audit_sample_scores_long.csv"),
      p("dominant_state_score_audit",
        "dominant_state_score_audit_calibrated_annotations.csv"),
      p("dominant_state_score_audit",
        "dominant_state_score_audit_state_summary.csv"),
      p("dominant_state_score_audit",
        "dominant_state_score_audit_winner_counts.csv"),
      p("dominant_state_score_audit",
        "dominant_state_score_audit_oxphos_summary.csv")
    ),
    state_activity_metadata = c(
      p("human_dominant_state_activity_metadata.fst"),
      p("human_dominant_state_activity_metadata_long.csv"),
      p("human_dominant_state_activity_metadata_manifest.csv")
    ),
    auxiliary_state_screen = c(
      p("dominant_state_auxiliary_modules",
        "dominant_state_auxiliary_module_genes.csv"),
      p("dominant_state_auxiliary_modules",
        "dominant_state_auxiliary_module_scores.csv"),
      p("dominant_state_auxiliary_modules",
        "dominant_state_auxiliary_module_scores_wide.csv"),
      p("dominant_state_auxiliary_modules",
        "dominant_state_auxiliary_existing_state_correlations.csv"),
      p("dominant_state_auxiliary_modules",
        "dominant_state_auxiliary_module_summary.csv"),
      p("dominant_state_auxiliary_modules",
        "dominant_state_auxiliary_metadata_enrichment.csv"),
      p("dominant_state_auxiliary_modules",
        "dominant_state_auxiliary_numeric_metadata_correlations.csv"),
      p("dominant_state_auxiliary_modules",
        "dominant_state_auxiliary_rdg_gene_links.csv"),
      p("dominant_state_auxiliary_modules",
        "dominant_state_auxiliary_summary_metrics.csv"),
      p("dominant_state_auxiliary_modules", "figures",
        "auxiliary_state_novelty_and_availability.png"),
      p("dominant_state_auxiliary_modules", "figures",
        "auxiliary_state_existing_state_correlations.png"),
      p("dominant_state_auxiliary_modules", "figures",
        "auxiliary_state_top_metadata_effects.png")
    ),
    adjusted_enrichment = c(
      p("dominant_state_clean_cds_enrichment_all_main.csv"),
      p("dominant_state_clean_cds_enrichment_interactions.csv")
    ),
    control_aware = c(
      p("dominant_state_control_aware_condition_contrasts.csv"),
      p("dominant_state_control_aware_condition_summary.csv"),
      p("dominant_state_control_aware_enrichment_review.csv")
    ),
    broad_metadata = c(
      p("dominant_state_broad_metadata_main_terms.csv"),
      p("dominant_state_broad_metadata_nwise_terms.csv"),
      p("dominant_state_broad_metadata_rf_performance.csv"),
      p("dominant_state_broad_metadata_rdg_terms.csv")
    ),
    glmnet = c(
      p("dominant_state_glmnet_metadata_performance.csv"),
      p("dominant_state_glmnet_metadata_selected_term_scores.csv"),
      p("dominant_state_glmnet_metadata_rdg_terms.csv")
    ),
    latent = c(
      p("dominant_state_latent_program_pca_variance.csv"),
      p("dominant_state_latent_program_compact_report.csv")
    ),
    known_uorf = c(
      p("curated_inputs", "known_uorf_effects_with_citations.csv")
    ),
    uorf_regulation = c(
      p("dominant_uorf_feature_expression.csv"),
      p("dominant_uorf_regulation_sample_model.csv"),
      p("dominant_uorf_regulation_term_model.csv"),
      p("dominant_uorf_regulation_edges.csv"),
      p("dominant_uorf_regulation_allocation_model.csv")
    ),
    score_semantics = c(
      p("dominant_state_score_semantics",
        "dominant_state_score_semantics_label_summary.csv"),
      p("dominant_state_score_semantics",
        "dominant_state_score_semantics_metadata_call_enrichment.csv"),
      p("dominant_state_score_semantics",
        "dominant_state_score_semantics_metadata_state_summary.csv"),
      p("dominant_state_score_semantics",
        "dominant_state_score_semantics_rdg_clean_cds_contrasts.csv"),
      p("dominant_state_score_semantics",
        "dominant_state_score_semantics_rdg_raw_vs_percentile.csv"),
      p("dominant_state_score_semantics",
        "dominant_state_score_semantics_rdg_state_summary.csv"),
      p("dominant_state_score_semantics",
        "dominant_state_score_semantics_summary_metrics.csv")
    ),
    score_semantics_design = c(
      p("dominant_state_score_semantics_design_aware",
        "dominant_state_score_semantics_label_design_summary.csv"),
      p("dominant_state_score_semantics_design_aware",
        "dominant_state_score_semantics_metadata_fixed_effects.csv"),
      p("dominant_state_score_semantics_design_aware",
        "dominant_state_score_semantics_rdg_study_clean_cds_contrasts.csv"),
      p("dominant_state_score_semantics_design_aware",
        "dominant_state_score_semantics_rdg_clean_cds_fixed_effects.csv"),
      p("dominant_state_score_semantics_design_aware",
        "dominant_state_score_semantics_design_summary_metrics.csv")
    ),
    rdg_graphs = c(
      p("dominant_rdg_outputs", "rdg_feature_annotation.csv"),
      p("dominant_rdg_outputs", "rdg_branch_feature_contrasts.csv"),
      p("dominant_rdg_outputs", "rdg_branch_terms.csv"),
      p("dominant_rdg_outputs", "relative_usage_matrices",
        "relative_usage_in_long.csv")
    ),
    uorf_structure = c(
      p("dominant_uorf_structure_rules", "uorf_structure_gene_geometry.csv"),
      p("dominant_uorf_structure_rules",
        "sample_gene_clean_cds_usage_with_structure.csv"),
      p("dominant_uorf_structure_rules",
        "branch_gene_clean_cds_usage_with_structure.csv"),
      p("dominant_uorf_structure_rules",
        "uorf_structure_multivariable_rule_models.csv")
    ),
    next_model = c(
      p("dominant_next_model", "dominant_next_model_candidate_rankings_v4.csv"),
      p("dominant_next_model",
        "dominant_next_model_cds_buffering_effect_sizes_relevant_rdgs.csv"),
      p("dominant_next_model",
        "dominant_next_model_matched_control_cds_buffering.csv"),
      p("dominant_next_model",
        "dominant_next_model_condition_design_family_audit.csv"),
      p("dominant_next_model",
        "dominant_next_model_generic_condition_evidence.csv"),
      p("dominant_next_model",
        "dominant_next_model_internal_orf_matched_bump_evidence.csv"),
      p("dominant_next_model",
        "dominant_next_model_internal_orf_initiation_evidence.csv"),
      p("dominant_next_model",
        "dominant_next_model_internal_orf_architecture.csv"),
      p("dominant_next_model",
        "dominant_next_model_feature_confidence.csv"),
      p("dominant_next_model",
        "dominant_next_model_browser_review_handoff.csv"),
      p("dominant_next_model", "dominant_next_model_summary_metrics.csv")
    ),
    branch_usage = c(
      p("dominant_rdg_branch_usage",
        "dominant_rdg_branch_usage_contrast_table.csv"),
      p("dominant_rdg_branch_usage",
        "dominant_rdg_branch_usage_feature_summary.csv"),
      p("dominant_rdg_branch_usage",
        "dominant_rdg_branch_usage_gene_design_summary.csv"),
      p("dominant_rdg_branch_usage",
        "dominant_rdg_branch_usage_review.csv"),
      p("dominant_rdg_branch_usage",
        "dominant_rdg_branch_usage_summary_metrics.csv"),
      p("dominant_rdg_branch_usage", "figures",
        "dominant_rdg_branch_usage_top_shifts.png"),
      p("dominant_rdg_branch_usage", "figures",
        "dominant_rdg_branch_usage_viral_review.png")
    ),
    branch_coverage_qc = c(
      p("dominant_rdg_branch_coverage_qc",
        "branch_coverage_metadata_match_summary.csv"),
      p("dominant_rdg_branch_coverage_qc",
        "branch_coverage_gate_summary.csv"),
      p("dominant_rdg_branch_coverage_qc",
        "branch_coverage_feature_summary.csv"),
      p("dominant_rdg_branch_coverage_qc",
        "branch_coverage_gene_summary.csv"),
      p("dominant_rdg_branch_coverage_qc",
        "branch_coverage_unit_summary.csv"),
      p("dominant_rdg_branch_coverage_qc",
        "branch_coverage_replicate_group_summary.csv"),
      p("dominant_rdg_branch_coverage_qc",
        "branch_coverage_small_uorf_single_vs_merged.csv"),
      p("dominant_rdg_branch_coverage_qc",
        "branch_coverage_recommendation.csv"),
      p("dominant_rdg_branch_coverage_qc",
        "branch_coverage_summary_metrics.csv"),
      p("dominant_rdg_branch_coverage_qc", "figures",
        "branch_coverage_count_gate_by_length.png"),
      p("dominant_rdg_branch_coverage_qc", "figures",
        "branch_coverage_small_uorf_top_features.png")
    ),
    grouped_branch_usage = c(
      p("dominant_rdg_grouped_branch_usage",
        "dominant_rdg_grouped_branch_usage_adapter.csv"),
      p("dominant_rdg_grouped_branch_usage",
        "dominant_rdg_grouped_branch_usage_contrast_table.csv"),
      p("dominant_rdg_grouped_branch_usage",
        "dominant_rdg_grouped_branch_usage_feature_summary.csv"),
      p("dominant_rdg_grouped_branch_usage",
        "dominant_rdg_grouped_branch_usage_gene_design_summary.csv"),
      p("dominant_rdg_grouped_branch_usage",
        "dominant_rdg_grouped_branch_usage_review.csv"),
      p("dominant_rdg_grouped_branch_usage",
        "dominant_rdg_grouped_branch_usage_summary_metrics.csv"),
      p("dominant_rdg_grouped_branch_usage", "figures",
        "dominant_rdg_grouped_branch_usage_top_shifts.png"),
      p("dominant_rdg_grouped_branch_usage", "figures",
        "dominant_rdg_grouped_branch_usage_viral_review.png")
    ),
    branch_allocation = c(
      p("dominant_rdg_branch_allocation",
        "dominant_rdg_branch_allocation_posterior_contrasts.csv"),
      p("dominant_rdg_branch_allocation",
        "dominant_rdg_branch_allocation_feature_summary.csv"),
      p("dominant_rdg_branch_allocation",
        "dominant_rdg_branch_allocation_gene_design_summary.csv"),
      p("dominant_rdg_branch_allocation",
        "dominant_rdg_branch_allocation_review.csv"),
      p("dominant_rdg_branch_allocation",
        "dominant_rdg_branch_allocation_summary_metrics.csv"),
      p("dominant_rdg_branch_allocation", "figures",
        "dominant_rdg_branch_allocation_top_shifts.png"),
      p("dominant_rdg_branch_allocation", "figures",
        "dominant_rdg_branch_allocation_viral_review.png")
    ),
    ouorf_cds_coupling = c(
      p("dominant_ouorf_cds_coupling",
        "dominant_ouorf_cds_coupling_gene_summary.csv"),
      p("dominant_ouorf_cds_coupling",
        "dominant_ouorf_cds_coupling_candidate_review.csv"),
      p("dominant_ouorf_cds_coupling",
        "dominant_ouorf_cds_coupling_matched_summary.csv"),
      p("dominant_ouorf_cds_coupling",
        "dominant_ouorf_cds_coupling_onoff_pairs.csv"),
      p("dominant_ouorf_cds_coupling",
        "dominant_ouorf_cds_coupling_onoff_summary.csv"),
      p("dominant_ouorf_cds_coupling",
        "dominant_ouorf_cds_coupling_sparse_onoff_design_summary.csv"),
      p("dominant_ouorf_cds_coupling",
        "dominant_ouorf_cds_coupling_sparse_onoff_candidates.csv"),
      p("dominant_ouorf_cds_coupling",
        "dominant_ouorf_cds_coupling_summary_metrics.csv")
    ),
    downstream_tis_rescue = c(
      p("dominant_downstream_tis_rescue",
        "downstream_tis_rescue_candidates.csv"),
      p("dominant_downstream_tis_rescue",
        "downstream_tis_rescue_best_candidates.csv"),
      p("dominant_downstream_tis_rescue",
        "downstream_tis_rescue_feature_summary.csv"),
      p("dominant_downstream_tis_rescue",
        "downstream_tis_rescue_summary_metrics.csv")
    ),
    fatigue_summary = c(
      p("postviral_fatigue_outputs", "postviral_fatigue_summary_metrics.csv"),
      p("postviral_fatigue_outputs",
        "postviral_fatigue_rdg_feature_shift_summary.csv"),
      p("postviral_fatigue_outputs",
        "postviral_fatigue_cds_buffering_effects.csv"),
      p("postviral_fatigue_outputs",
        "postviral_fatigue_submodule_summary.csv"),
      p("postviral_fatigue_outputs",
        "postviral_fatigue_design_family_terms.csv"),
      p("postviral_fatigue_outputs",
        "postviral_fatigue_matched_design_family_summary.csv"),
      p("postviral_fatigue_outputs",
        "postviral_fatigue_iorf_evidence_classifier.csv"),
      p("postviral_fatigue_outputs",
        "postviral_fatigue_counterfactual_predictions.csv"),
      p("postviral_fatigue_outputs",
        "postviral_fatigue_matched_viral_ifn_review.csv"),
      p("postviral_fatigue_outputs",
        "postviral_fatigue_matched_design_cards.csv"),
      p("postviral_fatigue_outputs",
        "postviral_fatigue_matched_gene_design_evidence.csv"),
      p("postviral_fatigue_outputs",
        "postviral_fatigue_matched_gene_design_loo.csv"),
      p("postviral_fatigue_outputs",
        "postviral_fatigue_matched_gene_design_summary.csv"),
      p("postviral_fatigue_outputs",
        "postviral_fatigue_validation_queue.csv")
    ),
    hierarchical_model = c(
      p("dominant_rdg_hierarchical",
        "dominant_rdg_hierarchical_gene_design_effects.csv"),
      p("dominant_rdg_hierarchical",
        "dominant_rdg_hierarchical_counterfactual_predictions.csv"),
      p("dominant_rdg_hierarchical",
        "dominant_rdg_hierarchical_leave_one_study_calibration.csv"),
      p("dominant_rdg_hierarchical",
        "dominant_rdg_hierarchical_calibration_summary.csv"),
      p("dominant_rdg_hierarchical",
        "dominant_rdg_hierarchical_interval_calibration.csv"),
      p("dominant_rdg_hierarchical",
        "dominant_rdg_hierarchical_summary_metrics.csv"),
      p("dominant_rdg_hierarchical", "figures",
        "hierarchical_gene_design_heatmap.png"),
      p("dominant_rdg_hierarchical", "figures",
        "hierarchical_leave_one_study_calibration.png")
    ),
    hierarchical_residual_audit = c(
      p("dominant_rdg_hierarchical_residual_audit",
        "dominant_rdg_hierarchical_residual_audit_loo_context.csv"),
      p("dominant_rdg_hierarchical_residual_audit",
        "dominant_rdg_hierarchical_residual_audit_gene_design.csv"),
      p("dominant_rdg_hierarchical_residual_audit",
        "dominant_rdg_hierarchical_residual_audit_context_summary.csv"),
      p("dominant_rdg_hierarchical_residual_audit",
        "dominant_rdg_hierarchical_residual_audit_viral_contexts.csv"),
      p("dominant_rdg_hierarchical_residual_audit",
        "dominant_rdg_hierarchical_residual_audit_viral_gene_design.csv"),
      p("dominant_rdg_hierarchical_residual_audit",
        "dominant_rdg_hierarchical_residual_audit_summary_metrics.csv"),
      p("dominant_rdg_hierarchical_residual_audit", "figures",
        "hierarchical_residual_gene_design_scatter.png"),
      p("dominant_rdg_hierarchical_residual_audit", "figures",
        "hierarchical_residual_viral_contexts.png")
    ),
    hierarchical_context_calibration = c(
      p("dominant_rdg_hierarchical_context_calibration",
        "dominant_rdg_hierarchical_context_calibration_loo.csv"),
      p("dominant_rdg_hierarchical_context_calibration",
        "dominant_rdg_hierarchical_context_calibration_gene_design.csv"),
      p("dominant_rdg_hierarchical_context_calibration",
        "dominant_rdg_hierarchical_context_calibration_summary.csv"),
      p("dominant_rdg_hierarchical_context_calibration",
        "dominant_rdg_hierarchical_context_calibration_offset_summary.csv"),
      p("dominant_rdg_hierarchical_context_calibration",
        "dominant_rdg_hierarchical_context_calibration_viral_gene_design.csv"),
      p("dominant_rdg_hierarchical_context_calibration", "figures",
        "hierarchical_context_calibration_rmse.png"),
      p("dominant_rdg_hierarchical_context_calibration", "figures",
        "hierarchical_context_calibration_viral_genes.png")
    ),
    hierarchical_context_model = c(
      p("dominant_rdg_hierarchical_context_model",
        "dominant_rdg_hierarchical_context_model_loo.csv"),
      p("dominant_rdg_hierarchical_context_model",
        "dominant_rdg_hierarchical_context_model_gene_design.csv"),
      p("dominant_rdg_hierarchical_context_model",
        "dominant_rdg_hierarchical_context_model_summary.csv"),
      p("dominant_rdg_hierarchical_context_model",
        "dominant_rdg_hierarchical_context_model_viral_gene_design.csv"),
      p("dominant_rdg_hierarchical_context_model",
        "dominant_rdg_hierarchical_context_model_coefficients.csv"),
      p("dominant_rdg_hierarchical_context_model",
        "dominant_rdg_hierarchical_context_model_lambda_selection.csv"),
      p("dominant_rdg_hierarchical_context_model",
        "dominant_rdg_hierarchical_context_model_metadata.csv"),
      p("dominant_rdg_hierarchical_context_model", "figures",
        "hierarchical_context_model_rmse.png"),
      p("dominant_rdg_hierarchical_context_model", "figures",
        "hierarchical_context_model_viral_genes.png")
    ),
    hierarchical_branch_covariates = c(
      p("dominant_rdg_hierarchical_branch_covariates",
        "dominant_rdg_hierarchical_branch_covariate_model_summary.csv"),
      p("dominant_rdg_hierarchical_branch_covariates",
        "dominant_rdg_hierarchical_branch_covariate_model_loo.csv"),
      p("dominant_rdg_hierarchical_branch_covariates",
        "dominant_rdg_hierarchical_branch_covariate_model_gene_design.csv"),
      p("dominant_rdg_hierarchical_branch_covariates",
        "dominant_rdg_hierarchical_branch_covariate_model_viral_gene_design.csv"),
      p("dominant_rdg_hierarchical_branch_covariates",
        "dominant_rdg_hierarchical_branch_covariate_model_coefficients.csv"),
      p("dominant_rdg_hierarchical_branch_covariates",
        "dominant_rdg_hierarchical_branch_covariate_model_lambda_selection.csv"),
      p("dominant_rdg_hierarchical_branch_covariates",
        "dominant_rdg_hierarchical_branch_covariate_model_metadata.csv"),
      p("dominant_rdg_hierarchical_branch_covariates",
        "dominant_rdg_hierarchical_branch_covariates_training_covariates.csv"),
      p("dominant_rdg_hierarchical_branch_covariates", "figures",
        "branch_covariate_model_rmse.png"),
      p("dominant_rdg_hierarchical_branch_covariates", "figures",
        "branch_covariate_model_viral_genes.png")
    ),
    graph_geometry_model = c(
      p("dominant_rdg_graph_geometry_model",
        "dominant_rdg_graph_geometry_model_summary.csv"),
      p("dominant_rdg_graph_geometry_model",
        "dominant_rdg_graph_geometry_model_loo.csv"),
      p("dominant_rdg_graph_geometry_model",
        "dominant_rdg_graph_geometry_model_gene_design.csv"),
      p("dominant_rdg_graph_geometry_model",
        "dominant_rdg_graph_geometry_model_viral_gene_design.csv"),
      p("dominant_rdg_graph_geometry_model",
        "dominant_rdg_graph_geometry_model_covariates.csv"),
      p("dominant_rdg_graph_geometry_model",
        "dominant_rdg_graph_geometry_model_coefficients.csv"),
      p("dominant_rdg_graph_geometry_model",
        "dominant_rdg_graph_geometry_model_lambda_selection.csv"),
      p("dominant_rdg_graph_geometry_model",
        "dominant_rdg_graph_geometry_model_metadata.csv"),
      p("dominant_rdg_graph_geometry_model", "figures",
        "graph_geometry_model_rmse.png"),
      p("dominant_rdg_graph_geometry_model", "figures",
        "graph_geometry_model_viral_genes.png"),
      p("dominant_rdg_graph_geometry_model", "figures",
        "graph_geometry_model_coefficients.png")
    ),
    mechanistic_motif_model = c(
      p("dominant_rdg_mechanistic_motif_model",
        "dominant_rdg_mechanistic_motif_model_summary.csv"),
      p("dominant_rdg_mechanistic_motif_model",
        "dominant_rdg_mechanistic_motif_model_loo.csv"),
      p("dominant_rdg_mechanistic_motif_model",
        "dominant_rdg_mechanistic_motif_model_gene_design.csv"),
      p("dominant_rdg_mechanistic_motif_model",
        "dominant_rdg_mechanistic_motif_model_viral_gene_design.csv"),
      p("dominant_rdg_mechanistic_motif_model",
        "dominant_rdg_mechanistic_motif_scores.csv"),
      p("dominant_rdg_mechanistic_motif_model",
        "dominant_rdg_mechanistic_motif_model_coefficients.csv"),
      p("dominant_rdg_mechanistic_motif_model",
        "dominant_rdg_mechanistic_motif_model_lambda_selection.csv"),
      p("dominant_rdg_mechanistic_motif_model",
        "dominant_rdg_mechanistic_motif_model_metadata.csv"),
      p("dominant_rdg_mechanistic_motif_model", "figures",
        "mechanistic_motif_model_rmse.png"),
      p("dominant_rdg_mechanistic_motif_model", "figures",
        "mechanistic_motif_model_viral_genes.png"),
      p("dominant_rdg_mechanistic_motif_model", "figures",
        "mechanistic_motif_model_coefficients.png")
    ),
    motif_prior_calibration = c(
      p("dominant_rdg_motif_prior_calibration",
        "dominant_rdg_motif_prior_calibration_summary.csv"),
      p("dominant_rdg_motif_prior_calibration",
        "dominant_rdg_motif_prior_branch_rows.csv"),
      p("dominant_rdg_motif_prior_calibration",
        "dominant_rdg_motif_prior_gene_design_calibration.csv"),
      p("dominant_rdg_motif_prior_calibration",
        "dominant_rdg_motif_prior_viral_review.csv"),
      p("dominant_rdg_motif_prior_calibration",
        "dominant_rdg_motif_prior_review_queue.csv"),
      p("dominant_rdg_motif_prior_calibration", "figures",
        "motif_prior_vs_branch_evidence.png"),
      p("dominant_rdg_motif_prior_calibration", "figures",
        "motif_prior_viral_alignment.png")
    ),
    joint_branch_allocation = c(
      p("dominant_rdg_joint_branch_allocation",
        "dominant_rdg_joint_branch_allocation_rows.csv"),
      p("dominant_rdg_joint_branch_allocation",
        "dominant_rdg_joint_branch_allocation_class_rows.csv"),
      p("dominant_rdg_joint_branch_allocation",
        "dominant_rdg_joint_branch_allocation_gene_design_summary.csv"),
      p("dominant_rdg_joint_branch_allocation",
        "dominant_rdg_joint_branch_allocation_review.csv"),
      p("dominant_rdg_joint_branch_allocation",
        "dominant_rdg_joint_branch_allocation_viral_review.csv"),
      p("dominant_rdg_joint_branch_allocation",
        "dominant_rdg_joint_branch_allocation_summary_metrics.csv"),
      p("dominant_rdg_joint_branch_allocation", "figures",
        "joint_branch_allocation_top_shifts.png"),
      p("dominant_rdg_joint_branch_allocation", "figures",
        "joint_branch_allocation_viral_review.png"),
      p("dominant_rdg_joint_branch_allocation", "figures",
        "joint_branch_allocation_delta_heatmap.png")
    ),
    hierarchical_joint_allocation = c(
      p("dominant_rdg_hierarchical_joint_allocation",
        "dominant_rdg_hierarchical_joint_allocation_branch_effects.csv"),
      p("dominant_rdg_hierarchical_joint_allocation",
        "dominant_rdg_hierarchical_joint_allocation_study_effects.csv"),
      p("dominant_rdg_hierarchical_joint_allocation",
        "dominant_rdg_hierarchical_joint_allocation_gene_design_summary.csv"),
      p("dominant_rdg_hierarchical_joint_allocation",
        "dominant_rdg_hierarchical_joint_allocation_review.csv"),
      p("dominant_rdg_hierarchical_joint_allocation",
        "dominant_rdg_hierarchical_joint_allocation_viral_review.csv"),
      p("dominant_rdg_hierarchical_joint_allocation",
        "dominant_rdg_hierarchical_joint_allocation_summary_metrics.csv"),
      p("dominant_rdg_hierarchical_joint_allocation", "figures",
        "hierarchical_joint_allocation_top_effects.png"),
      p("dominant_rdg_hierarchical_joint_allocation", "figures",
        "hierarchical_joint_allocation_viral_review.png"),
      p("dominant_rdg_hierarchical_joint_allocation", "figures",
        "hierarchical_joint_allocation_delta_heatmap.png")
    ),
    dirichlet_multinomial_allocation = c(
      p("dominant_rdg_dirichlet_multinomial",
        "dominant_rdg_dirichlet_multinomial_branch_rows.csv"),
      p("dominant_rdg_dirichlet_multinomial",
        "dominant_rdg_dirichlet_multinomial_study_effects.csv"),
      p("dominant_rdg_dirichlet_multinomial",
        "dominant_rdg_dirichlet_multinomial_branch_effects.csv"),
      p("dominant_rdg_dirichlet_multinomial",
        "dominant_rdg_dirichlet_multinomial_gene_design_summary.csv"),
      p("dominant_rdg_dirichlet_multinomial",
        "dominant_rdg_dirichlet_multinomial_review.csv"),
      p("dominant_rdg_dirichlet_multinomial",
        "dominant_rdg_dirichlet_multinomial_viral_review.csv"),
      p("dominant_rdg_dirichlet_multinomial",
        "dominant_rdg_dirichlet_multinomial_summary_metrics.csv"),
      p("dominant_rdg_dirichlet_multinomial", "figures",
        "dirichlet_multinomial_top_effects.png"),
      p("dominant_rdg_dirichlet_multinomial", "figures",
        "dirichlet_multinomial_viral_review.png"),
      p("dominant_rdg_dirichlet_multinomial", "figures",
        "dirichlet_multinomial_delta_heatmap.png")
    ),
    model_agreement = c(
      p("dominant_rdg_model_agreement",
        "dominant_rdg_model_agreement_branch_comparison.csv"),
      p("dominant_rdg_model_agreement",
        "dominant_rdg_model_agreement_branch_summary.csv"),
      p("dominant_rdg_model_agreement",
        "dominant_rdg_model_agreement_gene_design.csv"),
      p("dominant_rdg_model_agreement",
        "dominant_rdg_model_agreement_review.csv"),
      p("dominant_rdg_model_agreement",
        "dominant_rdg_model_agreement_viral_review.csv"),
      p("dominant_rdg_model_agreement",
        "dominant_rdg_model_agreement_summary_metrics.csv"),
      p("dominant_rdg_model_agreement", "figures",
        "model_agreement_branch_comparison.png"),
      p("dominant_rdg_model_agreement", "figures",
        "model_agreement_viral_matrix.png")
    ),
    model_agreement_loo = c(
      p("dominant_rdg_model_agreement_loo",
        "dominant_rdg_model_agreement_study_branch_evidence.csv"),
      p("dominant_rdg_model_agreement_loo",
        "dominant_rdg_model_agreement_study_evidence.csv"),
      p("dominant_rdg_model_agreement_loo",
        "dominant_rdg_model_agreement_loo_rows.csv"),
      p("dominant_rdg_model_agreement_loo",
        "dominant_rdg_model_agreement_loo_branch_summary.csv"),
      p("dominant_rdg_model_agreement_loo",
        "dominant_rdg_model_agreement_loo_gene_design.csv"),
      p("dominant_rdg_model_agreement_loo",
        "dominant_rdg_model_agreement_loo_viral_review.csv"),
      p("dominant_rdg_model_agreement_loo",
        "dominant_rdg_model_agreement_loo_fragile_review.csv"),
      p("dominant_rdg_model_agreement_loo",
        "dominant_rdg_model_agreement_loo_summary_metrics.csv"),
      p("dominant_rdg_model_agreement_loo", "figures",
        "model_agreement_loo_consensus_stability.png"),
      p("dominant_rdg_model_agreement_loo", "figures",
        "model_agreement_loo_viral_stability.png")
    ),
    model_agreement_transport = c(
      p("dominant_rdg_model_agreement_transport",
        "dominant_rdg_model_agreement_transport_study_metadata.csv"),
      p("dominant_rdg_model_agreement_transport",
        "dominant_rdg_model_agreement_transport_study_groups.csv"),
      p("dominant_rdg_model_agreement_transport",
        "dominant_rdg_model_agreement_transport_omission_rows.csv"),
      p("dominant_rdg_model_agreement_transport",
        "dominant_rdg_model_agreement_transport_axis_model_summary.csv"),
      p("dominant_rdg_model_agreement_transport",
        "dominant_rdg_model_agreement_transport_branch_axis_summary.csv"),
      p("dominant_rdg_model_agreement_transport",
        "dominant_rdg_model_agreement_transport_branch_summary.csv"),
      p("dominant_rdg_model_agreement_transport",
        "dominant_rdg_model_agreement_transport_gene_axis_summary.csv"),
      p("dominant_rdg_model_agreement_transport",
        "dominant_rdg_model_agreement_transport_gene_design.csv"),
      p("dominant_rdg_model_agreement_transport",
        "dominant_rdg_model_agreement_transport_viral_review.csv"),
      p("dominant_rdg_model_agreement_transport",
        "dominant_rdg_model_agreement_transport_fragile_review.csv"),
      p("dominant_rdg_model_agreement_transport",
        "dominant_rdg_model_agreement_transport_summary_metrics.csv"),
      p("dominant_rdg_model_agreement_transport", "figures",
        "model_agreement_transport_consensus.png"),
      p("dominant_rdg_model_agreement_transport", "figures",
        "model_agreement_transport_viral_axes.png")
    ),
    context_interactions = c(
      p("dominant_rdg_context_interactions",
        "dominant_rdg_context_interaction_study_context_effects.csv"),
      p("dominant_rdg_context_interactions",
        "dominant_rdg_context_interaction_branch_context_effects.csv"),
      p("dominant_rdg_context_interactions",
        "dominant_rdg_context_interaction_review.csv"),
      p("dominant_rdg_context_interactions",
        "dominant_rdg_context_interaction_gene_design_summary.csv"),
      p("dominant_rdg_context_interactions",
        "dominant_rdg_context_interaction_viral_consensus_summary.csv"),
      p("dominant_rdg_context_interactions",
        "dominant_rdg_context_interaction_summary_metrics.csv"),
      p("dominant_rdg_context_interactions", "figures",
        "context_interaction_viral_top_contexts.png"),
      p("dominant_rdg_context_interactions", "figures",
        "context_interaction_ifih1_tissue.png"),
      p("dominant_rdg_context_interactions", "figures",
        "context_interaction_viral_consensus_tissue_heatmap.png")
    ),
    context_partial_pooling = c(
      p("dominant_rdg_context_partial_pooling",
        "dominant_rdg_context_partial_pooling_rows.csv"),
      p("dominant_rdg_context_partial_pooling",
        "dominant_rdg_context_partial_pooling_gene_design_summary.csv"),
      p("dominant_rdg_context_partial_pooling",
        "dominant_rdg_context_partial_pooling_viral_review.csv"),
      p("dominant_rdg_context_partial_pooling",
        "dominant_rdg_context_partial_pooling_raw_vs_pooled.csv"),
      p("dominant_rdg_context_partial_pooling",
        "dominant_rdg_context_partial_pooling_summary_metrics.csv"),
      p("dominant_rdg_context_partial_pooling", "figures",
        "context_partial_pooling_viral_top_effects.png"),
      p("dominant_rdg_context_partial_pooling", "figures",
        "context_partial_pooling_viral_top_effects.html"),
      p("dominant_rdg_context_partial_pooling", "figures",
        "context_partial_pooling_ifih1_lung.png")
    ),
    pooling_bias_audit = c(
      p("dominant_rdg_pooling_bias_audit",
        "pooling_claim_language.csv"),
      p("dominant_rdg_pooling_bias_audit",
        "normalization_limits.csv"),
      p("dominant_rdg_pooling_bias_audit",
        "pooling_gain_loss_summary.csv"),
      p("dominant_rdg_pooling_bias_audit",
        "pooling_contrast_tier_summary.csv"),
      p("dominant_rdg_pooling_bias_audit",
        "grouped_review_support_summary.csv"),
      p("dominant_rdg_pooling_bias_audit",
        "design_family_pooling_bias_summary.csv"),
      p("dominant_rdg_pooling_bias_audit",
        "protocol_context_bias_summary.csv"),
      p("dominant_rdg_pooling_bias_audit",
        "partial_pooling_tradeoff_summary.csv"),
      p("dominant_rdg_pooling_bias_audit",
        "partial_pooling_pool_loss_review.csv"),
      p("dominant_rdg_pooling_bias_audit",
        "replicate_pooling_loss_summary.csv"),
      p("dominant_rdg_pooling_bias_audit",
        "protocol_bias_run_summary.csv"),
      p("dominant_rdg_pooling_bias_audit",
        "condition_protocol_hypotheses.csv"),
      p("dominant_rdg_pooling_bias_audit",
        "pooling_bias_summary_metrics.csv"),
      p("dominant_rdg_pooling_bias_audit",
        "pooling_bias_report.md")
    ),
    rdg_atlas_cards = c(
      p("dominant_rdg_atlas", "dominant_rdg_atlas_cards.csv"),
      p("dominant_rdg_atlas", "dominant_rdg_atlas_branch_contexts.csv"),
      p("dominant_rdg_atlas", "dominant_rdg_atlas_cards_review.csv"),
      p("dominant_rdg_atlas", "dominant_rdg_atlas_summary_metrics.csv"),
      p("dominant_rdg_atlas",
        "dominant_rdg_atlas_browser_validation_queue.csv"),
      p("dominant_rdg_atlas",
        "dominant_rdg_atlas_viral_infection_review.csv"),
      p("dominant_rdg_atlas", "figures",
        "dominant_rdg_atlas_buffering_heatmap.png"),
      p("dominant_rdg_atlas", "figures",
        "dominant_rdg_atlas_confidence_vs_effect.png"),
      p("dominant_rdg_atlas", "figures",
        "dominant_rdg_atlas_viral_review_priority.png")
    ),
    auxiliary_state_rdg_audit = c(
      p("dominant_rdg_auxiliary_state_audit",
        "dominant_rdg_auxiliary_state_stratum_scores.csv"),
      p("dominant_rdg_auxiliary_state_audit",
        "dominant_rdg_auxiliary_state_branch_rows.csv"),
      p("dominant_rdg_auxiliary_state_audit",
        "dominant_rdg_auxiliary_state_gene_design_summary.csv"),
      p("dominant_rdg_auxiliary_state_audit",
        "dominant_rdg_auxiliary_state_atlas_annotations.csv"),
      p("dominant_rdg_auxiliary_state_audit",
        "dominant_rdg_auxiliary_state_selected_state_matrix.csv"),
      p("dominant_rdg_auxiliary_state_audit",
        "dominant_rdg_auxiliary_state_summary_by_state.csv"),
      p("dominant_rdg_auxiliary_state_audit",
        "dominant_rdg_auxiliary_state_summary_metrics.csv"),
      p("dominant_rdg_auxiliary_state_audit", "figures",
        "auxiliary_state_rdg_overlap_heatmap.png"),
      p("dominant_rdg_auxiliary_state_audit", "figures",
        "auxiliary_state_overlap_counts.png"),
      p("dominant_rdg_auxiliary_state_audit", "figures",
        "auxiliary_state_vs_atlas_effect.png")
    ),
    allocation_output_disagreements = c(
      p("dominant_rdg_allocation_output_disagreements",
        "dominant_rdg_allocation_output_disagreement_audit.csv"),
      p("dominant_rdg_allocation_output_disagreements",
        "dominant_rdg_allocation_output_disagreement_review.csv"),
      p("dominant_rdg_allocation_output_disagreements",
        "dominant_rdg_allocation_output_disagreement_contexts.csv"),
      p("dominant_rdg_allocation_output_disagreements",
        "dominant_rdg_allocation_output_disagreement_study_evidence.csv"),
      p("dominant_rdg_allocation_output_disagreements",
        "dominant_rdg_allocation_output_disagreement_summary_metrics.csv"),
      p("dominant_rdg_allocation_output_disagreements", "figures",
        "allocation_output_disagreement_scatter.png"),
      p("dominant_rdg_allocation_output_disagreements", "figures",
        "allocation_output_disagreement_branch_heatmap.png")
    ),
    therapeutic_hypotheses = c(
      p("dominant_rdg_therapeutic_hypotheses",
        "dominant_rdg_therapeutic_hypothesis_gene_evidence.csv"),
      p("dominant_rdg_therapeutic_hypotheses",
        "dominant_rdg_therapeutic_hypothesis_scores.csv"),
      p("dominant_rdg_therapeutic_hypotheses",
        "dominant_rdg_therapeutic_hypothesis_review.csv"),
      p("dominant_rdg_therapeutic_hypotheses",
        "dominant_rdg_therapeutic_hypothesis_summary_metrics.csv"),
      p("dominant_rdg_therapeutic_hypotheses", "figures",
        "therapeutic_hypothesis_scores.png"),
      p("dominant_rdg_therapeutic_hypotheses", "figures",
        "therapeutic_hypothesis_evidence_matrix.png")
    ),
    therapeutic_perturbation_predictions = c(
      p("dominant_rdg_therapeutic_perturbation_predictions",
        "dominant_rdg_therapeutic_perturbation_gene_predictions.csv"),
      p("dominant_rdg_therapeutic_perturbation_predictions",
        "dominant_rdg_therapeutic_perturbation_hypothesis_summary.csv"),
      p("dominant_rdg_therapeutic_perturbation_predictions",
        "dominant_rdg_therapeutic_perturbation_review.csv"),
      p("dominant_rdg_therapeutic_perturbation_predictions",
        "dominant_rdg_therapeutic_perturbation_summary_metrics.csv"),
      p("dominant_rdg_therapeutic_perturbation_predictions", "figures",
        "therapeutic_perturbation_top_gene_effects.png"),
      p("dominant_rdg_therapeutic_perturbation_predictions", "figures",
        "therapeutic_perturbation_hypothesis_summary.png")
    ),
    rdg_consensus = c(
      p("dominant_rdg_consensus", "dominant_rdg_consensus_candidates.csv"),
      p("dominant_rdg_consensus", "dominant_rdg_consensus_review_queue.csv"),
      p("dominant_rdg_consensus", "dominant_rdg_consensus_evidence_matrix.csv"),
      p("dominant_rdg_consensus", "dominant_rdg_consensus_gene_summary.csv"),
      p("dominant_rdg_consensus", "dominant_rdg_consensus_summary_metrics.csv"),
      p("dominant_rdg_consensus", "figures",
        "dominant_rdg_consensus_top_candidates.png"),
      p("dominant_rdg_consensus", "figures",
        "dominant_rdg_consensus_evidence_heatmap.png"),
      p("dominant_rdg_consensus", "figures",
        "dominant_rdg_consensus_class_counts.png")
    ),
    rdg_validation_dossiers = c(
      p("dominant_rdg_validation_dossiers",
        "dominant_rdg_validation_dossier_candidates.csv"),
      p("dominant_rdg_validation_dossiers",
        "dominant_rdg_validation_dossier_branch_contexts.csv"),
      p("dominant_rdg_validation_dossiers",
        "dominant_rdg_validation_dossier_manual_review_template.csv"),
      p("dominant_rdg_validation_dossiers",
        "dominant_rdg_validation_dossier_gene_summary.csv"),
      p("dominant_rdg_validation_dossiers",
        "dominant_rdg_validation_dossier_summary_metrics.csv"),
      p("dominant_rdg_validation_dossiers",
        "dominant_rdg_validation_dossier_report.md"),
      p("dominant_rdg_validation_dossiers", "figures",
        "rdg_validation_dossier_priority.png"),
      p("dominant_rdg_validation_dossiers", "figures",
        "rdg_validation_dossier_evidence_heatmap.png"),
      p("dominant_rdg_validation_dossiers", "figures",
        "rdg_validation_dossier_workload.png")
    ),
    rdg_evidence_gaps = c(
      p("dominant_rdg_evidence_gaps",
        "dominant_rdg_evidence_gap_priorities.csv"),
      p("dominant_rdg_evidence_gaps",
        "dominant_rdg_evidence_gap_gene_summary.csv"),
      p("dominant_rdg_evidence_gaps",
        "dominant_rdg_evidence_gap_action_summary.csv"),
      p("dominant_rdg_evidence_gaps",
        "dominant_rdg_evidence_gap_negative_controls.csv"),
      p("dominant_rdg_evidence_gaps",
        "dominant_rdg_evidence_gap_summary_metrics.csv"),
      p("dominant_rdg_evidence_gaps",
        "dominant_rdg_evidence_gap_report.md"),
      p("dominant_rdg_evidence_gaps", "figures",
        "evidence_gap_priority_map.png"),
      p("dominant_rdg_evidence_gaps", "figures",
        "evidence_gap_action_counts.png"),
      p("dominant_rdg_evidence_gaps", "figures",
        "evidence_gap_gene_action_heatmap.png")
    ),
    merged_replicate_translon_shift = c(
      p("dominant_rdg_merged_replicate_translon_shift",
        "merged_replicate_translon_shift_rows.csv"),
      p("dominant_rdg_merged_replicate_translon_shift",
        "merged_replicate_translon_shift_gene_design_summary.csv"),
      p("dominant_rdg_merged_replicate_translon_shift",
        "merged_replicate_translon_shift_review_queue.csv"),
      p("dominant_rdg_merged_replicate_translon_shift",
        "merged_replicate_translon_shift_candidate_queue.csv"),
      p("dominant_rdg_merged_replicate_translon_shift",
        "merged_replicate_translon_shift_pooling_summary.csv"),
      p("dominant_rdg_merged_replicate_translon_shift",
        "merged_replicate_translon_shift_summary_metrics.csv"),
      p("dominant_rdg_merged_replicate_translon_shift",
        "merged_replicate_translon_shift_report.md"),
      p("dominant_rdg_merged_replicate_translon_shift", "figures",
        "merged_replicate_translon_shift_clean_uorf_scatter.png"),
      p("dominant_rdg_merged_replicate_translon_shift", "figures",
        "merged_replicate_translon_shift_priority_top.png"),
      p("dominant_rdg_merged_replicate_translon_shift", "figures",
        "merged_replicate_translon_shift_candidate_top.png")
    ),
    rdg_review_batches = c(
      p("dominant_rdg_review_batches",
        "dominant_rdg_review_batches.csv"),
      p("dominant_rdg_review_batches",
        "dominant_rdg_review_batch_template.csv"),
      p("dominant_rdg_review_batches",
        "dominant_rdg_review_batch_summary.csv"),
      p("dominant_rdg_review_batches",
        "dominant_rdg_review_batch_gene_balance.csv"),
      p("dominant_rdg_review_batches",
        "dominant_rdg_review_batch_metrics.csv"),
      p("dominant_rdg_review_batches",
        "dominant_rdg_review_batch_report.md"),
      p("dominant_rdg_review_batches", "figures",
        "review_batch_action_counts.png"),
      p("dominant_rdg_review_batches", "figures",
        "review_batch_priority_readiness.png"),
      p("dominant_rdg_review_batches", "figures",
        "review_batch_gene_matrix.png")
    ),
    rdg_review_links = c(
      p("dominant_rdg_review_links",
        "dominant_rdg_review_batch_links.csv"),
      p("dominant_rdg_review_links",
        "dominant_rdg_review_batch_template_with_links.csv"),
      p("dominant_rdg_review_links",
        "dominant_rdg_review_link_summary.csv"),
      p("dominant_rdg_review_links",
        "dominant_rdg_review_cards.html"),
      p("dominant_rdg_review_links",
        "dominant_rdg_review_cards.md"),
      p("dominant_rdg_review_links",
        "dominant_rdg_review_link_report.md")
    ),
    rdg_review_feedback = c(
      p("dominant_rdg_review_feedback",
        "dominant_rdg_review_feedback_labels.csv"),
      p("dominant_rdg_review_feedback",
        "dominant_rdg_review_feedback_training_set.csv"),
      p("dominant_rdg_review_feedback",
        "dominant_rdg_review_feedback_browser_validation_notes_candidate.csv"),
      p("dominant_rdg_review_feedback",
        "dominant_rdg_review_feedback_summary.csv"),
      p("dominant_rdg_review_feedback",
        "dominant_rdg_review_feedback_batch_summary.csv"),
      p("dominant_rdg_review_feedback",
        "dominant_rdg_review_feedback_report.md")
    ),
    rdg_validation_labels = c(
      p("dominant_rdg_validation_label_bridge",
        "dominant_rdg_validation_label_set.csv"),
      p("dominant_rdg_validation_label_bridge",
        "dominant_rdg_validation_label_gene_design.csv"),
      p("dominant_rdg_validation_label_bridge",
        "dominant_rdg_validation_label_motif_overlay.csv"),
      p("dominant_rdg_validation_label_bridge",
        "dominant_rdg_validation_label_gap_queue.csv"),
      p("dominant_rdg_validation_label_bridge",
        "dominant_rdg_validation_label_balance.csv"),
      p("dominant_rdg_validation_label_bridge",
        "dominant_rdg_validation_label_bridge_summary.csv"),
      p("dominant_rdg_validation_label_bridge",
        "dominant_rdg_validation_label_bridge_report.md")
    ),
    cds_frame_bumpiness_qc = c(
      p("dominant_cds_frame_bumpiness_qc",
        "cds_frame_bumpiness_run_metrics.csv"),
      p("dominant_cds_frame_bumpiness_qc",
        "housekeeping_baseline_run_qc.csv"),
      p("dominant_cds_frame_bumpiness_qc",
        "myc_rrm1_infected_library_context_summary.csv"),
      p("dominant_cds_frame_bumpiness_qc",
        "myc_rrm1_viral_design_housekeeping_contrast_summary.csv"),
      p("dominant_cds_frame_bumpiness_qc",
        "target_gene_frame_bumpiness_summary.csv")
    ),
    infected_off_frame_scan = c(
      p("dominant_infected_off_frame_cds_scan",
        "infected_off_frame_cds_candidate_gene_scores.csv"),
      p("dominant_infected_off_frame_cds_scan",
        "infected_off_frame_cds_candidate_contexts.csv"),
      p("dominant_infected_off_frame_cds_scan",
        "infected_off_frame_cds_timepoint_contexts.csv"),
      p("dominant_infected_off_frame_cds_scan",
        "infected_off_frame_cds_top_candidate_runs.csv"),
      p("dominant_infected_off_frame_cds_scan",
        "infected_off_frame_cds_frame_risk_runs.csv"),
      p("dominant_infected_off_frame_cds_scan",
        "infected_off_frame_cds_frame_risk_gene_scores.csv"),
      p("dominant_infected_off_frame_cds_scan",
        "infected_off_frame_cds_summary_metrics.csv"),
      p("dominant_infected_off_frame_cds_scan", "figures",
        "infected_off_frame_cds_top_gene_scores.png"),
      p("dominant_infected_off_frame_cds_scan", "figures",
        "infected_off_frame_cds_context_matrix.png"),
      p("dominant_infected_off_frame_cds_scan", "figures",
        "infected_off_frame_cds_top_candidate_runs.png")
    ),
    infected_off_frame_frame_risk = c(
      p("dominant_infected_off_frame_frame_risk",
        "infected_off_frame_frame_risk_runs.csv"),
      p("dominant_infected_off_frame_frame_risk",
        "infected_off_frame_frame_risk_run_positions.csv"),
      p("dominant_infected_off_frame_frame_risk",
        "infected_off_frame_frame_risk_hotspots.csv"),
      p("dominant_infected_off_frame_frame_risk",
        "infected_off_frame_frame_risk_gene_summary.csv"),
      p("dominant_infected_off_frame_frame_risk",
        "infected_off_frame_frame_risk_summary_metrics.csv"),
      p("dominant_infected_off_frame_frame_risk", "figures",
        "infected_off_frame_frame_risk_hotspots.png")
    ),
    infected_off_frame_hotspots = c(
      p("dominant_infected_off_frame_hotspots",
        "infected_off_frame_hotspot_gene_summary.csv"),
      p("dominant_infected_off_frame_hotspots",
        "infected_off_frame_hotspot_recurrence.csv"),
      p("dominant_infected_off_frame_hotspots",
        "infected_off_frame_hotspot_run_positions.csv"),
      p("dominant_infected_off_frame_hotspots",
        "infected_off_frame_hotspot_aggregate_profiles.csv"),
      p("dominant_infected_off_frame_hotspots",
        "infected_off_frame_hotspot_group_run_summary.csv"),
      p("dominant_infected_off_frame_hotspots",
        "infected_off_frame_hotspot_summary_metrics.csv"),
      p("dominant_infected_off_frame_hotspots", "figures",
        "infected_off_frame_hotspot_recurrence.png"),
      p("dominant_infected_off_frame_hotspots", "figures",
        "infected_off_frame_hotspot_aggregate_profiles.png")
    ),
    infected_off_frame_review = c(
      p("dominant_infected_off_frame_review",
        "infected_off_frame_hotspot_review_queue.csv"),
      p("dominant_infected_off_frame_review",
        "infected_off_frame_hotspot_review_run_map.csv"),
      p("dominant_infected_off_frame_review",
        "infected_off_frame_hotspot_review_summary.csv"),
      p("dominant_infected_off_frame_review",
        "infected_off_frame_hotspot_review_metrics.csv"),
      p("dominant_infected_off_frame_review", "figures",
        "infected_off_frame_hotspot_review_priority.png"),
      p("dominant_infected_off_frame_review", "figures",
        "infected_off_frame_hotspot_review_windows.png")
    ),
    infected_off_frame_frame_diagnostics = c(
      p("dominant_infected_off_frame_frame_diagnostics",
        "infected_off_frame_hotspot_frame_diagnostics_runs.csv"),
      p("dominant_infected_off_frame_frame_diagnostics",
        "infected_off_frame_hotspot_frame_diagnostics_hotspots.csv"),
      p("dominant_infected_off_frame_frame_diagnostics",
        "infected_off_frame_hotspot_frame_diagnostics_class_summary.csv"),
      p("dominant_infected_off_frame_frame_diagnostics",
        "infected_off_frame_hotspot_frame_diagnostics_summary_metrics.csv"),
      p("dominant_infected_off_frame_frame_diagnostics",
        "infected_off_frame_b2m_codon55_window_profiles.csv"),
      p("dominant_infected_off_frame_frame_diagnostics", "figures",
        "infected_off_frame_hotspot_frame_diagnostic_quadrants.png"),
      p("dominant_infected_off_frame_frame_diagnostics", "figures",
        "infected_off_frame_b2m_codon55_diagnostic_windows.png")
    ),
    infected_off_frame_readlength_diagnostics = c(
      p("dominant_infected_off_frame_readlength_diagnostics",
        "b2m_readlength_frame_diagnostics_summary_metrics.csv"),
      p("dominant_infected_off_frame_readlength_diagnostics",
        "b2m_readlength_run_summary.csv"),
      p("dominant_infected_off_frame_readlength_diagnostics",
        "b2m_readlength_psite_profile.csv"),
      p("dominant_infected_off_frame_readlength_diagnostics",
        "b2m_readlength_window_frame_summary.csv"),
      p("dominant_infected_off_frame_readlength_diagnostics",
        "b2m_codon55_readlength_summary.csv"),
      p("dominant_infected_off_frame_readlength_diagnostics",
        "b2m_codon55_psite_offset_sensitivity.csv"),
      p("dominant_infected_off_frame_readlength_diagnostics",
        "b2m_codon55_fst_vs_bam_pshift_counts.csv"),
      p("dominant_infected_off_frame_readlength_diagnostics", "figures",
        "b2m_codon55_readlength_heatmap.png"),
      p("dominant_infected_off_frame_readlength_diagnostics", "figures",
        "b2m_codon55_psite_offset_sensitivity.png")
    )
  )
}

pipeline_read_csv_paths <- function(path, columns) {
  if (!file.exists(path) || !requireNamespace("data.table", quietly = TRUE)) {
    return(character())
  }
  dt <- tryCatch(data.table::fread(path, showProgress = FALSE),
                 error = function(e) data.table::data.table())
  paths <- character()
  for (column in intersect(columns, names(dt))) {
    paths <- c(paths, dt[[column]])
  }
  paths[!is.na(paths) & nzchar(paths)]
}

pipeline_canonical_isoform_file <- function() {
  df <- tryCatch(
    ORFik::read.experiment("all_samples-Homo_sapiens", validate = FALSE),
    error = function(e) NULL
  )
  if (is.null(df)) return(character())
  file.path(ORFik::refFolder(df), "canonical_isoforms.txt")
}

resolve_metadata_file <- function() {
  candidates <- c(
    Sys.getenv("DOMINANT_METADATA_FILE", unset = NA_character_),
    "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv",
    path.expand("~/livemount/Bio_data/NGS_pipeline/metadata_done_samples_extended_qc.csv")
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) existing[[1]] else candidates[[1]]
}

pipeline_step_dependencies <- function(step_name, analysis_dir) {
  p <- function(...) file.path(analysis_dir, ...)
  rp <- function(...) file.path(analysis_results_dir(analysis_dir), ...)
  metadata_file <- resolve_metadata_file()
  outputs <- pipeline_step_outputs(analysis_dir)
  canonical_isoform_file <- pipeline_canonical_isoform_file()
  base <- switch(
    step_name,
    manual_translons = c(
      p("manual_translons", "translons_to_be_added_to_manual.txt")
    ),
    fst_pages = {
      required <- p("required_fst_pages.csv")
      c(
        p("scripts", "dominant_cell_state_required_fst_pages.R"),
        p("scripts", "dominant_state_module_definitions.R"),
        rp("curated_inputs", "expanded_oxphos_marker_genes.csv"),
        rp("curated_inputs", "postviral_fatigue_dominant_state_genes.csv"),
        rp("curated_inputs", "curated_additional_gene_priorities.csv"),
        canonical_isoform_file,
        outputs$manual_translons,
        required,
        pipeline_read_csv_paths(
          required,
          c("expected_local_path_after_download", "fst_file_from_index")
        )
      )
    },
    clean_cds = c(
      p("scripts", "dominant_state_module_definitions.R"),
      rp("curated_inputs", "expanded_oxphos_marker_genes.csv"),
      rp("curated_inputs", "postviral_fatigue_dominant_state_genes.csv"),
      rp("curated_inputs", "curated_additional_gene_priorities.csv"),
      outputs$fst_pages,
      canonical_isoform_file,
      outputs$manual_translons
    ),
    score_audit = c(
      outputs$clean_cds,
      p("scripts", "dominant_state_module_definitions.R")
    ),
    state_activity_metadata = c(
      outputs$clean_cds,
      outputs$score_audit,
      p("scripts", "dominant_state_module_definitions.R")
    ),
    auxiliary_state_screen = c(
      metadata_file,
      outputs$clean_cds,
      outputs$score_audit
    ),
    adjusted_enrichment = c(
      metadata_file,
      outputs$clean_cds
    ),
    control_aware = c(
      metadata_file,
      outputs$clean_cds
    ),
    broad_metadata = c(
      metadata_file,
      outputs$clean_cds
    ),
    glmnet = c(
      metadata_file,
      outputs$clean_cds
    ),
    latent = c(
      metadata_file,
      outputs$clean_cds,
      outputs$broad_metadata,
      outputs$glmnet
    ),
    known_uorf = character(),
    uorf_regulation = c(
      metadata_file,
      outputs$clean_cds,
      outputs$manual_translons
    ),
    score_semantics = c(
      metadata_file,
      outputs$clean_cds,
      outputs$score_audit,
      outputs$broad_metadata,
      outputs$uorf_regulation
    ),
    score_semantics_design = c(
      metadata_file,
      outputs$clean_cds,
      outputs$broad_metadata,
      outputs$uorf_regulation,
      outputs$score_semantics
    ),
    rdg_graphs = c(
      pipeline_runtime_option_file(analysis_dir, "rdg_plot_mode"),
      metadata_file,
      outputs$clean_cds,
      outputs$adjusted_enrichment,
      outputs$broad_metadata,
      outputs$glmnet,
      outputs$uorf_regulation
    ),
    uorf_structure = c(
      outputs$uorf_regulation,
      outputs$rdg_graphs
    ),
    next_model = c(
      metadata_file,
      outputs$clean_cds,
      outputs$control_aware,
      outputs$broad_metadata,
      outputs$glmnet,
      outputs$known_uorf,
      outputs$uorf_regulation,
      outputs$rdg_graphs,
      outputs$uorf_structure
    ),
    ouorf_cds_coupling = c(
      outputs$uorf_regulation,
      outputs$rdg_graphs,
      outputs$uorf_structure,
      outputs$next_model
    ),
    branch_usage = c(
      outputs$next_model
    ),
    branch_coverage_qc = c(
      metadata_file,
      outputs$clean_cds,
      outputs$uorf_regulation
    ),
    grouped_branch_usage = c(
      metadata_file,
      outputs$uorf_regulation,
      outputs$branch_coverage_qc
    ),
    branch_allocation = c(
      outputs$grouped_branch_usage
    ),
    downstream_tis_rescue = c(
      outputs$ouorf_cds_coupling
    ),
    fatigue_summary = c(
      metadata_file,
      outputs$clean_cds,
      outputs$adjusted_enrichment,
      outputs$broad_metadata,
      outputs$glmnet,
      outputs$rdg_graphs,
      outputs$next_model,
      outputs$ouorf_cds_coupling,
      outputs$downstream_tis_rescue
    ),
    hierarchical_model = c(
      outputs$next_model,
      outputs$branch_usage,
      outputs$branch_allocation,
      outputs$fatigue_summary
    ),
    hierarchical_residual_audit = c(
      outputs$next_model,
      outputs$hierarchical_model
    ),
    hierarchical_context_calibration = c(
      outputs$hierarchical_model,
      outputs$hierarchical_residual_audit
    ),
    hierarchical_context_model = c(
      outputs$hierarchical_model,
      outputs$hierarchical_residual_audit,
      outputs$hierarchical_context_calibration
    ),
    hierarchical_branch_covariates = c(
      outputs$hierarchical_residual_audit,
      outputs$branch_allocation
    ),
    graph_geometry_model = c(
      outputs$hierarchical_residual_audit,
      outputs$rdg_graphs,
      outputs$uorf_structure
    ),
    mechanistic_motif_model = c(
      outputs$hierarchical_residual_audit,
      outputs$graph_geometry_model,
      outputs$downstream_tis_rescue
    ),
    motif_prior_calibration = c(
      outputs$mechanistic_motif_model,
      outputs$branch_allocation
    ),
    joint_branch_allocation = c(
      outputs$grouped_branch_usage,
      outputs$branch_allocation,
      outputs$mechanistic_motif_model,
      outputs$motif_prior_calibration
    ),
    hierarchical_joint_allocation = c(
      outputs$joint_branch_allocation
    ),
    dirichlet_multinomial_allocation = c(
      outputs$joint_branch_allocation,
      outputs$hierarchical_joint_allocation
    ),
    model_agreement = c(
      outputs$hierarchical_model,
      outputs$motif_prior_calibration,
      outputs$joint_branch_allocation,
      outputs$hierarchical_joint_allocation,
      outputs$dirichlet_multinomial_allocation,
      p("dominant_rdg_atlas", "browser_validation_notes.csv"),
      rp("curated_inputs", "browser_validation_notes.csv")
    ),
    model_agreement_loo = c(
      outputs$model_agreement,
      outputs$hierarchical_joint_allocation,
      outputs$dirichlet_multinomial_allocation
    ),
    model_agreement_transport = c(
      outputs$model_agreement_loo,
      outputs$model_agreement,
      outputs$hierarchical_joint_allocation,
      outputs$dirichlet_multinomial_allocation,
      outputs$joint_branch_allocation,
      resolve_metadata_file()
    ),
    context_interactions = c(
      outputs$joint_branch_allocation,
      outputs$model_agreement_transport
    ),
    context_partial_pooling = c(
      outputs$context_interactions
    ),
    rdg_atlas_cards = c(
      outputs$next_model,
      outputs$branch_usage,
      outputs$grouped_branch_usage,
      outputs$branch_allocation,
      outputs$joint_branch_allocation,
      outputs$hierarchical_joint_allocation,
      outputs$dirichlet_multinomial_allocation,
      outputs$model_agreement,
      outputs$model_agreement_loo,
      outputs$model_agreement_transport,
      outputs$context_interactions,
      outputs$context_partial_pooling,
      outputs$ouorf_cds_coupling,
      outputs$downstream_tis_rescue,
      outputs$fatigue_summary,
      outputs$hierarchical_model,
      outputs$hierarchical_residual_audit,
      outputs$hierarchical_context_calibration,
      outputs$hierarchical_context_model,
      p("dominant_rdg_atlas", "browser_validation_notes.csv"),
      rp("curated_inputs", "browser_validation_notes.csv")
    ),
    pooling_bias_audit = c(
      outputs$branch_coverage_qc,
      outputs$grouped_branch_usage,
      outputs$context_partial_pooling,
      outputs$rdg_atlas_cards,
      p("dominant_rdg_qualified_inspection", "qualified_replicate_summary.csv"),
      p("dominant_cds_frame_bumpiness_qc", "cds_frame_bumpiness_run_metrics.csv")
    ),
    auxiliary_state_rdg_audit = c(
      metadata_file,
      outputs$auxiliary_state_screen,
      outputs$joint_branch_allocation,
      outputs$rdg_atlas_cards
    ),
    allocation_output_disagreements = c(
      outputs$rdg_atlas_cards,
      outputs$model_agreement,
      outputs$model_agreement_loo
    ),
    therapeutic_hypotheses = c(
      outputs$rdg_atlas_cards,
      outputs$model_agreement,
      outputs$model_agreement_transport,
      outputs$context_interactions
    ),
    therapeutic_perturbation_predictions = c(
      outputs$therapeutic_hypotheses,
      outputs$rdg_atlas_cards
    ),
    rdg_consensus = c(
      outputs$rdg_atlas_cards,
      outputs$model_agreement,
      outputs$model_agreement_loo,
      outputs$model_agreement_transport,
      outputs$context_interactions,
      outputs$context_partial_pooling,
      outputs$auxiliary_state_rdg_audit,
      outputs$therapeutic_hypotheses,
      outputs$therapeutic_perturbation_predictions
    ),
    rdg_validation_dossiers = c(
      outputs$rdg_consensus,
      outputs$rdg_atlas_cards,
      outputs$therapeutic_perturbation_predictions
    ),
    rdg_evidence_gaps = c(
      outputs$rdg_validation_dossiers,
      outputs$rdg_atlas_cards,
      outputs$model_agreement,
      outputs$model_agreement_loo,
      outputs$model_agreement_transport,
      outputs$context_interactions,
      outputs$context_partial_pooling,
      outputs$auxiliary_state_rdg_audit,
      outputs$therapeutic_hypotheses,
      outputs$therapeutic_perturbation_predictions,
      outputs$next_model
    ),
    merged_replicate_translon_shift = c(
      outputs$joint_branch_allocation,
      outputs$rdg_atlas_cards,
      outputs$rdg_evidence_gaps,
      outputs$allocation_output_disagreements
    ),
    rdg_review_batches = c(
      outputs$rdg_evidence_gaps,
      outputs$merged_replicate_translon_shift,
      outputs$rdg_validation_dossiers
    ),
    rdg_review_links = c(
      outputs$rdg_review_batches,
      resolve_metadata_file()
    ),
    rdg_review_feedback = c(
      outputs$rdg_review_links,
      p("dominant_rdg_review_links",
        "dominant_rdg_review_batch_template_with_links.csv")
    ),
    rdg_validation_labels = c(
      outputs$rdg_review_feedback,
      outputs$motif_prior_calibration,
      p("dominant_rdg_atlas", "browser_validation_notes.csv"),
      rp("curated_inputs", "browser_validation_notes.csv"),
      p("dominant_rdg_review_links",
        "dominant_rdg_review_batch_template_with_links.csv")
    ),
    cds_frame_bumpiness_qc = c(
      outputs$clean_cds,
      outputs$next_model
    ),
    infected_off_frame_scan = c(
      outputs$cds_frame_bumpiness_qc
    ),
    infected_off_frame_frame_risk = c(
      outputs$infected_off_frame_scan,
      outputs$clean_cds
    ),
    infected_off_frame_hotspots = c(
      outputs$infected_off_frame_scan,
      outputs$clean_cds
    ),
    infected_off_frame_review = c(
      outputs$infected_off_frame_hotspots,
      outputs$cds_frame_bumpiness_qc
    ),
    infected_off_frame_frame_diagnostics = c(
      outputs$infected_off_frame_review,
      outputs$infected_off_frame_frame_risk,
      outputs$clean_cds
    ),
    infected_off_frame_readlength_diagnostics = c(
      outputs$infected_off_frame_frame_diagnostics,
      p("readlength_inputs", "SRR2052924.bam"),
      p("readlength_inputs", "SRR27450688.bam"),
      p("readlength_inputs", "PRJNA285961_shifting_tables.rds"),
      p("readlength_inputs", "PRJNA1062064_shifting_tables.rds")
    ),
    character()
  )
  unique(base[!is.na(base) & nzchar(base)])
}

pipeline_scopes <- function() {
  list(
    validate_manual = c("manual_translons"),
    core = c(
      "manual_translons",
      "fst_pages",
      "clean_cds",
      "score_audit",
      "state_activity_metadata"
    ),
    metadata = c(
      "state_activity_metadata",
      "auxiliary_state_screen",
      "adjusted_enrichment",
      "control_aware",
      "broad_metadata",
      "glmnet",
      "latent"
    ),
    rdg = c(
      "uorf_regulation",
      "score_semantics",
      "score_semantics_design",
      "rdg_graphs",
      "uorf_structure",
      "next_model",
      "branch_usage",
      "branch_coverage_qc",
      "grouped_branch_usage",
      "branch_allocation",
      "ouorf_cds_coupling",
      "downstream_tis_rescue",
      "fatigue_summary",
      "hierarchical_model",
      "hierarchical_residual_audit",
      "hierarchical_context_calibration",
      "hierarchical_context_model",
      "hierarchical_branch_covariates",
      "graph_geometry_model",
      "mechanistic_motif_model",
      "motif_prior_calibration",
      "joint_branch_allocation",
      "hierarchical_joint_allocation",
      "dirichlet_multinomial_allocation",
      "model_agreement",
      "model_agreement_loo",
      "model_agreement_transport",
      "context_interactions",
      "context_partial_pooling",
      "rdg_atlas_cards",
      "pooling_bias_audit",
      "auxiliary_state_rdg_audit",
      "allocation_output_disagreements",
      "therapeutic_hypotheses",
      "therapeutic_perturbation_predictions",
      "rdg_consensus",
      "rdg_validation_dossiers",
      "rdg_evidence_gaps",
      "merged_replicate_translon_shift",
      "rdg_review_batches",
      "rdg_review_links",
      "rdg_review_feedback",
      "rdg_validation_labels"
    ),
    frame_qc = c(
      "cds_frame_bumpiness_qc",
      "infected_off_frame_scan",
      "infected_off_frame_frame_risk",
      "infected_off_frame_hotspots",
      "infected_off_frame_review",
      "infected_off_frame_frame_diagnostics",
      "infected_off_frame_readlength_diagnostics"
    ),
    manual = c(
      "manual_translons",
      "fst_pages",
      "clean_cds",
      "score_audit",
      "state_activity_metadata",
      "auxiliary_state_screen",
      "adjusted_enrichment",
      "control_aware",
      "broad_metadata",
      "glmnet",
      "latent",
      "uorf_regulation",
      "score_semantics",
      "score_semantics_design",
      "rdg_graphs",
      "uorf_structure",
      "next_model",
      "branch_usage",
      "branch_coverage_qc",
      "grouped_branch_usage",
      "branch_allocation",
      "ouorf_cds_coupling",
      "downstream_tis_rescue",
      "fatigue_summary",
      "hierarchical_model",
      "hierarchical_residual_audit",
      "hierarchical_context_calibration",
      "hierarchical_context_model",
      "hierarchical_branch_covariates",
      "graph_geometry_model",
      "mechanistic_motif_model",
      "motif_prior_calibration",
      "joint_branch_allocation",
      "hierarchical_joint_allocation",
      "dirichlet_multinomial_allocation",
      "model_agreement",
      "model_agreement_loo",
      "model_agreement_transport",
      "context_interactions",
      "context_partial_pooling",
      "rdg_atlas_cards",
      "pooling_bias_audit",
      "auxiliary_state_rdg_audit",
      "allocation_output_disagreements",
      "therapeutic_hypotheses",
      "therapeutic_perturbation_predictions",
      "rdg_consensus",
      "rdg_validation_dossiers",
      "rdg_evidence_gaps",
      "merged_replicate_translon_shift",
      "rdg_review_batches",
      "rdg_review_links",
      "rdg_review_feedback",
      "rdg_validation_labels"
    ),
    all = c(
      "manual_translons",
      "fst_pages",
      "clean_cds",
      "score_audit",
      "state_activity_metadata",
      "auxiliary_state_screen",
      "adjusted_enrichment",
      "control_aware",
      "broad_metadata",
      "glmnet",
      "latent",
      "known_uorf",
      "uorf_regulation",
      "score_semantics",
      "score_semantics_design",
      "rdg_graphs",
      "uorf_structure",
      "next_model",
      "branch_usage",
      "branch_coverage_qc",
      "grouped_branch_usage",
      "branch_allocation",
      "ouorf_cds_coupling",
      "downstream_tis_rescue",
      "fatigue_summary",
      "hierarchical_model",
      "hierarchical_residual_audit",
      "hierarchical_context_calibration",
      "hierarchical_context_model",
      "hierarchical_branch_covariates",
      "graph_geometry_model",
      "mechanistic_motif_model",
      "motif_prior_calibration",
      "joint_branch_allocation",
      "hierarchical_joint_allocation",
      "dirichlet_multinomial_allocation",
      "model_agreement",
      "model_agreement_loo",
      "model_agreement_transport",
      "context_interactions",
      "context_partial_pooling",
      "rdg_atlas_cards",
      "pooling_bias_audit",
      "auxiliary_state_rdg_audit",
      "allocation_output_disagreements",
      "therapeutic_hypotheses",
      "therapeutic_perturbation_predictions",
      "rdg_consensus",
      "rdg_validation_dossiers",
      "rdg_evidence_gaps",
      "merged_replicate_translon_shift",
      "rdg_review_batches",
      "rdg_review_links",
      "rdg_review_feedback",
      "rdg_validation_labels"
    )
  )
}

slice_steps <- function(selected, from = NULL, to = NULL) {
  if (!is.null(from)) {
    idx <- match(from, selected)
    if (is.na(idx)) stop("--from step is not in selected scope: ", from, call. = FALSE)
    selected <- selected[idx:length(selected)]
  }
  if (!is.null(to)) {
    idx <- match(to, selected)
    if (is.na(idx)) stop("--to step is not in selected scope: ", to, call. = FALSE)
    selected <- selected[seq_len(idx)]
  }
  selected
}

run_one_step <- function(step_name, script_file, analysis_dir,
                         rdg_plot_mode = "all") {
  if (!file.exists(script_file)) {
    stop("Missing pipeline script for step '", step_name, "': ", script_file,
         call. = FALSE)
  }

  old_wd <- setwd(analysis_dir)
  on.exit(setwd(old_wd), add = TRUE)
  old_rdg_plot_mode <- Sys.getenv("DOMINANT_RDG_PLOT_MODE", unset = NA_character_)
  old_analysis_dir <- Sys.getenv("DOMINANT_CELL_STATES_DIR", unset = NA_character_)
  old_results_dir <- Sys.getenv("DOMINANT_CELL_STATES_RESULTS_DIR",
                                unset = NA_character_)
  Sys.setenv(DOMINANT_RDG_PLOT_MODE = rdg_plot_mode)
  Sys.setenv(
    DOMINANT_CELL_STATES_DIR = analysis_dir,
    DOMINANT_CELL_STATES_RESULTS_DIR = analysis_results_dir(
      analysis_dir,
      create = TRUE
    )
  )
  on.exit({
    if (is.na(old_rdg_plot_mode)) {
      Sys.unsetenv("DOMINANT_RDG_PLOT_MODE")
    } else {
      Sys.setenv(DOMINANT_RDG_PLOT_MODE = old_rdg_plot_mode)
    }
    if (is.na(old_analysis_dir)) {
      Sys.unsetenv("DOMINANT_CELL_STATES_DIR")
    } else {
      Sys.setenv(DOMINANT_CELL_STATES_DIR = old_analysis_dir)
    }
    if (is.na(old_results_dir)) {
      Sys.unsetenv("DOMINANT_CELL_STATES_RESULTS_DIR")
    } else {
      Sys.setenv(DOMINANT_CELL_STATES_RESULTS_DIR = old_results_dir)
    }
  }, add = TRUE)

  step_env <- new.env(parent = globalenv())
  step_env$dominant_rdg_plot_mode <- rdg_plot_mode
  step_results_dir <- analysis_results_dir(analysis_dir, create = TRUE)
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
  is_analysis_root_arg <- function(x) {
    if (!is.character(x) || length(x) != 1L) return(FALSE)
    normalized <- normalizePath(x, mustWork = FALSE)
    x %in% c(".", analysis_dir) || identical(normalized, analysis_dir)
  }
  is_generated_root_entry <- function(x) {
    if (!is.character(x) || length(x) < 1L) return(FALSE)
    any(vapply(generated_root_names, grepl, logical(1), x = x[[1]]))
  }
  is_generated_manual_entry <- function(parts) {
    length(parts) >= 3L &&
      identical(parts[[2]], "manual_translons") &&
      is.character(parts[[3]]) &&
      parts[[3]][[1]] %in% generated_manual_files
  }
  step_env$file.path <- function(..., fsep = .Platform$file.sep) {
    parts <- list(...)
    if (length(parts) >= 2L &&
        is_analysis_root_arg(parts[[1]]) &&
        (is_generated_root_entry(parts[[2]]) ||
         is_generated_manual_entry(parts))) {
      parts[[1]] <- step_results_dir
    }
    do.call(base::file.path, c(parts, list(fsep = fsep)))
  }
  htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
  if (file.exists(htmlwidget_helper)) {
    source(htmlwidget_helper, local = step_env)
  }
  source(script_file, local = step_env)
}

run_dominant_cell_state_pipeline <- function(scope = "all",
                                             steps = NULL,
                                             dry_run = FALSE,
                                             seed_cache = FALSE,
                                             stop_on_error = TRUE,
                                             from = NULL,
                                             to = NULL,
                                             use_cache = TRUE,
                                             force = FALSE,
                                             force_steps = character(),
                                             rdg_plot_mode = Sys.getenv(
                                               "DOMINANT_RDG_PLOT_MODE",
                                               "all"
                                             ),
                                             cache_version = Sys.getenv(
                                               "DOMINANT_PIPELINE_CACHE_VERSION",
                                               "dominant-pipeline-cache-v1"
                                             ),
                                             repo_root = find_repo_root(),
                                             analysis_dir = find_analysis_dir(
                                               repo_root = repo_root
                                             )) {
  repo_root <- normalizePath(repo_root, mustWork = TRUE)
  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)
  if (!dir.exists(analysis_dir)) {
    stop("Missing analysis directory: ", analysis_dir, call. = FALSE)
  }

  setup_runtime_env(analysis_dir)
  results_dir <- analysis_results_dir(analysis_dir, create = TRUE)
  source_pipeline_cache(analysis_dir)
  load_ribocrypt(repo_root)
  Sys.setenv(RIBOCRYPT_REPO = repo_root)
  rdg_plot_mode <- normalize_rdg_plot_mode(rdg_plot_mode)
  write_pipeline_runtime_options(analysis_dir, rdg_plot_mode)

  all_steps <- pipeline_steps(analysis_dir)
  scopes <- pipeline_scopes()

  if (is.null(steps)) {
    if (!scope %in% names(scopes)) {
      stop("Unknown scope '", scope, "'. Available scopes: ",
           paste(names(scopes), collapse = ", "), call. = FALSE)
    }
    selected <- scopes[[scope]]
  } else {
    selected <- steps
  }

  unknown <- setdiff(selected, all_steps$step)
  if (length(unknown)) {
    stop("Unknown step(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  selected <- slice_steps(selected, from = from, to = to)
  step_table <- all_steps[match(selected, all_steps$step), , drop = FALSE]

  message("Dominant-cell-state pipeline")
  message("  repo:  ", repo_root)
  message("  analysis: ", analysis_dir)
  message("  scope: ", if (is.null(steps)) scope else "<custom steps>")
  message("  steps: ", paste(step_table$step, collapse = " -> "))
  message("  cache: ", if (use_cache) {
    paste0("on (version ", cache_version, ")")
  } else {
    "off"
  })
  message("  RDG plot mode: ", rdg_plot_mode)

  if (dry_run) {
    message("Dry run only; no scripts were executed.")
    return(invisible(step_table))
  }

  summary <- data.frame(
    step = character(),
    script = character(),
    status = character(),
    started = character(),
    ended = character(),
    seconds = numeric(),
    error = character(),
    cache_reason = character(),
    stringsAsFactors = FALSE
  )

  output_map <- pipeline_step_outputs(analysis_dir)
  source_helpers <- file.path(
    analysis_dir,
    "scripts",
    c(
      "dominant_state_module_definitions.R",
      "dominant_coverage_cache.R",
      "dominant_internal_orf_signal_helpers.R",
      "dominant_htmlwidgets.R"
    )
  )

  if (seed_cache) {
    for (i in seq_len(nrow(step_table))) {
      step_name <- step_table$step[[i]]
      script_file <- step_table$script[[i]]
      outputs <- output_map[[step_name]]
      dependencies <- pipeline_step_dependencies(step_name, analysis_dir)
      started <- Sys.time()
      missing_outputs <- pipeline_missing_outputs(outputs)
      ended <- Sys.time()

      if (length(missing_outputs)) {
        message("[", timestamp(), "] cannot seed ", step_name,
                " (missing output: ",
                paste(basename(missing_outputs), collapse = ", "), ")")
        summary <- rbind(
          summary,
          data.frame(
            step = step_name,
            script = basename(script_file),
            status = "seed_skipped_missing_outputs",
            started = format(started, "%Y-%m-%d %H:%M:%S"),
            ended = format(ended, "%Y-%m-%d %H:%M:%S"),
            seconds = as.numeric(difftime(ended, started, units = "secs")),
            error = "",
            cache_reason = paste0("missing output(s): ",
                                  paste(basename(missing_outputs),
                                        collapse = ", ")),
            stringsAsFactors = FALSE
          )
        )
        next
      }

      cache_file <- pipeline_cache_write(
        step_name = step_name,
        script_file = script_file,
        analysis_dir = analysis_dir,
        outputs = outputs,
        dependencies = dependencies,
        cache_version = cache_version,
        extra_source_files = source_helpers,
        seconds = NA_real_
      )
      message("[", timestamp(), "] seeded cache for ", step_name)
      summary <- rbind(
        summary,
        data.frame(
          step = step_name,
          script = basename(script_file),
          status = "seeded_cache",
          started = format(started, "%Y-%m-%d %H:%M:%S"),
          ended = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
          error = "",
          cache_reason = paste0("cache manifest written: ", basename(cache_file)),
          stringsAsFactors = FALSE
        )
      )
    }

    out_file <- file.path(results_dir, "dominant_pipeline_last_run_summary.csv")
    data.table::fwrite(summary, out_file)
    message("Seeded cache only; no scripts were executed.")
    message("Saved run summary: ", out_file)
    return(invisible(summary))
  }

  for (i in seq_len(nrow(step_table))) {
    step_name <- step_table$step[[i]]
    script_file <- step_table$script[[i]]
    outputs <- output_map[[step_name]]
    dependencies <- pipeline_step_dependencies(step_name, analysis_dir)
    started <- Sys.time()
    forced <- isTRUE(force) || step_name %in% force_steps

    if (use_cache && !forced) {
      cache_state <- pipeline_cache_status(
        step_name = step_name,
        script_file = script_file,
        analysis_dir = analysis_dir,
        outputs = outputs,
        dependencies = dependencies,
        cache_version = cache_version,
        extra_source_files = source_helpers
      )
      if (isTRUE(cache_state$valid)) {
        ended <- Sys.time()
        message("[", timestamp(), "] skipping ", step_name,
                " (", cache_state$reason, ")")
        summary <- rbind(
          summary,
          data.frame(
            step = step_name,
            script = basename(script_file),
            status = "skipped_cache",
            started = format(started, "%Y-%m-%d %H:%M:%S"),
            ended = format(ended, "%Y-%m-%d %H:%M:%S"),
            seconds = as.numeric(difftime(ended, started, units = "secs")),
            error = "",
            cache_reason = cache_state$reason,
            stringsAsFactors = FALSE
          )
        )
        next
      }
      message("[", timestamp(), "] starting ", step_name,
              " (cache miss: ", cache_state$reason, ")")
    } else {
      message("[", timestamp(), "] starting ", step_name,
              if (forced) " (forced)" else "")
    }

    err <- NULL
    tryCatch(
      run_one_step(step_name, script_file, analysis_dir,
                   rdg_plot_mode = rdg_plot_mode),
      error = function(e) {
        err <<- conditionMessage(e)
      }
    )

    ended <- Sys.time()
    ok <- is.null(err)
    elapsed <- as.numeric(difftime(ended, started, units = "secs"))
    cache_reason <- ""
    if (ok && use_cache && !is.null(outputs)) {
      cache_file <- pipeline_cache_write(
        step_name = step_name,
        script_file = script_file,
        analysis_dir = analysis_dir,
        outputs = outputs,
        dependencies = dependencies,
        cache_version = cache_version,
        extra_source_files = source_helpers,
        seconds = elapsed
      )
      cache_reason <- paste0("cache manifest written: ", basename(cache_file))
    }
    summary <- rbind(
      summary,
      data.frame(
        step = step_name,
        script = basename(script_file),
        status = if (ok) "ok" else "error",
        started = format(started, "%Y-%m-%d %H:%M:%S"),
        ended = format(ended, "%Y-%m-%d %H:%M:%S"),
        seconds = elapsed,
        error = if (ok) "" else err,
        cache_reason = if (ok) cache_reason else "",
        stringsAsFactors = FALSE
      )
    )

    if (ok) {
      message("[", timestamp(), "] finished ", step_name)
    } else {
      message("[", timestamp(), "] failed ", step_name, ": ", err)
      if (stop_on_error) {
        break
      }
    }
  }

  out_file <- file.path(results_dir, "dominant_pipeline_last_run_summary.csv")
  data.table::fwrite(summary, out_file)
  message("Saved run summary: ", out_file)

  if (any(summary$status == "error") && stop_on_error) {
    stop("Pipeline stopped after an error. See: ", out_file, call. = FALSE)
  }

  invisible(summary)
}

print_pipeline_help <- function() {
  scopes <- pipeline_scopes()
  cat(
    "Dominant-cell-state pipeline\n\n",
    "Usage:\n",
    "  Rscript dominant_cell_states/scripts/run_dominant_cell_state_pipeline.R --scope all\n",
    "  Rscript dominant_cell_states/scripts/run_dominant_cell_state_pipeline.R --scope manual\n",
    "  Rscript dominant_cell_states/scripts/run_dominant_cell_state_pipeline.R --scope validate_manual\n\n",
    "Scopes:\n",
    paste(sprintf("  %-15s %s", names(scopes), vapply(scopes, paste, character(1), collapse = " -> ")),
          collapse = "\n"),
    "\n\nOptions:\n",
    "  --scope NAME       One of the scopes above. Default: all\n",
    "  --from STEP        Start at this step inside the selected scope.\n",
    "  --to STEP          Stop after this step inside the selected scope.\n",
    "  --dry-run          Print selected steps without running them.\n",
    "  --seed-cache       Write cache manifests for existing outputs without running scripts.\n",
    "  --no-cache         Run selected steps without reading or writing cache manifests.\n",
    "  --force            Run selected steps even when cache manifests are valid.\n",
    "  --force-step STEP  Force one named step; repeat or use commas for several steps.\n",
    "  --skip-rdg-plots   Run RDG graph step in table-only mode.\n",
    "  --rdg-plot-mode X  RDG graph plotting mode: all or tables.\n",
    "  --cache-version X  Add a manual invalidation string to every cache key.\n",
    "  --keep-going       Continue after step errors and write them to the summary CSV.\n",
    "  --list             Print available steps and exit.\n",
    "  --help             Print this help and exit.\n",
    sep = ""
  )
}

parse_pipeline_args <- function(args) {
  cfg <- list(
    scope = "all",
    dry_run = FALSE,
    seed_cache = FALSE,
    stop_on_error = TRUE,
    from = NULL,
    to = NULL,
    use_cache = TRUE,
    force = FALSE,
    force_steps = character(),
    rdg_plot_mode = Sys.getenv("DOMINANT_RDG_PLOT_MODE", "all"),
    cache_version = Sys.getenv(
      "DOMINANT_PIPELINE_CACHE_VERSION",
      "dominant-pipeline-cache-v1"
    ),
    list = FALSE,
    help = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--help", "-h")) {
      cfg$help <- TRUE
    } else if (arg == "--list") {
      cfg$list <- TRUE
    } else if (arg == "--dry-run") {
      cfg$dry_run <- TRUE
    } else if (arg == "--seed-cache") {
      cfg$seed_cache <- TRUE
    } else if (arg == "--no-cache") {
      cfg$use_cache <- FALSE
    } else if (arg == "--force") {
      cfg$force <- TRUE
    } else if (arg == "--force-step") {
      i <- i + 1L
      if (i > length(args)) stop("--force-step requires a value", call. = FALSE)
      cfg$force_steps <- c(
        cfg$force_steps,
        unlist(strsplit(args[[i]], ",", fixed = TRUE), use.names = FALSE)
      )
    } else if (startsWith(arg, "--force-step=")) {
      value <- sub("^--force-step=", "", arg)
      cfg$force_steps <- c(
        cfg$force_steps,
        unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE)
      )
    } else if (arg == "--skip-rdg-plots") {
      cfg$rdg_plot_mode <- "tables"
    } else if (arg == "--rdg-plot-mode") {
      i <- i + 1L
      if (i > length(args)) stop("--rdg-plot-mode requires a value", call. = FALSE)
      cfg$rdg_plot_mode <- args[[i]]
    } else if (startsWith(arg, "--rdg-plot-mode=")) {
      cfg$rdg_plot_mode <- sub("^--rdg-plot-mode=", "", arg)
    } else if (arg == "--cache-version") {
      i <- i + 1L
      if (i > length(args)) stop("--cache-version requires a value", call. = FALSE)
      cfg$cache_version <- args[[i]]
    } else if (startsWith(arg, "--cache-version=")) {
      cfg$cache_version <- sub("^--cache-version=", "", arg)
    } else if (arg == "--keep-going") {
      cfg$stop_on_error <- FALSE
    } else if (arg == "--scope") {
      i <- i + 1L
      if (i > length(args)) stop("--scope requires a value", call. = FALSE)
      cfg$scope <- args[[i]]
    } else if (startsWith(arg, "--scope=")) {
      cfg$scope <- sub("^--scope=", "", arg)
    } else if (arg == "--from") {
      i <- i + 1L
      if (i > length(args)) stop("--from requires a value", call. = FALSE)
      cfg$from <- args[[i]]
    } else if (startsWith(arg, "--from=")) {
      cfg$from <- sub("^--from=", "", arg)
    } else if (arg == "--to") {
      i <- i + 1L
      if (i > length(args)) stop("--to requires a value", call. = FALSE)
      cfg$to <- args[[i]]
    } else if (startsWith(arg, "--to=")) {
      cfg$to <- sub("^--to=", "", arg)
    } else if (!startsWith(arg, "-")) {
      cfg$scope <- arg
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
    i <- i + 1L
  }

  cfg$force_steps <- unique(trimws(cfg$force_steps))
  cfg$force_steps <- cfg$force_steps[nzchar(cfg$force_steps)]
  cfg$rdg_plot_mode <- normalize_rdg_plot_mode(cfg$rdg_plot_mode)
  cfg
}

run_dominant_cell_state_pipeline_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- parse_pipeline_args(args)

  if (cfg$help) {
    print_pipeline_help()
    return(invisible(NULL))
  }

  repo_root <- find_repo_root()
  analysis_dir <- find_analysis_dir(repo_root = repo_root)

  if (cfg$list) {
    print(pipeline_steps(analysis_dir), row.names = FALSE)
    return(invisible(NULL))
  }

  run_dominant_cell_state_pipeline(
    scope = cfg$scope,
    dry_run = cfg$dry_run,
    seed_cache = cfg$seed_cache,
    stop_on_error = cfg$stop_on_error,
    from = cfg$from,
    to = cfg$to,
    use_cache = cfg$use_cache,
    force = cfg$force,
    force_steps = cfg$force_steps,
    rdg_plot_mode = cfg$rdg_plot_mode,
    cache_version = cfg$cache_version,
    repo_root = repo_root,
    analysis_dir = analysis_dir
  )
}

if (any(grepl("^--file=", commandArgs(trailingOnly = FALSE)))) {
  run_dominant_cell_state_pipeline_cli()
}
