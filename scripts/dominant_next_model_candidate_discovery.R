dominant_loader <- c(
  file.path("scripts", "dominant_ribocrypt_loader.R"),
  file.path("dominant_cell_states", "scripts", "dominant_ribocrypt_loader.R")
)
dominant_loader <- dominant_loader[file.exists(dominant_loader)][1]
if (!is.na(dominant_loader)) {
  source(dominant_loader)
  dominant_load_ribocrypt_if_available()
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(plotly)
  library(htmlwidgets)
})

htmlwidget_helper <- c(
  file.path("scripts", "dominant_htmlwidgets.R"),
  "dominant_htmlwidgets.R",
  file.path("dominant_cell_states", "scripts", "dominant_htmlwidgets.R")
)
htmlwidget_helper <- htmlwidget_helper[file.exists(htmlwidget_helper)][1]
if (!is.na(htmlwidget_helper)) source(htmlwidget_helper)

condition_family_helper <- c(
  file.path("scripts", "dominant_condition_families.R"),
  "dominant_condition_families.R",
  file.path("dominant_cell_states", "scripts",
            "dominant_condition_families.R")
)
condition_family_helper <-
  condition_family_helper[file.exists(condition_family_helper)][1]
if (!is.na(condition_family_helper)) source(condition_family_helper)

analysis_dir <- if (file.exists("dominant_rdg_outputs/rdg_feature_annotation.csv")) {
  "."
} else if (file.exists("dominant_cell_states/dominant_rdg_outputs/rdg_feature_annotation.csv")) {
  "dominant_cell_states"
} else {
  "."
}

results_dir <- Sys.getenv(
  "DOMINANT_CELL_STATES_RESULTS_DIR",
  unset = file.path(analysis_dir, "results")
)
output_dir <- file.path(analysis_dir, "dominant_next_model")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

message("Dominant-state next-model candidate discovery:")
message("  1. Build constrained RDG flow proxy tables")
message("  2. Compact branch-feature contrasts into gene/state/branch/edge rows")
message("  3. Combine structure, allocation, metadata confidence, and known-uORF status")
message("  4. Add compact latent reductions, matched-control uncertainty, and rank candidate discoveries")

feature_annotation_file <- file.path(analysis_dir, "dominant_rdg_outputs", "rdg_feature_annotation.csv")
relative_usage_file <- file.path(analysis_dir, "dominant_rdg_outputs", "relative_usage_matrices", "relative_usage_in_long.csv")
branch_contrast_file <- file.path(analysis_dir, "dominant_rdg_outputs", "rdg_branch_feature_contrasts.csv")
structure_file <- file.path(analysis_dir, "dominant_uorf_structure_rules", "uorf_structure_gene_rule_summary.csv")
allocation_file <- file.path(analysis_dir, "dominant_uorf_regulation_allocation_model.csv")
feature_expression_file <- file.path(analysis_dir, "dominant_uorf_feature_expression.csv")
control_file <- file.path(analysis_dir, "dominant_state_control_aware_enrichment_review.csv")
glmnet_selected_file <- file.path(analysis_dir, "dominant_state_glmnet_metadata_selected_term_scores.csv")
glmnet_stability_file <- file.path(analysis_dir, "dominant_state_glmnet_metadata_stability_selection.csv")
rf_performance_file <- file.path(analysis_dir, "dominant_state_broad_metadata_rf_performance.csv")
state_file <- file.path(analysis_dir, "human_dominant_cell_states_clean_cds.csv")
diagnostics_file <- file.path(
  analysis_dir,
  "human_dominant_cell_states_clean_cds_gene_diagnostics.csv"
)
known_uorf_file <- file.path(
  results_dir,
  "curated_inputs",
  "known_uorf_effects_with_citations.csv"
)

flow_output <- file.path(output_dir, "dominant_next_model_rdg_flow_proxy.csv")
compact_edge_output <- file.path(output_dir, "dominant_next_model_compact_gene_state_branch_edge.csv")
candidate_output <- file.path(output_dir, "dominant_next_model_candidate_rankings.csv")
novel_output <- file.path(output_dir, "dominant_next_model_novel_discovery_candidates.csv")
latent_output <- file.path(output_dir, "dominant_next_model_latent_reductions.csv")
summary_output <- file.path(output_dir, "dominant_next_model_summary_metrics.csv")
buffering_all_output <- file.path(output_dir, "dominant_next_model_cds_buffering_effect_sizes_all.csv")
buffering_relevant_output <- file.path(output_dir, "dominant_next_model_cds_buffering_effect_sizes_relevant_rdgs.csv")
buffering_condition_output <- file.path(output_dir, "dominant_next_model_cds_buffering_effect_sizes_condition_focus.csv")
buffering_condition_representative_output <- file.path(output_dir, "dominant_next_model_cds_buffering_effect_sizes_condition_representative.csv")
candidate_v2_output <- file.path(output_dir, "dominant_next_model_candidate_rankings_v2.csv")
candidate_v3_output <- file.path(output_dir, "dominant_next_model_candidate_rankings_v3.csv")
candidate_v4_output <- file.path(output_dir, "dominant_next_model_candidate_rankings_v4.csv")
matched_feature_contrast_output <- file.path(
  output_dir,
  "dominant_next_model_matched_control_feature_contrasts.csv"
)
condition_design_family_output <- file.path(
  output_dir,
  "dominant_next_model_condition_design_family_audit.csv"
)
matched_cds_buffering_output <- file.path(
  output_dir,
  "dominant_next_model_matched_control_cds_buffering.csv"
)
matched_cds_bootstrap_output <- file.path(
  output_dir,
  "dominant_next_model_matched_control_cds_bootstrap.csv"
)
matched_gene_summary_output <- file.path(
  output_dir,
  "dominant_next_model_matched_control_gene_summary.csv"
)
matched_design_summary_output <- file.path(
  output_dir,
  "dominant_next_model_matched_control_design_summary.csv"
)
generic_condition_evidence_output <- file.path(
  output_dir,
  "dominant_next_model_generic_condition_evidence.csv"
)
internal_orf_matched_bump_output <- file.path(
  output_dir,
  "dominant_next_model_internal_orf_matched_bump_evidence.csv"
)
internal_orf_initiation_output <- file.path(
  output_dir,
  "dominant_next_model_internal_orf_initiation_evidence.csv"
)
internal_orf_architecture_output <- file.path(
  output_dir,
  "dominant_next_model_internal_orf_architecture.csv"
)
browser_review_handoff_output <- file.path(
  output_dir,
  "dominant_next_model_browser_review_handoff.csv"
)
feature_confidence_output <- file.path(
  output_dir,
  "dominant_next_model_feature_confidence.csv"
)
metadata_file <- local({
  candidates <- c(
    Sys.getenv("DOMINANT_METADATA_FILE", unset = NA_character_),
    "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv",
    path.expand("~/livemount/Bio_data/NGS_pipeline/metadata_done_samples_extended_qc.csv")
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) existing[[1]] else candidates[[1]]
})
bootstrap_iterations <- 120L
min_matched_group_samples <- 2L
min_matched_feature_raw_counts_any <- 30L
min_matched_feature_raw_counts_both_strict <- 30L

safe_fread <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop("Missing required input: ", path)
    return(data.table())
  }
  fread(path)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(max(x))
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(min(x))
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x)
}

safe_div <- function(num, den, eps = 1e-9) {
  out <- rep(NA_real_, length(num))
  ok <- is.finite(num) & is.finite(den) & abs(den) > eps
  out[ok] <- num[ok] / den[ok]
  out
}

safe_quantile <- function(x, prob, default = NA_real_) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(default)
  as.numeric(quantile(x, prob, na.rm = TRUE, names = FALSE))
}

control_status_score <- function(status) {
  status <- as.character(status)
  fcase(
    status == "supported_by_matched_control", 2,
    status == "supported_by_matched_controls", 2,
    status == "matched_control_weak_or_underpowered", 1,
    status == "possible_global_false_positive", -1,
    status == "no_matched_control_available", 0,
    default = 0
  )
}

extract_level_field_value <- function(field, level, target_field) {
  field <- as.character(field)
  level <- as.character(level)
  if (is.na(field) || is.na(level) || !nzchar(field) || !nzchar(level)) {
    return(NA_character_)
  }
  if (!grepl(":", field, fixed = TRUE)) {
    if (identical(field, target_field)) return(level)
    return(NA_character_)
  }
  parts <- strsplit(level, " | ", fixed = TRUE)[[1]]
  for (part in parts) {
    eq <- regexpr("=", part, fixed = TRUE)[1]
    if (!is.finite(eq) || eq < 1L) next
    key <- substr(part, 1L, eq - 1L)
    value <- substr(part, eq + 1L, nchar(part))
    if (identical(key, target_field)) return(value)
  }
  NA_character_
}

norm01 <- function(x, inverse = FALSE) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  rng <- range(x[ok])
  if (rng[[1]] == rng[[2]]) {
    out[ok] <- 0
  } else {
    out[ok] <- (x[ok] - rng[[1]]) / diff(rng)
  }
  if (inverse) out[ok] <- 1 - out[ok]
  out
}

safe_setcolorder <- function(dt, cols) {
  present <- intersect(cols, names(dt))
  setcolorder(dt, c(present, setdiff(names(dt), present)))
  invisible(dt)
}

bootstrap_cds_contrast <- function(case_rel, control_rel,
                                   case_fpkm, control_fpkm,
                                   n_boot = bootstrap_iterations) {
  case_rel <- as.numeric(case_rel[is.finite(case_rel)])
  control_rel <- as.numeric(control_rel[is.finite(control_rel)])
  case_fpkm <- as.numeric(case_fpkm[is.finite(case_fpkm)])
  control_fpkm <- as.numeric(control_fpkm[is.finite(control_fpkm)])
  if (length(case_rel) < 2 || length(control_rel) < 2 ||
      length(case_fpkm) < 2 || length(control_fpkm) < 2) {
    return(data.table(
      bootstrap_iterations = 0L,
      bootstrap_relative_use_delta_sd = NA_real_,
      bootstrap_relative_use_delta_lower = NA_real_,
      bootstrap_relative_use_delta_upper = NA_real_,
      bootstrap_fpkm_log2_fc_sd = NA_real_,
      bootstrap_fpkm_log2_fc_lower = NA_real_,
      bootstrap_fpkm_log2_fc_upper = NA_real_
    ))
  }
  rel_delta <- numeric(n_boot)
  fpkm_log2_fc <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    case_rel_b <- case_rel[sample.int(length(case_rel), length(case_rel), replace = TRUE)]
    control_rel_b <- control_rel[sample.int(length(control_rel), length(control_rel), replace = TRUE)]
    case_fpkm_b <- case_fpkm[sample.int(length(case_fpkm), length(case_fpkm), replace = TRUE)]
    control_fpkm_b <- control_fpkm[sample.int(length(control_fpkm), length(control_fpkm), replace = TRUE)]
    rel_delta[[b]] <- mean(case_rel_b, na.rm = TRUE) - mean(control_rel_b, na.rm = TRUE)
    fpkm_log2_fc[[b]] <- log2((mean(case_fpkm_b, na.rm = TRUE) + 1) /
                                (mean(control_fpkm_b, na.rm = TRUE) + 1))
  }
  data.table(
    bootstrap_iterations = n_boot,
    bootstrap_relative_use_delta_sd = safe_sd(rel_delta),
    bootstrap_relative_use_delta_lower = safe_quantile(rel_delta, 0.025),
    bootstrap_relative_use_delta_upper = safe_quantile(rel_delta, 0.975),
    bootstrap_fpkm_log2_fc_sd = safe_sd(fpkm_log2_fc),
    bootstrap_fpkm_log2_fc_lower = safe_quantile(fpkm_log2_fc, 0.025),
    bootstrap_fpkm_log2_fc_upper = safe_quantile(fpkm_log2_fc, 0.975)
  )
}

file_safe <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  substr(x, 1L, 160L)
}

save_widget <- function(widget, output_file, title) {
  dominant_save_widget(widget, output_file = output_file, title = title)
}

write_plot <- function(plot, name, width = 9, height = 6) {
  png_file <- file.path(figure_dir, paste0(name, ".png"))
  pdf_file <- file.path(figure_dir, paste0(name, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, dpi = 220, bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, bg = "white")
}

state_short <- function(x) {
  x <- gsub(" / ", " ", x, fixed = TRUE)
  x <- gsub(" \\(.*\\)$", "", x)
  x <- gsub("mTOR Translation capacity.*", "mTOR/TOP", x)
  x <- gsub("Proliferation Cell cycle", "Proliferation", x)
  x <- gsub("Hypoxia HIF1A program", "Hypoxia", x)
  x <- gsub("Interferon antiviral", "Interferon", x)
  x <- gsub("EMT Mesenchymal shift", "EMT", x)
  x <- gsub("OXPHOS mitochondrial", "OXPHOS", x)
  x <- gsub("Post-viral fatigue ribosome stress", "Post-viral fatigue", x)
  x
}

truncate_label <- function(x, max_chars = 48L) {
  x <- as.character(x)
  too_long <- nchar(x) > max_chars
  x[too_long] <- paste0(substr(x[too_long], 1L, max_chars - 3L), "...")
  x
}

clean_level <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "NA", "N/A", "na", "n/a", "NULL", "null", "None", "none")] <- NA_character_
  x
}

normalise_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

is_control_condition <- function(x) {
  y <- normalise_label(x)
  y %chin% c(
    "wt", "wildtype", "wild_type", "control", "ctrl", "mock", "vehicle",
    "dmso", "untreated", "unperturbed", "baseline", "none", "empty_vector",
    "ev", "scramble", "scrambled", "scr", "shctrl", "sh_control",
    "non_targeting", "nontargeting", "nt", "ntc"
  ) |
    grepl("^(wt|ctrl|control|mock|vehicle|dmso)$", y) |
    grepl("(^|_)(wild_type|wildtype|untreated|non_targeting|nontargeting)($|_)", y)
}

classify_condition <- function(x) {
  y <- normalise_label(x)
  fifelse(
    grepl("(^|_)(ko|knockout|knock_out|kd|knockdown|knock_down|sirna|shrna|crispr|null)($|_)", y),
    "loss_of_function",
    fifelse(
      grepl("(^|_)(oe|overexpression|over_expression|rescue|reconstitution|addback)($|_)", y),
      "gain_of_function",
      fifelse(
        grepl("(infect|virus|viral|ifn|interferon|lps|tnf|stress|starv|hypox|treat|drug|stim)", y),
        "treatment_or_stimulus",
        "other_non_control"
      )
    )
  )
}

if (exists("dominant_is_control_condition", mode = "function")) {
  is_control_condition <- dominant_is_control_condition
}
if (exists("dominant_classify_perturbation_class", mode = "function")) {
  classify_condition <- function(x) dominant_classify_perturbation_class(x)
}
if (!exists("classify_design_family", mode = "function")) {
  classify_design_family <- function(condition = NA_character_,
                                     fraction = NA_character_,
                                     inhibitor = NA_character_,
                                     cell_line = NA_character_,
                                     tissue = NA_character_,
                                     gene = NA_character_,
                                     cancer_type = NA_character_) {
    if (exists("dominant_classify_design_family", mode = "function")) {
      dominant_classify_design_family(
        condition = condition,
        fraction = fraction,
        inhibitor = inhibitor,
        cell_line = cell_line,
        tissue = tissue,
        gene = gene,
        cancer_type = cancer_type
      )
    } else {
      classify_condition(condition)
    }
  }
}

safe_wilcox <- function(case_score, control_score, alternative = "two.sided", min_n = 2L) {
  case_score <- case_score[is.finite(case_score)]
  control_score <- control_score[is.finite(control_score)]
  if (length(case_score) < min_n || length(control_score) < min_n) return(NA_real_)
  if (length(unique(c(case_score, control_score))) < 2) return(NA_real_)
  tryCatch(
    wilcox.test(case_score, control_score, alternative = alternative,
                exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

compact_run_signature <- function(runs, max_chars = 240L) {
  runs <- sort(unique(as.character(runs)))
  if (length(runs) == 0) return("")
  joined <- paste(runs, collapse = ";")
  if (nchar(joined) > max_chars) {
    joined <- paste0(substr(joined, 1L, max_chars - 18L), "...;n=", length(runs))
  }
  joined
}

compact_level_label <- function(field, level) {
  field <- as.character(field)
  level <- as.character(level)
  if (grepl(":", field, fixed = TRUE)) {
    fields <- strsplit(field, ":", fixed = TRUE)[[1]]
    parts <- strsplit(level, " | ", fixed = TRUE)[[1]]
    values <- setNames(rep(NA_character_, length(fields)), fields)
    for (part in parts) {
      eq <- regexpr("=", part, fixed = TRUE)[1]
      if (!is.finite(eq) || eq < 1L) next
      key <- substr(part, 1L, eq - 1L)
      value <- substr(part, eq + 1L, nchar(part))
      if (key %chin% names(values)) values[[key]] <- value
    }
    keep <- names(values)[!is.na(values) & nzchar(values)]
    keep <- keep[seq_len(min(length(keep), 3L))]
    compact <- paste0(keep, "=", truncate_label(values[keep], 14L))
    extra <- sum(!is.na(values) & nzchar(values)) - length(keep)
    if (extra > 0L) compact <- c(compact, paste0("+", extra, " more"))
    return(paste(compact, collapse = " | "))
  }
  paste0(field, "=", truncate_label(sub(paste0("^", field, "="), "", level), 26L))
}

feature_annotation <- safe_fread(feature_annotation_file)
relative_usage <- safe_fread(relative_usage_file)
branch_contrasts <- safe_fread(branch_contrast_file)
structure <- safe_fread(structure_file)
feature_expression <- safe_fread(feature_expression_file, required = FALSE)
allocation <- safe_fread(allocation_file, required = FALSE)
control_review <- safe_fread(control_file, required = FALSE)
glmnet_selected <- safe_fread(glmnet_selected_file, required = FALSE)
glmnet_stability <- safe_fread(glmnet_stability_file, required = FALSE)
rf_performance <- safe_fread(rf_performance_file, required = FALSE)
state_dt <- safe_fread(state_file, required = FALSE)
diagnostics <- safe_fread(diagnostics_file, required = FALSE)
known_uorf <- safe_fread(known_uorf_file, required = FALSE)

feature_annotation[, gene_symbol := toupper(gene_symbol)]
relative_usage[, gene_symbol := toupper(gene_symbol)]
branch_contrasts[, gene_symbol := toupper(gene_symbol)]
structure[, gene_symbol := toupper(gene_symbol)]
if (nrow(feature_expression) > 0) feature_expression[, gene_symbol := toupper(gene_symbol)]
if (nrow(allocation) > 0) allocation[, gene_symbol := toupper(gene_symbol)]
if (nrow(diagnostics) > 0) diagnostics[, gene_symbol := toupper(gene_symbol)]

feature_annotation[, feature_key := paste(gene_symbol, tx_id, feature_id, sep = "|")]
relative_usage[, feature_key := paste(gene_symbol, tx_id, feature_id, sep = "|")]
branch_contrasts[, feature_key := paste(gene_symbol, tx_id, feature_id, sep = "|")]

flow <- merge(
  relative_usage,
  feature_annotation[, .(
    feature_key, category, translon_source, translon_sources_merged,
    tx_start, tx_end, feature_bases, cds_overlap_bases,
    start_codon, stop_codon, start_codon_class, frame_relative_to_cds,
    peptide_length_aa, kozak_strength, coding_potential_proxy
  )],
  by = "feature_key",
  all.x = TRUE,
  sort = FALSE
)
flow[, clean_cds_valid := passes_clean_cds_count_cutoff == TRUE &
       is.finite(clean_cds_raw_counts) & clean_cds_raw_counts >= 10]
flow[, feature_start_probability_proxy := pmin(pmax(initiation_probability_proxy, 0), 1)]
flow[, feature_start_probability_proxy := fifelse(
  feature_type == "clean_CDS",
  pmin(pmax(relative_translation, 0), 1),
  feature_start_probability_proxy
)]
setorder(flow, branch_id, gene_symbol, tx_id, feature_order, tx_start)
flow[, ribosome_flow_in_proxy := {
  p <- shift(feature_start_probability_proxy, fill = 0)
  p[!is.finite(p)] <- 0
  cumprod(pmax(0, 1 - p))
}, by = .(branch_id, gene_symbol, tx_id)]
flow[, feature_capture_proxy := fifelse(
  feature_type == "clean_CDS",
  pmin(pmax(relative_translation, 0), 1),
  ribosome_flow_in_proxy * pmin(pmax(feature_start_probability_proxy, 0), 1)
)]
flow[, post_feature_remaining_proxy := pmax(0, ribosome_flow_in_proxy - feature_capture_proxy)]
flow[, post_uorf_loss_proxy := fifelse(
  feature_type == "clean_CDS",
  pmax(0, ribosome_flow_in_proxy - pmin(pmax(relative_translation, 0), 1)),
  pmax(0, ribosome_flow_in_proxy - post_feature_remaining_proxy)
)]
fwrite(flow, flow_output)

compact_edges <- merge(
  branch_contrasts,
  feature_annotation[, .(
    feature_key, category, translon_source, translon_sources_merged,
    tx_start, tx_end, feature_bases, cds_overlap_bases,
    start_codon, start_codon_class, frame_relative_to_cds,
    peptide_length_aa, kozak_strength, coding_potential_proxy
  )],
  by = "feature_key",
  all.x = TRUE,
  sort = FALSE
)
compact_edges <- compact_edges[
  comparison_passes_clean_cds_count_cutoff == TRUE &
    is.finite(delta_relative_translation)
]
compact_edges[, abs_delta_relative_translation := abs(delta_relative_translation)]
setorder(compact_edges, gene_symbol, dominant_state, branch_source, branch_id,
         -abs_delta_relative_translation)
compact_edges <- compact_edges[
  ,
  head(.SD, 1L),
  by = .(gene_symbol, dominant_state, branch_source, branch_id, feature_id)
]
compact_edges[, edge_direction := fifelse(
  delta_relative_translation > 0,
  "IN higher relative use",
  "IN lower relative use"
)]
compact_edges <- compact_edges[, .(
  gene_symbol, tx_id, dominant_state, branch_source, field, level, branch_id,
  feature_id, feature_type, category, translon_sources_merged, tx_start, tx_end,
  feature_bases, cds_overlap_bases, start_codon, start_codon_class,
  frame_relative_to_cds, peptide_length_aa, kozak_strength,
  coding_potential_proxy, n_samples_in, n_samples_out,
  clean_cds_raw_counts_in, clean_cds_raw_counts_out,
  relative_translation_in, relative_translation_out,
  log2_mean_fpkm_like_in_vs_out,
  delta_relative_translation, abs_delta_relative_translation,
  delta_initiation_probability_proxy, edge_direction
)]
fwrite(compact_edges, compact_edge_output)

allocation_summary <- data.table(gene_symbol = unique(structure$gene_symbol))
if (nrow(allocation) > 0) {
  allocation_summary <- allocation[
    scope == "sample",
    .(
      min_sample_allocation_coefficient = safe_min(coefficient),
      max_sample_allocation_abs_coefficient = safe_max(abs(coefficient)),
      strongest_negative_allocation_feature = feature_id[which.min(coefficient)][1],
      strongest_negative_allocation_p_adj = p_adj[which.min(coefficient)][1]
    ),
    by = gene_symbol
  ]
}

edge_summary <- compact_edges[
  ,
  .(
    n_branch_edges = .N,
    max_abs_branch_delta = safe_max(abs_delta_relative_translation),
    max_clean_cds_branch_delta = safe_max(abs_delta_relative_translation[feature_type == "clean_CDS"]),
    max_overlapping_uorf_branch_delta = safe_max(abs_delta_relative_translation[feature_type == "overlapping_uorf"]),
    strongest_branch_state = dominant_state[which.max(abs_delta_relative_translation)][1],
    strongest_branch_field = field[which.max(abs_delta_relative_translation)][1],
    strongest_branch_level = level[which.max(abs_delta_relative_translation)][1],
    strongest_branch_feature = feature_id[which.max(abs_delta_relative_translation)][1],
    strongest_branch_feature_type = feature_type[which.max(abs_delta_relative_translation)][1]
  ),
  by = gene_symbol
]

metadata_confidence <- data.table()
if (nrow(control_review) > 0) {
  metadata_confidence <- copy(control_review)
  metadata_confidence[, confidence_score := fifelse(
    direct_control_available == TRUE &
      control_aware_status %chin% c("supported_by_matched_controls", "supported"),
    2,
    fifelse(direct_control_available == TRUE, 1, 0)
  )]
  metadata_confidence <- metadata_confidence[
    is.finite(adjusted_effect_z) & adjusted_effect_z > 0,
    .(
      max_control_aware_effect_z = safe_max(adjusted_effect_z),
      n_control_aware_supported_terms = sum(confidence_score >= 2, na.rm = TRUE),
      n_direct_control_terms = sum(direct_control_available == TRUE, na.rm = TRUE)
    ),
    by = dominant_state
  ]
}

glmnet_confidence <- data.table()
if (nrow(glmnet_selected) > 0) {
  glmnet_confidence <- copy(glmnet_selected)
  if (nrow(glmnet_stability) > 0 &&
      all(c("dominant_state", "feature_id", "stability_fraction") %chin% names(glmnet_stability))) {
    if ("stability_fraction" %chin% names(glmnet_confidence)) {
      glmnet_confidence[, stability_fraction := NULL]
    }
    glmnet_confidence <- merge(
      glmnet_confidence,
      glmnet_stability[, .(dominant_state, feature_id, stability_fraction)],
      by = c("dominant_state", "feature_id"),
      all.x = TRUE,
      sort = FALSE
    )
  } else if (!"stability_fraction" %chin% names(glmnet_confidence)) {
    glmnet_confidence[, stability_fraction := NA_real_]
  }
  glmnet_confidence <- glmnet_confidence[
    coefficient > 0 & is.finite(effect_z) & effect_z > 0,
    .(
      max_glmnet_effect_z = safe_max(effect_z),
      max_glmnet_stability_fraction = safe_max(stability_fraction),
      n_positive_regularized_terms = .N,
      top_regularized_term = feature_id[which.max(effect_z)][1]
    ),
    by = dominant_state
  ]
}

known_terms <- character()
if (nrow(known_uorf) > 0 && "gene_or_element" %chin% names(known_uorf)) {
  marker_genes <- unique(structure$gene_symbol)
  known_text <- toupper(known_uorf$gene_or_element)
  known_terms <- marker_genes[vapply(marker_genes, function(g) {
    any(grepl(paste0("\\b", g, "\\b"), known_text))
  }, logical(1))]
}

manual_manifest_file <- file.path(analysis_dir, "manual_translons", "manual_translons_manifest.csv")
manual_manifest <- safe_fread(manual_manifest_file, required = FALSE)
manual_genes <- if (nrow(manual_manifest) > 0 && "gene_symbol" %chin% names(manual_manifest)) {
  toupper(unique(manual_manifest$gene_symbol))
} else {
  character()
}

candidate <- merge(structure, allocation_summary, by = "gene_symbol", all.x = TRUE, sort = FALSE)
candidate <- merge(candidate, edge_summary, by = "gene_symbol", all.x = TRUE, sort = FALSE)
candidate[, known_uorf_regulation := gene_symbol %chin% known_terms]
candidate[, manual_m_source_gene := gene_symbol %chin% manual_genes |
            fifelse(is.na(n_uorfs_with_m_support), 0L, n_uorfs_with_m_support) > 0]
candidate[, hidden_or_manual_geometry := hidden_uorf_candidate == TRUE | manual_m_source_gene == TRUE]
candidate[, clean_cds_usage_component := norm01(primary_clean_cds_usage, inverse = TRUE)]
candidate[, branch_shift_component := norm01(max_abs_branch_delta)]
candidate[, allocation_component := norm01(abs(pmin(min_sample_allocation_coefficient, 0)))]
candidate[, overlap_component := as.numeric(has_overlapping_uorf == TRUE)]
candidate[, manual_component := as.numeric(manual_m_source_gene == TRUE)]
candidate[, hidden_component := as.numeric(hidden_uorf_candidate == TRUE)]
candidate[, unknown_component := as.numeric(!known_uorf_regulation)]
for (component_col in c(
  "clean_cds_usage_component", "branch_shift_component", "allocation_component",
  "overlap_component", "manual_component", "hidden_component", "unknown_component"
)) {
  candidate[!is.finite(get(component_col)), (component_col) := 0]
}
candidate[, low_count_penalty := fifelse(
  (!is.na(sample_rows) & sample_rows < 100) |
    (!is.na(branch_rows) & branch_rows < 20),
  1,
  0
)]
candidate[, novel_candidate_score :=
            1.50 * clean_cds_usage_component +
            1.25 * branch_shift_component +
            1.25 * allocation_component +
            1.00 * overlap_component +
            0.90 * manual_component +
            0.60 * hidden_component +
            0.50 * unknown_component -
            0.75 * fifelse(known_uorf_regulation == TRUE, 1, 0) -
            0.50 * low_count_penalty]
candidate[, novelty_class := fcase(
  known_uorf_regulation == TRUE, "known uORF benchmark",
  manual_m_source_gene == TRUE & has_overlapping_uorf == TRUE, "manual-overlap candidate",
  hidden_uorf_candidate == TRUE, "under-annotated candidate",
  has_overlapping_uorf == TRUE, "predicted-overlap candidate",
  default = "structure-associated candidate"
)]
setorder(candidate, -novel_candidate_score, primary_clean_cds_usage)
candidate[, candidate_rank := seq_len(.N)]
candidate <- candidate[, .(
  candidate_rank, gene_symbol, tx_id, novelty_class, known_uorf_regulation,
  manual_m_source_gene, hidden_uorf_candidate, has_overlapping_uorf,
  n_uorfs, n_overlapping_uorfs, n_uorfs_with_m_support,
  primary_clean_cds_usage, clean_cds_delta_vs_count_expected,
  max_abs_branch_delta, max_clean_cds_branch_delta,
  max_overlapping_uorf_branch_delta, min_sample_allocation_coefficient,
  strongest_negative_allocation_feature, strongest_negative_allocation_p_adj,
  strongest_branch_state, strongest_branch_field, strongest_branch_level,
  strongest_branch_feature, strongest_branch_feature_type,
  novel_candidate_score, structure_note, hidden_uorf_candidate_reason
)]
fwrite(candidate, candidate_output)

novel <- candidate[
  known_uorf_regulation == FALSE &
    (manual_m_source_gene == TRUE |
       hidden_uorf_candidate == TRUE |
       has_overlapping_uorf == TRUE |
       novel_candidate_score >= quantile(candidate$novel_candidate_score, 0.75, na.rm = TRUE))
]
setorder(novel, -novel_candidate_score, primary_clean_cds_usage)
fwrite(novel, novel_output)

buffering_effects <- data.table()
buffering_relevant <- data.table()
buffering_condition <- data.table()
buffering_condition_representative <- data.table()
if (nrow(branch_contrasts) > 0) {
  upstream_delta <- branch_contrasts[
    feature_type %chin% c("leader_uorf", "overlapping_uorf") &
      comparison_passes_clean_cds_count_cutoff == TRUE &
      is.finite(delta_relative_translation),
    .(
      n_upstream_features = .N,
      mean_upstream_delta_relative_translation = safe_mean(delta_relative_translation),
      mean_abs_upstream_delta_relative_translation = safe_mean(abs(delta_relative_translation)),
      max_abs_upstream_delta_relative_translation = safe_max(abs(delta_relative_translation))
    ),
    by = .(gene_symbol, tx_id, branch_id)
  ]

  buffering_effects <- branch_contrasts[
    feature_id == "clean_CDS" &
      comparison_passes_clean_cds_count_cutoff == TRUE &
      is.finite(mean_fpkm_like_in) &
      is.finite(mean_fpkm_like_out) &
      is.finite(relative_translation_in) &
      is.finite(relative_translation_out),
    .(
      gene_symbol, tx_id, dominant_state, branch_source, field, level, branch_id,
      n_samples_in, n_samples_out,
      clean_cds_mean_fpkm_like_in = mean_fpkm_like_in,
      clean_cds_mean_fpkm_like_out = mean_fpkm_like_out,
      clean_cds_raw_counts_in, clean_cds_raw_counts_out,
      clean_cds_relative_use_in = relative_translation_in,
      clean_cds_relative_use_out = relative_translation_out,
      clean_cds_relative_use_delta = delta_relative_translation,
      clean_cds_log2_mean_fpkm_like_in_vs_out = log2_mean_fpkm_like_in_vs_out
    )
  ]
  buffering_effects <- merge(
    buffering_effects,
    upstream_delta,
    by = c("gene_symbol", "tx_id", "branch_id"),
    all.x = TRUE,
    sort = FALSE
  )
  buffering_effects <- merge(
    buffering_effects,
    candidate[, .(
      gene_symbol, tx_id, candidate_rank, novelty_class,
      known_uorf_regulation, manual_m_source_gene, novel_candidate_score
    )],
    by = c("gene_symbol", "tx_id"),
    all.x = TRUE,
    sort = FALSE
  )

  buffering_effects[, context_complexity := lengths(strsplit(field, ":", fixed = TRUE))]
  buffering_effects[, condition_value := mapply(
    extract_level_field_value,
    field,
    level,
    MoreArgs = list(target_field = "CONDITION"),
    USE.NAMES = FALSE
  )]
  buffering_effects[, has_condition_context := !is.na(condition_value) & nzchar(condition_value)]

  if (nrow(control_review) > 0) {
    control_exact <- control_review[
      ,
      .(
        field, level, dominant_state,
        exact_direct_control_available = direct_control_available,
        exact_control_aware_status = control_aware_status,
        exact_weighted_score_delta = weighted_score_delta,
        exact_median_effect_z = median_effect_z,
        exact_n_case_samples = n_case_samples,
        exact_n_control_samples = n_control_samples,
        exact_stouffer_greater_p_adj = stouffer_greater_p_adj
      )
    ]
    control_exact[, exact_control_score := control_status_score(exact_control_aware_status)]
    buffering_effects <- merge(
      buffering_effects,
      control_exact,
      by = c("field", "level", "dominant_state"),
      all.x = TRUE,
      sort = FALSE
    )

    condition_control <- control_review[
      field == "CONDITION",
      .(
        condition_direct_control_available = any(direct_control_available == TRUE, na.rm = TRUE),
        condition_control_aware_status = control_aware_status[
          which.max(control_status_score(control_aware_status))
        ][1],
        condition_weighted_score_delta = safe_mean(weighted_score_delta),
        condition_median_effect_z = safe_mean(median_effect_z),
        condition_n_case_samples = safe_max(n_case_samples),
        condition_n_control_samples = safe_max(n_control_samples),
        condition_stouffer_greater_p_adj = safe_min(stouffer_greater_p_adj)
      ),
      by = .(condition_value = level, dominant_state)
    ]
    condition_control[, condition_control_score := control_status_score(condition_control_aware_status)]
    buffering_effects <- merge(
      buffering_effects,
      condition_control,
      by = c("condition_value", "dominant_state"),
      all.x = TRUE,
      sort = FALSE
    )
  }

  for (col in c("exact_control_score", "condition_control_score")) {
    if (!col %chin% names(buffering_effects)) buffering_effects[, (col) := 0]
    buffering_effects[is.na(get(col)), (col) := 0]
  }
  for (col in c("exact_direct_control_available", "condition_direct_control_available")) {
    if (!col %chin% names(buffering_effects)) buffering_effects[, (col) := FALSE]
    buffering_effects[is.na(get(col)), (col) := FALSE]
  }
  buffering_effects[, best_control_score := pmax(exact_control_score, condition_control_score)]
  buffering_effects[, best_control_source := fifelse(
    exact_control_score >= condition_control_score & exact_control_score != 0,
    "exact_branch",
    fifelse(condition_control_score != 0, "condition_level", "none")
  )]
  buffering_effects[, direct_control_available_any :=
                      exact_direct_control_available == TRUE |
                      condition_direct_control_available == TRUE]
  buffering_effects[, matched_control_weighted_score_delta := fifelse(
    best_control_source == "exact_branch",
    exact_weighted_score_delta,
    fifelse(best_control_source == "condition_level", condition_weighted_score_delta, NA_real_)
  )]
  buffering_effects[, matched_control_median_effect_z := fifelse(
    best_control_source == "exact_branch",
    exact_median_effect_z,
    fifelse(best_control_source == "condition_level", condition_median_effect_z, NA_real_)
  )]

  fpkm_denominator_floor <- max(
    5,
    safe_quantile(buffering_effects$clean_cds_mean_fpkm_like_out, 0.10, default = 5)
  )
  buffering_effects[, cds_fpkm_change := clean_cds_mean_fpkm_like_in - clean_cds_mean_fpkm_like_out]
  buffering_effects[, cds_fpkm_change_pct := 100 * safe_div(
    cds_fpkm_change,
    clean_cds_mean_fpkm_like_out
  )]
  buffering_effects[, cds_fpkm_change_pct_stabilized := 100 * safe_div(
    cds_fpkm_change,
    pmax(clean_cds_mean_fpkm_like_out, fpkm_denominator_floor)
  )]
  buffering_effects[, cds_fpkm_stabilized_denominator := pmax(
    clean_cds_mean_fpkm_like_out,
    fpkm_denominator_floor
  )]
  buffering_effects[, low_out_fpkm_like_flag :=
                      clean_cds_mean_fpkm_like_out < fpkm_denominator_floor]
  buffering_effects[, cds_fpkm_log2_fc := log2(
    (clean_cds_mean_fpkm_like_in + 1e-6) /
      (clean_cds_mean_fpkm_like_out + 1e-6)
  )]
  buffering_effects[, cds_fpkm_log2_fc_pseudocount1 := log2(
    (clean_cds_mean_fpkm_like_in + 1) /
      (clean_cds_mean_fpkm_like_out + 1)
  )]
  buffering_effects[, cds_log2_fc_poisson_se_approx := sqrt(
    1 / (pmax(clean_cds_raw_counts_in, 0) + 1) +
      1 / (pmax(clean_cds_raw_counts_out, 0) + 1)
  ) / log(2)]
  buffering_effects[, cds_log2_fc_pseudocount1_lower_approx :=
                      cds_fpkm_log2_fc_pseudocount1 -
                      1.96 * cds_log2_fc_poisson_se_approx]
  buffering_effects[, cds_log2_fc_pseudocount1_upper_approx :=
                      cds_fpkm_log2_fc_pseudocount1 +
                      1.96 * cds_log2_fc_poisson_se_approx]
  buffering_effects[, cds_relative_use_change_pct := 100 * safe_div(
    clean_cds_relative_use_in - clean_cds_relative_use_out,
    clean_cds_relative_use_out
  )]
  buffering_effects[, clean_cds_relative_use_delta_se_approx := sqrt(
    pmax(clean_cds_relative_use_in, 0) *
      pmax(1 - clean_cds_relative_use_in, 0) / pmax(n_samples_in, 1) +
      pmax(clean_cds_relative_use_out, 0) *
      pmax(1 - clean_cds_relative_use_out, 0) / pmax(n_samples_out, 1)
  )]
  buffering_effects[, clean_cds_relative_use_delta_lower_approx :=
                      clean_cds_relative_use_delta -
                      1.96 * clean_cds_relative_use_delta_se_approx]
  buffering_effects[, clean_cds_relative_use_delta_upper_approx :=
                      clean_cds_relative_use_delta +
                      1.96 * clean_cds_relative_use_delta_se_approx]
  buffering_effects[, buffering_ratio_cds_vs_upstream := safe_div(
    abs(clean_cds_relative_use_delta),
    mean_abs_upstream_delta_relative_translation
  )]
  buffering_effects[, buffering_score := pmax(
    0,
    pmin(1, 1 - pmin(buffering_ratio_cds_vs_upstream, 1))
  )]
  buffering_effects[, buffering_interpretation := fifelse(
    !is.finite(buffering_ratio_cds_vs_upstream),
    "insufficient_upstream_signal",
    fifelse(
      buffering_ratio_cds_vs_upstream <= 0.50,
      "buffered_cds_change",
      fifelse(
        buffering_ratio_cds_vs_upstream <= 1.00,
        "partial_buffering",
        "amplified_or_unbuffered_cds_change"
      )
    )
  )]
  buffering_effects[, cds_effect_direction := fifelse(
    cds_fpkm_change > 0,
    "IN higher clean_CDS proxy",
    "IN lower clean_CDS proxy"
  )]
  buffering_effects[, sample_confidence_component := pmin(
    1,
    log10(pmin(n_samples_in, n_samples_out) + 1) / 2
  )]
  buffering_effects[, coverage_confidence_component := pmin(
    1,
    log10(pmin(clean_cds_raw_counts_in, clean_cds_raw_counts_out) + 1) / 4
  )]
  buffering_effects[, control_confidence_component := pmin(
    1,
    pmax(best_control_score, 0) / 2
  )]
  buffering_effects[, effect_confidence_component := pmin(
    1,
    abs(cds_fpkm_log2_fc_pseudocount1) /
      pmax(3 * cds_log2_fc_poisson_se_approx, 1e-6)
  )]
  buffering_effects[, buffering_confidence_score :=
                      0.30 * sample_confidence_component +
                      0.30 * coverage_confidence_component +
                      0.25 * control_confidence_component +
                      0.15 * effect_confidence_component]
  buffering_effects[, effect_size_reliability := fcase(
    low_out_fpkm_like_flag == TRUE & direct_control_available_any == FALSE,
    "low_out_baseline_no_direct_control",
    low_out_fpkm_like_flag == TRUE,
    "low_out_baseline",
    pmin(n_samples_in, n_samples_out) < 10,
    "small_branch_group",
    direct_control_available_any == TRUE & best_control_score >= 1,
    "matched_control_supported",
    default = "descriptive_compendium"
  )]
  buffering_effects[, abs_cds_fpkm_change_pct := abs(cds_fpkm_change_pct)]
  buffering_effects[, abs_cds_fpkm_change_pct_stabilized := abs(cds_fpkm_change_pct_stabilized)]
  buffering_effects[, abs_clean_cds_relative_use_delta := abs(clean_cds_relative_use_delta)]
  buffering_effects[!is.finite(abs_cds_fpkm_change_pct), abs_cds_fpkm_change_pct := -Inf]
  buffering_effects[!is.finite(abs_cds_fpkm_change_pct_stabilized), abs_cds_fpkm_change_pct_stabilized := -Inf]
  buffering_effects[!is.finite(abs_clean_cds_relative_use_delta), abs_clean_cds_relative_use_delta := -Inf]
  setorder(buffering_effects,
           candidate_rank, -buffering_confidence_score,
           -abs_cds_fpkm_change_pct_stabilized, -abs_clean_cds_relative_use_delta)
  fwrite(buffering_effects, buffering_all_output)

  relevant_genes <- unique(c(
    head(candidate$gene_symbol, 40L),
    candidate[manual_m_source_gene == TRUE, gene_symbol],
    candidate[known_uorf_regulation == TRUE, gene_symbol],
    c("ATF4", "PPP1R15A", "MYC")
  ))
  buffering_relevant <- buffering_effects[
    gene_symbol %chin% relevant_genes &
      is.finite(cds_fpkm_change_pct)
  ]
  if (nrow(buffering_relevant) > 0) {
    setorder(buffering_relevant, gene_symbol, -buffering_confidence_score,
             -abs_cds_fpkm_change_pct_stabilized, -abs_clean_cds_relative_use_delta)
    buffering_relevant <- buffering_relevant[
      ,
      head(.SD, 12L),
      by = .(gene_symbol, tx_id)
    ]
    setorder(buffering_relevant, candidate_rank, gene_symbol,
             -buffering_confidence_score,
             -abs_cds_fpkm_change_pct_stabilized, -abs_clean_cds_relative_use_delta)
  }
  fwrite(buffering_relevant, buffering_relevant_output)

  buffering_condition <- buffering_effects[
    grepl("CONDITION", field, fixed = TRUE) &
      is.finite(cds_fpkm_change_pct_stabilized)
  ]
  buffering_condition_representative <- data.table()
  if (nrow(buffering_condition) > 0) {
    buffering_condition[, context_signature := paste(
      gene_symbol, tx_id, dominant_state, sign(cds_fpkm_change),
      n_samples_in, n_samples_out,
      round(clean_cds_mean_fpkm_like_in, 3),
      round(clean_cds_mean_fpkm_like_out, 3),
      round(clean_cds_relative_use_delta, 4),
      sep = "|"
    )]
    buffering_condition_representative <- copy(buffering_condition)
    setorder(
      buffering_condition_representative,
      candidate_rank, -best_control_score, context_complexity,
      -buffering_confidence_score, -abs_cds_fpkm_change_pct_stabilized,
      -abs_clean_cds_relative_use_delta
    )
    buffering_condition_representative <- buffering_condition_representative[
      ,
      head(.SD, 1L),
      by = context_signature
    ]
    setorder(buffering_condition_representative, candidate_rank,
             -buffering_confidence_score, -abs_cds_fpkm_change_pct_stabilized)

    setorder(buffering_condition, gene_symbol, -buffering_confidence_score,
             -abs_cds_fpkm_change_pct_stabilized, -abs_clean_cds_relative_use_delta)
    buffering_condition <- buffering_condition[
      ,
      head(.SD, 8L),
      by = .(gene_symbol, tx_id)
    ]
    setorder(buffering_condition, candidate_rank, gene_symbol,
             -buffering_confidence_score,
             -abs_cds_fpkm_change_pct_stabilized, -abs_clean_cds_relative_use_delta)
  }
  fwrite(buffering_condition, buffering_condition_output)
  fwrite(buffering_condition_representative, buffering_condition_representative_output)
}

matched_feature_contrasts <- data.table()
matched_cds_buffering <- data.table()
matched_cds_bootstrap <- data.table()
matched_gene_summary <- data.table()
matched_design_summary <- data.table()
internal_orf_matched_bump <- data.table()
internal_orf_matched_bump_summary <- data.table()
if (nrow(feature_expression) > 0 && file.exists(metadata_file)) {
  matched_metadata_fields <- c(
    "Run", "CELL_LINE", "TISSUE", "GENE", "CONDITION",
    "study", "AUTHOR", "Cancer_type", "Sex", "INHIBITOR", "FRACTION"
  )
  m <- fread(metadata_file)
  keep_metadata_fields <- intersect(matched_metadata_fields, names(m))
  m <- m[, ..keep_metadata_fields]
  for (field in setdiff(matched_metadata_fields, names(m))) {
    m[, (field) := NA_character_]
  }
  for (field in setdiff(matched_metadata_fields, "Run")) {
    m[, (field) := clean_level(get(field))]
  }

  stratum_fields <- c("study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE")
  for (field in stratum_fields) {
    m[is.na(get(field)), (field) := "MISSING"]
  }
  m[, CONDITION := clean_level(CONDITION)]
  m[, is_control_condition := is_control_condition(CONDITION)]
  m[, design_family := classify_design_family(
    condition = CONDITION,
    fraction = FRACTION,
    inhibitor = INHIBITOR,
    cell_line = CELL_LINE,
    tissue = TISSUE,
    gene = GENE,
    cancer_type = Cancer_type
  )]
  m[, perturbation_class :=
      classify_condition(CONDITION)]
  if (exists("dominant_classify_perturbation_class", mode = "function")) {
    m[, perturbation_class :=
        dominant_classify_perturbation_class(CONDITION, design_family)]
  }
  condition_design_family_audit <- m[
    ,
    .(
      n_samples = .N,
      n_studies = uniqueN(study[!is.na(study)]),
      n_cell_lines = uniqueN(CELL_LINE[!is.na(CELL_LINE)]),
      n_tissues = uniqueN(TISSUE[!is.na(TISSUE)]),
      example_studies =
        paste(head(sort(unique(study[!is.na(study)])), 5L),
              collapse = ";"),
      example_cell_lines =
        paste(head(sort(unique(CELL_LINE[!is.na(CELL_LINE)])), 5L),
              collapse = ";")
    ),
    by = .(CONDITION, is_control_condition, perturbation_class,
           design_family)
  ]
  setorder(condition_design_family_audit, design_family,
           -n_samples, CONDITION)
  fwrite(condition_design_family_audit, condition_design_family_output)
  m[, stratum_id := do.call(paste, c(.SD, sep = " | ")), .SDcols = stratum_fields]

  structure_key <- structure[, .(
    gene_symbol, tx_id,
    has_modeled_uorf = if ("has_modeled_uorf" %chin% names(structure)) {
      has_modeled_uorf
    } else {
      n_uorfs > 0
    },
    n_uorfs, n_overlapping_uorfs, has_overlapping_uorf,
    n_uorfs_with_m_support
  )]

  matched_fx <- feature_expression[
    feature_type %chin% c("leader_uorf", "overlapping_uorf", "clean_CDS")
  ]
  matched_fx <- merge(
    matched_fx,
    structure_key,
    by = c("gene_symbol", "tx_id"),
    all.x = TRUE,
    sort = FALSE
  )
  matched_fx <- matched_fx[has_modeled_uorf == TRUE]
  matched_fx <- merge(matched_fx, m, by = "Run", all.x = TRUE, sort = FALSE)
  matched_fx <- matched_fx[!is.na(CONDITION)]
  matched_fx[, raw_counts := as.numeric(raw_counts)]
  matched_fx[, fpkm_like := as.numeric(fpkm_like)]
  matched_fx[, total_uorf_model_fpkm := sum(fpkm_like, na.rm = TRUE),
             by = .(Run, gene_symbol, tx_id)]
  matched_fx[, relative_use_uorf_model := fifelse(
    total_uorf_model_fpkm > 0,
    fpkm_like / total_uorf_model_fpkm,
    NA_real_
  )]

  contrast_terms <- unique(
    matched_fx[
      is_control_condition != TRUE,
      .(stratum_id, case_condition = CONDITION, perturbation_class,
        design_family)
    ]
  )
  control_strata <- unique(matched_fx[
    is_control_condition == TRUE,
    .(stratum_id)
  ])
  contrast_terms <- merge(
    contrast_terms,
    control_strata[, has_control := TRUE],
    by = "stratum_id",
    all = FALSE,
    sort = FALSE
  )

  feature_cols <- c(
    "gene_symbol", "tx_id", "feature_id", "feature_type", "category",
    "translon_sources_merged", "feature_bases", "cds_overlap_bases",
    "n_uorfs", "n_overlapping_uorfs", "has_overlapping_uorf",
    "n_uorfs_with_m_support"
  )
  matched_feature_contrasts <- rbindlist(
    lapply(seq_len(nrow(contrast_terms)), function(i) {
      term <- contrast_terms[i]
      s <- matched_fx[
        stratum_id == term$stratum_id &
          (is_control_condition == TRUE | CONDITION == term$case_condition)
      ]
      if (nrow(s) == 0) return(NULL)
      s[, matched_group := fifelse(is_control_condition == TRUE, "control", "case")]
      case_runs <- unique(s[matched_group == "case", Run])
      control_runs <- unique(s[matched_group == "control", Run])
      if (length(case_runs) == 0 || length(control_runs) == 0) return(NULL)
      control_conditions <- paste(sort(unique(s[matched_group == "control", CONDITION])), collapse = ";")
      case_run_signature <- compact_run_signature(case_runs)
      control_run_signature <- compact_run_signature(control_runs)
      sample_set_signature <- paste(
        term$stratum_id, term$case_condition,
        paste(sort(case_runs), collapse = ";"),
        paste(sort(control_runs), collapse = ";"),
        sep = " || "
      )

      out <- s[
        ,
        {
          is_case <- matched_group == "case"
          case_rel <- relative_use_uorf_model[is_case]
          control_rel <- relative_use_uorf_model[!is_case]
          case_fpkm <- fpkm_like[is_case]
          control_fpkm <- fpkm_like[!is_case]
          case_raw <- raw_counts[is_case]
          control_raw <- raw_counts[!is_case]
          rel_delta <- mean(case_rel, na.rm = TRUE) - mean(control_rel, na.rm = TRUE)
          fpkm_log2_fc <- log2((mean(case_fpkm, na.rm = TRUE) + 1) /
                                 (mean(control_fpkm, na.rm = TRUE) + 1))
          pooled_rel_sd <- sd(c(case_rel, control_rel), na.rm = TRUE)
          rel_se <- sqrt(
            var(case_rel, na.rm = TRUE) / pmax(uniqueN(Run[is_case]), 1) +
              var(control_rel, na.rm = TRUE) / pmax(uniqueN(Run[!is_case]), 1)
          )
          if (!is.finite(rel_se)) rel_se <- NA_real_
          .(
            n_case_samples = uniqueN(Run[is_case]),
            n_control_samples = uniqueN(Run[!is_case]),
            case_mean_relative_use = mean(case_rel, na.rm = TRUE),
            control_mean_relative_use = mean(control_rel, na.rm = TRUE),
            matched_relative_use_delta = rel_delta,
            matched_relative_use_delta_se = rel_se,
            matched_relative_use_delta_lower = rel_delta - 1.96 * rel_se,
            matched_relative_use_delta_upper = rel_delta + 1.96 * rel_se,
            matched_relative_use_effect_z = if (is.finite(pooled_rel_sd) && pooled_rel_sd > 0) {
              rel_delta / pooled_rel_sd
            } else NA_real_,
            case_mean_fpkm_like = mean(case_fpkm, na.rm = TRUE),
            control_mean_fpkm_like = mean(control_fpkm, na.rm = TRUE),
            matched_fpkm_log2_fc_pseudocount1 = fpkm_log2_fc,
            case_sum_raw_counts = sum(case_raw, na.rm = TRUE),
            control_sum_raw_counts = sum(control_raw, na.rm = TRUE),
            matched_relative_use_wilcox_p = safe_wilcox(case_rel, control_rel, "two.sided"),
            matched_fpkm_wilcox_p = safe_wilcox(case_fpkm, control_fpkm, "two.sided")
          )
        },
        by = feature_cols
      ]
      if (nrow(out) == 0) return(NULL)
      out[, `:=`(
        stratum_id = term$stratum_id,
        case_condition = term$case_condition,
        perturbation_class = term$perturbation_class,
        design_family = term$design_family,
        control_conditions = control_conditions,
        sample_set_signature = sample_set_signature,
        case_run_signature = case_run_signature,
        control_run_signature = control_run_signature
      )]
      for (field in stratum_fields) {
        out[, (field) := s[[field]][1]]
      }
      out
    }),
    fill = TRUE
  )

  if (nrow(matched_feature_contrasts) > 0) {
    matched_feature_contrasts[, `:=`(
      matched_relative_use_wilcox_q = p.adjust(matched_relative_use_wilcox_p, "BH"),
      matched_fpkm_wilcox_q = p.adjust(matched_fpkm_wilcox_p, "BH"),
      min_matched_samples = pmin(n_case_samples, n_control_samples),
      min_feature_raw_counts = pmin(case_sum_raw_counts, control_sum_raw_counts),
      max_feature_raw_counts = pmax(case_sum_raw_counts, control_sum_raw_counts),
      abs_matched_relative_use_delta = abs(matched_relative_use_delta)
    )]
    matched_feature_contrasts[, `:=`(
      matched_feature_passes_sample_gate =
        min_matched_samples >= min_matched_group_samples,
      matched_feature_passes_count_gate_any =
        max_feature_raw_counts >= min_matched_feature_raw_counts_any,
      matched_feature_passes_count_gate_both =
        min_feature_raw_counts >= min_matched_feature_raw_counts_both_strict
    )]
    matched_feature_contrasts[, matched_feature_evaluable :=
                                matched_feature_passes_sample_gate &
                                matched_feature_passes_count_gate_any]
    matched_feature_contrasts[, matched_feature_confidence_score :=
                                fifelse(
                                  matched_feature_evaluable,
                                  0.35 * pmin(1, log10(min_matched_samples + 1) / 1.5) +
                                    0.30 * pmin(1, log10(max_feature_raw_counts + 1) / 3.5) +
                                    0.20 * pmin(1, abs(matched_relative_use_delta) /
                                                  pmax(3 * matched_relative_use_delta_se, 1e-6)) +
                                    0.15 * fifelse(
                                      is.finite(matched_relative_use_wilcox_q) &
                                        matched_relative_use_wilcox_q < 0.10,
                                      1, 0
                                    ),
                                  NA_real_
                                )]
    matched_feature_contrasts[, matched_feature_support := fcase(
      matched_feature_passes_sample_gate != TRUE,
      "single_sample_or_no_replicate",
      matched_feature_passes_count_gate_any != TRUE,
      "low_feature_counts",
      is.finite(matched_relative_use_wilcox_q) &
        matched_relative_use_wilcox_q < 0.10 &
        abs(matched_relative_use_delta) >= 0.10,
      "matched_feature_shift_supported",
      is.finite(matched_fpkm_wilcox_q) &
        matched_fpkm_wilcox_q < 0.10 &
        abs(matched_fpkm_log2_fc_pseudocount1) >= 0.5,
      "matched_feature_fpkm_shift_supported",
      default = "matched_descriptive"
    )]
    safe_setcolorder(
      matched_feature_contrasts,
      c(
        "gene_symbol", "tx_id", "feature_id", "feature_type", "category",
        "stratum_id", stratum_fields, "case_condition", "control_conditions",
        "perturbation_class", "design_family",
        "n_case_samples", "n_control_samples",
        "case_mean_relative_use", "control_mean_relative_use",
        "matched_relative_use_delta", "matched_relative_use_delta_se",
        "matched_relative_use_delta_lower", "matched_relative_use_delta_upper",
        "matched_relative_use_effect_z", "case_mean_fpkm_like",
        "control_mean_fpkm_like", "matched_fpkm_log2_fc_pseudocount1",
        "case_sum_raw_counts", "control_sum_raw_counts",
        "matched_relative_use_wilcox_p", "matched_relative_use_wilcox_q",
        "matched_fpkm_wilcox_p", "matched_fpkm_wilcox_q",
        "min_matched_samples", "min_feature_raw_counts",
        "max_feature_raw_counts", "matched_feature_evaluable",
        "matched_feature_passes_count_gate_any",
        "matched_feature_passes_count_gate_both",
        "matched_feature_confidence_score", "matched_feature_support",
        "sample_set_signature", "case_run_signature", "control_run_signature"
      )
    )
    setorder(
      matched_feature_contrasts,
      gene_symbol, tx_id, case_condition,
      -matched_feature_confidence_score,
      -abs_matched_relative_use_delta
    )
  }

  if (nrow(matched_feature_contrasts) > 0) {
    upstream_matched <- matched_feature_contrasts[
      feature_type != "clean_CDS" &
        matched_feature_evaluable == TRUE &
        is.finite(matched_relative_use_delta),
      .(
        n_upstream_features = .N,
        mean_abs_upstream_relative_use_delta = safe_mean(abs(matched_relative_use_delta)),
        max_abs_upstream_relative_use_delta = safe_max(abs(matched_relative_use_delta)),
        strongest_upstream_feature = feature_id[which.max(abs(matched_relative_use_delta))][1],
        strongest_upstream_feature_type = feature_type[which.max(abs(matched_relative_use_delta))][1]
      ),
      by = .(gene_symbol, tx_id, stratum_id, case_condition, sample_set_signature)
    ]
    matched_cds_buffering <- matched_feature_contrasts[
      feature_type == "clean_CDS"
    ]
    matched_cds_buffering <- merge(
      matched_cds_buffering,
      upstream_matched,
      by = c("gene_symbol", "tx_id", "stratum_id", "case_condition", "sample_set_signature"),
      all.x = TRUE,
      sort = FALSE
    )
    if (nrow(matched_cds_buffering) > 0) {
      fpkm_floor <- max(
        5,
        safe_quantile(matched_cds_buffering$control_mean_fpkm_like, 0.10, default = 5)
      )
      matched_cds_buffering[, matched_cds_fpkm_change :=
                              case_mean_fpkm_like - control_mean_fpkm_like]
      matched_cds_buffering[, matched_cds_fpkm_change_pct_stabilized :=
                              100 * safe_div(
                                matched_cds_fpkm_change,
                                pmax(control_mean_fpkm_like, fpkm_floor)
                              )]
      matched_cds_buffering[, matched_buffering_ratio_cds_vs_upstream :=
                              safe_div(
                                abs(matched_relative_use_delta),
                                mean_abs_upstream_relative_use_delta
                              )]
      matched_cds_buffering[, matched_buffering_interpretation := fifelse(
        !is.finite(matched_buffering_ratio_cds_vs_upstream),
        "insufficient_upstream_signal",
        fifelse(
          matched_buffering_ratio_cds_vs_upstream <= 0.50,
          "matched_buffered_cds_change",
          fifelse(
            matched_buffering_ratio_cds_vs_upstream <= 1.00,
            "matched_partial_buffering",
            "matched_amplified_or_unbuffered_cds_change"
          )
        )
      )]
      matched_cds_buffering[, matched_cds_support := fcase(
        matched_feature_passes_sample_gate != TRUE,
        "single_sample_or_no_replicate",
        matched_feature_passes_count_gate_any != TRUE,
        "low_clean_cds_counts",
        is.finite(matched_relative_use_wilcox_q) &
          matched_relative_use_wilcox_q < 0.10 &
          abs(matched_relative_use_delta) >= 0.10,
        "matched_clean_cds_relative_shift_supported",
        is.finite(matched_fpkm_wilcox_q) &
          matched_fpkm_wilcox_q < 0.10 &
          abs(matched_fpkm_log2_fc_pseudocount1) >= 0.5,
        "matched_clean_cds_fpkm_shift_supported",
        default = "matched_clean_cds_descriptive"
      )]
      matched_cds_buffering[, abs_matched_cds_fpkm_change_pct_stabilized :=
                              abs(matched_cds_fpkm_change_pct_stabilized)]
      matched_cds_buffering[, abs_matched_clean_cds_relative_use_delta :=
                              abs(matched_relative_use_delta)]
      setorder(
        matched_cds_buffering,
        gene_symbol, -matched_feature_confidence_score,
        -abs_matched_cds_fpkm_change_pct_stabilized,
        -abs_matched_clean_cds_relative_use_delta
      )
    }
  }

  internal_bump_fx <- feature_expression[
    feature_type == "internal_orf" &
      is.finite(feature_vs_cds_mean_ratio) &
      is.finite(cds_reference_raw_counts) &
      cds_reference_raw_counts >= 10
  ]
  if (nrow(internal_bump_fx) > 0) {
    internal_bump_fx <- merge(
      internal_bump_fx,
      m,
      by = "Run",
      all.x = TRUE,
      sort = FALSE
    )
    internal_bump_fx <- internal_bump_fx[!is.na(CONDITION)]
    internal_bump_fx[, `:=`(
      raw_counts = as.numeric(raw_counts),
      feature_vs_cds_mean_ratio = as.numeric(feature_vs_cds_mean_ratio),
      internal_log2_ratio = log2(pmax(feature_vs_cds_mean_ratio, 0) + 0.10),
      internal_ratio_has_bump = feature_vs_cds_mean_ratio > 1
    )]

    internal_bump_terms <- unique(
      internal_bump_fx[
        is_control_condition != TRUE,
        .(stratum_id, case_condition = CONDITION, perturbation_class,
          design_family)
      ]
    )
    internal_bump_controls <- unique(internal_bump_fx[
      is_control_condition == TRUE,
      .(stratum_id)
    ])
    internal_bump_terms <- merge(
      internal_bump_terms,
      internal_bump_controls[, has_control := TRUE],
      by = "stratum_id",
      all = FALSE,
      sort = FALSE
    )
    internal_bump_feature_cols <- c(
      "gene_symbol", "tx_id", "feature_id", "feature_type", "category",
      "translon_sources_merged", "feature_bases", "cds_overlap_bases",
      "measured_bases", "measurement_scope"
    )

    internal_orf_matched_bump <- rbindlist(
      lapply(seq_len(nrow(internal_bump_terms)), function(i) {
        term <- internal_bump_terms[i]
        s <- internal_bump_fx[
          stratum_id == term$stratum_id &
            (is_control_condition == TRUE | CONDITION == term$case_condition)
        ]
        if (nrow(s) == 0) return(NULL)
        s[, matched_group := fifelse(
          is_control_condition == TRUE,
          "control",
          "case"
        )]
        case_runs <- unique(s[matched_group == "case", Run])
        control_runs <- unique(s[matched_group == "control", Run])
        if (length(case_runs) == 0 || length(control_runs) == 0) return(NULL)
        control_conditions <- paste(
          sort(unique(s[matched_group == "control", CONDITION])),
          collapse = ";"
        )
        case_run_signature <- compact_run_signature(case_runs)
        control_run_signature <- compact_run_signature(control_runs)
        sample_set_signature <- paste(
          term$stratum_id, term$case_condition,
          paste(sort(case_runs), collapse = ";"),
          paste(sort(control_runs), collapse = ";"),
          sep = " || "
        )

        out <- s[
          ,
          {
            is_case <- matched_group == "case"
            case_ratio <- feature_vs_cds_mean_ratio[is_case]
            control_ratio <- feature_vs_cds_mean_ratio[!is_case]
            case_log_ratio <- internal_log2_ratio[is_case]
            control_log_ratio <- internal_log2_ratio[!is_case]
            case_bump <- internal_ratio_has_bump[is_case]
            control_bump <- internal_ratio_has_bump[!is_case]
            .(
              n_case_samples = uniqueN(Run[is_case]),
              n_control_samples = uniqueN(Run[!is_case]),
              case_mean_feature_vs_cds_ratio = safe_mean(case_ratio),
              control_mean_feature_vs_cds_ratio = safe_mean(control_ratio),
              case_median_feature_vs_cds_ratio = safe_quantile(case_ratio, 0.50),
              control_median_feature_vs_cds_ratio =
                safe_quantile(control_ratio, 0.50),
              matched_feature_vs_cds_ratio_delta =
                safe_mean(case_ratio) - safe_mean(control_ratio),
              case_bump_fraction = safe_mean(as.numeric(case_bump)),
              control_bump_fraction = safe_mean(as.numeric(control_bump)),
              matched_bump_fraction_delta =
                safe_mean(as.numeric(case_bump)) -
                safe_mean(as.numeric(control_bump)),
              matched_log2_feature_vs_cds_ratio_shift =
                safe_mean(case_log_ratio) - safe_mean(control_log_ratio),
              case_sum_raw_counts = sum(raw_counts[is_case], na.rm = TRUE),
              control_sum_raw_counts = sum(raw_counts[!is_case], na.rm = TRUE),
              matched_ratio_wilcox_p =
                safe_wilcox(case_log_ratio, control_log_ratio, "two.sided")
            )
          },
          by = internal_bump_feature_cols
        ]
        if (nrow(out) == 0) return(NULL)
        out[, `:=`(
          stratum_id = term$stratum_id,
          case_condition = term$case_condition,
          control_conditions = control_conditions,
          perturbation_class = term$perturbation_class,
          design_family = term$design_family,
          sample_set_signature = sample_set_signature,
          case_run_signature = case_run_signature,
          control_run_signature = control_run_signature
        )]
        for (field in stratum_fields) {
          out[, (field) := s[[field]][1]]
        }
        out
      }),
      fill = TRUE
    )
  }

  if (nrow(internal_orf_matched_bump) > 0) {
    internal_orf_matched_bump[
      ,
      matched_ratio_wilcox_q := p.adjust(matched_ratio_wilcox_p, "BH")
    ]
    internal_orf_matched_bump[, `:=`(
      min_matched_samples = pmin(n_case_samples, n_control_samples),
      min_internal_raw_counts = pmin(case_sum_raw_counts, control_sum_raw_counts),
      max_internal_raw_counts = pmax(case_sum_raw_counts, control_sum_raw_counts),
      internal_bump_passes_sample_gate =
        pmin(n_case_samples, n_control_samples) >= min_matched_group_samples,
      internal_bump_passes_count_gate_any =
        pmax(case_sum_raw_counts, control_sum_raw_counts) >=
          min_matched_feature_raw_counts_any,
      internal_bump_passes_count_gate_both =
        pmin(case_sum_raw_counts, control_sum_raw_counts) >=
          min_matched_feature_raw_counts_both_strict,
      internal_condition_bump_confidence_score =
        fifelse(
          pmin(n_case_samples, n_control_samples) >= min_matched_group_samples &
            pmax(case_sum_raw_counts, control_sum_raw_counts) >=
              min_matched_feature_raw_counts_any,
          0.30 * pmin(1, log10(pmin(n_case_samples, n_control_samples) + 1) / 1.5) +
            0.25 * pmin(1, log10(pmax(case_sum_raw_counts, control_sum_raw_counts) + 1) / 3) +
            0.20 * pmin(1, pmax(case_bump_fraction, 0)) +
            0.15 * pmin(
              1,
              pmax(matched_log2_feature_vs_cds_ratio_shift, 0) / 1
            ) +
            0.10 * fifelse(
              is.finite(matched_ratio_wilcox_q) & matched_ratio_wilcox_q < 0.10,
              1,
              0
            ),
          NA_real_
        )
    )]
    internal_orf_matched_bump[, internal_condition_bump_class := fcase(
      internal_bump_passes_sample_gate != TRUE,
      "single_sample_or_no_replicate",
      case_sum_raw_counts < min_matched_feature_raw_counts_any,
      "low_case_internal_counts",
      is.finite(case_mean_feature_vs_cds_ratio) &
        case_mean_feature_vs_cds_ratio > 1 &
        case_bump_fraction >= 0.50 &
        is.finite(matched_ratio_wilcox_q) &
        matched_ratio_wilcox_q < 0.10 &
        (matched_log2_feature_vs_cds_ratio_shift >= 0.50 |
           matched_bump_fraction_delta >= 0.20),
      "matched_case_cds_bump_supported",
      is.finite(case_mean_feature_vs_cds_ratio) &
        case_mean_feature_vs_cds_ratio > 1 &
        case_bump_fraction >= 0.25 &
        (matched_log2_feature_vs_cds_ratio_shift >= 0.25 |
           matched_bump_fraction_delta >= 0.10),
      "matched_case_cds_bump_review",
      is.finite(case_mean_feature_vs_cds_ratio) &
        case_mean_feature_vs_cds_ratio > 1,
      "case_cds_bump_descriptive",
      default = "no_case_cds_bump"
    )]
    safe_setcolorder(
      internal_orf_matched_bump,
      c(
        "gene_symbol", "tx_id", "feature_id", "feature_type", "category",
        "internal_condition_bump_class",
        "internal_condition_bump_confidence_score",
        "stratum_id", stratum_fields, "case_condition", "control_conditions",
        "perturbation_class", "design_family",
        "n_case_samples", "n_control_samples",
        "case_mean_feature_vs_cds_ratio", "control_mean_feature_vs_cds_ratio",
        "case_bump_fraction", "control_bump_fraction",
        "matched_bump_fraction_delta",
        "matched_log2_feature_vs_cds_ratio_shift",
        "matched_ratio_wilcox_p", "matched_ratio_wilcox_q",
        "case_sum_raw_counts", "control_sum_raw_counts",
        "min_matched_samples", "min_internal_raw_counts",
        "max_internal_raw_counts", "internal_bump_passes_count_gate_any",
        "internal_bump_passes_count_gate_both",
        "sample_set_signature", "case_run_signature", "control_run_signature"
      )
    )
    setorder(
      internal_orf_matched_bump,
      -internal_condition_bump_confidence_score,
      -case_bump_fraction,
      -matched_log2_feature_vs_cds_ratio_shift,
      gene_symbol,
      feature_id
    )

    internal_orf_matched_bump_summary <- internal_orf_matched_bump[
      ,
      {
        supported <- internal_condition_bump_class ==
          "matched_case_cds_bump_supported"
        review <- internal_condition_bump_class %chin% c(
          "matched_case_cds_bump_supported",
          "matched_case_cds_bump_review"
        )
        idx <- order(
          -as.integer(supported),
          -as.integer(review),
          -internal_condition_bump_confidence_score,
          -pmax(matched_log2_feature_vs_cds_ratio_shift, 0),
          na.last = TRUE
        )[1]
        .(
          n_internal_matched_bump_rows = .N,
          n_internal_condition_bump_review_rows = sum(review, na.rm = TRUE),
          n_internal_condition_bump_supported_rows =
            sum(supported, na.rm = TRUE),
          n_internal_condition_bump_review_studies =
            uniqueN(study[review]),
          n_internal_condition_bump_supported_studies =
            uniqueN(study[supported]),
          max_internal_condition_bump_case_ratio =
            safe_max(case_mean_feature_vs_cds_ratio),
          max_internal_condition_bump_log2_shift =
            safe_max(matched_log2_feature_vs_cds_ratio_shift),
          best_internal_condition_bump_class =
            internal_condition_bump_class[idx][1],
          best_internal_condition_bump_score =
            internal_condition_bump_confidence_score[idx][1],
          best_internal_condition_bump_case_condition =
            case_condition[idx][1],
          best_internal_condition_bump_context =
            paste(stratum_id[idx][1], case_condition[idx][1], sep = " || ")
        )
      },
      by = .(gene_symbol, tx_id, feature_id)
    ]
  }

  if (nrow(matched_cds_buffering) > 0) {
    set.seed(37)
    cds_sample <- matched_fx[
      feature_type == "clean_CDS",
      .(
        Run, gene_symbol, tx_id, feature_id, stratum_id, CONDITION,
        is_control_condition, relative_use_uorf_model, fpkm_like
      )
    ]
    setkey(cds_sample, gene_symbol, tx_id, stratum_id)
    bootstrap_keys <- unique(matched_cds_buffering[, .(
      gene_symbol, tx_id, feature_id, stratum_id, case_condition,
      sample_set_signature
    )])
    matched_cds_bootstrap <- rbindlist(
      lapply(seq_len(nrow(bootstrap_keys)), function(i) {
        key <- bootstrap_keys[i]
        s <- cds_sample[.(key$gene_symbol, key$tx_id, key$stratum_id)]
        if (nrow(s) == 0) return(NULL)
        is_case <- s$CONDITION == key$case_condition
        is_control <- s$is_control_condition == TRUE
        out <- bootstrap_cds_contrast(
          s$relative_use_uorf_model[is_case],
          s$relative_use_uorf_model[is_control],
          s$fpkm_like[is_case],
          s$fpkm_like[is_control]
        )
        out[, `:=`(
          gene_symbol = key$gene_symbol,
          tx_id = key$tx_id,
          feature_id = key$feature_id,
          stratum_id = key$stratum_id,
          case_condition = key$case_condition,
          sample_set_signature = key$sample_set_signature
        )]
        out
      }),
      fill = TRUE
    )
    if (nrow(matched_cds_bootstrap) > 0) {
      safe_setcolorder(
        matched_cds_bootstrap,
        c(
          "gene_symbol", "tx_id", "feature_id", "stratum_id",
          "case_condition", "sample_set_signature", "bootstrap_iterations"
        )
      )
      matched_cds_buffering <- merge(
        matched_cds_buffering,
        matched_cds_bootstrap,
        by = c(
          "gene_symbol", "tx_id", "feature_id", "stratum_id",
          "case_condition", "sample_set_signature"
        ),
        all.x = TRUE,
        sort = FALSE
      )
    }

    matched_cds_buffering[, bootstrap_relative_use_ci_excludes_zero :=
                            is.finite(bootstrap_relative_use_delta_lower) &
                            is.finite(bootstrap_relative_use_delta_upper) &
                            (bootstrap_relative_use_delta_lower > 0 |
                               bootstrap_relative_use_delta_upper < 0)]
    matched_cds_buffering[, bootstrap_fpkm_log2_fc_ci_excludes_zero :=
                            is.finite(bootstrap_fpkm_log2_fc_lower) &
                            is.finite(bootstrap_fpkm_log2_fc_upper) &
                            (bootstrap_fpkm_log2_fc_lower > 0 |
                               bootstrap_fpkm_log2_fc_upper < 0)]
    rel_tau2 <- var(matched_cds_buffering$matched_relative_use_delta, na.rm = TRUE)
    if (!is.finite(rel_tau2) || rel_tau2 <= 0) rel_tau2 <- 0.01
    rel_se2 <- fifelse(
      is.finite(matched_cds_buffering$bootstrap_relative_use_delta_sd),
      matched_cds_buffering$bootstrap_relative_use_delta_sd^2,
      matched_cds_buffering$matched_relative_use_delta_se^2
    )
    matched_cds_buffering[, relative_use_empirical_bayes_weight :=
                            rel_tau2 / (rel_tau2 + pmax(rel_se2, 1e-8))]
    matched_cds_buffering[, shrunk_matched_relative_use_delta :=
                            relative_use_empirical_bayes_weight *
                            matched_relative_use_delta]

    logfc_tau2 <- var(matched_cds_buffering$matched_fpkm_log2_fc_pseudocount1, na.rm = TRUE)
    if (!is.finite(logfc_tau2) || logfc_tau2 <= 0) logfc_tau2 <- 0.25
    logfc_se2 <- matched_cds_buffering$bootstrap_fpkm_log2_fc_sd^2
    matched_cds_buffering[, fpkm_log2_fc_empirical_bayes_weight :=
                            logfc_tau2 / (logfc_tau2 + pmax(logfc_se2, 1e-8))]
    matched_cds_buffering[!is.finite(fpkm_log2_fc_empirical_bayes_weight),
                          fpkm_log2_fc_empirical_bayes_weight := 0]
    matched_cds_buffering[, shrunk_matched_fpkm_log2_fc :=
                            fpkm_log2_fc_empirical_bayes_weight *
                            matched_fpkm_log2_fc_pseudocount1]
    matched_cds_buffering[, matched_bootstrap_certainty_score := pmax(
      pmin(1, abs(shrunk_matched_relative_use_delta) /
             pmax(3 * bootstrap_relative_use_delta_sd, 1e-6)),
      pmin(1, abs(shrunk_matched_fpkm_log2_fc) /
             pmax(3 * bootstrap_fpkm_log2_fc_sd, 1e-6)),
      fifelse(bootstrap_relative_use_ci_excludes_zero == TRUE |
                bootstrap_fpkm_log2_fc_ci_excludes_zero == TRUE, 0.75, 0),
      na.rm = TRUE
    )]
    matched_cds_buffering[!is.finite(matched_bootstrap_certainty_score),
                          matched_bootstrap_certainty_score := 0]

    design_fields <- c("study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE")
    matched_cds_buffering[, matched_design_signature := paste(
      paste0(design_fields, "=", unlist(.SD)),
      collapse = " | "
    ), by = seq_len(nrow(matched_cds_buffering)), .SDcols = design_fields]
    matched_cds_buffering[, matched_design_signature := paste(
      matched_design_signature,
      paste0("CONDITION=", case_condition),
      sep = " | "
    )]
    matched_cds_buffering[, matched_study_signature := paste(study, AUTHOR, sep = " | ")]
    supported_rows <- matched_cds_buffering[
      matched_cds_support %chin% c(
        "matched_clean_cds_relative_shift_supported",
        "matched_clean_cds_fpkm_shift_supported"
      )
    ]
    design_gene_counts <- supported_rows[
      ,
      .(supported_genes_in_design = uniqueN(gene_symbol)),
      by = matched_design_signature
    ]
    matched_cds_buffering <- merge(
      matched_cds_buffering,
      design_gene_counts,
      by = "matched_design_signature",
      all.x = TRUE,
      sort = FALSE
    )
    matched_cds_buffering[is.na(supported_genes_in_design),
                          supported_genes_in_design := 0L]
    matched_cds_buffering[, matched_design_specificity :=
                            fifelse(
                              supported_genes_in_design > 0,
                              1 / sqrt(supported_genes_in_design),
                              0
                            )]

    matched_gene_summary <- matched_cds_buffering[
      ,
      {
        supported <- matched_cds_support %chin% c(
          "matched_clean_cds_relative_shift_supported",
          "matched_clean_cds_fpkm_shift_supported"
        )
        support_rank <- fcase(
          supported == TRUE, 2,
          matched_cds_support == "matched_clean_cds_descriptive", 1,
          default = 0
        )
        best_idx <- order(
          -support_rank,
          -matched_feature_confidence_score,
          -abs_matched_cds_fpkm_change_pct_stabilized
        )[1]
        .(
          n_matched_control_cds_rows = .N,
          n_supported_matched_control_cds_rows = sum(supported, na.rm = TRUE),
          n_supported_matched_control_designs =
            uniqueN(matched_design_signature[supported]),
          n_supported_matched_control_studies =
            uniqueN(matched_study_signature[supported]),
          max_abs_matched_control_cds_fpkm_change_pct_stabilized =
            safe_max(abs_matched_cds_fpkm_change_pct_stabilized),
          max_abs_matched_control_clean_cds_relative_use_delta =
            safe_max(abs_matched_clean_cds_relative_use_delta),
          max_matched_control_confidence_score =
            safe_max(matched_feature_confidence_score),
          max_matched_control_bootstrap_certainty_score =
            safe_max(matched_bootstrap_certainty_score),
          max_abs_shrunk_matched_control_relative_use_delta =
            safe_max(abs(shrunk_matched_relative_use_delta)),
          max_abs_shrunk_matched_control_fpkm_log2_fc =
            safe_max(abs(shrunk_matched_fpkm_log2_fc)),
          max_supported_design_specificity =
            safe_max(matched_design_specificity[supported]),
          max_supported_design_supported_gene_count =
            safe_max(supported_genes_in_design[supported]),
          best_matched_control_context = paste(
            stratum_id[best_idx][1],
            case_condition[best_idx][1],
            sep = " | case="
          ),
          best_matched_control_support =
            matched_cds_support[best_idx][1],
          best_matched_buffering_interpretation =
            matched_buffering_interpretation[best_idx][1]
        )
      },
      by = .(gene_symbol, tx_id)
    ]

    matched_design_summary <- matched_cds_buffering[
      ,
      .(
        n_clean_cds_rows = .N,
        n_supported_clean_cds_rows = sum(
          matched_cds_support %chin% c(
            "matched_clean_cds_relative_shift_supported",
            "matched_clean_cds_fpkm_shift_supported"
          ),
          na.rm = TRUE
        ),
        n_genes = uniqueN(gene_symbol),
        n_supported_genes = uniqueN(
          gene_symbol[
            matched_cds_support %chin% c(
              "matched_clean_cds_relative_shift_supported",
              "matched_clean_cds_fpkm_shift_supported"
            )
          ]
        ),
        mean_abs_shrunk_clean_cds_relative_use_delta =
          safe_mean(abs(shrunk_matched_relative_use_delta)),
        mean_abs_shrunk_clean_cds_fpkm_log2_fc =
          safe_mean(abs(shrunk_matched_fpkm_log2_fc)),
        mean_bootstrap_certainty_score =
          safe_mean(matched_bootstrap_certainty_score),
        top_supported_genes = paste(
          head(unique(gene_symbol[
            matched_cds_support %chin% c(
              "matched_clean_cds_relative_shift_supported",
              "matched_clean_cds_fpkm_shift_supported"
            )
          ]), 12L),
          collapse = ";"
        )
      ),
      by = .(matched_design_signature, matched_study_signature)
    ]
    setorder(
      matched_design_summary,
      -n_supported_genes,
      -mean_bootstrap_certainty_score,
      -mean_abs_shrunk_clean_cds_relative_use_delta
    )
  }
}
fwrite(matched_feature_contrasts, matched_feature_contrast_output)
fwrite(matched_cds_buffering, matched_cds_buffering_output)
fwrite(matched_cds_bootstrap, matched_cds_bootstrap_output)
fwrite(matched_gene_summary, matched_gene_summary_output)
fwrite(matched_design_summary, matched_design_summary_output)
fwrite(internal_orf_matched_bump, internal_orf_matched_bump_output)

generic_condition_evidence <- data.table()
if (nrow(matched_cds_buffering) > 0) {
  generic_condition_evidence <- copy(matched_cds_buffering)
  generic_condition_evidence[, matched_support_row :=
                               matched_cds_support %chin% c(
                                 "matched_clean_cds_relative_shift_supported",
                                 "matched_clean_cds_fpkm_shift_supported"
                               )]
  generic_condition_evidence[, bootstrap_interval_support :=
                               bootstrap_relative_use_ci_excludes_zero == TRUE |
                               bootstrap_fpkm_log2_fc_ci_excludes_zero == TRUE]
  generic_condition_evidence[, strict_matched_support_row :=
                               matched_support_row == TRUE &
                               bootstrap_interval_support == TRUE]
  generic_condition_evidence <- generic_condition_evidence[
    ,
    {
      weights <- pmax(matched_bootstrap_certainty_score, 0.05)
      weighted_log2_fc <- if (all(!is.finite(shrunk_matched_fpkm_log2_fc))) {
        NA_real_
      } else {
        weighted.mean(
          shrunk_matched_fpkm_log2_fc[is.finite(shrunk_matched_fpkm_log2_fc)],
          weights[is.finite(shrunk_matched_fpkm_log2_fc)],
          na.rm = TRUE
        )
      }
      review_order <- order(
        -as.integer(strict_matched_support_row),
        -as.integer(matched_support_row),
        -matched_bootstrap_certainty_score,
        -abs(shrunk_matched_fpkm_log2_fc)
      )
      .(
        n_matched_rows = .N,
        n_matched_designs = uniqueN(matched_design_signature),
        n_matched_studies = uniqueN(matched_study_signature),
        n_supported_rows = sum(matched_support_row, na.rm = TRUE),
        n_supported_designs =
          uniqueN(matched_design_signature[matched_support_row == TRUE]),
        n_supported_studies =
          uniqueN(matched_study_signature[matched_support_row == TRUE]),
        n_strict_rows = sum(strict_matched_support_row, na.rm = TRUE),
        n_strict_designs =
          uniqueN(matched_design_signature[strict_matched_support_row == TRUE]),
        n_strict_studies =
          uniqueN(matched_study_signature[strict_matched_support_row == TRUE]),
        weighted_shrunk_clean_cds_log2_fc = weighted_log2_fc,
        weighted_shrunk_clean_cds_pct_change = if (is.finite(weighted_log2_fc)) {
          100 * (2^weighted_log2_fc - 1)
        } else {
          NA_real_
        },
        max_abs_shrunk_clean_cds_log2_fc =
          safe_max(abs(shrunk_matched_fpkm_log2_fc)),
        max_bootstrap_certainty_score =
          safe_max(matched_bootstrap_certainty_score),
        strongest_case_condition = case_condition[review_order][1],
        strongest_control_conditions = control_conditions[review_order][1],
        strongest_design_signature = matched_design_signature[review_order][1],
        strongest_matched_support = matched_cds_support[review_order][1],
        strongest_buffering_interpretation =
          matched_buffering_interpretation[review_order][1]
      )
    },
    by = .(gene_symbol, tx_id, perturbation_class)
  ]
  generic_condition_evidence[, generic_condition_evidence_gate := fcase(
    n_strict_designs >= 2 & n_strict_studies >= 2,
    "generic_replicated_strict_interval_review",
    n_supported_designs >= 2 & n_supported_studies >= 2,
    "generic_replicated_matched_review",
    n_supported_rows > 0,
    "generic_supported_single_context_review",
    n_matched_rows > 0,
    "generic_descriptive_matched_context",
    default = "no_generic_condition_evidence"
  )]
  generic_condition_evidence[, generic_condition_evidence_score :=
                               pmin(n_strict_designs, 3) +
                               0.75 * pmin(n_supported_designs, 3) +
                               0.50 * pmin(n_strict_studies, 2) +
                               pmin(max_bootstrap_certainty_score, 1) +
                               pmin(max_abs_shrunk_clean_cds_log2_fc, 2) / 2]
  setorder(
    generic_condition_evidence,
    -generic_condition_evidence_score,
    -n_strict_designs,
    -n_supported_designs,
    gene_symbol,
    perturbation_class
  )
}
fwrite(generic_condition_evidence, generic_condition_evidence_output)

candidate_v2 <- copy(candidate)
buffering_gene_summary <- data.table()
if (nrow(buffering_effects) > 0) {
  buffering_gene_summary <- buffering_effects[
    is.finite(cds_fpkm_change_pct_stabilized),
    .(
      n_buffering_rows = .N,
      n_condition_buffering_rows = sum(has_condition_context == TRUE, na.rm = TRUE),
      n_direct_control_buffering_rows = sum(direct_control_available_any == TRUE, na.rm = TRUE),
      n_supported_control_buffering_rows = sum(best_control_score >= 2, na.rm = TRUE),
      max_abs_cds_fpkm_change_pct_stabilized = safe_max(abs(cds_fpkm_change_pct_stabilized)),
      max_abs_clean_cds_relative_use_delta = safe_max(abs(clean_cds_relative_use_delta)),
      max_buffering_confidence_score = safe_max(buffering_confidence_score),
      best_buffering_ratio_cds_vs_upstream = buffering_ratio_cds_vs_upstream[
        which.max(buffering_confidence_score)
      ][1],
      best_buffering_interpretation = buffering_interpretation[
        which.max(buffering_confidence_score)
      ][1],
      best_buffering_effect_size_reliability = effect_size_reliability[
        which.max(buffering_confidence_score)
      ][1],
      best_buffering_context = paste(
        field[which.max(buffering_confidence_score)][1],
        level[which.max(buffering_confidence_score)][1],
        sep = "="
      )
    ),
    by = .(gene_symbol, tx_id)
  ]
}
candidate_v2 <- merge(candidate_v2, buffering_gene_summary,
                      by = c("gene_symbol", "tx_id"), all.x = TRUE, sort = FALSE)
if (nrow(matched_gene_summary) > 0) {
  candidate_v2 <- merge(candidate_v2, matched_gene_summary,
                        by = c("gene_symbol", "tx_id"), all.x = TRUE, sort = FALSE)
}
for (col in c(
  "n_buffering_rows", "n_condition_buffering_rows", "n_direct_control_buffering_rows",
  "n_supported_control_buffering_rows", "max_abs_cds_fpkm_change_pct_stabilized",
  "max_abs_clean_cds_relative_use_delta", "max_buffering_confidence_score",
  "n_matched_control_cds_rows", "n_supported_matched_control_cds_rows",
  "n_supported_matched_control_designs", "n_supported_matched_control_studies",
  "max_abs_matched_control_cds_fpkm_change_pct_stabilized",
  "max_abs_matched_control_clean_cds_relative_use_delta",
  "max_matched_control_confidence_score",
  "max_matched_control_bootstrap_certainty_score",
  "max_abs_shrunk_matched_control_relative_use_delta",
  "max_abs_shrunk_matched_control_fpkm_log2_fc",
  "max_supported_design_specificity",
  "max_supported_design_supported_gene_count"
)) {
  if (!col %chin% names(candidate_v2)) candidate_v2[, (col) := 0]
  candidate_v2[is.na(get(col)) | !is.finite(get(col)), (col) := 0]
}
candidate_v2[, preliminary_candidate_rank := candidate_rank]
candidate_v2[, buffering_effect_component := norm01(max_abs_cds_fpkm_change_pct_stabilized)]
candidate_v2[, buffering_confidence_component_v2 := pmin(pmax(max_buffering_confidence_score, 0), 1)]
candidate_v2[, direct_control_component := pmin(n_direct_control_buffering_rows, 1)]
candidate_v2[, supported_control_component := pmin(n_supported_control_buffering_rows, 1)]
for (component_col in c(
  "buffering_effect_component", "buffering_confidence_component_v2",
  "direct_control_component", "supported_control_component"
)) {
  candidate_v2[!is.finite(get(component_col)), (component_col) := 0]
}
candidate_v2[, next_iteration_candidate_score :=
               novel_candidate_score +
               0.75 * buffering_effect_component +
               0.65 * buffering_confidence_component_v2 +
               0.40 * direct_control_component +
               0.30 * supported_control_component]
candidate_v2[, next_iteration_candidate_score_v2 := next_iteration_candidate_score]
setorder(candidate_v2, -next_iteration_candidate_score, primary_clean_cds_usage)
candidate_v2[, candidate_rank_v2 := seq_len(.N)]
safe_setcolorder(candidate_v2, c("candidate_rank_v2", "preliminary_candidate_rank"))
candidate_v2_ranked <- copy(candidate_v2)
fwrite(candidate_v2_ranked, candidate_v2_output)

candidate_v3 <- copy(candidate_v2_ranked)
candidate_v3[, matched_control_effect_component :=
               norm01(max_abs_matched_control_cds_fpkm_change_pct_stabilized)]
candidate_v3[, matched_control_relative_component :=
               norm01(max_abs_matched_control_clean_cds_relative_use_delta)]
candidate_v3[, matched_control_confidence_component :=
               pmin(pmax(max_matched_control_confidence_score, 0), 1)]
candidate_v3[, matched_supported_control_component :=
               pmin(n_supported_matched_control_cds_rows, 1)]
for (component_col in c(
  "matched_control_effect_component", "matched_control_relative_component",
  "matched_control_confidence_component", "matched_supported_control_component"
)) {
  candidate_v3[!is.finite(get(component_col)), (component_col) := 0]
}
candidate_v3[, next_iteration_candidate_score_v3 :=
               next_iteration_candidate_score_v2 +
               0.55 * matched_control_effect_component +
               0.45 * matched_control_relative_component +
               0.55 * matched_control_confidence_component +
               0.45 * matched_supported_control_component]
candidate_v3[, next_iteration_candidate_score := next_iteration_candidate_score_v3]
setorder(candidate_v3, -next_iteration_candidate_score_v3, primary_clean_cds_usage)
candidate_v3[, candidate_rank_v3 := seq_len(.N)]
safe_setcolorder(
  candidate_v3,
  c("candidate_rank_v3", "candidate_rank_v2", "preliminary_candidate_rank")
)
fwrite(candidate_v3, candidate_v3_output)

candidate_v4 <- copy(candidate_v3)
candidate_v4[, design_replication_component :=
               pmin(1, log1p(n_supported_matched_control_designs) / log1p(3))]
candidate_v4[, study_replication_component :=
               pmin(1, log1p(n_supported_matched_control_studies) / log1p(2))]
candidate_v4[, design_specificity_component :=
               pmin(pmax(max_supported_design_specificity, 0), 1)]
candidate_v4[, bootstrap_certainty_component :=
               pmin(pmax(max_matched_control_bootstrap_certainty_score, 0), 1)]
candidate_v4[, shrunk_matched_effect_component :=
               norm01(max_abs_shrunk_matched_control_fpkm_log2_fc)]
candidate_v4[, broad_single_design_penalty := fifelse(
  n_supported_matched_control_designs <= 1 &
    max_supported_design_supported_gene_count >= 10,
  0.35,
  0
)]
for (component_col in c(
  "design_replication_component", "study_replication_component",
  "design_specificity_component", "bootstrap_certainty_component",
  "shrunk_matched_effect_component", "broad_single_design_penalty"
)) {
  candidate_v4[!is.finite(get(component_col)), (component_col) := 0]
}
candidate_v4[, next_iteration_candidate_score_v4 :=
               next_iteration_candidate_score_v3 +
               0.45 * design_replication_component +
               0.30 * study_replication_component +
               0.35 * design_specificity_component +
               0.35 * bootstrap_certainty_component +
               0.30 * shrunk_matched_effect_component -
               broad_single_design_penalty]
candidate_v4[, next_iteration_candidate_score := next_iteration_candidate_score_v4]
setorder(candidate_v4, -next_iteration_candidate_score_v4, primary_clean_cds_usage)
candidate_v4[, candidate_rank_v4 := seq_len(.N)]
candidate_v4[, candidate_rank := candidate_rank_v4]
safe_setcolorder(
  candidate_v4,
  c("candidate_rank_v4", "candidate_rank",
    "candidate_rank_v3",
    "candidate_rank_v2", "preliminary_candidate_rank")
)
candidate_v2 <- candidate_v4

novel <- candidate_v2[
  known_uorf_regulation == FALSE &
    (manual_m_source_gene == TRUE |
       hidden_uorf_candidate == TRUE |
       has_overlapping_uorf == TRUE |
       next_iteration_candidate_score >= quantile(next_iteration_candidate_score, 0.75, na.rm = TRUE))
]
setorder(novel, -next_iteration_candidate_score, primary_clean_cds_usage)
fwrite(candidate_v2, candidate_output)
fwrite(candidate_v2, candidate_v4_output)
fwrite(novel, novel_output)

internal_orf_initiation_evidence <- data.table()
internal_orf_initiation_summary <- data.table()
internal_signal_columns <- c(
  "internal_frame_relative_to_cds",
  "internal_phase_total_raw_counts",
  "internal_feature_phase0_raw_counts",
  "internal_canonical_phase_raw_counts",
  "internal_iorf_frame_fraction",
  "internal_canonical_frame_fraction",
  "internal_iorf_vs_canonical_frame_delta",
  "internal_start_window_bases",
  "internal_feature_core_bases",
  "internal_upstream_window_bases",
  "internal_start_window_raw_counts",
  "internal_feature_core_raw_counts",
  "internal_upstream_window_raw_counts",
  "internal_start_vs_body_density_ratio",
  "internal_start_vs_upstream_density_ratio"
)
if (nrow(feature_expression) > 0 &&
    "feature_type" %chin% names(feature_expression)) {
  for (column in setdiff(internal_signal_columns, names(feature_expression))) {
    feature_expression[, (column) := NA]
  }
  internal_signal_dt <- feature_expression[feature_type == "internal_orf"]
  if (nrow(internal_signal_dt) > 0) {
    internal_orf_initiation_summary <- internal_signal_dt[
      ,
      {
        frame_ok <- is.finite(internal_phase_total_raw_counts) &
          internal_phase_total_raw_counts >= 10 &
          is.finite(internal_iorf_vs_canonical_frame_delta)
        start_body_ok <- is.finite(internal_start_vs_body_density_ratio) &
          internal_start_window_raw_counts + internal_feature_core_raw_counts >= 10 &
          internal_feature_core_bases >= 12
        start_upstream_ok <- is.finite(internal_start_vs_upstream_density_ratio) &
          internal_start_window_raw_counts + internal_upstream_window_raw_counts >= 10 &
          internal_upstream_window_bases >= 12
        .(
          measured_internal_frame_relative_to_cds = {
            x <- internal_frame_relative_to_cds[
              is.finite(internal_frame_relative_to_cds)
            ]
            if (length(x) == 0) NA_integer_ else as.integer(x[[1]])
          },
          n_internal_frame_evaluable_runs = sum(frame_ok, na.rm = TRUE),
          total_internal_phase_raw_counts = sum(
            internal_phase_total_raw_counts,
            na.rm = TRUE
          ),
          weighted_internal_iorf_frame_fraction = safe_div(
            sum(internal_feature_phase0_raw_counts, na.rm = TRUE),
            sum(internal_phase_total_raw_counts, na.rm = TRUE)
          ),
          weighted_internal_canonical_frame_fraction = safe_div(
            sum(internal_canonical_phase_raw_counts, na.rm = TRUE),
            sum(internal_phase_total_raw_counts, na.rm = TRUE)
          ),
          median_internal_iorf_frame_fraction = safe_quantile(
            internal_iorf_frame_fraction[frame_ok],
            0.50
          ),
          median_internal_canonical_frame_fraction = safe_quantile(
            internal_canonical_frame_fraction[frame_ok],
            0.50
          ),
          median_internal_iorf_vs_canonical_frame_delta = safe_quantile(
            internal_iorf_vs_canonical_frame_delta[frame_ok],
            0.50
          ),
          n_internal_iorf_frame_advantage_runs = sum(
            frame_ok & internal_iorf_vs_canonical_frame_delta >= 0.10,
            na.rm = TRUE
          ),
          n_internal_canonical_frame_advantage_runs = sum(
            frame_ok & internal_iorf_vs_canonical_frame_delta <= -0.10,
            na.rm = TRUE
          ),
          n_internal_start_body_evaluable_runs = sum(start_body_ok, na.rm = TRUE),
          median_internal_start_vs_body_ratio = safe_quantile(
            internal_start_vs_body_density_ratio[start_body_ok],
            0.50
          ),
          q90_internal_start_vs_body_ratio = safe_quantile(
            internal_start_vs_body_density_ratio[start_body_ok],
            0.90
          ),
          n_internal_start_body_peak_runs = sum(
            start_body_ok & internal_start_vs_body_density_ratio >= 1.50,
            na.rm = TRUE
          ),
          n_internal_start_upstream_evaluable_runs = sum(
            start_upstream_ok,
            na.rm = TRUE
          ),
          median_internal_start_vs_upstream_ratio = safe_quantile(
            internal_start_vs_upstream_density_ratio[start_upstream_ok],
            0.50
          ),
          q90_internal_start_vs_upstream_ratio = safe_quantile(
            internal_start_vs_upstream_density_ratio[start_upstream_ok],
            0.90
          )
        )
      },
      by = .(gene_symbol, tx_id, feature_id)
    ]
    internal_orf_initiation_summary[
      ,
      internal_iorf_frame_advantage_run_fraction := safe_div(
        n_internal_iorf_frame_advantage_runs,
        n_internal_frame_evaluable_runs
      )
    ]
    internal_orf_initiation_summary[
      ,
      internal_canonical_frame_advantage_run_fraction := safe_div(
        n_internal_canonical_frame_advantage_runs,
        n_internal_frame_evaluable_runs
      )
    ]
    internal_orf_initiation_summary[
      ,
      internal_start_body_peak_run_fraction := safe_div(
        n_internal_start_body_peak_runs,
        n_internal_start_body_evaluable_runs
      )
    ]
    internal_orf_initiation_summary[, internal_frame_evidence_class := fcase(
      measured_internal_frame_relative_to_cds == 0,
      "same_frame_not_discriminating",
      !is.finite(total_internal_phase_raw_counts) |
        total_internal_phase_raw_counts < 100 |
        n_internal_frame_evaluable_runs < 3,
      "low_internal_frame_counts",
      is.finite(weighted_internal_iorf_frame_fraction) &
        is.finite(weighted_internal_canonical_frame_fraction) &
        weighted_internal_iorf_frame_fraction -
          weighted_internal_canonical_frame_fraction >= 0.10 &
        internal_iorf_frame_advantage_run_fraction >= 0.50,
      "iorf_frame_advantage_review",
      is.finite(weighted_internal_iorf_frame_fraction) &
        is.finite(weighted_internal_canonical_frame_fraction) &
        weighted_internal_canonical_frame_fraction -
          weighted_internal_iorf_frame_fraction >= 0.10 &
        internal_canonical_frame_advantage_run_fraction >= 0.50,
      "canonical_frame_like",
      default = "mixed_or_weak_frame_signal"
    )]
    internal_orf_initiation_summary[, internal_start_shape_class := fcase(
      n_internal_start_body_evaluable_runs < 3,
      "low_start_shape_counts",
      is.finite(median_internal_start_vs_body_ratio) &
        median_internal_start_vs_body_ratio >= 1.50 &
        internal_start_body_peak_run_fraction >= 0.50,
      "recurrent_start_proximal_peak_review",
      is.finite(q90_internal_start_vs_body_ratio) &
        q90_internal_start_vs_body_ratio >= 2 &
        internal_start_body_peak_run_fraction >= 0.20,
      "intermittent_start_proximal_peak_review",
      default = "no_start_proximal_peak"
    )]
    internal_orf_initiation_evidence <- copy(internal_orf_initiation_summary)
    if (nrow(feature_annotation) > 0) {
      internal_orf_initiation_evidence <- merge(
        feature_annotation[
          feature_type == "internal_orf",
          .(
            gene_symbol, tx_id, feature_id, feature_type, category,
            translon_source, translon_sources_merged,
            tx_start, tx_end, feature_bases, cds_overlap_bases,
            start_codon, start_codon_class, frame_relative_to_cds,
            peptide_length_aa, kozak_strength
          )
        ],
        internal_orf_initiation_evidence,
        by = c("gene_symbol", "tx_id", "feature_id"),
        all.y = TRUE,
        sort = FALSE
      )
    }
    safe_setcolorder(
      internal_orf_initiation_evidence,
      c(
        "gene_symbol", "tx_id", "feature_id", "feature_type", "category",
        "internal_frame_evidence_class", "internal_start_shape_class",
        "translon_sources_merged", "tx_start", "tx_end", "feature_bases",
        "start_codon", "start_codon_class", "frame_relative_to_cds",
        "weighted_internal_iorf_frame_fraction",
        "weighted_internal_canonical_frame_fraction",
        "median_internal_iorf_vs_canonical_frame_delta",
        "n_internal_frame_evaluable_runs",
        "internal_iorf_frame_advantage_run_fraction",
        "internal_canonical_frame_advantage_run_fraction",
        "median_internal_start_vs_body_ratio",
        "q90_internal_start_vs_body_ratio",
        "internal_start_body_peak_run_fraction"
      )
    )
  }
}
fwrite(internal_orf_initiation_evidence, internal_orf_initiation_output)

internal_orf_architecture <- data.table()
if (nrow(feature_annotation) > 0 &&
    "feature_type" %chin% names(feature_annotation)) {
  internal_orf_architecture <- feature_annotation[
    feature_type == "internal_orf",
    .(
      gene_symbol, tx_id, feature_id, feature_type, category,
      translon_source, translon_sources_merged,
      tx_start, tx_end, feature_bases, cds_overlap_bases,
      start_codon, stop_codon, start_codon_class,
      frame_relative_to_cds, peptide_length_aa,
      kozak_minus3, kozak_plus4, kozak_strength,
      coding_potential_proxy, feature_order
    )
  ]
}
if (nrow(internal_orf_architecture) > 0) {
  internal_branch_summary <- data.table()
  if (nrow(branch_contrasts) > 0) {
    internal_branch_summary <- branch_contrasts[
      feature_type == "internal_orf" &
        comparison_passes_clean_cds_count_cutoff == TRUE &
        is.finite(delta_relative_translation)
    ]
    if (nrow(internal_branch_summary) > 0) {
      internal_branch_summary[, abs_delta_relative_translation := abs(delta_relative_translation)]
      internal_branch_summary <- internal_branch_summary[
        ,
        {
          idx <- which.max(abs_delta_relative_translation)
          .(
            n_internal_branch_rows = .N,
            n_internal_branch_states = uniqueN(dominant_state),
            max_abs_internal_branch_delta = abs_delta_relative_translation[idx][1],
            strongest_internal_branch_delta = delta_relative_translation[idx][1],
            strongest_internal_branch_log2_feature_change =
              log2_mean_fpkm_like_in_vs_out[idx][1],
            strongest_internal_branch_state = dominant_state[idx][1],
            strongest_internal_branch_field = field[idx][1],
            strongest_internal_branch_level = level[idx][1],
            strongest_internal_branch_id = branch_id[idx][1],
            strongest_internal_branch_label = branch_label[idx][1]
          )
        },
        by = .(gene_symbol, tx_id, feature_id)
      ]
    }
  }

  internal_expression_summary <- data.table()
  if (nrow(feature_expression) > 0) {
    internal_expression_summary <- feature_expression[
      feature_type == "internal_orf",
      .(
        n_internal_expression_runs = uniqueN(Run),
        n_internal_runs_with_counts = sum(raw_counts > 0, na.rm = TRUE),
        total_internal_raw_counts = sum(raw_counts, na.rm = TRUE),
        median_internal_fpkm_like = median(fpkm_like, na.rm = TRUE),
        max_internal_fpkm_like = safe_max(fpkm_like),
        n_internal_cds_reference_runs = sum(
          is.finite(feature_vs_cds_mean_ratio) &
            is.finite(cds_reference_raw_counts) &
            cds_reference_raw_counts >= 10,
          na.rm = TRUE
        ),
        n_internal_cds_bump_runs = sum(
          is.finite(feature_vs_cds_mean_ratio) &
            feature_vs_cds_mean_ratio > 1 &
            is.finite(cds_reference_raw_counts) &
            cds_reference_raw_counts >= 10,
          na.rm = TRUE
        ),
        median_internal_vs_cds_mean_ratio = safe_quantile(
          feature_vs_cds_mean_ratio[
            is.finite(cds_reference_raw_counts) &
              cds_reference_raw_counts >= 10
          ],
          0.50
        ),
        q90_internal_vs_cds_mean_ratio = safe_quantile(
          feature_vs_cds_mean_ratio[
            is.finite(cds_reference_raw_counts) &
              cds_reference_raw_counts >= 10
          ],
          0.90
        ),
        max_internal_vs_cds_mean_ratio = safe_max(
          feature_vs_cds_mean_ratio[
            is.finite(cds_reference_raw_counts) &
              cds_reference_raw_counts >= 10
          ]
        )
      ),
      by = .(gene_symbol, tx_id, feature_id)
    ]
  }

  internal_orf_architecture <- merge(
    internal_orf_architecture,
    internal_branch_summary,
    by = c("gene_symbol", "tx_id", "feature_id"),
    all.x = TRUE,
    sort = FALSE
  )
  internal_orf_architecture <- merge(
    internal_orf_architecture,
    internal_expression_summary,
    by = c("gene_symbol", "tx_id", "feature_id"),
    all.x = TRUE,
    sort = FALSE
  )
  internal_orf_architecture <- merge(
    internal_orf_architecture,
    candidate_v2[, .(
      gene_symbol, tx_id, candidate_rank_v4, novelty_class,
      primary_clean_cds_usage, next_iteration_candidate_score_v4,
      n_uorfs, n_overlapping_uorfs, has_overlapping_uorf
    )],
    by = c("gene_symbol", "tx_id"),
    all.x = TRUE,
    sort = FALSE
  )
  if (nrow(internal_orf_matched_bump_summary) > 0) {
    internal_orf_architecture <- merge(
      internal_orf_architecture,
      internal_orf_matched_bump_summary,
      by = c("gene_symbol", "tx_id", "feature_id"),
      all.x = TRUE,
      sort = FALSE
    )
  }
  if (nrow(internal_orf_initiation_summary) > 0) {
    internal_orf_architecture <- merge(
      internal_orf_architecture,
      internal_orf_initiation_summary,
      by = c("gene_symbol", "tx_id", "feature_id"),
      all.x = TRUE,
      sort = FALSE
    )
  }
  internal_bump_count_cols <- c(
    "n_internal_matched_bump_rows",
    "n_internal_condition_bump_review_rows",
    "n_internal_condition_bump_supported_rows",
    "n_internal_condition_bump_review_studies",
    "n_internal_condition_bump_supported_studies"
  )
  for (column in internal_bump_count_cols) {
    if (!column %chin% names(internal_orf_architecture)) {
      internal_orf_architecture[, (column) := 0L]
    }
    internal_orf_architecture[is.na(get(column)), (column) := 0L]
  }

  internal_orf_architecture[, internal_orf_evidence_score := 0L]
  internal_orf_architecture[
    start_codon_class %chin% c("ATG", "near_cognate"),
    internal_orf_evidence_score := internal_orf_evidence_score + 1L
  ]
  internal_orf_architecture[
    kozak_strength %chin% c("strong", "moderate"),
    internal_orf_evidence_score := internal_orf_evidence_score + 1L
  ]
  internal_orf_architecture[
    is.finite(peptide_length_aa) & peptide_length_aa >= 10,
    internal_orf_evidence_score := internal_orf_evidence_score + 1L
  ]
  internal_orf_architecture[
    is.finite(frame_relative_to_cds) & frame_relative_to_cds != 0,
    internal_orf_evidence_score := internal_orf_evidence_score + 1L
  ]
  internal_orf_architecture[
    is.finite(max_abs_internal_branch_delta) &
      max_abs_internal_branch_delta >= 0.10,
    internal_orf_evidence_score := internal_orf_evidence_score + 1L
  ]
  internal_orf_architecture[
    is.finite(strongest_internal_branch_log2_feature_change) &
      abs(strongest_internal_branch_log2_feature_change) >= 0.50,
    internal_orf_evidence_score := internal_orf_evidence_score + 1L
  ]
  internal_orf_architecture[
    grepl("\\bM\\b", fifelse(
      is.na(translon_sources_merged),
      "",
      translon_sources_merged
    )),
    internal_orf_evidence_score := internal_orf_evidence_score + 1L
  ]
  internal_orf_architecture[, internal_cds_signal_bump_supported :=
                              is.finite(q90_internal_vs_cds_mean_ratio) &
                              q90_internal_vs_cds_mean_ratio > 1 &
                              n_internal_cds_bump_runs > 0]
  internal_orf_architecture[, internal_cds_signal_bump_fraction := safe_div(
    n_internal_cds_bump_runs,
    n_internal_cds_reference_runs
  )]
  internal_orf_architecture[, internal_compendium_bump_class := fcase(
    internal_cds_signal_bump_supported != TRUE,
    "no_compendium_cds_bump",
    is.finite(median_internal_vs_cds_mean_ratio) &
      median_internal_vs_cds_mean_ratio >= 1 &
      internal_cds_signal_bump_fraction >= 0.50,
    "recurrent_median_cds_bump",
    internal_cds_signal_bump_fraction >= 0.25,
    "intermittent_cds_bump",
    default = "sparse_upper_tail_cds_bump"
  )]
  internal_orf_architecture[
    internal_cds_signal_bump_supported == TRUE,
    internal_orf_evidence_score := internal_orf_evidence_score + 1L
  ]
  internal_orf_architecture[
    n_internal_condition_bump_supported_rows > 0,
    internal_orf_evidence_score := internal_orf_evidence_score + 1L
  ]
  internal_orf_architecture[
    internal_frame_evidence_class == "iorf_frame_advantage_review",
    internal_orf_evidence_score := internal_orf_evidence_score + 1L
  ]
  internal_orf_architecture[, internal_orf_tc_only :=
                              grepl("\\bTC\\b", fifelse(
                                is.na(translon_sources_merged),
                                "",
                                translon_sources_merged
                              )) &
                              !grepl("\\b(T|M)\\b", fifelse(
                                is.na(translon_sources_merged),
                                "",
                                translon_sources_merged
                              ))]
  internal_orf_architecture[, internal_orf_architecture_class := fifelse(
    is.finite(frame_relative_to_cds) & frame_relative_to_cds == 0,
    "in_frame_cds_shared_coverage_flag",
    fifelse(
      internal_frame_evidence_class == "canonical_frame_like",
      "canonical_frame_dominant_coverage_review",
      fifelse(
      internal_orf_tc_only == TRUE &
        internal_cds_signal_bump_supported != TRUE,
      "tc_iORF_without_cds_signal_bump",
      fifelse(
        internal_orf_tc_only == TRUE &
          internal_compendium_bump_class == "sparse_upper_tail_cds_bump" &
          !(n_internal_condition_bump_supported_rows > 0),
        "tc_iORF_sparse_compendium_bump_review",
        fifelse(
      internal_orf_evidence_score >= 5,
      "high_priority_iORF_review",
      fifelse(
        internal_orf_evidence_score >= 3,
        "medium_priority_iORF_review",
        "low_evidence_internal_coverage_flag"
      )
        )
      )
      )
    )
  )]
  internal_orf_architecture[, candidate_score_component :=
                              norm01(next_iteration_candidate_score_v4)]
  internal_orf_architecture[!is.finite(candidate_score_component),
                            candidate_score_component := 0]
  internal_orf_architecture[, internal_orf_priority_score :=
                              internal_orf_evidence_score +
                              2 * fifelse(
                                is.finite(max_abs_internal_branch_delta),
                                max_abs_internal_branch_delta,
                                0
                              ) +
                              candidate_score_component]
  setorder(
    internal_orf_architecture,
    -internal_orf_priority_score,
    candidate_rank_v4,
    gene_symbol,
    feature_order
  )
  safe_setcolorder(
    internal_orf_architecture,
    c(
      "gene_symbol", "tx_id", "feature_id", "internal_orf_architecture_class",
      "internal_orf_priority_score", "internal_orf_evidence_score",
      "candidate_rank_v4", "novelty_class", "feature_type", "category",
      "tx_start", "tx_end", "feature_bases", "cds_overlap_bases",
      "start_codon", "start_codon_class", "kozak_strength",
      "frame_relative_to_cds", "peptide_length_aa",
      "translon_sources_merged", "max_abs_internal_branch_delta",
      "strongest_internal_branch_state", "strongest_internal_branch_field",
      "strongest_internal_branch_level", "strongest_internal_branch_label",
      "n_internal_expression_runs", "n_internal_runs_with_counts",
      "total_internal_raw_counts", "median_internal_fpkm_like",
      "internal_orf_tc_only", "internal_cds_signal_bump_supported",
      "internal_frame_evidence_class", "internal_start_shape_class",
      "weighted_internal_iorf_frame_fraction",
      "weighted_internal_canonical_frame_fraction",
      "median_internal_iorf_vs_canonical_frame_delta",
      "n_internal_frame_evaluable_runs",
      "internal_iorf_frame_advantage_run_fraction",
      "internal_canonical_frame_advantage_run_fraction",
      "median_internal_start_vs_body_ratio",
      "q90_internal_start_vs_body_ratio",
      "internal_start_body_peak_run_fraction",
      "n_internal_cds_reference_runs", "n_internal_cds_bump_runs",
      "internal_cds_signal_bump_fraction", "internal_compendium_bump_class",
      "median_internal_vs_cds_mean_ratio", "q90_internal_vs_cds_mean_ratio",
      "n_internal_condition_bump_review_rows",
      "n_internal_condition_bump_supported_rows",
      "n_internal_condition_bump_review_studies",
      "n_internal_condition_bump_supported_studies",
      "best_internal_condition_bump_class",
      "best_internal_condition_bump_context"
    )
  )
}
fwrite(internal_orf_architecture, internal_orf_architecture_output)

feature_confidence <- copy(feature_annotation)
if (nrow(feature_expression) > 0) {
  ratio_cols <- c(
    "feature_vs_cds_mean_ratio", "cds_reference_raw_counts", "measurement_scope",
    "raw_counts", "measured_bases"
  )
  for (column in setdiff(ratio_cols, names(feature_expression))) {
    feature_expression[, (column) := NA]
  }
  feature_measurement_summary <- feature_expression[
    ,
    .(
      expression_runs = uniqueN(Run),
      runs_with_feature_counts = sum(raw_counts > 0, na.rm = TRUE),
      total_feature_raw_counts = sum(raw_counts, na.rm = TRUE),
      median_measured_bases = safe_quantile(measured_bases, 0.50),
      measurement_scope = {
        values <- unique(as.character(measurement_scope))
        values <- values[!is.na(values) & nzchar(values)]
        if (length(values) == 0) NA_character_ else paste(values, collapse = ";")
      },
      cds_reference_runs = sum(
        is.finite(cds_reference_raw_counts) & cds_reference_raw_counts >= 10,
        na.rm = TRUE
      ),
      feature_vs_cds_mean_ratio_median = safe_quantile(
        feature_vs_cds_mean_ratio[
          is.finite(cds_reference_raw_counts) & cds_reference_raw_counts >= 10
        ],
        0.50
      ),
      feature_vs_cds_mean_ratio_q90 = safe_quantile(
        feature_vs_cds_mean_ratio[
          is.finite(cds_reference_raw_counts) & cds_reference_raw_counts >= 10
        ],
        0.90
      ),
      runs_with_cds_signal_bump = sum(
        is.finite(feature_vs_cds_mean_ratio) &
          feature_vs_cds_mean_ratio > 1 &
          is.finite(cds_reference_raw_counts) &
          cds_reference_raw_counts >= 10,
        na.rm = TRUE
      )
    ),
    by = .(gene_symbol, tx_id, feature_id)
  ]
  feature_confidence <- merge(
    feature_confidence,
    feature_measurement_summary,
    by = c("gene_symbol", "tx_id", "feature_id"),
    all.x = TRUE,
    sort = FALSE
  )
}
if (nrow(diagnostics) > 0) {
  diag_columns <- intersect(
    c(
      "gene_symbol", "tx_id", "is_canonical_isoform",
      "canonical_isoform_qc_priority", "canonical_isoform_source",
      "cds_exon_qc_status", "cds_exon_qc_measurable",
      "cds_exon_qc_transcript_counts", "cds_exon_qc_n_exons",
      "cds_exon_qc_min_ratio", "cds_exon_qc_low_exon_ranks",
      "clean_cds_bases", "excluded_cds_bases", "expression_source"
    ),
    names(diagnostics)
  )
  feature_confidence <- merge(
    feature_confidence,
    unique(diagnostics[, ..diag_columns]),
    by = c("gene_symbol", "tx_id"),
    all.x = TRUE,
    sort = FALSE
  )
}
if (nrow(internal_orf_matched_bump_summary) > 0) {
  feature_confidence <- merge(
    feature_confidence,
    internal_orf_matched_bump_summary,
    by = c("gene_symbol", "tx_id", "feature_id"),
    all.x = TRUE,
    sort = FALSE
  )
}
if (nrow(internal_orf_initiation_summary) > 0) {
  feature_confidence <- merge(
    feature_confidence,
    internal_orf_initiation_summary,
    by = c("gene_symbol", "tx_id", "feature_id"),
    all.x = TRUE,
    sort = FALSE
  )
}
for (column in c(
  "n_internal_matched_bump_rows",
  "n_internal_condition_bump_review_rows",
  "n_internal_condition_bump_supported_rows",
  "n_internal_condition_bump_review_studies",
  "n_internal_condition_bump_supported_studies"
)) {
  if (!column %chin% names(feature_confidence)) {
    feature_confidence[, (column) := 0L]
  }
  feature_confidence[is.na(get(column)), (column) := 0L]
}
feature_confidence[, translon_sources_text := fifelse(
  is.na(translon_sources_merged),
  fifelse(is.na(translon_source), "", as.character(translon_source)),
  as.character(translon_sources_merged)
)]
feature_confidence[, internal_orf_tc_only :=
                     feature_type == "internal_orf" &
                     grepl("\\bTC\\b", translon_sources_text) &
                     !grepl("\\b(T|M)\\b", translon_sources_text)]
feature_confidence[, cds_signal_bump_supported :=
                     is.finite(feature_vs_cds_mean_ratio_q90) &
                     feature_vs_cds_mean_ratio_q90 > 1 &
                     runs_with_cds_signal_bump > 0]
feature_confidence[, cds_signal_bump_fraction := safe_div(
  runs_with_cds_signal_bump,
  cds_reference_runs
)]
feature_confidence[, internal_compendium_bump_class := fcase(
  feature_type != "internal_orf",
  NA_character_,
  cds_signal_bump_supported != TRUE,
  "no_compendium_cds_bump",
  is.finite(feature_vs_cds_mean_ratio_median) &
    feature_vs_cds_mean_ratio_median >= 1 &
    cds_signal_bump_fraction >= 0.50,
  "recurrent_median_cds_bump",
  cds_signal_bump_fraction >= 0.25,
  "intermittent_cds_bump",
  default = "sparse_upper_tail_cds_bump"
)]
feature_confidence[, feature_measurement_class := fcase(
  feature_id == "clean_CDS",
  "clean_cds",
  feature_type == "internal_orf" & internal_orf_tc_only == TRUE &
    cds_signal_bump_supported != TRUE,
  "tc_iORF_without_cds_signal_bump",
  feature_type == "internal_orf" & internal_orf_tc_only == TRUE &
    internal_compendium_bump_class == "sparse_upper_tail_cds_bump" &
    n_internal_condition_bump_supported_rows == 0,
  "tc_iORF_sparse_compendium_bump_review",
  feature_type == "internal_orf" &
    internal_frame_evidence_class == "canonical_frame_like",
  "diagnostic_iORF_canonical_frame_dominant",
  feature_type == "internal_orf" &
    internal_frame_evidence_class == "iorf_frame_advantage_review",
  "diagnostic_iORF_with_feature_frame_advantage",
  feature_type == "internal_orf" &
    n_internal_condition_bump_supported_rows > 0,
  "diagnostic_iORF_with_matched_condition_bump",
  feature_type == "internal_orf" & cds_signal_bump_supported == TRUE,
  "diagnostic_iORF_with_cds_signal_bump",
  feature_type == "internal_orf",
  "diagnostic_iORF_without_cds_signal_bump",
  feature_type == "nte_extension",
  "unique_nte_extension_branch_evidence",
  feature_type == "ntt_truncation",
  "cds_shared_ntt_branch_evidence",
  feature_type %chin% c("leader_uorf", "overlapping_uorf"),
  "upstream_or_overlapping_orf_measured",
  default = "review_feature_measurement"
)]
feature_confidence[, transcript_qc_class := fcase(
  cds_exon_qc_status == "low_aggregate_cds_exon",
  "suspect_low_cds_exon",
  is_canonical_isoform == TRUE & canonical_isoform_qc_priority == TRUE,
  "canonical_qc_priority",
  is_canonical_isoform == TRUE,
  "canonical_without_qc_priority",
  cds_exon_qc_status == "pass",
  "noncanonical_exon_qc_pass",
  default = "transcript_qc_unresolved"
)]
safe_setcolorder(
  feature_confidence,
  c(
    "gene_symbol", "tx_id", "feature_id", "feature_type", "category",
    "feature_measurement_class", "transcript_qc_class",
    "is_canonical_isoform", "canonical_isoform_qc_priority",
    "cds_exon_qc_status", "cds_exon_qc_min_ratio",
    "cds_exon_qc_low_exon_ranks", "translon_source",
    "translon_sources_merged", "tx_start", "tx_end", "feature_bases",
    "cds_overlap_bases", "measurement_scope", "total_feature_raw_counts",
    "cds_reference_runs", "feature_vs_cds_mean_ratio_median",
    "feature_vs_cds_mean_ratio_q90", "runs_with_cds_signal_bump",
    "cds_signal_bump_fraction", "cds_signal_bump_supported",
    "internal_compendium_bump_class",
    "internal_frame_evidence_class", "internal_start_shape_class",
    "weighted_internal_iorf_frame_fraction",
    "weighted_internal_canonical_frame_fraction",
    "median_internal_iorf_vs_canonical_frame_delta",
    "n_internal_frame_evaluable_runs",
    "internal_iorf_frame_advantage_run_fraction",
    "internal_canonical_frame_advantage_run_fraction",
    "median_internal_start_vs_body_ratio",
    "q90_internal_start_vs_body_ratio",
    "internal_start_body_peak_run_fraction",
    "n_internal_condition_bump_review_rows",
    "n_internal_condition_bump_supported_rows",
    "n_internal_condition_bump_review_studies",
    "n_internal_condition_bump_supported_studies",
    "best_internal_condition_bump_class",
    "best_internal_condition_bump_context"
  )
)
fwrite(feature_confidence, feature_confidence_output)

browser_review_handoff <- data.table()
if (nrow(candidate_v2) > 0) {
  browser_review_handoff <- copy(
    candidate_v2[seq_len(min(nrow(candidate_v2), 80L))]
  )
  browser_review_handoff[, review_feature_id := fifelse(
    !is.na(strongest_branch_feature) & nzchar(strongest_branch_feature),
    strongest_branch_feature,
    "clean_CDS"
  )]
  browser_review_handoff[, review_reason := fcase(
    strongest_branch_feature_type == "internal_orf",
    "candidate strongest branch feature is diagnostic internal ORF",
    strongest_branch_feature_type == "overlapping_uorf",
    "candidate strongest branch feature is overlapping uORF",
    strongest_branch_feature_type == "clean_CDS",
    "candidate strongest branch feature is clean CDS",
    has_overlapping_uorf == TRUE,
    "candidate has modeled overlapping-uORF architecture",
    manual_m_source_gene == TRUE,
    "candidate has manual M-source translon support",
    default = "ranked structural/RDG candidate"
  )]

  review_features <- feature_confidence[, .(
    gene_symbol, tx_id, review_feature_id = feature_id,
    review_feature_type = feature_type,
    review_feature_category = category,
    review_feature_measurement_class = feature_measurement_class,
    review_transcript_qc_class = transcript_qc_class,
    review_translon_sources = translon_sources_merged,
    review_tx_start = tx_start,
    review_tx_end = tx_end,
    review_feature_bases = feature_bases,
    review_cds_overlap_bases = cds_overlap_bases,
    review_start_codon = start_codon,
    review_start_codon_class = start_codon_class,
    review_frame_relative_to_cds = frame_relative_to_cds,
    review_peptide_length_aa = peptide_length_aa,
    review_kozak_strength = kozak_strength
  )]
  browser_review_handoff <- merge(
    browser_review_handoff,
    review_features,
    by = c("gene_symbol", "tx_id", "review_feature_id"),
    all.x = TRUE,
    sort = FALSE
  )

  review_edge <- data.table()
  if (nrow(compact_edges) > 0) {
    review_edge <- compact_edges[
      is.finite(abs_delta_relative_translation),
      {
        idx <- which.max(abs_delta_relative_translation)
        .(
          review_branch_state = dominant_state[idx][1],
          review_branch_field = field[idx][1],
          review_branch_level = level[idx][1],
          review_branch_id = branch_id[idx][1],
          review_branch_delta_relative_translation =
            delta_relative_translation[idx][1],
          review_branch_abs_delta_relative_translation =
            abs_delta_relative_translation[idx][1],
          review_branch_log2_feature_change =
            log2_mean_fpkm_like_in_vs_out[idx][1],
          review_n_samples_in = n_samples_in[idx][1],
          review_n_samples_out = n_samples_out[idx][1]
        )
      },
      by = .(gene_symbol, tx_id, review_feature_id = feature_id)
    ]
  }
  browser_review_handoff <- merge(
    browser_review_handoff,
    review_edge,
    by = c("gene_symbol", "tx_id", "review_feature_id"),
    all.x = TRUE,
    sort = FALSE
  )

  matched_feature_best <- data.table()
  if (nrow(matched_feature_contrasts) > 0) {
    matched_feature_best <- matched_feature_contrasts[
      ,
      {
        supported <- matched_feature_support %chin% c(
          "matched_feature_shift_supported",
          "matched_feature_fpkm_shift_supported"
        )
        idx <- order(
          -as.integer(supported),
          -matched_feature_confidence_score,
          -abs_matched_relative_use_delta
        )[1]
        .(
          review_best_matched_case_condition = case_condition[idx][1],
          review_best_matched_control_conditions =
            control_conditions[idx][1],
          review_best_matched_stratum = stratum_id[idx][1],
          review_best_matched_support = matched_feature_support[idx][1],
          review_best_matched_confidence =
            matched_feature_confidence_score[idx][1],
          review_best_matched_relative_use_delta =
            matched_relative_use_delta[idx][1],
          review_best_matched_log2_fc =
            matched_fpkm_log2_fc_pseudocount1[idx][1]
        )
      },
      by = .(gene_symbol, tx_id, review_feature_id = feature_id)
    ]
  }
  browser_review_handoff <- merge(
    browser_review_handoff,
    matched_feature_best,
    by = c("gene_symbol", "tx_id", "review_feature_id"),
    all.x = TRUE,
    sort = FALSE
  )

  browser_review_handoff[, review_state_for_file := fifelse(
    !is.na(review_branch_state) & nzchar(review_branch_state),
    review_branch_state,
    strongest_branch_state
  )]
  browser_review_handoff[, rdg_png_file := file.path(
    "dominant_rdg_outputs",
    "figures",
    paste0("rdg_", gene_symbol, ".png")
  )]
  browser_review_handoff[, rdg_pdf_file := file.path(
    "dominant_rdg_outputs",
    "figures",
    paste0("rdg_", gene_symbol, ".pdf")
  )]
  browser_review_handoff[, state_rdg_png_file := file.path(
    "dominant_rdg_outputs",
    "figures_by_state",
    file_safe(review_state_for_file),
    paste0("rdg_", gene_symbol, "_", file_safe(review_state_for_file), ".png")
  )]
  browser_review_handoff[, state_rdg_pdf_file := file.path(
    "dominant_rdg_outputs",
    "figures_by_state",
    file_safe(review_state_for_file),
    paste0("rdg_", gene_symbol, "_", file_safe(review_state_for_file), ".pdf")
  )]
  browser_review_handoff[, relative_usage_matrix_file := file.path(
    "dominant_rdg_outputs",
    "relative_usage_matrices",
    "relative_usage_in_long.csv"
  )]
  browser_review_handoff[, matched_control_table_file := file.path(
    "dominant_next_model",
    "dominant_next_model_matched_control_feature_contrasts.csv"
  )]
  browser_review_handoff[, rdg_png_exists :=
                            file.exists(file.path(analysis_dir, rdg_png_file))]
  browser_review_handoff[, state_rdg_png_exists :=
                            file.exists(file.path(analysis_dir, state_rdg_png_file))]
  browser_review_handoff[, browser_review_question := fcase(
    review_feature_type == "internal_orf",
    "Check whether local p-site/start evidence supports independent iORF initiation or shared CDS elongation coverage.",
    review_feature_type == "overlapping_uorf",
    "Check whether the overlapping uORF start/stop geometry and coverage explain reduced clean-CDS allocation.",
    review_feature_type == "leader_uorf",
    "Check whether leader-uORF spacing and branch shifts support reinitiation or scanning-bypass behavior.",
    review_feature_type == "clean_CDS",
    "Check whether clean-CDS change is buffered, amplified, or a denominator artifact against upstream feature use.",
    default = "Inspect RDG topology, matched controls, and browser coverage before mechanistic interpretation."
  )]
  browser_review_handoff[, browser_review_priority_score :=
                            next_iteration_candidate_score_v4 +
                            fifelse(
                              is.finite(review_branch_abs_delta_relative_translation),
                              review_branch_abs_delta_relative_translation,
                              0
                            ) +
                            fifelse(review_feature_type == "internal_orf", 0.35, 0) +
                            fifelse(review_best_matched_support %chin% c(
                              "matched_feature_shift_supported",
                              "matched_feature_fpkm_shift_supported"
                            ), 0.35, 0)]
  setorder(browser_review_handoff, -browser_review_priority_score)
  browser_review_handoff[, browser_review_rank := seq_len(.N)]
  safe_setcolorder(
    browser_review_handoff,
    c(
      "browser_review_rank", "gene_symbol", "tx_id", "review_feature_id",
      "review_feature_type", "review_reason", "browser_review_question",
      "review_feature_measurement_class", "review_transcript_qc_class",
      "candidate_rank_v4", "novelty_class", "known_uorf_regulation",
      "manual_m_source_gene", "has_overlapping_uorf", "n_uorfs",
      "n_overlapping_uorfs", "primary_clean_cds_usage",
      "next_iteration_candidate_score_v4", "review_tx_start",
      "review_tx_end", "review_feature_bases",
      "review_cds_overlap_bases", "review_translon_sources",
      "review_start_codon", "review_start_codon_class",
      "review_kozak_strength", "review_frame_relative_to_cds",
      "review_peptide_length_aa", "review_branch_state",
      "review_branch_field", "review_branch_level",
      "review_branch_delta_relative_translation",
      "review_branch_log2_feature_change",
      "best_matched_control_context", "best_matched_control_support",
      "review_best_matched_case_condition",
      "review_best_matched_control_conditions",
      "review_best_matched_support", "rdg_png_file",
      "state_rdg_png_file", "relative_usage_matrix_file",
      "matched_control_table_file"
    )
  )
}
fwrite(browser_review_handoff, browser_review_handoff_output)

latent_reductions <- data.table()
if (nrow(state_dt) > 0) {
  state_cols <- setdiff(names(state_dt), c("Run", "sample", "dominant_state",
                                           "dominant_state_score"))
  state_cols <- state_cols[state_cols != "Baseline (control)"]
  state_cols <- state_cols[vapply(state_dt[, ..state_cols], is.numeric, logical(1))]
  if (length(state_cols) >= 2) {
    X <- as.matrix(state_dt[, ..state_cols])
    X <- X[complete.cases(X), , drop = FALSE]
    if (nrow(X) >= 10) {
      pca <- prcomp(X, center = TRUE, scale. = TRUE)
      pca_loadings <- as.data.table(pca$rotation, keep.rownames = "dominant_state")
      pca_long <- melt(
        pca_loadings,
        id.vars = "dominant_state",
        variable.name = "component",
        value.name = "loading"
      )
      pca_long[, method := "PCA"]
      pca_long[, abs_loading := abs(loading)]
      pca_long[, sparse_loading := fifelse(
        abs_loading >= quantile(abs_loading, 0.75, na.rm = TRUE),
        loading,
        0
      ), by = component]
      pca_long <- pca_long[, .(
        method, component, dominant_state, loading, sparse_loading,
        component_weight = abs_loading
      )]

      X_nonneg <- X - min(X, na.rm = TRUE)
      X_nonneg <- X_nonneg + 1e-6
      set.seed(11)
      nmf_rows <- rbindlist(lapply(2:4, function(k) {
        W <- matrix(runif(nrow(X_nonneg) * k, 0.1, 1), nrow(X_nonneg), k)
        H <- matrix(runif(k * ncol(X_nonneg), 0.1, 1), k, ncol(X_nonneg))
        for (iter in seq_len(200L)) {
          H <- H * (crossprod(W, X_nonneg) / (crossprod(W, W) %*% H + 1e-9))
          W <- W * ((X_nonneg %*% t(H)) / (W %*% (H %*% t(H)) + 1e-9))
        }
        rss <- sum((X_nonneg - W %*% H)^2)
        out <- as.data.table(H)
        setnames(out, state_cols)
        out[, component := paste0("NMF", k, "_", seq_len(.N))]
        long <- melt(out, id.vars = "component",
                     variable.name = "dominant_state",
                     value.name = "loading")
        long[, `:=`(
          method = paste0("NMF_k", k),
          sparse_loading = fifelse(
            loading >= quantile(loading, 0.75, na.rm = TRUE),
            loading,
            0
          ),
          component_weight = loading,
          reconstruction_rss = rss
        )]
        long
      }), fill = TRUE)

      state_cor <- cor(X, use = "pairwise.complete.obs")
      hc <- hclust(as.dist(1 - state_cor), method = "average")
      hier <- data.table(
        method = "hierarchical_state_correlation",
        component = paste0("cluster_order_", seq_along(hc$order)),
        dominant_state = colnames(state_cor)[hc$order],
        loading = NA_real_,
        sparse_loading = NA_real_,
        component_weight = seq_along(hc$order),
        reconstruction_rss = NA_real_
      )
      latent_reductions <- rbind(pca_long, nmf_rows, hier, fill = TRUE)
    }
  }
}
fwrite(latent_reductions, latent_output)

top_candidate_plot <- candidate_v2[is.finite(next_iteration_candidate_score)]
top_candidate_plot <- top_candidate_plot[seq_len(min(nrow(top_candidate_plot), 25L))]
if (nrow(top_candidate_plot) > 0) {
  top_candidate_plot[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
  p_rank <- ggplot(
    top_candidate_plot,
    aes(next_iteration_candidate_score, gene_symbol, fill = novelty_class)
  ) +
    geom_col(width = 0.72) +
    geom_text(
      aes(label = sprintf("CDS %.0f%%", 100 * primary_clean_cds_usage)),
      hjust = -0.05,
      size = 2.7,
      color = "#111827"
    ) +
    scale_fill_manual(values = c(
      "known uORF benchmark" = "#6b7280",
      "manual-overlap candidate" = "#b2182b",
      "under-annotated candidate" = "#ef8a62",
      "predicted-overlap candidate" = "#2166ac",
      "structure-associated candidate" = "#4d9221"
    )) +
    coord_cartesian(xlim = c(0, max(top_candidate_plot$next_iteration_candidate_score, na.rm = TRUE) * 1.16)) +
    labs(
      title = "Next-model candidate ranking",
      subtitle = paste(
        "Score combines low clean-CDS use, branch shifts, allocation effects,",
        "buffering evidence, matched-control support, and known-uORF status."
      ),
      x = "Next-iteration candidate score",
      y = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.subtitle = element_text(size = 9.3, lineheight = 1.08),
      plot.title = element_text(face = "bold")
    )
  write_plot(p_rank, "next_model_candidate_ranking", width = 9.2, height = 7.2)
}

scatter_input <- candidate_v2[
  is.finite(primary_clean_cds_usage) & is.finite(max_abs_branch_delta)
]
if (nrow(scatter_input) > 0) {
  p_scatter <- ggplot(
    scatter_input,
    aes(primary_clean_cds_usage, max_abs_branch_delta,
        color = novelty_class, size = pmax(n_uorfs, 1))
  ) +
    geom_point(alpha = 0.82) +
    geom_text(
      data = scatter_input[next_iteration_candidate_score >= quantile(next_iteration_candidate_score, 0.80, na.rm = TRUE)],
      aes(label = gene_symbol),
      hjust = -0.08,
      vjust = 0.45,
      size = 2.8,
      show.legend = FALSE
    ) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_color_manual(values = c(
      "known uORF benchmark" = "#6b7280",
      "manual-overlap candidate" = "#b2182b",
      "under-annotated candidate" = "#ef8a62",
      "predicted-overlap candidate" = "#2166ac",
      "structure-associated candidate" = "#4d9221"
    )) +
    labs(
      title = "RDG evidence space",
      subtitle = "Lower clean-CDS usage plus larger branch shifts marks stronger transcript-decision candidates.",
      x = "Primary clean-CDS relative usage",
      y = "Maximum absolute branch shift",
      color = NULL,
      size = "uORFs"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
  write_plot(p_scatter, "next_model_rdg_evidence_space", width = 9, height = 6.2)
  save_widget(
    ggplotly(p_scatter),
    file.path(figure_dir, "next_model_rdg_evidence_space.html"),
    "Next-model RDG evidence space"
  )
}

flow_plot <- flow[
  clean_cds_valid == TRUE &
    feature_type %chin% c("overlapping_uorf", "internal_orf", "clean_CDS")
]
if (nrow(flow_plot) > 0) {
  candidate_genes <- candidate_v2[seq_len(min(nrow(candidate_v2), 12L)), gene_symbol]
  selected_branch_pairs <- data.table()
  if (nrow(compact_edges) > 0 &&
      all(c("gene_symbol", "branch_id", "abs_delta_relative_translation") %chin% names(compact_edges))) {
    branch_rank <- compact_edges[
      gene_symbol %chin% candidate_genes &
        is.finite(abs_delta_relative_translation),
      .(branch_score = max(abs_delta_relative_translation, na.rm = TRUE)),
      by = .(gene_symbol, branch_id)
    ]
    branch_rank[, candidate_order := match(gene_symbol, candidate_genes)]
    branch_rank <- branch_rank[order(candidate_order, -branch_score)]
    selected_branch_pairs <- branch_rank[
      ,
      head(.SD, 3L),
      by = gene_symbol
    ][, .(gene_symbol, branch_id)]
  }
  flow_plot <- flow_plot[
    gene_symbol %chin% candidate_genes
  ]
  if (nrow(selected_branch_pairs) > 0) {
    flow_plot <- merge(
      flow_plot,
      selected_branch_pairs,
      by = c("gene_symbol", "branch_id"),
      all = FALSE,
      sort = FALSE
    )
  }
  if (nrow(flow_plot) == 0) {
    flow_plot <- flow[
      clean_cds_valid == TRUE &
        feature_type %chin% c("overlapping_uorf", "internal_orf", "clean_CDS") &
        gene_symbol %chin% candidate_genes
    ]
  }
  flow_plot[, candidate_order := match(gene_symbol, candidate_genes)]
  if (exists("branch_rank") && nrow(branch_rank) > 0) {
    flow_plot <- merge(
      flow_plot,
      branch_rank[, .(gene_symbol, branch_id, branch_score)],
      by = c("gene_symbol", "branch_id"),
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    flow_plot[, branch_score := 0]
  }
  flow_plot[!is.finite(branch_score), branch_score := 0]
  flow_plot[, compact_level := mapply(
    compact_level_label,
    field,
    level,
    USE.NAMES = FALSE
  )]
  flow_plot[, row_label := paste(
    gene_symbol,
    state_short(dominant_state),
    compact_level,
    sep = " | "
  )]
  flow_plot <- flow_plot[
    ,
    .(
      relative_translation = safe_mean(relative_translation),
      feature_capture_proxy = safe_mean(feature_capture_proxy),
      post_uorf_loss_proxy = safe_mean(post_uorf_loss_proxy),
      candidate_order = min(candidate_order, na.rm = TRUE),
      branch_score = max(branch_score, na.rm = TRUE)
    ),
    by = .(row_label, feature_id, feature_type)
  ]
  flow_plot[, feature_label := fifelse(
    feature_id == "clean_CDS",
    "clean_CDS",
    paste0(feature_id, "\n", sub("_uorf", " uORF", feature_type, fixed = TRUE))
  )]
  row_levels <- unique(flow_plot[order(candidate_order, -branch_score, row_label), row_label])
  flow_plot[, row_label := factor(row_label, levels = rev(row_levels))]
  flow_height <- min(11.5, max(5.8, 2.0 + 0.25 * uniqueN(flow_plot$row_label)))
  p_flow <- ggplot(flow_plot, aes(feature_label, row_label, fill = relative_translation)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#8b0000",
      midpoint = 0.5,
      limits = c(0, 1),
      name = "relative use"
    ) +
    labs(
      title = "Relative feature use in candidate RDGs",
      subtitle = "Blue is low, white is 0.5, dark red is high relative use inside each branch/gene transcript.",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.text.y = element_text(size = 6.8),
      panel.grid = element_blank(),
      plot.subtitle = element_text(size = 9.3, lineheight = 1.08),
      plot.title = element_text(face = "bold")
    )
  write_plot(p_flow, "next_model_candidate_relative_use_matrix", width = 11.5, height = flow_height)
  save_widget(
    ggplotly(p_flow),
    file.path(figure_dir, "next_model_candidate_relative_use_matrix.html"),
    "Next-model candidate relative-use matrix"
  )
}

if (nrow(buffering_condition) > 0 || nrow(buffering_condition_representative) > 0) {
  buffer_plot <- if (nrow(buffering_condition_representative) > 0) {
    copy(buffering_condition_representative)
  } else {
    copy(buffering_condition)
  }
  buffer_plot <- merge(
    buffer_plot,
    candidate_v2[, .(gene_symbol, tx_id, candidate_rank_v2)],
    by = c("gene_symbol", "tx_id"),
    all.x = TRUE,
    sort = FALSE
  )
  buffer_plot[, candidate_rank_plot := fifelse(is.finite(candidate_rank_v2), candidate_rank_v2, 1e9)]
  setorder(buffer_plot, candidate_rank_plot, -buffering_confidence_score,
           -abs_cds_fpkm_change_pct_stabilized)
  keep_genes <- unique(buffer_plot$gene_symbol)[seq_len(min(uniqueN(buffer_plot$gene_symbol), 14L))]
  buffer_plot <- buffer_plot[gene_symbol %chin% keep_genes]
  if (nrow(buffer_plot) > 0) {
    buffer_plot[, branch_context := mapply(
      compact_level_label,
      field,
      level,
      USE.NAMES = FALSE
    )]
    buffer_plot[, gene_state := paste(gene_symbol, state_short(dominant_state), sep = " | ")]
    gene_state_levels <- unique(buffer_plot$gene_state)
    buffer_plot[, gene_state := factor(gene_state, levels = rev(gene_state_levels))]
    p_buffer <- ggplot(buffer_plot, aes(branch_context, gene_state, fill = cds_fpkm_change_pct_stabilized)) +
      geom_tile(color = "white", linewidth = 0.25) +
      scale_fill_gradient2(
        low = "#2166ac",
        mid = "white",
        high = "#b2182b",
        midpoint = 0,
        name = "stabilized CDS change %"
      ) +
      labs(
        title = "Condition-level CDS buffering effect size (relevant RDGs)",
        subtitle = paste(
          "Values are denominator-stabilized IN vs OUT clean_CDS proxy percent change;",
          "representative rows collapse near-duplicate nested condition contexts."
        ),
        x = "Condition branch context",
        y = NULL
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 40, hjust = 1, size = 7.8),
        axis.text.y = element_text(size = 7.6),
        panel.grid = element_blank(),
        plot.subtitle = element_text(size = 9.1, lineheight = 1.08),
        plot.title = element_text(face = "bold")
      )
    buffer_height <- min(10.8, max(4.8, 2.2 + 0.28 * uniqueN(buffer_plot$gene_state)))
    write_plot(p_buffer, "next_model_cds_buffering_condition_effects", width = 12, height = buffer_height)
    save_widget(
      ggplotly(p_buffer),
      file.path(figure_dir, "next_model_cds_buffering_condition_effects.html"),
      "Condition-level CDS buffering effect sizes"
    )
  }
}

if (nrow(internal_orf_matched_bump) > 0) {
  bump_plot <- copy(internal_orf_matched_bump[
    internal_condition_bump_class %chin% c(
      "matched_case_cds_bump_supported",
      "matched_case_cds_bump_review",
      "case_cds_bump_descriptive"
    ) &
      is.finite(matched_log2_feature_vs_cds_ratio_shift)
  ])
  if (nrow(bump_plot) > 0) {
    setorder(
      bump_plot,
      -internal_condition_bump_confidence_score,
      -case_bump_fraction,
      -matched_log2_feature_vs_cds_ratio_shift
    )
    bump_plot <- bump_plot[seq_len(min(.N, 36L))]
    bump_plot[, matched_context := paste0(
      gene_symbol, " ", feature_id, " | ", case_condition, " | ",
      CELL_LINE, " / ", TISSUE
    )]
    bump_plot[, matched_context := factor(
      matched_context,
      levels = rev(unique(matched_context))
    )]
    p_internal_bump <- ggplot(
      bump_plot,
      aes(
        matched_log2_feature_vs_cds_ratio_shift,
        matched_context,
        color = internal_condition_bump_class,
        size = case_bump_fraction
      )
    ) +
      geom_vline(xintercept = 0, color = "#cbd5e1", linewidth = 0.35) +
      geom_segment(
        data = bump_plot,
        aes(
          x = 0,
          xend = matched_log2_feature_vs_cds_ratio_shift,
          y = matched_context,
          yend = matched_context
        ),
        inherit.aes = FALSE,
        linewidth = 0.35,
        color = "#94a3b8",
        show.legend = FALSE
      ) +
      geom_point(alpha = 0.92) +
      scale_color_manual(
        values = c(
          matched_case_cds_bump_supported = "#8b0000",
          matched_case_cds_bump_review = "#c2410c",
          case_cds_bump_descriptive = "#2563eb"
        ),
        name = "matched bump class"
      ) +
      scale_size_continuous(
        range = c(1.8, 5.5),
        limits = c(0, 1),
        name = "case bump\nfraction"
      ) +
      labs(
        title = "Internal-ORF matched feature-vs-CDS bump evidence",
        subtitle = paste(
          "Positive x means the internal feature/CDS density ratio is higher in case than matched control;",
          "iORFs stay diagnostic until start and frame evidence is added."
        ),
        x = "matched log2 shift of feature-vs-CDS ratio",
        y = NULL
      ) +
      theme_minimal(base_size = 10.2) +
      theme(
        axis.text.y = element_text(size = 7.2),
        panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 9, lineheight = 1.08),
        plot.title = element_text(face = "bold"),
        legend.position = "right"
      )
    bump_height <- min(11.2, max(5.5, 2.3 + 0.23 * nrow(bump_plot)))
    write_plot(
      p_internal_bump,
      "next_model_internal_orf_matched_bump_evidence",
      width = 12,
      height = bump_height
    )
    save_widget(
      ggplotly(p_internal_bump),
      file.path(
        figure_dir,
        "next_model_internal_orf_matched_bump_evidence.html"
      ),
      "Matched internal ORF feature-vs-CDS bump evidence"
    )
  }
}

if (nrow(internal_orf_initiation_evidence) > 0) {
  frame_plot <- copy(internal_orf_initiation_evidence[
    is.finite(weighted_internal_iorf_frame_fraction) &
      is.finite(weighted_internal_canonical_frame_fraction)
  ])
  if (nrow(frame_plot) > 0) {
    frame_plot[, internal_feature_label := paste(gene_symbol, feature_id)]
    p_internal_frame <- ggplot(
      frame_plot,
      aes(
        weighted_internal_canonical_frame_fraction,
        weighted_internal_iorf_frame_fraction,
        color = internal_frame_evidence_class,
        shape = internal_start_shape_class,
        size = total_internal_phase_raw_counts
      )
    ) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linewidth = 0.35,
        color = "#94a3b8",
        linetype = "dashed"
      ) +
      geom_point(alpha = 0.92) +
      geom_text(
        aes(label = internal_feature_label),
        size = 2.25,
        nudge_y = 0.018,
        check_overlap = TRUE,
        show.legend = FALSE
      ) +
      scale_color_manual(
        values = c(
          iorf_frame_advantage_review = "#8b0000",
          canonical_frame_like = "#2563eb",
          mixed_or_weak_frame_signal = "#64748b",
          low_internal_frame_counts = "#94a3b8",
          same_frame_not_discriminating = "#0f766e"
        ),
        name = "frame evidence"
      ) +
      scale_size_continuous(
        trans = "log10",
        range = c(2, 7),
        name = "internal body\nraw counts"
      ) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(
        title = "Internal-ORF frame evidence from cached positional coverage",
        subtitle = paste(
          "Above the diagonal, the feature-start frame exceeds the canonical CDS frame.",
          "Start-window shape is a review cue, not initiation proof."
        ),
        x = "weighted canonical-CDS frame fraction inside internal ORF",
        y = "weighted feature-start frame fraction inside internal ORF",
        shape = "start-window shape"
      ) +
      theme_minimal(base_size = 10.2) +
      theme(
        panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 9, lineheight = 1.08),
        plot.title = element_text(face = "bold"),
        legend.position = "right"
      )
    write_plot(
      p_internal_frame,
      "next_model_internal_orf_initiation_evidence",
      width = 10.8,
      height = 7.2
    )
    save_widget(
      ggplotly(p_internal_frame),
      file.path(
        figure_dir,
        "next_model_internal_orf_initiation_evidence.html"
      ),
      "Internal ORF cached positional frame evidence"
    )
  }
}

summary_rows <- list(
  data.table(metric = "genes_ranked", value = nrow(candidate_v2)),
  data.table(metric = "known_uorf_benchmarks_in_marker_panel", value = sum(candidate_v2$known_uorf_regulation, na.rm = TRUE)),
  data.table(metric = "manual_m_source_genes", value = sum(candidate_v2$manual_m_source_gene, na.rm = TRUE)),
  data.table(metric = "overlapping_uorf_genes", value = sum(candidate_v2$has_overlapping_uorf, na.rm = TRUE)),
  data.table(metric = "novel_candidate_rows", value = nrow(novel)),
  data.table(metric = "compact_branch_edge_rows", value = nrow(compact_edges)),
  data.table(metric = "rdg_flow_rows", value = nrow(flow)),
  data.table(metric = "cds_buffering_effect_rows_all", value = nrow(buffering_effects)),
  data.table(metric = "cds_buffering_effect_rows_relevant", value = nrow(buffering_relevant)),
  data.table(metric = "cds_buffering_effect_rows_condition_focus", value = nrow(buffering_condition)),
  data.table(metric = "cds_buffering_effect_rows_condition_representative", value = nrow(buffering_condition_representative)),
  data.table(metric = "direct_control_buffering_rows", value = sum(buffering_effects$direct_control_available_any == TRUE, na.rm = TRUE)),
  data.table(metric = "supported_control_buffering_rows", value = sum(buffering_effects$best_control_score >= 2, na.rm = TRUE)),
  data.table(metric = "matched_control_feature_contrast_rows", value = nrow(matched_feature_contrasts)),
  data.table(metric = "matched_control_feature_evaluable_rows", value = sum(matched_feature_contrasts$matched_feature_evaluable == TRUE, na.rm = TRUE)),
  data.table(metric = "matched_control_feature_low_count_rows", value = sum(matched_feature_contrasts$matched_feature_support == "low_feature_counts", na.rm = TRUE)),
  data.table(metric = "matched_control_cds_buffering_rows", value = nrow(matched_cds_buffering)),
  data.table(metric = "matched_control_cds_evaluable_rows", value = sum(matched_cds_buffering$matched_feature_evaluable == TRUE, na.rm = TRUE)),
  data.table(metric = "matched_control_cds_low_count_rows", value = sum(matched_cds_buffering$matched_cds_support == "low_clean_cds_counts", na.rm = TRUE)),
  data.table(metric = "matched_control_cds_bootstrap_rows", value = nrow(matched_cds_bootstrap)),
  data.table(metric = "matched_control_design_rows", value = nrow(matched_design_summary)),
  data.table(metric = "internal_orf_matched_bump_rows", value = nrow(internal_orf_matched_bump)),
  data.table(metric = "internal_orf_matched_bump_evaluable_rows", value = sum(internal_orf_matched_bump$internal_bump_passes_sample_gate == TRUE & internal_orf_matched_bump$internal_bump_passes_count_gate_any == TRUE, na.rm = TRUE)),
  data.table(metric = "internal_orf_matched_bump_features", value = uniqueN(
    paste(
      internal_orf_matched_bump$gene_symbol,
      internal_orf_matched_bump$tx_id,
      internal_orf_matched_bump$feature_id,
      sep = "|"
    )
  )),
  data.table(metric = "internal_orf_matched_bump_supported_rows", value = sum(
    internal_orf_matched_bump$internal_condition_bump_class ==
      "matched_case_cds_bump_supported",
    na.rm = TRUE
  )),
  data.table(metric = "internal_orf_matched_bump_review_rows", value = sum(
    internal_orf_matched_bump$internal_condition_bump_class %chin% c(
      "matched_case_cds_bump_supported",
      "matched_case_cds_bump_review"
    ),
    na.rm = TRUE
  )),
  data.table(metric = "internal_orf_initiation_evidence_rows", value = nrow(
    internal_orf_initiation_evidence
  )),
  data.table(metric = "internal_orf_feature_frame_advantage_rows", value = sum(
    internal_orf_initiation_evidence$internal_frame_evidence_class ==
      "iorf_frame_advantage_review",
    na.rm = TRUE
  )),
  data.table(metric = "internal_orf_canonical_frame_like_rows", value = sum(
    internal_orf_initiation_evidence$internal_frame_evidence_class ==
      "canonical_frame_like",
    na.rm = TRUE
  )),
  data.table(metric = "supported_matched_control_cds_buffering_rows", value = sum(
    matched_cds_buffering$matched_cds_support %chin% c(
      "matched_clean_cds_relative_shift_supported",
      "matched_clean_cds_fpkm_shift_supported"
    ),
    na.rm = TRUE
  )),
  data.table(metric = "bootstrap_supported_matched_control_cds_rows", value = sum(
    matched_cds_buffering$bootstrap_relative_use_ci_excludes_zero == TRUE |
      matched_cds_buffering$bootstrap_fpkm_log2_fc_ci_excludes_zero == TRUE,
    na.rm = TRUE
  )),
  data.table(metric = "supported_matched_control_design_rows", value = sum(
    matched_design_summary$n_supported_genes > 0,
    na.rm = TRUE
  )),
  data.table(metric = "matched_control_genes", value = uniqueN(matched_gene_summary$gene_symbol)),
  data.table(metric = "internal_orf_architecture_rows", value = nrow(internal_orf_architecture)),
  data.table(metric = "browser_review_handoff_rows", value = nrow(browser_review_handoff))
)
if (nrow(rf_performance) > 0 && "variance_explained_final" %chin% names(rf_performance)) {
  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    metric = paste0("rf_variance_explained_", file_safe(rf_performance$dominant_state)),
    value = rf_performance$variance_explained_final
  )
}
if (nrow(glmnet_stability) > 0 && "stability_fraction" %chin% names(glmnet_stability)) {
  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    metric = "stable_positive_glmnet_terms",
    value = sum(glmnet_stability$stable_positive == TRUE, na.rm = TRUE)
  )
}
summary_dt <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_dt, summary_output)

message("Saved next-model outputs in ", output_dir)
message("Saved: ", buffering_all_output)
message("Saved: ", buffering_relevant_output)
message("Saved: ", buffering_condition_output)
message("Saved: ", buffering_condition_representative_output)
message("Saved: ", matched_feature_contrast_output)
message("Saved: ", matched_cds_buffering_output)
message("Saved: ", matched_cds_bootstrap_output)
message("Saved: ", matched_gene_summary_output)
message("Saved: ", matched_design_summary_output)
message("Saved: ", internal_orf_matched_bump_output)
message("Saved: ", internal_orf_initiation_output)
message("Saved: ", internal_orf_architecture_output)
message("Saved: ", browser_review_handoff_output)
message("Saved: ", candidate_v2_output)
message("Saved: ", candidate_v3_output)
message("Saved: ", candidate_v4_output)
message("Top next-iteration candidates: ",
        paste(head(candidate_v2$gene_symbol, 8L), collapse = ", "))
