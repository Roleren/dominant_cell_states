#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (basename(here) == "dominant_cell_states" &&
        dir.exists(file.path(here, "scripts"))) {
      return(here)
    }
    candidate <- file.path(here, "dominant_cell_states")
    if (dir.exists(candidate) &&
        dir.exists(file.path(candidate, "scripts"))) {
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
input_dir <- file.path(analysis_dir, "dominant_infected_off_frame_hotspots")
output_dir <- file.path(analysis_dir, "dominant_infected_off_frame_review")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)

message("Infected off-frame hotspot browser-review queue")
message("  1. Convert recurrent hotspot calls into exact review rows")
message("  2. Add hotspot-specific run maps with metadata")
message("  3. Plot compact windows around top recurrent hotspots")

hotspot_file <- file.path(input_dir, "infected_off_frame_hotspot_recurrence.csv")
run_position_file <- file.path(input_dir, "infected_off_frame_hotspot_run_positions.csv")
profile_file <- file.path(input_dir, "infected_off_frame_hotspot_aggregate_profiles.csv")
gene_summary_file <- file.path(input_dir, "infected_off_frame_hotspot_gene_summary.csv")
run_metrics_file <- file.path(
  analysis_dir, "dominant_cds_frame_bumpiness_qc",
  "cds_frame_bumpiness_run_metrics.csv"
)

for (file in c(hotspot_file, run_position_file, profile_file,
               gene_summary_file, run_metrics_file)) {
  if (!file.exists(file)) stop("Missing required file: ", file, call. = FALSE)
}

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

collapse_text <- function(x, max_n = 8L) {
  x <- unique(safe_chr(x))
  x <- x[nzchar(x)]
  if (!length(x)) return(NA_character_)
  suffix <- if (length(x) > max_n) ";..." else ""
  paste0(paste(head(x, max_n), collapse = ";"), suffix)
}

first_existing_cols <- function(file, wanted) {
  available <- names(fread(file, nrows = 0, showProgress = FALSE))
  intersect(wanted, available)
}

hotspots <- fread(hotspot_file, showProgress = FALSE)
run_pos <- fread(run_position_file, showProgress = FALSE)
profiles <- fread(profile_file, showProgress = FALSE)
gene_summary <- fread(gene_summary_file, showProgress = FALSE)

for (dt in list(hotspots, run_pos, profiles, gene_summary)) {
  if ("gene_symbol" %in% names(dt)) dt[, gene_symbol := toupper(gene_symbol)]
}

hotspots[, hotspot_id := sprintf(
  "%s|%s|txpos:%s|codon:%s|%s",
  gene_symbol, tx_id, transcript_position, codon_index, cds_frame
)]
run_pos[, hotspot_id := sprintf(
  "%s|%s|txpos:%s|codon:%s|%s",
  gene_symbol, tx_id, transcript_position, codon_index, cds_frame
)]

class_base_score <- c(
  matched_specific_recurrent_infected_hotspot = 100,
  matched_enriched_shared_hotspot = 65,
  unmatched_recurrent_infected_hotspot = 50,
  shared_or_control_hotspot = 20,
  weak_recurrent_infected_hotspot = 18,
  single_run_infected_hotspot = 8
)
hotspots[, class_base := unname(class_base_score[hotspot_review_class])]
hotspots[is.na(class_base), class_base := 0]
hotspots[, matched_control_status := fifelse(
  group_n_runs_matched_control >= 2,
  fifelse(n_runs_with_hotspot_matched_control > 0,
          "matched_controls_with_hotspot",
          "matched_controls_without_hotspot"),
  "no_or_too_few_matched_controls"
)]
hotspots[, hotspot_review_priority_score :=
           class_base +
           8 * log1p(n_runs_with_hotspot_candidate_infected) +
           12 * pmax(0, case_minus_control_fraction) +
           2 * pmax(0, hotspot_score) +
           fifelse(group_n_runs_matched_control >= 2, 5, -8) -
           12 * pmax(0, fraction_runs_with_hotspot_matched_control)]
hotspots[, review_tier := fifelse(
  hotspot_review_priority_score >= 80 &
    hotspot_review_class %chin% c(
      "matched_specific_recurrent_infected_hotspot",
      "matched_enriched_shared_hotspot"
    ),
  "tier1_matched_review",
  fifelse(
    hotspot_review_priority_score >= 55,
    "tier2_recurrent_review",
    fifelse(
      hotspot_review_priority_score >= 35,
      "tier3_context_review",
      "tier4_low_priority"
    )
  )
)]
hotspots[, review_question := fifelse(
  hotspot_review_class == "matched_specific_recurrent_infected_hotspot",
  "Does this exact off-frame CDS hotspot appear only in infected samples?",
  fifelse(
    hotspot_review_class == "matched_enriched_shared_hotspot",
    "This hotspot is enriched in infected samples but also appears in controls; is infection amplifying a baseline peak?",
    fifelse(
      hotspot_review_class == "unmatched_recurrent_infected_hotspot",
      "Recurrent in infected candidates but lacks matched controls; find controls or validate in browser first.",
      fifelse(
        hotspot_review_class == "shared_or_control_hotspot",
        "Hotspot is shared with or stronger in controls; likely baseline or technical unless browser shows otherwise.",
        "Single-run hotspot; only review if gene/context is otherwise compelling."
      )
    )
  )
)]

gene_context <- gene_summary[
  ,
  .(gene_symbol, candidate_score, review_tier_gene = review_tier,
    total_case_n_off_frame_peak, top_context_label)
]
queue <- merge(hotspots, gene_context, by = "gene_symbol", all.x = TRUE,
               sort = FALSE)
setorder(queue, -hotspot_review_priority_score, -hotspot_score,
         gene_symbol, codon_index)
queue[, browser_review_rank := seq_len(.N)]
queue[, browser_review_label := paste0(
  browser_review_rank, ". ", gene_symbol, " codon ", codon_index,
  " ", cds_frame
)]
queue[, suggested_browser_action := fifelse(
  hotspot_review_class %chin% c("matched_specific_recurrent_infected_hotspot",
                                "matched_enriched_shared_hotspot"),
  "Open candidate infected runs and matched controls at the hotspot; compare frame-colored peaks.",
  fifelse(
    hotspot_review_class == "unmatched_recurrent_infected_hotspot",
    "Open recurrent infected runs first; then search the same study for usable controls or noncandidate infected runs.",
    "Low priority unless the same hotspot is visually striking or overlaps known biology."
  )
)]

queue_out <- queue[
  ,
  .(
    browser_review_rank,
    browser_review_label,
    review_tier,
    hotspot_review_priority_score,
    gene_symbol,
    tx_id,
    codon_index,
    cds_frame,
    cds_relative_position,
    cds_percent,
    transcript_position,
    hotspot_review_class,
    matched_control_status,
    n_runs_with_hotspot_candidate_infected,
    group_n_runs_candidate_infected,
    n_runs_with_hotspot_matched_control,
    group_n_runs_matched_control,
    n_runs_with_hotspot_same_context_infected_noncandidate,
    group_n_runs_same_context_infected_noncandidate,
    fraction_runs_with_hotspot_candidate_infected,
    fraction_runs_with_hotspot_matched_control,
    fraction_runs_with_hotspot_same_context_infected_noncandidate,
    case_minus_control_fraction,
    case_minus_same_context_infected_fraction,
    hotspot_score,
    fisher_p_matched_control,
    q_matched_control,
    median_position_count_share_candidate_infected,
    max_position_count_share_candidate_infected,
    best_rank_in_run_candidate_infected,
    run_examples_candidate_infected,
    run_examples_matched_control,
    run_examples_same_context_infected_noncandidate,
    contexts_candidate_infected,
    candidate_score,
    review_tier_gene,
    total_case_n_off_frame_peak,
    top_context_label,
    review_question,
    suggested_browser_action
  )
]
fwrite(queue_out, file.path(output_dir,
                            "infected_off_frame_hotspot_review_queue.csv"))

metadata_cols <- c(
  "Run", "BioProject", "Study_Pubmed_id", "AUTHOR", "sample_source",
  "sample_title", "LIBRARYTYPE", "REPLICATE", "CONDITION", "INHIBITOR",
  "TIMEPOINT", "TISSUE", "CELL_LINE", "GENE", "FRACTION", "Cancer_type",
  "Cell_model", "Cell_type", "Organ_system", "Sex", "Life_stage", "study",
  "Frame_usage_0", "Frame_usage_1", "Frame_usage_2",
  "metadata_canonical_frame_fraction", "metadata_dominant_frame",
  "metadata_possible_misshift", "global_canonical_frame_fraction",
  "global_possible_misshift"
)
metadata_cols <- first_existing_cols(run_metrics_file, metadata_cols)
message("Reading compact run metadata for hotspot run map.")
run_metadata <- fread(run_metrics_file, select = metadata_cols,
                      showProgress = FALSE)
run_metadata <- unique(run_metadata, by = "Run")

run_map <- merge(
  run_pos[hotspot_id %chin% queue$hotspot_id],
  queue_out[
    ,
    .(hotspot_id = sprintf("%s|%s|txpos:%s|codon:%s|%s",
                           gene_symbol, tx_id, transcript_position,
                           codon_index, cds_frame),
      browser_review_rank, review_tier, hotspot_review_priority_score,
      hotspot_review_class, matched_control_status)
  ],
  by = "hotspot_id",
  all.x = TRUE,
  sort = FALSE
)
run_map <- merge(run_map, run_metadata, by = "Run", all.x = TRUE,
                 sort = FALSE)
setorder(run_map, browser_review_rank, hotspot_run_group,
         hotspot_rank_in_run, Run)
fwrite(run_map, file.path(output_dir,
                          "infected_off_frame_hotspot_review_run_map.csv"))

queue_summary <- queue_out[
  ,
  .(
    n_hotspots = .N,
    n_genes = uniqueN(gene_symbol),
    n_candidate_run_hits =
      sum(n_runs_with_hotspot_candidate_infected, na.rm = TRUE),
    n_control_run_hits =
      sum(n_runs_with_hotspot_matched_control, na.rm = TRUE),
    top_examples = collapse_text(browser_review_label, max_n = 8L)
  ),
  by = .(review_tier, hotspot_review_class, matched_control_status)
]
setorder(queue_summary, review_tier, hotspot_review_class)
fwrite(queue_summary, file.path(output_dir,
                                "infected_off_frame_hotspot_review_summary.csv"))

summary_metrics <- data.table(
  metric = c(
    "n_review_hotspots",
    "n_review_genes",
    "n_tier1_matched_review",
    "n_tier2_recurrent_review",
    "n_tier3_context_review",
    "n_review_run_map_rows"
  ),
  value = as.character(c(
    nrow(queue_out),
    uniqueN(queue_out$gene_symbol),
    queue_out[review_tier == "tier1_matched_review", .N],
    queue_out[review_tier == "tier2_recurrent_review", .N],
    queue_out[review_tier == "tier3_context_review", .N],
    nrow(run_map)
  ))
)
fwrite(summary_metrics, file.path(output_dir,
                                  "infected_off_frame_hotspot_review_metrics.csv"))

top_bar_dt <- head(queue_out, 30L)
if (nrow(top_bar_dt) > 0) {
  top_bar_dt[, browser_review_label := factor(
    browser_review_label, levels = rev(browser_review_label)
  )]
  p <- ggplot(
    top_bar_dt,
    aes(
      x = browser_review_label,
      y = hotspot_review_priority_score,
      fill = hotspot_review_class,
      text = paste0(
        as.character(browser_review_label),
        "<br>priority=", round(hotspot_review_priority_score, 2),
        "<br>class=", hotspot_review_class,
        "<br>case hits=", n_runs_with_hotspot_candidate_infected,
        " / ", group_n_runs_candidate_infected,
        "<br>control hits=", n_runs_with_hotspot_matched_control,
        " / ", group_n_runs_matched_control,
        "<br>question=", review_question
      )
    )
  ) +
    geom_col(width = 0.72) +
    coord_flip() +
    scale_fill_manual(
      values = c(
        matched_specific_recurrent_infected_hotspot = "#8b0000",
        matched_enriched_shared_hotspot = "#d95f02",
        unmatched_recurrent_infected_hotspot = "#3182bd",
        shared_or_control_hotspot = "#9e9e9e",
        weak_recurrent_infected_hotspot = "#bdbdbd",
        single_run_infected_hotspot = "#d9d9d9"
      ),
      drop = FALSE
    ) +
    labs(
      title = "Off-Frame CDS Hotspot Browser-Review Queue",
      subtitle = "Priority favors matched specificity first, then recurrent unmatched hotspots.",
      x = NULL,
      y = "Review priority score",
      fill = "Hotspot class"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.major.y = element_blank()
    )
  ggsave(file.path(figure_dir,
                   "infected_off_frame_hotspot_review_priority.png"),
         p, width = 11, height = 8, dpi = 150)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p,
      file.path(figure_dir,
                "infected_off_frame_hotspot_review_priority.html"),
      title = "infected_off_frame_hotspot_review_priority"
    )
  }
}

top_windows <- head(
  queue_out[hotspot_review_class != "single_run_infected_hotspot"],
  16L
)
if (nrow(top_windows) == 0) {
  top_windows <- head(queue_out, 16L)
}
window_radius_codons <- 30L
window_dt <- rbindlist(lapply(seq_len(nrow(top_windows)), function(i) {
  row <- top_windows[i]
  prof <- profiles[
    gene_symbol == row$gene_symbol &
      tx_id == row$tx_id &
      abs(codon_index - row$codon_index) <= window_radius_codons
  ]
  if (!nrow(prof)) return(data.table())
  prof[, browser_review_rank := row$browser_review_rank]
  prof[, hotspot_label := paste0(
    row$browser_review_rank, ". ", row$gene_symbol,
    " codon ", row$codon_index, " ", row$cds_frame
  )]
  prof[, hotspot_codon_index := row$codon_index]
  prof[, codon_offset := codon_index - row$codon_index]
  prof[, hotspot_review_class := row$hotspot_review_class]
  prof[, hotspot_review_priority_score := row$hotspot_review_priority_score]
  prof
}), fill = TRUE)

if (nrow(window_dt) > 0) {
  window_dt[, hotspot_label := factor(
    hotspot_label,
    levels = top_windows$browser_review_label
  )]
  window_dt[, profile_group := factor(
    profile_group,
    levels = c("candidate_infected", "matched_control",
               "same_context_infected_noncandidate")
  )]
  p2 <- ggplot(
    window_dt,
    aes(
      x = codon_offset,
      y = normalized_density,
      color = profile_group,
      group = profile_group,
      text = paste0(
        as.character(hotspot_label),
        "<br>group=", profile_group,
        "<br>codon offset=", codon_offset,
        "<br>codon=", codon_index,
        "<br>frame=", cds_frame,
        "<br>density=", round(normalized_density, 3),
        "<br>n runs=", n_runs
      )
    )
  ) +
    geom_vline(xintercept = 0, color = "#222222", linewidth = 0.35) +
    geom_line(linewidth = 0.35, alpha = 0.85, na.rm = TRUE) +
    facet_wrap(~ hotspot_label, scales = "free_y", ncol = 2) +
    scale_color_manual(
      values = c(
        candidate_infected = "#8b0000",
        matched_control = "#2166ac",
        same_context_infected_noncandidate = "#7f7f7f"
      ),
      drop = FALSE
    ) +
    labs(
      title = "Zoomed Windows Around Prioritized Off-Frame Hotspots",
      subtitle = paste0(
        "Vertical line marks the candidate hotspot; window = +/- ",
        window_radius_codons, " codons."
      ),
      x = "Codon offset from hotspot",
      y = "Normalized density per 1000 CDS reads",
      color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  ggsave(file.path(figure_dir,
                   "infected_off_frame_hotspot_review_windows.png"),
         p2, width = 11, height = 12, dpi = 150)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p2,
      file.path(figure_dir,
                "infected_off_frame_hotspot_review_windows.html"),
      title = "infected_off_frame_hotspot_review_windows"
    )
  }
}

message("Saved hotspot review queue outputs to: ", output_dir)
message("Top review rows:")
print(head(queue_out[
  ,
  .(browser_review_rank, gene_symbol, codon_index, cds_frame, review_tier,
    hotspot_review_class, hotspot_review_priority_score,
    n_runs_with_hotspot_candidate_infected,
    n_runs_with_hotspot_matched_control, group_n_runs_matched_control)
], 20L))
