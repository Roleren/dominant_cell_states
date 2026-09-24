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
})

analysis_dir <- if (dir.exists("dominant_cell_states")) "dominant_cell_states" else "."
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
labels_file <- file.path(
  analysis_dir,
  "dominant_state_score_semantics",
  "dominant_state_score_semantics_labels.csv"
)
feature_file <- file.path(analysis_dir, "dominant_uorf_feature_expression.csv")
state_file <- file.path(analysis_dir, "human_dominant_cell_states_clean_cds.csv")
selected_fields_file <- file.path(
  analysis_dir,
  "dominant_state_broad_metadata_selected_fields.csv"
)
output_dir <- file.path(analysis_dir, "dominant_state_score_semantics_design_aware")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- c(
  file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R"),
  file.path("scripts", "dominant_htmlwidgets.R"),
  "dominant_htmlwidgets.R"
)
htmlwidget_helper <- htmlwidget_helper[file.exists(htmlwidget_helper)][1]
if (!is.na(htmlwidget_helper)) source(htmlwidget_helper)

required_files <- c(metadata_file, labels_file, feature_file, state_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Design-aware score-semantics step is missing files: ",
       paste(missing_files, collapse = "; "), call. = FALSE)
}

oxphos_arm_semantics <- "oxphos_arm_top_decile"
min_term_samples <- 8L
min_study_branch_n <- 2L
min_clean_cds_raw_counts <- 30L
design_field_exclusions <- c(
  "study", "AUTHOR", "LibrarySelection", "FRACTION", "Platform",
  "LIBRARYTYPE", "mode", "BATCH"
)
fallback_metadata_fields <- c(
  "CELL_LINE", "TISSUE", "Cell_type", "Cancer_type", "Organ_system",
  "Sex", "Life_stage", "CONDITION", "GENE", "INHIBITOR", "Cell_model"
)

clean_level <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "NA", "N/A", "na", "n/a", "NULL", "null",
               "None", "none", "missing", "Missing", "unknown",
               "Unknown")] <- NA_character_
  x
}

file_safe <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  substr(x, 1L, 160L)
}

safe_sum <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  sum(x, na.rm = TRUE)
}

safe_weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  weighted.mean(x[ok], w[ok])
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || uniqueN(x[ok]) < 2 || uniqueN(y[ok]) < 2) {
    return(NA_real_)
  }
  cor(x[ok], y[ok], method = "spearman")
}

read_metadata_fields <- function(metadata) {
  fields <- character()
  if (file.exists(selected_fields_file)) {
    profile <- fread(selected_fields_file, showProgress = FALSE)
    if (all(c("field", "used_for_broad_model") %chin% names(profile))) {
      fields <- profile[used_for_broad_model == TRUE, field]
    }
  }
  if (length(fields) == 0) fields <- fallback_metadata_fields
  fields <- intersect(fields, names(metadata))
  setdiff(unique(fields), design_field_exclusions)
}

save_plot_pair <- function(plot, name, width, height, html = TRUE) {
  ggsave(
    file.path(figure_dir, paste0(name, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 220,
    bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(name, ".pdf")),
    plot,
    width = width,
    height = height,
    bg = "white"
  )
  if (html && exists("dominant_save_ggplotly", mode = "function") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      plot,
      file.path(figure_dir, paste0(name, ".html")),
      title = name
    )
  }
}

fixed_effect_summary <- function(strata, by_cols, beta_col, variance_col,
                                 sign_col = beta_col) {
  strata[
    is.finite(get(beta_col)) & is.finite(get(variance_col)) &
      get(variance_col) > 0,
    {
      beta <- get(beta_col)
      variance <- get(variance_col)
      weight <- 1 / variance
      beta_fixed <- sum(beta * weight) / sum(weight)
      se_fixed <- sqrt(1 / sum(weight))
      z_fixed <- beta_fixed / se_fixed
      top_weight <- max(weight) / sum(weight)
      .(
        n_studies = uniqueN(study_id),
        n_study_contrasts = .N,
        total_in_samples = sum(n_in, na.rm = TRUE),
        total_out_samples = sum(n_out, na.rm = TRUE),
        fixed_effect = beta_fixed,
        fixed_effect_se = se_fixed,
        fixed_effect_z = z_fixed,
        fixed_effect_greater_p = pnorm(z_fixed, lower.tail = FALSE),
        median_effect = median(beta, na.rm = TRUE),
        weighted_effect = safe_weighted_mean(beta, weight),
        fraction_positive_studies = mean(get(sign_col) > 0, na.rm = TRUE),
        top_study_weight_fraction = top_weight,
        max_abs_study_effect = max(abs(beta), na.rm = TRUE)
      )
    },
    by = by_cols
  ][
    ,
    fixed_effect_greater_q := p.adjust(fixed_effect_greater_p, "BH")
  ][
    order(fixed_effect_greater_q, -fixed_effect, -n_studies)
  ][]
}

metadata_study_contrasts_one_label <- function(context, label, label_runs,
                                               metadata_fields) {
  dt <- copy(context)
  dt[, branch_in := Run %chin% label_runs]
  rbindlist(lapply(metadata_fields, function(field) {
    levels_dt <- dt[
      !is.na(get(field)),
      .(term_n_global = .N),
      by = .(level = get(field))
    ][term_n_global >= min_term_samples]
    if (nrow(levels_dt) == 0) return(NULL)

    valid_dt <- dt[
      !is.na(get(field)) & get(field) %chin% levels_dt$level,
      .(Run, study_id, branch_in, level = get(field))
    ]
    study_total <- valid_dt[
      ,
      .(
        study_n = uniqueN(Run),
        label_n = uniqueN(Run[branch_in == TRUE])
      ),
      by = study_id
    ]
    term_counts <- valid_dt[
      ,
      .(
        term_n = uniqueN(Run),
        label_term_n = uniqueN(Run[branch_in == TRUE])
      ),
      by = .(study_id, level)
    ]
    term_counts <- merge(
      term_counts,
      study_total,
      by = "study_id",
      all.x = TRUE,
      sort = FALSE
    )
    term_counts[, `:=`(
      other_n = study_n - term_n,
      label_other_n = label_n - label_term_n
    )]
    term_counts <- term_counts[
      term_n > 0 & other_n > 0 &
        label_n > 0 & label_n < study_n
    ]
    if (nrow(term_counts) == 0) return(NULL)

    term_counts[, `:=`(
      no_label_term_n = term_n - label_term_n,
      no_label_other_n = other_n - label_other_n
    )]
    term_counts[, log_or_ln := log(
      ((label_term_n + 0.5) * (no_label_other_n + 0.5)) /
        ((no_label_term_n + 0.5) * (label_other_n + 0.5))
    )]
    term_counts[, log_or_variance := 1 / (label_term_n + 0.5) +
                  1 / (no_label_term_n + 0.5) +
                  1 / (label_other_n + 0.5) +
                  1 / (no_label_other_n + 0.5)]
    term_counts[, `:=`(
      field = field,
      label_semantics = label$label_semantics,
      label_name = label$label_name,
      dominant_state = label$dominant_state,
      score_measure = label$score_measure,
      n_in = label_n,
      n_out = study_n - label_n
    )]
    term_counts[
      ,
      .(
        label_semantics, label_name, dominant_state, score_measure,
        field, level, study_id, term_n, other_n, label_term_n,
        no_label_term_n, label_other_n, no_label_other_n, n_in, n_out,
        log_or_ln, log_or_variance
      )
    ]
  }), fill = TRUE)
}

metadata_study_contrasts <- function(context, label_keys, labels,
                                     metadata_fields) {
  rbindlist(lapply(seq_len(nrow(label_keys)), function(i) {
    label <- label_keys[i]
    label_runs <- labels[
      label_semantics == label$label_semantics &
        label_name == label$label_name &
        score_measure == label$score_measure,
      unique(Run)
    ]
    metadata_study_contrasts_one_label(
      context,
      label,
      label_runs,
      metadata_fields
    )
  }), fill = TRUE)
}

aggregate_feature_scores_one_label <- function(feature_context, label,
                                               label_runs, mixed_studies) {
  fx <- copy(feature_context[study_id %chin% mixed_studies])
  if (nrow(fx) == 0) return(data.table())
  fx[, branch_group := fifelse(Run %chin% label_runs, "in", "out")]
  scores <- fx[
    ,
    .(
      n_samples = uniqueN(Run),
      sum_raw_counts = safe_sum(raw_counts),
      mean_fpkm_like = mean(fpkm_like, na.rm = TRUE)
    ),
    by = .(
      study_id, branch_group, gene_symbol, tx_id, feature_id, feature_type
    )
  ]
  scores[, `:=`(
    label_semantics = label$label_semantics,
    label_name = label$label_name,
    dominant_state = label$dominant_state,
    score_measure = label$score_measure
  )]
  scores[
    ,
    allocation_feature :=
      feature_type %chin% c("leader_uorf", "overlapping_uorf", "clean_CDS")
  ][
    ,
    total_mean_fpkm_like := sum(mean_fpkm_like, na.rm = TRUE),
    by = .(label_semantics, label_name, score_measure, study_id, branch_group,
           gene_symbol)
  ][
    ,
    allocation_total_mean_fpkm_like :=
      sum(mean_fpkm_like[allocation_feature], na.rm = TRUE),
    by = .(label_semantics, label_name, score_measure, study_id, branch_group,
           gene_symbol)
  ][
    ,
    clean_cds_raw_counts := sum_raw_counts[feature_id == "clean_CDS"][1],
    by = .(label_semantics, label_name, score_measure, study_id, branch_group,
           gene_symbol)
  ][
    ,
    passes_clean_cds_count_cutoff :=
      is.finite(clean_cds_raw_counts) &
        clean_cds_raw_counts >= min_clean_cds_raw_counts
  ][
    ,
    relative_translation := fifelse(
      passes_clean_cds_count_cutoff &
        is.finite(allocation_total_mean_fpkm_like) &
        allocation_total_mean_fpkm_like > 0,
      mean_fpkm_like / allocation_total_mean_fpkm_like,
      NA_real_
    )
  ][]
}

contrast_feature_scores <- function(scores) {
  key_cols <- c(
    "label_semantics", "label_name", "dominant_state", "score_measure",
    "study_id", "gene_symbol", "tx_id", "feature_id", "feature_type"
  )
  in_scores <- scores[branch_group == "in"]
  out_scores <- scores[branch_group == "out"]
  out <- merge(
    in_scores[
      ,
      c(
        key_cols, "n_samples", "sum_raw_counts", "mean_fpkm_like",
        "clean_cds_raw_counts", "passes_clean_cds_count_cutoff",
        "relative_translation"
      ),
      with = FALSE
    ],
    out_scores[
      ,
      c(
        key_cols, "n_samples", "sum_raw_counts", "mean_fpkm_like",
        "clean_cds_raw_counts", "passes_clean_cds_count_cutoff",
        "relative_translation"
      ),
      with = FALSE
    ],
    by = key_cols,
    suffixes = c("_in", "_out"),
    all = FALSE,
    sort = FALSE
  )
  out[, `:=`(
    comparison_passes_clean_cds_count_cutoff =
      passes_clean_cds_count_cutoff_in & passes_clean_cds_count_cutoff_out,
    delta_relative_translation =
      relative_translation_in - relative_translation_out,
    log2_mean_fpkm_like_in_vs_out =
      log2((mean_fpkm_like_in + 1) / (mean_fpkm_like_out + 1))
  )]
  out[
    comparison_passes_clean_cds_count_cutoff == FALSE,
    `:=`(
      delta_relative_translation = NA_real_,
      log2_mean_fpkm_like_in_vs_out = NA_real_
    )
  ]
  out[]
}

rdg_study_scores <- function(feature_context, label_keys, labels,
                             mixed_study_counts) {
  rbindlist(lapply(seq_len(nrow(label_keys)), function(i) {
    label <- label_keys[i]
    label_runs <- labels[
      label_semantics == label$label_semantics &
        label_name == label$label_name &
        score_measure == label$score_measure,
      unique(Run)
    ]
    mixed_studies <- mixed_study_counts[
      label_semantics == label$label_semantics &
        label_name == label$label_name &
        score_measure == label$score_measure &
        in_n >= min_study_branch_n & out_n >= min_study_branch_n,
      study_id
    ]
    aggregate_feature_scores_one_label(
      feature_context,
      label,
      label_runs,
      mixed_studies
    )
  }), fill = TRUE)
}

summarize_clean_cds_studies <- function(clean_cds) {
  d <- clean_cds[
    is.finite(delta_relative_translation) &
      n_samples_in >= min_study_branch_n &
      n_samples_out >= min_study_branch_n
  ]
  d[, `:=`(
    n_in = n_samples_in,
    n_out = n_samples_out,
    min_group_n = pmin(n_samples_in, n_samples_out)
  )]
  d[, delta_variance_proxy := 1 / pmax(min_group_n, 1)]
  fixed_effect_summary(
    d,
    c("label_semantics", "label_name", "dominant_state",
      "score_measure", "gene_symbol", "tx_id"),
    "delta_relative_translation",
    "delta_variance_proxy"
  )[
    ,
    .(
      label_semantics, label_name, dominant_state, score_measure,
      gene_symbol, tx_id, n_studies, n_study_contrasts,
      total_in_samples, total_out_samples,
      weighted_delta_relative_translation = fixed_effect,
      weighted_delta_se_proxy = fixed_effect_se,
      weighted_delta_z_proxy = fixed_effect_z,
      weighted_delta_greater_p_proxy = fixed_effect_greater_p,
      weighted_delta_greater_q_proxy = fixed_effect_greater_q,
      median_delta_relative_translation = median_effect,
      fraction_positive_studies,
      top_study_weight_fraction,
      max_abs_study_delta = max_abs_study_effect
    )
  ][order(-n_studies, -abs(weighted_delta_relative_translation))][]
}

labels <- fread(labels_file, showProgress = FALSE)
labels <- unique(labels[, .(
  Run, label_semantics, label_name, dominant_state, score_measure,
  score_question
)])
state_runs <- fread(state_file, select = "Run", showProgress = FALSE)
metadata <- fread(metadata_file, showProgress = FALSE)
if (metadata[Run %chin% state_runs$Run, anyDuplicated(Run)] > 0) {
  stop("Metadata has duplicate Runs matching the dominant-state table.",
       call. = FALSE)
}
metadata_fields <- read_metadata_fields(metadata)

context <- merge(
  state_runs[, .(Run, input_order = .I)],
  metadata,
  by = "Run",
  all.x = TRUE,
  sort = FALSE
)
setorder(context, input_order)
context[, input_order := NULL]
context[, study_clean := clean_level(study)]
context[is.na(study_clean), study_clean := "MISSING_STUDY"]
context[, study_id := study_clean]
for (field in metadata_fields) {
  context[, (field) := clean_level(get(field))]
}

label_keys <- unique(labels[, .(
  label_semantics, label_name, dominant_state, score_measure, score_question
)])
label_context <- merge(
  labels,
  context[, .(Run, study_id, study, AUTHOR)],
  by = "Run",
  all.x = TRUE,
  sort = FALSE
)
study_totals <- context[, .(study_n = uniqueN(Run)), by = study_id]
label_study_counts <- label_context[
  ,
  .(
    in_n = uniqueN(Run),
    example_study = study[1],
    example_author = AUTHOR[1]
  ),
  by = .(label_semantics, label_name, dominant_state, score_measure, study_id)
]
label_study_counts <- merge(
  label_study_counts,
  study_totals,
  by = "study_id",
  all.x = TRUE,
  sort = FALSE
)
label_study_counts[, out_n := study_n - in_n]
label_study_counts[, mixed_study :=
                     in_n >= min_study_branch_n &
                       out_n >= min_study_branch_n]

label_design_summary <- label_study_counts[
  ,
  {
    idx <- which.max(in_n)
    label_n <- sum(in_n)
    .(
      label_n = label_n,
      label_studies = uniqueN(study_id),
      mixed_studies = sum(mixed_study),
      mixed_in_samples = sum(in_n[mixed_study]),
      mixed_out_samples = sum(out_n[mixed_study]),
      top_study = study_id[idx][1],
      top_study_label_n = in_n[idx][1],
      top_study_fraction = in_n[idx][1] / label_n,
      studies_with_at_least_5_labels = sum(in_n >= 5)
    )
  },
  by = .(label_semantics, label_name, dominant_state, score_measure)
][order(label_semantics, -label_n, label_name)][]

metadata_strata <- metadata_study_contrasts(
  context,
  label_keys,
  labels,
  metadata_fields
)
metadata_summary <- fixed_effect_summary(
  metadata_strata,
  c("label_semantics", "label_name", "dominant_state", "score_measure",
    "field", "level"),
  "log_or_ln",
  "log_or_variance"
)
metadata_summary[, `:=`(
  fixed_log2_odds = fixed_effect / log(2),
  fixed_log2_odds_se = fixed_effect_se / log(2),
  median_study_log2_odds = median_effect / log(2),
  max_abs_study_log2_odds = max_abs_study_effect / log(2)
)]
metadata_summary <- metadata_summary[
  ,
  .(
    label_semantics, label_name, dominant_state, score_measure,
    field, level, n_studies, n_study_contrasts, total_in_samples,
    total_out_samples, fixed_log2_odds, fixed_log2_odds_se,
    fixed_effect_z, fixed_effect_greater_p, fixed_effect_greater_q,
    median_study_log2_odds, fraction_positive_studies,
    top_study_weight_fraction, max_abs_study_log2_odds
  )
][order(fixed_effect_greater_q, -fixed_log2_odds, -n_studies)][]

feature_dt <- fread(feature_file, showProgress = FALSE)
required_feature_cols <- c(
  "Run", "gene_symbol", "tx_id", "feature_id", "feature_type",
  "raw_counts", "fpkm_like"
)
if (!all(required_feature_cols %chin% names(feature_dt))) {
  stop("Feature table lacks design-aware RDG columns: ",
       paste(setdiff(required_feature_cols, names(feature_dt)), collapse = ", "),
       call. = FALSE)
}
feature_dt <- feature_dt[
  feature_type %chin% c("leader_uorf", "overlapping_uorf", "internal_orf",
                        "clean_CDS") &
    is.finite(fpkm_like)
]
feature_dt[, raw_counts := as.numeric(raw_counts)]
feature_context <- merge(
  feature_dt,
  context[, .(Run, study_id)],
  by = "Run",
  all.x = TRUE,
  sort = FALSE
)
feature_scores <- rdg_study_scores(
  feature_context,
  label_keys,
  labels,
  label_study_counts
)
feature_contrasts <- contrast_feature_scores(feature_scores)
clean_cds_study_contrasts <- feature_contrasts[feature_id == "clean_CDS"]
clean_cds_summary <- summarize_clean_cds_studies(clean_cds_study_contrasts)

summary_metrics <- data.table(
  metric = c(
    "metadata_fields_design_aware",
    "score_semantics_labels",
    "label_study_rows",
    "mixed_label_study_rows",
    "metadata_study_term_contrast_rows",
    "metadata_fixed_effect_rows",
    "rdg_study_feature_contrast_rows",
    "rdg_study_clean_cds_contrast_rows",
    "rdg_study_clean_cds_valid_rows",
    "rdg_clean_cds_fixed_effect_rows"
  ),
  value = c(
    length(metadata_fields),
    nrow(label_design_summary),
    nrow(label_study_counts),
    label_study_counts[mixed_study == TRUE, .N],
    nrow(metadata_strata),
    nrow(metadata_summary),
    nrow(feature_contrasts),
    nrow(clean_cds_study_contrasts),
    clean_cds_study_contrasts[is.finite(delta_relative_translation), .N],
    nrow(clean_cds_summary)
  )
)

fwrite(
  label_study_counts,
  file.path(output_dir, "dominant_state_score_semantics_label_study_counts.csv")
)
fwrite(
  label_design_summary,
  file.path(output_dir, "dominant_state_score_semantics_label_design_summary.csv")
)
fwrite(
  metadata_strata,
  file.path(output_dir, "dominant_state_score_semantics_metadata_study_contrasts.csv")
)
fwrite(
  metadata_summary,
  file.path(output_dir, "dominant_state_score_semantics_metadata_fixed_effects.csv")
)
fwrite(
  feature_contrasts,
  file.path(output_dir, "dominant_state_score_semantics_rdg_study_feature_contrasts.csv")
)
fwrite(
  clean_cds_study_contrasts,
  file.path(output_dir, "dominant_state_score_semantics_rdg_study_clean_cds_contrasts.csv")
)
fwrite(
  clean_cds_summary,
  file.path(output_dir, "dominant_state_score_semantics_rdg_clean_cds_fixed_effects.csv")
)
fwrite(
  summary_metrics,
  file.path(output_dir, "dominant_state_score_semantics_design_summary_metrics.csv")
)

diversity_plot <- copy(label_design_summary)
diversity_plot[, label_text := paste0(label_name, "\n", label_n, " samples")]
p_diversity <- ggplot(
  diversity_plot,
  aes(mixed_studies, top_study_fraction, color = label_semantics)
) +
  geom_point(aes(size = label_n), alpha = 0.78) +
  geom_text(
    aes(label = label_text),
    check_overlap = TRUE,
    hjust = -0.08,
    size = 2.8,
    show.legend = FALSE
  ) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  scale_color_manual(
    values = c(
      "raw_signed_winner" = "#345B8C",
      "within_state_percentile_winner" = "#B24A3C",
      "oxphos_arm_top_decile" = "#D89A36"
    )
  ) +
  scale_size_continuous(range = c(2.2, 7.8)) +
  labs(
    title = "Score-semantics labels differ in design diversity",
    subtitle = paste(
      "Mixed studies have at least", min_study_branch_n,
      "IN and OUT samples for within-study RDG contrasts."
    ),
    x = "Mixed studies",
    y = "Top study share of label samples",
    color = NULL,
    size = "Label samples"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot"
  )
save_plot_pair(
  p_diversity,
  "dominant_state_score_semantics_design_diversity",
  width = 12,
  height = 7.6
)

ox_metadata <- metadata_summary[
  label_semantics == oxphos_arm_semantics &
    n_studies >= 2 &
    is.finite(fixed_log2_odds) &
    is.finite(fixed_effect_greater_q)
][
  order(fixed_effect_greater_q, -fixed_log2_odds),
  head(.SD, 10L),
  by = label_name
]
if (nrow(ox_metadata) > 0) {
  ox_metadata[, term := paste(field, level, sep = " :: ")]
  term_order <- ox_metadata[
    ,
    .(
      best_q = min(fixed_effect_greater_q, na.rm = TRUE),
      best_log2 = max(fixed_log2_odds, na.rm = TRUE)
    ),
    by = term
  ][order(best_q, -best_log2), term]
  ox_metadata[, term := factor(term, levels = rev(term_order))]
  p_ox_metadata <- ggplot(
    ox_metadata,
    aes(label_name, term, fill = fixed_log2_odds)
  ) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = n_studies), size = 3) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0
    ) +
    labs(
      title = "Within-study OXPHOS arm metadata support",
      subtitle = "Fill is fixed-effect log2 odds; tile text is supporting study strata.",
      x = NULL,
      y = NULL,
      fill = "log2 odds"
    ) +
    theme_minimal(base_size = 10.3) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1),
      plot.title.position = "plot"
    )
  save_plot_pair(
    p_ox_metadata,
    "dominant_state_score_semantics_oxphos_arm_metadata_fixed_effects",
    width = 11.7,
    height = 8.6
  )
}

ox_cds <- clean_cds_summary[
  label_semantics == oxphos_arm_semantics &
    n_studies >= 2 &
    is.finite(weighted_delta_relative_translation)
]
if (nrow(ox_cds) > 0) {
  top_genes <- ox_cds[
    ,
    .(max_abs = max(abs(weighted_delta_relative_translation), na.rm = TRUE)),
    by = gene_symbol
  ][order(-max_abs), head(gene_symbol, 24L)]
  ox_cds_plot <- ox_cds[gene_symbol %chin% top_genes]
  ox_cds_plot[, gene_symbol := factor(gene_symbol, levels = rev(top_genes))]
  p_ox_cds <- ggplot(
    ox_cds_plot,
    aes(label_name, gene_symbol, fill = weighted_delta_relative_translation)
  ) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = n_studies), size = 3) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0
    ) +
    labs(
      title = "Within-study RDG clean-CDS shifts for OXPHOS arms",
      subtitle = "Fill is study-weighted IN minus OUT clean-CDS relative use; text is studies.",
      x = NULL,
      y = NULL,
      fill = "delta use"
    ) +
    theme_minimal(base_size = 10.3) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1),
      plot.title.position = "plot"
    )
  save_plot_pair(
    p_ox_cds,
    "dominant_state_score_semantics_oxphos_arm_rdg_fixed_effects",
    width = 11.2,
    height = 8.2
  )
}

stopifnot(nrow(label_design_summary) == nrow(label_keys))
stopifnot(all(label_study_counts$out_n >= 0))
stopifnot(all(clean_cds_study_contrasts$feature_id == "clean_CDS"))

message("Design-aware score-semantics comparison:")
message("  metadata fields: ", length(metadata_fields))
message("  mixed label-study rows: ",
        label_study_counts[mixed_study == TRUE, .N])
message("  valid study clean-CDS contrasts: ",
        clean_cds_study_contrasts[is.finite(delta_relative_translation), .N])
message("  saved outputs: ", output_dir)
print(label_design_summary)
print(head(metadata_summary, 30))
print(head(clean_cds_summary, 30))
print(summary_metrics)
