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
state_file <- file.path(analysis_dir, "human_dominant_cell_states_clean_cds.csv")
audit_score_file <- file.path(
  analysis_dir,
  "dominant_state_score_audit",
  "dominant_state_score_audit_sample_scores_long.csv"
)
audit_annotation_file <- file.path(
  analysis_dir,
  "dominant_state_score_audit",
  "dominant_state_score_audit_calibrated_annotations.csv"
)
feature_file <- file.path(analysis_dir, "dominant_uorf_feature_expression.csv")
selected_fields_file <- file.path(
  analysis_dir,
  "dominant_state_broad_metadata_selected_fields.csv"
)
output_dir <- file.path(analysis_dir, "dominant_state_score_semantics")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- c(
  file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R"),
  file.path("scripts", "dominant_htmlwidgets.R"),
  "dominant_htmlwidgets.R"
)
htmlwidget_helper <- htmlwidget_helper[file.exists(htmlwidget_helper)][1]
if (!is.na(htmlwidget_helper)) source(htmlwidget_helper)

required_files <- c(
  metadata_file, state_file, audit_score_file, audit_annotation_file, feature_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Score-semantics comparison is missing files: ",
       paste(missing_files, collapse = "; "), call. = FALSE)
}

raw_semantics <- "raw_signed_winner"
percentile_semantics <- "within_state_percentile_winner"
oxphos_arm_semantics <- "oxphos_arm_top_decile"
oxphos_state <- "OXPHOS (mitochondrial)"
min_metadata_level_n <- 5L
min_branch_samples <- 3L
min_clean_cds_raw_counts <- 30
oxphos_arm_top_fraction <- 0.10
fallback_metadata_fields <- c("CELL_LINE", "TISSUE", "GENE", "CONDITION")

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

safe_fisher_greater <- function(call_in, no_call_in, call_out, no_call_out) {
  mat <- matrix(c(call_in, no_call_in, call_out, no_call_out),
                nrow = 2, byrow = TRUE)
  if (sum(mat[1, ]) < min_metadata_level_n ||
      sum(mat[2, ]) < min_metadata_level_n ||
      sum(mat[, 1]) < min_branch_samples) {
    return(list(p = NA_real_, odds_ratio = NA_real_))
  }
  test <- tryCatch(
    fisher.test(mat, alternative = "greater"),
    error = function(e) NULL
  )
  if (is.null(test)) return(list(p = NA_real_, odds_ratio = NA_real_))
  list(p = test$p.value, odds_ratio = unname(test$estimate))
}

log2_odds_pseudocount <- function(call_in, no_call_in, call_out, no_call_out) {
  log2(((call_in + 0.5) / (no_call_in + 0.5)) /
         ((call_out + 0.5) / (no_call_out + 0.5)))
}

percentile_rank <- function(x) {
  out <- rep(NA_real_, length(x))
  valid <- is.finite(x)
  if (!any(valid)) return(out)
  out[valid] <- (rank(x[valid], ties.method = "average") - 0.5) / sum(valid)
  out
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

read_metadata_fields <- function(metadata) {
  if (file.exists(selected_fields_file)) {
    fields <- fread(selected_fields_file, showProgress = FALSE)
    if (all(c("field", "used_for_broad_model") %chin% names(fields))) {
      selected <- fields[used_for_broad_model == TRUE, field]
      selected <- selected[selected %chin% names(metadata)]
      if (length(selected) > 0) return(unique(selected))
    }
  }
  intersect(fallback_metadata_fields, names(metadata))
}

make_winner_labels <- function(annotations) {
  raw <- annotations[
    !is.na(raw_winner) & nzchar(raw_winner),
    .(
      Run,
      label_semantics = raw_semantics,
      label_name = raw_winner,
      dominant_state = raw_winner,
      score_measure = "signed_score",
      score_question = "largest raw signed module score"
    )
  ]
  percentile <- annotations[
    !is.na(calibrated_winner) & nzchar(calibrated_winner),
    .(
      Run,
      label_semantics = percentile_semantics,
      label_name = calibrated_winner,
      dominant_state = calibrated_winner,
      score_measure = "within_state_percentile_z",
      score_question = "largest within-state percentile z score"
    )
  ]
  rbind(raw, percentile, fill = TRUE)
}

make_oxphos_arm_labels <- function(audit_scores) {
  ox <- copy(audit_scores[dominant_state == oxphos_state])
  if (nrow(ox) == 0) stop("Score audit lacks OXPHOS rows.", call. = FALSE)
  arm_scores <- melt(
    ox[
      ,
      .(
        Run,
        oxphos_capacity = up_arm_score,
        glycolysis_hypoxia_low = down_low_score,
        signed_oxphos_balance = signed_score
      )
    ],
    id.vars = "Run",
    variable.name = "score_measure",
    value.name = "arm_score",
    variable.factor = FALSE
  )
  arm_scores[, arm_percentile := percentile_rank(arm_score), by = score_measure]
  arm_scores[, branch_in := arm_percentile >= 1 - oxphos_arm_top_fraction]
  arm_names <- c(
    oxphos_capacity = "OXPHOS capacity top decile",
    glycolysis_hypoxia_low = "Glycolysis/hypoxia low top decile",
    signed_oxphos_balance = "Signed OXPHOS balance top decile"
  )
  labels <- arm_scores[
    branch_in == TRUE,
    .(
      Run,
      label_semantics = oxphos_arm_semantics,
      label_name = unname(arm_names[score_measure]),
      dominant_state = oxphos_state,
      score_measure,
      score_question = paste0(
        "top ", round(100 * oxphos_arm_top_fraction),
        "% OXPHOS arm score"
      )
    )
  ]
  list(labels = labels, arm_scores = arm_scores)
}

metadata_call_enrichment <- function(context, labels, metadata_fields) {
  label_summary <- labels[, .(label_n = uniqueN(Run)), by = .(
    label_semantics, label_name, dominant_state, score_measure, score_question
  )]
  rbindlist(lapply(metadata_fields, function(field) {
    values <- clean_level(context[[field]])
    valid <- !is.na(values)
    levels <- sort(unique(values[valid]))
    levels <- levels[vapply(levels, function(level) {
      sum(valid & values == level) >= min_metadata_level_n
    }, logical(1))]
    if (length(levels) == 0) return(NULL)
    rbindlist(lapply(seq_len(nrow(label_summary)), function(i) {
      label <- label_summary[i]
      label_runs <- labels[
        label_semantics == label$label_semantics &
          label_name == label$label_name,
        Run
      ]
      call <- context$Run %chin% label_runs
      rbindlist(lapply(levels, function(level) {
        in_level <- valid & values == level
        out_level <- valid & values != level
        n_in <- sum(in_level)
        n_out <- sum(out_level)
        call_in <- sum(call & in_level)
        call_out <- sum(call & out_level)
        no_call_in <- n_in - call_in
        no_call_out <- n_out - call_out
        fish <- safe_fisher_greater(call_in, no_call_in, call_out, no_call_out)
        data.table(
          label_semantics = label$label_semantics,
          label_name = label$label_name,
          dominant_state = label$dominant_state,
          score_measure = label$score_measure,
          score_question = label$score_question,
          label_n = label$label_n,
          field = field,
          level = level,
          n_in = n_in,
          n_out = n_out,
          calls_in = call_in,
          calls_out = call_out,
          call_fraction_in = call_in / n_in,
          call_fraction_out = call_out / n_out,
          call_fraction_delta = (call_in / n_in) - (call_out / n_out),
          log2_odds_pseudocount = log2_odds_pseudocount(
            call_in, no_call_in, call_out, no_call_out
          ),
          fisher_odds_ratio = fish$odds_ratio,
          fisher_greater_p = fish$p
        )
      }), fill = TRUE)
    }), fill = TRUE)
  }), fill = TRUE)[
    ,
    fisher_greater_p_adj := p.adjust(fisher_greater_p, "BH")
  ][
    order(fisher_greater_p_adj, -log2_odds_pseudocount, label_semantics,
          label_name, field, level)
  ][]
}

aggregate_one_feature_split <- function(feature_dt, label, in_runs, out_runs) {
  if (length(in_runs) < min_branch_samples || length(out_runs) < min_branch_samples) {
    return(data.table())
  }
  rbindlist(lapply(list("in" = in_runs, "out" = out_runs), function(group_runs) {
    feature_dt[
      Run %chin% group_runs,
      .(
        n_samples = uniqueN(Run),
        sum_raw_counts = safe_sum(raw_counts),
        mean_fpkm_like = mean(fpkm_like, na.rm = TRUE)
      ),
      by = .(gene_symbol, tx_id, feature_id, feature_type)
    ]
  }), idcol = "branch_group", fill = TRUE)[
    ,
    `:=`(
      label_semantics = label$label_semantics,
      label_name = label$label_name,
      dominant_state = label$dominant_state,
      score_measure = label$score_measure,
      branch_id = paste(
        label$label_semantics,
        file_safe(label$label_name),
        file_safe(label$score_measure),
        sep = "|"
      )
    )
  ][
    ,
    allocation_feature :=
      feature_type %chin% c("leader_uorf", "overlapping_uorf", "clean_CDS")
  ][
    ,
    total_mean_fpkm_like := sum(mean_fpkm_like, na.rm = TRUE),
    by = .(branch_id, branch_group, gene_symbol)
  ][
    ,
    allocation_total_mean_fpkm_like :=
      sum(mean_fpkm_like[allocation_feature], na.rm = TRUE),
    by = .(branch_id, branch_group, gene_symbol)
  ][
    ,
    clean_cds_raw_counts := sum_raw_counts[feature_id == "clean_CDS"][1],
    by = .(branch_id, branch_group, gene_symbol)
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

aggregate_feature_splits <- function(feature_dt, context, labels) {
  label_summary <- unique(labels[, .(
    label_semantics, label_name, dominant_state, score_measure
  )])
  rbindlist(lapply(seq_len(nrow(label_summary)), function(i) {
    label <- label_summary[i]
    in_runs <- labels[
      label_semantics == label$label_semantics &
        label_name == label$label_name,
      unique(Run)
    ]
    out_runs <- setdiff(context$Run, in_runs)
    aggregate_one_feature_split(feature_dt, label, in_runs, out_runs)
  }), fill = TRUE)
}

contrast_feature_splits <- function(scores) {
  key_cols <- c(
    "branch_id", "label_semantics", "label_name", "dominant_state",
    "score_measure", "gene_symbol", "tx_id", "feature_id", "feature_type"
  )
  in_scores <- scores[branch_group == "in"]
  out_scores <- scores[branch_group == "out"]
  contrasts <- merge(
    in_scores[
      ,
      c(
        key_cols, "n_samples", "mean_fpkm_like", "sum_raw_counts",
        "clean_cds_raw_counts", "passes_clean_cds_count_cutoff",
        "relative_translation"
      ),
      with = FALSE
    ],
    out_scores[
      ,
      c(
        key_cols, "n_samples", "mean_fpkm_like", "sum_raw_counts",
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
  contrasts[, `:=`(
    comparison_passes_clean_cds_count_cutoff =
      passes_clean_cds_count_cutoff_in & passes_clean_cds_count_cutoff_out,
    log2_mean_fpkm_like_in_vs_out =
      log2((mean_fpkm_like_in + 1) / (mean_fpkm_like_out + 1)),
    delta_relative_translation =
      relative_translation_in - relative_translation_out
  )]
  contrasts[
    comparison_passes_clean_cds_count_cutoff == FALSE,
    `:=`(
      log2_mean_fpkm_like_in_vs_out = NA_real_,
      delta_relative_translation = NA_real_
    )
  ]
  contrasts[]
}

compare_state_metadata_terms <- function(enrichment) {
  raw <- enrichment[
    label_semantics == raw_semantics,
    .(
      field, level, dominant_state,
      raw_label_n = label_n,
      raw_calls_in = calls_in,
      raw_log2_odds = log2_odds_pseudocount,
      raw_q = fisher_greater_p_adj
    )
  ]
  percentile <- enrichment[
    label_semantics == percentile_semantics,
    .(
      field, level, dominant_state,
      percentile_label_n = label_n,
      percentile_calls_in = calls_in,
      percentile_log2_odds = log2_odds_pseudocount,
      percentile_q = fisher_greater_p_adj
    )
  ]
  merge(raw, percentile, by = c("field", "level", "dominant_state"),
        all = TRUE, sort = FALSE)[]
}

compare_state_rdg_clean_cds <- function(clean_cds) {
  raw <- clean_cds[
    label_semantics == raw_semantics,
    .(
      gene_symbol, tx_id, dominant_state,
      raw_n_samples_in = n_samples_in,
      raw_delta_relative_translation = delta_relative_translation,
      raw_log2_mean_fpkm_like_in_vs_out = log2_mean_fpkm_like_in_vs_out
    )
  ]
  percentile <- clean_cds[
    label_semantics == percentile_semantics,
    .(
      gene_symbol, tx_id, dominant_state,
      percentile_n_samples_in = n_samples_in,
      percentile_delta_relative_translation = delta_relative_translation,
      percentile_log2_mean_fpkm_like_in_vs_out =
        log2_mean_fpkm_like_in_vs_out
    )
  ]
  merge(raw, percentile, by = c("gene_symbol", "tx_id", "dominant_state"),
        all = TRUE, sort = FALSE)[]
}

summarize_state_metadata_comparison <- function(comparison) {
  comparison[
    raw_label_n >= 10 & percentile_label_n >= 10 &
      is.finite(raw_log2_odds) & is.finite(percentile_log2_odds),
    .(
      metadata_terms = .N,
      raw_label_n = first(raw_label_n),
      percentile_label_n = first(percentile_label_n),
      log2_odds_spearman = cor(raw_log2_odds, percentile_log2_odds,
                               method = "spearman"),
      median_abs_log2_odds_change =
        median(abs(percentile_log2_odds - raw_log2_odds)),
      q95_abs_log2_odds_change =
        quantile(abs(percentile_log2_odds - raw_log2_odds), 0.95),
      raw_q05_terms = sum(is.finite(raw_q) & raw_q < 0.05),
      percentile_q05_terms =
        sum(is.finite(percentile_q) & percentile_q < 0.05)
    ),
    by = dominant_state
  ][order(log2_odds_spearman)][]
}

summarize_state_rdg_comparison <- function(comparison) {
  comparison[
    raw_n_samples_in >= min_branch_samples &
      percentile_n_samples_in >= min_branch_samples &
      is.finite(raw_delta_relative_translation) &
      is.finite(percentile_delta_relative_translation),
    .(
      clean_cds_genes = .N,
      raw_n_samples_in = first(raw_n_samples_in),
      percentile_n_samples_in = first(percentile_n_samples_in),
      clean_cds_delta_spearman =
        cor(raw_delta_relative_translation,
            percentile_delta_relative_translation, method = "spearman"),
      median_abs_clean_cds_delta_change =
        median(abs(percentile_delta_relative_translation -
                     raw_delta_relative_translation)),
      q95_abs_clean_cds_delta_change =
        quantile(abs(percentile_delta_relative_translation -
                       raw_delta_relative_translation), 0.95)
    ),
    by = dominant_state
  ][order(-q95_abs_clean_cds_delta_change)][]
}

annotations <- fread(audit_annotation_file, showProgress = FALSE)
audit_scores <- fread(audit_score_file, showProgress = FALSE)
state_dt <- fread(state_file, showProgress = FALSE)
metadata <- fread(metadata_file, showProgress = FALSE)
feature_dt <- fread(feature_file, showProgress = FALSE)

if (metadata[Run %chin% state_dt$Run, anyDuplicated(Run)] > 0) {
  stop("Metadata has duplicated Runs matching score-semantics context.",
       call. = FALSE)
}
if (!all(c("raw_counts", "fpkm_like", "Run", "feature_id", "feature_type",
           "gene_symbol", "tx_id") %chin% names(feature_dt))) {
  stop("uORF feature-expression table lacks score-semantics columns.",
       call. = FALSE)
}
feature_dt <- feature_dt[
  is.finite(fpkm_like) & feature_type %chin%
    c("leader_uorf", "overlapping_uorf", "clean_CDS", "internal_orf")
]
feature_dt[, raw_counts := as.numeric(raw_counts)]

winner_labels <- make_winner_labels(annotations)
oxphos_arm <- make_oxphos_arm_labels(audit_scores)
labels <- rbind(winner_labels, oxphos_arm$labels, fill = TRUE)
labels[, branch_label_n := uniqueN(Run), by = .(
  label_semantics, label_name, dominant_state, score_measure
)]
setorder(labels, label_semantics, label_name, Run)

context <- merge(
  state_dt[, .(Run, sample)],
  metadata,
  by = "Run",
  all.x = TRUE,
  sort = FALSE
)
context <- merge(
  state_dt[, .(Run, input_order = .I)],
  context,
  by = "Run",
  all.x = TRUE,
  sort = FALSE
)
setorder(context, input_order)
context[, input_order := NULL]

metadata_fields <- read_metadata_fields(metadata)
call_enrichment <- metadata_call_enrichment(context, labels, metadata_fields)
metadata_term_comparison <- compare_state_metadata_terms(call_enrichment)

feature_scores <- aggregate_feature_splits(feature_dt, context, labels)
feature_contrasts <- contrast_feature_splits(feature_scores)
clean_cds_contrasts <- feature_contrasts[feature_id == "clean_CDS"]
rdg_state_comparison <- compare_state_rdg_clean_cds(clean_cds_contrasts)
metadata_state_summary <- summarize_state_metadata_comparison(
  metadata_term_comparison
)
rdg_state_summary <- summarize_state_rdg_comparison(rdg_state_comparison)

label_summary <- unique(labels[, .(
  label_semantics, label_name, dominant_state, score_measure, score_question,
  branch_label_n
)])
label_summary <- merge(
  label_summary,
  call_enrichment[
    is.finite(fisher_greater_p_adj) & fisher_greater_p_adj < 0.05,
    .(metadata_q05_terms = .N),
    by = .(label_semantics, label_name)
  ],
  by = c("label_semantics", "label_name"),
  all.x = TRUE,
  sort = FALSE
)
label_summary <- merge(
  label_summary,
  clean_cds_contrasts[
    is.finite(delta_relative_translation),
    .(rdg_clean_cds_genes = uniqueN(gene_symbol)),
    by = .(label_semantics, label_name)
  ],
  by = c("label_semantics", "label_name"),
  all.x = TRUE,
  sort = FALSE
)
label_summary[is.na(metadata_q05_terms), metadata_q05_terms := 0L]
label_summary[is.na(rdg_clean_cds_genes), rdg_clean_cds_genes := 0L]
setorder(label_summary, label_semantics, -branch_label_n, label_name)

summary_metrics <- data.table(
  metric = c(
    "metadata_fields_tested",
    "branch_labels",
    "metadata_call_enrichment_rows",
    "metadata_q05_rows",
    "rdg_feature_contrast_rows",
    "rdg_clean_cds_contrast_rows",
    "rdg_clean_cds_valid_rows",
    "raw_percentile_metadata_term_rows",
    "raw_percentile_rdg_clean_cds_rows"
  ),
  value = c(
    length(metadata_fields),
    nrow(label_summary),
    nrow(call_enrichment),
    call_enrichment[
      is.finite(fisher_greater_p_adj) & fisher_greater_p_adj < 0.05, .N
    ],
    nrow(feature_contrasts),
    nrow(clean_cds_contrasts),
    clean_cds_contrasts[is.finite(delta_relative_translation), .N],
    nrow(metadata_term_comparison),
    nrow(rdg_state_comparison)
  )
)

fwrite(
  labels,
  file.path(output_dir, "dominant_state_score_semantics_labels.csv")
)
fwrite(
  oxphos_arm$arm_scores,
  file.path(output_dir, "dominant_state_score_semantics_oxphos_arm_scores.csv")
)
fwrite(
  label_summary,
  file.path(output_dir, "dominant_state_score_semantics_label_summary.csv")
)
fwrite(
  call_enrichment,
  file.path(output_dir, "dominant_state_score_semantics_metadata_call_enrichment.csv")
)
fwrite(
  metadata_term_comparison,
  file.path(output_dir, "dominant_state_score_semantics_metadata_raw_vs_percentile.csv")
)
fwrite(
  metadata_state_summary,
  file.path(output_dir, "dominant_state_score_semantics_metadata_state_summary.csv")
)
fwrite(
  feature_scores,
  file.path(output_dir, "dominant_state_score_semantics_rdg_feature_scores.csv")
)
fwrite(
  feature_contrasts,
  file.path(output_dir, "dominant_state_score_semantics_rdg_feature_contrasts.csv")
)
fwrite(
  clean_cds_contrasts,
  file.path(output_dir, "dominant_state_score_semantics_rdg_clean_cds_contrasts.csv")
)
fwrite(
  rdg_state_comparison,
  file.path(output_dir, "dominant_state_score_semantics_rdg_raw_vs_percentile.csv")
)
fwrite(
  rdg_state_summary,
  file.path(output_dir, "dominant_state_score_semantics_rdg_state_summary.csv")
)
fwrite(
  summary_metrics,
  file.path(output_dir, "dominant_state_score_semantics_summary_metrics.csv")
)

metadata_plot <- metadata_term_comparison[
  raw_label_n >= 10 & percentile_label_n >= 10 &
    is.finite(raw_log2_odds) & is.finite(percentile_log2_odds)
]
if (nrow(metadata_plot) > 0) {
  metadata_plot[, either_q05 :=
    (is.finite(raw_q) & raw_q < 0.05) |
      (is.finite(percentile_q) & percentile_q < 0.05)]
  p_metadata <- ggplot(
    metadata_plot,
    aes(raw_log2_odds, percentile_log2_odds, color = either_q05)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "#808080") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "#808080") +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.45,
                linetype = "dashed", color = "#303030") +
    geom_point(alpha = 0.45, size = 0.9) +
    facet_wrap(~ dominant_state, scales = "free", ncol = 3) +
    scale_color_manual(
      values = c("FALSE" = "#8293A5", "TRUE" = "#B24A3C"),
      labels = c("FALSE" = "q >= 0.05", "TRUE" = "q < 0.05")
    ) +
    labs(
      title = "Metadata call enrichment moves under winner-score semantics",
      subtitle = paste(
        "Each point is one metadata term. Axes are Fisher log2 odds",
        "for raw signed winners and within-state percentile winners."
      ),
      x = "Raw signed winner metadata log2 odds",
      y = "Percentile winner metadata log2 odds",
      color = "Either label"
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      plot.title.position = "plot"
    )
  save_plot_pair(
    p_metadata,
    "dominant_state_score_semantics_metadata_raw_vs_percentile",
    width = 12.2,
    height = 9.2,
    html = FALSE
  )
}

oxphos_metadata <- call_enrichment[
  label_semantics == oxphos_arm_semantics &
    is.finite(fisher_greater_p_adj) &
    is.finite(log2_odds_pseudocount)
][
  order(fisher_greater_p_adj, -log2_odds_pseudocount),
  head(.SD, 10L),
  by = label_name
]
if (nrow(oxphos_metadata) > 0) {
  oxphos_metadata[, term := paste(field, level, sep = " :: ")]
  term_order <- oxphos_metadata[
    ,
    .(best_q = min(fisher_greater_p_adj, na.rm = TRUE),
      best_odds = max(log2_odds_pseudocount, na.rm = TRUE)),
    by = term
  ][order(best_q, -best_odds), term]
  oxphos_metadata[, term := factor(term, levels = rev(term_order))]
  p_ox_metadata <- ggplot(
    oxphos_metadata,
    aes(label_name, term, fill = log2_odds_pseudocount)
  ) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = calls_in), size = 3) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0
    ) +
    labs(
      title = "OXPHOS arm labels expose different metadata contexts",
      subtitle = "Fill is Fisher log2 odds; tile text is top-decile label calls in the term.",
      x = NULL,
      y = NULL,
      fill = "log2 odds"
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1),
      plot.title.position = "plot"
    )
  save_plot_pair(
    p_ox_metadata,
    "dominant_state_score_semantics_oxphos_arm_metadata",
    width = 11.4,
    height = 8.5
  )
}

rdg_plot <- rdg_state_comparison[
  raw_n_samples_in >= min_branch_samples &
    percentile_n_samples_in >= min_branch_samples &
    is.finite(raw_delta_relative_translation) &
    is.finite(percentile_delta_relative_translation)
]
if (nrow(rdg_plot) > 0) {
  p_rdg <- ggplot(
    rdg_plot,
    aes(
      raw_delta_relative_translation,
      percentile_delta_relative_translation,
      color = dominant_state
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "#808080") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "#808080") +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.45,
                linetype = "dashed", color = "#303030") +
    geom_point(alpha = 0.56, size = 1.2) +
    labs(
      title = "Winner semantics change RDG clean-CDS split estimates",
      subtitle = "Each point is one gene-state clean-CDS relative-use contrast.",
      x = "Raw winner IN minus OUT clean-CDS relative use",
      y = "Percentile winner IN minus OUT clean-CDS relative use",
      color = "State"
    ) +
    theme_minimal(base_size = 10.8) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title.position = "plot"
    )
  save_plot_pair(
    p_rdg,
    "dominant_state_score_semantics_rdg_raw_vs_percentile",
    width = 11.6,
    height = 7.6
  )
}

oxphos_clean_cds <- clean_cds_contrasts[
  label_semantics == oxphos_arm_semantics &
    is.finite(delta_relative_translation)
]
if (nrow(oxphos_clean_cds) > 0) {
  top_oxphos_genes <- oxphos_clean_cds[
    ,
    .(max_abs_delta = max(abs(delta_relative_translation), na.rm = TRUE)),
    by = gene_symbol
  ][order(-max_abs_delta), head(gene_symbol, 24L)]
  plot_ox_cds <- oxphos_clean_cds[gene_symbol %chin% top_oxphos_genes]
  plot_ox_cds[, gene_symbol := factor(
    gene_symbol,
    levels = rev(top_oxphos_genes)
  )]
  p_ox_cds <- ggplot(
    plot_ox_cds,
    aes(label_name, gene_symbol, fill = delta_relative_translation)
  ) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0
    ) +
    labs(
      title = "RDG clean-CDS shifts across OXPHOS arm labels",
      subtitle = "IN minus OUT clean-CDS relative use for the strongest arm-split genes.",
      x = NULL,
      y = NULL,
      fill = "delta use"
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1),
      plot.title.position = "plot"
    )
  save_plot_pair(
    p_ox_cds,
    "dominant_state_score_semantics_oxphos_arm_rdg_clean_cds",
    width = 10.8,
    height = 8
  )
}

stopifnot(uniqueN(labels[label_semantics == raw_semantics, Run]) == nrow(state_dt))
stopifnot(uniqueN(labels[label_semantics == percentile_semantics, Run]) == nrow(state_dt))
stopifnot(all(clean_cds_contrasts$feature_id == "clean_CDS"))

message("Dominant-state score-semantics comparison:")
message("  metadata fields tested: ", length(metadata_fields))
message("  branch labels: ", nrow(label_summary))
message("  RDG clean-CDS valid contrasts: ",
        clean_cds_contrasts[is.finite(delta_relative_translation), .N])
message("  saved outputs: ", output_dir)
print(label_summary)
print(summary_metrics)
print(metadata_state_summary)
print(rdg_state_summary)
