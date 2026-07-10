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

analysis_dir <- if (file.exists("dominant_uorf_feature_expression.csv")) {
  "."
} else if (file.exists(file.path("dominant_cell_states",
                                 "dominant_uorf_feature_expression.csv"))) {
  "dominant_cell_states"
} else {
  "."
}

feature_expression_file <- file.path(analysis_dir, "dominant_uorf_feature_expression.csv")
feature_annotation_file <- file.path(analysis_dir, "dominant_rdg_outputs",
                                     "rdg_feature_annotation.csv")
matched_feature_file <- file.path(analysis_dir, "dominant_next_model",
                                  "dominant_next_model_matched_control_feature_contrasts.csv")
output_dir <- file.path(analysis_dir, "dominant_ouorf_cds_coupling")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

min_matched_group_samples <- 2L
min_matched_feature_raw_counts_any <- 30L
min_matched_feature_raw_counts_both <- 30L

message("ouORF/CDS coupling analysis:")
message("  1. Pair each overlapping uORF with clean_CDS per sample")
message("  2. Separate raw absolute co-expression from total-route residual coupling")
message("  3. Summarize matched-control co-movement when next-model contrasts exist")
message("  4. Save review tables and figures in ", output_dir)

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
}
stop_if_missing(feature_expression_file)
stop_if_missing(feature_annotation_file)

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || uniqueN(x[ok]) < 2L || uniqueN(y[ok]) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(median(x))
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(mean(x))
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(max(x))
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(min(x))
}

collapse_top <- function(x, n = 5L) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  paste(head(sort(x), n), collapse = ";")
}

safe_div <- function(x, y) {
  fifelse(is.finite(x) & is.finite(y) & y != 0, x / y, NA_real_)
}

sign_consistency_p <- function(pos, neg) {
  n <- pos + neg
  fifelse(n > 0, pbinom(pmax(pos, neg) - 1L, n, 0.5,
                         lower.tail = FALSE), NA_real_)
}

raw_count_log2_ratio <- function(case_counts, control_counts,
                                 pseudocount = 0.5) {
  log2((case_counts + pseudocount) / (control_counts + pseudocount))
}

raw_count_log2_ratio_se <- function(case_counts, control_counts,
                                    pseudocount = 0.5) {
  sqrt(1 / (case_counts + pseudocount) +
         1 / (control_counts + pseudocount)) / log(2)
}

raw_count_log2_ratio_z <- function(case_counts, control_counts,
                                   pseudocount = 0.5) {
  ratio <- raw_count_log2_ratio(case_counts, control_counts, pseudocount)
  se <- raw_count_log2_ratio_se(case_counts, control_counts, pseudocount)
  safe_div(ratio, se)
}

residualize_one <- function(x, total) {
  ok <- is.finite(x) & is.finite(total)
  out <- rep(NA_real_, length(x))
  if (sum(ok) < 10L || uniqueN(total[ok]) < 2L || uniqueN(x[ok]) < 2L) {
    return(out)
  }
  out[ok] <- residuals(lm(x[ok] ~ total[ok]))
  out
}

write_plot_pair <- function(plot, name, width = 9, height = 6) {
  png_file <- file.path(figure_dir, paste0(name, ".png"))
  pdf_file <- file.path(figure_dir, paste0(name, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, dpi = 220, bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, bg = "white")
}

write_plotly <- function(plot, name, title = name) {
  if (exists("dominant_save_ggplotly", mode = "function")) {
    dominant_save_ggplotly(
      plot,
      file.path(figure_dir, paste0(name, ".html")),
      title = title
    )
  }
}

annotation <- fread(feature_annotation_file)
ouorf_annotation <- annotation[
  feature_type == "overlapping_uorf" | category == "uoORF",
  .(
    gene_symbol, tx_id, feature_id, feature_type, category,
    translon_sources_merged, tx_start, tx_end, feature_bases,
    cds_overlap_bases, start_codon, start_codon_class,
    frame_relative_to_cds, peptide_length_aa, kozak_strength
  )
]
clean_annotation <- annotation[
  feature_type == "clean_CDS",
  .(
    gene_symbol, tx_id,
    clean_cds_tx_start = tx_start,
    clean_cds_tx_end = tx_end,
    clean_cds_bases = feature_bases
  )
]
if (nrow(ouorf_annotation) == 0) {
  stop("No overlapping_uorf/uoORF annotations found in ", feature_annotation_file,
       call. = FALSE)
}

ouorf_annotation <- merge(
  ouorf_annotation,
  clean_annotation,
  by = c("gene_symbol", "tx_id"),
  all.x = TRUE,
  sort = FALSE
)
ouorf_annotation[, `:=`(
  cds_overlap_fraction_of_clean_cds =
    safe_div(cds_overlap_bases, clean_cds_bases),
  ouorf_start_relative_to_cds_start = tx_start - clean_cds_tx_start,
  ouorf_stop_relative_to_cds_start = tx_end - clean_cds_tx_start,
  downstream_clean_cds_after_ouorf_bases =
    pmax(0, clean_cds_tx_end - pmax(tx_end, clean_cds_tx_start - 1L)),
  stop_inside_cds_or_later = tx_end >= clean_cds_tx_start
)]

feature_cols <- c(
  "gene_symbol", "tx_id", "feature_id", "feature_type", "category",
  "translon_sources_merged", "feature_bases", "cds_overlap_bases",
  "Run", "raw_counts", "fpkm_like", "cds_reference_raw_counts",
  "cds_reference_fpkm_like", "feature_vs_cds_mean_ratio"
)
feature_expression <- fread(feature_expression_file, select = feature_cols)
feature_expression[, raw_counts := as.numeric(raw_counts)]
feature_expression[, fpkm_like := as.numeric(fpkm_like)]

clean_expr <- feature_expression[
  feature_type == "clean_CDS",
  .(
    gene_symbol, tx_id, Run,
    clean_raw_counts = raw_counts,
    clean_fpkm_like = fpkm_like
  )
]
route_totals <- feature_expression[
  feature_type %chin% c("leader_uorf", "overlapping_uorf", "clean_CDS"),
  .(
    route_total_fpkm_like = sum(fpkm_like, na.rm = TRUE),
    route_total_raw_counts = sum(raw_counts, na.rm = TRUE),
    n_route_features = uniqueN(feature_id)
  ),
  by = .(gene_symbol, tx_id, Run)
]
ouorf_expr <- feature_expression[
  feature_type == "overlapping_uorf" | category == "uoORF",
  .(
    gene_symbol, tx_id, feature_id, Run,
    ouorf_raw_counts = raw_counts,
    ouorf_fpkm_like = fpkm_like,
    ouorf_feature_vs_cds_mean_ratio = feature_vs_cds_mean_ratio
  )
]

sample_pairs <- merge(
  ouorf_expr,
  clean_expr,
  by = c("gene_symbol", "tx_id", "Run"),
  all.x = TRUE,
  sort = FALSE
)
sample_pairs <- merge(
  sample_pairs,
  route_totals,
  by = c("gene_symbol", "tx_id", "Run"),
  all.x = TRUE,
  sort = FALSE
)
sample_pairs <- merge(
  sample_pairs,
  ouorf_annotation,
  by = c("gene_symbol", "tx_id", "feature_id"),
  all.x = TRUE,
  sort = FALSE
)

sample_pairs[, `:=`(
  log2_clean_fpkm_like = log2(clean_fpkm_like + 1),
  log2_ouorf_fpkm_like = log2(ouorf_fpkm_like + 1),
  log2_route_total_fpkm_like = log2(route_total_fpkm_like + 1),
  clean_relative_route_use = safe_div(clean_fpkm_like, route_total_fpkm_like),
  ouorf_relative_route_use = safe_div(ouorf_fpkm_like, route_total_fpkm_like),
  clean_raw_share = safe_div(clean_raw_counts, route_total_raw_counts),
  ouorf_raw_share = safe_div(ouorf_raw_counts, route_total_raw_counts),
  pair_evaluable = clean_raw_counts >= 10 & ouorf_raw_counts >= 1 &
    is.finite(clean_fpkm_like) & is.finite(ouorf_fpkm_like) &
    is.finite(route_total_fpkm_like) & route_total_fpkm_like > 0
)]

sample_pairs[
  ,
  clean_total_route_residual := residualize_one(
    log2_clean_fpkm_like,
    log2_route_total_fpkm_like
  ),
  by = .(gene_symbol, tx_id, feature_id)
]
sample_pairs[
  ,
  ouorf_total_route_residual := residualize_one(
    log2_ouorf_fpkm_like,
    log2_route_total_fpkm_like
  ),
  by = .(gene_symbol, tx_id, feature_id)
]

gene_summary <- sample_pairs[
  pair_evaluable == TRUE,
  {
    raw_rho <- safe_cor(log2_clean_fpkm_like, log2_ouorf_fpkm_like)
    residual_rho <- safe_cor(clean_total_route_residual,
                             ouorf_total_route_residual)
    relative_rho <- safe_cor(clean_relative_route_use,
                             ouorf_relative_route_use)
    raw_pearson <- safe_cor(log2_clean_fpkm_like, log2_ouorf_fpkm_like,
                            method = "pearson")
    .(
      n_evaluable_samples = .N,
      n_detected_ouorf_samples = sum(ouorf_raw_counts > 0, na.rm = TRUE),
      total_clean_raw_counts = sum(clean_raw_counts, na.rm = TRUE),
      total_ouorf_raw_counts = sum(ouorf_raw_counts, na.rm = TRUE),
      median_clean_fpkm_like = safe_median(clean_fpkm_like),
      median_ouorf_fpkm_like = safe_median(ouorf_fpkm_like),
      median_route_total_fpkm_like = safe_median(route_total_fpkm_like),
      median_clean_relative_route_use = safe_median(clean_relative_route_use),
      median_ouorf_relative_route_use = safe_median(ouorf_relative_route_use),
      median_ouorf_feature_vs_cds_mean_ratio =
        safe_median(ouorf_feature_vs_cds_mean_ratio),
      raw_log2_fpkm_spearman_rho = raw_rho,
      raw_log2_fpkm_pearson_r = raw_pearson,
      total_route_residual_spearman_rho = residual_rho,
      relative_route_use_spearman_rho = relative_rho,
      raw_minus_residual_rho = raw_rho - residual_rho
    )
  },
  by = .(gene_symbol, tx_id, feature_id)
]
gene_summary <- merge(
  gene_summary,
  unique(ouorf_annotation),
  by = c("gene_symbol", "tx_id", "feature_id"),
  all.x = TRUE,
  sort = FALSE
)
gene_summary[, coupling_class := fcase(
  n_evaluable_samples < 30 | total_ouorf_raw_counts < 100,
  "low_count_review",
  raw_log2_fpkm_spearman_rho >= 0.30 &
    total_route_residual_spearman_rho <= 0.05,
  "raw_positive_shared_route_signal",
  raw_log2_fpkm_spearman_rho >= 0.30 &
    total_route_residual_spearman_rho > 0.05 &
    total_route_residual_spearman_rho < 0.20,
  "raw_positive_weak_residual_review",
  raw_log2_fpkm_spearman_rho >= 0.30 &
    total_route_residual_spearman_rho >= 0.20,
  "positive_residual_coupling_review",
  raw_log2_fpkm_spearman_rho < 0 &
    total_route_residual_spearman_rho < 0,
  "absolute_and_residual_tradeoff",
  default = "mixed_or_weak_coupling"
)]
gene_summary[, truncated_cds_rescue_review := coupling_class ==
               "positive_residual_coupling_review" &
               is.finite(relative_route_use_spearman_rho) &
               relative_route_use_spearman_rho >= 0.20 &
               is.finite(cds_overlap_bases) &
               cds_overlap_bases >= 30 &
               is.finite(frame_relative_to_cds) &
               frame_relative_to_cds != 0 &
               downstream_clean_cds_after_ouorf_bases >= 90]
gene_summary[, interpretation := fcase(
  truncated_cds_rescue_review == TRUE,
  paste(
    "Positive ouORF/CDS residual coupling remains after total-route adjustment;",
    "because the out-of-frame ouORF overlaps the CDS and leaves a downstream",
    "clean-CDS segment, inspect for downstream TIS or N-terminal truncation."
  ),
  coupling_class == "raw_positive_shared_route_signal",
  paste(
    "Raw ouORF and clean-CDS levels rise together, but the positive relation",
    "largely disappears after total route signal is removed; likely shared",
    "transcript abundance, global ribosome loading, or sample/state expression."
  ),
  coupling_class == "positive_residual_coupling_review",
  paste(
    "ouORF and clean-CDS remain positively coupled after total route adjustment;",
    "review browser coverage, downstream start sites, and isoform support."
  ),
  coupling_class == "absolute_and_residual_tradeoff",
  "Both absolute and residual coupling are negative, consistent with a direct tradeoff.",
  coupling_class == "low_count_review",
  "Too few ouORF reads or evaluable samples for a stable coupling interpretation.",
  default = "Mixed or weak coupling; treat as descriptive."
)]

matched_pairs <- data.table()
matched_summary <- data.table()
matched_onoff_pairs <- data.table()
matched_onoff_summary <- data.table()
matched_onoff_design_summary <- data.table()
matched_onoff_sparse_candidates <- data.table()
if (file.exists(matched_feature_file)) {
  matched_cols <- c(
    "gene_symbol", "tx_id", "feature_id", "feature_type", "category",
    "stratum_id", "study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE",
    "case_condition", "control_conditions", "perturbation_class",
    "design_family",
    "n_case_samples", "n_control_samples", "case_mean_relative_use",
    "control_mean_relative_use", "matched_relative_use_delta",
    "matched_relative_use_effect_z", "case_mean_fpkm_like",
    "control_mean_fpkm_like", "matched_fpkm_log2_fc_pseudocount1",
    "case_sum_raw_counts", "control_sum_raw_counts",
    "matched_feature_support"
  )
  matched_available_cols <- names(fread(matched_feature_file, nrows = 0))
  matched <- fread(
    matched_feature_file,
    select = intersect(matched_cols, matched_available_cols)
  )
  if (!"design_family" %chin% names(matched)) {
    matched[, design_family := NA_character_]
  }
  if (!"perturbation_class" %chin% names(matched)) {
    matched[, perturbation_class := NA_character_]
  }
  clean_matched <- matched[
    feature_type == "clean_CDS",
    .(
      gene_symbol, tx_id, stratum_id, study, AUTHOR, CELL_LINE, TISSUE, GENE,
      case_condition, control_conditions, perturbation_class, design_family,
      n_case_samples, n_control_samples,
      clean_case_relative_use = case_mean_relative_use,
      clean_control_relative_use = control_mean_relative_use,
      clean_matched_relative_use_delta = matched_relative_use_delta,
      clean_matched_relative_use_effect_z = matched_relative_use_effect_z,
      clean_case_mean_fpkm_like = case_mean_fpkm_like,
      clean_control_mean_fpkm_like = control_mean_fpkm_like,
      clean_matched_fpkm_log2_fc = matched_fpkm_log2_fc_pseudocount1,
      clean_case_sum_raw_counts = case_sum_raw_counts,
      clean_control_sum_raw_counts = control_sum_raw_counts,
      clean_matched_feature_support = matched_feature_support
    )
  ]
  ouorf_matched <- matched[
    feature_type == "overlapping_uorf" | category == "uoORF",
    .(
      gene_symbol, tx_id, feature_id, stratum_id,
      ouorf_case_relative_use = case_mean_relative_use,
      ouorf_control_relative_use = control_mean_relative_use,
      ouorf_matched_relative_use_delta = matched_relative_use_delta,
      ouorf_matched_relative_use_effect_z = matched_relative_use_effect_z,
      ouorf_case_mean_fpkm_like = case_mean_fpkm_like,
      ouorf_control_mean_fpkm_like = control_mean_fpkm_like,
      ouorf_matched_fpkm_log2_fc = matched_fpkm_log2_fc_pseudocount1,
      ouorf_case_sum_raw_counts = case_sum_raw_counts,
      ouorf_control_sum_raw_counts = control_sum_raw_counts,
      ouorf_matched_feature_support = matched_feature_support
    )
  ]
  matched_pairs <- merge(
    ouorf_matched,
    clean_matched,
    by = c("gene_symbol", "tx_id", "stratum_id"),
    all.x = TRUE,
    sort = FALSE
  )
  matched_pairs <- merge(
    matched_pairs,
    ouorf_annotation,
    by = c("gene_symbol", "tx_id", "feature_id"),
    all.x = TRUE,
    sort = FALSE
  )
  matched_pairs[, min_case_control_samples := pmin(n_case_samples,
                                                   n_control_samples)]
  matched_pairs[, `:=`(
    ouorf_max_pair_raw_counts =
      pmax(ouorf_case_sum_raw_counts, ouorf_control_sum_raw_counts),
    ouorf_min_pair_raw_counts =
      pmin(ouorf_case_sum_raw_counts, ouorf_control_sum_raw_counts),
    clean_max_pair_raw_counts =
      pmax(clean_case_sum_raw_counts, clean_control_sum_raw_counts),
    clean_min_pair_raw_counts =
      pmin(clean_case_sum_raw_counts, clean_control_sum_raw_counts)
  )]
  matched_pairs[, `:=`(
    matched_pair_passes_sample_gate =
      min_case_control_samples >= min_matched_group_samples,
    matched_pair_passes_ouorf_count_gate_any =
      ouorf_max_pair_raw_counts >= min_matched_feature_raw_counts_any,
    matched_pair_passes_clean_count_gate_any =
      clean_max_pair_raw_counts >= min_matched_feature_raw_counts_any,
    matched_pair_passes_ouorf_count_gate_both =
      ouorf_min_pair_raw_counts >= min_matched_feature_raw_counts_both,
    matched_pair_passes_clean_count_gate_both =
      clean_min_pair_raw_counts >= min_matched_feature_raw_counts_both
  )]
  matched_pairs[, `:=`(
    matched_pair_passes_ouorf_count_gate =
      matched_pair_passes_ouorf_count_gate_both,
    matched_pair_passes_clean_count_gate =
      matched_pair_passes_clean_count_gate_both
  )]
  matched_pairs[, matched_pair_evaluable_any :=
                  matched_pair_passes_sample_gate &
                  matched_pair_passes_ouorf_count_gate_any &
                  matched_pair_passes_clean_count_gate_any &
                  is.finite(clean_matched_fpkm_log2_fc) &
                  is.finite(ouorf_matched_fpkm_log2_fc)]
  matched_pairs[, matched_pair_evaluable_strict :=
                  matched_pair_passes_sample_gate &
                  matched_pair_passes_ouorf_count_gate_both &
                  matched_pair_passes_clean_count_gate_both &
                  is.finite(clean_matched_fpkm_log2_fc) &
                  is.finite(ouorf_matched_fpkm_log2_fc)]
  matched_pairs[, matched_pair_evaluable := matched_pair_evaluable_strict]
  matched_pairs[, matched_coupling_pattern := fcase(
    matched_pair_passes_sample_gate != TRUE,
    "low_matched_sample_count",
    matched_pair_passes_ouorf_count_gate_any != TRUE,
    "low_ouorf_counts",
    matched_pair_passes_clean_count_gate_any != TRUE,
    "low_clean_cds_counts",
    matched_pair_passes_ouorf_count_gate_both != TRUE,
    "one_sided_low_ouorf_counts",
    matched_pair_passes_clean_count_gate_both != TRUE,
    "one_sided_low_clean_cds_counts",
    clean_matched_fpkm_log2_fc > 0 & ouorf_matched_fpkm_log2_fc > 0 &
      clean_matched_relative_use_delta < 0 &
      ouorf_matched_relative_use_delta > 0,
    "absolute_coinduction_but_ouorf_allocation_gain",
    clean_matched_fpkm_log2_fc > 0 & ouorf_matched_fpkm_log2_fc > 0,
    "absolute_coinduction",
    clean_matched_fpkm_log2_fc < 0 & ouorf_matched_fpkm_log2_fc > 0,
    "ouorf_up_cds_down",
    clean_matched_fpkm_log2_fc > 0 & ouorf_matched_fpkm_log2_fc < 0,
    "cds_up_ouorf_down",
    clean_matched_fpkm_log2_fc < 0 & ouorf_matched_fpkm_log2_fc < 0,
    "absolute_codepletion",
    default = "mixed_or_near_zero"
  )]

  matched_summary <- matched_pairs[
    matched_pair_evaluable == TRUE,
    {
      abs_coind <- clean_matched_fpkm_log2_fc > 0 &
        ouorf_matched_fpkm_log2_fc > 0
      rel_trade <- clean_matched_relative_use_delta < 0 &
        ouorf_matched_relative_use_delta > 0
      strong_idx <- order(
        -(as.numeric(abs_coind) + as.numeric(rel_trade)),
        -abs(ouorf_matched_fpkm_log2_fc),
        -abs(clean_matched_fpkm_log2_fc)
      )[1]
      .(
        n_matched_pairs = .N,
        n_absolute_coinduction_pairs = sum(abs_coind, na.rm = TRUE),
        fraction_absolute_coinduction = mean(abs_coind, na.rm = TRUE),
        n_ouorf_up_cds_down_pairs = sum(
          ouorf_matched_fpkm_log2_fc > 0 &
            clean_matched_fpkm_log2_fc < 0,
          na.rm = TRUE
        ),
        n_absolute_coinduction_with_ouorf_allocation_gain =
          sum(abs_coind & rel_trade, na.rm = TRUE),
        fraction_absolute_coinduction_with_ouorf_allocation_gain =
          mean(abs_coind & rel_trade, na.rm = TRUE),
        n_relative_ouorf_gain_cds_loss_pairs = sum(rel_trade, na.rm = TRUE),
        fraction_relative_ouorf_gain_cds_loss =
          mean(rel_trade, na.rm = TRUE),
        median_clean_matched_fpkm_log2_fc =
          safe_median(clean_matched_fpkm_log2_fc),
        median_ouorf_matched_fpkm_log2_fc =
          safe_median(ouorf_matched_fpkm_log2_fc),
        median_clean_matched_relative_use_delta =
          safe_median(clean_matched_relative_use_delta),
        median_ouorf_matched_relative_use_delta =
          safe_median(ouorf_matched_relative_use_delta),
        matched_fpkm_log2_fc_spearman_rho =
          safe_cor(clean_matched_fpkm_log2_fc,
                   ouorf_matched_fpkm_log2_fc),
        matched_relative_delta_spearman_rho =
          safe_cor(clean_matched_relative_use_delta,
                   ouorf_matched_relative_use_delta),
        strongest_review_study = study[strong_idx],
        strongest_review_context = paste(
          "CELL_LINE", CELL_LINE[strong_idx],
          "TISSUE", TISSUE[strong_idx],
          "GENE", GENE[strong_idx],
          "CONDITION", case_condition[strong_idx],
          "vs", control_conditions[strong_idx],
          sep = "="
        ),
        strongest_review_pattern = matched_coupling_pattern[strong_idx],
        strongest_review_clean_log2_fc =
          clean_matched_fpkm_log2_fc[strong_idx],
        strongest_review_ouorf_log2_fc =
          ouorf_matched_fpkm_log2_fc[strong_idx],
        strongest_review_clean_relative_delta =
          clean_matched_relative_use_delta[strong_idx],
        strongest_review_ouorf_relative_delta =
          ouorf_matched_relative_use_delta[strong_idx]
      )
    },
    by = .(gene_symbol, tx_id, feature_id)
  ]

  matched_onoff_pairs <- matched_pairs[
    matched_pair_evaluable_any == TRUE &
      matched_pair_evaluable_strict != TRUE &
      is.finite(clean_matched_fpkm_log2_fc) &
      is.finite(ouorf_matched_fpkm_log2_fc)
  ]
  if (nrow(matched_onoff_pairs) > 0) {
    matched_onoff_pairs[, `:=`(
      ouorf_case_on = ouorf_case_sum_raw_counts >=
        min_matched_feature_raw_counts_any &
        ouorf_control_sum_raw_counts < min_matched_feature_raw_counts_both,
      ouorf_control_on = ouorf_control_sum_raw_counts >=
        min_matched_feature_raw_counts_any &
        ouorf_case_sum_raw_counts < min_matched_feature_raw_counts_both,
      clean_case_on = clean_case_sum_raw_counts >=
        min_matched_feature_raw_counts_any &
        clean_control_sum_raw_counts < min_matched_feature_raw_counts_both,
      clean_control_on = clean_control_sum_raw_counts >=
        min_matched_feature_raw_counts_any &
        clean_case_sum_raw_counts < min_matched_feature_raw_counts_both
    )]
    matched_onoff_pairs[, `:=`(
      ouorf_raw_count_log2_ratio =
        raw_count_log2_ratio(ouorf_case_sum_raw_counts,
                             ouorf_control_sum_raw_counts),
      clean_raw_count_log2_ratio =
        raw_count_log2_ratio(clean_case_sum_raw_counts,
                             clean_control_sum_raw_counts),
      ouorf_raw_count_log2_ratio_se =
        raw_count_log2_ratio_se(ouorf_case_sum_raw_counts,
                                ouorf_control_sum_raw_counts),
      clean_raw_count_log2_ratio_se =
        raw_count_log2_ratio_se(clean_case_sum_raw_counts,
                                clean_control_sum_raw_counts)
    )]
    matched_onoff_pairs[, `:=`(
      ouorf_raw_count_z =
        raw_count_log2_ratio_z(ouorf_case_sum_raw_counts,
                               ouorf_control_sum_raw_counts),
      clean_raw_count_z =
        raw_count_log2_ratio_z(clean_case_sum_raw_counts,
                               clean_control_sum_raw_counts)
    )]
    matched_onoff_pairs[, `:=`(
      ouorf_raw_count_p =
        2 * pnorm(-abs(ouorf_raw_count_z)),
      clean_raw_count_p =
        2 * pnorm(-abs(clean_raw_count_z))
    )]
    matched_onoff_pairs[, `:=`(
      ouorf_raw_count_q = p.adjust(ouorf_raw_count_p, "BH"),
      clean_raw_count_q = p.adjust(clean_raw_count_p, "BH")
    )]
    matched_onoff_pairs[, onoff_pattern := fcase(
      ouorf_case_on & clean_case_on,
      "case_on_ouorf_and_clean_cds",
      ouorf_control_on & clean_control_on,
      "control_on_ouorf_and_clean_cds",
      ouorf_case_on & clean_control_on,
      "case_on_ouorf_control_on_clean_cds",
      ouorf_control_on & clean_case_on,
      "case_on_clean_cds_control_on_ouorf",
      ouorf_case_on,
      "case_on_ouorf_only",
      ouorf_control_on,
      "control_on_ouorf_only",
      clean_case_on,
      "case_on_clean_cds_only",
      clean_control_on,
      "control_on_clean_cds_only",
      default = "one_sided_near_threshold"
    )]
    matched_onoff_pairs[, onoff_review_class := fcase(
      (ouorf_case_on | ouorf_control_on) &
        (clean_case_on | clean_control_on),
      "paired_onoff_review",
      ouorf_case_on | ouorf_control_on,
      "ouorf_onoff_review",
      clean_case_on | clean_control_on,
      "clean_cds_onoff_review",
      default = "near_threshold_review"
    )]
    matched_onoff_pairs[, onoff_evidence_score :=
                          0.25 * pmin(1, log10(ouorf_max_pair_raw_counts + 1) / 3) +
                          0.20 * pmin(1, log10(clean_max_pair_raw_counts + 1) / 3.5) +
                          0.20 * pmin(1, min_case_control_samples / 6) +
                          0.20 * pmin(1, abs(ouorf_raw_count_z) / 5) +
                          0.15 * pmin(1, abs(clean_raw_count_z) / 5)]
    matched_onoff_pairs[, onoff_interpretation := fcase(
      onoff_review_class == "paired_onoff_review",
      paste(
        "Both ouORF and clean-CDS have one-sided count support.",
        "Review as a possible on/off routing or transcript-abundance event,"
      ),
      onoff_review_class == "ouorf_onoff_review",
      paste(
        "The ouORF has one-sided count support while clean-CDS passes only",
        "the any-side gate. Treat as ouORF on/off triage, not a default",
        "matched log2FC result."
      ),
      onoff_review_class == "clean_cds_onoff_review",
      paste(
        "The clean-CDS feature has one-sided count support while the ouORF",
        "passes only the any-side gate. Treat as clean-CDS on/off triage."
      ),
      default = "Near-threshold one-sided row retained for audit."
    )]
    matched_onoff_pairs[, `:=`(
      abs_ouorf_raw_count_z = abs(ouorf_raw_count_z),
      abs_clean_raw_count_z = abs(clean_raw_count_z)
    )]
    setorder(
      matched_onoff_pairs,
      -onoff_evidence_score,
      -abs_ouorf_raw_count_z,
      -abs_clean_raw_count_z
    )

    matched_onoff_summary <- matched_onoff_pairs[
      ,
      {
        best_idx <- order(
          -onoff_evidence_score,
          -abs(ouorf_matched_fpkm_log2_fc),
          -abs(clean_matched_fpkm_log2_fc)
        )[1]
        .(
          n_onoff_rows = .N,
          n_paired_onoff_rows =
            sum(onoff_review_class == "paired_onoff_review", na.rm = TRUE),
          n_ouorf_onoff_rows =
            sum(onoff_review_class == "ouorf_onoff_review", na.rm = TRUE),
          n_clean_cds_onoff_rows =
            sum(onoff_review_class == "clean_cds_onoff_review", na.rm = TRUE),
          fraction_case_on_ouorf = mean(ouorf_case_on, na.rm = TRUE),
          fraction_control_on_ouorf = mean(ouorf_control_on, na.rm = TRUE),
          fraction_case_on_clean_cds = mean(clean_case_on, na.rm = TRUE),
          fraction_control_on_clean_cds = mean(clean_control_on, na.rm = TRUE),
          median_ouorf_raw_count_log2_ratio =
            safe_median(ouorf_raw_count_log2_ratio),
          median_clean_raw_count_log2_ratio =
            safe_median(clean_raw_count_log2_ratio),
          max_onoff_evidence_score = safe_max(onoff_evidence_score),
          strongest_onoff_pattern = onoff_pattern[best_idx],
          strongest_onoff_review_class = onoff_review_class[best_idx],
          strongest_onoff_context = paste(
            "STUDY", study[best_idx],
            "CELL_LINE", CELL_LINE[best_idx],
            "TISSUE", TISSUE[best_idx],
            "GENE", GENE[best_idx],
            "CONDITION", case_condition[best_idx],
            "vs", control_conditions[best_idx],
            sep = "="
          ),
          strongest_onoff_ouorf_counts = paste0(
            ouorf_case_sum_raw_counts[best_idx], "/",
            ouorf_control_sum_raw_counts[best_idx]
          ),
          strongest_onoff_clean_counts = paste0(
            clean_case_sum_raw_counts[best_idx], "/",
            clean_control_sum_raw_counts[best_idx]
          )
        )
      },
      by = .(gene_symbol, tx_id, feature_id)
    ]
    setorder(
      matched_onoff_summary,
      -max_onoff_evidence_score,
      -n_paired_onoff_rows,
      gene_symbol,
      feature_id
    )

    matched_onoff_pairs[, `:=`(
      onoff_design_family = fifelse(
        !is.na(design_family) & nzchar(design_family),
        design_family,
        fifelse(
          !is.na(perturbation_class) & nzchar(perturbation_class),
          perturbation_class,
          "unknown"
        )
      ),
      ouorf_onoff_direction = fcase(
        ouorf_case_on, "case_on",
        ouorf_control_on, "control_on",
        default = "near_threshold"
      ),
      clean_onoff_direction = fcase(
        clean_case_on, "case_on",
        clean_control_on, "control_on",
        default = "near_threshold"
      ),
      sparse_branch_switch =
        onoff_pattern %chin% c(
          "case_on_ouorf_control_on_clean_cds",
          "case_on_clean_cds_control_on_ouorf"
        ),
      sparse_shared_onoff =
        onoff_pattern %chin% c(
          "case_on_ouorf_and_clean_cds",
          "control_on_ouorf_and_clean_cds"
        )
    )]

    matched_onoff_design_summary <- matched_onoff_pairs[
      ,
      {
        ouorf_pos <- sum(ouorf_raw_count_log2_ratio > 0, na.rm = TRUE)
        ouorf_neg <- sum(ouorf_raw_count_log2_ratio < 0, na.rm = TRUE)
        clean_pos <- sum(clean_raw_count_log2_ratio > 0, na.rm = TRUE)
        clean_neg <- sum(clean_raw_count_log2_ratio < 0, na.rm = TRUE)
        best_idx <- order(
          -onoff_evidence_score,
          -abs_ouorf_raw_count_z,
          -abs_clean_raw_count_z
        )[1]
        .(
          n_onoff_rows = .N,
          n_studies = uniqueN(study),
          n_cell_lines = uniqueN(CELL_LINE[!is.na(CELL_LINE)]),
          n_tissues = uniqueN(TISSUE[!is.na(TISSUE)]),
          n_case_conditions = uniqueN(case_condition[!is.na(case_condition)]),
          n_control_conditions =
            uniqueN(control_conditions[!is.na(control_conditions)]),
          total_case_samples = sum(n_case_samples, na.rm = TRUE),
          total_control_samples = sum(n_control_samples, na.rm = TRUE),
          median_min_case_control_samples =
            safe_median(min_case_control_samples),
          median_onoff_evidence_score = safe_median(onoff_evidence_score),
          max_onoff_evidence_score = safe_max(onoff_evidence_score),
          median_ouorf_raw_count_log2_ratio =
            safe_median(ouorf_raw_count_log2_ratio),
          median_clean_raw_count_log2_ratio =
            safe_median(clean_raw_count_log2_ratio),
          median_ouorf_raw_count_z = safe_median(ouorf_raw_count_z),
          median_clean_raw_count_z = safe_median(clean_raw_count_z),
          min_ouorf_raw_count_q = safe_min(ouorf_raw_count_q),
          min_clean_raw_count_q = safe_min(clean_raw_count_q),
          max_ouorf_pair_raw_counts = safe_max(ouorf_max_pair_raw_counts),
          max_clean_pair_raw_counts = safe_max(clean_max_pair_raw_counts),
          ouorf_case_on_rows = sum(ouorf_case_on, na.rm = TRUE),
          ouorf_control_on_rows = sum(ouorf_control_on, na.rm = TRUE),
          clean_case_on_rows = sum(clean_case_on, na.rm = TRUE),
          clean_control_on_rows = sum(clean_control_on, na.rm = TRUE),
          ouorf_positive_rows = ouorf_pos,
          ouorf_negative_rows = ouorf_neg,
          clean_positive_rows = clean_pos,
          clean_negative_rows = clean_neg,
          ouorf_direction_consistency =
            safe_div(max(ouorf_pos, ouorf_neg), ouorf_pos + ouorf_neg),
          clean_direction_consistency =
            safe_div(max(clean_pos, clean_neg), clean_pos + clean_neg),
          example_studies = collapse_top(study, 4L),
          example_cell_lines = collapse_top(CELL_LINE, 4L),
          example_tissues = collapse_top(TISSUE, 4L),
          example_case_conditions = collapse_top(case_condition, 4L),
          example_control_conditions = collapse_top(control_conditions, 4L),
          strongest_onoff_study = study[best_idx],
          strongest_onoff_context = paste(
            "CELL_LINE", CELL_LINE[best_idx],
            "TISSUE", TISSUE[best_idx],
            "GENE", GENE[best_idx],
            "CONDITION", case_condition[best_idx],
            "vs", control_conditions[best_idx],
            sep = "="
          ),
          strongest_ouorf_counts = paste0(
            ouorf_case_sum_raw_counts[best_idx], "/",
            ouorf_control_sum_raw_counts[best_idx]
          ),
          strongest_clean_counts = paste0(
            clean_case_sum_raw_counts[best_idx], "/",
            clean_control_sum_raw_counts[best_idx]
          )
        )
      },
      by = .(
        gene_symbol, tx_id, feature_id, onoff_design_family,
        onoff_pattern, onoff_review_class
      )
    ]
    matched_onoff_design_summary[, `:=`(
      ouorf_sign_consistency_p =
        sign_consistency_p(ouorf_positive_rows, ouorf_negative_rows),
      clean_sign_consistency_p =
        sign_consistency_p(clean_positive_rows, clean_negative_rows)
    )]
    matched_onoff_design_summary[, `:=`(
      ouorf_sign_consistency_q = p.adjust(ouorf_sign_consistency_p, "BH"),
      clean_sign_consistency_q = p.adjust(clean_sign_consistency_p, "BH")
    )]
    matched_onoff_design_summary[, `:=`(
      sparse_replication_score =
        0.22 * safe_div(n_studies, n_studies + 4) +
        0.08 * safe_div(n_cell_lines, n_cell_lines + 4) +
        0.10 * safe_div(n_onoff_rows, n_onoff_rows + 15),
      sparse_count_score =
        0.12 * safe_div(log10(max_ouorf_pair_raw_counts + 1),
                        log10(max_ouorf_pair_raw_counts + 1) + 2.5) +
        0.08 * safe_div(log10(max_clean_pair_raw_counts + 1),
                        log10(max_clean_pair_raw_counts + 1) + 3),
      sparse_direction_score =
        0.15 * pmax(
          0,
          2 * pmax(ouorf_direction_consistency,
                   clean_direction_consistency,
                   na.rm = TRUE) - 1
        ),
      sparse_evidence_score =
        0.18 * pmin(1, median_onoff_evidence_score),
      sparse_specificity_score =
        0.03 * fifelse(onoff_design_family == "other_non_control", 0, 1),
      sparse_pairing_score =
        0.04 * fifelse(onoff_review_class == "paired_onoff_review", 1, 0) +
        0.02 * fifelse(
          onoff_pattern %chin% c(
            "case_on_ouorf_control_on_clean_cds",
            "case_on_clean_cds_control_on_ouorf"
          ),
          1,
          0
        )
    )]
    matched_onoff_design_summary[, sparse_onoff_score :=
                                   pmin(
                                     1,
                                     sparse_replication_score +
                                       sparse_count_score +
                                       sparse_direction_score +
                                       sparse_evidence_score +
                                       sparse_specificity_score +
                                       sparse_pairing_score
                                   )]
    matched_onoff_design_summary[
      !is.finite(sparse_onoff_score),
      sparse_onoff_score := 0
    ]
    matched_onoff_design_summary[, sparse_onoff_class := fcase(
      onoff_review_class == "paired_onoff_review" &
        onoff_pattern %chin% c(
          "case_on_ouorf_control_on_clean_cds",
          "case_on_clean_cds_control_on_ouorf"
        ) &
        n_studies >= 2L,
      "sparse_branch_switch_review",
      onoff_review_class == "paired_onoff_review" &
        n_studies >= 3L,
      "replicated_sparse_shared_onoff_review",
      onoff_review_class == "ouorf_onoff_review" &
        onoff_design_family != "other_non_control" &
        n_studies >= 4L &
        n_onoff_rows >= 6L &
        ouorf_direction_consistency >= 0.75,
      "replicated_specific_ouorf_onoff_review",
      onoff_review_class == "ouorf_onoff_review" &
        n_studies >= 6L &
        n_onoff_rows >= 10L &
        ouorf_direction_consistency >= 0.75,
      "broad_replicated_ouorf_onoff_review",
      sparse_onoff_score >= 0.72 &
        n_studies >= 2L,
      "strong_sparse_onoff_audit",
      default = "sparse_onoff_audit"
    )]
    matched_onoff_sparse_candidates <- matched_onoff_design_summary[
      sparse_onoff_class != "sparse_onoff_audit"
    ]
    setorder(
      matched_onoff_design_summary,
      -sparse_onoff_score,
      -n_studies,
      -n_onoff_rows,
      gene_symbol,
      feature_id
    )
    setorder(
      matched_onoff_sparse_candidates,
      -sparse_onoff_score,
      -n_studies,
      -n_onoff_rows,
      gene_symbol,
      feature_id
    )
  }
}

if (nrow(matched_summary) > 0) {
  gene_summary <- merge(
    gene_summary,
    matched_summary,
    by = c("gene_symbol", "tx_id", "feature_id"),
    all.x = TRUE,
    sort = FALSE
  )
} else {
  gene_summary[, `:=`(
    n_matched_pairs = NA_integer_,
    fraction_absolute_coinduction = NA_real_,
    fraction_absolute_coinduction_with_ouorf_allocation_gain = NA_real_,
    fraction_relative_ouorf_gain_cds_loss = NA_real_,
    matched_fpkm_log2_fc_spearman_rho = NA_real_,
    matched_relative_delta_spearman_rho = NA_real_,
    strongest_review_context = NA_character_
  )]
}

gene_summary[, review_priority_score :=
               fifelse(truncated_cds_rescue_review == TRUE, 3, 0) +
               fifelse(coupling_class == "positive_residual_coupling_review", 2, 0) +
               fifelse(coupling_class == "raw_positive_weak_residual_review", 1, 0) +
               fifelse(fraction_absolute_coinduction_with_ouorf_allocation_gain >= 0.20, 1, 0) +
               fifelse(fraction_relative_ouorf_gain_cds_loss >= 0.30, 1, 0) +
               fifelse(translon_sources_merged %like% "M", 1, 0)]
gene_summary[!is.finite(review_priority_score), review_priority_score := 0]

candidate_review <- gene_summary[
  coupling_class %chin% c(
    "positive_residual_coupling_review",
    "raw_positive_weak_residual_review"
  ) |
    truncated_cds_rescue_review == TRUE |
    fraction_absolute_coinduction_with_ouorf_allocation_gain >= 0.20 |
    fraction_relative_ouorf_gain_cds_loss >= 0.30
]
setorder(candidate_review, -review_priority_score,
         -total_route_residual_spearman_rho,
         -fraction_absolute_coinduction_with_ouorf_allocation_gain,
         gene_symbol, feature_id)

setorder(gene_summary, coupling_class, gene_symbol, feature_id)
fwrite(
  gene_summary,
  file.path(output_dir, "dominant_ouorf_cds_coupling_gene_summary.csv")
)
fwrite(
  sample_pairs,
  file.path(output_dir, "dominant_ouorf_cds_coupling_sample_pairs.csv")
)
fwrite(
  matched_pairs,
  file.path(output_dir, "dominant_ouorf_cds_coupling_matched_pairs.csv")
)
fwrite(
  matched_summary,
  file.path(output_dir, "dominant_ouorf_cds_coupling_matched_summary.csv")
)
fwrite(
  matched_onoff_pairs,
  file.path(output_dir, "dominant_ouorf_cds_coupling_onoff_pairs.csv")
)
fwrite(
  matched_onoff_summary,
  file.path(output_dir, "dominant_ouorf_cds_coupling_onoff_summary.csv")
)
fwrite(
  matched_onoff_design_summary,
  file.path(output_dir,
            "dominant_ouorf_cds_coupling_sparse_onoff_design_summary.csv")
)
fwrite(
  matched_onoff_sparse_candidates,
  file.path(output_dir,
            "dominant_ouorf_cds_coupling_sparse_onoff_candidates.csv")
)
fwrite(
  candidate_review,
  file.path(output_dir, "dominant_ouorf_cds_coupling_candidate_review.csv")
)

summary_metrics <- rbindlist(list(
  data.table(metric = "ouorf_features", value = nrow(gene_summary)),
  data.table(metric = "ouorf_features_raw_positive",
             value = sum(gene_summary$raw_log2_fpkm_spearman_rho >= 0.30,
                         na.rm = TRUE)),
  data.table(metric = "ouorf_features_raw_positive_shared_route_signal",
             value = sum(gene_summary$coupling_class ==
                           "raw_positive_shared_route_signal",
                         na.rm = TRUE)),
  data.table(metric = "ouorf_features_positive_residual_review",
             value = sum(gene_summary$coupling_class ==
                           "positive_residual_coupling_review",
                         na.rm = TRUE)),
  data.table(metric = "ouorf_features_truncated_cds_rescue_review",
             value = sum(gene_summary$truncated_cds_rescue_review == TRUE,
                         na.rm = TRUE)),
  data.table(metric = "matched_ouorf_cds_pairs",
             value = nrow(matched_pairs)),
  data.table(metric = "matched_ouorf_cds_sample_gate_pairs",
             value = sum(matched_pairs$matched_pair_passes_sample_gate == TRUE,
                         na.rm = TRUE)),
  data.table(metric = "matched_ouorf_cds_ouorf_count_gate_any_pairs",
             value = sum(matched_pairs$matched_pair_passes_ouorf_count_gate_any == TRUE,
                         na.rm = TRUE)),
  data.table(metric = "matched_ouorf_cds_ouorf_count_gate_both_pairs",
             value = sum(matched_pairs$matched_pair_passes_ouorf_count_gate_both == TRUE,
                         na.rm = TRUE)),
  data.table(metric = "matched_ouorf_cds_clean_count_gate_any_pairs",
             value = sum(matched_pairs$matched_pair_passes_clean_count_gate_any == TRUE,
                         na.rm = TRUE)),
  data.table(metric = "matched_ouorf_cds_clean_count_gate_both_pairs",
             value = sum(matched_pairs$matched_pair_passes_clean_count_gate_both == TRUE,
                         na.rm = TRUE)),
  data.table(metric = "matched_ouorf_cds_evaluable_any_pairs",
             value = sum(matched_pairs$matched_pair_evaluable_any == TRUE,
                         na.rm = TRUE)),
  data.table(metric = "matched_ouorf_cds_evaluable_pairs",
             value = sum(matched_pairs$matched_pair_evaluable == TRUE,
                         na.rm = TRUE)),
  data.table(metric = "matched_ouorf_cds_onoff_review_pairs",
             value = nrow(matched_onoff_pairs)),
  data.table(metric = "matched_ouorf_cds_onoff_review_features",
             value = nrow(matched_onoff_summary)),
  data.table(metric = "matched_ouorf_cds_paired_onoff_review_pairs",
             value = sum(
               matched_onoff_pairs$onoff_review_class == "paired_onoff_review",
               na.rm = TRUE
             )),
  data.table(metric = "matched_ouorf_cds_sparse_onoff_designs",
             value = nrow(matched_onoff_design_summary)),
  data.table(metric = "matched_ouorf_cds_sparse_onoff_candidates",
             value = nrow(matched_onoff_sparse_candidates)),
  data.table(metric = "matched_ouorf_cds_sparse_branch_switch_candidates",
             value = sum(
               matched_onoff_sparse_candidates$sparse_onoff_class ==
                 "sparse_branch_switch_review",
               na.rm = TRUE
             )),
  data.table(metric = "candidate_review_rows",
             value = nrow(candidate_review))
), fill = TRUE)
fwrite(summary_metrics, file.path(output_dir,
                                  "dominant_ouorf_cds_coupling_summary_metrics.csv"))

plot_dt <- gene_summary[
  is.finite(raw_log2_fpkm_spearman_rho) &
    is.finite(total_route_residual_spearman_rho)
]
if (nrow(plot_dt) > 0) {
  plot_dt[, label_gene := fifelse(
    coupling_class == "positive_residual_coupling_review" |
      truncated_cds_rescue_review == TRUE |
      review_priority_score >= 3,
    paste(gene_symbol, feature_id),
    ""
  )]
  p_raw_resid <- ggplot(
    plot_dt,
    aes(
      raw_log2_fpkm_spearman_rho,
      total_route_residual_spearman_rho,
      color = coupling_class,
      size = log10(total_ouorf_raw_counts + 1),
      text = paste0(
        gene_symbol, " ", feature_id,
        "<br>raw rho=", round(raw_log2_fpkm_spearman_rho, 3),
        "<br>residual rho=", round(total_route_residual_spearman_rho, 3),
        "<br>relative-use rho=", round(relative_route_use_spearman_rho, 3),
        "<br>class=", coupling_class
      )
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "#64748b") +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "#64748b") +
    geom_vline(xintercept = 0.30, linewidth = 0.35, linetype = "dashed",
               color = "#64748b") +
    geom_hline(yintercept = 0.20, linewidth = 0.35, linetype = "dashed",
               color = "#64748b") +
    geom_point(alpha = 0.82) +
    geom_text(
      aes(label = label_gene),
      color = "#111827",
      size = 3,
      vjust = -0.8,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        absolute_and_residual_tradeoff = "#0f766e",
        low_count_review = "#94a3b8",
        mixed_or_weak_coupling = "#64748b",
        positive_residual_coupling_review = "#b91c1c",
        raw_positive_shared_route_signal = "#2563eb",
        raw_positive_weak_residual_review = "#f59e0b"
      ),
      drop = FALSE
    ) +
    scale_size_continuous(range = c(2.2, 7.5), name = "log10 ouORF reads") +
    coord_cartesian(xlim = c(-0.25, 1), ylim = c(-1, 0.85)) +
    labs(
      title = "ouORF/CDS coupling separates shared expression from routing effects",
      subtitle = paste(
        "x: raw absolute log-FPKM correlation across samples;",
        "y: correlation after removing total leader+ouORF+clean-CDS route signal."
      ),
      x = "raw ouORF vs clean-CDS Spearman rho",
      y = "total-route residual Spearman rho",
      color = "coupling class"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
  write_plot_pair(p_raw_resid, "ouorf_cds_coupling_raw_vs_residual",
                  width = 10.5, height = 6.6)
  write_plotly(p_raw_resid, "ouorf_cds_coupling_raw_vs_residual",
               "ouORF/CDS raw vs residual coupling")
}

if (nrow(matched_summary) > 0) {
  matched_plot <- copy(gene_summary[
    is.finite(fraction_absolute_coinduction) &
      is.finite(fraction_relative_ouorf_gain_cds_loss)
  ])
  if (nrow(matched_plot) > 0) {
    matched_plot[, label_gene := fifelse(
      truncated_cds_rescue_review == TRUE |
        fraction_absolute_coinduction_with_ouorf_allocation_gain >= 0.20 |
        coupling_class == "positive_residual_coupling_review",
      paste(gene_symbol, feature_id),
      ""
    )]
    p_matched <- ggplot(
      matched_plot,
      aes(
        fraction_absolute_coinduction,
        fraction_relative_ouorf_gain_cds_loss,
        color = coupling_class,
        size = pmax(n_matched_pairs, 1),
        text = paste0(
          gene_symbol, " ", feature_id,
          "<br>absolute coinduction fraction=",
          round(fraction_absolute_coinduction, 3),
          "<br>ouORF relative gain/CDS loss fraction=",
          round(fraction_relative_ouorf_gain_cds_loss, 3),
          "<br>median clean log2FC=",
          round(median_clean_matched_fpkm_log2_fc, 3),
          "<br>median ouORF log2FC=",
          round(median_ouorf_matched_fpkm_log2_fc, 3)
        )
      )
    ) +
      geom_vline(xintercept = 0.5, linewidth = 0.35, linetype = "dashed",
                 color = "#64748b") +
      geom_hline(yintercept = 0.3, linewidth = 0.35, linetype = "dashed",
                 color = "#64748b") +
      geom_point(alpha = 0.82) +
      geom_text(
        aes(label = label_gene),
        color = "#111827",
        size = 3,
        vjust = -0.8,
        check_overlap = TRUE,
        show.legend = FALSE
      ) +
      scale_color_manual(
        values = c(
          absolute_and_residual_tradeoff = "#0f766e",
          low_count_review = "#94a3b8",
          mixed_or_weak_coupling = "#64748b",
          positive_residual_coupling_review = "#b91c1c",
          raw_positive_shared_route_signal = "#2563eb",
          raw_positive_weak_residual_review = "#f59e0b"
        ),
        drop = FALSE
      ) +
      scale_size_continuous(range = c(2.2, 8), name = "matched pairs") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(
        title = "Matched-control ouORF/CDS coupling",
        subtitle = paste(
          "Absolute coinduction can coexist with allocation tradeoff;",
          "upper-right points are browser-review candidates."
        ),
        x = "fraction of matched contrasts where clean CDS and ouORF both rise",
        y = "fraction where ouORF relative use rises while clean-CDS relative use falls",
        color = "sample-level class"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "right",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold")
      )
    write_plot_pair(p_matched, "ouorf_cds_coupling_matched_control_summary",
                    width = 10.5, height = 6.6)
    write_plotly(p_matched, "ouorf_cds_coupling_matched_control_summary",
                 "Matched-control ouORF/CDS coupling")
  }

  matched_row_plot <- matched_pairs[
    matched_pair_evaluable == TRUE &
      is.finite(clean_matched_fpkm_log2_fc) &
      is.finite(ouorf_matched_fpkm_log2_fc)
  ]
  if (nrow(matched_row_plot) > 0) {
    matched_row_plot[, abs_sum_fc := abs(clean_matched_fpkm_log2_fc) +
                       abs(ouorf_matched_fpkm_log2_fc)]
    setorder(matched_row_plot, -abs_sum_fc)
    matched_row_plot <- matched_row_plot[seq_len(min(.N, 900L))]
    p_quadrant <- ggplot(
      matched_row_plot,
      aes(
        clean_matched_fpkm_log2_fc,
        ouorf_matched_fpkm_log2_fc,
        color = matched_coupling_pattern,
        size = min_case_control_samples,
        text = paste0(
          gene_symbol, " ", feature_id,
          "<br>", study,
          "<br>", case_condition, " vs ", control_conditions,
          "<br>clean log2FC=", round(clean_matched_fpkm_log2_fc, 3),
          "<br>ouORF log2FC=", round(ouorf_matched_fpkm_log2_fc, 3),
          "<br>clean relative delta=",
          round(clean_matched_relative_use_delta, 3),
          "<br>ouORF relative delta=",
          round(ouorf_matched_relative_use_delta, 3),
          "<br>ouORF raw case/control=",
          ouorf_case_sum_raw_counts, "/", ouorf_control_sum_raw_counts,
          "<br>clean raw case/control=",
          clean_case_sum_raw_counts, "/", clean_control_sum_raw_counts,
          "<br>strict count gate=", matched_pair_evaluable_strict
        )
      )
    ) +
      geom_hline(yintercept = 0, linewidth = 0.35, color = "#64748b") +
      geom_vline(xintercept = 0, linewidth = 0.35, color = "#64748b") +
      geom_point(alpha = 0.72) +
      scale_color_manual(
        values = c(
          absolute_coinduction = "#2563eb",
          absolute_coinduction_but_ouorf_allocation_gain = "#b91c1c",
          absolute_codepletion = "#0f766e",
          cds_up_ouorf_down = "#7c3aed",
          low_clean_cds_counts = "#cbd5e1",
          low_matched_sample_count = "#94a3b8",
          low_ouorf_counts = "#fca5a5",
          mixed_or_near_zero = "#64748b",
          one_sided_low_clean_cds_counts = "#e2e8f0",
          one_sided_low_ouorf_counts = "#fecaca",
          ouorf_up_cds_down = "#f59e0b"
        ),
        drop = FALSE
      ) +
      scale_size_continuous(range = c(1.4, 6), name = "min n") +
      coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4)) +
      labs(
        title = "Largest matched-control ouORF and clean-CDS absolute shifts",
        subtitle = paste(
          "Shown rows: n>=2 per side; each feature has >=30 raw reads in both case and control.\n",
          "Hover shows case/control raw counts and strict-gate status."
        ),
        x = "clean-CDS matched log2FC",
        y = "ouORF matched log2FC",
        color = "matched pattern"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "right",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold")
      )
    write_plot_pair(p_quadrant, "ouorf_cds_coupling_matched_shift_quadrants",
                    width = 11.5, height = 7)
    write_plotly(p_quadrant, "ouorf_cds_coupling_matched_shift_quadrants",
                 "Matched-control ouORF/CDS shift quadrants")
  }
}

if (nrow(matched_onoff_pairs) > 0) {
  onoff_plot <- copy(matched_onoff_pairs)
  onoff_plot[, abs_count_z_sum := abs(ouorf_raw_count_z) +
               abs(clean_raw_count_z)]
  setorder(onoff_plot, -onoff_evidence_score, -abs_count_z_sum)
  onoff_plot <- onoff_plot[seq_len(min(.N, 500L))]
  p_onoff <- ggplot(
    onoff_plot,
    aes(
      clean_raw_count_log2_ratio,
      ouorf_raw_count_log2_ratio,
      color = onoff_pattern,
      size = min_case_control_samples,
      text = paste0(
        gene_symbol, " ", feature_id,
        "<br>", study,
        "<br>", case_condition, " vs ", control_conditions,
        "<br>on/off class=", onoff_review_class,
        "<br>ouORF raw case/control=",
        ouorf_case_sum_raw_counts, "/", ouorf_control_sum_raw_counts,
        "<br>clean raw case/control=",
        clean_case_sum_raw_counts, "/", clean_control_sum_raw_counts,
        "<br>ouORF count log2 ratio=",
        round(ouorf_raw_count_log2_ratio, 3),
        "<br>clean count log2 ratio=",
        round(clean_raw_count_log2_ratio, 3),
        "<br>evidence score=", round(onoff_evidence_score, 3)
      )
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "#64748b") +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "#64748b") +
    geom_point(alpha = 0.76) +
    scale_color_manual(
      values = c(
        case_on_clean_cds_control_on_ouorf = "#7c3aed",
        case_on_clean_cds_only = "#a78bfa",
        case_on_ouorf_and_clean_cds = "#2563eb",
        case_on_ouorf_control_on_clean_cds = "#b91c1c",
        case_on_ouorf_only = "#f59e0b",
        control_on_clean_cds_only = "#c4b5fd",
        control_on_ouorf_and_clean_cds = "#0f766e",
        control_on_ouorf_only = "#f97316",
        one_sided_near_threshold = "#64748b"
      ),
      drop = FALSE
    ) +
    scale_size_continuous(range = c(1.5, 6), name = "min n") +
    coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8)) +
    labs(
      title = "One-sided matched ouORF/CDS count-support review",
      subtitle = paste(
        "Rows fail the default both-side count gate but pass any-side >=30 counts;",
        "use as on/off triage, not as ordinary matched log2FC evidence."
      ),
      x = "clean-CDS raw count log2(case/control)",
      y = "ouORF raw count log2(case/control)",
      color = "on/off pattern"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
  write_plot_pair(p_onoff, "ouorf_cds_coupling_onoff_review",
                  width = 11.5, height = 7)
  write_plotly(p_onoff, "ouorf_cds_coupling_onoff_review",
               "One-sided matched ouORF/CDS count-support review")
}

if (nrow(matched_onoff_sparse_candidates) > 0) {
  sparse_score_plot <- copy(matched_onoff_sparse_candidates)
  sparse_score_plot[, candidate_label := paste(
    gene_symbol, feature_id, onoff_design_family, onoff_pattern,
    sep = " | "
  )]
  sparse_score_plot <- sparse_score_plot[seq_len(min(.N, 35L))]
  p_sparse_score <- ggplot(
    sparse_score_plot,
    aes(
      sparse_onoff_score,
      reorder(candidate_label, sparse_onoff_score),
      color = sparse_onoff_class,
      size = n_studies,
      text = paste0(
        gene_symbol, " ", feature_id,
        "<br>design family=", onoff_design_family,
        "<br>pattern=", onoff_pattern,
        "<br>class=", sparse_onoff_class,
        "<br>rows=", n_onoff_rows,
        "<br>studies=", n_studies,
        "<br>cell lines=", n_cell_lines,
        "<br>ouORF median count log2 ratio=",
        round(median_ouorf_raw_count_log2_ratio, 3),
        "<br>clean median count log2 ratio=",
        round(median_clean_raw_count_log2_ratio, 3),
        "<br>strongest ouORF counts=", strongest_ouorf_counts,
        "<br>strongest clean counts=", strongest_clean_counts
      )
    )
  ) +
    geom_segment(
      data = sparse_score_plot,
      aes(
        x = 0,
        xend = sparse_onoff_score,
        y = reorder(candidate_label, sparse_onoff_score),
        yend = reorder(candidate_label, sparse_onoff_score)
      ),
      inherit.aes = FALSE,
      linewidth = 0.4,
      color = "#cbd5e1"
    ) +
    geom_point(alpha = 0.88) +
    scale_color_manual(
      values = c(
        broad_replicated_ouorf_onoff_review = "#f97316",
        replicated_sparse_shared_onoff_review = "#2563eb",
        replicated_specific_ouorf_onoff_review = "#b91c1c",
        sparse_branch_switch_review = "#7c3aed",
        strong_sparse_onoff_audit = "#64748b"
      ),
      drop = FALSE
    ) +
    scale_size_continuous(range = c(2, 7), name = "studies") +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      title = "Sparse on/off ouORF/CDS candidate score",
      subtitle = paste(
        "Aggregated from one-sided count-supported matched rows;",
        "not part of the default strict log2FC review."
      ),
      x = "sparse on/off score",
      y = NULL,
      color = "sparse class"
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 8.5)
    )
  write_plot_pair(p_sparse_score,
                  "ouorf_cds_sparse_onoff_candidate_scores",
                  width = 12, height = 8.5)
  write_plotly(p_sparse_score,
               "ouorf_cds_sparse_onoff_candidate_scores",
               "Sparse on/off ouORF/CDS candidate score")
}

if (nrow(matched_onoff_design_summary) > 0) {
  sparse_heat <- copy(matched_onoff_design_summary)
  sparse_heat <- sparse_heat[
    sparse_onoff_class != "sparse_onoff_audit" |
      sparse_onoff_score >= 0.70
  ]
  if (nrow(sparse_heat) > 0) {
    sparse_heat[, gene_feature := paste(gene_symbol, feature_id)]
    setorder(sparse_heat, gene_feature, onoff_design_family,
             -sparse_onoff_score)
    sparse_heat <- sparse_heat[
      ,
      .SD[1],
      by = .(gene_feature, onoff_design_family)
    ]
    top_features <- sparse_heat[
      ,
      .(max_score = max(sparse_onoff_score, na.rm = TRUE)),
      by = gene_feature
    ][order(-max_score)][seq_len(min(.N, 30L)), gene_feature]
    sparse_heat <- sparse_heat[gene_feature %chin% top_features]
    sparse_heat[, gene_feature := factor(
      gene_feature,
      levels = rev(top_features)
    )]
    sparse_heat[, heat_label := paste0("S", n_studies, "\nN", n_onoff_rows)]
    heat_height <- max(5, min(12, 0.34 * length(top_features) + 2))
    p_sparse_heat <- ggplot(
      sparse_heat,
      aes(
        onoff_design_family,
        gene_feature,
        fill = median_ouorf_raw_count_log2_ratio,
        text = paste0(
          gene_feature,
          "<br>design family=", onoff_design_family,
          "<br>pattern=", onoff_pattern,
          "<br>class=", sparse_onoff_class,
          "<br>score=", round(sparse_onoff_score, 3),
          "<br>rows=", n_onoff_rows,
          "<br>studies=", n_studies,
          "<br>median ouORF count log2 ratio=",
          round(median_ouorf_raw_count_log2_ratio, 3),
          "<br>median clean count log2 ratio=",
          round(median_clean_raw_count_log2_ratio, 3)
        )
      )
    ) +
      geom_tile(color = "white", linewidth = 0.35) +
      geom_text(aes(label = heat_label), size = 2.8, lineheight = 0.9) +
      scale_fill_gradient2(
        low = "#2563eb",
        mid = "white",
        high = "#991b1b",
        midpoint = 0,
        name = "median ouORF\ncount log2 ratio"
      ) +
      labs(
        title = "Sparse on/off ouORF support by design family",
        subtitle = "Tile label: S = studies, N = one-sided rows. Fill is median ouORF raw count log2(case/control).",
        x = NULL,
        y = NULL
      ) +
      theme_minimal(base_size = 10.5) +
      theme(
        legend.position = "right",
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(size = 8.5)
      )
    write_plot_pair(p_sparse_heat,
                    "ouorf_cds_sparse_onoff_design_heatmap",
                    width = 10.5, height = heat_height)
    write_plotly(p_sparse_heat,
                 "ouorf_cds_sparse_onoff_design_heatmap",
                 "Sparse on/off ouORF support by design family")
  }
}

message("Saved: ", file.path(output_dir,
                             "dominant_ouorf_cds_coupling_gene_summary.csv"))
message("Saved: ", file.path(output_dir,
                             "dominant_ouorf_cds_coupling_candidate_review.csv"))
message("Saved: ", file.path(output_dir,
                             "dominant_ouorf_cds_coupling_summary_metrics.csv"))
