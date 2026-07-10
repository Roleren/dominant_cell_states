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
  library(randomForest)
})

htmlwidget_helper <- c(
  file.path("scripts", "dominant_htmlwidgets.R"),
  "dominant_htmlwidgets.R",
  file.path("dominant_cell_states", "scripts", "dominant_htmlwidgets.R")
)
htmlwidget_helper <- htmlwidget_helper[file.exists(htmlwidget_helper)][1]
if (!is.na(htmlwidget_helper)) source(htmlwidget_helper)

analysis_dir <- if (file.exists("dominant_rdg_outputs/rdg_feature_annotation.csv")) {
  "."
} else if (file.exists("dominant_cell_states/dominant_rdg_outputs/rdg_feature_annotation.csv")) {
  "dominant_cell_states"
} else {
  "."
}

feature_annotation_file <- file.path(analysis_dir, "dominant_rdg_outputs/rdg_feature_annotation.csv")
branch_usage_file <- file.path(analysis_dir, "dominant_rdg_outputs/relative_usage_matrices/relative_usage_in_wide.csv")
feature_expression_file <- file.path(analysis_dir, "dominant_uorf_feature_expression.csv")
output_dir <- file.path(analysis_dir, "dominant_uorf_structure_rules")
min_clean_cds_raw_counts <- 30

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

message("uORF structure rule inference:")
message("  1. Build transcript-geometry predictors from RDG feature annotation")
message("  2. Summarize clean-CDS relative usage at sample, branch, and gene level")
message("  3. Fit descriptive rule models for uORF count, spacing, length, and overlap")
message("  4. Save CSV tables and figures in ", output_dir)

logit_clip <- function(x, eps = 1e-4) {
  qlogis(pmin(pmax(x, eps), 1 - eps))
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

safe_min <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

finite_complete_rows <- function(dt, cols) {
  if (length(cols) == 0) return(rep(TRUE, nrow(dt)))
  keep <- rep(TRUE, nrow(dt))
  for (col in cols) {
    v <- dt[[col]]
    if (is.numeric(v)) {
      keep <- keep & is.finite(v)
    } else {
      keep <- keep & !is.na(v)
    }
  }
  keep
}

union_width <- function(starts, ends) {
  starts <- as.integer(starts)
  ends <- as.integer(ends)
  valid <- is.finite(starts) & is.finite(ends)
  if (!any(valid)) return(0L)
  starts <- starts[valid]
  ends <- ends[valid]
  length(unique(unlist(Map(seq.int, starts, ends), use.names = FALSE)))
}

clean_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  substr(x, 1L, 160L)
}

write_plot <- function(plot, name, width = 9, height = 6) {
  png_file <- file.path(output_dir, "figures", paste0(name, ".png"))
  pdf_file <- file.path(output_dir, "figures", paste0(name, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, dpi = 220, bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, bg = "white")
}

write_plotly <- function(plot, name) {
  html_file <- file.path(output_dir, "figures", paste0(name, ".html"))
  dominant_save_ggplotly(plot, html_file, title = name)
}

weighted_mean_or_na <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  weighted.mean(x[ok], w[ok])
}

annotation <- fread(feature_annotation_file)
branch_usage <- fread(branch_usage_file)
feature_expression <- fread(feature_expression_file)

uorf_annotation <- annotation[feature_type != "clean_CDS" & category %chin% c("uORF", "uoORF")]
clean_cds_annotation <- annotation[feature_type == "clean_CDS"]

if (nrow(uorf_annotation) == 0) {
  stop("No modeled uORF features found in ", feature_annotation_file)
}
if (nrow(clean_cds_annotation) == 0) {
  stop("No clean_CDS features found in ", feature_annotation_file)
}

make_structure_one <- function(gene, tx) {
  u <- copy(uorf_annotation[gene_symbol == gene & tx_id == tx])
  cds <- clean_cds_annotation[gene_symbol == gene & tx_id == tx][1]
  setorder(u, tx_start, tx_end)

  n_u <- nrow(u)
  leader <- u[feature_type == "leader_uorf"]
  overlap <- u[feature_type == "overlapping_uorf" | cds_overlap_bases > 0]
  gaps <- if (n_u >= 2) u$tx_start[-1] - u$tx_end[-n_u] - 1 else numeric()
  start_spacing <- if (n_u >= 2) diff(u$tx_start) else numeric()
  leader_length <- max(cds$tx_start - 1L, 0L)
  uorf_union <- union_width(u$tx_start, u$tx_end)
  leader_union <- union_width(leader$tx_start, leader$tx_end)
  overlap_union <- union_width(overlap$tx_start, overlap$tx_end)
  last_leader_end <- if (nrow(leader) > 0) max(leader$tx_end, na.rm = TRUE) else NA_real_
  last_any_end <- if (n_u > 0) max(u$tx_end, na.rm = TRUE) else NA_real_
  first_uorf_bases <- if (n_u > 0) u$feature_bases[1] else NA_real_
  first_uorf_to_cds_gap <- if (n_u > 0) cds$tx_start - u$tx_end[1] - 1 else NA_real_
  leader_stop_to_cds_gaps <- if (nrow(leader) > 0) {
    cds$tx_start - leader$tx_end - 1
  } else {
    numeric()
  }
  last_leader_stop_to_clean_cds_gap <- if (is.finite(last_leader_end)) {
    cds$tx_start - last_leader_end - 1
  } else {
    NA_real_
  }
  last_any_stop_to_clean_cds_gap <- if (is.finite(last_any_end)) {
    cds$tx_start - last_any_end - 1
  } else {
    NA_real_
  }
  sources <- as.character(u$translon_sources_merged)
  sources[is.na(sources)] <- "unknown"
  max_cds_overlap_bases <- if (nrow(overlap) > 0) max(overlap$cds_overlap_bases, na.rm = TRUE) else 0
  total_cds_overlap_bases_double_counted <- if (nrow(overlap) > 0) {
    sum(overlap$cds_overlap_bases, na.rm = TRUE)
  } else {
    0
  }
  max_cds_overlap_fraction_of_clean_cds <- if (cds$feature_bases > 0) {
    max_cds_overlap_bases / cds$feature_bases
  } else {
    NA_real_
  }

  data.table(
    gene_symbol = gene,
    tx_id = tx,
    n_uorfs = n_u,
    n_leader_uorfs = nrow(leader),
    n_overlapping_uorfs = nrow(overlap),
    n_uorfs_t_only = sum(sources == "T"),
    n_uorfs_tc_only = sum(sources == "TC"),
    n_uorfs_m_only = sum(sources == "M"),
    n_uorfs_t_tc_shared = sum(sources == "T+TC"),
    n_uorfs_with_tc_support = sum(grepl("TC", sources, fixed = TRUE)),
    n_uorfs_with_m_support = sum(grepl("M", sources, fixed = TRUE)),
    uorf_source_summary = paste(names(table(sources)), as.integer(table(sources)),
                                sep = "=", collapse = ";"),
    has_overlapping_uorf = nrow(overlap) > 0,
    n_overlapping_or_nested_pairs = sum(gaps < 0, na.rm = TRUE),
    n_abutting_pairs = sum(gaps == 0, na.rm = TRUE),
    n_close_pairs_10bp = sum(gaps <= 10, na.rm = TRUE),
    n_close_pairs_30bp = sum(gaps <= 30, na.rm = TRUE),
    n_close_pairs_60bp = sum(gaps <= 60, na.rm = TRUE),
    min_inter_uorf_gap = safe_min(gaps),
    median_inter_uorf_gap = safe_median(gaps),
    mean_inter_uorf_gap = safe_mean(gaps),
    min_start_spacing = safe_min(start_spacing),
    leader_length = leader_length,
    clean_cds_bases = cds$feature_bases,
    total_uorf_bases_double_counted = sum(u$feature_bases, na.rm = TRUE),
    uorf_union_bases = uorf_union,
    leader_uorf_union_bases = leader_union,
    overlapping_uorf_union_bases = overlap_union,
    uorf_union_fraction_of_leader = ifelse(leader_length > 0, uorf_union / leader_length, NA_real_),
    leader_uorf_union_fraction_of_leader = ifelse(leader_length > 0, leader_union / leader_length, NA_real_),
    longest_uorf_bases = safe_max(u$feature_bases),
    mean_uorf_bases = safe_mean(u$feature_bases),
    first_uorf_bases = first_uorf_bases,
    first_uorf_is_short_30bp = is.finite(first_uorf_bases) && first_uorf_bases <= 30,
    first_uorf_is_short_90bp = is.finite(first_uorf_bases) && first_uorf_bases <= 90,
    first_uorf_to_cds_gap = first_uorf_to_cds_gap,
    min_leader_stop_to_clean_cds_gap = safe_min(leader_stop_to_cds_gaps),
    n_leader_uorfs_stop_within_30bp_cds = sum(leader_stop_to_cds_gaps <= 30, na.rm = TRUE),
    n_leader_uorfs_stop_within_60bp_cds = sum(leader_stop_to_cds_gaps <= 60, na.rm = TRUE),
    n_leader_uorfs_stop_within_100bp_cds = sum(leader_stop_to_cds_gaps <= 100, na.rm = TRUE),
    last_leader_stop_to_clean_cds_gap = last_leader_stop_to_clean_cds_gap,
    last_any_stop_to_clean_cds_gap = last_any_stop_to_clean_cds_gap,
    max_cds_overlap_bases = max_cds_overlap_bases,
    total_cds_overlap_bases_double_counted = total_cds_overlap_bases_double_counted,
    max_cds_overlap_fraction_of_clean_cds = max_cds_overlap_fraction_of_clean_cds
  )
}

gene_tx <- unique(annotation[, .(gene_symbol, tx_id)])
structure_summary <- rbindlist(
  lapply(seq_len(nrow(gene_tx)), function(i) {
    make_structure_one(gene_tx$gene_symbol[i], gene_tx$tx_id[i])
  }),
  fill = TRUE
)

structure_summary[, uorf_count_bin := fifelse(
  n_uorfs == 0, "0",
  fifelse(n_uorfs == 1, "1",
          fifelse(n_uorfs == 2, "2",
                  fifelse(n_uorfs <= 4, "3-4", "5+")))
)]
structure_summary[, uorf_count_bin := factor(uorf_count_bin, levels = c("0", "1", "2", "3-4", "5+"))]
structure_summary[, overlap_class := fifelse(
  n_uorfs == 0,
  "no modeled uORF",
  fifelse(has_overlapping_uorf, "has uoORF/CDS overlap", "leader-only uORFs")
)]
structure_summary[, overlap_class := factor(
  overlap_class,
  levels = c("no modeled uORF", "leader-only uORFs", "has uoORF/CDS overlap")
)]
structure_summary[, spacing_class := fifelse(
  n_uorfs == 0, "no uORF",
  fifelse(n_uorfs == 1, "single uORF",
          fifelse(min_inter_uorf_gap < 0, "overlapping/nested uORFs",
                  fifelse(min_inter_uorf_gap <= 10, "0-10 bp between uORFs",
                          fifelse(min_inter_uorf_gap <= 30, "11-30 bp between uORFs",
                                  fifelse(min_inter_uorf_gap <= 60, "31-60 bp between uORFs",
                                          ">60 bp between uORFs")))))
)]
structure_summary[, spacing_class := factor(
  spacing_class,
  levels = c(
    "no uORF",
    "single uORF",
    "overlapping/nested uORFs",
    "0-10 bp between uORFs",
    "11-30 bp between uORFs",
    "31-60 bp between uORFs",
    ">60 bp between uORFs"
  )
)]
structure_summary[, last_leader_gap_class := fifelse(
  is.na(last_leader_stop_to_clean_cds_gap), "no leader uORF",
  fifelse(last_leader_stop_to_clean_cds_gap <= 10, "<=10 bp to clean CDS",
          fifelse(last_leader_stop_to_clean_cds_gap <= 50, "11-50 bp to clean CDS",
                  fifelse(last_leader_stop_to_clean_cds_gap <= 150, "51-150 bp to clean CDS",
                          ">150 bp to clean CDS")))
)]
structure_summary[, first_uorf_length_class := fifelse(
  n_uorfs == 0 | is.na(first_uorf_bases), "no uORF",
  fifelse(first_uorf_bases <= 30, "<=30 bp",
          fifelse(first_uorf_bases <= 90, "31-90 bp", ">90 bp"))
)]
structure_summary[, has_modeled_uorf := n_uorfs > 0]

feature_expression[, raw_counts := as.numeric(raw_counts)]
feature_expression[, fpkm_like := as.numeric(fpkm_like)]
feature_expression[, is_clean_cds := feature_type == "clean_CDS"]

sample_usage <- feature_expression[
  feature_type %chin% c("leader_uorf", "overlapping_uorf", "clean_CDS"),
  .(
    total_feature_fpkm = sum(fpkm_like, na.rm = TRUE),
    clean_cds_fpkm = sum(fpkm_like[is_clean_cds], na.rm = TRUE),
    total_feature_raw_counts = sum(raw_counts, na.rm = TRUE),
    clean_cds_raw_counts = sum(raw_counts[is_clean_cds], na.rm = TRUE),
    n_features_measured = uniqueN(feature_id)
  ),
  by = .(Run, gene_symbol, tx_id)
]
sample_usage[, clean_cds_relative_density := fifelse(
  total_feature_fpkm > 0,
  clean_cds_fpkm / total_feature_fpkm,
  NA_real_
)]
sample_usage[, clean_cds_raw_share := fifelse(
  total_feature_raw_counts > 0,
  clean_cds_raw_counts / total_feature_raw_counts,
  NA_real_
)]
sample_usage[, passes_clean_cds_count_cutoff := clean_cds_raw_counts >= min_clean_cds_raw_counts]
sample_usage <- merge(sample_usage, structure_summary, by = c("gene_symbol", "tx_id"), all.x = TRUE, sort = FALSE)

branch_usage[, clean_CDS := as.numeric(clean_CDS)]
branch_usage[, n_samples := as.numeric(n_samples)]
branch_usage[, clean_cds_raw_counts := as.numeric(clean_cds_raw_counts)]
branch_usage <- merge(branch_usage, structure_summary, by = c("gene_symbol", "tx_id"), all.x = TRUE, sort = FALSE)
branch_usage[, valid_clean_cds_usage :=
               has_modeled_uorf == TRUE &
               is.finite(clean_CDS) &
               passes_clean_cds_count_cutoff == TRUE]

sample_gene_summary <- sample_usage[
  has_modeled_uorf == TRUE &
    passes_clean_cds_count_cutoff == TRUE &
    is.finite(clean_cds_relative_density) &
    is.finite(clean_cds_raw_share),
  .(
    sample_rows = .N,
    sample_mean_clean_cds_relative_density = mean(clean_cds_relative_density),
    sample_median_clean_cds_relative_density = median(clean_cds_relative_density),
    sample_mean_clean_cds_raw_share = mean(clean_cds_raw_share),
    sample_median_clean_cds_raw_share = median(clean_cds_raw_share),
    sample_weighted_clean_cds_relative_density = weighted_mean_or_na(
      clean_cds_relative_density,
      sqrt(clean_cds_raw_counts + 1)
    )
  ),
  by = .(gene_symbol, tx_id)
]

branch_gene_summary <- branch_usage[
  valid_clean_cds_usage == TRUE,
  .(
    branch_rows = .N,
    branch_mean_clean_cds_usage = mean(clean_CDS),
    branch_median_clean_cds_usage = median(clean_CDS),
    branch_weighted_clean_cds_usage = weighted_mean_or_na(clean_CDS, sqrt(n_samples)),
    branch_sd_clean_cds_usage = sd(clean_CDS),
    branch_min_clean_cds_usage = min(clean_CDS),
    branch_max_clean_cds_usage = max(clean_CDS)
  ),
  by = .(gene_symbol, tx_id)
]

source_gene_summary <- branch_usage[
  valid_clean_cds_usage == TRUE,
  .(mean_clean_cds_usage = mean(clean_CDS)),
  by = .(gene_symbol, tx_id, branch_source)
]
source_gene_summary[, source_col := paste0("branch_source_", branch_source)]
source_gene_wide <- dcast(
  source_gene_summary,
  gene_symbol + tx_id ~ source_col,
  value.var = "mean_clean_cds_usage"
)

gene_rule_table <- Reduce(
  function(x, y) merge(x, y, by = c("gene_symbol", "tx_id"), all.x = TRUE, sort = FALSE),
  list(structure_summary, branch_gene_summary, sample_gene_summary, source_gene_wide)
)
numeric_gene_cols <- names(gene_rule_table)[vapply(gene_rule_table, is.numeric, logical(1))]
for (col in numeric_gene_cols) {
  set(gene_rule_table, i = which(!is.finite(gene_rule_table[[col]])), j = col, value = NA_real_)
}
gene_rule_table[, primary_clean_cds_usage := branch_weighted_clean_cds_usage]
gene_rule_table[is.na(primary_clean_cds_usage), primary_clean_cds_usage := sample_weighted_clean_cds_relative_density]
gene_rule_table[, primary_clean_cds_usage_logit := logit_clip(primary_clean_cds_usage)]
gene_rule_table[, clean_cds_usage_rank := frank(-primary_clean_cds_usage, ties.method = "average")]

count_prior_fit <- lm(
  primary_clean_cds_usage_logit ~ n_uorfs,
  data = gene_rule_table[is.finite(primary_clean_cds_usage_logit) & is.finite(n_uorfs)]
)
gene_rule_table[, expected_clean_cds_usage_from_uorf_count := plogis(
  as.numeric(predict(count_prior_fit, newdata = gene_rule_table))
)]
gene_rule_table[, clean_cds_delta_vs_count_expected := primary_clean_cds_usage - expected_clean_cds_usage_from_uorf_count]
gene_rule_table[, clean_cds_logit_residual_vs_count_expected := primary_clean_cds_usage_logit - as.numeric(predict(count_prior_fit, newdata = gene_rule_table))]
gene_rule_table[has_modeled_uorf != TRUE, `:=`(
  expected_clean_cds_usage_from_uorf_count = NA_real_,
  clean_cds_delta_vs_count_expected = NA_real_,
  clean_cds_logit_residual_vs_count_expected = NA_real_
)]
gene_rule_table[, hidden_uorf_candidate := has_modeled_uorf == TRUE &
                  !has_overlapping_uorf &
                  is.finite(clean_cds_delta_vs_count_expected) &
                  primary_clean_cds_usage <= 0.30 &
                  clean_cds_delta_vs_count_expected <= -0.15]
gene_rule_table[, hidden_uorf_candidate_priority := fifelse(
  hidden_uorf_candidate == TRUE &
    (primary_clean_cds_usage <= 0.15 | clean_cds_delta_vs_count_expected <= -0.25),
  "high",
  fifelse(hidden_uorf_candidate == TRUE, "medium", "not_flagged")
)]
gene_rule_table[, hidden_uorf_candidate_reason := fifelse(
  hidden_uorf_candidate == TRUE,
  paste0(
    "No predicted overlapping uoORF, but clean-CDS usage ",
    round(100 * primary_clean_cds_usage, 1),
    "% is far below the count-only expectation ",
    round(100 * expected_clean_cds_usage_from_uorf_count, 1),
    "%."
  ),
  ""
)]
gene_rule_table[, uorf_structure_flag := fifelse(
  has_overlapping_uorf,
  "predicted_uoORF",
  fifelse(hidden_uorf_candidate_priority == "high", "hidden_candidate_high",
          fifelse(hidden_uorf_candidate_priority == "medium", "hidden_candidate_medium", "not_flagged"))
)]

manual_note <- data.table(
  gene_symbol = c("LMNB1", "SNAI2", "MMP2"),
  structure_note = c(
    paste(
      "Manual M source adds the user-specified most-upstream active CTG",
      "overlapping uORF at transcript coordinates 313-795."
    ),
    paste(
      "Manual M source adds the strongest start-supported missing overlapping uORF",
      "from the ORFik findORFs plus categorize_and_filter_ORFs candidate space."
    ),
    paste(
      "Manual M source adds the strongest start-supported missing overlapping uORF",
      "from the ORFik findORFs plus categorize_and_filter_ORFs candidate space."
    )
  )
)
gene_rule_table <- merge(gene_rule_table, manual_note, by = "gene_symbol", all.x = TRUE, sort = FALSE)
gene_rule_table[is.na(structure_note), structure_note := ""]

predictor_specs <- list(
  n_uorfs = list(label = "Additional predicted uORF", unit = 1),
  n_leader_uorfs = list(label = "Additional leader-only uORF", unit = 1),
  n_overlapping_uorfs = list(label = "Additional overlapping uoORF", unit = 1),
  has_overlapping_uorf = list(label = "Has overlapping uoORF", unit = 1),
  n_overlapping_or_nested_pairs = list(label = "Additional nested/overlapping pair", unit = 1),
  n_close_pairs_10bp = list(label = "Additional pair <=10 bp apart", unit = 1),
  n_close_pairs_30bp = list(label = "Additional pair <=30 bp apart", unit = 1),
  min_inter_uorf_gap = list(label = "Min inter-uORF gap, per 30 bp", unit = 30),
  min_start_spacing = list(label = "Min start-to-start spacing, per 30 bp", unit = 30),
  n_leader_uorfs_stop_within_30bp_cds = list(label = "Leader uORFs ending <=30 bp before CDS", unit = 1),
  n_leader_uorfs_stop_within_60bp_cds = list(label = "Leader uORFs ending <=60 bp before CDS", unit = 1),
  n_leader_uorfs_stop_within_100bp_cds = list(label = "Leader uORFs ending <=100 bp before CDS", unit = 1),
  min_leader_stop_to_clean_cds_gap = list(label = "Closest leader uORF stop-to-CDS gap, per 30 bp", unit = 30),
  last_leader_stop_to_clean_cds_gap = list(label = "Last leader uORF stop-to-CDS gap, per 50 bp", unit = 50),
  first_uorf_to_cds_gap = list(label = "First uORF stop-to-CDS gap, per 100 bp", unit = 100),
  longest_uorf_bases = list(label = "Longest uORF length, per 30 bp", unit = 30),
  first_uorf_bases = list(label = "First uORF length, per 30 bp", unit = 30),
  uorf_union_fraction_of_leader = list(label = "Leader occupied by uORFs, per 10%", unit = 0.10),
  max_cds_overlap_fraction_of_clean_cds = list(label = "Max CDS-overlap fraction, per 10%", unit = 0.10)
)

fit_univariate <- function(predictor) {
  d <- gene_rule_table[
    is.finite(primary_clean_cds_usage_logit) &
      is.finite(get(predictor))
  ]
  if (nrow(d) < 6 || uniqueN(d[[predictor]]) < 2) return(NULL)
  unit <- predictor_specs[[predictor]]$unit
  d[, predictor_scaled := get(predictor) / unit]
  fit <- lm(primary_clean_cds_usage_logit ~ predictor_scaled, data = d)
  coefs <- summary(fit)$coefficients
  x_med <- median(d$predictor_scaled, na.rm = TRUE)
  pred0 <- plogis(unname(coef(fit)[1] + coef(fit)[2] * x_med))
  pred1 <- plogis(unname(coef(fit)[1] + coef(fit)[2] * (x_med + 1)))
  q <- quantile(d$predictor_scaled, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
  pred_q1 <- plogis(unname(coef(fit)[1] + coef(fit)[2] * q[1]))
  pred_q3 <- plogis(unname(coef(fit)[1] + coef(fit)[2] * q[2]))
  x_range <- range(d$predictor_scaled, na.rm = TRUE)
  pred_min <- plogis(unname(coef(fit)[1] + coef(fit)[2] * x_range[1]))
  pred_max <- plogis(unname(coef(fit)[1] + coef(fit)[2] * x_range[2]))
  effect_iqr <- 100 * (pred_q3 - pred_q1)
  effect_range <- 100 * (pred_max - pred_min)
  use_range <- isTRUE(all.equal(q[1], q[2])) || abs(effect_iqr) < 1e-9
  data.table(
    predictor = predictor,
    predictor_label = predictor_specs[[predictor]]$label,
    unit = unit,
    n_genes = nrow(d),
    coefficient_logit_per_unit = unname(coef(fit)[2]),
    p_value = coefs[2, 4],
    r_squared = summary(fit)$r.squared,
    adjusted_r_squared = summary(fit)$adj.r.squared,
    effect_per_unit_percentage_points = 100 * (pred1 - pred0),
    effect_iqr_percentage_points = effect_iqr,
    effect_range_percentage_points = effect_range,
    effect_display_percentage_points = ifelse(use_range, effect_range, effect_iqr),
    effect_display_basis = ifelse(use_range, "min-to-max", "IQR"),
    median_predictor_value = median(d[[predictor]], na.rm = TRUE),
    predictor_q25 = quantile(d[[predictor]], 0.25, na.rm = TRUE, names = FALSE),
    predictor_q75 = quantile(d[[predictor]], 0.75, na.rm = TRUE, names = FALSE),
    predictor_min = min(d[[predictor]], na.rm = TRUE),
    predictor_max = max(d[[predictor]], na.rm = TRUE)
  )
}

structure_rule_effects <- rbindlist(lapply(names(predictor_specs), fit_univariate), fill = TRUE)
if (nrow(structure_rule_effects) > 0) {
  structure_rule_effects[, q_value := p.adjust(p_value, method = "BH")]
  structure_rule_effects[, abs_effect_display_percentage_points := abs(effect_display_percentage_points)]
  setorder(structure_rule_effects, -abs_effect_display_percentage_points)
}

summarize_category <- function(category_col) {
  gene_rule_table[
    is.finite(primary_clean_cds_usage),
    .(
      n_genes = .N,
      genes = paste(gene_symbol[order(primary_clean_cds_usage)], collapse = ";"),
      mean_clean_cds_usage = mean(primary_clean_cds_usage),
      median_clean_cds_usage = median(primary_clean_cds_usage),
      sd_clean_cds_usage = sd(primary_clean_cds_usage),
      min_clean_cds_usage = min(primary_clean_cds_usage),
      max_clean_cds_usage = max(primary_clean_cds_usage)
    ),
    by = .(category = get(category_col))
  ][, category_variable := category_col][]
}

structure_category_effects <- rbindlist(
  lapply(
    c("uorf_count_bin", "overlap_class", "spacing_class", "last_leader_gap_class", "first_uorf_length_class"),
    summarize_category
  ),
  fill = TRUE
)
setcolorder(structure_category_effects, c("category_variable", "category"))
structure_category_effects[
  ,
  delta_vs_overall_percentage_points := 100 * (
    mean_clean_cds_usage - mean(gene_rule_table$primary_clean_cds_usage, na.rm = TRUE)
  )
]
setorder(structure_category_effects, category_variable, category)

model_formulas <- list(
  count_only = primary_clean_cds_usage_logit ~ n_uorfs,
  count_overlap = primary_clean_cds_usage_logit ~ n_uorfs + has_overlapping_uorf,
  count_close_pairs = primary_clean_cds_usage_logit ~ n_uorfs + n_close_pairs_30bp,
  count_overlap_gap = primary_clean_cds_usage_logit ~ n_uorfs + has_overlapping_uorf + last_leader_stop_to_clean_cds_gap,
  compact_rule = primary_clean_cds_usage_logit ~ n_uorfs + has_overlapping_uorf + n_close_pairs_30bp + longest_uorf_bases
)

loo_metrics <- function(formula, data) {
  response <- all.vars(formula)[1]
  predictors <- setdiff(all.vars(formula), response)
  cols <- c(response, predictors)
  d <- data[finite_complete_rows(data, cols)]
  if (nrow(d) < length(predictors) + 5) return(NULL)
  predicted <- rep(NA_real_, nrow(d))
  for (i in seq_len(nrow(d))) {
    train <- d[-i]
    test <- d[i]
    fit_i <- tryCatch(lm(formula, data = train), error = function(e) NULL)
    if (is.null(fit_i)) next
    predicted[i] <- tryCatch(as.numeric(predict(fit_i, newdata = test)), error = function(e) NA_real_)
  }
  ok <- is.finite(predicted)
  if (!any(ok)) return(NULL)
  observed <- d[[response]][ok]
  pred_prob <- plogis(predicted[ok])
  obs_prob <- plogis(observed)
  data.table(
    loo_n = sum(ok),
    loo_rmse_probability = sqrt(mean((pred_prob - obs_prob)^2)),
    loo_correlation = ifelse(length(unique(pred_prob)) > 1, cor(pred_prob, obs_prob), NA_real_)
  )
}

fit_model_summary <- function(name, formula) {
  response <- all.vars(formula)[1]
  predictors <- setdiff(all.vars(formula), response)
  cols <- c(response, predictors)
  d <- gene_rule_table[finite_complete_rows(gene_rule_table, cols)]
  if (nrow(d) < length(predictors) + 5) return(NULL)
  fit <- lm(formula, data = d)
  loo <- loo_metrics(formula, gene_rule_table)
  out <- data.table(
    model = name,
    formula = paste(deparse(formula), collapse = " "),
    n_genes = nrow(d),
    r_squared = summary(fit)$r.squared,
    adjusted_r_squared = summary(fit)$adj.r.squared,
    aic = AIC(fit),
    residual_sd_logit = sigma(fit)
  )
  if (!is.null(loo)) out <- cbind(out, loo)
  out
}

structure_rule_models <- rbindlist(
  Map(fit_model_summary, names(model_formulas), model_formulas),
  fill = TRUE
)
if (nrow(structure_rule_models) > 0) setorder(structure_rule_models, aic)

model_coefs <- rbindlist(lapply(names(model_formulas), function(name) {
  formula <- model_formulas[[name]]
  response <- all.vars(formula)[1]
  predictors <- setdiff(all.vars(formula), response)
  cols <- c(response, predictors)
  d <- gene_rule_table[finite_complete_rows(gene_rule_table, cols)]
  if (nrow(d) < length(predictors) + 5) return(NULL)
  fit <- lm(formula, data = d)
  out <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
  setnames(out, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
           c("estimate", "std_error", "t_value", "p_value"))
  out[, model := name]
  out
}), fill = TRUE)
if (nrow(model_coefs) > 0) {
  model_coefs[, q_value_within_model := p.adjust(p_value, method = "BH"), by = model]
  setcolorder(model_coefs, c("model", "term"))
}

branch_structure_rf <- branch_usage[
  valid_clean_cds_usage == TRUE,
  .(
    clean_CDS,
    n_uorfs,
    n_leader_uorfs,
    n_overlapping_uorfs,
    has_overlapping_uorf = as.factor(has_overlapping_uorf),
    n_close_pairs_10bp,
    n_close_pairs_30bp,
    min_inter_uorf_gap,
    last_leader_stop_to_clean_cds_gap,
    longest_uorf_bases,
    uorf_union_fraction_of_leader,
    max_cds_overlap_fraction_of_clean_cds,
    dominant_state = as.factor(dominant_state),
    branch_source = as.factor(branch_source)
  )
]
branch_structure_rf <- branch_structure_rf[complete.cases(branch_structure_rf)]
branch_rf_importance <- data.table()
branch_rf_performance <- data.table()
if (nrow(branch_structure_rf) >= 100) {
  set.seed(1)
  rf_fit <- randomForest(
    clean_CDS ~ .,
    data = branch_structure_rf,
    ntree = 500,
    importance = TRUE
  )
  branch_rf_importance <- as.data.table(importance(rf_fit), keep.rownames = "feature")
  setorder(branch_rf_importance, -`%IncMSE`)
  branch_rf_performance <- data.table(
    model = "branch_level_random_forest_structure_plus_branch_labels",
    n_rows = nrow(branch_structure_rf),
    n_unique_genes = uniqueN(branch_usage[valid_clean_cds_usage == TRUE, gene_symbol]),
    percent_variance_explained = tail(rf_fit$rsq, 1),
    mean_squared_residual = tail(rf_fit$mse, 1),
    warning = paste(
      "Descriptive only: branch rows repeat the same gene structures many times.",
      "Use gene-level effect tables for rule interpretation."
    )
  )
}

sample_rule_summary <- sample_usage[
  has_modeled_uorf == TRUE &
    passes_clean_cds_count_cutoff == TRUE &
    is.finite(clean_cds_relative_density),
  .(
    n_sample_gene_rows = .N,
    n_genes = uniqueN(gene_symbol),
    mean_clean_cds_relative_density = mean(clean_cds_relative_density),
    median_clean_cds_relative_density = median(clean_cds_relative_density)
  ),
  by = .(uorf_count_bin, overlap_class, spacing_class)
]

hidden_uorf_candidates <- gene_rule_table[
  has_modeled_uorf == TRUE & !has_overlapping_uorf & is.finite(primary_clean_cds_usage),
  .(
    gene_symbol,
    tx_id,
    hidden_uorf_candidate,
    hidden_uorf_candidate_priority,
    uorf_structure_flag,
    hidden_uorf_candidate_reason,
    n_uorfs,
    n_leader_uorfs,
    n_close_pairs_10bp,
    n_close_pairs_30bp,
    min_inter_uorf_gap,
    last_leader_stop_to_clean_cds_gap,
    primary_clean_cds_usage,
    expected_clean_cds_usage_from_uorf_count,
    clean_cds_delta_vs_count_expected,
    clean_cds_logit_residual_vs_count_expected,
    sample_weighted_clean_cds_relative_density,
    branch_rows,
    structure_note
  )
]
hidden_uorf_candidates[, priority_order := fifelse(
  hidden_uorf_candidate_priority == "high", 1L,
  fifelse(hidden_uorf_candidate_priority == "medium", 2L, 3L)
)]
setorder(hidden_uorf_candidates, priority_order, clean_cds_delta_vs_count_expected)
hidden_uorf_candidates[, priority_order := NULL]

rule_interpretation <- data.table(
  rule = c(
    "Additional predicted uORFs",
    "Overlapping uoORFs",
    "False-negative uoORF candidates",
    "Short spacing between uORFs",
    "Distance from the last leader uORF to clean CDS",
    "Long uORFs / leader occupancy",
    "Manual structure overrides"
  ),
  current_inference = c(
    "More predicted uORFs are associated with lower density-normalized clean-CDS usage across the 20 modeled genes.",
    "Genes with overlapping uoORFs, including manually added M-source overlaps, have low clean-CDS usage and strong negative allocation effects.",
    "Low clean-CDS usage in genes without predicted uoORFs should be treated as an annotation-quality signal. Remaining high-priority hidden-structure candidates need browser-level review before being promoted to M.",
    "Nested or very close uORF starts/stops track with lower clean-CDS usage, consistent with limited time/space for scanning complexes to reacquire initiation factors.",
    "A larger gap after the last leader uORF is not a strong simple rescue rule in this small panel; ATF4/MYC violate simple distance logic because overlapping uoORFs reset the clean-CDS coordinate.",
    "High leader occupancy and longer uORFs tend to reduce clean-CDS usage, but this is collinear with uORF count and overlap.",
    "Curated M-source geometry now includes several missing overlapping uoORFs, including MKI67, DDIT4, HSPA5, IGFBP5, SNAI2, MMP2, and LMNB1."
  ),
  caution = c(
    "n = 20 genes, so this is a triage rule rather than a genome-wide estimate.",
    "Only two predicted overlapping-uoORF genes are present, so the percentage effect is a candidate prior, not a calibrated constant.",
    "A low clean-CDS residual is not proof of a missed uoORF; it could also reflect RNA structure, poor uORF prediction, alternate transcript choice, or other translated leader features.",
    "Pair spacing is measured between predicted ORF intervals; overlapped/nested ORFs can make the minimum gap negative.",
    "The effective reinitiation distance probably depends on ORF length, termination context, Kozak strength, and RNA structure, not distance alone.",
    "Feature relative usage is density-normalized; short highly occupied uORFs can dominate the denominator.",
    "Manual browser review should be converted into explicit transcript-coordinate features before recalibrating global uORF rules."
  )
)

fwrite(structure_summary, file.path(output_dir, "uorf_structure_gene_geometry.csv"))
fwrite(sample_usage, file.path(output_dir, "sample_gene_clean_cds_usage_with_structure.csv"))
fwrite(branch_usage, file.path(output_dir, "branch_gene_clean_cds_usage_with_structure.csv"))
fwrite(gene_rule_table, file.path(output_dir, "uorf_structure_gene_rule_summary.csv"))
fwrite(structure_rule_effects, file.path(output_dir, "uorf_structure_univariate_rule_effects.csv"))
fwrite(structure_category_effects, file.path(output_dir, "uorf_structure_category_effects.csv"))
fwrite(structure_rule_models, file.path(output_dir, "uorf_structure_multivariable_rule_models.csv"))
fwrite(model_coefs, file.path(output_dir, "uorf_structure_multivariable_coefficients.csv"))
fwrite(branch_rf_importance, file.path(output_dir, "uorf_structure_branch_rf_importance.csv"))
fwrite(branch_rf_performance, file.path(output_dir, "uorf_structure_branch_rf_performance.csv"))
fwrite(sample_rule_summary, file.path(output_dir, "uorf_structure_sample_group_summary.csv"))
fwrite(hidden_uorf_candidates, file.path(output_dir, "uorf_structure_hidden_uorf_candidates.csv"))
fwrite(rule_interpretation, file.path(output_dir, "uorf_structure_rule_interpretation.csv"))

plot_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    legend.position = "right"
  )

label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggrepel::geom_text_repel(aes(label = gene_symbol), size = 3,
                           max.overlaps = Inf, show.legend = FALSE)
} else {
  geom_text(aes(label = gene_symbol), size = 2.8, vjust = -0.7,
            check_overlap = TRUE, show.legend = FALSE)
}

p_count <- ggplot(
  gene_rule_table[is.finite(primary_clean_cds_usage)],
  aes(n_uorfs, primary_clean_cds_usage, color = overlap_class, shape = spacing_class)
) +
  geom_point(aes(size = pmax(branch_rows, 1)), alpha = 0.9) +
  geom_smooth(
    data = gene_rule_table[is.finite(primary_clean_cds_usage)],
    aes(n_uorfs, primary_clean_cds_usage),
    inherit.aes = FALSE,
    method = "lm",
    se = TRUE,
    color = "grey25",
    linewidth = 0.7
  ) +
  label_layer +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  scale_x_continuous(breaks = sort(unique(gene_rule_table[is.finite(primary_clean_cds_usage), n_uorfs]))) +
  scale_size_continuous(name = "Branch rows", range = c(2.5, 7)) +
  labs(
    title = "Clean-CDS relative usage falls as predicted uORF count increases",
    subtitle = "Primary response is branch-weighted clean_CDS relative usage from RDG matrices; points are genes.",
    x = "Number of predicted uORFs/uoORFs",
    y = "Clean-CDS relative usage",
    color = NULL,
    shape = NULL
  ) +
  plot_theme
write_plot(p_count, "clean_cds_usage_vs_uorf_count", width = 10.5, height = 6.4)
write_plotly(p_count, "clean_cds_usage_vs_uorf_count")

p_hidden <- ggplot(
  gene_rule_table[is.finite(primary_clean_cds_usage) & is.finite(expected_clean_cds_usage_from_uorf_count)],
  aes(expected_clean_cds_usage_from_uorf_count, primary_clean_cds_usage)
) +
  geom_abline(slope = 1, intercept = 0, color = "grey45", linetype = "dashed") +
  geom_hline(yintercept = c(0.15, 0.30), color = "grey82", linetype = "dotted") +
  geom_point(aes(
    color = uorf_structure_flag,
    shape = overlap_class,
    size = pmax(n_uorfs, 1)
  ), alpha = 0.9) +
  label_layer +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  scale_color_manual(
    values = c(
      predicted_uoORF = "#762a83",
      hidden_candidate_high = "#b2182b",
      hidden_candidate_medium = "#ef8a62",
      not_flagged = "#2166ac"
    ),
    breaks = c("predicted_uoORF", "hidden_candidate_high", "hidden_candidate_medium", "not_flagged"),
    labels = c("predicted uoORF", "hidden candidate, high", "hidden candidate, medium", "not flagged"),
    name = "Structure flag"
  ) +
  scale_size_continuous(name = "Predicted uORFs", range = c(2.8, 7)) +
  labs(
    title = "Genes below their uORF-count expectation flag missing structure",
    subtitle = "Expected value comes from clean-CDS usage ~ predicted uORF count. No-overlap genes far below the diagonal are hidden uoORF candidates.",
    x = "Expected clean-CDS usage from predicted uORF count",
    y = "Observed clean-CDS relative usage",
    shape = NULL
  ) +
  plot_theme
write_plot(p_hidden, "hidden_uorf_candidate_residuals", width = 10.5, height = 6.4)
write_plotly(p_hidden, "hidden_uorf_candidate_residuals")

p_count_bin <- ggplot(
  gene_rule_table[is.finite(primary_clean_cds_usage)],
  aes(uorf_count_bin, primary_clean_cds_usage, color = overlap_class)
) +
  geom_boxplot(outlier.shape = NA, width = 0.55, color = "grey35", fill = "grey95") +
  geom_point(aes(size = pmax(branch_rows, 1)), position = position_jitter(width = 0.08, height = 0), alpha = 0.9) +
  geom_text(aes(label = gene_symbol), position = position_jitter(width = 0.08, height = 0), size = 2.7, vjust = -0.8, check_overlap = TRUE, show.legend = FALSE) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  scale_size_continuous(name = "Branch rows", range = c(2.5, 6)) +
  labs(
    title = "Count-bin rule: genes with fewer uORFs retain more clean-CDS usage",
    x = "Predicted uORF count bin",
    y = "Clean-CDS relative usage",
    color = NULL
  ) +
  plot_theme
write_plot(p_count_bin, "clean_cds_usage_by_uorf_count_bin", width = 9.5, height = 6)
write_plotly(p_count_bin, "clean_cds_usage_by_uorf_count_bin")

p_spacing <- ggplot(
  gene_rule_table[is.finite(primary_clean_cds_usage) & n_uorfs >= 2],
  aes(min_inter_uorf_gap, primary_clean_cds_usage, color = overlap_class)
) +
  geom_vline(xintercept = c(0, 10, 30, 60), linetype = c("solid", "dashed", "dashed", "dashed"), color = "grey65") +
  geom_point(aes(size = n_uorfs), alpha = 0.9) +
  label_layer +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  scale_size_continuous(name = "uORFs", range = c(3, 7)) +
  labs(
    title = "Inter-uORF spacing and clean-CDS usage",
    subtitle = "Gap = next uORF start minus previous uORF end minus 1; negative values mean nested/overlapping uORFs.",
    x = "Minimum inter-uORF gap (bp)",
    y = "Clean-CDS relative usage",
    color = NULL
  ) +
  plot_theme
write_plot(p_spacing, "clean_cds_usage_vs_min_inter_uorf_gap", width = 10.5, height = 6.2)
write_plotly(p_spacing, "clean_cds_usage_vs_min_inter_uorf_gap")

p_last_gap <- ggplot(
  gene_rule_table[is.finite(primary_clean_cds_usage) & is.finite(last_leader_stop_to_clean_cds_gap)],
  aes(last_leader_stop_to_clean_cds_gap, primary_clean_cds_usage, color = overlap_class)
) +
  geom_vline(xintercept = c(10, 50, 150), linetype = "dashed", color = "grey70") +
  geom_point(aes(size = n_uorfs), alpha = 0.9) +
  geom_smooth(
    data = gene_rule_table[is.finite(primary_clean_cds_usage) & is.finite(last_leader_stop_to_clean_cds_gap)],
    aes(last_leader_stop_to_clean_cds_gap, primary_clean_cds_usage),
    inherit.aes = FALSE,
    method = "lm",
    se = TRUE,
    color = "grey25",
    linewidth = 0.7
  ) +
  label_layer +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  scale_size_continuous(name = "uORFs", range = c(3, 7)) +
  labs(
    title = "Distance from last leader uORF to clean CDS is not enough by itself",
    subtitle = "Overlap class and uORF count dominate simple stop-to-CDS spacing in this small panel.",
    x = "Last leader uORF stop to clean-CDS start (bp)",
    y = "Clean-CDS relative usage",
    color = NULL
  ) +
  plot_theme
write_plot(p_last_gap, "clean_cds_usage_vs_last_leader_gap", width = 10.5, height = 6.2)
write_plotly(p_last_gap, "clean_cds_usage_vs_last_leader_gap")

p_overlap <- ggplot(
  gene_rule_table[is.finite(primary_clean_cds_usage)],
  aes(overlap_class, primary_clean_cds_usage, color = overlap_class)
) +
  geom_boxplot(outlier.shape = NA, width = 0.45, color = "grey35", fill = "grey95") +
  geom_point(aes(size = n_uorfs), position = position_jitter(width = 0.08, height = 0), alpha = 0.9) +
  geom_text(aes(label = gene_symbol), position = position_jitter(width = 0.08, height = 0), size = 2.7, vjust = -0.8, check_overlap = TRUE, show.legend = FALSE) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  scale_size_continuous(name = "uORFs", range = c(2.5, 7)) +
  labs(
    title = "Overlapping uoORFs are the strongest low-CDS geometry class",
    subtitle = "ATF4/MYC are predicted benchmarks; curated M-source uoORFs add missing overlaps for MKI67, DDIT4, HSPA5, IGFBP5, SNAI2, MMP2, and LMNB1.",
    x = NULL,
    y = "Clean-CDS relative usage",
    color = NULL
  ) +
  plot_theme +
  theme(legend.position = "none")
write_plot(p_overlap, "clean_cds_usage_by_overlap_class", width = 8, height = 5.8)
write_plotly(p_overlap, "clean_cds_usage_by_overlap_class")

effect_plot_data <- structure_rule_effects[
  is.finite(effect_display_percentage_points)
][
  order(effect_display_percentage_points)
]
effect_plot_data[, predictor_label := factor(predictor_label, levels = predictor_label)]
p_effects <- ggplot(
  effect_plot_data,
  aes(effect_display_percentage_points, predictor_label, fill = effect_display_percentage_points < 0)
) +
  geom_col(width = 0.65) +
  geom_vline(xintercept = 0, color = "grey35") +
  scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"), guide = "none") +
  labs(
    title = "Descriptive rule effects on clean-CDS usage",
    subtitle = "Bars show IQR effect; rare/binary predictors use min-to-max effect. Gene-level summaries, n = 20.",
    x = "Effect on clean-CDS usage (percentage points)",
    y = NULL
  ) +
  plot_theme
write_plot(p_effects, "uorf_structure_univariate_effects", width = 10, height = 6.8)

heat_features <- c(
  "primary_clean_cds_usage",
  "n_uorfs",
  "n_uorfs_with_tc_support",
  "n_uorfs_with_m_support",
  "n_overlapping_uorfs",
  "n_close_pairs_10bp",
  "n_close_pairs_30bp",
  "n_leader_uorfs_stop_within_30bp_cds",
  "min_leader_stop_to_clean_cds_gap",
  "min_inter_uorf_gap",
  "last_leader_stop_to_clean_cds_gap",
  "longest_uorf_bases",
  "uorf_union_fraction_of_leader"
)
heat_labels <- c(
  primary_clean_cds_usage = "clean_CDS usage",
  n_uorfs = "uORF count",
  n_uorfs_with_tc_support = "uORFs with TC support",
  n_uorfs_with_m_support = "uORFs with M support",
  n_overlapping_uorfs = "uoORF count",
  n_close_pairs_10bp = "pairs <=10 bp",
  n_close_pairs_30bp = "pairs <=30 bp",
  n_leader_uorfs_stop_within_30bp_cds = "stops <=30 bp to CDS",
  min_leader_stop_to_clean_cds_gap = "closest stop-to-CDS",
  min_inter_uorf_gap = "min gap",
  last_leader_stop_to_clean_cds_gap = "last gap to CDS",
  longest_uorf_bases = "longest uORF",
  uorf_union_fraction_of_leader = "leader occupied"
)
heat_input <- copy(gene_rule_table[is.finite(primary_clean_cds_usage)][order(primary_clean_cds_usage)])
heat_input[, (heat_features) := lapply(.SD, as.numeric), .SDcols = heat_features]
heat <- melt(
  heat_input,
  id.vars = c("gene_symbol", "primary_clean_cds_usage"),
  measure.vars = heat_features,
  variable.name = "feature",
  value.name = "value"
)
heat[, value_scaled := {
  if (all(is.na(value)) || sd(value, na.rm = TRUE) == 0) {
    rep(0, .N)
  } else {
    as.numeric(scale(value))
  }
}, by = feature]
heat[, feature_label := factor(heat_labels[as.character(feature)], levels = heat_labels[heat_features])]
heat[, gene_symbol := factor(gene_symbol, levels = unique(gene_rule_table[is.finite(primary_clean_cds_usage)][order(primary_clean_cds_usage), gene_symbol]))]
p_heat <- ggplot(heat, aes(feature_label, gene_symbol, fill = value_scaled)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "z") +
  labs(
    title = "Gene-level uORF structure rule matrix",
    subtitle = "Rows ordered by clean-CDS relative usage; red means high relative to this 20-gene panel.",
    x = NULL,
    y = NULL
  ) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
write_plot(p_heat, "uorf_structure_rule_heatmap", width = 10.5, height = 7.2)
write_plotly(p_heat, "uorf_structure_rule_heatmap")

timeline <- copy(annotation)
timeline[, feature_class := fifelse(
  feature_type == "clean_CDS", "clean CDS",
  fifelse(feature_type == "overlapping_uorf", "overlapping uoORF", "leader uORF")
)]
timeline <- merge(
  timeline,
  gene_rule_table[, .(gene_symbol, tx_id, primary_clean_cds_usage, n_uorfs, structure_note)],
  by = c("gene_symbol", "tx_id"),
  all.x = TRUE,
  sort = FALSE
)
timeline[, gene_label := paste0(gene_symbol, "  CDS=", round(100 * primary_clean_cds_usage), "%")]
timeline[, gene_label := factor(gene_label, levels = unique(timeline[order(primary_clean_cds_usage), gene_label]))]
p_timeline <- ggplot(timeline, aes(xmin = tx_start, xmax = tx_end, ymin = -0.25, ymax = 0.25, fill = feature_class)) +
  geom_rect(color = NA) +
  facet_grid(gene_label ~ ., scales = "free_x", space = "free_y") +
  scale_fill_manual(values = c("leader uORF" = "#fdae61", "overlapping uoORF" = "#b2182b", "clean CDS" = "#2166ac")) +
  labs(
    title = "Modeled transcript feature geometry",
    subtitle = "Rectangles have no outlines so very short uORFs remain visible; curated M-source uoORFs are included where available.",
    x = "Transcript coordinate in modeled leader+CDS region",
    y = NULL,
    fill = NULL
  ) +
  plot_theme +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey95", color = NA)
  )
write_plot(p_timeline, "uorf_structure_transcript_timeline", width = 11, height = 12)

message("Wrote:")
message("  ", file.path(output_dir, "uorf_structure_gene_rule_summary.csv"))
message("  ", file.path(output_dir, "uorf_structure_univariate_rule_effects.csv"))
message("  ", file.path(output_dir, "uorf_structure_category_effects.csv"))
message("  ", file.path(output_dir, "uorf_structure_multivariable_rule_models.csv"))
message("  ", file.path(output_dir, "figures"))
