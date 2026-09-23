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
results_dir <- file.path(analysis_dir, "results")
output_dir <- file.path(
  results_dir, "dominant_rdg_eif4g2_evidence_packet"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

specification <- data.table(
  comparison = c(
    "tRNA-Glu(UUC) OE", "tRNA-Arg(CCG) OE\nnegative control",
    "MDA-LM2\noriginal", "MDA-LM2\nindependent 2023",
    "CN34-LM1a\nvs parental",
    "UV", "ABCE1 KO"
  ),
  short_label = c(
    "tRNA-Glu OE", "tRNA-Arg OE", "MDA-LM2 original",
    "MDA-LM2 independent", "CN34-LM1a", "UV", "ABCE1 KO"
  ),
  evidence_role = c(
    "causal_tRNA_perturbation", "matched_tRNA_specificity_control",
    "metastatic_state_original", "metastatic_state_nonreplication",
    "metastatic_state_supporting",
    "collision_stress_recurrence", "recycling_loss_recurrence"
  ),
  source_dir = c(
    "dominant_rdg_eif4g2_trna_glu_lane_collapsed",
    "dominant_rdg_eif4g2_trna_arg_negative_control",
    "dominant_rdg_eif4g2_mda_lm2_validation",
    "dominant_rdg_eif4g2_mda_lm2_independent_validation",
    "dominant_rdg_eif4g2_cn34_lm1a_validation",
    "dominant_rdg_eif4g2_recycling_stress_evidence",
    "dominant_rdg_eif4g2_recycling_stress_evidence"
  ),
  study = c(
    "PRJNA310107-homo_sapiens", "PRJNA310096-homo_sapiens",
    "PRJNA310015-homo_sapiens", "PRJNA774895-homo_sapiens",
    "PRJNA309995-homo_sapiens",
    "PRJNA593553-homo_sapiens", "PRJNA602917-homo_sapiens"
  )
)
specification[, comparison_order := seq_len(.N)]

read_evidence_file <- function(file_name) {
  rbindlist(lapply(seq_len(nrow(specification)), function(i) {
    spec <- specification[i]
    path <- file.path(results_dir, spec$source_dir, file_name)
    if (!file.exists(path)) stop("Missing evidence input: ", path,
                                call. = FALSE)
    x <- fread(path)
    if ("study" %in% names(x)) x <- x[study == spec$study]
    if (!nrow(x)) stop("No matching study row in: ", path, call. = FALSE)
    x[, `:=`(
      comparison = spec$comparison,
      short_label = spec$short_label,
      evidence_role = spec$evidence_role,
      comparison_order = spec$comparison_order
    )]
    x
  }), fill = TRUE)
}

summary_dt <- read_evidence_file("union_exact_contrast_summary.csv")
unit_dt <- read_evidence_file("union_biological_unit_metrics.csv")
feature_dt <- read_evidence_file("uorf_feature_biological_unit_metrics.csv")
leaveout_dt <- read_evidence_file("union_position_leaveout_sensitivity.csv")

summary_columns <- c(
  "comparison_order", "comparison", "short_label", "evidence_role",
  "exact_union_validation_class", "study", "case_condition",
  "control_conditions", "n_case_runs", "n_control_runs",
  "n_case_biological_units", "n_control_biological_units",
  "case_union_uorf_counts", "control_union_uorf_counts",
  "case_unique_clean_cds_counts", "control_unique_clean_cds_counts",
  "union_clean_cds_allocation_delta", "uorf_cpm_log2_ratio",
  "clean_cds_cpm_log2_ratio", "biological_unit_range_separation",
  "case_uorf_phase0_fraction", "control_uorf_phase0_fraction",
  "frame0_mean_delta", "coordinate_shift_min_abs_delta",
  "leave_one_uorf_position_min_abs_delta", "pooled_top3_uorf_positions",
  "pooled_top3_removed_allocation_delta", "peak_robustness_gate"
)
summary_columns <- intersect(summary_columns, names(summary_dt))
setorder(summary_dt, comparison_order)
fwrite(
  summary_dt[, ..summary_columns],
  file.path(output_dir, "eif4g2_contrast_evidence.csv")
)

unit_columns <- c(
  "comparison_order", "comparison", "short_label", "evidence_role",
  "study", "biological_unit", "biological_unit_raw", "group", "run_count",
  "runs", "uorf_counts", "clean_cds_counts", "clean_cds_allocation",
  "uorf_phase0_fraction", "uorf_nonzero_position_fraction", "frame0",
  "top_readlength"
)
unit_columns <- intersect(unit_columns, names(unit_dt))
setorder(unit_dt, comparison_order, group, biological_unit)
fwrite(
  unit_dt[, ..unit_columns],
  file.path(output_dir, "eif4g2_biological_unit_evidence.csv")
)

uorf1_dt <- feature_dt[feature_id == "uORF1"]
feature_columns <- c(
  "comparison_order", "comparison", "short_label", "evidence_role",
  "study", "biological_unit", "group", "run_count", "runs", "feature_id",
  "tx_start", "tx_end", "start_codon", "kozak_strength", "counts",
  "density", "phase0_fraction", "start_9nt_fraction", "stop_6nt_fraction",
  "max_position_fraction", "nonzero_position_fraction"
)
feature_columns <- intersect(feature_columns, names(uorf1_dt))
setorder(uorf1_dt, comparison_order, group, biological_unit)
fwrite(
  uorf1_dt[, ..feature_columns],
  file.path(output_dir, "eif4g2_uorf1_biological_unit_evidence.csv")
)

leaveout_columns <- c(
  "comparison_order", "comparison", "short_label", "evidence_role",
  "study", "removal_type", "removed_positions", "n_removed_positions",
  "removed_pooled_counts", "clean_cds_allocation_delta"
)
leaveout_columns <- intersect(leaveout_columns, names(leaveout_dt))
setorder(leaveout_dt, comparison_order, removal_type, removed_positions)
fwrite(
  leaveout_dt[, ..leaveout_columns],
  file.path(output_dir, "eif4g2_position_leaveout_evidence.csv")
)

trna_length_dir <- file.path(
  results_dir, "dominant_rdg_eif4g2_trna_length_resolved"
)
trna_length_tables <- c(
  "run_manifest.csv",
  "calibration_design.csv",
  "psite_offset_qc.csv",
  "study_shared_well_phased_lengths.csv",
  "lane_shared_well_phased_lengths.csv",
  "lane_category_shared_length_coverage.csv",
  "eif4g2_unit_length_metrics.csv",
  "eif4g2_length_contrast_summary.csv",
  "eif4g2_length_specificity_summary.csv",
  "eif4g2_position_profile.csv",
  "eif4g2_coordinate_shift_sensitivity.csv",
  "eif4g2_position_leaveout_sensitivity.csv"
)
trna_length_files <- file.path(trna_length_dir, trna_length_tables)
if (!all(file.exists(trna_length_files))) {
  stop(
    "Run the raw tRNA length-resolved validator before the evidence packet: ",
    paste(trna_length_files[!file.exists(trna_length_files)], collapse = ", "),
    call. = FALSE
  )
}
for (i in seq_along(trna_length_files)) {
  fwrite(
    fread(trna_length_files[[i]]),
    file.path(output_dir, paste0("eif4g2_trna_", trna_length_tables[[i]]))
  )
}

uv_length_dir <- file.path(
  results_dir, "dominant_rdg_eif4g2_uv_length_resolved"
)
uv_length_files <- file.path(uv_length_dir, c(
  "eif4g2_length_contrast_summary.csv",
  "eif4g2_run_length_metrics.csv",
  "eif4g2_position_profile.csv",
  "eif4g2_position_leaveout_sensitivity.csv",
  "calibration_design.csv"
))
if (!all(file.exists(uv_length_files))) {
  stop(
    "Run the raw UV length-resolved validator before the evidence packet: ",
    paste(uv_length_files[!file.exists(uv_length_files)], collapse = ", "),
    call. = FALSE
  )
}
uv_length_summary <- fread(uv_length_files[[1L]])
uv_length_runs <- fread(uv_length_files[[2L]])
uv_position_profile <- fread(uv_length_files[[3L]])
uv_position_leaveout <- fread(uv_length_files[[4L]])
uv_calibration_design <- fread(uv_length_files[[5L]])
fwrite(
  uv_length_summary,
  file.path(output_dir, "eif4g2_uv_length_contrast_evidence.csv")
)
fwrite(
  uv_length_runs,
  file.path(output_dir, "eif4g2_uv_run_length_evidence.csv")
)
fwrite(
  uv_position_leaveout,
  file.path(output_dir, "eif4g2_uv_position_leaveout_evidence.csv")
)

plot_uv_length_falsification <- function(device) {
  device()
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(
    mfrow = c(2, 2), mar = c(4.7, 4.6, 2.5, 0.8),
    oma = c(0, 0, 1.8, 0), las = 1
  )
  passing_categories <- c(
    "all_shared_well_phased",
    "canonical_28_32nt_shared_well_phased",
    "short_20_23nt_shared_well_phased"
  )
  category_labels <- c("All shared", "Canonical\n28-32 nt", "Short\n20-23 nt")
  colours <- c(control = "#4C78A8", UV = "#E45756")

  run_plot <- uv_length_runs[category %chin% passing_categories]
  run_plot[, category_order := match(category, passing_categories)]
  plot(
    NA, xlim = c(0.5, 3.5), ylim = c(0.3, 1), xaxt = "n",
    xlab = "", ylab = "Clean-CDS allocation",
    main = "Replicate separation within length strata"
  )
  abline(h = seq(0.4, 1, by = 0.1), col = "grey92")
  axis(1, at = 1:3, labels = category_labels)
  for (condition_value in c("control", "UV")) {
    x <- run_plot[condition == condition_value]
    offset <- if (condition_value == "control") -0.12 else 0.12
    points(
      x$category_order + offset,
      x$clean_cds_allocation,
      pch = ifelse(x$replicate == 1L, 21, 24),
      bg = colours[[condition_value]], col = "white", cex = 1.25
    )
  }
  legend(
    "bottomleft", legend = c("Control", "UV"), pch = 21,
    pt.bg = colours, col = "white", bty = "n", horiz = TRUE, cex = 0.8
  )

  exact_lengths <- uv_length_summary[grepl("^[0-9]+nt$", category)]
  exact_lengths[, read_length := as.integer(sub("nt$", "", category))]
  shared_lengths <- as.integer(strsplit(
    uv_calibration_design$shared_well_phased_lengths[[1L]], ";",
    fixed = TRUE
  )[[1L]])
  plot(
    NA, xlim = c(14.5, 34.5), ylim = c(-0.5, 0.1),
    xlab = "Footprint length (nt)",
    ylab = "UV - control clean-CDS allocation",
    main = "Exact-length direction"
  )
  rect(19.5, -1, 23.5, 1, col = rgb(0.90, 0.45, 0.35, 0.10), border = NA)
  rect(27.5, -1, 32.5, 1, col = rgb(0.25, 0.45, 0.75, 0.10), border = NA)
  abline(h = 0, lty = 2, col = "grey45")
  points(
    exact_lengths$read_length,
    exact_lengths$clean_cds_allocation_delta,
    pch = ifelse(exact_lengths$read_length %in% shared_lengths, 21, 1),
    bg = ifelse(
      exact_lengths$read_length %in% shared_lengths, "#6A51A3", "white"
    ),
    col = "#6A51A3", cex = 1.15
  )
  legend(
    "bottomright", legend = c("Shared well-phased", "Not shared"),
    pch = c(21, 1), pt.bg = c("#6A51A3", "white"),
    col = "#6A51A3", bty = "n", cex = 0.72
  )

  shared_profile <- uv_position_profile[
    shared_well_phased == TRUE &
      feature %chin% c("uORF1", "uORF2", "clean_CDS")
  ]
  profile_totals <- shared_profile[, .(target_total = sum(N)), by = condition]
  leader_profile <- shared_profile[tx_position <= 400L,
    .(counts = sum(N)), by = .(condition, tx_position)
  ]
  leader_grid <- CJ(
    condition = c("control", "UV"), tx_position = 1:400,
    unique = TRUE
  )
  leader_profile <- merge(
    leader_grid, leader_profile,
    by = c("condition", "tx_position"), all.x = TRUE, sort = FALSE
  )
  leader_profile[is.na(counts), counts := 0]
  leader_profile <- merge(
    leader_profile, profile_totals, by = "condition", all.x = TRUE,
    sort = FALSE
  )
  leader_profile[, target_per_thousand := counts / target_total * 1000]
  y_max <- max(leader_profile$target_per_thousand) * 1.04
  plot(
    NA, xlim = c(1, 400), ylim = c(0, y_max),
    xlab = "EIF4G2 transcript position",
    ylab = "P-sites per 1,000 target P-sites",
    main = "Shared-length position profile"
  )
  abline(v = c(23, 73, 128, 220, 309), col = "grey88", lty = 3)
  for (condition_value in c("control", "UV")) {
    x <- leader_profile[condition == condition_value]
    lines(
      x$tx_position, x$target_per_thousand,
      type = "h", col = grDevices::adjustcolor(
        colours[[condition_value]], alpha.f = 0.62
      ), lwd = 1.1
    )
  }
  text(c(48, 174, 340), y_max * 0.96, c("uORF1", "uORF2", "CDS"),
       cex = 0.72)
  legend(
    "topright", legend = c("Control", "UV"),
    col = colours, lwd = 2, bty = "n", cex = 0.78
  )

  robustness <- uv_length_summary[
    match(passing_categories, category),
    .(
      category,
      observed = clean_cds_allocation_delta,
      worst_coordinate = sign(clean_cds_allocation_delta) *
        coordinate_shift_min_abs_delta,
      worst_single = sign(clean_cds_allocation_delta) *
        single_position_min_abs_delta,
      worst_three = greedy_final_clean_cds_allocation_delta
    )
  ]
  mat <- t(as.matrix(robustness[, .(
    observed, worst_coordinate, worst_single, worst_three
  )]))
  colnames(mat) <- category_labels
  barplot(
    mat, beside = TRUE, ylim = c(-0.4, 0),
    col = c("#0072B2", "#CC79A7", "#E69F00", "#009E73"),
    border = NA, las = 1, ylab = "Clean-CDS allocation delta",
    main = "Worst-case robustness", cex.names = 0.78
  )
  abline(h = 0, col = "grey45")
  legend(
    "bottomright",
    legend = c("Observed", "Worst +/-2 nt", "Worst one position",
               "Worst three positions"),
    fill = c("#0072B2", "#CC79A7", "#E69F00", "#009E73"),
    border = NA, bty = "n", cex = 0.68
  )
  mtext(
    "Raw EIF4G2 UV falsification: the switch survives length, phase, replicate, and peak gates",
    outer = TRUE, side = 3, line = 0.35, font = 2
  )
}

plot_uv_length_falsification(function() {
  png(
    file.path(output_dir, "eif4g2_uv_length_falsification.png"),
    width = 2600, height = 1900, res = 220
  )
})
plot_uv_length_falsification(function() {
  pdf(
    file.path(output_dir, "eif4g2_uv_length_falsification.pdf"),
    width = 12, height = 8.6, useDingbats = FALSE
  )
})

plot_evidence <- function(device) {
  device()
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(
    mfrow = c(1, 2), mar = c(8.2, 4.6, 2.4, 0.8),
    oma = c(0, 0, 1.6, 0), las = 1
  )

  n_comparisons <- nrow(specification)
  plot(
    NA, xlim = c(0.5, n_comparisons + 0.5), ylim = c(0.45, 1),
    xaxt = "n", xlab = "", ylab = "Clean-CDS allocation",
    main = "Biological units"
  )
  abline(h = seq(0.5, 1, by = 0.1), col = "grey92", lwd = 0.8)
  axis(
    1, at = seq_len(n_comparisons), labels = specification$comparison,
    las = 2, cex.axis = 0.72
  )
  group_colours <- c(control = "#666666", case = "#D55E00")
  for (i in seq_len(n_comparisons)) {
    x <- unit_dt[comparison_order == i]
    for (group_name in c("control", "case")) {
      values <- x[group == group_name, clean_cds_allocation]
      offset <- if (group_name == "control") -0.13 else 0.13
      jitter <- if (length(values) > 1L) {
        seq(-0.025, 0.025, length.out = length(values))
      } else 0
      points(
        i + offset + jitter, values, pch = 21,
        bg = group_colours[[group_name]], col = "white", cex = 1.15
      )
      segments(
        i + offset - 0.08, mean(values), i + offset + 0.08, mean(values),
        col = group_colours[[group_name]], lwd = 2.2
      )
    }
  }
  legend(
    "bottomleft", inset = 0.01, legend = c("Control", "Case"), pch = 21,
    pt.bg = group_colours[c("control", "case")], col = "white", bty = "n",
    horiz = TRUE, cex = 0.8
  )

  robustness <- summary_dt[, .(
    comparison_order,
    observed = union_clean_cds_allocation_delta,
    worst_single = sign(union_clean_cds_allocation_delta) *
      leave_one_uorf_position_min_abs_delta,
    top3_removed = pooled_top3_removed_allocation_delta
  )]
  y <- rev(seq_len(n_comparisons))
  plot(
    NA, xlim = c(-0.36, 0.08), ylim = c(0.5, n_comparisons + 0.5),
    yaxt = "n", xlab = "Case - control clean-CDS allocation",
    ylab = "", main = "Peak-removal robustness"
  )
  abline(v = 0, col = "grey45", lty = 2)
  abline(v = seq(-0.3, 0, by = 0.1), col = "grey92", lwd = 0.8)
  axis(
    2, at = y, labels = specification$comparison,
    las = 1, cex.axis = 0.72
  )
  for (i in seq_len(n_comparisons)) {
    row <- robustness[comparison_order == i]
    points(row$observed, y[[i]], pch = 16, col = "#0072B2", cex = 1.05)
    points(row$worst_single, y[[i]], pch = 17, col = "#E69F00", cex = 1.0)
    points(row$top3_removed, y[[i]], pch = 15, col = "#009E73", cex = 0.9)
  }
  legend(
    "bottomleft", inset = 0.01,
    legend = c("Observed", "Worst single deletion", "Pooled top 3 removed"),
    pch = c(16, 17, 15), col = c("#0072B2", "#E69F00", "#009E73"),
    bty = "n", cex = 0.72
  )
  mtext(
    "EIF4G2 uORF-to-CDS redistribution: specificity, recurrence, and falsification",
    outer = TRUE, side = 3, line = 0.2, font = 2
  )
}

plot_evidence(function() {
  png(
    file.path(output_dir, "eif4g2_evidence_overview.png"),
    width = 2600, height = 1450, res = 220
  )
})
plot_evidence(function() {
  pdf(
    file.path(output_dir, "eif4g2_evidence_overview.pdf"),
    width = 12.2, height = 6.8, useDingbats = FALSE
  )
})

message("Saved EIF4G2 evidence packet: ", output_dir)
