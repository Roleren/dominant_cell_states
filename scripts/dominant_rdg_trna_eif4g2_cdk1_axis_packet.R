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
  results_dir, "dominant_rdg_trna_eif4g2_cdk1_axis_packet"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

specification <- data.table(
  gene_symbol = c("EIF4G2", "EIF4G2", "CDK1", "CDK1"),
  perturbation = c("tRNA-Glu", "tRNA-Arg", "tRNA-Glu", "tRNA-Arg"),
  short_label = c(
    "EIF4G2\ntRNA-Glu", "EIF4G2\ntRNA-Arg",
    "CDK1\ntRNA-Glu", "CDK1\ntRNA-Arg"
  ),
  architecture = c(
    "exact non-overlapping ATG uORF1",
    "exact non-overlapping ATG uORF1",
    "six overlapping near-cognate/other uORFs",
    "six overlapping near-cognate/other uORFs"
  ),
  source_dir = c(
    "dominant_rdg_eif4g2_trna_glu_lane_collapsed",
    "dominant_rdg_eif4g2_trna_arg_negative_control",
    "dominant_rdg_cdk1_trna_glu_validation",
    "dominant_rdg_cdk1_trna_arg_validation"
  )
)
specification[, comparison_order := seq_len(.N)]

read_specified <- function(file_name) {
  rbindlist(lapply(seq_len(nrow(specification)), function(i) {
    spec <- specification[i]
    path <- file.path(results_dir, spec$source_dir, file_name)
    if (!file.exists(path)) stop("Missing axis input: ", path, call. = FALSE)
    x <- fread(path)
    x <- x[gene_symbol == spec$gene_symbol]
    if (!nrow(x)) stop("No target row in: ", path, call. = FALSE)
    x[, `:=`(
      comparison_order = spec$comparison_order,
      perturbation = spec$perturbation,
      short_label = spec$short_label,
      architecture = spec$architecture
    )]
    x
  }), fill = TRUE)
}

summary <- read_specified("union_exact_contrast_summary.csv")
units <- read_specified("union_biological_unit_metrics.csv")
setorder(summary, comparison_order)
setorder(units, comparison_order, group, biological_unit)

summary_columns <- c(
  "comparison_order", "short_label", "gene_symbol", "tx_id",
  "perturbation", "architecture", "exact_union_validation_class",
  "n_case_biological_units", "n_control_biological_units",
  "case_union_uorf_counts", "control_union_uorf_counts",
  "case_unique_clean_cds_counts", "control_unique_clean_cds_counts",
  "union_clean_cds_allocation_delta", "biological_unit_range_separation",
  "case_uorf_phase0_fraction", "control_uorf_phase0_fraction",
  "coordinate_shift_min_abs_delta",
  "leave_one_uorf_position_min_abs_delta",
  "pooled_top3_removed_allocation_delta", "peak_robustness_gate"
)
fwrite(
  summary[, ..summary_columns],
  file.path(output_dir, "trna_eif4g2_cdk1_axis_contrasts.csv")
)

unit_columns <- c(
  "comparison_order", "short_label", "gene_symbol", "tx_id",
  "perturbation", "architecture", "biological_unit", "group", "run_count",
  "runs", "uorf_counts", "clean_cds_counts", "clean_cds_allocation",
  "uorf_phase0_fraction", "frame0", "top_readlength"
)
unit_columns <- intersect(unit_columns, names(units))
fwrite(
  units[, ..unit_columns],
  file.path(output_dir, "trna_eif4g2_cdk1_axis_biological_units.csv")
)

plot_axis <- function(device) {
  device()
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(
    mfrow = c(1, 2), mar = c(7.2, 4.7, 2.5, 0.8),
    oma = c(0, 0, 1.7, 0), las = 1
  )

  plot(
    NA, xlim = c(0.5, 4.5), ylim = c(0.35, 0.92), xaxt = "n",
    xlab = "", ylab = "Clean-CDS allocation",
    main = "Biological units (n = 2 per arm)"
  )
  abline(h = seq(0.4, 0.9, by = 0.1), col = "grey92")
  axis(1, at = 1:4, labels = specification$short_label,
       las = 2, cex.axis = 0.78)
  colours <- c(control = "#666666", case = "#D55E00")
  for (i in 1:4) {
    x <- units[comparison_order == i]
    for (group_name in c("control", "case")) {
      values <- x[group == group_name, clean_cds_allocation]
      offset <- if (group_name == "control") -0.13 else 0.13
      points(
        i + offset + seq(-0.025, 0.025, length.out = length(values)),
        values, pch = 21, bg = colours[[group_name]], col = "white",
        cex = 1.2
      )
      segments(
        i + offset - 0.08, mean(values), i + offset + 0.08, mean(values),
        col = colours[[group_name]], lwd = 2.2
      )
    }
  }
  legend(
    "bottomleft", legend = c("Control", "Perturbation"), pch = 21,
    pt.bg = colours, col = "white", bty = "n", horiz = TRUE, cex = 0.8
  )

  robustness <- summary[, .(
    comparison_order,
    observed = union_clean_cds_allocation_delta,
    worst_single = sign(union_clean_cds_allocation_delta) *
      leave_one_uorf_position_min_abs_delta,
    top3_removed = pooled_top3_removed_allocation_delta
  )]
  y <- rev(1:4)
  plot(
    NA, xlim = c(-0.24, 0.08), ylim = c(0.5, 4.5), yaxt = "n",
    xlab = "Perturbation - control clean-CDS allocation", ylab = "",
    main = "Position-deletion falsification"
  )
  abline(v = 0, col = "grey45", lty = 2)
  abline(v = seq(-0.2, 0, by = 0.05), col = "grey92")
  axis(2, at = y, labels = specification$short_label,
       las = 1, cex.axis = 0.78)
  for (i in 1:4) {
    row <- robustness[comparison_order == i]
    points(row$observed, y[[i]], pch = 16, col = "#0072B2", cex = 1.1)
    points(row$worst_single, y[[i]], pch = 17, col = "#E69F00", cex = 1.0)
    points(row$top3_removed, y[[i]], pch = 15, col = "#009E73", cex = 0.95)
  }
  legend(
    "bottomleft",
    legend = c("Observed", "Worst single deletion", "Top 3 removed"),
    pch = c(16, 17, 15), col = c("#0072B2", "#E69F00", "#009E73"),
    bty = "n", cex = 0.72
  )
  mtext(
    "tRNA-Glu-selective EIF4G2-CDK1 leader axis (CDK1 architecture remains ambiguous)",
    outer = TRUE, side = 3, line = 0.2, font = 2
  )
}

plot_axis(function() {
  png(
    file.path(output_dir, "trna_eif4g2_cdk1_axis.png"),
    width = 2400, height = 1350, res = 220
  )
})
plot_axis(function() {
  pdf(
    file.path(output_dir, "trna_eif4g2_cdk1_axis.pdf"),
    width = 11.2, height = 6.2, useDingbats = FALSE
  )
})

message("Saved tRNA EIF4G2-CDK1 axis packet: ", output_dir)
