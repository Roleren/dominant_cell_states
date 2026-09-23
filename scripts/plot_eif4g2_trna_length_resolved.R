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
input_dir <- file.path(
  analysis_dir, "results", "dominant_rdg_eif4g2_trna_length_resolved"
)
required <- file.path(input_dir, c(
  "eif4g2_unit_length_metrics.csv",
  "eif4g2_length_contrast_summary.csv",
  "eif4g2_length_specificity_summary.csv",
  "lane_shared_well_phased_lengths.csv"
))
if (!all(file.exists(required))) {
  stop(
    "Run dominant_rdg_eif4g2_trna_length_resolved.R first. Missing: ",
    paste(basename(required[!file.exists(required)]), collapse = ", "),
    call. = FALSE
  )
}

unit_metrics <- fread(required[[1L]])
contrast <- fread(required[[2L]])
specificity <- fread(required[[3L]])
lane_lengths <- fread(required[[4L]])
main_category <- "all_lane_matched_well_phased"

main_units <- unit_metrics[category == main_category]
if (!all(c("tRNA-Glu", "tRNA-Arg") %in% main_units$perturbation)) {
  stop("The lane-matched category is incomplete.", call. = FALSE)
}

exact <- contrast[grepl("^[0-9]+nt$", category)]
exact[, read_length := as.integer(sub("nt$", "", category))]
exact_plot <- exact[read_length >= 18L & read_length <= 32L]
lanes <- contrast[
  perturbation == "tRNA-Glu" &
    category %chin% paste0("lane", 1:3, "_shared_well_phased")
]
lanes[, lane := as.integer(sub("lane([0-9]+).*", "\\1", category))]
robust <- contrast[
  perturbation == "tRNA-Glu" & category == main_category
]
if (nrow(robust) != 1L) stop("Missing main tRNA-Glu robustness row.")

source_data <- rbindlist(list(
  main_units[, .(
    panel = "biological_units", perturbation, condition, replicate,
    category, clean_cds_allocation,
    clean_cds_vs_uorf1_allocation,
    uorf1_phase0_fraction = uorf_phase0_fraction
  )],
  exact[, .(
    panel = "exact_lengths", perturbation, category, read_length,
    clean_cds_allocation_delta,
    clean_cds_vs_uorf1_allocation_delta,
    all_runs_phase_gate, measurement_gate, evidence_gate
  )],
  lanes[, .(
    panel = "matched_lanes", perturbation, category, lane,
    clean_cds_allocation_delta,
    clean_cds_vs_uorf1_allocation_delta,
    case_uorf_phase0_fraction, control_uorf_phase0_fraction
  )],
  robust[, .(
    panel = "robustness", perturbation, category,
    clean_cds_allocation_delta,
    coordinate_shift_min_abs_delta,
    single_position_min_abs_delta,
    greedy_final_clean_cds_allocation_delta,
    evidence_gate
  )],
  specificity[category == main_category]
    [, c(list(panel = "specificity"), as.list(.SD))]
), fill = TRUE, use.names = TRUE)
fwrite(
  source_data,
  file.path(input_dir, "eif4g2_trna_length_figure_source_data.csv")
)

plot_falsification <- function(device) {
  device()
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(
    mfrow = c(2, 2), mar = c(4.8, 4.7, 2.7, 0.8),
    oma = c(0, 0, 1.8, 0), las = 1
  )
  colours <- c(
    "tRNA-Glu" = "#0072B2", "tRNA-Arg" = "#E69F00",
    control = "#4D4D4D", case = "#D55E00"
  )

  group_order <- data.table(
    perturbation = rep(c("tRNA-Glu", "tRNA-Arg"), each = 2L),
    condition = rep(c("control", "case"), 2L),
    x = 1:4,
    label = c("Glu\ncontrol", "Glu\nOE", "Arg\ncontrol", "Arg\nOE")
  )
  unit_plot <- merge(
    main_units, group_order,
    by = c("perturbation", "condition"), all.x = TRUE, sort = FALSE
  )
  plot(
    NA, xlim = c(0.5, 4.5), ylim = range(
      c(unit_plot$clean_cds_allocation,
        unit_plot$clean_cds_vs_uorf1_allocation),
      finite = TRUE
    ) + c(-0.025, 0.025),
    xaxt = "n", xlab = "", ylab = "Clean-CDS allocation",
    main = "True biological-unit separation"
  )
  abline(h = seq(0.5, 1, by = 0.05), col = "grey93")
  axis(1, at = 1:4, labels = group_order$label)
  jitter <- ifelse(unit_plot$replicate == 1L, -0.06, 0.06)
  points(
    unit_plot$x + jitter, unit_plot$clean_cds_allocation,
    pch = 21, bg = colours[unit_plot$condition], col = "white", cex = 1.3
  )
  points(
    unit_plot$x + jitter, unit_plot$clean_cds_vs_uorf1_allocation,
    pch = 2, col = colours[unit_plot$condition], cex = 0.9
  )
  legend(
    "bottomright", legend = c("uORF union", "exact uORF1"),
    pch = c(21, 2), pt.bg = c("grey55", NA), col = "grey35",
    bty = "n", cex = 0.75
  )

  plot(
    NA, xlim = range(exact_plot$read_length, finite = TRUE) + c(-0.5, 0.5),
    ylim = range(exact_plot$clean_cds_allocation_delta, finite = TRUE) +
      c(-0.025, 0.025),
    xlab = "Footprint length (nt)",
    ylab = "Case - control allocation",
    main = "18-32-nt exact-length direction"
  )
  abline(h = 0, lty = 2, col = "grey45")
  for (perturbation_value in c("tRNA-Glu", "tRNA-Arg")) {
    z <- exact_plot[perturbation == perturbation_value][order(read_length)]
    lines(
      z$read_length, z$clean_cds_allocation_delta,
      col = grDevices::adjustcolor(colours[[perturbation_value]], 0.55)
    )
    points(
      z$read_length, z$clean_cds_allocation_delta,
      pch = ifelse(z$all_runs_phase_gate, 21, 1),
      bg = ifelse(z$all_runs_phase_gate,
                  colours[[perturbation_value]], "white"),
      col = colours[[perturbation_value]], cex = 1.05
    )
  }
  legend(
    "bottomright", legend = c("tRNA-Glu", "tRNA-Arg"),
    col = colours[c("tRNA-Glu", "tRNA-Arg")], pch = 21,
    pt.bg = colours[c("tRNA-Glu", "tRNA-Arg")], bty = "n", cex = 0.75
  )

  setorder(lanes, lane)
  plot(
    NA, xlim = c(0.5, 3.5),
    ylim = range(c(
      lanes$clean_cds_allocation_delta,
      lanes$clean_cds_vs_uorf1_allocation_delta, 0
    ), finite = TRUE) + c(-0.02, 0.02),
    xaxt = "n", xlab = "Matched technical lane",
    ylab = "tRNA-Glu - control allocation",
    main = "Independent lane consistency"
  )
  axis(1, at = 1:3, labels = paste("Lane", 1:3))
  abline(h = 0, lty = 2, col = "grey45")
  points(
    lanes$lane - 0.05, lanes$clean_cds_allocation_delta,
    pch = 21, bg = colours[["tRNA-Glu"]], col = "white", cex = 1.35
  )
  points(
    lanes$lane + 0.05, lanes$clean_cds_vs_uorf1_allocation_delta,
    pch = 2, col = colours[["tRNA-Glu"]], cex = 1.05
  )

  robustness_values <- c(
    Observed = robust$clean_cds_allocation_delta,
    `Worst +/-2 nt` = sign(robust$clean_cds_allocation_delta) *
      robust$coordinate_shift_min_abs_delta,
    `Worst one position` = sign(robust$clean_cds_allocation_delta) *
      robust$single_position_min_abs_delta,
    `Worst three positions` = robust$greedy_final_clean_cds_allocation_delta
  )
  barplot(
    robustness_values, ylim = range(c(robustness_values, 0)) + c(-0.02, 0),
    col = c("#0072B2", "#CC79A7", "#E69F00", "#009E73"),
    border = NA, las = 2, cex.names = 0.73,
    ylab = "Clean-CDS allocation delta",
    main = "Worst-case raw robustness"
  )
  abline(h = 0, col = "grey45")
  mtext(
    "Raw EIF4G2 tRNA falsification: exact-uORF1, length, lane and peak gates",
    outer = TRUE, side = 3, line = 0.35, font = 2
  )
}

plot_falsification(function() {
  png(
    file.path(input_dir, "eif4g2_trna_length_falsification.png"),
    width = 2600, height = 1900, res = 220
  )
})
plot_falsification(function() {
  pdf(
    file.path(input_dir, "eif4g2_trna_length_falsification.pdf"),
    width = 12, height = 8.6, useDingbats = FALSE
  )
})

message("Saved raw tRNA EIF4G2 figure and source data: ", input_dir)
