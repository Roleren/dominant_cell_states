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
state_file <- file.path(analysis_dir, "human_dominant_cell_states_clean_cds.csv")
expression_file <- file.path(
  analysis_dir,
  "human_dominant_cell_states_clean_cds_gene_expression.csv"
)
output_dir <- file.path(analysis_dir, "dominant_state_score_audit")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

module_definition_file <- file.path(
  analysis_dir,
  "scripts",
  "dominant_state_module_definitions.R"
)
if (!file.exists(module_definition_file)) {
  stop("Missing dominant-state module definition file: ",
       module_definition_file, call. = FALSE)
}
source(module_definition_file)
modules <- load_dominant_state_modules(analysis_dir)

htmlwidget_helper <- c(
  file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R"),
  file.path("scripts", "dominant_htmlwidgets.R"),
  "dominant_htmlwidgets.R"
)
htmlwidget_helper <- htmlwidget_helper[file.exists(htmlwidget_helper)][1]
if (!is.na(htmlwidget_helper)) source(htmlwidget_helper)

if (!file.exists(state_file) || !file.exists(expression_file)) {
  stop("Score audit requires clean-CDS state and gene-expression files.",
       call. = FALSE)
}

baseline_state <- "Baseline (control)"
oxphos_state <- "OXPHOS (mitochondrial)"

read_gene_expression_matrix <- function(path) {
  dt <- fread(path, showProgress = FALSE)
  if (!"gene_symbol" %chin% names(dt)) {
    stop("Gene-expression file lacks gene_symbol: ", path, call. = FALSE)
  }
  gene_symbols <- toupper(as.character(dt$gene_symbol))
  if (anyDuplicated(gene_symbols)) {
    stop("Gene-expression file contains duplicated gene symbols.", call. = FALSE)
  }
  expr <- as.matrix(dt[, setdiff(names(dt), "gene_symbol"), with = FALSE])
  storage.mode(expr) <- "numeric"
  rownames(expr) <- gene_symbols
  expr
}

mean_present_genes <- function(gene_z, genes, absent_value = NA_real_) {
  genes <- unique(toupper(genes))
  present <- intersect(genes, rownames(gene_z))
  score <- if (length(present) == 0) {
    rep(absent_value, ncol(gene_z))
  } else {
    colMeans(gene_z[present, , drop = FALSE], na.rm = TRUE)
  }
  list(
    score = as.numeric(score),
    present = length(present),
    total = length(genes)
  )
}

score_module_arms <- function(module, gene_z, run_ids, sample_names) {
  up <- mean_present_genes(gene_z, module$up_genes[[1]])
  has_down_arm <- length(module$down_genes[[1]]) > 0
  down <- mean_present_genes(
    gene_z,
    module$down_genes[[1]],
    absent_value = if (has_down_arm) NA_real_ else 0
  )
  score <- up$score - down$score
  data.table(
    Run = run_ids,
    sample = sample_names,
    dominant_state = module$dominant_state,
    signed_score = score,
    up_arm_score = up$score,
    down_arm_score = if (has_down_arm) down$score else NA_real_,
    down_low_score = if (has_down_arm) -down$score else NA_real_,
    has_down_arm = has_down_arm,
    up_genes_present = up$present,
    up_genes_total = up$total,
    down_genes_present = down$present,
    down_genes_total = down$total
  )
}

mad_or_sd <- function(x) {
  spread <- mad(x, center = median(x, na.rm = TRUE), na.rm = TRUE)
  if (!is.finite(spread) || spread <= 0) {
    spread <- sd(x, na.rm = TRUE)
  }
  if (!is.finite(spread) || spread <= 0) 1 else spread
}

percentile_z <- function(x) {
  out <- rep(NA_real_, length(x))
  valid <- is.finite(x)
  if (!any(valid)) return(out)
  p <- (rank(x[valid], ties.method = "average") - 0.5) / sum(valid)
  out[valid] <- qnorm(pmin(pmax(p, 1e-6), 1 - 1e-6))
  out
}

save_plot <- function(plot, name, width, height) {
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
  if (exists("dominant_save_ggplotly", mode = "function") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      plot,
      file.path(figure_dir, paste0(name, ".html")),
      title = name
    )
  }
}

padded_quantile_range <- function(x, probs = c(0.005, 0.995), include = 0) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(-1, 1))
  limits <- as.numeric(quantile(x, probs = probs, na.rm = TRUE))
  limits <- range(c(limits, include), na.rm = TRUE)
  span <- diff(limits)
  if (!is.finite(span) || span <= 0) span <- 1
  limits + c(-0.06, 0.06) * span
}

gene_expr <- read_gene_expression_matrix(expression_file)
state_dt <- fread(state_file, showProgress = FALSE)
if (!all(c("Run", "sample", "dominant_state", "dominant_state_score") %chin%
         names(state_dt))) {
  stop("State table lacks raw dominant-state annotation columns.", call. = FALSE)
}
if (!identical(colnames(gene_expr), state_dt$Run)) {
  missing_runs <- setdiff(state_dt$Run, colnames(gene_expr))
  if (length(missing_runs) > 0) {
    stop("Gene-expression matrix is missing state-table runs: ",
         paste(head(missing_runs, 8), collapse = ", "), call. = FALSE)
  }
  gene_expr <- gene_expr[, state_dt$Run, drop = FALSE]
}

expr_norm <- log2(gene_expr + 1)
baseline_genes <- modules[dominant_state == baseline_state, up_genes[[1]]]
baseline_present <- intersect(baseline_genes, rownames(expr_norm))
if (length(baseline_present) < 5) {
  stop("Too few baseline genes found for score audit: ",
       length(baseline_present), call. = FALSE)
}
baseline <- colMeans(expr_norm[baseline_present, , drop = FALSE], na.rm = TRUE)
gene_z <- sweep(expr_norm, 2, baseline, "-")

score_long <- rbindlist(lapply(seq_len(nrow(modules)), function(i) {
  score_module_arms(modules[i], gene_z, state_dt$Run, state_dt$sample)
}), fill = TRUE)

saved_score_states <- intersect(modules$dominant_state, names(state_dt))
saved_scores <- melt(
  state_dt[, c("Run", saved_score_states), with = FALSE],
  id.vars = "Run",
  variable.name = "dominant_state",
  value.name = "saved_signed_score",
  variable.factor = FALSE
)
score_long <- merge(
  score_long,
  saved_scores,
  by = c("Run", "dominant_state"),
  all.x = TRUE,
  sort = FALSE
)
score_long[, score_reconstruction_delta := signed_score - saved_signed_score]
max_reconstruction_delta <- score_long[
  is.finite(score_reconstruction_delta),
  max(abs(score_reconstruction_delta))
]
if (!is.finite(max_reconstruction_delta) || max_reconstruction_delta > 1e-8) {
  stop("Score audit does not reconstruct saved clean-CDS scores; max delta = ",
       signif(max_reconstruction_delta, 4), call. = FALSE)
}

audit_scores <- score_long[dominant_state != baseline_state]
audit_scores[
  ,
  `:=`(
    state_median = median(signed_score, na.rm = TRUE),
    state_mad = mad_or_sd(signed_score)
  ),
  by = dominant_state
]
audit_scores[, robust_z := (signed_score - state_median) / state_mad]
audit_scores[, percentile_z := percentile_z(signed_score), by = dominant_state]
audit_scores[, percentile := pnorm(percentile_z)]

raw_winners <- audit_scores[
  ,
  .SD[which.max(signed_score)],
  by = Run
][
  ,
  .(
    Run,
    sample,
    raw_winner = dominant_state,
    raw_winner_score = signed_score,
    raw_winner_percentile_z = percentile_z
  )
]
calibrated_winners <- audit_scores[
  ,
  .SD[which.max(percentile_z)],
  by = Run
][
  ,
  .(
    Run,
    calibrated_winner = dominant_state,
    calibrated_winner_raw_score = signed_score,
    calibrated_winner_percentile_z = percentile_z
  )
]
winner_comparison <- merge(
  raw_winners,
  calibrated_winners,
  by = "Run",
  all = TRUE,
  sort = FALSE
)
winner_comparison <- merge(
  state_dt[, .(Run, run_order = .I)],
  winner_comparison,
  by = "Run",
  all.x = TRUE,
  sort = FALSE
)
setorder(winner_comparison, run_order)
winner_comparison[, `:=`(
  run_order = NULL,
  winner_changed = raw_winner != calibrated_winner
)]

percentile_wide <- dcast(
  audit_scores[, .(Run, dominant_state, percentile_z)],
  Run ~ dominant_state,
  value.var = "percentile_z"
)
setnames(
  percentile_wide,
  setdiff(names(percentile_wide), "Run"),
  paste0("percentile_z::", setdiff(names(percentile_wide), "Run"))
)
calibrated_annotations <- merge(
  winner_comparison,
  percentile_wide,
  by = "Run",
  all.x = TRUE,
  sort = FALSE
)
calibrated_annotations <- merge(
  state_dt[, .(Run, run_order = .I)],
  calibrated_annotations,
  by = "Run",
  all.x = TRUE,
  sort = FALSE
)
setorder(calibrated_annotations, run_order)
calibrated_annotations[, run_order := NULL]

winner_counts <- rbindlist(list(
  raw_winners[, .(winner_type = "Raw signed score", N = .N),
              by = .(dominant_state = raw_winner)],
  calibrated_winners[
    ,
    .(winner_type = "Within-state percentile", N = .N),
    by = .(dominant_state = calibrated_winner)
  ]
))
winner_counts[
  ,
  fraction := N / sum(N),
  by = winner_type
]
winner_transitions <- winner_comparison[
  ,
  .(N = .N),
  by = .(raw_winner, calibrated_winner)
][
  ,
  raw_winner_fraction := N / sum(N),
  by = raw_winner
]

state_summary <- audit_scores[
  ,
  .(
    has_down_arm = first(has_down_arm),
    up_genes_present = first(up_genes_present),
    up_genes_total = first(up_genes_total),
    down_genes_present = first(down_genes_present),
    down_genes_total = first(down_genes_total),
    score_median = median(signed_score, na.rm = TRUE),
    score_mad = first(state_mad),
    score_q95 = quantile(signed_score, 0.95, na.rm = TRUE),
    score_max = max(signed_score, na.rm = TRUE),
    up_arm_median = median(up_arm_score, na.rm = TRUE),
    down_arm_median = if (first(has_down_arm)) {
      median(down_arm_score, na.rm = TRUE)
    } else {
      NA_real_
    }
  ),
  by = dominant_state
]
state_summary <- merge(
  state_summary,
  dcast(winner_counts, dominant_state ~ winner_type, value.var = "N", fill = 0),
  by = "dominant_state",
  all.x = TRUE,
  sort = FALSE
)

oxphos_scores <- audit_scores[dominant_state == oxphos_state]
if (nrow(oxphos_scores) == 0) {
  stop("Missing OXPHOS score rows in score audit.", call. = FALSE)
}
oxphos_scores <- merge(
  oxphos_scores,
  winner_comparison[
    ,
    .(Run, raw_winner, calibrated_winner, winner_changed)
  ],
  by = "Run",
  all.x = TRUE,
  sort = FALSE
)
oxphos_scores[, `:=`(
  raw_oxphos_winner = raw_winner == oxphos_state,
  calibrated_oxphos_winner = calibrated_winner == oxphos_state
)]
oxphos_scores[, oxphos_winner_class := fcase(
  raw_oxphos_winner & calibrated_oxphos_winner,
  "Raw + percentile OXPHOS",
  raw_oxphos_winner,
  "Raw OXPHOS only",
  calibrated_oxphos_winner,
  "Percentile OXPHOS only",
  default = "Not OXPHOS winner"
)]
oxphos_summary <- oxphos_scores[
  ,
  .(
    n_samples = .N,
    capacity_median = median(up_arm_score, na.rm = TRUE),
    capacity_q95 = quantile(up_arm_score, 0.95, na.rm = TRUE),
    fraction_capacity_above_baseline = mean(up_arm_score > 0, na.rm = TRUE),
    glycolysis_hypoxia_arm_median = median(down_arm_score, na.rm = TRUE),
    glycolysis_hypoxia_low_median = median(down_low_score, na.rm = TRUE),
    balance_median = median(signed_score, na.rm = TRUE),
    balance_q95 = quantile(signed_score, 0.95, na.rm = TRUE)
  ),
  by = .(raw_oxphos_winner, calibrated_oxphos_winner, oxphos_winner_class)
]

fwrite(
  audit_scores,
  file.path(output_dir, "dominant_state_score_audit_sample_scores_long.csv")
)
fwrite(
  state_summary,
  file.path(output_dir, "dominant_state_score_audit_state_summary.csv")
)
fwrite(
  calibrated_annotations,
  file.path(output_dir, "dominant_state_score_audit_calibrated_annotations.csv")
)
fwrite(
  winner_counts,
  file.path(output_dir, "dominant_state_score_audit_winner_counts.csv")
)
fwrite(
  winner_transitions,
  file.path(output_dir, "dominant_state_score_audit_winner_transitions.csv")
)
fwrite(
  winner_comparison,
  file.path(output_dir, "dominant_state_score_audit_winner_comparison.csv")
)
fwrite(
  oxphos_scores,
  file.path(output_dir, "dominant_state_score_audit_oxphos_components.csv")
)
fwrite(
  oxphos_summary,
  file.path(output_dir, "dominant_state_score_audit_oxphos_summary.csv")
)

state_levels <- state_summary[order(-`Raw signed score`, dominant_state),
                              dominant_state]
winner_counts[, dominant_state := factor(dominant_state, levels = rev(state_levels))]
p_counts <- ggplot(winner_counts, aes(N, dominant_state, fill = winner_type)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.76)) +
  geom_text(
    aes(label = N),
    position = position_dodge(width = 0.76),
    hjust = -0.12,
    size = 3.4
  ) +
  scale_fill_manual(
    values = c(
      "Raw signed score" = "#345B8C",
      "Within-state percentile" = "#B24A3C"
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.11))) +
  labs(
    title = "Dominant-state winners depend on score calibration",
    subtitle = paste(
      "Raw winners compare signed module scores directly.",
      "Percentile winners compare each score against its own compendium background."
    ),
    x = "Samples",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 11.5) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    plot.title.position = "plot"
  )
save_plot(p_counts, "dominant_state_score_audit_winner_counts", 10.4, 6.2)

transition_plot <- copy(winner_transitions)
transition_plot[, raw_winner := factor(raw_winner, levels = rev(state_levels))]
transition_plot[, calibrated_winner := factor(calibrated_winner, levels = state_levels)]
p_transition <- ggplot(
  transition_plot,
  aes(calibrated_winner, raw_winner, fill = raw_winner_fraction)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = N), size = 3.15) +
  scale_fill_gradient(
    low = "#EFF4FA",
    high = "#345B8C",
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  labs(
    title = "Raw winner to percentile-calibrated winner transitions",
    subtitle = "Fill is the fraction within each raw winner row; labels are sample counts.",
    x = "Within-state percentile winner",
    y = "Raw signed-score winner",
    fill = "Row share"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title.position = "plot"
  )
save_plot(p_transition, "dominant_state_score_audit_winner_transitions", 10.8, 8.2)

oxphos_scores[, oxphos_winner_class := factor(
  oxphos_winner_class,
  levels = c(
    "Not OXPHOS winner",
    "Raw OXPHOS only",
    "Percentile OXPHOS only",
    "Raw + percentile OXPHOS"
  )
)]
p_oxphos <- ggplot(
  oxphos_scores,
  aes(up_arm_score, down_low_score, color = oxphos_winner_class)
) +
  geom_hline(yintercept = 0, color = "#8C8C8C", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "#8C8C8C", linewidth = 0.35) +
  geom_point(alpha = 0.55, size = 1.25) +
  scale_color_manual(
    values = c(
      "Not OXPHOS winner" = "#8A939D",
      "Raw OXPHOS only" = "#345B8C",
      "Percentile OXPHOS only" = "#D89A36",
      "Raw + percentile OXPHOS" = "#B24A3C"
    )
  ) +
  coord_cartesian(
    xlim = padded_quantile_range(oxphos_scores$up_arm_score),
    ylim = padded_quantile_range(oxphos_scores$down_low_score)
  ) +
  labs(
    title = "OXPHOS balance score separates capacity from a low down arm",
    subtitle = paste(
      "x = OXPHOS-capacity arm; y = negative glycolysis/hypoxia arm.",
      "Signed balance score is x + y."
    ),
    x = "OXPHOS capacity arm relative to baseline",
    y = "Glycolysis/hypoxia low arm relative to baseline",
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot"
  )
save_plot(p_oxphos, "dominant_state_score_audit_oxphos_components", 9.6, 7)

stopifnot(nrow(calibrated_annotations) == nrow(state_dt))
stopifnot(identical(calibrated_annotations$Run, state_dt$Run))
stopifnot(all(
  oxphos_scores[
    is.finite(signed_score) & is.finite(up_arm_score) &
      is.finite(down_arm_score),
    abs(signed_score - (up_arm_score - down_arm_score))
  ] < 1e-8
))

message("Dominant-state score audit:")
message("  score reconstruction max delta: ", signif(max_reconstruction_delta, 4))
message("  winner changes after within-state percentile calibration: ",
        sum(winner_comparison$winner_changed, na.rm = TRUE), " / ",
        nrow(winner_comparison))
message("  saved outputs: ", output_dir)

print(state_summary[order(-`Raw signed score`)])
print(oxphos_summary[order(-n_samples)])
