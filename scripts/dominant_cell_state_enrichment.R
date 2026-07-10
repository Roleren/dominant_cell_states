library(data.table)
library(plotly)
library(htmlwidgets)

htmlwidget_helper <- c(
  file.path("scripts", "dominant_htmlwidgets.R"),
  "dominant_htmlwidgets.R",
  file.path("dominant_cell_states", "scripts", "dominant_htmlwidgets.R")
)
htmlwidget_helper <- htmlwidget_helper[file.exists(htmlwidget_helper)][1]
if (!is.na(htmlwidget_helper)) source(htmlwidget_helper)

analysis_dir <- if (dir.exists("dominant_cell_states")) "dominant_cell_states" else "."
results_dir <- Sys.getenv(
  "DOMINANT_CELL_STATES_RESULTS_DIR",
  unset = file.path(analysis_dir, "results")
)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- function(...) file.path(results_dir, ...)

state_file <- result_file("human_dominant_cell_states.csv")
metadata_file <- "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv"
metadata_fields <- c("CELL_LINE", "TISSUE", "GENE", "CONDITION")
min_level_n <- 3L

output_prefix <- "dominant_state_enrichment"

clean_level <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "NA", "N/A", "na", "n/a", "NULL", "null", "None", "none")] <- NA_character_
  x
}

safe_wilcox_greater <- function(in_score, out_score) {
  in_score <- in_score[is.finite(in_score)]
  out_score <- out_score[is.finite(out_score)]
  if (length(in_score) < min_level_n || length(out_score) < min_level_n) {
    return(NA_real_)
  }
  if (length(unique(c(in_score, out_score))) < 2) return(NA_real_)

  tryCatch(
    wilcox.test(in_score, out_score, alternative = "greater",
                exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

safe_fisher_greater <- function(call_in, no_call_in, call_out, no_call_out) {
  mat <- matrix(c(call_in, no_call_in, call_out, no_call_out),
                nrow = 2, byrow = TRUE)
  if (sum(mat[1, ]) < min_level_n || sum(mat[2, ]) < min_level_n) {
    return(list(p = NA_real_, odds_ratio = NA_real_))
  }

  test <- tryCatch(
    fisher.test(mat, alternative = "greater"),
    error = function(e) NULL
  )
  if (is.null(test)) return(list(p = NA_real_, odds_ratio = NA_real_))

  list(
    p = test$p.value,
    odds_ratio = unname(test$estimate)
  )
}

enrich_one_field <- function(dt, field, state_cols) {
  field_dt <- copy(dt)
  field_dt[, level := clean_level(get(field))]
  field_dt <- field_dt[!is.na(level)]

  if (nrow(field_dt) == 0) {
    warning("No non-missing metadata levels for ", field)
    return(data.table())
  }

  levels <- sort(unique(field_dt$level))

  res <- rbindlist(lapply(levels, function(current_level) {
    in_level <- field_dt$level == current_level
    n_in <- sum(in_level)
    n_out <- sum(!in_level)

    rbindlist(lapply(state_cols, function(state) {
      all_scores <- field_dt[[state]]
      in_score <- all_scores[in_level]
      out_score <- all_scores[!in_level]
      score_sd <- sd(all_scores, na.rm = TRUE)
      score_delta <- mean(in_score, na.rm = TRUE) - mean(out_score, na.rm = TRUE)

      call_in <- sum(field_dt$dominant_state[in_level] == state, na.rm = TRUE)
      call_out <- sum(field_dt$dominant_state[!in_level] == state, na.rm = TRUE)
      no_call_in <- n_in - call_in
      no_call_out <- n_out - call_out
      fisher <- safe_fisher_greater(call_in, no_call_in, call_out, no_call_out)

      data.table(
        field = field,
        level = current_level,
        dominant_state = state,
        n_in = n_in,
        n_out = n_out,
        mean_score_in = mean(in_score, na.rm = TRUE),
        mean_score_out = mean(out_score, na.rm = TRUE),
        median_score_in = median(in_score, na.rm = TRUE),
        median_score_out = median(out_score, na.rm = TRUE),
        score_delta = score_delta,
        effect_z = if (is.finite(score_sd) && score_sd > 0) score_delta / score_sd else NA_real_,
        wilcox_greater_p = safe_wilcox_greater(in_score, out_score),
        dominant_calls_in = call_in,
        dominant_calls_out = call_out,
        dominant_call_fraction_in = call_in / n_in,
        dominant_call_fraction_out = call_out / n_out,
        dominant_call_fraction_delta = (call_in / n_in) - (call_out / n_out),
        fisher_greater_p = fisher$p,
        fisher_odds_ratio = fisher$odds_ratio
      )
    }))
  }))

  res[, wilcox_greater_p_adj := p.adjust(wilcox_greater_p, method = "BH")]
  res[, fisher_greater_p_adj := p.adjust(fisher_greater_p, method = "BH")]
  setcolorder(
    res,
    c("field", "level", "dominant_state", "n_in", "n_out",
      "effect_z", "score_delta", "wilcox_greater_p", "wilcox_greater_p_adj",
      "mean_score_in", "mean_score_out", "median_score_in", "median_score_out",
      "dominant_calls_in", "dominant_calls_out",
      "dominant_call_fraction_in", "dominant_call_fraction_out",
      "dominant_call_fraction_delta", "fisher_odds_ratio",
      "fisher_greater_p", "fisher_greater_p_adj")
  )
  res[, p_missing := is.na(wilcox_greater_p_adj)]
  setorder(res, p_missing, wilcox_greater_p_adj, -effect_z, level, dominant_state)
  res[, p_missing := NULL]
  res
}

save_enrichment_heatmap <- function(enrich_dt, field, state_cols, output_file) {
  if (nrow(enrich_dt) == 0) return(invisible(NULL))

  plot_dt <- copy(enrich_dt)
  plot_dt[, hover := paste0(
    field, ": ", level,
    "<br>State: ", dominant_state,
    "<br>n in/out: ", n_in, " / ", n_out,
    "<br>effect z: ", signif(effect_z, 3),
    "<br>score delta: ", signif(score_delta, 3),
    "<br>Wilcoxon adj p: ", signif(wilcox_greater_p_adj, 3),
    "<br>call fraction in/out: ",
    signif(dominant_call_fraction_in, 3), " / ",
    signif(dominant_call_fraction_out, 3),
    "<br>Fisher adj p: ", signif(fisher_greater_p_adj, 3)
  )]

  level_order <- plot_dt[
    ,
    .(
      best_effect_z = max(effect_z, na.rm = TRUE),
      best_p = suppressWarnings(min(wilcox_greater_p_adj, na.rm = TRUE)),
      n = max(n_in, na.rm = TRUE)
    ),
    by = level
  ]
  level_order[!is.finite(best_effect_z), best_effect_z := NA_real_]
  level_order[!is.finite(best_p), best_p := NA_real_]
  level_order[, p_missing := is.na(best_p)]
  setorder(level_order, p_missing, best_p, -best_effect_z, -n, level)
  level_order[, p_missing := NULL]

  z_wide <- dcast(plot_dt, level ~ dominant_state, value.var = "effect_z")
  text_wide <- dcast(plot_dt, level ~ dominant_state, value.var = "hover")
  z_wide <- z_wide[level_order, on = "level"]
  text_wide <- text_wide[level_order, on = "level"]

  missing_states <- setdiff(state_cols, names(z_wide))
  if (length(missing_states) > 0) {
    z_wide[, (missing_states) := NA_real_]
    text_wide[, (missing_states) := ""]
  }

  z_matrix <- as.matrix(z_wide[, ..state_cols])
  text_matrix <- as.matrix(text_wide[, ..state_cols])
  row_labels <- z_wide$level
  max_abs <- suppressWarnings(quantile(abs(z_matrix), 0.98, na.rm = TRUE))
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1

  p <- plot_ly(
    x = state_cols,
    y = row_labels,
    z = z_matrix,
    text = text_matrix,
    hoverinfo = "text",
    type = "heatmap",
    zmin = -max_abs,
    zmax = max_abs,
    colorscale = list(
      c(0, "#2166AC"),
      c(0.5, "#F7F7F7"),
      c(1, "#B2182B")
    ),
    colorbar = list(title = "effect z")
  ) %>%
    layout(
      title = paste("Dominant-state score enrichment by", field),
      xaxis = list(title = "", tickangle = -35, automargin = TRUE),
      yaxis = list(title = "", automargin = TRUE),
      margin = list(l = 220, b = 150, t = 70, r = 30)
    )

  dominant_save_widget(
    p,
    output_file = output_file,
    title = paste("Dominant-state enrichment", field)
  )
}

dt <- fread(state_file)
m <- fread(metadata_file)
m <- m[, .(Run, CELL_LINE, TISSUE, CONDITION, GENE)]

matching_duplicate_runs <- m[
  Run %chin% dt$Run & (duplicated(Run) | duplicated(Run, fromLast = TRUE))
]
if (nrow(matching_duplicate_runs) > 0) {
  stop("Metadata has duplicate Run rows matching the dominant-state table: ",
       paste(unique(matching_duplicate_runs$Run), collapse = ", "))
}

dt[, input_order := .I]
dt <- merge.data.table(dt, m, all.x = TRUE, sort = FALSE)
setorder(dt, input_order)
dt[, input_order := NULL]

score_cols <- setdiff(
  names(dt),
  c("Run", "sample", "dominant_state", "dominant_state_score",
    "Baseline (control)", metadata_fields)
)

all_enrichment <- rbindlist(lapply(metadata_fields, function(field) {
  message("Testing enrichment for ", field)
  enrichment <- enrich_one_field(dt, field, score_cols)
  csv_file <- result_file(paste0(output_prefix, "_", tolower(field), ".csv"))
  heatmap_file <- result_file(paste0(output_prefix, "_", tolower(field), "_heatmap.html"))
  fwrite(enrichment, csv_file)
  save_enrichment_heatmap(enrichment, field, score_cols, heatmap_file)
  enrichment
}), fill = TRUE)

fwrite(all_enrichment, result_file(paste0(output_prefix, "_all.csv")))

message("Saved enrichment CSVs and heatmaps for: ",
        paste(metadata_fields, collapse = ", "))

invisible(all_enrichment)
