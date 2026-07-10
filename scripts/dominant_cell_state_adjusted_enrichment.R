dominant_loader <- c(
  file.path("scripts", "dominant_ribocrypt_loader.R"),
  file.path("dominant_cell_states", "scripts", "dominant_ribocrypt_loader.R")
)
dominant_loader <- dominant_loader[file.exists(dominant_loader)][1]
if (!is.na(dominant_loader)) {
  source(dominant_loader)
  dominant_load_ribocrypt_if_available()
}

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

state_file <- result_file("human_dominant_cell_states_clean_cds.csv")
metadata_file <- "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv"
metadata_fields <- c("CELL_LINE", "TISSUE", "GENE", "CONDITION")
min_level_n <- 3L
min_covariate_n <- 5L
max_covariate_levels <- 80L
output_prefix <- "dominant_state_clean_cds_enrichment"

clean_level <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "NA", "N/A", "na", "n/a", "NULL", "null", "None", "none")] <- NA_character_
  x
}

collapse_covariate <- function(x, min_n = min_covariate_n,
                               max_levels = max_covariate_levels) {
  x <- clean_level(x)
  x[is.na(x)] <- "MISSING"
  tab <- sort(table(x), decreasing = TRUE)
  keep <- names(tab)[tab >= min_n]
  if (length(keep) > max_levels) keep <- keep[seq_len(max_levels)]
  x[!x %chin% keep] <- "OTHER"
  factor(x)
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
  test <- tryCatch(fisher.test(mat, alternative = "greater"),
                   error = function(e) NULL)
  if (is.null(test)) return(list(p = NA_real_, odds_ratio = NA_real_))
  list(p = test$p.value, odds_ratio = unname(test$estimate))
}

adjusted_residuals <- function(dt, score_col, adjust_fields) {
  y <- as.numeric(dt[[score_col]])
  out <- rep(NA_real_, length(y))
  valid <- is.finite(y)
  if (!any(valid)) return(out)

  if (length(adjust_fields) == 0) {
    out[valid] <- y[valid] - mean(y[valid], na.rm = TRUE)
    return(out)
  }

  model_dt <- data.table(score = y[valid])
  for (field in adjust_fields) {
    model_dt[, (field) := collapse_covariate(dt[[field]][valid])]
  }
  rhs <- paste(adjust_fields, collapse = " + ")
  mm <- tryCatch(
    model.matrix(as.formula(paste("score ~", rhs)), model_dt),
    error = function(e) NULL
  )
  if (is.null(mm) || ncol(mm) == 0) {
    out[valid] <- y[valid] - mean(y[valid], na.rm = TRUE)
    return(out)
  }

  fit <- tryCatch(lm.fit(mm, model_dt$score), error = function(e) NULL)
  if (is.null(fit)) {
    out[valid] <- y[valid] - mean(y[valid], na.rm = TRUE)
  } else {
    out[valid] <- fit$residuals
  }
  out
}

enrich_one_field_adjusted <- function(dt, field, state_cols) {
  adjust_fields <- setdiff(metadata_fields, field)
  field_levels <- clean_level(dt[[field]])
  valid_field <- !is.na(field_levels)
  levels <- sort(unique(field_levels[valid_field]))

  rbindlist(lapply(state_cols, function(state) {
    residual <- adjusted_residuals(dt, state, adjust_fields)
    raw_score <- as.numeric(dt[[state]])

    state_res <- rbindlist(lapply(levels, function(current_level) {
      in_level <- valid_field & field_levels == current_level
      out_level <- valid_field & field_levels != current_level
      n_in <- sum(in_level)
      n_out <- sum(out_level)

      residual_sd <- sd(residual[valid_field], na.rm = TRUE)
      residual_delta <- mean(residual[in_level], na.rm = TRUE) -
        mean(residual[out_level], na.rm = TRUE)
      raw_sd <- sd(raw_score[valid_field], na.rm = TRUE)
      raw_delta <- mean(raw_score[in_level], na.rm = TRUE) -
        mean(raw_score[out_level], na.rm = TRUE)

      call_in <- sum(dt$dominant_state[in_level] == state, na.rm = TRUE)
      call_out <- sum(dt$dominant_state[out_level] == state, na.rm = TRUE)
      no_call_in <- n_in - call_in
      no_call_out <- n_out - call_out
      fisher <- safe_fisher_greater(call_in, no_call_in, call_out, no_call_out)

      data.table(
        field = field,
        level = current_level,
        dominant_state = state,
        n_in = n_in,
        n_out = n_out,
        adjusted_effect_z = if (is.finite(residual_sd) && residual_sd > 0) {
          residual_delta / residual_sd
        } else NA_real_,
        adjusted_score_delta = residual_delta,
        adjusted_wilcox_greater_p = safe_wilcox_greater(
          residual[in_level], residual[out_level]
        ),
        raw_effect_z = if (is.finite(raw_sd) && raw_sd > 0) raw_delta / raw_sd else NA_real_,
        raw_score_delta = raw_delta,
        raw_wilcox_greater_p = safe_wilcox_greater(
          raw_score[in_level], raw_score[out_level]
        ),
        dominant_calls_in = call_in,
        dominant_calls_out = call_out,
        dominant_call_fraction_in = call_in / n_in,
        dominant_call_fraction_out = call_out / n_out,
        dominant_call_fraction_delta = (call_in / n_in) - (call_out / n_out),
        fisher_odds_ratio = fisher$odds_ratio,
        fisher_greater_p = fisher$p
      )
    }))
    state_res
  }), fill = TRUE)[
    ,
    `:=`(
      adjusted_wilcox_greater_p_adj = p.adjust(adjusted_wilcox_greater_p, "BH"),
      raw_wilcox_greater_p_adj = p.adjust(raw_wilcox_greater_p, "BH"),
      fisher_greater_p_adj = p.adjust(fisher_greater_p, "BH")
    )
  ][
    ,
    p_missing := is.na(adjusted_wilcox_greater_p_adj)
  ][
    order(p_missing, adjusted_wilcox_greater_p_adj, -adjusted_effect_z,
          level, dominant_state)
  ][
    ,
    p_missing := NULL
  ][]
}

enrich_interactions_adjusted <- function(dt, state_cols) {
  pairs <- combn(metadata_fields, 2, simplify = FALSE)

  rbindlist(lapply(state_cols, function(state) {
    residual <- adjusted_residuals(dt, state, metadata_fields)

    rbindlist(lapply(pairs, function(pair) {
      f1 <- pair[[1]]
      f2 <- pair[[2]]
      v1 <- clean_level(dt[[f1]])
      v2 <- clean_level(dt[[f2]])
      valid <- !is.na(v1) & !is.na(v2)
      interaction_level <- paste0(f1, "=", v1, " | ", f2, "=", v2)
      levels <- sort(unique(interaction_level[valid]))

      rbindlist(lapply(levels, function(level) {
        in_level <- valid & interaction_level == level
        out_level <- valid & interaction_level != level
        n_in <- sum(in_level)
        n_out <- sum(out_level)
        if (n_in < min_level_n || n_out < min_level_n) return(NULL)

        residual_sd <- sd(residual[valid], na.rm = TRUE)
        residual_delta <- mean(residual[in_level], na.rm = TRUE) -
          mean(residual[out_level], na.rm = TRUE)

        data.table(
          field = paste(f1, f2, sep = ":"),
          level = level,
          dominant_state = state,
          n_in = n_in,
          n_out = n_out,
          adjusted_effect_z = if (is.finite(residual_sd) && residual_sd > 0) {
            residual_delta / residual_sd
          } else NA_real_,
          adjusted_score_delta = residual_delta,
          adjusted_wilcox_greater_p = safe_wilcox_greater(
            residual[in_level], residual[out_level]
          )
        )
      }))
    }), fill = TRUE)
  }), fill = TRUE)[
    ,
    adjusted_wilcox_greater_p_adj := p.adjust(adjusted_wilcox_greater_p, "BH")
  ][
    ,
    p_missing := is.na(adjusted_wilcox_greater_p_adj)
  ][
    order(p_missing, adjusted_wilcox_greater_p_adj, -adjusted_effect_z,
          field, level, dominant_state)
  ][
    ,
    p_missing := NULL
  ][]
}

save_adjusted_heatmap <- function(enrich_dt, title_field, state_cols, output_file,
                                  max_rows = Inf) {
  if (nrow(enrich_dt) == 0) return(invisible(NULL))
  plot_dt <- copy(enrich_dt)
  plot_dt[, row_label := ifelse(field == title_field, level,
                                paste(field, level, sep = " :: "))]
  plot_dt[, hover := paste0(
    "Term: ", row_label,
    "<br>State: ", dominant_state,
    "<br>n in/out: ", n_in, " / ", n_out,
    "<br>adjusted effect z: ", signif(adjusted_effect_z, 3),
    "<br>adjusted delta: ", signif(adjusted_score_delta, 3),
    "<br>adjusted Wilcoxon adj p: ",
    signif(adjusted_wilcox_greater_p_adj, 3)
  )]

  level_order <- plot_dt[
    ,
    .(
      best_effect_z = max(adjusted_effect_z, na.rm = TRUE),
      best_p = suppressWarnings(min(adjusted_wilcox_greater_p_adj, na.rm = TRUE)),
      n = max(n_in, na.rm = TRUE)
    ),
    by = row_label
  ]
  level_order[!is.finite(best_effect_z), best_effect_z := NA_real_]
  level_order[!is.finite(best_p), best_p := NA_real_]
  level_order[, p_missing := is.na(best_p)]
  setorder(level_order, p_missing, best_p, -best_effect_z, -n, row_label)
  level_order[, p_missing := NULL]
  if (is.finite(max_rows) && nrow(level_order) > max_rows) {
    level_order <- level_order[seq_len(max_rows)]
  }

  z_wide <- dcast(plot_dt, row_label ~ dominant_state,
                  value.var = "adjusted_effect_z")
  text_wide <- dcast(plot_dt, row_label ~ dominant_state, value.var = "hover")
  z_wide <- z_wide[level_order, on = "row_label"]
  text_wide <- text_wide[level_order, on = "row_label"]

  missing_states <- setdiff(state_cols, names(z_wide))
  if (length(missing_states) > 0) {
    z_wide[, (missing_states) := NA_real_]
    text_wide[, (missing_states) := ""]
  }

  z_matrix <- as.matrix(z_wide[, ..state_cols])
  text_matrix <- as.matrix(text_wide[, ..state_cols])
  max_abs <- suppressWarnings(quantile(abs(z_matrix), 0.98, na.rm = TRUE))
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1

  p <- plot_ly(
    x = state_cols,
    y = z_wide$row_label,
    z = z_matrix,
    text = text_matrix,
    hoverinfo = "text",
    type = "heatmap",
    zmin = -max_abs,
    zmax = max_abs,
    colorscale = list(c(0, "#2166AC"), c(0.5, "#F7F7F7"), c(1, "#B2182B")),
    colorbar = list(title = "adjusted effect z")
  ) %>%
    layout(
      title = paste("Adjusted dominant-state enrichment:", title_field),
      xaxis = list(title = "", tickangle = -35, automargin = TRUE),
      yaxis = list(title = "", automargin = TRUE),
      margin = list(l = 260, b = 150, t = 70, r = 30)
    )

  dominant_save_widget(
    p,
    output_file = output_file,
    title = paste("Adjusted enrichment", title_field)
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

state_cols <- setdiff(
  names(dt),
  c("Run", "sample", "dominant_state", "dominant_state_score",
    "Baseline (control)", metadata_fields)
)

all_main <- rbindlist(lapply(metadata_fields, function(field) {
  message("Adjusted enrichment for ", field)
  enrichment <- enrich_one_field_adjusted(dt, field, state_cols)
  csv_file <- result_file(paste0(output_prefix, "_", tolower(field), ".csv"))
  heatmap_file <- result_file(paste0(output_prefix, "_", tolower(field), "_heatmap.html"))
  fwrite(enrichment, csv_file)
  save_adjusted_heatmap(enrichment, field, state_cols, heatmap_file)
  enrichment
}), fill = TRUE)

fwrite(all_main, result_file(paste0(output_prefix, "_all_main.csv")))

message("Adjusted interaction enrichment")
interactions <- enrich_interactions_adjusted(dt, state_cols)
fwrite(interactions, result_file(paste0(output_prefix, "_interactions.csv")))
save_adjusted_heatmap(
  interactions,
  "interactions",
  state_cols,
  result_file(paste0(output_prefix, "_interactions_heatmap.html")),
  max_rows = 200L
)

message("Saved adjusted enrichment CSVs and heatmaps.")

invisible(list(main = all_main, interactions = interactions))
