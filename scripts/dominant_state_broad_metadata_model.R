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

analysis_dir <- if (dir.exists("dominant_cell_states")) "dominant_cell_states" else "."
results_dir <- Sys.getenv(
  "DOMINANT_CELL_STATES_RESULTS_DIR",
  unset = file.path(analysis_dir, "results")
)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- function(...) file.path(results_dir, ...)

state_file <- result_file("human_dominant_cell_states_clean_cds.csv")
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
output_prefix <- "dominant_state_broad_metadata"

selected_fields_output <- result_file(paste0(output_prefix, "_selected_fields.csv"))
numeric_fields_output <- result_file(paste0(output_prefix, "_numeric_fields.csv"))
main_terms_output <- result_file(paste0(output_prefix, "_main_terms.csv"))
nwise_terms_output <- result_file(paste0(output_prefix, "_nwise_terms.csv"))
rf_importance_output <- result_file(paste0(output_prefix, "_rf_importance.csv"))
rf_performance_output <- result_file(paste0(output_prefix, "_rf_performance.csv"))
rdg_terms_output <- result_file(paste0(output_prefix, "_rdg_terms.csv"))
compact_summary_output <- result_file(paste0(output_prefix, "_compact_summary.csv"))
nwise_heatmap_output <- result_file(paste0(output_prefix, "_nwise_heatmap.html"))
rf_heatmap_output <- result_file(paste0(output_prefix, "_rf_importance_heatmap.html"))

min_field_nonmissing_fraction <- 0.05
min_level_n <- 8L
min_nwise_n <- 8L
min_nwise_out <- 30L
max_categorical_levels <- 60L
max_rf_factor_levels <- 32L
max_interaction_order <- 4L
max_nwise_fields <- 14L
max_heatmap_rows <- 120L
rf_ntree <- 150L

priority_fields <- c(
  "Cancer_type", "Cell_model", "Cell_type", "Organ_system", "Sex",
  "Life_stage", "CONDITION", "TISSUE", "CELL_LINE", "GENE", "INHIBITOR",
  "FRACTION", "LibrarySelection", "LibrarySource", "Model", "YEAR", "MONTH",
  "study", "AUTHOR", "mode"
)
nwise_priority_fields <- c(
  "Cancer_type", "Cell_model", "Cell_type", "Organ_system", "Sex",
  "Life_stage", "CONDITION", "TISSUE", "CELL_LINE", "GENE", "INHIBITOR",
  "FRACTION", "study", "AUTHOR"
)

identifier_patterns <- paste(
  c(
    "accession", "^Run$", "^Experiment$", "^LibraryName$", "^SRAStudy$",
    "^BioProject$", "^ProjectID$", "^Sample$", "^BioSample$", "^SampleName$",
    "^sample_id$", "^sample$", "^Submission$", "^CenterName$",
    "^Study_Pubmed_id$", "^TaxID$", "^ScientificName$", "title", "source$"
  ),
  collapse = "|"
)

clean_level <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "NA", "N/A", "na", "n/a", "NULL", "null", "None", "none",
               "missing", "Missing", "unknown", "Unknown")] <- NA_character_
  x
}

file_safe <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  substr(x, 1L, 160L)
}

safe_name <- function(x) {
  out <- gsub("[^A-Za-z0-9_]+", "_", x)
  out <- gsub("^([0-9])", "X\\1", out)
  make.unique(out)
}

collapse_levels <- function(x, min_n = min_level_n,
                            max_levels = max_categorical_levels) {
  x <- clean_level(x)
  x[is.na(x)] <- "MISSING"
  tab <- sort(table(x), decreasing = TRUE)
  keep <- names(tab)[tab >= min_n]
  keep <- setdiff(keep, "MISSING")
  if (length(keep) > max_levels) keep <- keep[seq_len(max_levels)]
  x[!x %chin% c(keep, "MISSING")] <- "OTHER"
  x
}

is_artificial_level <- function(x) {
  x %chin% c("MISSING", "OTHER")
}

score_group_terms <- function(values_dt, group_cols, score_matrix, state_cols,
                              score_type, min_in, min_out,
                              complete_terms_only = FALSE) {
  if (length(group_cols) == 0) return(data.table())
  group_values <- values_dt[, ..group_cols]
  artificial <- Reduce(`|`, lapply(group_cols, function(field) {
    is_artificial_level(group_values[[field]])
  }))
  term_level <- do.call(paste, c(lapply(group_cols, function(field) {
    paste0(field, "=", group_values[[field]])
  }), sep = " | "))
  valid_term <- !is.na(term_level)
  if (complete_terms_only) valid_term <- valid_term & !artificial
  if (!any(valid_term)) return(data.table())

  term_level <- term_level[valid_term]
  artificial <- artificial[valid_term]
  field_name <- paste(group_cols, collapse = ":")
  term_order <- length(group_cols)

  rbindlist(lapply(state_cols, function(state) {
    score <- as.numeric(score_matrix[valid_term, state])
    finite <- is.finite(score)
    if (sum(finite) < min_in + min_out) return(NULL)
    tmp <- data.table(
      level = term_level[finite],
      term_has_artificial_level = artificial[finite],
      score = score[finite]
    )
    total_n <- nrow(tmp)
    total_sum <- sum(tmp$score)
    total_sumsq <- sum(tmp$score^2)
    total_sd <- sd(tmp$score)
    if (!is.finite(total_sd) || total_sd <= 0) return(NULL)

    agg <- tmp[
      ,
      .(
        n_in = .N,
        term_has_artificial_level = any(term_has_artificial_level),
        sum_in = sum(score),
        sumsq_in = sum(score^2)
      ),
      by = level
    ]
    agg[, `:=`(
      field = field_name,
      term_order = term_order,
      dominant_state = state,
      score_type = score_type,
      n_out = total_n - n_in
    )]
    agg <- agg[n_in >= min_in & n_out >= min_out]
    if (nrow(agg) == 0) return(NULL)
    agg[, `:=`(
      sum_out = total_sum - sum_in,
      sumsq_out = total_sumsq - sumsq_in,
      mean_in = sum_in / n_in,
      mean_out = (total_sum - sum_in) / n_out
    )]
    agg[, `:=`(
      var_in = pmax((sumsq_in - (sum_in^2 / n_in)) / pmax(n_in - 1, 1), 0),
      var_out = pmax((sumsq_out - (sum_out^2 / n_out)) / pmax(n_out - 1, 1), 0)
    )]
    agg[, se_delta := sqrt(var_in / n_in + var_out / n_out)]
    agg[, effect_z := (mean_in - mean_out) / total_sd]
    agg[, t_stat := fifelse(is.finite(se_delta) & se_delta > 0,
                             (mean_in - mean_out) / se_delta, NA_real_)]
    agg[, df := {
      numerator <- (var_in / n_in + var_out / n_out)^2
      denominator <- ((var_in / n_in)^2 / pmax(n_in - 1, 1)) +
        ((var_out / n_out)^2 / pmax(n_out - 1, 1))
      fifelse(is.finite(denominator) & denominator > 0,
              numerator / denominator, NA_real_)
    }]
    agg[, greater_p := pt(t_stat, df = df, lower.tail = FALSE)]
    agg[, less_p := pt(t_stat, df = df, lower.tail = TRUE)]
    agg[, .(
      field, level, term_order, dominant_state, score_type,
      n_in, n_out, mean_in, mean_out,
      score_delta = mean_in - mean_out,
      effect_z, t_stat, df, greater_p, less_p,
      term_has_artificial_level
    )]
  }), fill = TRUE)
}

select_metadata_fields <- function(dt, metadata_cols) {
  n <- nrow(dt)
  candidates <- setdiff(metadata_cols, "Run")
  candidates <- intersect(candidates, names(dt))
  candidates <- candidates[!grepl(identifier_patterns, candidates, ignore.case = TRUE)]
  candidates <- candidates[!vapply(dt[, ..candidates], is.numeric, logical(1))]
  candidates <- unique(c(intersect(priority_fields, candidates), candidates))

  rbindlist(lapply(candidates, function(field) {
    x <- clean_level(dt[[field]])
    non_missing <- sum(!is.na(x))
    unique_n <- uniqueN(x, na.rm = TRUE)
    data.table(
      field = field,
      non_missing = non_missing,
      non_missing_fraction = non_missing / n,
      raw_unique_n = unique_n,
      selected = non_missing / n >= min_field_nonmissing_fraction &&
        unique_n >= 2 &&
        unique_n <= 700
    )
  }), fill = TRUE)[
    ,
    selected := selected | field %chin% priority_fields
  ][
    raw_unique_n < 2,
    selected := FALSE
  ][
    order(!field %chin% priority_fields, -non_missing_fraction, raw_unique_n)
  ][]
}

select_numeric_fields <- function(dt, metadata_cols) {
  candidates <- setdiff(metadata_cols, "Run")
  candidates <- intersect(candidates, names(dt))
  candidates <- candidates[vapply(dt[, ..candidates], is.numeric, logical(1))]
  candidates <- candidates[!grepl(identifier_patterns, candidates, ignore.case = TRUE)]
  n <- nrow(dt)
  rbindlist(lapply(candidates, function(field) {
    x <- as.numeric(dt[[field]])
    finite <- is.finite(x)
    data.table(
      field = field,
      non_missing = sum(finite),
      non_missing_fraction = sum(finite) / n,
      unique_n = uniqueN(x[finite]),
      selected = sum(finite) / n >= 0.25 && uniqueN(x[finite]) >= 5
    )
  }), fill = TRUE)[order(-selected, -non_missing_fraction, field)][]
}

build_residual_matrix <- function(dt, state_cols, categorical_fields,
                                  numeric_fields) {
  safe_fields <- safe_name(c(categorical_fields, numeric_fields))
  names(safe_fields) <- c(categorical_fields, numeric_fields)
  model_dt <- data.table(row_id = seq_len(nrow(dt)))
  for (field in categorical_fields) {
    model_dt[, (safe_fields[[field]]) := factor(dt[[field]])]
  }
  for (field in numeric_fields) {
    x <- as.numeric(dt[[field]])
    med <- median(x[is.finite(x)], na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    x[!is.finite(x)] <- med
    sx <- sd(x)
    if (is.finite(sx) && sx > 0) x <- (x - mean(x)) / sx else x <- 0
    model_dt[, (safe_fields[[field]]) := x]
  }

  rhs <- setdiff(names(model_dt), "row_id")
  if (length(rhs) == 0) {
    out <- as.matrix(dt[, ..state_cols])
    storage.mode(out) <- "numeric"
    return(out)
  }

  mm <- model.matrix(~ . , data = model_dt[, ..rhs])
  out <- matrix(NA_real_, nrow = nrow(dt), ncol = length(state_cols),
                dimnames = list(NULL, state_cols))
  for (state in state_cols) {
    y <- as.numeric(dt[[state]])
    valid <- is.finite(y)
    if (sum(valid) <= ncol(mm)) {
      out[valid, state] <- y[valid] - mean(y[valid], na.rm = TRUE)
      next
    }
    fit <- tryCatch(lm.fit(mm[valid, , drop = FALSE], y[valid]),
                    error = function(e) NULL)
    if (is.null(fit)) {
      out[valid, state] <- y[valid] - mean(y[valid], na.rm = TRUE)
    } else {
      out[valid, state] <- fit$residuals
    }
  }
  out
}

save_state_heatmap <- function(term_dt, output_file, value_col,
                               title, max_rows = max_heatmap_rows) {
  plot_dt <- copy(term_dt)
  plot_dt <- plot_dt[
    is.finite(get(value_col)) &
      is.finite(greater_p_adj) &
      effect_z > 0
  ]
  if (nrow(plot_dt) == 0) return(invisible(NULL))
  plot_dt[, row_label := paste(field, level, sep = " :: ")]
  row_order <- plot_dt[
    ,
    .(
      best_q = min(greater_p_adj, na.rm = TRUE),
      best_effect = max(effect_z, na.rm = TRUE),
      n = max(n_in, na.rm = TRUE)
    ),
    by = row_label
  ][order(best_q, -best_effect, -n)]
  row_order <- row_order[seq_len(min(.N, max_rows))]
  plot_dt <- plot_dt[row_label %chin% row_order$row_label]
  plot_dt[, hover := paste0(
    "Term: ", row_label,
    "<br>State: ", dominant_state,
    "<br>order: ", term_order,
    "<br>n in/out: ", n_in, " / ", n_out,
    "<br>effect z: ", signif(effect_z, 3),
    "<br>delta: ", signif(score_delta, 3),
    "<br>q: ", signif(greater_p_adj, 3)
  )]
  z_wide <- dcast(plot_dt, row_label ~ dominant_state, value.var = value_col)
  text_wide <- dcast(plot_dt, row_label ~ dominant_state, value.var = "hover")
  z_wide <- z_wide[row_order, on = "row_label"]
  text_wide <- text_wide[row_order, on = "row_label"]
  state_cols <- setdiff(names(z_wide), "row_label")
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
    colorbar = list(title = value_col)
  ) %>%
    layout(
      title = title,
      xaxis = list(title = "", tickangle = -35, automargin = TRUE),
      yaxis = list(title = "", automargin = TRUE),
      margin = list(l = 300, b = 150, t = 70, r = 30)
    )
  dominant_save_widget(p, output_file = output_file, title = title)
}

run_random_forest <- function(dt, state_cols, categorical_fields, numeric_fields) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    warning("randomForest package not installed; skipping RF model.")
    return(list(importance = data.table(), performance = data.table()))
  }
  rf_dt <- data.frame(row.names = seq_len(nrow(dt)))
  for (field in categorical_fields) {
    x <- collapse_levels(dt[[paste0(field, "_raw")]],
                         min_n = min_level_n,
                         max_levels = max_rf_factor_levels)
    rf_dt[[safe_name(field)]] <- factor(x)
  }
  for (field in numeric_fields) {
    x <- as.numeric(dt[[field]])
    med <- median(x[is.finite(x)], na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    x[!is.finite(x)] <- med
    rf_dt[[safe_name(field)]] <- x
  }
  keep <- vapply(rf_dt, function(x) uniqueN(x) > 1, logical(1))
  rf_dt <- rf_dt[, keep, drop = FALSE]
  if (ncol(rf_dt) == 0) {
    return(list(importance = data.table(), performance = data.table()))
  }

  rf_results <- lapply(state_cols, function(state) {
    y <- as.numeric(dt[[state]])
    valid <- is.finite(y)
    if (sum(valid) < 100) return(NULL)
    set.seed(100 + match(state, state_cols))
    fit <- randomForest::randomForest(
      x = rf_dt[valid, , drop = FALSE],
      y = y[valid],
      ntree = rf_ntree,
      importance = TRUE
    )
    imp <- as.data.table(randomForest::importance(fit), keep.rownames = "predictor")
    imp[, dominant_state := state]
    list(
      importance = imp,
      performance = data.table(
        dominant_state = state,
        n = sum(valid),
        ntree = rf_ntree,
        mse_final = tail(fit$mse, 1),
        variance_explained_final = tail(fit$rsq, 1)
      )
    )
  })
  rf_results <- Filter(Negate(is.null), rf_results)
  importance <- rbindlist(lapply(rf_results, `[[`, "importance"), fill = TRUE)
  performance <- rbindlist(lapply(rf_results, `[[`, "performance"), fill = TRUE)

  if (nrow(importance) > 0) {
    reverse_map <- data.table(
      predictor = safe_name(c(categorical_fields, numeric_fields)),
      field = c(categorical_fields, numeric_fields),
      predictor_type = c(rep("categorical", length(categorical_fields)),
                         rep("numeric", length(numeric_fields)))
    )
    importance <- merge(importance, reverse_map, by = "predictor",
                        all.x = TRUE, sort = FALSE)
    setcolorder(importance, c("dominant_state", "field", "predictor_type",
                              "predictor", setdiff(names(importance),
                                                    c("dominant_state", "field",
                                                      "predictor_type", "predictor"))))
  }
  list(importance = importance, performance = performance)
}

save_rf_heatmap <- function(rf_importance) {
  if (nrow(rf_importance) == 0 || !"%IncMSE" %chin% names(rf_importance)) {
    return(invisible(NULL))
  }
  plot_dt <- copy(rf_importance)
  plot_dt <- plot_dt[is.finite(`%IncMSE`)]
  field_order <- plot_dt[
    ,
    .(max_importance = max(`%IncMSE`, na.rm = TRUE)),
    by = field
  ][order(-max_importance)][seq_len(min(.N, 60L))]
  plot_dt <- plot_dt[field %chin% field_order$field]
  wide <- dcast(plot_dt, field ~ dominant_state, value.var = "%IncMSE")
  wide <- wide[field_order, on = "field"]
  states <- setdiff(names(wide), "field")
  z <- as.matrix(wide[, ..states])
  p <- plot_ly(
    x = states,
    y = wide$field,
    z = z,
    type = "heatmap",
    colorscale = "Viridis",
    colorbar = list(title = "%IncMSE")
  ) %>%
    layout(
      title = "Random-forest metadata importance for dominant-state scores",
      xaxis = list(title = "", tickangle = -35, automargin = TRUE),
      yaxis = list(title = "", automargin = TRUE),
      margin = list(l = 220, b = 150, t = 70, r = 30)
    )
  dominant_save_widget(
    p,
    output_file = rf_heatmap_output,
    title = "RF metadata importance"
  )
}

if (!file.exists(state_file)) stop("Missing ", state_file)
if (!file.exists(metadata_file)) stop("Missing ", metadata_file)

state_dt <- fread(state_file)
metadata <- fread(metadata_file)
matching_duplicate_runs <- metadata[
  Run %chin% state_dt$Run & (duplicated(Run) | duplicated(Run, fromLast = TRUE))
]
if (nrow(matching_duplicate_runs) > 0) {
  stop("Metadata has duplicate Run rows matching the dominant-state table: ",
       paste(unique(matching_duplicate_runs$Run), collapse = ", "))
}

state_dt[, input_order := .I]
dt <- merge(state_dt, metadata, by = "Run", all.x = TRUE, sort = FALSE)
setorder(dt, input_order)
dt[, input_order := NULL]

state_cols <- setdiff(
  names(state_dt),
  c("Run", "sample", "dominant_state", "dominant_state_score",
    "Baseline (control)")
)
state_cols <- state_cols[state_cols %chin% names(dt)]

metadata_cols <- setdiff(names(metadata), "Run")
field_profile <- select_metadata_fields(dt, names(metadata))
selected_fields <- field_profile[selected == TRUE, field]
selected_fields <- selected_fields[selected_fields %chin% names(dt)]
if (length(selected_fields) == 0) stop("No usable metadata fields selected.")

numeric_profile <- select_numeric_fields(dt, names(metadata))
numeric_fields <- numeric_profile[selected == TRUE, field]
numeric_fields <- numeric_fields[numeric_fields %chin% names(dt)]

for (field in selected_fields) {
  dt[, (paste0(field, "_raw")) := clean_level(get(field))]
  dt[, (field) := collapse_levels(get(field))]
}

field_profile[, used_for_broad_model := field %chin% selected_fields]
fwrite(field_profile, selected_fields_output)
fwrite(numeric_profile, numeric_fields_output)

raw_score_matrix <- as.matrix(dt[, ..state_cols])
storage.mode(raw_score_matrix) <- "numeric"

message("Selected categorical metadata fields: ",
        paste(selected_fields, collapse = ", "))
message("Selected numeric metadata fields: ",
        paste(numeric_fields, collapse = ", "))

main_terms <- rbindlist(lapply(selected_fields, function(field) {
  message("Scoring broad main terms for ", field)
  score_group_terms(
    dt,
    group_cols = field,
    score_matrix = raw_score_matrix,
    state_cols = state_cols,
    score_type = "raw_score",
    min_in = min_level_n,
    min_out = min_nwise_out,
    complete_terms_only = FALSE
  )
}), fill = TRUE)
if (nrow(main_terms) > 0) {
  main_terms[, `:=`(
    greater_p_adj = p.adjust(greater_p, "BH"),
    less_p_adj = p.adjust(less_p, "BH")
  )]
  setorder(main_terms, greater_p_adj, -effect_z, field, level)
}
fwrite(main_terms, main_terms_output)

residual_matrix <- build_residual_matrix(
  dt,
  state_cols = state_cols,
  categorical_fields = selected_fields,
  numeric_fields = numeric_fields
)

nwise_fields <- intersect(nwise_priority_fields, selected_fields)
if (length(nwise_fields) > max_nwise_fields) {
  nwise_fields <- nwise_fields[seq_len(max_nwise_fields)]
}
message("n-wise fields: ", paste(nwise_fields, collapse = ", "))
nwise_combos <- unlist(
  lapply(seq_len(max_interaction_order), function(k) {
    if (k < 2 || length(nwise_fields) < k) return(list())
    combn(nwise_fields, k, simplify = FALSE)
  }),
  recursive = FALSE
)

nwise_terms <- rbindlist(lapply(seq_along(nwise_combos), function(i) {
  fields <- nwise_combos[[i]]
  message("Scoring n-wise term ", i, "/", length(nwise_combos),
          ": ", paste(fields, collapse = " + "))
  score_group_terms(
    dt,
    group_cols = fields,
    score_matrix = residual_matrix,
    state_cols = state_cols,
    score_type = "main_effect_residual",
    min_in = min_nwise_n,
    min_out = min_nwise_out,
    complete_terms_only = TRUE
  )
}), fill = TRUE)
if (nrow(nwise_terms) > 0) {
  nwise_terms[, `:=`(
    greater_p_adj = p.adjust(greater_p, "BH"),
    less_p_adj = p.adjust(less_p, "BH")
  )]
  setorder(nwise_terms, greater_p_adj, -effect_z, -n_in, field, level)
}
fwrite(nwise_terms, nwise_terms_output)
save_state_heatmap(
  nwise_terms,
  nwise_heatmap_output,
  value_col = "effect_z",
  title = paste0("Broad n-wise metadata terms up to order ", max_interaction_order)
)

rf <- run_random_forest(
  dt,
  state_cols = state_cols,
  categorical_fields = selected_fields,
  numeric_fields = numeric_fields
)
fwrite(rf$importance, rf_importance_output)
fwrite(rf$performance, rf_performance_output)
save_rf_heatmap(rf$importance)

rdg_terms <- nwise_terms[
  is.finite(greater_p_adj) &
    greater_p_adj < 0.05 &
    effect_z > 0 &
    term_has_artificial_level == FALSE
][
  ,
  head(.SD, 4L),
  by = dominant_state
]
if (nrow(rdg_terms) > 0) {
  rdg_terms[, `:=`(
    adjusted_effect_z = effect_z,
    adjusted_score_delta = score_delta,
    adjusted_wilcox_greater_p_adj = greater_p_adj,
    broad_term_order = term_order,
    branch_source = "broad_metadata_nwise"
  )]
  rdg_terms <- rdg_terms[, .(
    field, level, dominant_state, n_in, n_out,
    adjusted_effect_z, adjusted_score_delta,
    adjusted_wilcox_greater_p_adj, broad_term_order, branch_source
  )]
}
fwrite(rdg_terms, rdg_terms_output)

compact_summary <- rbindlist(list(
  main_terms[
    is.finite(greater_p_adj) & greater_p_adj < 0.05 & effect_z > 0
  ][
    ,
    head(.SD, 10L),
    by = dominant_state
  ][
    ,
    .(summary_type = "main_term", dominant_state, field, level, term_order,
      n_in, effect_z, q = greater_p_adj)
  ],
  nwise_terms[
    is.finite(greater_p_adj) & greater_p_adj < 0.05 & effect_z > 0
  ][
    ,
    head(.SD, 10L),
    by = dominant_state
  ][
    ,
    .(summary_type = "nwise_residual", dominant_state, field, level, term_order,
      n_in, effect_z, q = greater_p_adj)
  ]
), fill = TRUE)
setorder(compact_summary, summary_type, dominant_state, q, -effect_z)
fwrite(compact_summary, compact_summary_output)

message("Saved: ", selected_fields_output)
message("Saved: ", numeric_fields_output)
message("Saved: ", main_terms_output)
message("Saved: ", nwise_terms_output)
message("Saved: ", rf_importance_output)
message("Saved: ", rf_performance_output)
message("Saved: ", rdg_terms_output)
message("Saved: ", compact_summary_output)

print(field_profile[used_for_broad_model == TRUE, .(
  field, non_missing, non_missing_fraction, raw_unique_n
)])
print(head(compact_summary, 40))
print(rf$performance)

invisible(list(
  selected_fields = field_profile,
  numeric_fields = numeric_profile,
  main_terms = main_terms,
  nwise_terms = nwise_terms,
  rf_importance = rf$importance,
  rf_performance = rf$performance,
  rdg_terms = rdg_terms,
  compact_summary = compact_summary
))
