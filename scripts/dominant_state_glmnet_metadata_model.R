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
  library(Matrix)
  library(glmnet)
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
metadata_file <- "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv"
selected_fields_file <- result_file("dominant_state_broad_metadata_selected_fields.csv")
numeric_fields_file <- result_file("dominant_state_broad_metadata_numeric_fields.csv")
output_prefix <- "dominant_state_glmnet_metadata"

feature_metadata_output <- result_file(paste0(output_prefix, "_feature_metadata.csv"))
coefficients_output <- result_file(paste0(output_prefix, "_coefficients.csv"))
performance_output <- result_file(paste0(output_prefix, "_performance.csv"))
stability_output <- result_file(paste0(output_prefix, "_stability_selection.csv"))
selected_term_scores_output <- result_file(paste0(output_prefix, "_selected_term_scores.csv"))
rdg_terms_output <- result_file(paste0(output_prefix, "_rdg_terms.csv"))
compact_summary_output <- result_file(paste0(output_prefix, "_compact_summary.csv"))
coefficient_heatmap_output <- result_file(paste0(output_prefix, "_coefficient_heatmap.html"))

min_level_n <- 8L
min_interaction_n <- 8L
min_out_n <- 30L
max_main_levels <- 50L
max_interaction_levels <- 80L
max_numeric_fields <- 50L
nfolds <- 5L
alpha <- 1
stability_repeats <- 12L
stability_subsample_fraction <- 0.80
max_heatmap_rows <- 120L
max_rdg_terms_per_state <- 4L
min_rdg_effect_z <- 0.5

biological_fields <- c(
  "Cancer_type", "Cell_model", "Cell_type", "Organ_system", "Sex",
  "Life_stage", "TISSUE", "CELL_LINE"
)
perturbation_fields <- c("CONDITION", "GENE", "INHIBITOR", "FRACTION")
nuisance_fields <- c(
  "study", "AUTHOR", "LibrarySelection", "Model", "LibraryLayout",
  "Library_Strandedness", "LibraryStrategy", "BATCH", "TIMEPOINT"
)

interaction_pairs <- list(
  c("CONDITION", "CELL_LINE"),
  c("CONDITION", "TISSUE"),
  c("CONDITION", "GENE"),
  c("CONDITION", "Cancer_type"),
  c("CONDITION", "Sex"),
  c("CONDITION", "Cell_model"),
  c("CONDITION", "Cell_type"),
  c("CONDITION", "Organ_system"),
  c("CONDITION", "INHIBITOR"),
  c("CELL_LINE", "TISSUE"),
  c("CELL_LINE", "GENE"),
  c("CELL_LINE", "Cancer_type"),
  c("CELL_LINE", "Sex"),
  c("CELL_LINE", "Cell_model"),
  c("TISSUE", "Cancer_type"),
  c("TISSUE", "Sex"),
  c("TISSUE", "Cell_type"),
  c("Cancer_type", "Sex"),
  c("Cancer_type", "Cell_model"),
  c("Cancer_type", "Organ_system"),
  c("Sex", "Cell_model"),
  c("Sex", "Cell_type"),
  c("GENE", "INHIBITOR")
)

clean_level <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "NA", "N/A", "na", "n/a", "NULL", "null", "None",
               "none", "missing", "Missing", "unknown", "Unknown")] <- NA_character_
  x
}

file_safe <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  substr(x, 1L, 160L)
}

collapse_levels <- function(x, min_n = min_level_n, max_levels = max_main_levels) {
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

field_role <- function(fields) {
  if (any(fields %chin% nuisance_fields)) return("nuisance_or_design")
  if (any(fields %chin% perturbation_fields) && any(fields %chin% biological_fields)) {
    return("biology_perturbation")
  }
  if (all(fields %chin% perturbation_fields)) return("perturbation")
  if (all(fields %chin% biological_fields)) return("biology")
  "metadata_other"
}

penalty_for_feature <- function(feature_role, feature_type, term_order) {
  ifelse(
    feature_role == "nuisance_or_design",
    0,
    ifelse(
      feature_type == "numeric_qc",
      0.25,
      ifelse(term_order >= 2L, 1.25, 1)
    )
  )
}

save_widget <- function(widget, output_file, title) {
  dominant_save_widget(widget, output_file = output_file, title = title)
}

read_selected_fields <- function(metadata) {
  if (file.exists(selected_fields_file)) {
    profile <- fread(selected_fields_file)
    if (all(c("field", "used_for_broad_model") %chin% names(profile))) {
      fields <- profile[used_for_broad_model == TRUE, field]
      fields <- fields[fields %chin% names(metadata)]
      if (length(fields) > 0) return(fields)
    }
  }
  intersect(c(biological_fields, perturbation_fields, nuisance_fields), names(metadata))
}

read_numeric_fields <- function(metadata) {
  if (!file.exists(numeric_fields_file)) return(character())
  profile <- fread(numeric_fields_file)
  if (!all(c("field", "selected") %chin% names(profile))) return(character())
  fields <- profile[selected == TRUE, field]
  fields <- fields[fields %chin% names(metadata)]
  fields[seq_len(min(length(fields), max_numeric_fields))]
}

make_feature_id <- function(type, field, level) {
  paste(type, file_safe(field), file_safe(level), sep = "|")
}

build_sparse_features <- function(dt, categorical_fields, numeric_fields) {
  n <- nrow(dt)
  model_dt <- data.table(row_id = seq_len(n))
  for (field in categorical_fields) {
    max_levels <- if (field %chin% nuisance_fields) 40L else max_main_levels
    model_dt[, (field) := collapse_levels(dt[[field]], max_levels = max_levels)]
  }

  row_idx <- integer()
  col_idx <- integer()
  feature_rows <- list()
  feature_counter <- 0L

  add_binary_feature <- function(rows, feature_id, feature_type, fields,
                                 level, term_order, role,
                                 artificial = FALSE) {
    if (length(rows) < min_level_n || n - length(rows) < min_out_n) return(NULL)
    feature_counter <<- feature_counter + 1L
    row_idx <<- c(row_idx, rows)
    col_idx <<- c(col_idx, rep(feature_counter, length(rows)))
    feature_rows[[feature_counter]] <<- data.table(
      feature_id = feature_id,
      feature_index = feature_counter,
      feature_type = feature_type,
      feature_role = role,
      term_order = term_order,
      field = paste(fields, collapse = ":"),
      level = level,
      n_in = length(rows),
      n_out = n - length(rows),
      term_has_artificial_level = artificial
    )
    NULL
  }

  for (field in categorical_fields) {
    values <- model_dt[[field]]
    levels <- sort(unique(values))
    for (level in levels) {
      rows <- which(values == level)
      add_binary_feature(
        rows = rows,
        feature_id = make_feature_id("main", field, level),
        feature_type = "categorical_main",
        fields = field,
        level = level,
        term_order = 1L,
        role = field_role(field),
        artificial = is_artificial_level(level)
      )
    }
  }

  valid_pairs <- Filter(function(pair) all(pair %chin% categorical_fields),
                        interaction_pairs)
  for (pair in valid_pairs) {
    f1 <- pair[[1]]
    f2 <- pair[[2]]
    v1 <- model_dt[[f1]]
    v2 <- model_dt[[f2]]
    usable <- !is_artificial_level(v1) & !is_artificial_level(v2)
    if (!any(usable)) next
    term <- paste0(f1, "=", v1, " | ", f2, "=", v2)
    tab <- sort(table(term[usable]), decreasing = TRUE)
    keep <- names(tab)[tab >= min_interaction_n]
    if (length(keep) > max_interaction_levels) keep <- keep[seq_len(max_interaction_levels)]
    for (level in keep) {
      rows <- which(usable & term == level)
      add_binary_feature(
        rows = rows,
        feature_id = make_feature_id("interaction", paste(pair, collapse = ":"), level),
        feature_type = "categorical_interaction",
        fields = pair,
        level = level,
        term_order = 2L,
        role = field_role(pair),
        artificial = FALSE
      )
    }
  }

  feature_meta <- rbindlist(feature_rows, fill = TRUE)
  if (nrow(feature_meta) == 0) {
    stop("No categorical features passed filters.")
  }
  X_cat <- sparseMatrix(
    i = row_idx,
    j = col_idx,
    x = 1,
    dims = c(n, nrow(feature_meta))
  )
  colnames(X_cat) <- feature_meta$feature_id

  numeric_meta <- data.table()
  if (length(numeric_fields) > 0) {
    numeric_values <- lapply(numeric_fields, function(field) {
      x <- as.numeric(dt[[field]])
      finite <- is.finite(x)
      if (sum(finite) < min_level_n || uniqueN(x[finite]) < 3) return(NULL)
      med <- median(x[finite], na.rm = TRUE)
      x[!finite] <- med
      x
    })
    names(numeric_values) <- numeric_fields
    keep_numeric <- !vapply(numeric_values, is.null, logical(1))
    numeric_values <- numeric_values[keep_numeric]
    if (length(numeric_values) > 0) {
      X_num <- Matrix(do.call(cbind, numeric_values), sparse = TRUE)
      numeric_ids <- make_feature_id("numeric", names(numeric_values), names(numeric_values))
      colnames(X_num) <- numeric_ids
      numeric_meta <- data.table(
        feature_id = numeric_ids,
        feature_index = nrow(feature_meta) + seq_along(numeric_ids),
        feature_type = "numeric_qc",
        feature_role = "technical_qc",
        term_order = 1L,
        field = names(numeric_values),
        level = names(numeric_values),
        n_in = sum(is.finite(unlist(numeric_values[[1]], use.names = FALSE))),
        n_out = NA_integer_,
        term_has_artificial_level = FALSE
      )
      feature_meta <- rbind(feature_meta, numeric_meta, fill = TRUE)
      X <- cbind(X_cat, X_num)
    } else {
      X <- X_cat
    }
  } else {
    X <- X_cat
  }

  if (anyDuplicated(feature_meta$feature_id)) {
    feature_meta[, feature_id := make.unique(feature_id)]
    colnames(X) <- feature_meta$feature_id
  }
  feature_meta[, penalty_factor := penalty_for_feature(
    feature_role, feature_type, term_order
  )]
  list(X = X, feature_meta = feature_meta, model_dt = model_dt)
}

coefs_to_dt <- function(fit, state, lambda_choice) {
  mat <- coef(fit, s = lambda_choice)
  vals <- as.numeric(mat)
  names(vals) <- rownames(mat)
  vals <- vals[vals != 0]
  vals <- vals[names(vals) != "(Intercept)"]
  if (length(vals) == 0) return(data.table())
  data.table(
    dominant_state = state,
    lambda_choice = lambda_choice,
    feature_id = names(vals),
    coefficient = unname(vals)
  )
}

score_selected_feature <- function(model_dt, y, feature_row) {
  if (feature_row$feature_type == "numeric_qc") return(NULL)
  if (feature_row$term_order == 1L) {
    in_group <- model_dt[[feature_row$field]] == feature_row$level
  } else {
    fields <- strsplit(feature_row$field, ":", fixed = TRUE)[[1]]
    parts <- strsplit(feature_row$level, " | ", fixed = TRUE)[[1]]
    values <- rep(NA_character_, length(fields))
    names(values) <- fields
    for (part in parts) {
      eq <- regexpr("=", part, fixed = TRUE)[1]
      if (!is.finite(eq) || eq < 1L) next
      key <- substr(part, 1L, eq - 1L)
      value <- substr(part, eq + 1L, nchar(part))
      if (key %chin% fields) values[[key]] <- value
    }
    if (anyNA(values) || !all(fields %chin% names(model_dt))) return(NULL)
    in_group <- rep(TRUE, nrow(model_dt))
    for (field in fields) {
      in_group <- in_group & model_dt[[field]] == values[[field]]
    }
  }
  in_group[is.na(in_group)] <- FALSE
  valid <- is.finite(y)
  n_in <- sum(valid & in_group)
  n_out <- sum(valid & !in_group)
  if (n_in < min_level_n || n_out < min_out_n) return(NULL)
  y_in <- y[valid & in_group]
  y_out <- y[valid & !in_group]
  total_sd <- sd(y[valid])
  if (!is.finite(total_sd) || total_sd <= 0) return(NULL)
  var_in <- var(y_in)
  var_out <- var(y_out)
  se <- sqrt(var_in / n_in + var_out / n_out)
  t_stat <- if (is.finite(se) && se > 0) {
    (mean(y_in) - mean(y_out)) / se
  } else {
    NA_real_
  }
  df <- {
    numerator <- (var_in / n_in + var_out / n_out)^2
    denominator <- ((var_in / n_in)^2 / pmax(n_in - 1, 1)) +
      ((var_out / n_out)^2 / pmax(n_out - 1, 1))
    if (is.finite(denominator) && denominator > 0) numerator / denominator else NA_real_
  }
  data.table(
    n_in = n_in,
    n_out = n_out,
    mean_in = mean(y_in),
    mean_out = mean(y_out),
    score_delta = mean(y_in) - mean(y_out),
    effect_z = (mean(y_in) - mean(y_out)) / total_sd,
    t_stat = t_stat,
    df = df,
    greater_p = pt(t_stat, df = df, lower.tail = FALSE),
    less_p = pt(t_stat, df = df, lower.tail = TRUE)
  )
}

make_coefficient_heatmap <- function(coefs) {
  plot_dt <- coefs[
    lambda_choice == "lambda.1se" &
      feature_type != "numeric_qc" &
      feature_role != "nuisance_or_design"
  ]
  if (nrow(plot_dt) == 0) return(invisible(NULL))
  plot_dt[, abs_coef := abs(coefficient)]
  row_order <- plot_dt[
    ,
    .(max_abs_coef = max(abs_coef, na.rm = TRUE),
      n_states = uniqueN(dominant_state)),
    by = .(field, level)
  ][order(-max_abs_coef, -n_states)]
  row_order <- row_order[seq_len(min(.N, max_heatmap_rows))]
  plot_dt <- plot_dt[row_order, on = c("field", "level"), nomatch = 0]
  plot_dt[, row_label := paste(field, level, sep = " :: ")]
  row_order[, row_label := paste(field, level, sep = " :: ")]
  plot_dt[, hover := paste0(
    "Feature: ", row_label,
    "<br>State: ", dominant_state,
    "<br>role: ", feature_role,
    "<br>coef: ", signif(coefficient, 3),
    "<br>n: ", n_in
  )]
  z_wide <- dcast(plot_dt, row_label ~ dominant_state, value.var = "coefficient")
  text_wide <- dcast(plot_dt, row_label ~ dominant_state, value.var = "hover")
  z_wide <- z_wide[row_order[, .(row_label)], on = "row_label"]
  text_wide <- text_wide[row_order[, .(row_label)], on = "row_label"]
  state_cols <- setdiff(names(z_wide), "row_label")
  z <- as.matrix(z_wide[, ..state_cols])
  text <- as.matrix(text_wide[, ..state_cols])
  max_abs <- suppressWarnings(quantile(abs(z), 0.98, na.rm = TRUE))
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
  p <- plot_ly(
    x = state_cols,
    y = z_wide$row_label,
    z = z,
    text = text,
    hoverinfo = "text",
    type = "heatmap",
    zmin = -max_abs,
    zmax = max_abs,
    colorscale = list(c(0, "#2166AC"), c(0.5, "#F7F7F7"), c(1, "#B2182B")),
    colorbar = list(title = "coef")
  ) %>%
    layout(
      title = "Elastic-net metadata coefficients at lambda.1se",
      xaxis = list(title = "", tickangle = -35, automargin = TRUE),
      yaxis = list(title = "", automargin = TRUE),
      margin = list(l = 320, b = 150, t = 70, r = 30)
    )
  save_widget(p, coefficient_heatmap_output, "glmnet metadata coefficient heatmap")
}

run_stability_selection <- function(fits, X, dt, state_cols, feature_meta,
                                    penalty_factor) {
  if (stability_repeats <= 0) return(data.table())
  rbindlist(lapply(state_cols, function(state) {
    fit <- fits[[state]]
    if (is.null(fit)) return(NULL)
    y <- as.numeric(dt[[state]])
    valid <- which(is.finite(y))
    if (length(valid) < 100) return(NULL)
    lambda <- fit$lambda.1se
    repeat_rows <- lapply(seq_len(stability_repeats), function(rep_id) {
      set.seed(10000 + rep_id * 37 + match(state, state_cols))
      n_sub <- max(50L, floor(length(valid) * stability_subsample_fraction))
      sub <- sort(sample(valid, n_sub))
      sub_fit <- tryCatch(
        glmnet(
          x = X[sub, , drop = FALSE],
          y = y[sub],
          family = "gaussian",
          alpha = alpha,
          lambda = lambda,
          standardize = TRUE,
          intercept = TRUE,
          penalty.factor = penalty_factor
        ),
        error = function(e) NULL
      )
      if (is.null(sub_fit)) return(NULL)
      mat <- coef(sub_fit, s = lambda)
      vals <- as.numeric(mat)
      names(vals) <- rownames(mat)
      vals <- vals[names(vals) != "(Intercept)"]
      data.table(
        dominant_state = state,
        repeat_id = rep_id,
        feature_id = names(vals),
        coefficient = vals,
        selected = vals != 0
      )
    })
    rep_dt <- rbindlist(repeat_rows, fill = TRUE)
    if (nrow(rep_dt) == 0) return(NULL)
    rep_dt[
      ,
      .(
        stability_fraction = mean(selected),
        median_coefficient = median(coefficient),
        mean_abs_coefficient = mean(abs(coefficient))
      ),
      by = .(dominant_state, feature_id)
    ]
  }), fill = TRUE)[
    feature_meta,
    on = "feature_id"
  ][
    ,
    stable_positive := stability_fraction >= 0.5 & median_coefficient > 0
  ][
    order(dominant_state, -stability_fraction, -mean_abs_coefficient)
  ][]
}

if (!file.exists(state_file)) stop("Missing ", state_file)
if (!file.exists(metadata_file)) stop("Missing ", metadata_file)

state_dt <- fread(state_file)
metadata <- fread(metadata_file)
categorical_fields <- read_selected_fields(metadata)
numeric_fields <- read_numeric_fields(metadata)

state_dt[, input_order := .I]
keep_metadata <- intersect(c("Run", categorical_fields, numeric_fields), names(metadata))
dt <- merge(state_dt, metadata[, ..keep_metadata], by = "Run", all.x = TRUE, sort = FALSE)
setorder(dt, input_order)
dt[, input_order := NULL]

state_cols <- setdiff(
  names(state_dt),
  c("Run", "sample", "dominant_state", "dominant_state_score",
    "Baseline (control)")
)
state_cols <- state_cols[state_cols %chin% names(dt)]
if (length(state_cols) < 2) stop("Need at least two state score columns.")

message("Building sparse metadata design:")
message("  categorical fields: ", paste(categorical_fields, collapse = ", "))
message("  numeric fields: ", paste(numeric_fields, collapse = ", "))
design <- build_sparse_features(dt, categorical_fields, numeric_fields)
X <- design$X
feature_meta <- design$feature_meta
model_dt <- design$model_dt
fwrite(feature_meta, feature_metadata_output)
message("  samples: ", nrow(X), ", features: ", ncol(X))
penalty_factor <- feature_meta$penalty_factor

set.seed(7)
foldid <- sample(rep(seq_len(nfolds), length.out = nrow(dt)))

fits <- vector("list", length(state_cols))
names(fits) <- state_cols
performance <- rbindlist(lapply(state_cols, function(state) {
  message("Fitting elastic-net model for ", state)
  y <- as.numeric(dt[[state]])
  valid <- is.finite(y)
  if (sum(valid) < 100) return(NULL)
  fit <- cv.glmnet(
    x = X[valid, , drop = FALSE],
    y = y[valid],
    family = "gaussian",
    alpha = alpha,
    nfolds = nfolds,
    foldid = foldid[valid],
    standardize = TRUE,
    intercept = TRUE,
    penalty.factor = penalty_factor
  )
  fits[[state]] <<- fit
  null_mse <- mean((y[valid] - mean(y[valid]))^2)
  idx_min <- match(fit$lambda.min, fit$lambda)
  idx_1se <- match(fit$lambda.1se, fit$lambda)
  coef_min <- coefs_to_dt(fit, state, "lambda.min")
  coef_1se <- coefs_to_dt(fit, state, "lambda.1se")
  data.table(
    dominant_state = state,
    n = sum(valid),
    alpha = alpha,
    nfolds = nfolds,
    n_features = ncol(X),
    lambda_min = fit$lambda.min,
    lambda_1se = fit$lambda.1se,
    cv_mse_min = fit$cvm[idx_min],
    cv_mse_1se = fit$cvm[idx_1se],
    cv_se_min = fit$cvsd[idx_min],
    cv_se_1se = fit$cvsd[idx_1se],
    null_mse = null_mse,
    cv_r2_min = 1 - fit$cvm[idx_min] / null_mse,
    cv_r2_1se = 1 - fit$cvm[idx_1se] / null_mse,
    n_nonzero_min = nrow(coef_min),
    n_nonzero_1se = nrow(coef_1se)
  )
}), fill = TRUE)
fwrite(performance, performance_output)

coefficients <- rbindlist(lapply(names(fits), function(state) {
  fit <- fits[[state]]
  if (is.null(fit)) return(NULL)
  rbind(
    coefs_to_dt(fit, state, "lambda.min"),
    coefs_to_dt(fit, state, "lambda.1se"),
    fill = TRUE
  )
}), fill = TRUE)
if (nrow(coefficients) > 0) {
  coefficients <- merge(coefficients, feature_meta, by = "feature_id",
                        all.x = TRUE, sort = FALSE)
  coefficients[, `:=`(
    coefficient_abs = abs(coefficient),
    coefficient_direction = fifelse(coefficient > 0, "positive", "negative")
  )]
  setcolorder(coefficients, c(
    "dominant_state", "lambda_choice", "feature_id", "feature_type",
    "feature_role", "term_order", "field", "level", "coefficient",
    "coefficient_abs", "coefficient_direction",
    setdiff(names(coefficients), c(
      "dominant_state", "lambda_choice", "feature_id", "feature_type",
      "feature_role", "term_order", "field", "level", "coefficient",
      "coefficient_abs", "coefficient_direction"
    ))
  ))
  setorder(coefficients, lambda_choice, dominant_state, -coefficient_abs)
}
fwrite(coefficients, coefficients_output)
make_coefficient_heatmap(coefficients)

stability <- run_stability_selection(
  fits, X, dt, state_cols, feature_meta, penalty_factor
)
fwrite(stability, stability_output)

selected <- coefficients[
  lambda_choice == "lambda.1se" &
    feature_type != "numeric_qc"
]
selected_scores <- rbindlist(lapply(seq_len(nrow(selected)), function(i) {
  row <- copy(selected[i])
  y <- as.numeric(dt[[row$dominant_state]])
  score <- score_selected_feature(model_dt, y, row)
  if (is.null(score)) return(NULL)
  if (all(c("n_in", "n_out") %chin% names(row))) {
    row[, `:=`(
      feature_n_in = n_in,
      feature_n_out = n_out,
      n_in = NULL,
      n_out = NULL
    )]
  }
  cbind(row, score)
}), fill = TRUE)
if (nrow(selected_scores) > 0) {
  if (nrow(stability) > 0) {
    selected_scores <- merge(
      selected_scores,
      stability[, .(dominant_state, feature_id, stability_fraction,
                    median_coefficient, stable_positive)],
      by = c("dominant_state", "feature_id"),
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    selected_scores[, `:=`(
      stability_fraction = NA_real_,
      median_coefficient = NA_real_,
      stable_positive = NA
    )]
  }
  selected_scores[, `:=`(
    greater_p_adj = p.adjust(greater_p, "BH"),
    less_p_adj = p.adjust(less_p, "BH")
  )]
  setorder(selected_scores, dominant_state, -coefficient_abs, greater_p_adj)
}
fwrite(selected_scores, selected_term_scores_output)

rdg_terms <- selected_scores[
  coefficient > 0 &
    effect_z >= min_rdg_effect_z &
    (is.na(stability_fraction) | stability_fraction >= 0.4) &
    is.finite(greater_p_adj) &
    greater_p_adj < 0.05 &
    feature_role %chin% c("biology", "perturbation", "biology_perturbation") &
    term_has_artificial_level == FALSE
][
  order(dominant_state, -coefficient_abs, greater_p_adj),
  head(.SD, max_rdg_terms_per_state),
  by = dominant_state
]
if (nrow(rdg_terms) > 0) {
  rdg_terms[, `:=`(
    adjusted_effect_z = effect_z,
    adjusted_score_delta = score_delta,
    adjusted_wilcox_greater_p_adj = greater_p_adj,
    glmnet_lambda_choice = lambda_choice,
    glmnet_coefficient = coefficient,
    stability_fraction = stability_fraction,
    branch_source = "glmnet_regularized_metadata"
  )]
  rdg_terms <- rdg_terms[, .(
    field, level, dominant_state, n_in, n_out,
    adjusted_effect_z, adjusted_score_delta, adjusted_wilcox_greater_p_adj,
    glmnet_lambda_choice, glmnet_coefficient, stability_fraction,
    feature_role, branch_source
  )]
}
fwrite(rdg_terms, rdg_terms_output)

compact_summary <- rbindlist(list(
  performance[
    ,
    .(
      summary_type = "performance",
      dominant_state,
      field = NA_character_,
      level = NA_character_,
      feature_role = NA_character_,
      lambda_choice = "lambda.1se",
      coefficient = NA_real_,
      effect_z = NA_real_,
      q = NA_real_,
      metric = paste0(
        "cv_r2_1se=", sprintf("%.3f", cv_r2_1se),
        "; nonzero=", n_nonzero_1se
      )
    )
  ],
  selected_scores[
    lambda_choice == "lambda.1se" &
      feature_role != "nuisance_or_design"
  ][
    order(dominant_state, -coefficient_abs),
    head(.SD, 12L),
    by = dominant_state
  ][
    ,
    .(
      summary_type = "selected_term",
      dominant_state,
      field,
      level,
      feature_role,
      lambda_choice,
      coefficient,
      effect_z,
      q = greater_p_adj,
      stability_fraction,
      metric = paste0("n=", n_in, "; delta=", sprintf("%.3f", score_delta))
    )
  ]
), fill = TRUE)
fwrite(compact_summary, compact_summary_output)

message("Saved: ", feature_metadata_output)
message("Saved: ", coefficients_output)
message("Saved: ", performance_output)
message("Saved: ", stability_output)
message("Saved: ", selected_term_scores_output)
message("Saved: ", rdg_terms_output)
message("Saved: ", compact_summary_output)

print(performance)
print(
  coefficients[
    lambda_choice == "lambda.1se",
    .N,
    by = .(dominant_state, feature_type, feature_role)
  ][order(dominant_state, feature_type, feature_role)]
)
print(head(compact_summary[summary_type == "selected_term"], 60))

invisible(list(
  feature_metadata = feature_meta,
  coefficients = coefficients,
  performance = performance,
  selected_term_scores = selected_scores,
  rdg_terms = rdg_terms,
  compact_summary = compact_summary
))
