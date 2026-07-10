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
metadata_file <- "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv"
selected_fields_file <- result_file("dominant_state_broad_metadata_selected_fields.csv")
numeric_fields_file <- result_file("dominant_state_broad_metadata_numeric_fields.csv")

pca_samples_output <- result_file("dominant_state_latent_program_pca_samples.csv")
pca_loadings_output <- result_file("dominant_state_latent_program_pca_loadings.csv")
pca_top_loadings_output <- result_file("dominant_state_latent_program_pca_top_loadings.csv")
pca_variance_output <- result_file("dominant_state_latent_program_pca_variance.csv")
pca_extreme_samples_output <- result_file("dominant_state_latent_program_pca_extreme_samples.csv")
pc_metadata_output <- result_file("dominant_state_latent_program_pc_metadata_associations.csv")
pc_numeric_output <- result_file("dominant_state_latent_program_pc_numeric_correlations.csv")
state_correlations_output <- result_file("dominant_state_latent_program_state_correlations.csv")
cluster_samples_output <- result_file("dominant_state_latent_program_clusters.csv")
cluster_state_summary_output <- result_file("dominant_state_latent_program_cluster_state_summary.csv")
cluster_metadata_output <- result_file("dominant_state_latent_program_cluster_metadata_enrichment.csv")
compact_report_output <- result_file("dominant_state_latent_program_compact_report.csv")
pc_metadata_heatmap_output <- result_file("dominant_state_latent_program_pc_metadata_heatmap.html")
state_correlation_heatmap_output <- result_file("dominant_state_latent_program_state_correlation_heatmap.html")
pca_scatter_output <- result_file("dominant_state_latent_program_pca_scatter.html")

n_pcs_to_save <- 8L
n_clusters <- 10L
min_level_n <- 8L
min_out_n <- 30L
min_cluster_level_n <- 5L
top_loadings_per_pc <- 4L
top_extreme_samples_per_pc <- 25L
max_pc_heatmap_rows <- 120L

fallback_metadata_fields <- c("CELL_LINE", "TISSUE", "GENE", "CONDITION")
priority_metadata_fields <- c(
  "Cancer_type", "Cell_model", "Cell_type", "Organ_system", "Sex",
  "Life_stage", "CONDITION", "TISSUE", "CELL_LINE", "GENE", "INHIBITOR",
  "FRACTION", "study", "AUTHOR"
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

safe_cor_test <- function(x, y, method) {
  valid <- is.finite(x) & is.finite(y)
  if (sum(valid) < 8L || uniqueN(x[valid]) < 3L || uniqueN(y[valid]) < 3L) {
    return(list(estimate = NA_real_, p = NA_real_, n = sum(valid)))
  }
  test <- tryCatch(
    suppressWarnings(cor.test(x[valid], y[valid], method = method)),
    error = function(e) NULL
  )
  if (is.null(test)) {
    return(list(estimate = NA_real_, p = NA_real_, n = sum(valid)))
  }
  list(estimate = unname(test$estimate), p = test$p.value, n = sum(valid))
}

safe_fisher_greater <- function(in_level, out_level, in_other, out_other) {
  mat <- matrix(c(in_level, in_other, out_level, out_other), nrow = 2,
                byrow = TRUE)
  if (sum(mat[1, ]) < min_cluster_level_n || sum(mat[, 1]) < min_cluster_level_n) {
    return(list(p = NA_real_, odds_ratio = NA_real_))
  }
  test <- tryCatch(fisher.test(mat, alternative = "greater"),
                   error = function(e) NULL)
  if (is.null(test)) return(list(p = NA_real_, odds_ratio = NA_real_))
  list(p = test$p.value, odds_ratio = unname(test$estimate))
}

save_widget <- function(widget, output_file, title) {
  dominant_save_widget(widget, output_file = output_file, title = title)
}

read_selected_metadata_fields <- function(metadata) {
  if (file.exists(selected_fields_file)) {
    profile <- fread(selected_fields_file)
    if (all(c("field", "used_for_broad_model") %chin% names(profile))) {
      fields <- profile[used_for_broad_model == TRUE, field]
      fields <- unique(c(intersect(priority_metadata_fields, fields), fields))
      fields <- fields[fields %chin% names(metadata)]
      if (length(fields) > 0) return(fields)
    }
  }
  intersect(fallback_metadata_fields, names(metadata))
}

read_selected_numeric_fields <- function(metadata) {
  if (!file.exists(numeric_fields_file)) return(character())
  profile <- fread(numeric_fields_file)
  if (!all(c("field", "selected") %chin% names(profile))) return(character())
  fields <- profile[selected == TRUE, field]
  fields[fields %chin% names(metadata)]
}

score_pc_metadata <- function(dt, metadata_fields, pc_cols) {
  rbindlist(lapply(metadata_fields, function(field) {
    values <- clean_level(dt[[field]])
    levels <- sort(unique(values[!is.na(values)]))
    levels <- levels[vapply(levels, function(level) {
      sum(values == level, na.rm = TRUE) >= min_level_n
    }, logical(1))]
    if (length(levels) == 0) return(NULL)
    rbindlist(lapply(levels, function(level) {
      in_level <- values == level
      in_level[is.na(in_level)] <- FALSE
      rbindlist(lapply(pc_cols, function(pc) {
        y <- as.numeric(dt[[pc]])
        valid <- is.finite(y)
        n_in <- sum(valid & in_level)
        n_out <- sum(valid & !in_level)
        if (n_in < min_level_n || n_out < min_out_n) return(NULL)
        y_in <- y[valid & in_level]
        y_out <- y[valid & !in_level]
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
          if (is.finite(denominator) && denominator > 0) {
            numerator / denominator
          } else {
            NA_real_
          }
        }
        data.table(
          field = field,
          level = level,
          PC = pc,
          n_in = n_in,
          n_out = n_out,
          mean_in = mean(y_in),
          mean_out = mean(y_out),
          pc_delta = mean(y_in) - mean(y_out),
          effect_z = (mean(y_in) - mean(y_out)) / total_sd,
          t_stat = t_stat,
          df = df,
          greater_p = pt(t_stat, df = df, lower.tail = FALSE),
          less_p = pt(t_stat, df = df, lower.tail = TRUE)
        )
      }), fill = TRUE)
    }), fill = TRUE)
  }), fill = TRUE)
}

score_cluster_metadata <- function(clusters, metadata_fields) {
  rbindlist(lapply(metadata_fields, function(field) {
    field_values <- clean_level(clusters[[field]])
    levels <- sort(unique(field_values[!is.na(field_values)]))
    levels <- levels[vapply(levels, function(level) {
      sum(field_values == level, na.rm = TRUE) >= min_cluster_level_n
    }, logical(1))]
    if (length(levels) == 0) return(NULL)
    rbindlist(lapply(levels, function(level) {
      in_level_all <- field_values == level
      in_level_all[is.na(in_level_all)] <- FALSE
      rbindlist(lapply(sort(unique(clusters$latent_cluster)), function(cluster) {
        in_cluster <- clusters$latent_cluster == cluster
        level_in <- sum(in_cluster & in_level_all)
        other_in <- sum(in_cluster & !in_level_all)
        level_out <- sum(!in_cluster & in_level_all)
        other_out <- sum(!in_cluster & !in_level_all)
        fish <- safe_fisher_greater(level_in, level_out, other_in, other_out)
        data.table(
          latent_cluster = cluster,
          field = field,
          level = level,
          n_level_in_cluster = level_in,
          n_cluster = sum(in_cluster),
          n_level_total = sum(in_level_all),
          fraction_in_cluster = level_in / sum(in_cluster),
          fraction_global = sum(in_level_all) / length(in_level_all),
          odds_ratio = fish$odds_ratio,
          fisher_greater_p = fish$p
        )
      }), fill = TRUE)
    }), fill = TRUE)
  }), fill = TRUE)
}

make_pc_metadata_heatmap <- function(pc_metadata) {
  plot_dt <- pc_metadata[
    is.finite(greater_p_adj) &
      greater_p_adj < 0.1 &
      is.finite(effect_z)
  ]
  if (nrow(plot_dt) == 0) return(invisible(NULL))
  plot_dt[, row_label := paste(field, level, sep = " :: ")]
  row_order <- plot_dt[
    ,
    .(
      best_q = min(greater_p_adj, na.rm = TRUE),
      best_abs = max(abs(effect_z), na.rm = TRUE),
      n = max(n_in, na.rm = TRUE)
    ),
    by = row_label
  ][order(best_q, -best_abs, -n)]
  row_order <- row_order[seq_len(min(.N, max_pc_heatmap_rows))]
  plot_dt <- plot_dt[row_label %chin% row_order$row_label]
  plot_dt[, hover := paste0(
    "Term: ", row_label,
    "<br>PC: ", PC,
    "<br>n in/out: ", n_in, " / ", n_out,
    "<br>effect z: ", signif(effect_z, 3),
    "<br>delta: ", signif(pc_delta, 3),
    "<br>q: ", signif(greater_p_adj, 3)
  )]
  z_wide <- dcast(plot_dt, row_label ~ PC, value.var = "effect_z")
  text_wide <- dcast(plot_dt, row_label ~ PC, value.var = "hover")
  z_wide <- z_wide[row_order, on = "row_label"]
  text_wide <- text_wide[row_order, on = "row_label"]
  pcs <- setdiff(names(z_wide), "row_label")
  z <- as.matrix(z_wide[, ..pcs])
  text <- as.matrix(text_wide[, ..pcs])
  max_abs <- suppressWarnings(quantile(abs(z), 0.98, na.rm = TRUE))
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
  p <- plot_ly(
    x = pcs,
    y = z_wide$row_label,
    z = z,
    text = text,
    hoverinfo = "text",
    type = "heatmap",
    zmin = -max_abs,
    zmax = max_abs,
    colorscale = list(c(0, "#2166AC"), c(0.5, "#F7F7F7"), c(1, "#B2182B")),
    colorbar = list(title = "effect z")
  ) %>%
    layout(
      title = "Metadata enrichments across dominant-state PCs",
      xaxis = list(title = "", automargin = TRUE),
      yaxis = list(title = "", automargin = TRUE),
      margin = list(l = 300, b = 80, t = 70, r = 30)
    )
  save_widget(p, pc_metadata_heatmap_output, "PC metadata heatmap")
}

make_state_correlation_heatmap <- function(cor_long, state_cols) {
  wide <- dcast(cor_long, state_1 ~ state_2, value.var = "pearson_r")
  wide <- wide[match(state_cols, state_1)]
  z <- as.matrix(wide[, ..state_cols])
  p <- plot_ly(
    x = state_cols,
    y = wide$state_1,
    z = z,
    type = "heatmap",
    zmin = -1,
    zmax = 1,
    colorscale = list(c(0, "#2166AC"), c(0.5, "#F7F7F7"), c(1, "#B2182B")),
    colorbar = list(title = "r")
  ) %>%
    layout(
      title = "Dominant-state score correlations",
      xaxis = list(title = "", tickangle = -35, automargin = TRUE),
      yaxis = list(title = "", automargin = TRUE),
      margin = list(l = 230, b = 150, t = 70, r = 30)
    )
  save_widget(p, state_correlation_heatmap_output, "State correlation heatmap")
}

make_pca_scatter <- function(clusters, pca_variance, metadata_fields) {
  if (!all(c("PC1", "PC2") %chin% names(clusters))) return(invisible(NULL))
  hover_fields <- intersect(
    c("Run", "sample", "dominant_state", "latent_cluster",
      priority_metadata_fields, metadata_fields),
    names(clusters)
  )
  hover_fields <- unique(hover_fields)
  hover_text <- apply(clusters[, ..hover_fields], 1, function(row) {
    paste(paste0(names(row), ": ", row), collapse = "<br>")
  })
  x_title <- paste0(
    "PC1 (", sprintf("%.1f", 100 * pca_variance[PC == "PC1", fraction_variance]), "%)"
  )
  y_title <- paste0(
    "PC2 (", sprintf("%.1f", 100 * pca_variance[PC == "PC2", fraction_variance]), "%)"
  )
  p <- plot_ly(
    clusters,
    x = ~PC1,
    y = ~PC2,
    color = ~latent_cluster,
    text = hover_text,
    hoverinfo = "text",
    type = "scatter",
    mode = "markers",
    marker = list(size = 6, opacity = 0.75)
  ) %>%
    layout(
      title = "Dominant-state PCA sample map",
      xaxis = list(title = x_title, zeroline = TRUE),
      yaxis = list(title = y_title, zeroline = TRUE),
      margin = list(l = 70, b = 70, t = 70, r = 20)
    )
  save_widget(p, pca_scatter_output, "Dominant-state PCA")
}

if (!file.exists(state_file)) stop("Missing ", state_file)
if (!file.exists(metadata_file)) stop("Missing ", metadata_file)

state_dt <- fread(state_file)
metadata <- fread(metadata_file)
metadata_fields <- read_selected_metadata_fields(metadata)
numeric_fields <- read_selected_numeric_fields(metadata)

state_dt[, input_order := .I]
keep_metadata <- intersect(c("Run", metadata_fields, numeric_fields), names(metadata))
dt <- merge(state_dt, metadata[, ..keep_metadata], by = "Run", all.x = TRUE, sort = FALSE)
setorder(dt, input_order)
dt[, input_order := NULL]

for (field in metadata_fields) {
  dt[, (field) := clean_level(get(field))]
}

state_cols <- setdiff(
  names(state_dt),
  c("Run", "sample", "dominant_state", "dominant_state_score",
    "Baseline (control)")
)
state_cols <- state_cols[state_cols %chin% names(dt)]
if (length(state_cols) < 2) stop("Need at least two state score columns.")

x <- as.matrix(dt[, ..state_cols])
storage.mode(x) <- "numeric"
col_medians <- apply(x, 2, median, na.rm = TRUE)
for (j in seq_len(ncol(x))) {
  x[!is.finite(x[, j]), j] <- col_medians[[j]]
}
x_scaled <- scale(x)
x_scaled[!is.finite(x_scaled)] <- 0

pca <- prcomp(x_scaled, center = FALSE, scale. = FALSE)
n_pc <- min(n_pcs_to_save, ncol(pca$x))
pc_cols <- paste0("PC", seq_len(n_pc))

pca_variance <- data.table(
  PC = paste0("PC", seq_along(pca$sdev)),
  variance = pca$sdev^2,
  fraction_variance = (pca$sdev^2) / sum(pca$sdev^2),
  cumulative_fraction_variance = cumsum((pca$sdev^2) / sum(pca$sdev^2))
)
fwrite(pca_variance, pca_variance_output)

pca_sample_scores <- as.data.table(pca$x[, seq_len(n_pc), drop = FALSE])
setnames(pca_sample_scores, pc_cols)
pca_samples <- cbind(
  data.table(
    Run = dt$Run,
    sample = if ("sample" %chin% names(dt)) dt$sample else NA_character_,
    dominant_state = dt$dominant_state
  ),
  pca_sample_scores,
  dt[, ..metadata_fields]
)
fwrite(pca_samples, pca_samples_output)

pca_loading_scores <- as.data.table(pca$rotation[, seq_len(n_pc), drop = FALSE])
setnames(pca_loading_scores, pc_cols)
pca_loadings <- cbind(data.table(dominant_state_score_name = rownames(pca$rotation)),
                      pca_loading_scores)
fwrite(pca_loadings, pca_loadings_output)

pca_top_loadings <- rbindlist(lapply(pc_cols, function(pc) {
  vals <- pca_loadings[, .(dominant_state_score_name, loading = get(pc))]
  rbind(
    vals[order(-loading)][seq_len(min(.N, top_loadings_per_pc))][
      ,
      `:=`(PC = pc, direction = "positive", rank = seq_len(.N))
    ],
    vals[order(loading)][seq_len(min(.N, top_loadings_per_pc))][
      ,
      `:=`(PC = pc, direction = "negative", rank = seq_len(.N))
    ],
    fill = TRUE
  )
}), fill = TRUE)
setcolorder(pca_top_loadings, c("PC", "direction", "rank",
                                "dominant_state_score_name", "loading"))
fwrite(pca_top_loadings, pca_top_loadings_output)

pca_extreme_samples <- rbindlist(lapply(pc_cols, function(pc) {
  vals <- pca_samples[, c("Run", "sample", "dominant_state", pc, metadata_fields),
                      with = FALSE]
  setnames(vals, pc, "PC_score")
  rbind(
    vals[order(-PC_score)][seq_len(min(.N, top_extreme_samples_per_pc))][
      ,
      `:=`(PC = pc, side = "positive", rank = seq_len(.N))
    ],
    vals[order(PC_score)][seq_len(min(.N, top_extreme_samples_per_pc))][
      ,
      `:=`(PC = pc, side = "negative", rank = seq_len(.N))
    ],
    fill = TRUE
  )
}), fill = TRUE)
setcolorder(pca_extreme_samples, c("PC", "side", "rank", "PC_score",
                                   setdiff(names(pca_extreme_samples),
                                           c("PC", "side", "rank", "PC_score"))))
fwrite(pca_extreme_samples, pca_extreme_samples_output)

pc_metadata_dt <- copy(pca_samples)
pc_metadata <- score_pc_metadata(pc_metadata_dt, metadata_fields, pc_cols)
if (nrow(pc_metadata) > 0) {
  pc_metadata[, `:=`(
    greater_p_adj = p.adjust(greater_p, "BH"),
    less_p_adj = p.adjust(less_p, "BH"),
    abs_effect_z = abs(effect_z)
  )]
  setorder(pc_metadata, greater_p_adj, -abs_effect_z, PC, field, level)
}
fwrite(pc_metadata, pc_metadata_output)
make_pc_metadata_heatmap(pc_metadata)

pc_numeric <- rbindlist(lapply(numeric_fields, function(field) {
  xnum <- as.numeric(dt[[field]])
  rbindlist(lapply(pc_cols, function(pc) {
    y <- as.numeric(pca_samples[[pc]])
    spearman <- safe_cor_test(xnum, y, "spearman")
    pearson <- safe_cor_test(xnum, y, "pearson")
    data.table(
      field = field,
      PC = pc,
      n = spearman$n,
      spearman_rho = spearman$estimate,
      spearman_p = spearman$p,
      pearson_r = pearson$estimate,
      pearson_p = pearson$p
    )
  }), fill = TRUE)
}), fill = TRUE)
if (nrow(pc_numeric) > 0) {
  pc_numeric[, `:=`(
    spearman_p_adj = p.adjust(spearman_p, "BH"),
    pearson_p_adj = p.adjust(pearson_p, "BH"),
    abs_spearman_rho = abs(spearman_rho)
  )]
  setorder(pc_numeric, spearman_p_adj, -abs_spearman_rho, field, PC)
}
fwrite(pc_numeric, pc_numeric_output)

state_cor_mat <- cor(x, use = "pairwise.complete.obs", method = "pearson")
state_correlations <- as.data.table(as.table(state_cor_mat))
setnames(state_correlations, c("state_1", "state_2", "pearson_r"))
fwrite(state_correlations, state_correlations_output)
make_state_correlation_heatmap(state_correlations, state_cols)

set.seed(1)
cluster_pc_count <- min(6L, ncol(pca$x))
first_90 <- which(pca_variance$cumulative_fraction_variance >= 0.90)[1]
if (is.finite(first_90)) cluster_pc_count <- min(cluster_pc_count, first_90)
k <- min(n_clusters, max(2L, floor(nrow(dt) / 20)))
km <- kmeans(pca$x[, seq_len(cluster_pc_count), drop = FALSE],
             centers = k, nstart = 80, iter.max = 150)

cluster_scores <- as.data.table(pca$x[, seq_len(n_pc), drop = FALSE])
setnames(cluster_scores, pc_cols)
clusters <- cbind(data.table(
  Run = dt$Run,
  sample = if ("sample" %chin% names(dt)) dt$sample else NA_character_,
  latent_cluster = paste0("LC", km$cluster),
  dominant_state = dt$dominant_state
), cluster_scores, dt[, ..metadata_fields])
fwrite(clusters, cluster_samples_output)
make_pca_scatter(clusters, pca_variance, metadata_fields)

state_summary <- rbindlist(lapply(state_cols, function(state) {
  data.table(
    latent_cluster = clusters$latent_cluster,
    dominant_state_score_name = state,
    score = x[, state]
  )[
    ,
    .(
      n_samples = .N,
      mean_score = mean(score, na.rm = TRUE),
      median_score = median(score, na.rm = TRUE)
    ),
    by = .(latent_cluster, dominant_state_score_name)
  ]
}), fill = TRUE)
state_summary[, cluster_rank := frank(-mean_score, ties.method = "dense"),
              by = latent_cluster]
setorder(state_summary, latent_cluster, cluster_rank,
         -mean_score, dominant_state_score_name)
fwrite(state_summary, cluster_state_summary_output)

cluster_metadata <- score_cluster_metadata(clusters, metadata_fields)
if (nrow(cluster_metadata) > 0) {
  cluster_metadata[, fisher_greater_p_adj := p.adjust(fisher_greater_p, "BH")]
  setorder(cluster_metadata, fisher_greater_p_adj, -odds_ratio,
           latent_cluster, field, level)
}
fwrite(cluster_metadata, cluster_metadata_output)

top_states <- state_summary[
  cluster_rank <= 3,
  .(
    top_states = paste(
      paste0(dominant_state_score_name, "=", sprintf("%.2f", mean_score)),
      collapse = "; "
    )
  ),
  by = latent_cluster
]
top_metadata <- cluster_metadata[
  is.finite(fisher_greater_p_adj) &
    fisher_greater_p_adj < 0.1 &
    n_level_in_cluster >= min_cluster_level_n
][
  ,
  head(.SD, 6L),
  by = latent_cluster
][
  ,
  .(
    top_metadata = paste(
      paste0(field, "=", level, " OR=", signif(odds_ratio, 3),
             " q=", signif(fisher_greater_p_adj, 3)),
      collapse = "; "
    )
  ),
  by = latent_cluster
]
cluster_sizes <- clusters[, .(n_samples = .N), by = latent_cluster]
compact_report <- merge(cluster_sizes, top_states, by = "latent_cluster",
                        all.x = TRUE, sort = FALSE)
compact_report <- merge(compact_report, top_metadata, by = "latent_cluster",
                        all.x = TRUE, sort = FALSE)
setorder(compact_report, -n_samples, latent_cluster)
fwrite(compact_report, compact_report_output)

message("Selected categorical metadata fields: ",
        paste(metadata_fields, collapse = ", "))
message("Selected numeric metadata fields: ",
        paste(numeric_fields, collapse = ", "))
message("Saved: ", pca_samples_output)
message("Saved: ", pca_loadings_output)
message("Saved: ", pca_top_loadings_output)
message("Saved: ", pca_variance_output)
message("Saved: ", pca_extreme_samples_output)
message("Saved: ", pc_metadata_output)
message("Saved: ", pc_numeric_output)
message("Saved: ", state_correlations_output)
message("Saved: ", cluster_samples_output)
message("Saved: ", cluster_state_summary_output)
message("Saved: ", cluster_metadata_output)
message("Saved: ", compact_report_output)

print(head(pca_variance, n_pc))
print(pca_top_loadings[PC %chin% pc_cols[seq_len(min(4L, length(pc_cols)))]])
print(compact_report)

invisible(list(
  pca_samples = pca_samples,
  pca_loadings = pca_loadings,
  pca_variance = pca_variance,
  pc_metadata = pc_metadata,
  pc_numeric = pc_numeric,
  state_correlations = state_correlations,
  clusters = clusters,
  state_summary = state_summary,
  cluster_metadata = cluster_metadata,
  compact_report = compact_report
))
