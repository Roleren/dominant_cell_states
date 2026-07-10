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
global_enrichment_file <- result_file("dominant_state_clean_cds_enrichment_all_main.csv")

contrast_output <- result_file("dominant_state_control_aware_condition_contrasts.csv")
condition_summary_output <- result_file("dominant_state_control_aware_condition_summary.csv")
gene_condition_summary_output <- result_file("dominant_state_control_aware_gene_condition_summary.csv")
gene_summary_output <- result_file("dominant_state_control_aware_gene_summary.csv")
review_output <- result_file("dominant_state_control_aware_enrichment_review.csv")
heatmap_output <- result_file("dominant_state_control_aware_condition_summary_heatmap.html")

metadata_fields <- c("CELL_LINE", "TISSUE", "GENE", "CONDITION")
stratum_fields <- c("CELL_LINE", "TISSUE", "GENE")
min_direct_n <- 2L

clean_level <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "NA", "N/A", "na", "n/a", "NULL", "null", "None", "none")] <- NA_character_
  x
}

normalise_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

is_control_condition <- function(x) {
  y <- normalise_label(x)
  y %chin% c(
    "wt", "wildtype", "wild_type", "control", "ctrl", "mock", "vehicle",
    "dmso", "untreated", "unperturbed", "baseline", "none", "empty_vector",
    "ev", "scramble", "scrambled", "scr", "shctrl", "sh_control",
    "non_targeting", "nontargeting", "nt", "ntc"
  ) |
    grepl("^(wt|ctrl|control|mock|vehicle|dmso)$", y) |
    grepl("(^|_)(wild_type|wildtype|untreated|non_targeting|nontargeting)($|_)", y)
}

classify_condition <- function(x) {
  y <- normalise_label(x)
  fifelse(
    grepl("(^|_)(ko|knockout|knock_out|kd|knockdown|knock_down|sirna|shrna|crispr|null)($|_)", y),
    "loss_of_function",
    fifelse(
      grepl("(^|_)(oe|overexpression|over_expression|rescue|reconstitution|addback)($|_)", y),
      "gain_of_function",
      fifelse(
        grepl("(infect|virus|viral|ifn|interferon|lps|tnf|stress|starv|hypox|treat|drug|stim)", y),
        "treatment_or_stimulus",
        "other_non_control"
      )
    )
  )
}

condition_family_helper <- c(
  file.path("scripts", "dominant_condition_families.R"),
  "dominant_condition_families.R",
  file.path("dominant_cell_states", "scripts",
            "dominant_condition_families.R")
)
condition_family_helper <-
  condition_family_helper[file.exists(condition_family_helper)][1]
if (!is.na(condition_family_helper)) source(condition_family_helper)
if (exists("dominant_is_control_condition", mode = "function")) {
  is_control_condition <- dominant_is_control_condition
}
if (exists("dominant_classify_perturbation_class", mode = "function")) {
  classify_condition <- function(x) dominant_classify_perturbation_class(x)
}
classify_design_family <- function(condition = NA_character_,
                                   cell_line = NA_character_,
                                   tissue = NA_character_,
                                   gene = NA_character_) {
  if (exists("dominant_classify_design_family", mode = "function")) {
    dominant_classify_design_family(
      condition = condition,
      cell_line = cell_line,
      tissue = tissue,
      gene = gene
    )
  } else {
    classify_condition(condition)
  }
}

safe_wilcox <- function(case_score, control_score, alternative = "greater") {
  case_score <- case_score[is.finite(case_score)]
  control_score <- control_score[is.finite(control_score)]
  if (length(case_score) < min_direct_n || length(control_score) < min_direct_n) {
    return(NA_real_)
  }
  if (length(unique(c(case_score, control_score))) < 2) return(NA_real_)
  tryCatch(
    wilcox.test(case_score, control_score, alternative = alternative,
                exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

stouffer_greater <- function(p, weights) {
  ok <- is.finite(p) & p > 0 & p <= 1 & is.finite(weights) & weights > 0
  if (!any(ok)) return(list(z = NA_real_, p = NA_real_))
  p <- pmin(pmax(p[ok], 1e-300), 1 - 1e-16)
  weights <- sqrt(weights[ok])
  z <- sum(qnorm(1 - p) * weights) / sqrt(sum(weights^2))
  list(z = z, p = pnorm(z, lower.tail = FALSE))
}

summarise_direct <- function(contrasts, by_cols) {
  out <- contrasts[
    ,
    {
      combined <- stouffer_greater(wilcox_greater_p, n_case + n_control)
      .(
        n_strata = .N,
        n_case_samples = sum(n_case),
        n_control_samples = sum(n_control),
        n_case_conditions = uniqueN(case_condition),
        median_score_delta = median(score_delta, na.rm = TRUE),
        weighted_score_delta = weighted.mean(
          score_delta,
          w = pmax(n_case + n_control, 1),
          na.rm = TRUE
        ),
        median_effect_z = median(effect_z, na.rm = TRUE),
        fraction_positive_delta = mean(score_delta > 0, na.rm = TRUE),
        stouffer_z = combined$z,
        stouffer_greater_p = combined$p,
        min_wilcox_greater_p = suppressWarnings(min(wilcox_greater_p, na.rm = TRUE))
      )
    },
    by = by_cols
  ]
  out[!is.finite(min_wilcox_greater_p), min_wilcox_greater_p := NA_real_]
  out[, stouffer_greater_p_adj := p.adjust(stouffer_greater_p, "BH")]
  setorder(out, stouffer_greater_p_adj, -weighted_score_delta)
  out[]
}

save_condition_heatmap <- function(summary_dt, state_cols) {
  plot_dt <- summary_dt[
    is.finite(stouffer_greater_p_adj) &
      stouffer_greater_p_adj < 0.2 &
      n_case_samples >= min_direct_n &
      n_control_samples >= min_direct_n
  ]
  if (nrow(plot_dt) == 0) return(invisible(NULL))

  plot_dt[, row_label := paste0(case_condition, " (", perturbation_class, ")")]
  row_order <- plot_dt[
    ,
    .(
      best_p = min(stouffer_greater_p_adj, na.rm = TRUE),
      best_delta = max(weighted_score_delta, na.rm = TRUE),
      n = max(n_case_samples, na.rm = TRUE)
    ),
    by = row_label
  ][order(best_p, -best_delta, -n)]
  row_order <- row_order[seq_len(min(.N, 80L))]

  plot_dt <- plot_dt[row_label %chin% row_order$row_label]
  plot_dt[, hover := paste0(
    "Condition: ", case_condition,
    "<br>Class: ", perturbation_class,
    "<br>State: ", dominant_state,
    "<br>Strata: ", n_strata,
    "<br>Case/control samples: ", n_case_samples, " / ", n_control_samples,
    "<br>Weighted score delta: ", signif(weighted_score_delta, 3),
    "<br>Stouffer z: ", signif(stouffer_z, 3),
    "<br>BH q: ", signif(stouffer_greater_p_adj, 3),
    "<br>Positive strata fraction: ", signif(fraction_positive_delta, 3)
  )]

  z_wide <- dcast(plot_dt, row_label ~ dominant_state,
                  value.var = "weighted_score_delta")
  text_wide <- dcast(plot_dt, row_label ~ dominant_state, value.var = "hover")
  z_wide <- z_wide[row_order, on = "row_label"]
  text_wide <- text_wide[row_order, on = "row_label"]

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
    type = "heatmap",
    colorscale = list(
      list(0, "#2166ac"),
      list(0.5, "#ffffff"),
      list(1, "#8b0000")
    ),
    zmid = 0,
    zmin = -max_abs,
    zmax = max_abs,
    text = text_matrix,
    hoverinfo = "text",
    colorbar = list(title = "direct delta")
  ) %>%
    layout(
      title = paste0(
        "Control-aware condition effects: direct case vs WT/control within ",
        paste(stratum_fields, collapse = " + ")
      ),
      xaxis = list(title = "", tickangle = -35),
      yaxis = list(title = "", automargin = TRUE),
      margin = list(l = 220, b = 120)
    )
  dominant_save_widget(
    p,
    output_file = heatmap_output,
    title = "Control-aware condition effects"
  )
}

if (!file.exists(state_file)) stop("Missing ", state_file)
if (!file.exists(metadata_file)) stop("Missing ", metadata_file)

dt <- fread(state_file)
m <- fread(metadata_file)
m <- m[, intersect(c("Run", metadata_fields), names(m)), with = FALSE]
dt <- merge(dt, m, by = "Run", all.x = TRUE, sort = FALSE)

state_cols <- setdiff(names(dt), c("Run", "sample", "dominant_state",
                                   "dominant_state_score", metadata_fields))
state_cols <- state_cols[state_cols != "Baseline (control)"]

for (field in metadata_fields) {
  if (!field %chin% names(dt)) dt[, (field) := NA_character_]
  dt[, (paste0(field, "_clean")) := clean_level(get(field))]
}
for (field in stratum_fields) {
  dt[is.na(get(paste0(field, "_clean"))), (paste0(field, "_clean")) := "MISSING"]
}
dt[, CONDITION_clean := clean_level(CONDITION)]
dt[, is_control_condition := is_control_condition(CONDITION_clean)]
dt[, design_family := classify_design_family(
  condition = CONDITION_clean,
  cell_line = CELL_LINE_clean,
  tissue = TISSUE_clean,
  gene = GENE_clean
)]
dt[, perturbation_class := classify_condition(CONDITION_clean)]
if (exists("dominant_classify_perturbation_class", mode = "function")) {
  dt[, perturbation_class :=
       dominant_classify_perturbation_class(CONDITION_clean, design_family)]
}
dt[, stratum_id := do.call(paste, c(.SD, sep = " | ")),
   .SDcols = paste0(stratum_fields, "_clean")]

contrasts <- rbindlist(lapply(unique(dt$stratum_id), function(stratum) {
  sdt <- dt[stratum_id == stratum & !is.na(CONDITION_clean)]
  controls <- sdt[is_control_condition == TRUE]
  cases <- sdt[is_control_condition != TRUE]
  if (nrow(controls) == 0 || nrow(cases) == 0) return(NULL)

  rbindlist(lapply(sort(unique(cases$CONDITION_clean)), function(case_condition) {
    case_dt <- cases[CONDITION_clean == case_condition]
    if (nrow(case_dt) == 0) return(NULL)

    rbindlist(lapply(state_cols, function(state) {
      case_score <- as.numeric(case_dt[[state]])
      control_score <- as.numeric(controls[[state]])
      pooled <- c(case_score, control_score)
      pooled_sd <- sd(pooled, na.rm = TRUE)
      score_delta <- mean(case_score, na.rm = TRUE) -
        mean(control_score, na.rm = TRUE)

      data.table(
        stratum_id = stratum,
        CELL_LINE = case_dt$CELL_LINE_clean[[1]],
        TISSUE = case_dt$TISSUE_clean[[1]],
        GENE = case_dt$GENE_clean[[1]],
        case_condition = case_condition,
        perturbation_class = case_dt$perturbation_class[[1]],
        design_family = case_dt$design_family[[1]],
        control_conditions = paste(sort(unique(controls$CONDITION_clean)),
                                   collapse = ";"),
        dominant_state = state,
        n_case = uniqueN(case_dt$Run),
        n_control = uniqueN(controls$Run),
        case_mean_score = mean(case_score, na.rm = TRUE),
        control_mean_score = mean(control_score, na.rm = TRUE),
        score_delta = score_delta,
        effect_z = if (is.finite(pooled_sd) && pooled_sd > 0) {
          score_delta / pooled_sd
        } else NA_real_,
        case_dominant_call_fraction =
          mean(case_dt$dominant_state == state, na.rm = TRUE),
        control_dominant_call_fraction =
          mean(controls$dominant_state == state, na.rm = TRUE),
        wilcox_greater_p = safe_wilcox(case_score, control_score, "greater"),
        wilcox_less_p = safe_wilcox(case_score, control_score, "less")
      )
    }), fill = TRUE)
  }), fill = TRUE)
}), fill = TRUE)

if (is.null(contrasts)) contrasts <- data.table()
if (nrow(contrasts) > 0) {
  contrasts[, `:=`(
    wilcox_greater_p_adj = p.adjust(wilcox_greater_p, "BH"),
    wilcox_less_p_adj = p.adjust(wilcox_less_p, "BH")
  )]
  setorder(contrasts, wilcox_greater_p_adj, -effect_z)
}
fwrite(contrasts, contrast_output)

condition_summary <- if (nrow(contrasts) > 0) {
  summarise_direct(
    contrasts,
    c("case_condition", "perturbation_class", "dominant_state")
  )
} else data.table()
fwrite(condition_summary, condition_summary_output)

gene_condition_summary <- if (nrow(contrasts) > 0) {
  summarise_direct(
    contrasts,
    c("GENE", "case_condition", "perturbation_class", "dominant_state")
  )
} else data.table()
fwrite(gene_condition_summary, gene_condition_summary_output)

gene_summary <- if (nrow(contrasts) > 0) {
  summarise_direct(contrasts, c("GENE", "dominant_state"))
} else data.table()
fwrite(gene_summary, gene_summary_output)

review <- data.table()
if (file.exists(global_enrichment_file) && nrow(condition_summary) > 0) {
  enrichment <- fread(global_enrichment_file)
  condition_review <- merge(
    enrichment[field == "CONDITION"],
    condition_summary,
    by.x = c("level", "dominant_state"),
    by.y = c("case_condition", "dominant_state"),
    all.x = TRUE,
    sort = FALSE
  )
  condition_review[, review_scope := "CONDITION"]

  gene_review <- merge(
    enrichment[field == "GENE"],
    gene_summary,
    by.x = c("level", "dominant_state"),
    by.y = c("GENE", "dominant_state"),
    all.x = TRUE,
    sort = FALSE
  )
  gene_review[, review_scope := "GENE"]

  review <- rbind(condition_review, gene_review, fill = TRUE)
  review[, direct_control_available := !is.na(n_strata)]
  review[, control_aware_status := fifelse(
    !direct_control_available,
    "no_matched_control_available",
    fifelse(
      is.finite(stouffer_greater_p_adj) &
        stouffer_greater_p_adj < 0.05 &
        weighted_score_delta > 0 &
        fraction_positive_delta >= 0.5,
      "supported_by_matched_control",
      fifelse(
        is.finite(weighted_score_delta) &
          (weighted_score_delta <= 0 | fraction_positive_delta < 0.5),
        "possible_global_false_positive",
        "matched_control_weak_or_underpowered"
      )
    )
  )]
  setorder(review, review_scope, direct_control_available,
           stouffer_greater_p_adj, adjusted_wilcox_greater_p_adj)
}
fwrite(review, review_output)

if (nrow(condition_summary) > 0) {
  save_condition_heatmap(condition_summary, state_cols)
}

message("Saved: ", contrast_output)
message("Saved: ", condition_summary_output)
message("Saved: ", gene_condition_summary_output)
message("Saved: ", gene_summary_output)
message("Saved: ", review_output)
if (file.exists(heatmap_output)) message("Saved: ", heatmap_output)

print(head(contrasts, 30))
print(head(condition_summary, 30))
print(head(review, 30))

invisible(list(
  contrasts = contrasts,
  condition_summary = condition_summary,
  gene_condition_summary = gene_condition_summary,
  gene_summary = gene_summary,
  review = review
))
