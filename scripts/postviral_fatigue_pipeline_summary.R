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

analysis_dir <- if (file.exists("human_dominant_cell_states_clean_cds.csv")) {
  "."
} else if (file.exists("dominant_cell_states/human_dominant_cell_states_clean_cds.csv")) {
  "dominant_cell_states"
} else {
  "."
}

fatigue_state <- "Post-viral fatigue / ribosome stress"
metadata_file <- "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv"
output_dir <- file.path(analysis_dir, "postviral_fatigue_outputs")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

state_file <- file.path(analysis_dir, "human_dominant_cell_states_clean_cds.csv")
if (!file.exists(state_file)) stop("Missing state file: ", state_file)

message("Post-viral fatigue pipeline summary:")
message("  1. Summarize sample-level fatigue scores")
message("  2. Collect fatigue metadata terms from enrichment/model outputs")
message("  3. Summarize RDG feature shifts, including internal ORF diagnostics")
message("  4. Save focused CSVs and inspection figures")

safe_fread <- function(path) {
  if (!file.exists(path)) return(data.table())
  fread(path)
}

clean_level <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "NA", "N/A", "na", "n/a", "NULL", "null", "None", "none")] <- NA_character_
  x
}

z_score <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

first_existing_col <- function(dt, cols, default = NA_character_) {
  col <- cols[cols %chin% names(dt)][1]
  if (is.na(col)) return(rep(default, nrow(dt)))
  dt[[col]]
}

classify_design_family <- function(condition = NA_character_,
                                   fraction = NA_character_,
                                   inhibitor = NA_character_,
                                   cell_line = NA_character_,
                                   tissue = NA_character_,
                                   gene = NA_character_,
                                   cancer_type = NA_character_) {
  text <- paste(condition, fraction, inhibitor, cell_line, tissue, gene,
                cancer_type, sep = " ")
  text <- tolower(text)
  cancer_clean <- clean_level(cancer_type)
  cancer_present <- !is.na(cancer_clean) & nzchar(cancer_clean)
  fifelse(
    grepl("ifn|interferon|ifng|cytokine|tnf|poly\\s*\\(?i:?c\\)?", text),
    "acute_interferon_cytokine",
    fifelse(
      grepl("infect|virus|viral|sars|cov|ebv|hsv|hiv|hcv|iav|influenza|adenovirus|reovirus", text),
      "viral_infection",
      fifelse(
        grepl("arsen|tunicamycin|thapsigargin|dtt|er stress|\\bhs\\b|heat|starv|stress|harringtonine|lactimidomycin", text),
        "isr_or_translation_stress",
        fifelse(
          grepl("senesc|irradiat|radiation|doxo|etoposide|dna damage|p53", text),
          "senescence_damage",
          fifelse(
            grepl("lcl|lymphoblast|macrophage|monocyte|cd14|pbmc|blood|b cell|t cell|immune", text),
            "immune_cell_or_lcl_context",
            fifelse(
              grepl("ko|kd|knock|sirna|shrna|crispr|oe|overexpress|mutant|deplet", text),
              "genetic_perturbation",
              fifelse(
                grepl("inhib|torin|rapa|bortez|drug|treated|vehicle|dmso", text),
                "drug_or_inhibitor",
                fifelse(
                  cancer_present,
                  "cancer_tissue_context",
                  "other_or_uncurated"
                )
              )
            )
          )
        )
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
if (exists("dominant_classify_design_family", mode = "function")) {
  classify_design_family <- dominant_classify_design_family
}

confidence01 <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(!is.finite(x))) return(x)
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0, length(x)))
  pmin(pmax((x - rng[1]) / diff(rng), 0), 1)
}

counterfactual_confidence <- function(n_design_rows, n_studies, n_supported_rows) {
  row_support <- 1 - exp(-pmax(0, as.numeric(n_design_rows)) / 10)
  study_support <- 1 - exp(-pmax(0, as.numeric(n_studies)) / 4)
  support_bonus <- 0.35 + 0.65 * (1 - exp(-pmax(0, as.numeric(n_supported_rows)) / 3))
  pmin(pmax(row_support * study_support * support_bonus, 0), 1)
}

safe_weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  weighted.mean(x[ok], w[ok])
}

random_effects_log2fc_interval <- function(effect, se, study, design) {
  empty <- list(
    rows = 0L,
    designs = 0L,
    studies = 0L,
    log2fc = NA_real_,
    log2fc_lower = NA_real_,
    log2fc_upper = NA_real_,
    pct_change = NA_real_,
    pct_change_lower = NA_real_,
    pct_change_upper = NA_real_,
    tau2 = NA_real_,
    i2 = NA_real_,
    q = NA_real_,
    effect_direction_fraction = NA_real_,
    ci_excludes_zero = FALSE
  )
  d <- data.table(
    effect = as.numeric(effect),
    se = as.numeric(se),
    study = clean_level(study),
    design = clean_level(design)
  )
  d <- d[is.finite(effect) & is.finite(se) & se > 0]
  if (nrow(d) == 0) return(empty)
  d[is.na(study), study := "missing_study"]
  d[is.na(design), design := paste0("row_", seq_len(.N))]
  d[, inverse_variance := 1 / (se^2)]

  # One clean-CDS row is expected per gene/design. Collapse defensively if a
  # future input repeats the same design signature.
  by_design <- d[
    ,
    .(
      effect = weighted.mean(effect, inverse_variance),
      se = sqrt(1 / sum(inverse_variance)),
      study = study[1],
      source_rows = .N
    ),
    by = design
  ]
  by_design[, inverse_variance := 1 / (se^2)]
  fixed_mean <- weighted.mean(
    by_design$effect,
    by_design$inverse_variance
  )
  q <- sum(
    by_design$inverse_variance * (by_design$effect - fixed_mean)^2,
    na.rm = TRUE
  )
  k <- nrow(by_design)
  df_q <- max(k - 1L, 0L)
  c_term <- sum(by_design$inverse_variance) -
    sum(by_design$inverse_variance^2) / sum(by_design$inverse_variance)
  tau2 <- if (k >= 2L && is.finite(c_term) && c_term > 0) {
    max(0, (q - df_q) / c_term)
  } else {
    0
  }
  by_design[, random_weight := 1 / (se^2 + tau2)]
  random_mean <- weighted.mean(by_design$effect, by_design$random_weight)
  random_se <- sqrt(1 / sum(by_design$random_weight))
  critical <- if (k >= 2L) {
    qt(0.975, df = max(1L, k - 1L))
  } else {
    qnorm(0.975)
  }
  lower <- random_mean - critical * random_se
  upper <- random_mean + critical * random_se
  same_direction <- if (is.finite(random_mean) && random_mean != 0) {
    mean(sign(by_design$effect[by_design$effect != 0]) == sign(random_mean))
  } else {
    NA_real_
  }
  list(
    rows = sum(by_design$source_rows),
    designs = k,
    studies = uniqueN(by_design$study),
    log2fc = random_mean,
    log2fc_lower = lower,
    log2fc_upper = upper,
    pct_change = 100 * (2^random_mean - 1),
    pct_change_lower = 100 * (2^lower - 1),
    pct_change_upper = 100 * (2^upper - 1),
    tau2 = tau2,
    i2 = fifelse(is.finite(q) & q > 0 & df_q > 0,
                 pmax(0, (q - df_q) / q), 0),
    q = q,
    effect_direction_fraction = same_direction,
    ci_excludes_zero = lower > 0 | upper < 0
  )
}

prefix_list <- function(x, prefix) {
  setNames(x, paste0(prefix, names(x)))
}

collapse_examples <- function(x, n = 6L) {
  x <- unique(clean_level(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  paste(head(x, n), collapse = ", ")
}

strict_clean_cds_support_labels <- c(
  "supported",
  "matched_clean_cds_relative_shift_supported",
  "matched_clean_cds_fpkm_shift_supported"
)

file_safe <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  substr(x, 1L, 160L)
}

`%||%` <- function(x, y) {
  if (length(x) == 0 || is.na(x)) y else x
}

short_state <- function(x) {
  x <- gsub(" / ", " ", x, fixed = TRUE)
  x <- gsub(" \\(.*\\)$", "", x)
  x <- gsub("mTOR Translation capacity.*", "mTOR/TOP", x)
  x <- gsub("Proliferation Cell cycle", "Proliferation", x)
  x <- gsub("Hypoxia HIF1A program", "Hypoxia", x)
  x <- gsub("Interferon antiviral", "Interferon", x)
  x <- gsub("EMT Mesenchymal shift", "EMT", x)
  x <- gsub("OXPHOS mitochondrial", "OXPHOS", x)
  x <- gsub("Post-viral fatigue ribosome stress", "Post-viral fatigue", x)
  x
}

wrap_label <- function(x, width = 58L) {
  vapply(as.character(x), function(one) {
    paste(strwrap(one, width = width), collapse = "\n")
  }, character(1))
}

save_plot_pair <- function(plot, name, width = 9, height = 6) {
  png_file <- file.path(figure_dir, paste0(name, ".png"))
  pdf_file <- file.path(figure_dir, paste0(name, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, dpi = 220, bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, bg = "white")
  c(png = png_file, pdf = pdf_file)
}

save_plotly <- function(plot, name) {
  html_file <- file.path(figure_dir, paste0(name, ".html"))
  dominant_save_ggplotly(plot, html_file, title = name)
  html_file
}

state_dt <- fread(state_file)
if (!fatigue_state %chin% names(state_dt)) {
  stop("State file does not contain the fatigue score column: ", fatigue_state)
}

sample_scores_file <- file.path(output_dir, "postviral_fatigue_sample_scores.csv")
metadata_fields <- character()
if (file.exists(metadata_file)) {
  metadata <- fread(metadata_file)
  metadata_fields <- setdiff(names(metadata), "Run")
  state_dt <- merge(state_dt, metadata, by = "Run", all.x = TRUE, sort = FALSE)
  if ("sample.x" %chin% names(state_dt) && !"sample" %chin% names(state_dt)) {
    setnames(state_dt, "sample.x", "sample")
  }
  if ("sample.y" %chin% names(state_dt)) {
    setnames(state_dt, "sample.y", "metadata_sample")
  }
} else if (file.exists(sample_scores_file)) {
  old_scores <- fread(sample_scores_file)
  fallback_meta <- setdiff(
    names(old_scores),
    c("Run", "sample", "dominant_state", "dominant_state_score",
      "postviral_fatigue_score", "fatigue_rank", "fatigue_quantile")
  )
  fallback_meta <- fallback_meta[fallback_meta %chin% names(old_scores)]
  if (length(fallback_meta) > 0) {
    message("Metadata file unavailable; reusing metadata columns from previous fatigue sample-score output.")
    state_dt <- merge(
      state_dt,
      unique(old_scores[, c("Run", fallback_meta), with = FALSE], by = "Run"),
      by = "Run",
      all.x = TRUE,
      sort = FALSE
    )
  }
}

score_cols <- setdiff(
  names(fread(state_file, nrows = 0)),
  c("Run", "sample", "dominant_state", "dominant_state_score")
)
score_cols <- score_cols[score_cols != "Baseline (control)"]
score_cols <- score_cols[score_cols %chin% names(state_dt)]

fatigue_scores <- state_dt[, .(
  Run,
  sample = if ("sample" %chin% names(state_dt)) sample else NA_character_,
  dominant_state,
  dominant_state_score,
  postviral_fatigue_score = get(fatigue_state)
)]
keep_meta <- intersect(c("CELL_LINE", "TISSUE", "CONDITION", "GENE", "Cancer_type",
                         "Cell_model", "Cell_type", "Organ_system", "Sex", "study",
                         "AUTHOR", "FRACTION", "INHIBITOR", "LibrarySelection",
                         "LibraryLayout", "LibraryStrategy"),
                       names(state_dt))
if (length(keep_meta) > 0) {
  fatigue_scores <- cbind(
    fatigue_scores,
    state_dt[, ..keep_meta]
  )
}
setorder(fatigue_scores, -postviral_fatigue_score)
fatigue_scores[, fatigue_rank := seq_len(.N)]
fatigue_scores[, fatigue_quantile := frank(postviral_fatigue_score, ties.method = "average") / .N]
fwrite(fatigue_scores, file.path(output_dir, "postviral_fatigue_sample_scores.csv"))
fwrite(fatigue_scores[seq_len(min(.N, 200L))],
       file.path(output_dir, "postviral_fatigue_top_samples.csv"))

submodule_states <- list(
  interferon_antiviral = "Interferon / antiviral",
  isr_ribosome_stress = "ISR / ER stress (ATF4 axis)",
  senescence_sasp = "Senescence / SASP",
  hypoxia_inflammation = "Hypoxia / HIF1A program",
  suppressed_growth_translation = "mTOR / Translation capacity (TOP program)",
  suppressed_oxphos = "OXPHOS (mitochondrial)"
)
submodule_scores <- copy(state_dt[, .(Run)])
submodule_scores[, postviral_fatigue_z := z_score(state_dt[[fatigue_state]])]
for (nm in names(submodule_states)) {
  state_name <- submodule_states[[nm]]
  if (state_name %chin% names(state_dt)) {
    direction <- if (grepl("^suppressed_", nm)) -1 else 1
    submodule_scores[, (paste0(nm, "_z")) :=
                       postviral_fatigue_z + direction * z_score(state_dt[[state_name]])]
  }
}
submodule_cols <- setdiff(names(submodule_scores), c("Run", "postviral_fatigue_z"))
if (length(submodule_cols) > 0) {
  submodule_scores[, fatigue_submodule := submodule_cols[
    max.col(as.matrix(.SD), ties.method = "first")
  ], .SDcols = submodule_cols]
  submodule_scores[, fatigue_submodule_score := do.call(pmax, c(.SD, na.rm = TRUE)),
                   .SDcols = submodule_cols]
  submodule_scores[, fatigue_submodule := sub("_z$", "", fatigue_submodule)]
  submodule_scores <- merge(
    submodule_scores,
    fatigue_scores,
    by = "Run",
    all.x = TRUE,
    sort = FALSE
  )
  fwrite(submodule_scores, file.path(output_dir, "postviral_fatigue_submodule_sample_scores.csv"))

  submodule_summary <- submodule_scores[
    ,
    .(
      n = .N,
      mean_fatigue_score = mean(postviral_fatigue_score, na.rm = TRUE),
      median_fatigue_score = median(postviral_fatigue_score, na.rm = TRUE),
      top_quantile_samples = sum(fatigue_quantile >= 0.95, na.rm = TRUE)
    ),
    by = fatigue_submodule
  ][order(-top_quantile_samples, -mean_fatigue_score)]
  fwrite(submodule_summary, file.path(output_dir, "postviral_fatigue_submodule_summary.csv"))
}

winner_counts <- state_dt[, .N, by = dominant_state][order(-N)]
score_summary <- rbindlist(lapply(score_cols, function(state) {
  x <- as.numeric(state_dt[[state]])
  data.table(
    dominant_state_score = state,
    n = sum(is.finite(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    q05 = quantile(x, 0.05, na.rm = TRUE),
    q95 = quantile(x, 0.95, na.rm = TRUE),
    winner_count = winner_counts[dominant_state == state, N][1] %||% 0L
  )
}), fill = TRUE)
fwrite(score_summary, file.path(output_dir, "postviral_fatigue_state_score_summary.csv"))

score_matrix <- as.matrix(state_dt[, ..score_cols])
storage.mode(score_matrix) <- "numeric"
state_cor <- cor(score_matrix, use = "pairwise.complete.obs")
cor_dt <- as.data.table(as.table(state_cor))
setnames(cor_dt, c("state_x", "state_y", "correlation"))
cor_dt[, `:=`(state_x_short = short_state(state_x), state_y_short = short_state(state_y))]
fwrite(cor_dt, file.path(output_dir, "postviral_fatigue_state_correlations.csv"))

p_density <- ggplot(state_dt, aes(x = get(fatigue_state))) +
  geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "#6b7280",
                 color = "white", linewidth = 0.15) +
  geom_density(color = "#8b0000", linewidth = 0.8) +
  geom_vline(xintercept = quantile(state_dt[[fatigue_state]], c(0.9, 0.95), na.rm = TRUE),
             linetype = c("dashed", "solid"), color = c("#1f77b4", "#8b0000"),
             linewidth = 0.45) +
  labs(
    title = "Post-viral fatigue / ribosome stress score distribution",
    subtitle = "Dashed and solid lines mark the 90th and 95th percentiles",
    x = "signed dominant-state score",
    y = "density"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
save_plot_pair(p_density, "postviral_fatigue_score_distribution", width = 8.2, height = 5.2)

cor_plot <- copy(cor_dt)
cor_plot[, state_x_short := factor(state_x_short, levels = short_state(score_cols))]
cor_plot[, state_y_short := factor(state_y_short, levels = rev(short_state(score_cols)))]
p_cor <- ggplot(cor_plot, aes(state_x_short, state_y_short, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 2.8) +
  scale_fill_gradient2(low = "#2166ac", mid = "#ffffff", high = "#8b0000",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  labs(
    title = "Dominant-state score correlations",
    subtitle = "Correlation structure after adding the post-viral fatigue score",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )
save_plot_pair(p_cor, "postviral_fatigue_state_correlation_heatmap", width = 8.8, height = 7.2)
save_plotly(p_cor, "postviral_fatigue_state_correlation_heatmap")

collect_term_table <- function(path, source_name, effect_cols, q_cols) {
  dt <- safe_fread(path)
  if (nrow(dt) == 0 || !"dominant_state" %chin% names(dt)) return(data.table())
  dt <- dt[dominant_state == fatigue_state]
  if (nrow(dt) == 0) return(data.table())
  effect_col <- effect_cols[effect_cols %chin% names(dt)][1]
  q_col <- q_cols[q_cols %chin% names(dt)][1]
  if (is.na(effect_col)) effect_col <- NA_character_
  if (is.na(q_col)) q_col <- NA_character_
  dt[, source := source_name]
  dt[, effect_value := if (!is.na(effect_col)) as.numeric(get(effect_col)) else NA_real_]
  dt[, q_value := if (!is.na(q_col)) as.numeric(get(q_col)) else NA_real_]
  if (!"field" %chin% names(dt)) dt[, field := NA_character_]
  if (!"level" %chin% names(dt)) dt[, level := NA_character_]
  if (!"n_in" %chin% names(dt)) dt[, n_in := NA_integer_]
  dt[, term_label := paste0(source, " | ", field, "=", level)]
  dt[, .(source, field, level, term_label, dominant_state, n_in, effect_value, q_value)]
}

term_tables <- rbindlist(list(
  collect_term_table(
    file.path(analysis_dir, "dominant_state_clean_cds_enrichment_all_main.csv"),
    "adjusted_enrichment",
    c("adjusted_effect_z", "effect_z"),
    c("adjusted_wilcox_greater_p_adj", "raw_wilcox_greater_p_adj",
      "fisher_greater_p_adj", "q")
  ),
  collect_term_table(
    file.path(analysis_dir, "dominant_state_broad_metadata_rdg_terms.csv"),
    "broad_metadata",
    c("adjusted_effect_z", "effect_z"),
    c("adjusted_wilcox_greater_p_adj", "wilcox_greater_p_adj", "q")
  ),
  collect_term_table(
    file.path(analysis_dir, "dominant_state_glmnet_metadata_selected_term_scores.csv"),
    "glmnet",
    c("coefficient", "coefficient_abs", "effect_z", "adjusted_effect_z"),
    c("greater_p_adj", "q", "p_adj")
  ),
  collect_term_table(
    file.path(analysis_dir, "dominant_state_control_aware_enrichment_review.csv"),
    "matched_control_review",
    c("direct_control_delta_mean_score", "direct_control_delta_call_fraction",
      "adjusted_effect_z", "effect_z"),
    c("stouffer_greater_p_adj", "adjusted_wilcox_greater_p_adj",
      "raw_wilcox_greater_p_adj", "fisher_greater_p_adj", "q")
  )
), fill = TRUE)

if (nrow(term_tables) > 0) {
  term_tables[, abs_effect := abs(effect_value)]
  term_tables[, design_family := classify_design_family(
    condition = fifelse(field == "CONDITION", level, NA_character_),
    fraction = fifelse(field == "FRACTION", level, NA_character_),
    inhibitor = fifelse(field == "INHIBITOR", level, NA_character_),
    cell_line = fifelse(field %chin% c("CELL_LINE", "Cell_model", "Cell_type"), level, NA_character_),
    tissue = fifelse(field %chin% c("TISSUE", "Organ_system"), level, NA_character_),
    gene = fifelse(field == "GENE", level, NA_character_),
    cancer_type = fifelse(field == "Cancer_type", level, NA_character_)
  )]
  term_tables[, q_sort := fifelse(is.finite(q_value), q_value, Inf)]
  setorder(term_tables, q_sort, -abs_effect)
  term_tables[, q_sort := NULL]
}
fwrite(term_tables, file.path(output_dir, "postviral_fatigue_metadata_terms.csv"))

if (nrow(term_tables) > 0 && "design_family" %chin% names(term_tables)) {
  term_family_summary <- term_tables[
    ,
    .(
      n_terms = .N,
      best_q = suppressWarnings(min(q_value, na.rm = TRUE)),
      max_abs_effect = suppressWarnings(max(abs_effect, na.rm = TRUE)),
      top_term = term_label[which.max(abs_effect)][1]
    ),
    by = design_family
  ][order(best_q, -max_abs_effect)]
  term_family_summary[!is.finite(best_q), best_q := NA_real_]
  term_family_summary[!is.finite(max_abs_effect), max_abs_effect := NA_real_]
} else {
  term_family_summary <- data.table()
}
fwrite(term_family_summary, file.path(output_dir, "postviral_fatigue_design_family_terms.csv"))

plot_terms <- term_tables[is.finite(effect_value)]
if (nrow(plot_terms) > 0) {
  plot_terms <- plot_terms[order(fifelse(is.finite(q_value), q_value, Inf), -abs(effect_value))]
  plot_terms <- plot_terms[seq_len(min(.N, 45L))]
  plot_terms[, source_short := fifelse(
    source == "matched_control_review", "matched-control",
    fifelse(source == "adjusted_enrichment", "adjusted",
            fifelse(source == "broad_metadata", "broad", source))
  )]
  plot_terms[, term_label_plot := wrap_label(paste0(source_short, " | ", field, "=", level), 62L)]
  plot_terms[, term_label_plot := factor(term_label_plot, levels = rev(unique(term_label_plot)))]
  plot_terms[, q_plot := fifelse(is.finite(q_value), pmax(q_value, 1e-300), 1)]
  p_terms <- ggplot(plot_terms, aes(effect_value, term_label_plot, color = source, size = -log10(q_plot))) +
    geom_vline(xintercept = 0, color = "#d1d5db", linewidth = 0.35) +
    geom_point(alpha = 0.88) +
    scale_size_continuous(name = "-log10(q)", range = c(1.7, 6)) +
    labs(
      title = "Post-viral fatigue metadata terms",
      subtitle = "Positive x means higher fatigue score or positive model coefficient",
      x = "effect / model coefficient",
      y = NULL,
      color = "source"
    ) +
    theme_minimal(base_size = 10.5) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
  save_plot_pair(p_terms, "postviral_fatigue_metadata_terms", width = 13.2, height = 10.2)
  save_plotly(p_terms, "postviral_fatigue_metadata_terms")
}

branch_contrasts <- safe_fread(file.path(analysis_dir, "dominant_rdg_outputs", "rdg_branch_feature_contrasts.csv"))
feature_annotation <- safe_fread(file.path(analysis_dir, "dominant_rdg_outputs", "rdg_feature_annotation.csv"))
feature_confidence <- safe_fread(file.path(
  analysis_dir,
  "dominant_next_model",
  "dominant_next_model_feature_confidence.csv"
))
fatigue_contrasts <- branch_contrasts[dominant_state == fatigue_state]
if (nrow(fatigue_contrasts) > 0) {
  fatigue_contrasts[, abs_delta_relative_translation := abs(delta_relative_translation)]
  setorder(fatigue_contrasts, -abs_delta_relative_translation)
}
fwrite(fatigue_contrasts, file.path(output_dir, "postviral_fatigue_rdg_feature_contrasts.csv"))

iorf_candidates <- fatigue_contrasts[
  feature_type == "internal_orf" &
    is.finite(delta_relative_translation)
]
if (nrow(iorf_candidates) > 0 && nrow(feature_annotation) > 0) {
  iorf_candidates <- merge(
    iorf_candidates,
    feature_annotation[, .(
      gene_symbol, tx_id, feature_id, tx_start, tx_end, start_codon,
      start_codon_class, kozak_strength, frame_relative_to_cds,
      peptide_length_aa, coding_potential_proxy, translon_sources_merged
    )],
    by = c("gene_symbol", "tx_id", "feature_id"),
    all.x = TRUE,
    sort = FALSE
  )
}
if (nrow(iorf_candidates) > 0 && nrow(feature_confidence) > 0) {
  for (column in c(
    "internal_frame_evidence_class",
    "internal_start_shape_class",
    "weighted_internal_iorf_frame_fraction",
    "weighted_internal_canonical_frame_fraction",
    "median_internal_iorf_vs_canonical_frame_delta",
    "n_internal_frame_evaluable_runs",
    "internal_iorf_frame_advantage_run_fraction",
    "median_internal_start_vs_body_ratio"
  )) {
    if (!column %chin% names(feature_confidence)) {
      feature_confidence[, (column) := NA]
    }
  }
  iorf_candidates <- merge(
    iorf_candidates,
    feature_confidence[, .(
      gene_symbol, tx_id, feature_id, internal_orf_tc_only,
      cds_signal_bump_supported, feature_vs_cds_mean_ratio_median,
      feature_vs_cds_mean_ratio_q90, runs_with_cds_signal_bump,
      cds_signal_bump_fraction, internal_compendium_bump_class,
      internal_frame_evidence_class, internal_start_shape_class,
      weighted_internal_iorf_frame_fraction,
      weighted_internal_canonical_frame_fraction,
      median_internal_iorf_vs_canonical_frame_delta,
      n_internal_frame_evaluable_runs,
      internal_iorf_frame_advantage_run_fraction,
      median_internal_start_vs_body_ratio,
      n_internal_condition_bump_review_rows,
      n_internal_condition_bump_supported_rows,
      n_internal_condition_bump_review_studies,
      n_internal_condition_bump_supported_studies,
      best_internal_condition_bump_class,
      best_internal_condition_bump_context, feature_measurement_class
    )],
    by = c("gene_symbol", "tx_id", "feature_id"),
    all.x = TRUE,
    sort = FALSE
  )
}
if (nrow(iorf_candidates) > 0) {
  iorf_candidates[, iorf_evidence_score := 0]
  iorf_candidates[start_codon_class %chin% c("ATG", "near_cognate"),
                  iorf_evidence_score := iorf_evidence_score + 1]
  iorf_candidates[kozak_strength %chin% c("strong", "moderate"),
                  iorf_evidence_score := iorf_evidence_score + 1]
  iorf_candidates[is.finite(peptide_length_aa) & peptide_length_aa >= 10,
                  iorf_evidence_score := iorf_evidence_score + 1]
  iorf_candidates[is.finite(frame_relative_to_cds) & frame_relative_to_cds != 0,
                  iorf_evidence_score := iorf_evidence_score + 1]
  iorf_candidates[is.finite(abs_delta_relative_translation) &
                    abs_delta_relative_translation >= 0.10,
                  iorf_evidence_score := iorf_evidence_score + 1]
  iorf_candidates[is.finite(log2_mean_fpkm_like_in_vs_out) &
                    abs(log2_mean_fpkm_like_in_vs_out) >= 0.50,
                  iorf_evidence_score := iorf_evidence_score + 1]
  translon_source_col <- c(
    "translon_sources_merged",
    "translon_sources_merged.x",
    "translon_sources_merged.y"
  )
  translon_source_col <- translon_source_col[
    translon_source_col %chin% names(iorf_candidates)
  ][1]
  if (is.na(translon_source_col)) {
    iorf_candidates[, translon_sources_for_evidence := ""]
  } else {
    iorf_candidates[, translon_sources_for_evidence := fifelse(
      is.na(get(translon_source_col)),
      "",
      as.character(get(translon_source_col))
    )]
  }
  iorf_candidates[grepl("M", translon_sources_for_evidence),
                  iorf_evidence_score := iorf_evidence_score + 1]
  iorf_candidates[cds_signal_bump_supported == TRUE,
                  iorf_evidence_score := iorf_evidence_score + 1]
  iorf_candidates[n_internal_condition_bump_supported_rows > 0,
                  iorf_evidence_score := iorf_evidence_score + 1]
  iorf_candidates[internal_frame_evidence_class ==
                    "iorf_frame_advantage_review",
                  iorf_evidence_score := iorf_evidence_score + 1]
  iorf_candidates[, iorf_evidence_class := fifelse(
    is.finite(frame_relative_to_cds) & frame_relative_to_cds == 0,
    "in_frame_cds_shared_coverage_flag",
    fifelse(
      internal_frame_evidence_class == "canonical_frame_like",
      "canonical_frame_dominant_coverage_review",
      fifelse(
      internal_orf_tc_only == TRUE & cds_signal_bump_supported != TRUE,
      "tc_iORF_without_cds_signal_bump",
      fifelse(
        internal_orf_tc_only == TRUE &
          internal_compendium_bump_class == "sparse_upper_tail_cds_bump" &
          !(n_internal_condition_bump_supported_rows > 0),
        "tc_iORF_sparse_compendium_bump_review",
        fifelse(
      iorf_evidence_score >= 5,
      "high_priority_iORF_review",
      fifelse(iorf_evidence_score >= 3,
              "medium_priority_iORF_review",
              "low_evidence_coverage_routing_flag")
        )
      )
      )
    )
  )]
  setorder(iorf_candidates, -abs_delta_relative_translation)
}
fwrite(iorf_candidates, file.path(output_dir, "postviral_fatigue_iorf_candidates.csv"))
fwrite(iorf_candidates, file.path(output_dir, "postviral_fatigue_iorf_evidence_classifier.csv"))

feature_shift <- fatigue_contrasts[
  feature_type %chin% c("overlapping_uorf", "internal_orf", "leader_uorf",
                        "nte_extension", "ntt_truncation", "clean_CDS") &
    is.finite(delta_relative_translation),
  .(
    max_abs_delta_relative_translation = max(abs(delta_relative_translation), na.rm = TRUE),
    strongest_delta_relative_translation = delta_relative_translation[which.max(abs(delta_relative_translation))],
    strongest_log2_feature_change = log2_mean_fpkm_like_in_vs_out[which.max(abs(delta_relative_translation))],
    strongest_branch = branch_label[which.max(abs(delta_relative_translation))][1],
    n_branches = uniqueN(branch_id)
  ),
  by = .(gene_symbol, tx_id, feature_id, feature_type)
]
setorder(feature_shift, -max_abs_delta_relative_translation)
fwrite(feature_shift, file.path(output_dir, "postviral_fatigue_rdg_feature_shift_summary.csv"))

if (nrow(feature_shift) > 0) {
  top_genes <- feature_shift[
    ,
    .(gene_score = max(max_abs_delta_relative_translation, na.rm = TRUE)),
    by = gene_symbol
  ][order(-gene_score)][seq_len(min(.N, 35L)), gene_symbol]
  plot_shift <- feature_shift[gene_symbol %chin% top_genes]
  plot_shift[, gene_symbol := factor(gene_symbol, levels = rev(top_genes))]
  plot_shift[, feature_label := paste0(feature_id, " (", feature_type, ")")]
  p_shift <- ggplot(plot_shift, aes(feature_label, gene_symbol,
                                    fill = strongest_delta_relative_translation)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(low = "#2166ac", mid = "#ffffff", high = "#8b0000",
                         midpoint = 0, name = "delta\nrelative use") +
    labs(
      title = "Post-viral fatigue RDG feature shifts",
      subtitle = "IN minus OUT relative-use change; includes internal ORF diagnostics",
      x = "RDG feature",
      y = NULL
    ) +
    theme_minimal(base_size = 10.2) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
  save_plot_pair(p_shift, "postviral_fatigue_rdg_feature_shift_heatmap",
                 width = 12.5, height = 8.8)
  save_plotly(p_shift, "postviral_fatigue_rdg_feature_shift_heatmap")

  top_feature_shift <- copy(feature_shift[seq_len(min(.N, 45L))])
  top_feature_shift[, feature_label_plot := wrap_label(
    paste0(gene_symbol, " | ", feature_id, " (", feature_type, ")"),
    48L
  )]
  top_feature_shift[, feature_label_plot := factor(
    feature_label_plot,
    levels = rev(unique(feature_label_plot))
  )]
  p_top_shift <- ggplot(top_feature_shift, aes(y = feature_label_plot)) +
    geom_vline(xintercept = 0, color = "#cbd5e1", linewidth = 0.35) +
    geom_segment(aes(x = 0, xend = strongest_delta_relative_translation,
                     yend = feature_label_plot),
                 color = "#cbd5e1", linewidth = 0.5) +
    geom_point(aes(x = strongest_delta_relative_translation,
                   color = feature_type,
                   size = max_abs_delta_relative_translation),
               alpha = 0.9) +
    scale_size_continuous(name = "abs(delta)", range = c(1.8, 5.5)) +
    labs(
      title = "Strongest post-viral RDG feature shifts",
      subtitle = "IN minus OUT relative-use delta, ranked by absolute shift",
      x = "delta relative use",
      y = NULL,
      color = "feature type"
    ) +
    theme_minimal(base_size = 10.5) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
  save_plot_pair(p_top_shift, "postviral_fatigue_top_rdg_feature_shifts",
                 width = 11.8, height = 10.2)
  save_plotly(p_top_shift, "postviral_fatigue_top_rdg_feature_shifts")
}

cds_buffering <- safe_fread(file.path(
  analysis_dir,
  "dominant_next_model",
  "dominant_next_model_cds_buffering_effect_sizes_relevant_rdgs.csv"
))
fatigue_buffering <- cds_buffering[dominant_state == fatigue_state]
if (nrow(fatigue_buffering) > 0) {
  fatigue_buffering[, abs_cds_buffering_effect := abs(cds_fpkm_change_pct_stabilized)]
  setorder(fatigue_buffering, -buffering_confidence_score, -abs_cds_buffering_effect)
}
fwrite(fatigue_buffering, file.path(output_dir, "postviral_fatigue_cds_buffering_effects.csv"))

matched_cds <- safe_fread(file.path(
  analysis_dir,
  "dominant_next_model",
  "dominant_next_model_matched_control_cds_buffering.csv"
))
matched_viral_ifn_review <- data.table()
if (nrow(matched_cds) > 0) {
  matched_cds[, design_family := classify_design_family(
    condition = case_condition,
    cell_line = CELL_LINE,
    tissue = TISSUE,
    gene = GENE
  )]
  matched_cds[, clean_cds_strict_support :=
                matched_cds_support %chin% strict_clean_cds_support_labels]
  primary_postviral_families <- c(
    "viral_infection",
    "acute_interferon_cytokine"
  )
  adjacent_postviral_families <- c(
    "isr_or_translation_stress",
    "isr_er_translation_stress",
    "senescence_damage",
    "dna_damage_senescence",
    "immune_cell_or_lcl_context",
    "immune_or_infection_model_context",
    "hypoxia_oxidative_stress"
  )
  matched_cds[, context_priority := fifelse(
    design_family %chin% primary_postviral_families,
    "postviral_primary",
    fifelse(design_family %chin% adjacent_postviral_families,
            "postviral_adjacent",
            "background_context")
  )]
  matched_cds[, matched_prediction_weight := pmax(0, n_case_samples + n_control_samples) *
                (1 + pmax(0, matched_bootstrap_certainty_score, na.rm = TRUE)) *
                (1 + pmax(0, matched_design_specificity, na.rm = TRUE))]
  matched_cds[!is.finite(matched_prediction_weight) | matched_prediction_weight <= 0,
              matched_prediction_weight := pmax(1, n_case_samples + n_control_samples)]
  matched_cds[, prediction_effect := fifelse(
    is.finite(shrunk_matched_fpkm_log2_fc) & shrunk_matched_fpkm_log2_fc != 0,
    shrunk_matched_fpkm_log2_fc,
    matched_fpkm_log2_fc_pseudocount1
  )]

  design_family_summary <- matched_cds[
    feature_id == "clean_CDS",
    .(
      n_rows = .N,
      n_genes = uniqueN(gene_symbol),
      supported_rows = sum(clean_cds_strict_support == TRUE, na.rm = TRUE),
      median_abs_clean_cds_delta = median(abs(matched_relative_use_delta), na.rm = TRUE),
      median_abs_log2fc = median(abs(prediction_effect), na.rm = TRUE)
    ),
    by = .(design_family, context_priority)
  ][order(context_priority, -supported_rows, -median_abs_log2fc)]
  fwrite(design_family_summary, file.path(output_dir, "postviral_fatigue_matched_design_family_summary.csv"))

  counterfactual <- matched_cds[
    feature_id == "clean_CDS" &
      design_family %chin% c(primary_postviral_families,
                             adjacent_postviral_families) &
      is.finite(prediction_effect),
    .(
      predicted_clean_cds_log2fc = weighted.mean(prediction_effect, matched_prediction_weight, na.rm = TRUE),
      predicted_clean_cds_log2fc_lower = quantile(prediction_effect, 0.10, na.rm = TRUE),
      predicted_clean_cds_log2fc_upper = quantile(prediction_effect, 0.90, na.rm = TRUE),
      n_design_rows = .N,
      n_studies = uniqueN(study),
      n_designs = uniqueN(matched_design_signature),
      n_supported_rows = sum(clean_cds_strict_support == TRUE, na.rm = TRUE),
      nearest_designs = paste(head(unique(matched_design_signature), 5L), collapse = " || ")
    ),
    by = .(gene_symbol, tx_id, design_family)
  ]
  counterfactual[, predicted_clean_cds_pct_change :=
                   100 * (2^predicted_clean_cds_log2fc - 1)]
  counterfactual[, prediction_confidence := counterfactual_confidence(
    n_design_rows, n_studies, n_supported_rows
  )]
  counterfactual[, abs_predicted_clean_cds_pct_change :=
                   abs(predicted_clean_cds_pct_change)]
  setorder(counterfactual, -prediction_confidence,
           -abs_predicted_clean_cds_pct_change)
  fwrite(counterfactual, file.path(output_dir, "postviral_fatigue_counterfactual_predictions.csv"))

  if (nrow(counterfactual) > 0) {
    plot_counter <- counterfactual[
      ,
      head(.SD, 40L),
      by = design_family
    ]
    top_counter_genes <- plot_counter[
      ,
      .(score = max(abs(predicted_clean_cds_pct_change), na.rm = TRUE)),
      by = gene_symbol
    ][order(-score)][seq_len(min(.N, 35L)), gene_symbol]
    plot_counter <- plot_counter[gene_symbol %chin% top_counter_genes]
    plot_counter[, gene_symbol := factor(gene_symbol, levels = rev(top_counter_genes))]
    p_counter <- ggplot(plot_counter, aes(design_family, gene_symbol,
                                          fill = predicted_clean_cds_pct_change)) +
      geom_tile(color = "white", linewidth = 0.25) +
      scale_fill_gradient2(low = "#2166ac", mid = "#ffffff", high = "#8b0000",
                           midpoint = 0, name = "predicted\nCDS %") +
      labs(
        title = "Counterfactual clean-CDS predictions from matched designs",
        subtitle = "Prototype: gene x design-family estimates from observed matched-control RDG contrasts",
        x = "design family",
        y = NULL
      ) +
      theme_minimal(base_size = 10.2) +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold")
      )
    save_plot_pair(p_counter, "postviral_fatigue_counterfactual_predictions",
                   width = 10.8, height = 8.4)
    save_plotly(p_counter, "postviral_fatigue_counterfactual_predictions")
  }
} else {
  design_family_summary <- data.table()
  counterfactual <- data.table()
}

candidate_rankings <- safe_fread(file.path(
  analysis_dir,
  "dominant_next_model",
  "dominant_next_model_candidate_rankings_v4.csv"
))
if (nrow(matched_cds) > 0) {
  matched_viral_ifn_review <- matched_cds[
    feature_id == "clean_CDS" &
      context_priority %chin% c("postviral_primary", "postviral_adjacent") &
      is.finite(prediction_effect)
  ]
  if (nrow(matched_viral_ifn_review) > 0) {
    matched_viral_ifn_review[, matched_clean_cds_pct_change :=
                               100 * (2^prediction_effect - 1)]
    matched_viral_ifn_review[, matched_postviral_review_score :=
                               0.25 * pmin(abs(matched_clean_cds_pct_change) / 100, 1) +
                               0.20 * pmin(abs(shrunk_matched_relative_use_delta) / 0.25, 1) +
                               0.20 * pmin(pmax(matched_feature_confidence_score, 0), 1) +
                               0.20 * pmin(pmax(matched_bootstrap_certainty_score, 0), 1) +
                               0.30 * as.numeric(clean_cds_strict_support == TRUE) +
                               0.15 * as.numeric(context_priority == "postviral_primary") +
                               0.10 * pmin(pmax(matched_design_specificity, 0), 1)]
    matched_viral_ifn_review[
      !is.finite(matched_postviral_review_score),
      matched_postviral_review_score := 0
    ]
    if (nrow(candidate_rankings) > 0) {
      candidate_subset <- candidate_rankings[
        ,
        .SD[1],
        by = .(gene_symbol, tx_id)
      ][, .(
        gene_symbol, tx_id, candidate_rank_v4, novelty_class,
        known_uorf_regulation, manual_m_source_gene,
        has_overlapping_uorf, n_uorfs, n_overlapping_uorfs,
        primary_clean_cds_usage, next_iteration_candidate_score_v4
      )]
      matched_viral_ifn_review <- merge(
        matched_viral_ifn_review,
        candidate_subset,
        by = c("gene_symbol", "tx_id"),
        all.x = TRUE,
        sort = FALSE
      )
    }
    matched_viral_ifn_review[, review_context := paste(
      design_family,
      paste0("case=", case_condition),
      paste0("control=", control_conditions),
      sep = " | "
    )]
    matched_viral_ifn_review[, abs_matched_clean_cds_pct_change :=
                               abs(matched_clean_cds_pct_change)]
    setorder(
      matched_viral_ifn_review,
      -matched_postviral_review_score,
      -abs_matched_clean_cds_pct_change
    )
    matched_viral_ifn_review[, matched_postviral_review_rank := seq_len(.N)]
    front_cols <- c(
      "matched_postviral_review_rank", "gene_symbol", "tx_id",
      "design_family", "context_priority", "case_condition",
      "control_conditions", "CELL_LINE", "TISSUE", "GENE",
      "study", "AUTHOR", "n_case_samples", "n_control_samples",
      "matched_clean_cds_pct_change", "prediction_effect",
      "matched_relative_use_delta", "shrunk_matched_relative_use_delta",
      "matched_cds_support", "clean_cds_strict_support",
      "matched_feature_confidence_score", "matched_bootstrap_certainty_score",
      "matched_design_specificity", "matched_postviral_review_score",
      "candidate_rank_v4", "novelty_class", "primary_clean_cds_usage",
      "review_context", "matched_design_signature"
    )
    setcolorder(
      matched_viral_ifn_review,
      c(intersect(front_cols, names(matched_viral_ifn_review)),
        setdiff(names(matched_viral_ifn_review), front_cols))
    )
  }
}
fwrite(
  matched_viral_ifn_review,
  file.path(output_dir, "postviral_fatigue_matched_viral_ifn_review.csv")
)

matched_design_cards <- data.table()
matched_gene_design_evidence <- data.table()
matched_gene_design_loo <- data.table()
matched_gene_design_study_loo <- data.table()
matched_gene_design_summary <- data.table()
if (nrow(matched_viral_ifn_review) > 0) {
  matched_design_cards <- matched_viral_ifn_review[
    ,
    {
      ranked <- .SD[order(-matched_postviral_review_score,
                          -abs_matched_clean_cds_pct_change)]
      strict_ranked <- ranked[clean_cds_strict_support == TRUE]
      .(
        n_gene_rows = .N,
        n_genes = uniqueN(gene_symbol),
        n_strict_supported_rows = sum(clean_cds_strict_support == TRUE,
                                      na.rm = TRUE),
        n_strict_supported_genes =
          uniqueN(gene_symbol[clean_cds_strict_support == TRUE]),
        n_manual_m_source_genes =
          uniqueN(gene_symbol[manual_m_source_gene == TRUE]),
        n_known_uorf_genes =
          uniqueN(gene_symbol[known_uorf_regulation == TRUE]),
        median_clean_cds_pct_change =
          median(matched_clean_cds_pct_change, na.rm = TRUE),
        median_abs_clean_cds_pct_change =
          median(abs(matched_clean_cds_pct_change), na.rm = TRUE),
        fraction_positive_effects =
          mean(prediction_effect > 0, na.rm = TRUE),
        max_review_score = max(matched_postviral_review_score, na.rm = TRUE),
        top_gene = ranked$gene_symbol[1],
        top_gene_clean_cds_pct_change =
          ranked$matched_clean_cds_pct_change[1],
        top_review_genes = collapse_examples(ranked$gene_symbol, 8L),
        strict_supported_genes =
          collapse_examples(strict_ranked$gene_symbol, 8L)
      )
    },
    by = .(
      matched_design_signature, design_family, context_priority, study,
      AUTHOR, CELL_LINE, TISSUE, GENE, case_condition, control_conditions,
      n_case_samples, n_control_samples
    )
  ]
  matched_design_cards[, design_card_score :=
                         0.40 * pmin(n_strict_supported_genes / 6, 1) +
                         0.15 * pmin(n_genes / 20, 1) +
                         0.20 * pmin(max_review_score, 1) +
                         0.15 * as.numeric(context_priority == "postviral_primary") +
                         0.10 * pmin(median_abs_clean_cds_pct_change / 100, 1)]
  matched_design_cards[, design_card_class := fifelse(
    n_strict_supported_genes >= 2,
    "strict_multi_gene_design",
    fifelse(n_strict_supported_genes == 1,
            "strict_single_gene_design",
            "descriptive_design")
  )]
  setorder(matched_design_cards, -design_card_score,
           -n_strict_supported_genes, -n_genes)

  matched_gene_design_evidence <- matched_viral_ifn_review[
    ,
    {
      best_idx <- which.max(matched_postviral_review_score)
      effect_idx <- which.max(abs_matched_clean_cds_pct_change)
      random_effects_all <- random_effects_log2fc_interval(
        prediction_effect,
        bootstrap_fpkm_log2_fc_sd,
        study,
        matched_design_signature
      )
      random_effects_strict <- random_effects_log2fc_interval(
        prediction_effect[clean_cds_strict_support == TRUE],
        bootstrap_fpkm_log2_fc_sd[clean_cds_strict_support == TRUE],
        study[clean_cds_strict_support == TRUE],
        matched_design_signature[clean_cds_strict_support == TRUE]
      )
      weighted_log2fc <- safe_weighted_mean(
        prediction_effect,
        matched_prediction_weight
      )
      c(list(
        n_design_rows = .N,
        n_designs = uniqueN(matched_design_signature),
        n_studies = uniqueN(study),
        strict_supported_rows = sum(clean_cds_strict_support == TRUE,
                                    na.rm = TRUE),
        strict_supported_designs =
          uniqueN(matched_design_signature[clean_cds_strict_support == TRUE]),
        strict_supported_studies =
          uniqueN(study[clean_cds_strict_support == TRUE]),
        weighted_clean_cds_log2fc = weighted_log2fc,
        weighted_clean_cds_pct_change = 100 * (2^weighted_log2fc - 1),
        median_clean_cds_pct_change =
          median(matched_clean_cds_pct_change, na.rm = TRUE),
        median_abs_clean_cds_pct_change =
          median(abs_matched_clean_cds_pct_change, na.rm = TRUE),
        weighted_shrunk_relative_use_delta = safe_weighted_mean(
          shrunk_matched_relative_use_delta,
          matched_prediction_weight
        ),
        fraction_positive_effects =
          mean(prediction_effect > 0, na.rm = TRUE),
        best_review_score =
          max(matched_postviral_review_score, na.rm = TRUE),
        best_condition = case_condition[best_idx],
        best_control = control_conditions[best_idx],
        best_study = study[best_idx],
        strongest_effect_condition = case_condition[effect_idx],
        strongest_effect_study = study[effect_idx],
        strongest_abs_clean_cds_pct_change =
          abs_matched_clean_cds_pct_change[effect_idx],
        nearest_designs =
          collapse_examples(matched_design_signature, 4L)
      ),
      prefix_list(random_effects_all, "random_effects_all_"),
      prefix_list(random_effects_strict, "random_effects_strict_"))
    },
    by = .(gene_symbol, tx_id, design_family, context_priority)
  ]
  matched_gene_design_evidence[, gene_design_evidence_class := fifelse(
    strict_supported_studies >= 2 & strict_supported_designs >= 2,
    "replicated_strict_matched",
    fifelse(strict_supported_designs >= 1,
            "single_design_strict_matched",
            fifelse(n_studies >= 2,
                    "replicated_descriptive_matched",
                    "single_design_descriptive"))
  )]
  matched_gene_design_evidence[, matched_gene_design_score :=
                                 0.35 * pmin(strict_supported_designs / 2, 1) +
                                 0.30 * pmin(strict_supported_studies / 2, 1) +
                                 0.20 * pmin(n_studies / 4, 1) +
                                 0.15 * pmin(best_review_score, 1)]
  matched_gene_design_evidence[, abs_weighted_clean_cds_pct_change :=
                                 abs(weighted_clean_cds_pct_change)]
  setorder(matched_gene_design_evidence, -matched_gene_design_score,
           -abs_weighted_clean_cds_pct_change)

  matched_gene_design_loo <- matched_viral_ifn_review[
    ,
    {
      design_ids <- unique(matched_design_signature)
      full_log2fc <- safe_weighted_mean(
        prediction_effect,
        matched_prediction_weight
      )
      full_weight <- sum(
        matched_prediction_weight[
          is.finite(matched_prediction_weight) & matched_prediction_weight > 0
        ],
        na.rm = TRUE
      )
      if (length(design_ids) < 2L || !is.finite(full_log2fc)) {
        NULL
      } else {
        rbindlist(lapply(design_ids, function(excluded_design) {
          excluded <- .SD[matched_design_signature == excluded_design]
          remaining <- .SD[matched_design_signature != excluded_design]
          loo_log2fc <- safe_weighted_mean(
            remaining$prediction_effect,
            remaining$matched_prediction_weight
          )
          excluded_weight <- sum(
            excluded$matched_prediction_weight[
              is.finite(excluded$matched_prediction_weight) &
                excluded$matched_prediction_weight > 0
            ],
            na.rm = TRUE
          )
          data.table(
            excluded_design_signature = excluded_design,
            excluded_study = excluded$study[1],
            excluded_case_condition = excluded$case_condition[1],
            excluded_control_conditions = excluded$control_conditions[1],
            excluded_clean_cds_strict_support =
              any(excluded$clean_cds_strict_support == TRUE, na.rm = TRUE),
            excluded_design_weight_fraction =
              fifelse(full_weight > 0, excluded_weight / full_weight, NA_real_),
            full_weighted_clean_cds_log2fc = full_log2fc,
            full_weighted_clean_cds_pct_change =
              100 * (2^full_log2fc - 1),
            loo_weighted_clean_cds_log2fc = loo_log2fc,
            loo_weighted_clean_cds_pct_change =
              100 * (2^loo_log2fc - 1),
            loo_minus_full_log2fc = loo_log2fc - full_log2fc,
            loo_same_effect_direction = fifelse(
              is.finite(loo_log2fc) & full_log2fc != 0 & loo_log2fc != 0,
              sign(full_log2fc) == sign(loo_log2fc),
              NA
            ),
            remaining_designs = uniqueN(remaining$matched_design_signature),
            remaining_studies = uniqueN(remaining$study),
            remaining_strict_supported_designs =
              uniqueN(remaining$matched_design_signature[
                remaining$clean_cds_strict_support == TRUE
              ]),
            remaining_strict_supported_studies =
              uniqueN(remaining$study[
                remaining$clean_cds_strict_support == TRUE
              ])
          )
        }), fill = TRUE)
      }
    },
    by = .(gene_symbol, tx_id, design_family, context_priority)
  ]
  if (nrow(matched_gene_design_loo) > 0) {
    loo_summary <- matched_gene_design_loo[
      ,
      .(
        loo_excluded_designs = .N,
        loo_min_clean_cds_log2fc =
          min(loo_weighted_clean_cds_log2fc, na.rm = TRUE),
        loo_max_clean_cds_log2fc =
          max(loo_weighted_clean_cds_log2fc, na.rm = TRUE),
        loo_min_clean_cds_pct_change =
          min(loo_weighted_clean_cds_pct_change, na.rm = TRUE),
        loo_max_clean_cds_pct_change =
          max(loo_weighted_clean_cds_pct_change, na.rm = TRUE),
        loo_fraction_same_effect_direction =
          mean(loo_same_effect_direction == TRUE, na.rm = TRUE),
        loo_max_abs_log2fc_shift =
          max(abs(loo_minus_full_log2fc), na.rm = TRUE),
        loo_top_design_weight_fraction =
          max(excluded_design_weight_fraction, na.rm = TRUE),
        loo_min_remaining_strict_designs =
          min(remaining_strict_supported_designs, na.rm = TRUE),
        loo_min_remaining_strict_studies =
          min(remaining_strict_supported_studies, na.rm = TRUE),
        loo_most_influential_design =
          excluded_design_signature[which.max(abs(loo_minus_full_log2fc))][1],
        loo_most_influential_study =
          excluded_study[which.max(abs(loo_minus_full_log2fc))][1]
      ),
      by = .(gene_symbol, tx_id, design_family, context_priority)
    ]
    matched_gene_design_evidence <- merge(
      matched_gene_design_evidence,
      loo_summary,
      by = c("gene_symbol", "tx_id", "design_family", "context_priority"),
      all.x = TRUE,
      sort = FALSE
    )
  }
  matched_gene_design_study_loo <- matched_viral_ifn_review[
    ,
    {
      study_ids <- unique(study)
      full_log2fc <- safe_weighted_mean(
        prediction_effect,
        matched_prediction_weight
      )
      full_weight <- sum(
        matched_prediction_weight[
          is.finite(matched_prediction_weight) & matched_prediction_weight > 0
        ],
        na.rm = TRUE
      )
      if (length(study_ids) < 2L || !is.finite(full_log2fc)) {
        NULL
      } else {
        rbindlist(lapply(study_ids, function(excluded_study) {
          excluded <- .SD[study == excluded_study]
          remaining <- .SD[study != excluded_study]
          loo_log2fc <- safe_weighted_mean(
            remaining$prediction_effect,
            remaining$matched_prediction_weight
          )
          excluded_weight <- sum(
            excluded$matched_prediction_weight[
              is.finite(excluded$matched_prediction_weight) &
                excluded$matched_prediction_weight > 0
            ],
            na.rm = TRUE
          )
          data.table(
            excluded_study = excluded_study,
            excluded_designs = uniqueN(excluded$matched_design_signature),
            excluded_design_signatures =
              collapse_examples(excluded$matched_design_signature, 4L),
            excluded_case_conditions =
              collapse_examples(excluded$case_condition, 4L),
            excluded_clean_cds_strict_support =
              any(excluded$clean_cds_strict_support == TRUE, na.rm = TRUE),
            excluded_study_weight_fraction =
              fifelse(full_weight > 0, excluded_weight / full_weight, NA_real_),
            full_weighted_clean_cds_log2fc = full_log2fc,
            full_weighted_clean_cds_pct_change =
              100 * (2^full_log2fc - 1),
            loo_weighted_clean_cds_log2fc = loo_log2fc,
            loo_weighted_clean_cds_pct_change =
              100 * (2^loo_log2fc - 1),
            loo_minus_full_log2fc = loo_log2fc - full_log2fc,
            loo_same_effect_direction = fifelse(
              is.finite(loo_log2fc) & full_log2fc != 0 & loo_log2fc != 0,
              sign(full_log2fc) == sign(loo_log2fc),
              NA
            ),
            remaining_designs = uniqueN(remaining$matched_design_signature),
            remaining_studies = uniqueN(remaining$study),
            remaining_strict_supported_designs =
              uniqueN(remaining$matched_design_signature[
                remaining$clean_cds_strict_support == TRUE
              ]),
            remaining_strict_supported_studies =
              uniqueN(remaining$study[
                remaining$clean_cds_strict_support == TRUE
              ])
          )
        }), fill = TRUE)
      }
    },
    by = .(gene_symbol, tx_id, design_family, context_priority)
  ]
  if (nrow(matched_gene_design_study_loo) > 0) {
    study_loo_summary <- matched_gene_design_study_loo[
      ,
      .(
        study_loo_excluded_studies = .N,
        study_loo_min_clean_cds_log2fc =
          min(loo_weighted_clean_cds_log2fc, na.rm = TRUE),
        study_loo_max_clean_cds_log2fc =
          max(loo_weighted_clean_cds_log2fc, na.rm = TRUE),
        study_loo_min_clean_cds_pct_change =
          min(loo_weighted_clean_cds_pct_change, na.rm = TRUE),
        study_loo_max_clean_cds_pct_change =
          max(loo_weighted_clean_cds_pct_change, na.rm = TRUE),
        study_loo_fraction_same_effect_direction =
          mean(loo_same_effect_direction == TRUE, na.rm = TRUE),
        study_loo_max_abs_log2fc_shift =
          max(abs(loo_minus_full_log2fc), na.rm = TRUE),
        study_loo_top_study_weight_fraction =
          max(excluded_study_weight_fraction, na.rm = TRUE),
        study_loo_min_remaining_strict_designs =
          min(remaining_strict_supported_designs, na.rm = TRUE),
        study_loo_min_remaining_strict_studies =
          min(remaining_strict_supported_studies, na.rm = TRUE),
        study_loo_most_influential_study =
          excluded_study[which.max(abs(loo_minus_full_log2fc))][1]
      ),
      by = .(gene_symbol, tx_id, design_family, context_priority)
    ]
    matched_gene_design_evidence <- merge(
      matched_gene_design_evidence,
      study_loo_summary,
      by = c("gene_symbol", "tx_id", "design_family", "context_priority"),
      all.x = TRUE,
      sort = FALSE
    )
  }
  matched_gene_design_evidence[, loo_effect_direction_stable :=
                                 is.finite(loo_fraction_same_effect_direction) &
                                   loo_fraction_same_effect_direction >= 0.999 &
                                   abs(weighted_clean_cds_log2fc) >= 0.10]
  matched_gene_design_evidence[, study_loo_effect_direction_stable :=
                                 is.finite(study_loo_fraction_same_effect_direction) &
                                   study_loo_fraction_same_effect_direction >= 0.999 &
                                   abs(weighted_clean_cds_log2fc) >= 0.10]
  matched_gene_design_evidence[, strict_random_effects_interval_supported :=
                                 random_effects_strict_designs >= 2 &
                                   random_effects_strict_studies >= 2 &
                                   random_effects_strict_ci_excludes_zero == TRUE]
  matched_gene_design_evidence[, atlas_condition_evidence_gate := fifelse(
    strict_supported_designs >= 3 & strict_supported_studies >= 3 &
      loo_effect_direction_stable == TRUE &
      study_loo_effect_direction_stable == TRUE &
      strict_random_effects_interval_supported == TRUE &
      loo_min_remaining_strict_designs >= 2 &
      loo_min_remaining_strict_studies >= 2 &
      study_loo_min_remaining_strict_designs >= 2 &
      study_loo_min_remaining_strict_studies >= 2,
    "atlas_ready_strict_design_study_interval",
    fifelse(
      strict_supported_designs >= 2 & strict_supported_studies >= 2 &
        loo_effect_direction_stable == TRUE &
        study_loo_effect_direction_stable == TRUE,
      "replicated_strict_interval_review",
      fifelse(
        strict_supported_designs >= 1,
        "single_design_strict_review",
        fifelse(
          n_studies >= 2 & loo_effect_direction_stable == TRUE,
          "descriptive_loo_stable",
          "descriptive_review"
        )
      )
    )
  )]
  matched_gene_design_evidence[, atlas_condition_evidence_score :=
                                 matched_gene_design_score +
                                   0.14 * as.numeric(loo_effect_direction_stable == TRUE) +
                                   0.14 * as.numeric(study_loo_effect_direction_stable == TRUE) +
                                   0.16 * as.numeric(
                                     strict_random_effects_interval_supported == TRUE
                                   ) +
                                   0.10 * pmax(
                                     0,
                                     1 - fifelse(
                                       is.finite(loo_top_design_weight_fraction),
                                       loo_top_design_weight_fraction,
                                       1
                                     )
                                   ) +
                                   0.30 * as.numeric(
                                     atlas_condition_evidence_gate ==
                                       "atlas_ready_strict_design_study_interval"
                                   ) +
                                   0.12 * as.numeric(
                                     atlas_condition_evidence_gate ==
                                       "replicated_strict_interval_review"
                                   )]
  setorder(matched_gene_design_evidence, -atlas_condition_evidence_score,
           -abs_weighted_clean_cds_pct_change)

  matched_gene_design_summary <- matched_viral_ifn_review[
    ,
    {
      best_idx <- which.max(matched_postviral_review_score)
      weighted_log2fc <- safe_weighted_mean(
        prediction_effect,
        matched_prediction_weight
      )
      .(
        postviral_context_rows = .N,
        postviral_design_families = uniqueN(design_family),
        postviral_primary_design_families =
          uniqueN(design_family[context_priority == "postviral_primary"]),
        postviral_designs = uniqueN(matched_design_signature),
        postviral_studies = uniqueN(study),
        postviral_strict_rows = sum(clean_cds_strict_support == TRUE,
                                    na.rm = TRUE),
        postviral_strict_designs =
          uniqueN(matched_design_signature[clean_cds_strict_support == TRUE]),
        postviral_strict_studies =
          uniqueN(study[clean_cds_strict_support == TRUE]),
        postviral_weighted_clean_cds_log2fc = weighted_log2fc,
        postviral_weighted_clean_cds_pct_change =
          100 * (2^weighted_log2fc - 1),
        postviral_median_clean_cds_pct_change =
          median(matched_clean_cds_pct_change, na.rm = TRUE),
        postviral_fraction_positive_effects =
          mean(prediction_effect > 0, na.rm = TRUE),
        postviral_best_family = design_family[best_idx],
        postviral_best_condition = case_condition[best_idx],
        postviral_best_study = study[best_idx],
        postviral_best_review_score =
          matched_postviral_review_score[best_idx]
      )
    },
    by = .(gene_symbol, tx_id)
  ]
  matched_gene_design_summary[, postviral_design_replication_score :=
                                0.45 * pmin(postviral_strict_studies / 2, 1) +
                                0.25 * pmin(postviral_strict_designs / 2, 1) +
                                0.15 * pmin(postviral_studies / 4, 1) +
                                0.15 * pmin(postviral_best_review_score, 1)]
  matched_gene_design_summary[, postviral_design_evidence_class := fifelse(
    postviral_strict_studies >= 2 & postviral_strict_designs >= 2,
    "replicated_strict_matched",
    fifelse(postviral_strict_designs >= 1,
            "single_design_strict_matched",
            fifelse(postviral_studies >= 2,
                    "replicated_descriptive_matched",
                    "single_design_descriptive"))
  )]
  matched_gene_design_summary[, abs_postviral_weighted_clean_cds_pct_change :=
                                abs(postviral_weighted_clean_cds_pct_change)]
  if (nrow(matched_gene_design_evidence) > 0) {
    best_gene_design_gate <- matched_gene_design_evidence[
      ,
      .SD[which.max(atlas_condition_evidence_score)],
      by = .(gene_symbol, tx_id)
    ][
      ,
      .(
        gene_symbol, tx_id,
        postviral_best_gene_design_family = design_family,
        postviral_best_gene_design_context = context_priority,
        postviral_best_gene_design_gate = atlas_condition_evidence_gate,
        postviral_best_gene_design_score = atlas_condition_evidence_score,
        postviral_best_gene_design_loo_direction_stable =
          loo_effect_direction_stable,
        postviral_best_gene_design_loo_same_direction_fraction =
          loo_fraction_same_effect_direction,
        postviral_best_gene_design_loo_top_weight_fraction =
          loo_top_design_weight_fraction,
        postviral_best_gene_design_loo_min_strict_designs =
          loo_min_remaining_strict_designs,
        postviral_best_gene_design_study_loo_direction_stable =
          study_loo_effect_direction_stable,
        postviral_best_gene_design_study_loo_same_direction_fraction =
          study_loo_fraction_same_effect_direction,
        postviral_best_gene_design_study_loo_top_weight_fraction =
          study_loo_top_study_weight_fraction,
        postviral_best_gene_design_strict_interval_supported =
          strict_random_effects_interval_supported,
        postviral_best_gene_design_strict_random_effects_pct_change =
          random_effects_strict_pct_change,
        postviral_best_gene_design_strict_random_effects_pct_change_lower =
          random_effects_strict_pct_change_lower,
        postviral_best_gene_design_strict_random_effects_pct_change_upper =
          random_effects_strict_pct_change_upper
      )
    ]
    matched_gene_design_summary <- merge(
      matched_gene_design_summary,
      best_gene_design_gate,
      by = c("gene_symbol", "tx_id"),
      all.x = TRUE,
      sort = FALSE
    )
  }
  setorder(matched_gene_design_summary, -postviral_design_replication_score,
           -abs_postviral_weighted_clean_cds_pct_change)
}
fwrite(
  matched_design_cards,
  file.path(output_dir, "postviral_fatigue_matched_design_cards.csv")
)
fwrite(
  matched_gene_design_evidence,
  file.path(output_dir, "postviral_fatigue_matched_gene_design_evidence.csv")
)
fwrite(
  matched_gene_design_loo,
  file.path(output_dir, "postviral_fatigue_matched_gene_design_loo.csv")
)
fwrite(
  matched_gene_design_study_loo,
  file.path(output_dir, "postviral_fatigue_matched_gene_design_study_loo.csv")
)
fwrite(
  matched_gene_design_summary,
  file.path(output_dir, "postviral_fatigue_matched_gene_design_summary.csv")
)

if (nrow(matched_viral_ifn_review) > 0) {
  plot_review <- matched_viral_ifn_review[
    seq_len(min(nrow(matched_viral_ifn_review), 70L))
  ]
  top_review_genes <- plot_review[
    ,
    .(score = max(matched_postviral_review_score, na.rm = TRUE)),
    by = gene_symbol
  ][order(-score)][seq_len(min(.N, 35L)), gene_symbol]
  plot_review <- plot_review[gene_symbol %chin% top_review_genes]
  plot_review[, gene_symbol := factor(gene_symbol, levels = rev(top_review_genes))]
  plot_review[, design_family_plot := gsub("_", " ", design_family, fixed = TRUE)]
  plot_review[, review_context_plot := wrap_label(
    paste0(design_family_plot, " | ", case_condition),
    38L
  )]
  p_matched_review <- ggplot(
    plot_review,
    aes(review_context_plot, gene_symbol, fill = matched_clean_cds_pct_change)
  ) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_point(
      aes(size = matched_postviral_review_score,
          shape = clean_cds_strict_support),
      color = "#111827",
      alpha = 0.72
    ) +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "#ffffff",
      high = "#8b0000",
      midpoint = 0,
      limits = c(-100, 400),
      oob = scales::squish,
      name = "matched\nCDS %\n(capped)"
    ) +
    scale_size_continuous(name = "review score", range = c(1.2, 4.8)) +
    scale_shape_manual(
      name = "strict support",
      values = c("FALSE" = 1, "TRUE" = 16)
    ) +
    labs(
      title = "Matched viral/IFN clean-CDS review",
      subtitle = "Matched-control clean-CDS proxy changes in viral, IFN, ISR, senescence, and immune/LCL contexts",
      x = "design family and case condition",
      y = NULL
    ) +
    theme_minimal(base_size = 10.2) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 7.4),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
  save_plot_pair(
    p_matched_review,
    "postviral_fatigue_matched_viral_ifn_review",
    width = 12.4,
    height = 8.8
  )
  save_plotly(p_matched_review, "postviral_fatigue_matched_viral_ifn_review")
}

if (nrow(matched_design_cards) > 0) {
  plot_cards <- matched_design_cards[seq_len(min(.N, 24L))]
  plot_cards[, design_label := wrap_label(
    paste(
      study,
      paste0(CELL_LINE, "/", TISSUE),
      paste0(case_condition, " vs ", control_conditions),
      gsub("_", " ", design_family, fixed = TRUE),
      sep = " | "
    ),
    74L
  )]
  plot_cards[, design_label := factor(design_label,
                                      levels = rev(unique(design_label)))]
  p_design_cards <- ggplot(
    plot_cards,
    aes(n_strict_supported_genes, design_label)
  ) +
    geom_segment(
      aes(x = 0, xend = n_strict_supported_genes,
          yend = design_label),
      color = "#cbd5e1",
      linewidth = 0.65
    ) +
    geom_point(
      aes(size = n_genes, color = median_abs_clean_cds_pct_change,
          shape = design_card_class),
      alpha = 0.92
    ) +
    scale_color_gradient(
      low = "#9ccae5",
      high = "#8b0000",
      name = "median abs\nCDS %"
    ) +
    scale_size_continuous(name = "reviewed genes", range = c(2.1, 6.7)) +
    scale_shape_manual(
      name = "design card",
      values = c(
        "strict_multi_gene_design" = 16,
        "strict_single_gene_design" = 17,
        "descriptive_design" = 1
      )
    ) +
    labs(
      title = "Post-viral matched design cards",
      subtitle = "Each row keeps one study and control contrast visible; x is genes with strict clean-CDS support.",
      x = "strict-supported genes in design",
      y = NULL
    ) +
    theme_minimal(base_size = 10.1) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
  save_plot_pair(
    p_design_cards,
    "postviral_fatigue_matched_design_cards",
    width = 13.2,
    height = 10.0
  )
  save_plotly(p_design_cards, "postviral_fatigue_matched_design_cards")
}

if (nrow(matched_gene_design_evidence) > 0) {
  top_gene_designs <- matched_gene_design_summary[
    seq_len(min(.N, 32L)),
    gene_symbol
  ]
  plot_gene_design <- matched_gene_design_evidence[
    gene_symbol %chin% top_gene_designs
  ]
  plot_gene_design[, gene_symbol := factor(gene_symbol,
                                           levels = rev(top_gene_designs))]
  plot_gene_design[, design_family_plot :=
                     gsub("_", " ", design_family, fixed = TRUE)]
  plot_gene_design[, support_label := paste0(strict_supported_studies,
                                             "/", n_studies)]
  p_gene_design <- ggplot(
    plot_gene_design,
    aes(design_family_plot, gene_symbol,
        fill = weighted_clean_cds_pct_change)
  ) +
    geom_tile(color = "white", linewidth = 0.28) +
    geom_text(aes(label = support_label), size = 2.7) +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "#ffffff",
      high = "#8b0000",
      midpoint = 0,
      limits = c(-100, 400),
      oob = scales::squish,
      name = "weighted\nCDS %\n(capped)"
    ) +
    labs(
      title = "Post-viral matched gene evidence by design family",
      subtitle = "Fill is weighted matched clean-CDS proxy change; text is strict-supported studies / observed studies.",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 10.0) +
    theme(
      axis.text.x = element_text(angle = 28, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
  save_plot_pair(
    p_gene_design,
    "postviral_fatigue_matched_gene_design_evidence",
    width = 11.4,
    height = 9.5
  )
  save_plotly(
    p_gene_design,
    "postviral_fatigue_matched_gene_design_evidence"
  )
}

if (nrow(matched_gene_design_evidence) > 0 &&
    "loo_min_clean_cds_pct_change" %chin% names(matched_gene_design_evidence)) {
  plot_loo <- matched_gene_design_evidence[
    is.finite(loo_min_clean_cds_pct_change) &
      is.finite(loo_max_clean_cds_pct_change)
  ][seq_len(min(.N, 34L))]
  if (nrow(plot_loo) > 0) {
    plot_loo[, evidence_label := wrap_label(
      paste(
        gene_symbol,
        gsub("_", " ", design_family, fixed = TRUE),
        atlas_condition_evidence_gate,
        sep = " | "
      ),
      78L
    )]
    plot_loo[, evidence_label := factor(
      evidence_label,
      levels = rev(unique(evidence_label))
    )]
    plot_loo[, `:=`(
      full_pct_capped = pmin(pmax(weighted_clean_cds_pct_change, -100), 400),
      loo_min_pct_capped =
        pmin(pmax(loo_min_clean_cds_pct_change, -100), 400),
      loo_max_pct_capped =
        pmin(pmax(loo_max_clean_cds_pct_change, -100), 400)
    )]
    p_loo <- ggplot(plot_loo, aes(y = evidence_label)) +
      geom_vline(xintercept = 0, color = "#94a3b8", linewidth = 0.38) +
      geom_segment(
        aes(x = loo_min_pct_capped, xend = loo_max_pct_capped,
            yend = evidence_label, color = atlas_condition_evidence_gate),
        linewidth = 1.15,
        alpha = 0.84
      ) +
      geom_point(
        aes(x = full_pct_capped, color = atlas_condition_evidence_gate),
        size = 2.5
      ) +
      labs(
        title = "Post-viral matched clean-CDS leave-one-design-out stability",
        subtitle = "Point is full weighted effect; segment is the effect range after removing one matched design. Percent changes are capped for display.",
        x = "weighted clean-CDS proxy change (%)",
        y = NULL,
        color = "evidence gate"
      ) +
      theme_minimal(base_size = 9.9) +
      theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold")
      )
    save_plot_pair(
      p_loo,
      "postviral_fatigue_matched_gene_design_loo_stability",
      width = 13.2,
      height = 10.2
    )
    save_plotly(
      p_loo,
      "postviral_fatigue_matched_gene_design_loo_stability"
    )
  }
}

if (nrow(matched_gene_design_evidence) > 0 &&
    "random_effects_strict_pct_change_lower" %chin%
      names(matched_gene_design_evidence)) {
  plot_interval <- matched_gene_design_evidence[
    is.finite(random_effects_strict_pct_change_lower) &
      is.finite(random_effects_strict_pct_change_upper)
  ][seq_len(min(.N, 34L))]
  if (nrow(plot_interval) > 0) {
    plot_interval[, evidence_label := wrap_label(
      paste(
        gene_symbol,
        gsub("_", " ", design_family, fixed = TRUE),
        atlas_condition_evidence_gate,
        sep = " | "
      ),
      78L
    )]
    plot_interval[, evidence_label := factor(
      evidence_label,
      levels = rev(unique(evidence_label))
    )]
    plot_interval[, `:=`(
      interval_pct_capped =
        pmin(pmax(random_effects_strict_pct_change, -100), 400),
      interval_lower_pct_capped =
        pmin(pmax(random_effects_strict_pct_change_lower, -100), 400),
      interval_upper_pct_capped =
        pmin(pmax(random_effects_strict_pct_change_upper, -100), 400)
    )]
    p_interval <- ggplot(plot_interval, aes(y = evidence_label)) +
      geom_vline(xintercept = 0, color = "#94a3b8", linewidth = 0.38) +
      geom_segment(
        aes(x = interval_lower_pct_capped,
            xend = interval_upper_pct_capped,
            yend = evidence_label,
            color = atlas_condition_evidence_gate),
        linewidth = 1.15,
        alpha = 0.84
      ) +
      geom_point(
        aes(x = interval_pct_capped,
            fill = strict_random_effects_interval_supported),
        shape = 21,
        color = "#111827",
        size = 2.65,
        stroke = 0.3
      ) +
      scale_fill_manual(
        name = "strict interval\nsupported",
        values = c("FALSE" = "#ffffff", "TRUE" = "#111827")
      ) +
      labs(
        title = "Post-viral strict matched clean-CDS interval gate",
        subtitle = "Random-effects proxy intervals use strict matched bootstrap log2FC SEs; percent changes are capped for display.",
        x = "strict random-effects clean-CDS proxy change (%)",
        y = NULL,
        color = "evidence gate"
      ) +
      theme_minimal(base_size = 9.9) +
      theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold")
      )
    save_plot_pair(
      p_interval,
      "postviral_fatigue_matched_gene_design_strict_interval",
      width = 13.2,
      height = 10.2
    )
    save_plotly(
      p_interval,
      "postviral_fatigue_matched_gene_design_strict_interval"
    )
  }
}

validation_queue <- data.table()
if (nrow(candidate_rankings) > 0) {
  validation_queue <- candidate_rankings[, .SD[1], by = gene_symbol]
  if (nrow(feature_shift) > 0) {
    validation_queue <- merge(
      validation_queue,
      feature_shift[
        ,
        .(
          fatigue_max_feature_shift = max(max_abs_delta_relative_translation, na.rm = TRUE),
          fatigue_top_shift_feature = feature_id[which.max(max_abs_delta_relative_translation)][1],
          fatigue_top_shift_type = feature_type[which.max(max_abs_delta_relative_translation)][1]
        ),
        by = gene_symbol
      ],
      by = "gene_symbol",
      all.x = TRUE,
      sort = FALSE
    )
  }
  if (exists("counterfactual") && nrow(counterfactual) > 0) {
    validation_queue <- merge(
      validation_queue,
      counterfactual[
        ,
        .(
          top_counterfactual_family = design_family[which.max(abs(predicted_clean_cds_pct_change))][1],
          top_counterfactual_cds_pct = predicted_clean_cds_pct_change[which.max(abs(predicted_clean_cds_pct_change))][1],
          top_counterfactual_confidence = prediction_confidence[which.max(abs(predicted_clean_cds_pct_change))][1]
        ),
        by = gene_symbol
      ],
      by = "gene_symbol",
      all.x = TRUE,
      sort = FALSE
    )
  }
  if (exists("matched_viral_ifn_review") && nrow(matched_viral_ifn_review) > 0) {
    validation_queue <- merge(
      validation_queue,
      matched_viral_ifn_review[
        ,
        .(
          top_matched_postviral_family =
            design_family[which.max(matched_postviral_review_score)][1],
          top_matched_postviral_condition =
            case_condition[which.max(matched_postviral_review_score)][1],
          top_matched_postviral_cds_pct =
            matched_clean_cds_pct_change[which.max(matched_postviral_review_score)][1],
          top_matched_postviral_support =
            matched_cds_support[which.max(matched_postviral_review_score)][1],
          top_matched_postviral_review_score =
            max(matched_postviral_review_score, na.rm = TRUE)
        ),
        by = gene_symbol
      ],
      by = "gene_symbol",
      all.x = TRUE,
      sort = FALSE
    )
  }
  if (exists("matched_gene_design_summary") &&
      nrow(matched_gene_design_summary) > 0) {
    validation_queue <- merge(
      validation_queue,
      matched_gene_design_summary[
        ,
        .SD[1],
        by = gene_symbol
      ][
        ,
        .(
          gene_symbol,
          postviral_design_evidence_class,
          postviral_design_families,
          postviral_primary_design_families,
          postviral_designs,
          postviral_studies,
          postviral_strict_designs,
          postviral_strict_studies,
          postviral_weighted_clean_cds_pct_change,
          postviral_fraction_positive_effects,
          postviral_best_family,
          postviral_best_condition,
          postviral_best_study,
          postviral_design_replication_score,
          postviral_best_gene_design_family,
          postviral_best_gene_design_gate,
          postviral_best_gene_design_score,
          postviral_best_gene_design_loo_direction_stable,
          postviral_best_gene_design_loo_same_direction_fraction,
          postviral_best_gene_design_loo_top_weight_fraction,
          postviral_best_gene_design_loo_min_strict_designs,
          postviral_best_gene_design_study_loo_direction_stable,
          postviral_best_gene_design_study_loo_same_direction_fraction,
          postviral_best_gene_design_study_loo_top_weight_fraction,
          postviral_best_gene_design_strict_interval_supported,
          postviral_best_gene_design_strict_random_effects_pct_change,
          postviral_best_gene_design_strict_random_effects_pct_change_lower,
          postviral_best_gene_design_strict_random_effects_pct_change_upper
        )
      ],
      by = "gene_symbol",
      all.x = TRUE,
      sort = FALSE
    )
  }
  if (exists("iorf_candidates") && nrow(iorf_candidates) > 0) {
    validation_queue <- merge(
      validation_queue,
      iorf_candidates[
        ,
        .(
          best_iorf_evidence_class = iorf_evidence_class[which.max(iorf_evidence_score)][1],
          best_iorf_evidence_score = max(iorf_evidence_score, na.rm = TRUE)
        ),
        by = gene_symbol
      ],
      by = "gene_symbol",
      all.x = TRUE,
      sort = FALSE
    )
  }
  for (missing_col in c("fatigue_max_feature_shift",
                        "top_counterfactual_confidence",
                        "top_matched_postviral_review_score",
                        "postviral_design_replication_score",
                        "postviral_best_gene_design_score",
                        "best_iorf_evidence_class")) {
    if (!missing_col %chin% names(validation_queue)) {
      validation_queue[, (missing_col) := NA]
    }
  }
  validation_queue[, postviral_validation_score :=
                     fifelse(is.finite(fatigue_max_feature_shift), fatigue_max_feature_shift, 0) +
                     fifelse(is.finite(top_counterfactual_confidence), top_counterfactual_confidence, 0) +
                     fifelse(is.finite(top_matched_postviral_review_score),
                             pmin(top_matched_postviral_review_score, 1), 0) +
                     0.5 * fifelse(is.finite(postviral_design_replication_score),
                                    pmin(postviral_design_replication_score, 1), 0) +
                     0.35 * fifelse(is.finite(postviral_best_gene_design_score),
                                     pmin(postviral_best_gene_design_score, 1), 0) +
                     fifelse(best_iorf_evidence_class %chin% c("high_priority_iORF_review",
                                                               "medium_priority_iORF_review"),
                             0.5, 0)]
  setorder(validation_queue, -postviral_validation_score)
}
fwrite(validation_queue, file.path(output_dir, "postviral_fatigue_validation_queue.csv"))

if (nrow(fatigue_buffering) > 0) {
  top_buffering <- fatigue_buffering[
    is.finite(cds_fpkm_change_pct_stabilized)
  ][seq_len(min(.N, 40L))]
  if (nrow(top_buffering) > 0) {
    top_buffering[, buffering_label := wrap_label(
      paste0(gene_symbol, " | ", field, "=", level),
      62L
    )]
    top_buffering[, buffering_label := factor(buffering_label, levels = rev(unique(buffering_label)))]
    p_buffering <- ggplot(top_buffering, aes(cds_fpkm_change_pct_stabilized, buffering_label)) +
      geom_vline(xintercept = 0, color = "#cbd5e1", linewidth = 0.35) +
      geom_point(aes(color = buffering_interpretation,
                     size = buffering_confidence_score),
                 alpha = 0.9) +
      scale_size_continuous(name = "confidence", range = c(1.8, 5.5), limits = c(0, 1)) +
      labs(
        title = "Post-viral CDS buffering effects",
        subtitle = "Clean-CDS percent change in fatigue-linked RDG branches; stabilized denominator",
        x = "clean-CDS proxy change (%)",
        y = NULL,
        color = "buffering class"
      ) +
      theme_minimal(base_size = 10.5) +
      theme(plot.title = element_text(face = "bold"), legend.position = "right")
    save_plot_pair(p_buffering, "postviral_fatigue_cds_buffering_effects",
                   width = 12.6, height = 10.0)
    save_plotly(p_buffering, "postviral_fatigue_cds_buffering_effects")
  }
}

summary_metrics <- rbindlist(list(
  data.table(metric = "samples", value = nrow(state_dt)),
  data.table(metric = "fatigue_winner_samples",
             value = winner_counts[dominant_state == fatigue_state, N][1] %||% 0L),
  data.table(metric = "fatigue_score_q95",
             value = as.numeric(quantile(state_dt[[fatigue_state]], 0.95, na.rm = TRUE))),
  data.table(metric = "fatigue_metadata_terms",
             value = nrow(term_tables)),
  data.table(metric = "fatigue_rdg_feature_contrasts",
             value = nrow(fatigue_contrasts)),
  data.table(metric = "fatigue_internal_orf_candidate_rows",
             value = nrow(iorf_candidates)),
  data.table(metric = "fatigue_feature_shift_gene_count",
             value = uniqueN(feature_shift$gene_symbol)),
  data.table(metric = "fatigue_cds_buffering_rows",
             value = nrow(fatigue_buffering)),
  data.table(metric = "fatigue_submodule_count",
             value = if (exists("submodule_summary")) nrow(submodule_summary) else 0L),
  data.table(metric = "fatigue_counterfactual_prediction_rows",
             value = if (exists("counterfactual")) nrow(counterfactual) else 0L),
  data.table(metric = "fatigue_matched_viral_ifn_review_rows",
             value = nrow(matched_viral_ifn_review)),
  data.table(metric = "fatigue_matched_viral_ifn_strict_supported_rows",
             value = sum(matched_viral_ifn_review$clean_cds_strict_support == TRUE,
                         na.rm = TRUE)),
  data.table(metric = "fatigue_matched_design_cards",
             value = nrow(matched_design_cards)),
  data.table(metric = "fatigue_matched_design_cards_with_strict_support",
             value = matched_design_cards[n_strict_supported_genes > 0, .N]),
  data.table(metric = "fatigue_matched_gene_design_evidence_rows",
             value = nrow(matched_gene_design_evidence)),
  data.table(metric = "fatigue_matched_gene_design_replicated_strict_rows",
             value = matched_gene_design_evidence[
               gene_design_evidence_class == "replicated_strict_matched",
               .N
             ]),
  data.table(metric = "fatigue_matched_gene_design_loo_rows",
             value = nrow(matched_gene_design_loo)),
  data.table(metric = "fatigue_matched_gene_design_study_loo_rows",
             value = nrow(matched_gene_design_study_loo)),
  data.table(metric = "fatigue_matched_gene_design_strict_interval_supported_rows",
             value = matched_gene_design_evidence[
               strict_random_effects_interval_supported == TRUE,
               .N
             ]),
  data.table(metric = "fatigue_matched_gene_design_atlas_ready_rows",
             value = matched_gene_design_evidence[
               atlas_condition_evidence_gate ==
                 "atlas_ready_strict_design_study_interval",
               .N
             ]),
  data.table(metric = "fatigue_matched_gene_design_replicated_interval_review_rows",
             value = matched_gene_design_evidence[
               atlas_condition_evidence_gate ==
                 "replicated_strict_interval_review",
               .N
             ]),
  data.table(metric = "fatigue_validation_queue_rows",
             value = nrow(validation_queue))
), fill = TRUE)
fwrite(summary_metrics, file.path(output_dir, "postviral_fatigue_summary_metrics.csv"))

writeLines(c(
  "Post-viral fatigue pipeline summary",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste0("Samples: ", nrow(state_dt)),
  paste0("Fatigue winner samples: ", summary_metrics[metric == "fatigue_winner_samples", value]),
  paste0("Fatigue metadata terms: ", nrow(term_tables)),
  paste0("Fatigue RDG feature contrast rows: ", nrow(fatigue_contrasts)),
  paste0("Fatigue internal ORF candidate rows: ", nrow(iorf_candidates)),
  paste0("Fatigue CDS buffering rows: ", nrow(fatigue_buffering)),
  paste0("Fatigue counterfactual prediction rows: ",
         if (exists("counterfactual")) nrow(counterfactual) else 0L),
  paste0("Matched viral/IFN review rows: ", nrow(matched_viral_ifn_review)),
  paste0("Matched viral/IFN strict supported rows: ",
         sum(matched_viral_ifn_review$clean_cds_strict_support == TRUE,
             na.rm = TRUE)),
  paste0("Matched post-viral design cards: ", nrow(matched_design_cards)),
  paste0("Matched post-viral design cards with strict support: ",
         matched_design_cards[n_strict_supported_genes > 0, .N]),
  paste0("Matched gene x design-family evidence rows: ",
         nrow(matched_gene_design_evidence)),
  paste0("Replicated strict gene x design-family rows: ",
         matched_gene_design_evidence[
           gene_design_evidence_class == "replicated_strict_matched",
           .N
         ]),
  paste0("Matched gene x design-family leave-one-design rows: ",
         nrow(matched_gene_design_loo)),
  paste0("Matched gene x design-family leave-one-study rows: ",
         nrow(matched_gene_design_study_loo)),
  paste0("Strict random-effects interval-supported gene x design-family rows: ",
         matched_gene_design_evidence[
           strict_random_effects_interval_supported == TRUE,
           .N
         ]),
  paste0("Atlas-ready strict design/study/interval gene x design-family rows: ",
         matched_gene_design_evidence[
           atlas_condition_evidence_gate ==
             "atlas_ready_strict_design_study_interval",
           .N
         ]),
  paste0("Replicated strict interval-review gene x design-family rows: ",
         matched_gene_design_evidence[
           atlas_condition_evidence_gate ==
             "replicated_strict_interval_review",
           .N
         ]),
  paste0("Fatigue validation queue rows: ", nrow(validation_queue)),
  "",
  "Important interpretation caveat:",
  "Internal ORF rows are diagnostic features, not proof of independent internal initiation.",
  "For in-frame internal ORFs, canonical CDS translation can share the same downstream coverage.",
  "Candidate iORFs need frame/start-peak evidence or browser review before mechanistic claims."
), file.path(output_dir, "postviral_fatigue_summary.txt"))

message("Saved post-viral fatigue outputs in: ", output_dir)
print(summary_metrics)

invisible(list(
  fatigue_scores = fatigue_scores,
  metadata_terms = term_tables,
  feature_shift = feature_shift,
  iorf_candidates = iorf_candidates,
  matched_viral_ifn_review = matched_viral_ifn_review,
  matched_design_cards = matched_design_cards,
  matched_gene_design_evidence = matched_gene_design_evidence,
  matched_gene_design_study_loo = matched_gene_design_study_loo,
  summary_metrics = summary_metrics
))
