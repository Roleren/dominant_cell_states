#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (basename(here) == "dominant_cell_states" &&
        dir.exists(file.path(here, "scripts"))) {
      return(here)
    }
    candidate <- file.path(here, "dominant_cell_states")
    if (dir.exists(candidate) &&
        dir.exists(file.path(candidate, "scripts"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find dominant_cell_states analysis directory from: ",
       start, call. = FALSE)
}

analysis_dir <- find_analysis_dir()
output_dir <- file.path(analysis_dir, "dominant_rdg_model_agreement_transport")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)

metadata_file <- "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv"

read_dt <- function(path, required = FALSE, ...) {
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("Missing required input: ", path, call. = FALSE)
    return(data.table())
  }
  fread(path, showProgress = FALSE, ...)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

finite_or_zero <- function(x) {
  x <- safe_num(x)
  x[!is.finite(x)] <- 0
  x
}

clip01 <- function(x) pmin(pmax(finite_or_zero(x), 0), 1)

min_safe <- function(x, default = NA_real_) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) default else min(x)
}

max_safe <- function(x, default = NA_real_) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) default else max(x)
}

median_safe <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else median(x)
}

collapse_unique <- function(x, n = 8L) {
  x <- as.character(x[!is.na(x) & nzchar(as.character(x))])
  x <- trimws(unlist(strsplit(x, ";\\s*", perl = TRUE), use.names = FALSE))
  x <- sort(unique(x[nzchar(x)]))
  if (!length(x)) return("")
  paste(head(x, n), collapse = "; ")
}

mode_value <- function(x, default = "unknown") {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x) & !tolower(x) %chin%
           c("na", "n/a", "unknown")]
  if (!length(x)) return(default)
  counts <- sort(table(x), decreasing = TRUE)
  names(counts)[1]
}

slug_value <- function(x, default = "unknown") {
  x <- tolower(trimws(as.character(x)))
  x[is.na(x) | !nzchar(x) | x %chin% c("na", "n/a")] <- default
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[!nzchar(x)] <- default
  x
}

random_effects_summary <- function(effect, se) {
  effect <- safe_num(effect)
  se <- safe_num(se)
  keep <- is.finite(effect) & is.finite(se) & se > 0
  if (!any(keep)) {
    return(data.table(
      n = 0L, effect = NA_real_, se = NA_real_, lower = NA_real_,
      upper = NA_real_, p = NA_real_, tau2 = NA_real_, i2 = NA_real_,
      direction_fraction = NA_real_, ci_excludes_zero = FALSE
    ))
  }
  yi <- effect[keep]
  vi <- pmax(se[keep], 1e-6)^2
  wi <- 1 / vi
  fixed <- sum(wi * yi) / sum(wi)
  q_stat <- sum(wi * (yi - fixed)^2)
  df <- length(yi) - 1
  c_val <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (df > 0 && c_val > 0) max(0, (q_stat - df) / c_val) else 0
  wi_re <- 1 / (vi + tau2)
  mu <- sum(wi_re * yi) / sum(wi_re)
  se_mu <- sqrt(1 / sum(wi_re))
  lower <- mu - 1.96 * se_mu
  upper <- mu + 1.96 * se_mu
  direction_fraction <- if (sign(mu) == 0) {
    mean(yi == 0)
  } else {
    mean(sign(yi) == sign(mu))
  }
  data.table(
    n = length(yi),
    effect = mu,
    se = se_mu,
    lower = lower,
    upper = upper,
    p = 2 * pnorm(-abs(mu / se_mu)),
    tau2 = tau2,
    i2 = if (q_stat > 0 && df > 0) max(0, (q_stat - df) / q_stat) else 0,
    direction_fraction = direction_fraction,
    ci_excludes_zero = isTRUE(lower > 0 | upper < 0)
  )
}

classify_transport_support <- function(model, n, effect, p, ci,
                                       direction_fraction) {
  if (n < 2 || !is.finite(effect) || !is.finite(p) ||
      !is.finite(direction_fraction)) {
    return("not_evaluable")
  }
  if (model == "hierarchical_joint") {
    return(fcase(
      ci & abs(effect) >= 0.04 & direction_fraction >= 0.70,
      "replicated_shift",
      p <= 0.20 & abs(effect) >= 0.03 & direction_fraction >= 0.65,
      "replicated_review",
      default = "neutral"
    ))
  }
  fcase(
    ci & p <= 0.20 & abs(effect) >= 0.035 & direction_fraction >= 0.70,
    "replicated_shift",
    (ci | p <= 0.30) & abs(effect) >= 0.025 &
      direction_fraction >= 0.65,
    "replicated_review",
    default = "neutral"
  )
}

standardize_hierarchical <- function(dt) {
  out <- copy(dt)
  setnames(out, c(
    "study_delta", "study_delta_se", "study_delta_p",
    "study_context_rows", "study_strict_context_rows",
    "study_max_context_priority", "study_top_context",
    "study_top_case_condition", "study_top_control_conditions",
    "study_top_cell_line", "study_top_tissue", "study_top_match_scope"
  ), c(
    "model_delta", "model_delta_se", "model_delta_p",
    "model_context_rows", "model_strict_context_rows",
    "model_max_context_priority", "model_top_context",
    "model_top_case_condition", "model_top_control_conditions",
    "model_top_cell_line", "model_top_tissue", "model_top_match_scope"
  ), skip_absent = TRUE)
  out[, model := "hierarchical_joint"]
  out
}

standardize_dm <- function(dt) {
  out <- copy(dt)
  setnames(out, c(
    "dm_study_delta", "dm_study_delta_se", "dm_study_delta_p",
    "dm_study_context_rows", "dm_study_strict_context_rows",
    "dm_study_max_priority", "dm_study_top_context",
    "dm_study_top_match_scope", "dm_study_top_prior_source"
  ), c(
    "model_delta", "model_delta_se", "model_delta_p",
    "model_context_rows", "model_strict_context_rows",
    "model_max_context_priority", "model_top_context",
    "model_top_match_scope", "model_top_prior_source"
  ), skip_absent = TRUE)
  out[, model := "dirichlet_multinomial"]
  out
}

message("Grouped-omission RDG transportability analysis:")
message("  1. Build publication/lab/protocol/context families per study")
message("  2. Omit groups containing at least two studies")
message("  3. Require at least two independent studies to remain")

agreement <- read_dt(
  file.path(analysis_dir, "dominant_rdg_model_agreement_loo",
            "dominant_rdg_model_agreement_loo_gene_design.csv"),
  required = TRUE
)
branch_comparison <- read_dt(
  file.path(analysis_dir, "dominant_rdg_model_agreement",
            "dominant_rdg_model_agreement_branch_comparison.csv"),
  required = TRUE
)
hierarchical_studies <- standardize_hierarchical(read_dt(
  file.path(analysis_dir, "dominant_rdg_hierarchical_joint_allocation",
            "dominant_rdg_hierarchical_joint_allocation_study_effects.csv"),
  required = TRUE
))
dm_studies <- standardize_dm(read_dt(
  file.path(analysis_dir, "dominant_rdg_dirichlet_multinomial",
            "dominant_rdg_dirichlet_multinomial_study_effects.csv"),
  required = TRUE
))
joint_rows <- read_dt(
  file.path(analysis_dir, "dominant_rdg_joint_branch_allocation",
            "dominant_rdg_joint_branch_allocation_rows.csv"),
  required = TRUE,
  select = c("study", "AUTHOR", "CELL_LINE", "TISSUE")
)

key <- c("gene_symbol", "tx_id", "design_family")
branch_key <- c(key, "branch_class")
axis_names <- c(
  "publication_family", "author_family", "protocol_family",
  "inhibitor_family", "tissue_context", "cell_line_context",
  "footprint_length_family"
)

joint_study <- joint_rows[, .(
  joint_author = mode_value(AUTHOR),
  joint_cell_line = mode_value(CELL_LINE),
  joint_tissue = mode_value(TISSUE)
), by = study]

metadata_columns <- c(
  "study", "Study_Pubmed_id", "YEAR", "AUTHOR", "CELL_LINE", "TISSUE",
  "INHIBITOR", "LibrarySelection", "LibraryLayout", "top_readlength"
)
metadata_header <- read_dt(metadata_file, required = TRUE, nrows = 0)
metadata <- read_dt(
  metadata_file,
  required = TRUE,
  select = intersect(metadata_columns, names(metadata_header))
)
for (column in setdiff(metadata_columns, names(metadata))) {
  metadata[, (column) := NA]
}

study_metadata <- metadata[, {
  pmid <- safe_num(Study_Pubmed_id)
  pmid <- pmid[is.finite(pmid) & pmid >= 1e6]
  top_len <- safe_num(top_readlength)
  top_len <- top_len[is.finite(top_len) & top_len > 0]
  data.table(
    study_pubmed_id = if (length(pmid)) as.character(names(
      sort(table(pmid), decreasing = TRUE)
    )[1]) else "",
    study_year = mode_value(YEAR),
    metadata_author = mode_value(AUTHOR),
    metadata_cell_line = mode_value(CELL_LINE),
    metadata_tissue = mode_value(TISSUE),
    inhibitor = mode_value(INHIBITOR),
    library_selection = mode_value(LibrarySelection),
    library_layout = mode_value(LibraryLayout),
    median_top_readlength = if (length(top_len)) median(top_len) else NA_real_,
    metadata_runs = .N,
    metadata_author_values = uniqueN(AUTHOR[!is.na(AUTHOR) & nzchar(AUTHOR)]),
    metadata_cell_line_values =
      uniqueN(CELL_LINE[!is.na(CELL_LINE) & nzchar(CELL_LINE)]),
    metadata_tissue_values =
      uniqueN(TISSUE[!is.na(TISSUE) & nzchar(TISSUE)]),
    metadata_inhibitor_values =
      uniqueN(INHIBITOR[!is.na(INHIBITOR) & nzchar(INHIBITOR)])
  )
}, by = study]
study_metadata <- merge(study_metadata, joint_study, by = "study", all = TRUE,
                        sort = FALSE)

study_metadata[, author := fifelse(
  metadata_author != "unknown", metadata_author, joint_author
)]
study_metadata[, cell_line := fifelse(
  metadata_cell_line != "unknown", metadata_cell_line, joint_cell_line
)]
study_metadata[, tissue := fifelse(
  metadata_tissue != "unknown", metadata_tissue, joint_tissue
)]
study_metadata[, publication_family := fifelse(
  nzchar(study_pubmed_id), paste0("pubmed:", study_pubmed_id), "unknown"
)]
study_metadata[, author_family := paste0("author:", slug_value(author))]
study_metadata[author == "unknown", author_family := "unknown"]
study_metadata[, inhibitor_family := fcase(
  grepl("chx|cycloheximide", inhibitor, ignore.case = TRUE),
  "elongation_chx",
  grepl("harr|lactim|ltm", inhibitor, ignore.case = TRUE),
  "initiation_stall",
  grepl("frozen|none|unknown", inhibitor, ignore.case = TRUE),
  "no_chx_or_frozen",
  default = "other_inhibitor"
)]
study_metadata[, selection_family := fcase(
  grepl("size|mnase", library_selection, ignore.case = TRUE),
  "size_fractionation",
  tolower(library_selection) %chin%
    c("other", "cdna", "unknown", "unspecified"),
  "standard_or_unspecified",
  default = "specialized_selection"
)]
study_metadata[, layout_family := fifelse(
  tolower(library_layout) == "paired", "paired", "single_or_unknown"
)]
study_metadata[, protocol_family := paste(
  inhibitor_family, selection_family, layout_family, sep = "|"
)]
study_metadata[, tissue_context := paste0("tissue:", slug_value(tissue))]
study_metadata[tissue == "unknown", tissue_context := "unknown"]
study_metadata[, cell_line_context := paste0(
  "cell_line:", slug_value(cell_line)
)]
study_metadata[cell_line == "unknown", cell_line_context := "unknown"]
study_metadata[, footprint_length_family := fcase(
  !is.finite(median_top_readlength), "unknown",
  median_top_readlength <= 29, "footprint_short_le29",
  median_top_readlength <= 32, "footprint_standard_30_32",
  default = "footprint_long_ge33"
)]

study_axis_long <- melt(
  study_metadata,
  id.vars = c(
    "study", "study_pubmed_id", "study_year", "author", "cell_line", "tissue",
    "inhibitor", "library_selection", "library_layout",
    "median_top_readlength", "metadata_runs", "metadata_author_values",
    "metadata_cell_line_values", "metadata_tissue_values",
    "metadata_inhibitor_values"
  ),
  measure.vars = axis_names,
  variable.name = "transport_axis",
  value.name = "transport_group"
)
study_axis_long <- study_axis_long[
  !is.na(transport_group) & nzchar(transport_group) &
    transport_group != "unknown"
]
study_axis_long[, transport_group_studies := uniqueN(study),
                by = .(transport_axis, transport_group)]
setorder(study_axis_long, transport_axis, -transport_group_studies,
         transport_group, study)

studies <- rbindlist(
  list(hierarchical_studies, dm_studies),
  use.names = TRUE,
  fill = TRUE
)
for (column in c("model_delta", "model_delta_se", "model_delta_p")) {
  studies[, (column) := safe_num(get(column))]
}
studies <- studies[
  is.finite(model_delta) & is.finite(model_delta_se) & model_delta_se > 0
]
studies <- merge(
  studies,
  study_metadata[, c("study", axis_names), with = FALSE],
  by = "study",
  all.x = TRUE,
  sort = FALSE
)
for (axis in axis_names) {
  studies[is.na(get(axis)) | !nzchar(get(axis)), (axis) := "unknown"]
}
setorder(studies, gene_symbol, tx_id, design_family, branch_class, model, study)

run_transport_omissions <- function(dt) {
  dt[, {
    full <- random_effects_summary(model_delta, model_delta_se)
    pieces <- lapply(axis_names, function(axis) {
      values <- as.character(.SD[[axis]])
      groups <- sort(unique(values[values != "unknown"]))
      if (!length(groups)) return(NULL)
      rbindlist(lapply(groups, function(group) {
        omit <- values == group
        omitted_studies <- uniqueN(study[omit])
        remaining <- .SD[!omit]
        remaining_studies <- uniqueN(remaining$study)
        if (omitted_studies < 2 || remaining_studies < 2) return(NULL)
        re <- random_effects_summary(
          remaining$model_delta,
          remaining$model_delta_se
        )
        support_class <- classify_transport_support(
          model[1], re$n, re$effect, re$p, re$ci_excludes_zero,
          re$direction_fraction
        )
        data.table(
          transport_axis = axis,
          omitted_group = group,
          model_studies = uniqueN(study),
          omitted_studies = omitted_studies,
          remaining_studies = remaining_studies,
          omitted_study_names = collapse_unique(study[omit], 12L),
          full_delta = full$effect,
          transport_delta = re$effect,
          transport_delta_se = re$se,
          transport_delta_lower = re$lower,
          transport_delta_upper = re$upper,
          transport_delta_p = re$p,
          transport_direction_fraction = re$direction_fraction,
          transport_ci_excludes_zero = re$ci_excludes_zero,
          transport_same_direction =
            is.finite(re$effect) & sign(re$effect) == sign(full$effect) &
              sign(full$effect) != 0,
          transport_abs_change = abs(re$effect - full$effect),
          transport_support_class = support_class,
          transport_support_retained =
            support_class %chin% c("replicated_shift", "replicated_review")
        )
      }), use.names = TRUE, fill = TRUE)
    })
    rbindlist(pieces, use.names = TRUE, fill = TRUE)
  }, by = c(branch_key, "model")]
}

transport_rows <- run_transport_omissions(studies)
setorder(transport_rows, gene_symbol, design_family, branch_class, model,
         transport_axis, -transport_abs_change)

transport_axis_model <- transport_rows[, .(
  transport_group_tests = .N,
  transport_model_studies = max_safe(model_studies, 0),
  transport_min_remaining_studies = min_safe(remaining_studies, 0),
  transport_fraction_same_direction =
    mean(transport_same_direction, na.rm = TRUE),
  transport_fraction_support_retained =
    mean(transport_support_retained, na.rm = TRUE),
  transport_max_abs_change = max_safe(transport_abs_change, 0),
  transport_worst_omitted_group =
    omitted_group[which.max(transport_abs_change)][1],
  transport_worst_omitted_studies =
    omitted_study_names[which.max(transport_abs_change)][1]
), by = c(branch_key, "model", "transport_axis")]
transport_axis_model[, transport_model_axis_class := fcase(
  transport_fraction_same_direction >= 0.80 &
    transport_fraction_support_retained >= 0.67,
  "transport_support_robust",
  transport_fraction_same_direction >= 0.80,
  "transport_direction_robust_support_fragile",
  default = "transport_direction_fragile"
)]

support_pairs <- rbindlist(list(
  branch_comparison[hierarchical_joint_supported == TRUE,
                    ..branch_key][, model := "hierarchical_joint"],
  branch_comparison[dm_supported == TRUE,
                    ..branch_key][, model := "dirichlet_multinomial"]
), use.names = TRUE, fill = TRUE)
support_pairs <- unique(support_pairs)
supported_axis_grid <- support_pairs[, .(transport_axis = axis_names),
                                     by = c(branch_key, "model")]
supported_axis_grid <- merge(
  supported_axis_grid,
  transport_axis_model,
  by = c(branch_key, "model", "transport_axis"),
  all.x = TRUE,
  sort = FALSE
)
supported_axis_grid[, transport_model_axis_evaluable :=
  is.finite(transport_fraction_same_direction)]
supported_axis_grid[
  transport_model_axis_evaluable == FALSE,
  transport_model_axis_class := "transport_not_evaluable"
]

transport_branch_axis <- supported_axis_grid[, .(
  transport_supported_models = .N,
  transport_evaluable_models = sum(transport_model_axis_evaluable),
  transport_group_tests = sum(transport_group_tests, na.rm = TRUE),
  transport_min_direction_fraction =
    min_safe(transport_fraction_same_direction),
  transport_min_support_retention =
    min_safe(transport_fraction_support_retained),
  transport_min_remaining_studies =
    min_safe(transport_min_remaining_studies),
  transport_max_abs_change = max_safe(transport_max_abs_change, 0),
  transport_worst_omitted_groups =
    collapse_unique(transport_worst_omitted_group, 4L),
  transport_worst_omitted_studies =
    collapse_unique(transport_worst_omitted_studies, 8L),
  transport_direction_fragile_models =
    sum(transport_model_axis_class == "transport_direction_fragile"),
  transport_support_robust_models =
    sum(transport_model_axis_class == "transport_support_robust")
), by = c(branch_key, "transport_axis")]
transport_branch_axis[, transport_branch_axis_class := fcase(
  transport_evaluable_models == 0, "transport_not_evaluable",
  transport_direction_fragile_models > 0, "transport_direction_fragile",
  transport_evaluable_models < transport_supported_models,
  "transport_partially_evaluable",
  transport_support_robust_models == transport_supported_models,
  "transport_support_robust",
  default = "transport_direction_robust_support_fragile"
)]
setorder(transport_branch_axis, gene_symbol, design_family, branch_class,
         transport_axis)

transport_branch <- transport_branch_axis[, .(
  transport_evaluable_axes =
    sum(transport_branch_axis_class != "transport_not_evaluable"),
  transport_support_robust_axes =
    sum(transport_branch_axis_class == "transport_support_robust"),
  transport_direction_robust_support_fragile_axes =
    sum(transport_branch_axis_class ==
          "transport_direction_robust_support_fragile"),
  transport_partially_evaluable_axes =
    sum(transport_branch_axis_class == "transport_partially_evaluable"),
  transport_direction_fragile_axes =
    sum(transport_branch_axis_class == "transport_direction_fragile"),
  transport_min_direction_fraction =
    min_safe(transport_min_direction_fraction),
  transport_min_support_retention =
    min_safe(transport_min_support_retention),
  transport_min_remaining_studies =
    min_safe(transport_min_remaining_studies),
  transport_evaluable_axis_names = collapse_unique(
    transport_axis[transport_branch_axis_class != "transport_not_evaluable"],
    7L
  ),
  transport_support_robust_axis_names = collapse_unique(
    transport_axis[transport_branch_axis_class == "transport_support_robust"],
    7L
  ),
  transport_fragile_axis_names = collapse_unique(
    transport_axis[transport_branch_axis_class %chin% c(
      "transport_direction_fragile",
      "transport_direction_robust_support_fragile"
    )],
    7L
  ),
  transport_worst_omitted_groups =
    collapse_unique(transport_worst_omitted_groups, 8L),
  transport_worst_omitted_studies =
    collapse_unique(transport_worst_omitted_studies, 8L)
), by = branch_key]
transport_branch[, transport_branch_class := fcase(
  transport_evaluable_axes == 0, "transport_not_evaluable",
  transport_direction_fragile_axes > 0, "transport_direction_fragile",
  transport_support_robust_axes == transport_evaluable_axes,
  "transport_support_robust",
  default = "transport_direction_robust_support_fragile"
)]

relevant_branches <- branch_comparison[
  hierarchical_joint_supported == TRUE | dm_supported == TRUE,
  ..branch_key
]
relevant_branches <- unique(relevant_branches)
transport_branch <- merge(
  relevant_branches,
  transport_branch,
  by = branch_key,
  all.x = TRUE,
  sort = FALSE
)
for (column in c(
  "transport_evaluable_axes", "transport_support_robust_axes",
  "transport_direction_robust_support_fragile_axes",
  "transport_partially_evaluable_axes", "transport_direction_fragile_axes"
)) {
  transport_branch[is.na(get(column)), (column) := 0]
}
transport_branch[
  is.na(transport_branch_class),
  transport_branch_class := "transport_not_evaluable"
]

gene_transport <- transport_branch[, .(
  transport_relevant_supported_branches = .N,
  transport_evaluable_supported_branches =
    sum(transport_branch_class != "transport_not_evaluable"),
  transport_evaluable_axis_tests = sum(transport_evaluable_axes),
  transport_possible_axis_tests = .N * length(axis_names),
  transport_min_evaluable_axes = min_safe(transport_evaluable_axes, 0),
  transport_support_robust_branches =
    sum(transport_branch_class == "transport_support_robust"),
  transport_direction_robust_support_fragile_branches =
    sum(transport_branch_class ==
          "transport_direction_robust_support_fragile"),
  transport_direction_fragile_branches =
    sum(transport_branch_class == "transport_direction_fragile"),
  transport_not_evaluable_branches =
    sum(transport_branch_class == "transport_not_evaluable"),
  transport_min_direction_fraction =
    min_safe(transport_min_direction_fraction, 0),
  transport_min_support_retention =
    min_safe(transport_min_support_retention, 0),
  transport_min_remaining_studies =
    min_safe(transport_min_remaining_studies, 0),
  transport_evaluable_axis_names =
    collapse_unique(transport_evaluable_axis_names, 7L),
  transport_support_robust_axis_names =
    collapse_unique(transport_support_robust_axis_names, 7L),
  transport_fragile_axis_names =
    collapse_unique(transport_fragile_axis_names, 7L),
  transport_worst_omitted_groups =
    collapse_unique(transport_worst_omitted_groups, 8L),
  transport_worst_omitted_studies =
    collapse_unique(transport_worst_omitted_studies, 8L)
), by = key]

agreement_transport <- merge(
  agreement,
  gene_transport,
  by = key,
  all.x = TRUE,
  sort = FALSE
)
for (column in c(
  "transport_relevant_supported_branches",
  "transport_evaluable_supported_branches",
  "transport_evaluable_axis_tests",
  "transport_possible_axis_tests",
  "transport_min_evaluable_axes",
  "transport_support_robust_branches",
  "transport_direction_robust_support_fragile_branches",
  "transport_direction_fragile_branches",
  "transport_not_evaluable_branches",
  "transport_min_direction_fraction",
  "transport_min_support_retention",
  "transport_min_remaining_studies"
)) {
  agreement_transport[is.na(get(column)), (column) := 0]
}
for (column in c(
  "transport_evaluable_axis_names", "transport_support_robust_axis_names",
  "transport_fragile_axis_names", "transport_worst_omitted_groups",
  "transport_worst_omitted_studies"
)) {
  agreement_transport[is.na(get(column)), (column) := ""]
}
agreement_transport[, transport_coverage_fraction := fifelse(
  transport_relevant_supported_branches > 0,
  transport_evaluable_supported_branches /
    transport_relevant_supported_branches,
  0
)]
agreement_transport[, transport_axis_coverage_fraction := fifelse(
  transport_possible_axis_tests > 0,
  transport_evaluable_axis_tests / transport_possible_axis_tests,
  0
)]
agreement_transport[, transport_agreement_class := fcase(
  transport_evaluable_supported_branches == 0, "transport_not_evaluable",
  transport_direction_fragile_branches > 0, "transport_direction_fragile",
  transport_evaluable_supported_branches <
    transport_relevant_supported_branches,
  "transport_partially_evaluable",
  transport_support_robust_branches ==
    transport_evaluable_supported_branches,
  "transport_support_robust",
  default = "transport_direction_robust_support_fragile"
)]
consensus_tiers <- c(
  "A_browser_validated_consensus",
  "B_replicated_branch_consensus",
  "C_multi_model_consensus"
)
agreement_transport[, transport_consensus_class := fcase(
  model_agreement_tier %chin% consensus_tiers &
    transport_agreement_class == "transport_support_robust",
  "consensus_transport_support_robust",
  model_agreement_tier %chin% consensus_tiers &
    transport_agreement_class == "transport_direction_robust_support_fragile",
  "consensus_transport_direction_robust",
  model_agreement_tier %chin% consensus_tiers &
    transport_agreement_class == "transport_partially_evaluable",
  "consensus_transport_partially_evaluable",
  model_agreement_tier %chin% consensus_tiers &
    transport_agreement_class == "transport_not_evaluable",
  "consensus_transport_not_evaluable",
  model_agreement_tier %chin% consensus_tiers,
  "consensus_transport_fragile",
  model_agreement_tier == "D_model_disagreement_review",
  "disagreement_review",
  default = "non_consensus_review"
)]
agreement_transport[, transport_stability_score := clip01(
  0.45 * transport_min_direction_fraction +
    0.30 * transport_min_support_retention +
    0.08 * transport_coverage_fraction +
    0.07 * transport_axis_coverage_fraction +
    0.10 * pmin(transport_min_remaining_studies / 4, 1)
)]
agreement_transport[
  transport_agreement_class == "transport_not_evaluable",
  transport_stability_score := 0
]
agreement_transport[, transport_review_priority := pmin(
  1.5,
  loo_review_priority +
    0.12 * transport_stability_score +
    0.10 * as.numeric(transport_consensus_class ==
                        "consensus_transport_fragile") +
    0.06 * as.numeric(transport_consensus_class ==
                        "consensus_transport_partially_evaluable")
)]
setorder(agreement_transport, -transport_review_priority,
         -transport_stability_score, gene_symbol, design_family)

gene_axis <- transport_branch_axis[, .(
  transport_axis_evaluable_branches =
    sum(transport_branch_axis_class != "transport_not_evaluable"),
  transport_axis_support_robust_branches =
    sum(transport_branch_axis_class == "transport_support_robust"),
  transport_axis_direction_fragile_branches =
    sum(transport_branch_axis_class == "transport_direction_fragile"),
  transport_axis_min_direction_fraction =
    min_safe(transport_min_direction_fraction),
  transport_axis_min_support_retention =
    min_safe(transport_min_support_retention),
  transport_axis_min_remaining_studies =
    min_safe(transport_min_remaining_studies),
  transport_axis_worst_omitted_groups =
    collapse_unique(transport_worst_omitted_groups, 5L)
), by = c(key, "transport_axis")]
gene_axis[, transport_axis_class := fcase(
  transport_axis_evaluable_branches == 0, "transport_not_evaluable",
  transport_axis_direction_fragile_branches > 0,
  "transport_direction_fragile",
  transport_axis_support_robust_branches ==
    transport_axis_evaluable_branches,
  "transport_support_robust",
  default = "transport_direction_robust_support_fragile"
)]

viral_review <- agreement_transport[design_family == "viral_infection"]
fragile_review <- agreement_transport[
  transport_agreement_class %chin% c(
    "transport_direction_fragile",
    "transport_direction_robust_support_fragile",
    "transport_partially_evaluable"
  )
]

summarize_scope <- function(dt, scope) {
  consensus <- dt$model_agreement_tier %chin% consensus_tiers
  evaluable <- dt$transport_evaluable_supported_branches > 0
  data.table(
    scope = scope,
    transport_rule =
      "omit_group_ge2_studies_require_ge2_remaining_within_branch_p_proxy",
    gene_designs = nrow(dt),
    consensus_gene_designs = sum(consensus),
    transport_evaluable_gene_designs = sum(evaluable),
    consensus_transport_evaluable_gene_designs = sum(consensus & evaluable),
    consensus_transport_support_robust =
      sum(dt$transport_consensus_class ==
            "consensus_transport_support_robust"),
    consensus_transport_direction_robust =
      sum(dt$transport_consensus_class ==
            "consensus_transport_direction_robust"),
    consensus_transport_partially_evaluable =
      sum(dt$transport_consensus_class ==
            "consensus_transport_partially_evaluable"),
    consensus_transport_not_evaluable =
      sum(dt$transport_consensus_class ==
            "consensus_transport_not_evaluable"),
    consensus_transport_fragile =
      sum(dt$transport_consensus_class == "consensus_transport_fragile"),
    median_evaluable_transport_stability_score =
      median_safe(dt$transport_stability_score[evaluable]),
    median_evaluable_transport_min_direction_fraction =
      median_safe(dt$transport_min_direction_fraction[evaluable]),
    median_evaluable_transport_min_support_retention =
      median_safe(dt$transport_min_support_retention[evaluable]),
    median_evaluable_transport_axis_coverage =
      median_safe(dt$transport_axis_coverage_fraction[evaluable])
  )
}
summary_metrics <- rbindlist(list(
  summarize_scope(agreement_transport, "all"),
  summarize_scope(viral_review, "design_family:viral_infection")
))

fwrite(study_metadata, file.path(
  output_dir, "dominant_rdg_model_agreement_transport_study_metadata.csv"
))
fwrite(study_axis_long, file.path(
  output_dir, "dominant_rdg_model_agreement_transport_study_groups.csv"
))
fwrite(transport_rows, file.path(
  output_dir, "dominant_rdg_model_agreement_transport_omission_rows.csv"
))
fwrite(transport_axis_model, file.path(
  output_dir, "dominant_rdg_model_agreement_transport_axis_model_summary.csv"
))
fwrite(transport_branch_axis, file.path(
  output_dir, "dominant_rdg_model_agreement_transport_branch_axis_summary.csv"
))
fwrite(transport_branch, file.path(
  output_dir, "dominant_rdg_model_agreement_transport_branch_summary.csv"
))
fwrite(gene_axis, file.path(
  output_dir, "dominant_rdg_model_agreement_transport_gene_axis_summary.csv"
))
fwrite(agreement_transport, file.path(
  output_dir, "dominant_rdg_model_agreement_transport_gene_design.csv"
))
fwrite(viral_review, file.path(
  output_dir, "dominant_rdg_model_agreement_transport_viral_review.csv"
))
fwrite(fragile_review, file.path(
  output_dir, "dominant_rdg_model_agreement_transport_fragile_review.csv"
))
fwrite(summary_metrics, file.path(
  output_dir, "dominant_rdg_model_agreement_transport_summary_metrics.csv"
))

viral_consensus <- viral_review[model_agreement_tier %chin% consensus_tiers]
viral_axis_plot <- merge(
  viral_consensus[, ..key],
  gene_axis,
  by = key,
  all.x = TRUE,
  sort = FALSE
)
viral_axis_plot <- viral_axis_plot[is.finite(
  transport_axis_min_direction_fraction
)]
if (nrow(viral_axis_plot)) {
  viral_axis_plot[, gene_symbol := factor(
    gene_symbol,
    levels = rev(unique(viral_consensus$gene_symbol))
  )]
  viral_axis_plot[, transport_axis := factor(
    transport_axis,
    levels = axis_names
  )]
  p_axis <- ggplot(
    viral_axis_plot,
    aes(
      x = transport_axis,
      y = gene_symbol,
      fill = transport_axis_min_support_retention,
      text = paste0(
        "Gene: ", gene_symbol,
        "<br>Omission axis: ", transport_axis,
        "<br>Class: ", transport_axis_class,
        "<br>Minimum direction retention: ",
        round(transport_axis_min_direction_fraction, 3),
        "<br>Minimum support retention: ",
        round(transport_axis_min_support_retention, 3),
        "<br>Minimum remaining studies: ",
        transport_axis_min_remaining_studies,
        "<br>Worst omitted groups: ", transport_axis_worst_omitted_groups
      )
    )
  ) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(
      aes(label = sprintf("%.0f%%", 100 * transport_axis_min_direction_fraction)),
      size = 3.2
    ) +
    scale_fill_gradient2(
      low = "#b2473e", mid = "white", high = "#245b9e",
      midpoint = 0.67, limits = c(0, 1),
      name = "Minimum support\nretention"
    ) +
    labs(
      title = "Viral Consensus Transportability Across Grouped Omissions",
      subtitle = "Tile fill is support retention; label is branch-direction retention",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  ggsave(file.path(figure_dir, "model_agreement_transport_viral_axes.png"),
         p_axis, width = 11.5, height = 5.8, dpi = 180)
  ggsave(file.path(figure_dir, "model_agreement_transport_viral_axes.pdf"),
         p_axis, width = 11.5, height = 5.8)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_axis,
      file.path(figure_dir, "model_agreement_transport_viral_axes.html"),
      title = "Viral consensus grouped-omission transportability"
    )
  }
}

consensus_plot <- agreement_transport[
  model_agreement_tier %chin% c(consensus_tiers, "D_model_disagreement_review")
]
consensus_plot <- consensus_plot[transport_evaluable_supported_branches > 0]
if (nrow(consensus_plot)) {
  label_dt <- consensus_plot[
    design_family == "viral_infection" |
      transport_consensus_class == "consensus_transport_fragile"
  ][order(-transport_review_priority)][1:min(.N, 24)]
  label_dt[, transport_label := fifelse(
    design_family == "viral_infection",
    gene_symbol,
    paste(gene_symbol, design_family, sep = " / ")
  )]
  p_consensus <- ggplot(
    consensus_plot,
    aes(
      x = transport_min_direction_fraction,
      y = transport_min_support_retention,
      color = transport_consensus_class,
      size = pmax(model_agreement_score, 0.05),
      text = paste0(
        "Gene/design: ", gene_symbol, " / ", design_family,
        "<br>Transport class: ", transport_consensus_class,
        "<br>Minimum direction retention: ",
        round(transport_min_direction_fraction, 3),
        "<br>Minimum support retention: ",
        round(transport_min_support_retention, 3),
        "<br>Coverage fraction: ", round(transport_coverage_fraction, 3),
        "<br>Axis coverage fraction: ",
        round(transport_axis_coverage_fraction, 3),
        "<br>Worst omitted groups: ", transport_worst_omitted_groups
      )
    )
  ) +
    geom_vline(xintercept = 0.80, color = "#98a2b3", linetype = "dashed") +
    geom_hline(yintercept = 0.67, color = "#98a2b3", linetype = "dashed") +
    geom_point(alpha = 0.78) +
    ggrepel::geom_text_repel(
      data = label_dt,
      aes(label = transport_label),
      size = 2.6,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0,
      max.overlaps = 30,
      show.legend = FALSE
    ) +
    coord_cartesian(xlim = c(-0.04, 1.04), ylim = c(-0.04, 1.04),
                    clip = "off") +
    labs(
      title = "RDG Consensus Transportability Under Grouped Omission",
      subtitle = "Only grouped omissions with >=2 omitted and >=2 remaining studies are evaluated",
      x = "Minimum branch-direction retention",
      y = "Minimum thresholded-support retention",
      color = "Transport class",
      size = "Agreement score"
    ) +
    guides(
      size = "none",
      color = guide_legend(nrow = 2, byrow = TRUE)
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
  ggsave(file.path(figure_dir, "model_agreement_transport_consensus.png"),
         p_consensus, width = 11, height = 7.4, dpi = 180)
  ggsave(file.path(figure_dir, "model_agreement_transport_consensus.pdf"),
         p_consensus, width = 11, height = 7.4)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_consensus,
      file.path(figure_dir, "model_agreement_transport_consensus.html"),
      title = "RDG consensus grouped-omission transportability"
    )
  }
}

message("Saved grouped-omission transportability outputs in: ", output_dir)
print(summary_metrics)
