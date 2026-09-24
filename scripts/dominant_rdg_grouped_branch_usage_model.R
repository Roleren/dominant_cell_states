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
output_dir <- file.path(analysis_dir, "dominant_rdg_grouped_branch_usage")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)

condition_helper <- file.path(analysis_dir, "scripts",
                              "dominant_condition_families.R")
if (file.exists(condition_helper)) source(condition_helper)

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
feature_file <- file.path(analysis_dir, "dominant_uorf_feature_expression.csv")

message("Merged-replicate RDG branch-usage model")
message("  1. Build Run x gene x RDG feature adapter from FST-derived feature counts")
message("  2. Collapse runs into same-condition replicate groups")
message("  3. Compare grouped case/control branch allocation with count gates")
message("  4. Summarize gene x design-family branch shifts")

if (!file.exists(feature_file)) {
  stop("Missing feature-expression input: ", feature_file, call. = FALSE)
}
if (!file.exists(metadata_file)) {
  stop("Missing metadata input: ", metadata_file, call. = FALSE)
}

stage_message <- function(label) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", label)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

clean_level <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "NA", "N/A", "na", "n/a", "NULL", "null",
               "None", "none")] <- NA_character_
  x
}

collapse_examples <- function(x, n = 6L) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  if (!length(x)) return("")
  paste(head(x, n), collapse = ";")
}

first_nonempty <- function(x, default = "") {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) x[[1]] else default
}

median_safe <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else median(x)
}

max_safe <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else max(x)
}

weighted_mean_safe <- function(x, w = NULL) {
  x <- safe_num(x)
  if (is.null(w)) w <- rep(1, length(x))
  w <- safe_num(w)
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

logit <- function(p) log(p / (1 - p))

branch_feature_class <- function(feature_type, category) {
  fifelse(
    feature_type == "clean_CDS" | category == "clean_CDS",
    "clean_CDS",
    fifelse(
      feature_type == "overlapping_uorf" |
        category %chin% c("uoORF", "ouORF"),
      "overlapping_uORF",
      fifelse(
        feature_type == "leader_uorf" | category == "uORF",
        "leader_uORF",
        fifelse(
          grepl("internal|iorf", feature_type, ignore.case = TRUE) |
            category == "internal",
          "iORF",
          fifelse(
            grepl("nte|ntt", feature_type, ignore.case = TRUE) |
              category %chin% c("NTE", "NTT"),
            "NTE_NTT",
            "other_feature"
          )
        )
      )
    )
  )
}

random_effects_summary <- function(effect, se, weight = NULL) {
  effect <- safe_num(effect)
  se <- safe_num(se)
  if (is.null(weight)) weight <- rep(1, length(effect))
  weight <- safe_num(weight)
  keep <- is.finite(effect) & is.finite(se) & se > 0 &
    is.finite(weight) & weight > 0
  if (!any(keep)) {
    return(data.table(
      n = 0L,
      effect = NA_real_,
      se = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      z = NA_real_,
      p = NA_real_,
      tau2 = NA_real_,
      i2 = NA_real_,
      direction_fraction = NA_real_,
      ci_excludes_zero = FALSE
    ))
  }

  yi <- effect[keep]
  vi <- se[keep]^2
  external_weight <- weight[keep] / mean(weight[keep])
  wi <- external_weight / vi
  fixed <- sum(wi * yi) / sum(wi)
  q <- sum(wi * (yi - fixed)^2)
  df <- length(yi) - 1L
  c_val <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (df > 0 && c_val > 0) max(0, (q - df) / c_val) else 0
  wi_re <- external_weight / (vi + tau2)
  mu <- sum(wi_re * yi) / sum(wi_re)
  se_mu <- sqrt(1 / sum(wi_re))
  lower <- mu - 1.96 * se_mu
  upper <- mu + 1.96 * se_mu
  z <- mu / se_mu
  p <- 2 * pnorm(-abs(z))
  sign_ref <- sign(mu)
  direction_fraction <- if (sign_ref == 0) mean(yi == 0) else mean(sign(yi) == sign_ref)
  i2 <- if (q > 0 && df > 0) max(0, (q - df) / q) else 0
  data.table(
    n = length(yi),
    effect = mu,
    se = se_mu,
    lower = lower,
    upper = upper,
    z = z,
    p = p,
    tau2 = tau2,
    i2 = i2,
    direction_fraction = direction_fraction,
    ci_excludes_zero = isTRUE(lower > 0 | upper < 0)
  )
}

read_metadata <- function() {
  fields <- c(
    "Run", "study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE", "CONDITION",
    "INHIBITOR", "FRACTION", "Cancer_type", "Sex", "TIMEPOINT"
  )
  m <- fread(metadata_file, showProgress = FALSE)
  keep <- intersect(fields, names(m))
  m <- m[, ..keep]
  for (field in setdiff(fields, names(m))) m[, (field) := NA_character_]
  for (field in setdiff(fields, "Run")) {
    m[, (field) := clean_level(get(field))]
  }
  for (field in setdiff(fields, c("Run", "CONDITION"))) {
    m[is.na(get(field)) | !nzchar(as.character(get(field))),
      (field) := "MISSING"]
  }
  m[, CONDITION := clean_level(CONDITION)]
  m[is.na(CONDITION) | !nzchar(CONDITION), CONDITION := "MISSING"]
  if (exists("dominant_is_control_condition", mode = "function")) {
    m[, is_control_condition := dominant_is_control_condition(CONDITION)]
  } else {
    m[, is_control_condition :=
        tolower(CONDITION) %chin% c("wt", "control", "ctrl", "mock",
                                    "vehicle", "dmso", "untreated")]
  }
  if (exists("dominant_classify_design_family", mode = "function")) {
    m[, design_family := dominant_classify_design_family(
      condition = CONDITION,
      fraction = FRACTION,
      inhibitor = INHIBITOR,
      cell_line = CELL_LINE,
      tissue = TISSUE,
      gene = GENE,
      cancer_type = Cancer_type
    )]
    m[, perturbation_class :=
        dominant_classify_perturbation_class(CONDITION, design_family)]
  } else {
    m[, design_family := fifelse(is_control_condition, "control",
                                 "other_non_control")]
    m[, perturbation_class := fifelse(is_control_condition, "control",
                                      "other_non_control")]
  }
  unique(m, by = "Run")
}

metadata <- read_metadata()
stage_message("Loaded metadata")

feature_cols <- c(
  "gene_symbol", "tx_id", "feature_id", "feature_type", "category",
  "translon_sources_merged", "feature_bases", "cds_overlap_bases",
  "measured_bases", "measurement_scope", "Run", "raw_counts", "fpkm_like"
)
feature_cols <- intersect(feature_cols, names(fread(feature_file, nrows = 0)))
fx <- fread(feature_file, select = feature_cols, showProgress = FALSE)
stage_message("Loaded feature-expression table")
for (column in intersect(c("raw_counts", "fpkm_like", "feature_bases",
                           "cds_overlap_bases", "measured_bases"),
                         names(fx))) {
  fx[, (column) := safe_num(get(column))]
}
fx[, branch_feature_class := branch_feature_class(feature_type, category)]
fx <- fx[branch_feature_class %chin% c("clean_CDS", "leader_uORF",
                                       "overlapping_uORF")]

setkey(metadata, Run)
setkey(fx, Run)
fx <- metadata[fx]
setkey(fx, NULL)
stage_message("Joined metadata to allocation features")

context_fields <- c(
  "study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE", "CONDITION",
  "INHIBITOR", "FRACTION", "Cancer_type", "Sex", "TIMEPOINT"
)
for (field in context_fields) {
  fx[is.na(get(field)) | !nzchar(as.character(get(field))),
     (field) := "MISSING"]
}

fx[, merged_replicate_group_id := do.call(
  paste,
  c(.SD, sep = " | ")
), .SDcols = context_fields]

feature_id_cols <- c(
  "gene_symbol", "tx_id", "feature_id", "feature_type", "category",
  "branch_feature_class", "translon_sources_merged", "feature_bases",
  "cds_overlap_bases", "measured_bases", "measurement_scope"
)

replicate_group_runs <- unique(
  fx[, .(merged_replicate_group_id, Run)]
)[, .(
  n_runs_merged = .N,
  example_runs = collapse_examples(Run, 8L)
), by = merged_replicate_group_id]

group_features <- fx[, .(
  raw_counts = sum(raw_counts, na.rm = TRUE),
  fpkm_like_sum = sum(fpkm_like, na.rm = TRUE),
  fpkm_like_mean = mean(fpkm_like, na.rm = TRUE),
  study = study[1],
  AUTHOR = AUTHOR[1],
  CELL_LINE = CELL_LINE[1],
  TISSUE = TISSUE[1],
  GENE = GENE[1],
  CONDITION = CONDITION[1],
  INHIBITOR = INHIBITOR[1],
  FRACTION = FRACTION[1],
  Cancer_type = Cancer_type[1],
  Sex = Sex[1],
  TIMEPOINT = TIMEPOINT[1],
  is_control_condition = is_control_condition[1],
  design_family = design_family[1],
  perturbation_class = perturbation_class[1]
), by = c("merged_replicate_group_id", feature_id_cols)]

group_features <- merge(
  group_features,
  replicate_group_runs,
  by = "merged_replicate_group_id",
  all.x = TRUE,
  sort = FALSE
)
setnames(group_features, "merged_replicate_group_id", "group_id")

group_totals <- group_features[, .(
  branch_total_raw_counts = sum(raw_counts, na.rm = TRUE),
  clean_cds_raw_counts = sum(raw_counts[branch_feature_class == "clean_CDS"],
                             na.rm = TRUE),
  n_branch_features = .N,
  n_uorf_features = sum(branch_feature_class %chin% c("leader_uORF",
                                                      "overlapping_uORF"))
), by = .(group_id, gene_symbol, tx_id)]
group_features <- merge(
  group_features,
  group_totals,
  by = c("group_id", "gene_symbol", "tx_id"),
  all.x = TRUE,
  sort = FALSE
)
group_features[, relative_use_counts := fifelse(
  branch_total_raw_counts > 0,
  raw_counts / branch_total_raw_counts,
  NA_real_
)]
group_features[, replicate_gate_n2 := n_runs_merged >= 2]
group_features[, feature_count_ge_30 := raw_counts >= 30]
group_features[, branch_total_ge_100 := branch_total_raw_counts >= 100]
group_features[, feature_key := paste(gene_symbol, tx_id, feature_id,
                                      sep = " | ")]
feature_meta <- unique(group_features[, c("feature_key", feature_id_cols),
                                      with = FALSE])
stage_message("Built merged-replicate feature adapter")

adapter_file <- file.path(output_dir,
                          "dominant_rdg_grouped_branch_usage_adapter.csv")
fwrite(group_features, adapter_file)

make_stratum_id <- function(dt, fields) {
  fields <- intersect(fields, names(dt))
  do.call(paste, c(dt[, ..fields], sep = " | "))
}

build_contrasts_for_scope <- function(group_features, match_scope,
                                      stratum_fields) {
  dt <- copy(group_features)
  dt[, stratum_id := make_stratum_id(dt, stratum_fields)]
  dt <- dt[replicate_gate_n2 == TRUE]
  if (nrow(dt) == 0) return(data.table())

  case_feature <- dt[
    is_control_condition != TRUE,
    .(
      n_case_groups = uniqueN(group_id),
      n_case_runs_merged = sum(n_runs_merged, na.rm = TRUE),
      case_branch_total_counts = sum(branch_total_raw_counts, na.rm = TRUE),
      case_sum_raw_counts = sum(raw_counts, na.rm = TRUE),
      study = study[1],
      AUTHOR = AUTHOR[1],
      CELL_LINE = CELL_LINE[1],
      TISSUE = TISSUE[1],
      GENE = GENE[1],
      INHIBITOR = INHIBITOR[1],
      FRACTION = FRACTION[1],
      Cancer_type = Cancer_type[1],
      Sex = Sex[1],
      TIMEPOINT = TIMEPOINT[1]
    ),
    by = c(
      "stratum_id",
      "case_condition" = "CONDITION",
      "perturbation_class",
      "design_family",
      "feature_key"
    )
  ]
  control_feature <- dt[
    is_control_condition == TRUE,
    .(
      n_control_groups = uniqueN(group_id),
      n_control_runs_merged = sum(n_runs_merged, na.rm = TRUE),
      control_branch_total_counts =
        sum(branch_total_raw_counts, na.rm = TRUE),
      control_sum_raw_counts = sum(raw_counts, na.rm = TRUE),
      control_conditions = paste(sort(unique(CONDITION)), collapse = ";")
    ),
    by = c("stratum_id", "feature_key")
  ]

  if (nrow(case_feature) == 0 || nrow(control_feature) == 0) {
    return(data.table())
  }

  join_cols <- c("stratum_id", "feature_key")
  out <- merge(
    case_feature,
    control_feature,
    by = join_cols,
    all = FALSE,
    sort = FALSE
  )
  if (nrow(out) == 0) return(data.table())

  alpha <- 0.5
  out[, case_other_branch_counts :=
        pmax(case_branch_total_counts - case_sum_raw_counts, 0)]
  out[, control_other_branch_counts :=
        pmax(control_branch_total_counts - control_sum_raw_counts, 0)]
  out[, case_branch_usage_proportion :=
        (case_sum_raw_counts + alpha) /
        (case_branch_total_counts + 2 * alpha)]
  out[, control_branch_usage_proportion :=
        (control_sum_raw_counts + alpha) /
        (control_branch_total_counts + 2 * alpha)]
  out[, grouped_usage_delta :=
        case_branch_usage_proportion - control_branch_usage_proportion]
  out[, grouped_branch_usage_log_or :=
        logit(case_branch_usage_proportion) -
        logit(control_branch_usage_proportion)]
  out[, grouped_branch_usage_log_or_se := sqrt(
    1 / (case_sum_raw_counts + alpha) +
      1 / (case_other_branch_counts + alpha) +
      1 / (control_sum_raw_counts + alpha) +
      1 / (control_other_branch_counts + alpha)
  )]
  out[, grouped_branch_usage_z :=
        grouped_branch_usage_log_or / grouped_branch_usage_log_or_se]
  out[, grouped_branch_usage_p :=
        2 * pnorm(-abs(grouped_branch_usage_z))]
  out[, `:=`(
    match_scope = match_scope,
    sample_set_signature = paste(
      match_scope,
      stratum_id,
      case_condition,
      sep = " || "
    )
  )]
  out <- merge(out, feature_meta, by = "feature_key", all.x = TRUE,
               sort = FALSE)
  out[, feature_key := NULL]
  out
}

exact_fields <- c(
  "study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE", "INHIBITOR",
  "FRACTION", "Cancer_type", "Sex", "TIMEPOINT"
)
relaxed_fields <- c("study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE")

stage_message("Building exact-context grouped contrasts")
exact_contrasts <- build_contrasts_for_scope(
  group_features,
  "exact_context",
  exact_fields
)
stage_message("Building relaxed study/cell/gene grouped contrasts")
relaxed_contrasts <- build_contrasts_for_scope(
  group_features,
  "study_cell_gene_context",
  relaxed_fields
)
if (nrow(exact_contrasts) > 0 && nrow(relaxed_contrasts) > 0) {
  fallback_key_cols <- c(
    "study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE",
    "case_condition", "design_family", "gene_symbol", "tx_id", "feature_id"
  )
  exact_contrasts[, fallback_key := do.call(
    paste,
    c(.SD, sep = " | ")
  ), .SDcols = fallback_key_cols]
  relaxed_contrasts[, fallback_key := do.call(
    paste,
    c(.SD, sep = " | ")
  ), .SDcols = fallback_key_cols]
  relaxed_contrasts <- relaxed_contrasts[
    !fallback_key %chin% exact_contrasts$fallback_key
  ]
  exact_contrasts[, fallback_key := NULL]
  relaxed_contrasts[, fallback_key := NULL]
}
contrast_dt <- rbindlist(list(exact_contrasts, relaxed_contrasts), fill = TRUE)
stage_message("Built grouped case/control branch contrasts")

if (nrow(contrast_dt) == 0) {
  warning("No grouped case/control branch contrasts could be built.")
  fwrite(data.table(), file.path(output_dir,
                                 "dominant_rdg_grouped_branch_usage_contrast_table.csv"))
  fwrite(data.table(), file.path(output_dir,
                                 "dominant_rdg_grouped_branch_usage_feature_summary.csv"))
  fwrite(data.table(), file.path(output_dir,
                                 "dominant_rdg_grouped_branch_usage_gene_design_summary.csv"))
  fwrite(data.table(), file.path(output_dir,
                                 "dominant_rdg_grouped_branch_usage_review.csv"))
  fwrite(data.table(metric = "contrast_rows", value = 0),
         file.path(output_dir,
                   "dominant_rdg_grouped_branch_usage_summary_metrics.csv"))
  quit(save = "no", status = 0)
}

contrast_dt[, grouped_branch_usage_q :=
              p.adjust(grouped_branch_usage_p, method = "BH")]
contrast_dt[, `:=`(
  grouped_min_runs_merged = pmin(n_case_runs_merged, n_control_runs_merged),
  grouped_min_branch_total_counts =
    pmin(case_branch_total_counts, control_branch_total_counts),
  grouped_max_feature_counts =
    pmax(case_sum_raw_counts, control_sum_raw_counts),
  grouped_min_feature_counts =
    pmin(case_sum_raw_counts, control_sum_raw_counts)
)]
contrast_dt[, `:=`(
  grouped_passes_replicate_gate =
    n_case_runs_merged >= 2 & n_control_runs_merged >= 2,
  grouped_passes_total_count_gate =
    grouped_min_branch_total_counts >= 100,
  grouped_passes_feature_count_gate_any =
    grouped_max_feature_counts >= 30,
  grouped_passes_feature_count_gate_both =
    grouped_min_feature_counts >= 30
)]
contrast_dt[, grouped_branch_usage_evaluable :=
              grouped_passes_replicate_gate &
              grouped_passes_total_count_gate &
              grouped_passes_feature_count_gate_any &
              is.finite(grouped_branch_usage_log_or) &
              is.finite(grouped_branch_usage_log_or_se)]
contrast_dt[, grouped_branch_usage_strict :=
              grouped_branch_usage_evaluable &
              grouped_passes_feature_count_gate_both &
              match_scope == "exact_context"]
contrast_dt[, grouped_branch_usage_relaxed :=
              grouped_branch_usage_evaluable &
              grouped_passes_feature_count_gate_both &
              match_scope != "exact_context"]
contrast_dt[, grouped_row_weight := pmax(
  1,
  log1p(grouped_max_feature_counts) *
    sqrt(pmax(grouped_min_branch_total_counts, 1)) *
    fifelse(match_scope == "exact_context", 1, 0.55)
)]
contrast_dt[!is.finite(grouped_row_weight), grouped_row_weight := 1]

contrast_file <- file.path(output_dir,
                           "dominant_rdg_grouped_branch_usage_contrast_table.csv")
fwrite(contrast_dt, contrast_file)

feature_summary <- contrast_dt[, {
  re_strict <- random_effects_summary(
    grouped_branch_usage_log_or[grouped_branch_usage_strict == TRUE],
    grouped_branch_usage_log_or_se[grouped_branch_usage_strict == TRUE],
    grouped_row_weight[grouped_branch_usage_strict == TRUE]
  )
  re_all <- random_effects_summary(
    grouped_branch_usage_log_or[grouped_branch_usage_evaluable == TRUE],
    grouped_branch_usage_log_or_se[grouped_branch_usage_evaluable == TRUE],
    grouped_row_weight[grouped_branch_usage_evaluable == TRUE]
  )
  use_re <- if (re_strict$n[[1]] >= 2L) re_strict else re_all
  data.table(
    branch_feature_class = first_nonempty(branch_feature_class),
    feature_type = first_nonempty(feature_type),
    category = first_nonempty(category),
    n_grouped_contrast_rows = .N,
    n_grouped_evaluable_rows = sum(grouped_branch_usage_evaluable,
                                   na.rm = TRUE),
    n_grouped_strict_rows = sum(grouped_branch_usage_strict, na.rm = TRUE),
    n_grouped_relaxed_rows = sum(grouped_branch_usage_relaxed, na.rm = TRUE),
    n_grouped_strict_studies =
      uniqueN(study[grouped_branch_usage_strict == TRUE]),
    n_grouped_strict_cell_lines =
      uniqueN(CELL_LINE[grouped_branch_usage_strict == TRUE]),
    n_grouped_evaluable_studies =
      uniqueN(study[grouped_branch_usage_evaluable == TRUE]),
    n_grouped_case_conditions = uniqueN(case_condition),
    total_case_branch_counts =
      sum(case_branch_total_counts[grouped_branch_usage_evaluable == TRUE],
          na.rm = TRUE),
    total_control_branch_counts =
      sum(control_branch_total_counts[grouped_branch_usage_evaluable == TRUE],
          na.rm = TRUE),
    max_case_feature_count = max_safe(case_sum_raw_counts),
    max_control_feature_count = max_safe(control_sum_raw_counts),
    grouped_weighted_usage_delta = weighted_mean_safe(
      grouped_usage_delta[grouped_branch_usage_evaluable == TRUE],
      grouped_row_weight[grouped_branch_usage_evaluable == TRUE]
    ),
    grouped_median_usage_delta =
      median_safe(grouped_usage_delta[grouped_branch_usage_evaluable == TRUE]),
    grouped_max_abs_usage_delta =
      max_safe(abs(grouped_usage_delta[grouped_branch_usage_evaluable == TRUE])),
    grouped_fraction_positive_usage_delta = mean(
      grouped_usage_delta[grouped_branch_usage_evaluable == TRUE] > 0,
      na.rm = TRUE
    ),
    grouped_branch_usage_log_or = use_re$effect,
    grouped_branch_usage_log_or_se = use_re$se,
    grouped_branch_usage_log_or_lower = use_re$lower,
    grouped_branch_usage_log_or_upper = use_re$upper,
    grouped_branch_usage_log_or_p = use_re$p,
    grouped_branch_usage_log_or_i2 = use_re$i2,
    grouped_branch_usage_direction_fraction = use_re$direction_fraction,
    grouped_branch_usage_interval_supported = use_re$ci_excludes_zero,
    grouped_branch_usage_interval_source =
      if (re_strict$n[[1]] >= 2L) "strict_exact_random_effects" else "all_evaluable_random_effects",
    strongest_grouped_case_condition =
      case_condition[which.max(abs(grouped_usage_delta))][1],
    strongest_grouped_control_conditions =
      control_conditions[which.max(abs(grouped_usage_delta))][1],
    strongest_grouped_study = study[which.max(abs(grouped_usage_delta))][1],
    nearest_grouped_contexts = collapse_examples(sample_set_signature, 5L)
  )
}, by = .(gene_symbol, tx_id, design_family, feature_id)]

feature_summary[, grouped_branch_usage_log_or_q :=
                  p.adjust(grouped_branch_usage_log_or_p, method = "BH")]
feature_summary[, grouped_branch_usage_score :=
                  0.22 * pmin(n_grouped_strict_rows / 3, 1) +
                  0.16 * pmin(n_grouped_strict_studies / 2, 1) +
                  0.16 * pmin(n_grouped_evaluable_rows / 4, 1) +
                  0.18 * pmin(abs(grouped_weighted_usage_delta) / 0.20, 1) +
                  0.16 * fifelse(is.finite(grouped_branch_usage_log_or_q),
                                  pmin(-log10(pmax(
                                    grouped_branch_usage_log_or_q,
                                    1e-300
                                  )) / 4, 1),
                                  0) +
                  0.08 * as.numeric(
                    grouped_branch_usage_interval_supported == TRUE
                  ) +
                  0.04 * pmin(log1p(pmax(max_case_feature_count,
                                         max_control_feature_count,
                                         na.rm = TRUE)) / log1p(1000), 1)]
feature_summary[!is.finite(grouped_branch_usage_score),
                grouped_branch_usage_score := 0]

feature_summary[, grouped_branch_usage_support_class := fifelse(
  n_grouped_strict_rows >= 2 &
    n_grouped_strict_studies >= 2 &
    grouped_branch_usage_interval_supported == TRUE &
    is.finite(grouped_branch_usage_log_or_q) &
    grouped_branch_usage_log_or_q <= 0.10 &
    abs(grouped_weighted_usage_delta) >= 0.05 &
    grouped_branch_usage_direction_fraction >= 0.70,
  "replicated_exact_group_branch_shift",
  fifelse(
    n_grouped_strict_rows >= 1 &
      is.finite(grouped_branch_usage_log_or_p) &
      grouped_branch_usage_log_or_p <= 0.01 &
      abs(grouped_weighted_usage_delta) >= 0.10,
    "single_exact_group_branch_shift",
    fifelse(
      n_grouped_evaluable_rows >= 2 &
        n_grouped_relaxed_rows >= 2 &
        grouped_branch_usage_interval_supported == TRUE &
        is.finite(grouped_branch_usage_log_or_q) &
        grouped_branch_usage_log_or_q <= 0.10 &
        abs(grouped_weighted_usage_delta) >= 0.05,
      "relaxed_group_branch_shift_review",
      fifelse(
        n_grouped_evaluable_rows > 0,
        "group_evaluable_no_strong_shift",
        "low_count_or_no_grouped_control"
      )
    )
  )
)]

setorder(feature_summary, -grouped_branch_usage_score,
         grouped_branch_usage_log_or_q, gene_symbol, design_family, feature_id)
feature_file_out <- file.path(output_dir,
                              "dominant_rdg_grouped_branch_usage_feature_summary.csv")
fwrite(feature_summary, feature_file_out)

top_feature <- feature_summary[, .SD[1], by = .(gene_symbol, tx_id, design_family)]
class_wide <- dcast(
  feature_summary[
    branch_feature_class %chin% c("clean_CDS", "leader_uORF",
                                  "overlapping_uORF")
  ],
  gene_symbol + tx_id + design_family ~ branch_feature_class,
  value.var = c("feature_id", "grouped_weighted_usage_delta",
                "grouped_branch_usage_log_or",
                "grouped_branch_usage_log_or_q",
                "grouped_branch_usage_score"),
  fun.aggregate = function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) {
      if (is.numeric(x)) return(NA_real_)
      return(NA_character_)
    }
    if (is.numeric(x)) x[which.max(abs(x))] else as.character(x[[1]])
  },
  fill = NA
)

gene_design_summary <- merge(
  top_feature[, .(
    gene_symbol, tx_id, design_family,
    grouped_branch_usage_top_feature_id = feature_id,
    grouped_branch_usage_top_feature_class = branch_feature_class,
    grouped_branch_usage_top_feature_type = feature_type,
    grouped_branch_usage_top_support_class =
      grouped_branch_usage_support_class,
    grouped_branch_usage_top_score = grouped_branch_usage_score,
    grouped_branch_usage_top_weighted_usage_delta =
      grouped_weighted_usage_delta,
    grouped_branch_usage_top_log_or = grouped_branch_usage_log_or,
    grouped_branch_usage_top_q = grouped_branch_usage_log_or_q,
    grouped_branch_usage_top_i2 = grouped_branch_usage_log_or_i2,
    grouped_branch_usage_top_direction_fraction =
      grouped_branch_usage_direction_fraction,
    grouped_branch_usage_top_strict_rows = n_grouped_strict_rows,
    grouped_branch_usage_top_strict_studies = n_grouped_strict_studies,
    grouped_branch_usage_top_evaluable_rows = n_grouped_evaluable_rows,
    grouped_branch_usage_top_context = strongest_grouped_case_condition,
    grouped_branch_usage_top_control = strongest_grouped_control_conditions,
    grouped_branch_usage_top_study = strongest_grouped_study
  )],
  class_wide,
  by = c("gene_symbol", "tx_id", "design_family"),
  all.x = TRUE,
  sort = FALSE
)
setorder(gene_design_summary, -grouped_branch_usage_top_score,
         grouped_branch_usage_top_q, gene_symbol, design_family)
gene_design_file <- file.path(output_dir,
                              "dominant_rdg_grouped_branch_usage_gene_design_summary.csv")
fwrite(gene_design_summary, gene_design_file)

review <- gene_design_summary[
  grouped_branch_usage_top_support_class != "low_count_or_no_grouped_control" &
    (grouped_branch_usage_top_support_class !=
       "group_evaluable_no_strong_shift" |
       grouped_branch_usage_top_score >= 0.55)
]
setorder(review, -grouped_branch_usage_top_score,
         grouped_branch_usage_top_q, gene_symbol, design_family)
review_file <- file.path(output_dir,
                         "dominant_rdg_grouped_branch_usage_review.csv")
fwrite(review, review_file)

summary_metrics <- data.table(
  metric = c(
    "adapter_rows",
    "adapter_groups",
    "adapter_groups_n2plus",
    "contrast_rows",
    "exact_contrast_rows",
    "relaxed_contrast_rows",
    "feature_summary_rows",
    "gene_design_rows",
    "review_rows",
    "replicated_exact_group_branch_shift_rows",
    "single_exact_group_branch_shift_rows",
    "relaxed_group_branch_shift_review_rows",
    "viral_infection_review_rows",
    "genes",
    "design_families"
  ),
  value = c(
    nrow(group_features),
    uniqueN(group_features$group_id),
    uniqueN(group_features[replicate_gate_n2 == TRUE, group_id]),
    nrow(contrast_dt),
    nrow(contrast_dt[match_scope == "exact_context"]),
    nrow(contrast_dt[match_scope == "study_cell_gene_context"]),
    nrow(feature_summary),
    nrow(gene_design_summary),
    nrow(review),
    sum(feature_summary$grouped_branch_usage_support_class ==
          "replicated_exact_group_branch_shift", na.rm = TRUE),
    sum(feature_summary$grouped_branch_usage_support_class ==
          "single_exact_group_branch_shift", na.rm = TRUE),
    sum(feature_summary$grouped_branch_usage_support_class ==
          "relaxed_group_branch_shift_review", na.rm = TRUE),
    nrow(review[design_family == "viral_infection"]),
    uniqueN(contrast_dt$gene_symbol),
    uniqueN(contrast_dt$design_family)
  )
)
metrics_file <- file.path(output_dir,
                          "dominant_rdg_grouped_branch_usage_summary_metrics.csv")
fwrite(summary_metrics, metrics_file)

plot_dt <- feature_summary[
  grouped_branch_usage_support_class != "low_count_or_no_grouped_control" &
    grouped_branch_usage_score >= 0.45
][order(-grouped_branch_usage_score)][seq_len(min(.N, 45))]
if (nrow(plot_dt) > 0) {
  plot_dt[, plot_label := paste0(gene_symbol, " / ", design_family,
                                 " / ", feature_id)]
  plot_dt[, plot_label := factor(plot_label, levels = rev(plot_label))]
  p_top <- ggplot(
    plot_dt,
    aes(
      x = grouped_weighted_usage_delta,
      y = plot_label,
      color = branch_feature_class,
      shape = grouped_branch_usage_support_class,
      text = paste0(
        "Gene: ", gene_symbol,
        "<br>Design family: ", design_family,
        "<br>Feature: ", feature_id,
        "<br>Class: ", branch_feature_class,
        "<br>Usage delta: ", round(grouped_weighted_usage_delta, 3),
        "<br>logOR: ", round(grouped_branch_usage_log_or, 3),
        "<br>q: ", signif(grouped_branch_usage_log_or_q, 3),
        "<br>Score: ", round(grouped_branch_usage_score, 3),
        "<br>Strict rows: ", n_grouped_strict_rows,
        "<br>Strict studies: ", n_grouped_strict_studies
      )
    )
  ) +
    geom_vline(xintercept = 0, color = "#7b8494", linewidth = 0.35) +
    geom_point(size = 2.6, alpha = 0.9) +
    scale_color_manual(values = c(
      clean_CDS = "#1f5aa6",
      leader_uORF = "#c57b00",
      overlapping_uORF = "#9b1c31"
    )) +
    labs(
      title = "Merged-replicate RDG branch-usage shifts",
      subtitle = "Case/control feature-vs-other branch allocation; only clean-CDS, leader-uORF, and overlapping-uORF are in the denominator",
      x = "Weighted case-control branch usage delta",
      y = "",
      color = "Feature class",
      shape = "Support"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave(file.path(figure_dir,
                   "dominant_rdg_grouped_branch_usage_top_shifts.png"),
         p_top, width = 11.2, height = 8.8, dpi = 180)
  ggsave(file.path(figure_dir,
                   "dominant_rdg_grouped_branch_usage_top_shifts.pdf"),
         p_top, width = 11.2, height = 8.8)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p_top,
      file.path(figure_dir,
                "dominant_rdg_grouped_branch_usage_top_shifts.html"),
      title = "Merged-replicate RDG branch-usage shifts"
    )
  }
}

viral_dt <- feature_summary[
  design_family == "viral_infection" &
    n_grouped_evaluable_rows > 0
][order(-grouped_branch_usage_score)][seq_len(min(.N, 45))]
if (nrow(viral_dt) > 0) {
  viral_dt[, plot_label := paste0(gene_symbol, " / ", feature_id)]
  viral_dt[, plot_label := factor(plot_label, levels = rev(plot_label))]
  p_viral <- ggplot(
    viral_dt,
    aes(
      x = grouped_weighted_usage_delta,
      y = plot_label,
      color = branch_feature_class,
      size = pmin(grouped_branch_usage_score, 1),
      text = paste0(
        "Gene: ", gene_symbol,
        "<br>Feature: ", feature_id,
        "<br>Class: ", branch_feature_class,
        "<br>Usage delta: ", round(grouped_weighted_usage_delta, 3),
        "<br>q: ", signif(grouped_branch_usage_log_or_q, 3),
        "<br>Support: ", grouped_branch_usage_support_class,
        "<br>Top context: ", strongest_grouped_case_condition,
        "<br>Study: ", strongest_grouped_study
      )
    )
  ) +
    geom_vline(xintercept = 0, color = "#7b8494", linewidth = 0.35) +
    geom_point(alpha = 0.9) +
    scale_size(range = c(1.8, 5), guide = "none") +
    scale_color_manual(values = c(
      clean_CDS = "#1f5aa6",
      leader_uORF = "#c57b00",
      overlapping_uORF = "#9b1c31"
    )) +
    labs(
      title = "Viral-infection merged-replicate branch-usage review",
      subtitle = "Positive x means higher feature share in infected/case groups",
      x = "Weighted case-control branch usage delta",
      y = "",
      color = "Feature class"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave(file.path(figure_dir,
                   "dominant_rdg_grouped_branch_usage_viral_review.png"),
         p_viral, width = 10.6, height = 8.2, dpi = 180)
  ggsave(file.path(figure_dir,
                   "dominant_rdg_grouped_branch_usage_viral_review.pdf"),
         p_viral, width = 10.6, height = 8.2)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p_viral,
      file.path(figure_dir,
                "dominant_rdg_grouped_branch_usage_viral_review.html"),
      title = "Viral-infection merged-replicate branch-usage review"
    )
  }
}

message("Saved grouped branch adapter: ", adapter_file)
message("Saved grouped branch contrast table: ", contrast_file)
message("Saved grouped branch feature summary: ", feature_file_out)
message("Saved grouped branch gene-design summary: ", gene_design_file)
message("Saved grouped branch review: ", review_file)
message("Saved grouped branch metrics: ", metrics_file)
