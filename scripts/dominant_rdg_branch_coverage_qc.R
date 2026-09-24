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
    if (dir.exists(candidate) && dir.exists(file.path(candidate, "scripts"))) {
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
output_dir <- file.path(analysis_dir, "dominant_rdg_branch_coverage_qc")
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
state_file <- file.path(analysis_dir, "human_dominant_cell_states_clean_cds.csv")

message("RDG branch coverage QC:")
message("  1. Compare single-run feature coverage to merged same-condition groups")
message("  2. Quantify small-uORF count gates before sample-resolution modeling")
message("  3. Recommend whether branch modeling should use single runs or merged replicates")

stage_message <- function(label) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", label)
}

if (!file.exists(feature_file)) {
  stop("Missing feature-expression input: ", feature_file, call. = FALSE)
}
if (!file.exists(metadata_file)) {
  stop("Missing metadata input: ", metadata_file, call. = FALSE)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

clean_level <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "NA", "N/A", "na", "n/a", "NULL", "null",
               "None", "none")] <- NA_character_
  x
}

median_safe <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else median(x)
}

quantile_safe <- function(x, prob = 0.95) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE,
                             type = 8))
}

collapse_examples <- function(x, n = 5L) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  if (!length(x)) return("")
  paste(head(x, n), collapse = ";")
}

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

feature_length_bin <- function(width) {
  width <- safe_num(width)
  cut(
    width,
    breaks = c(-Inf, 15, 30, 60, 120, Inf),
    labels = c("<=15 nt", "16-30 nt", "31-60 nt", "61-120 nt", ">120 nt"),
    right = TRUE
  )
}

read_metadata <- function() {
  fields <- c(
    "Run", "study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE", "CONDITION",
    "INHIBITOR", "FRACTION", "Cancer_type", "Sex", "TIMEPOINT",
    "Frame_usage_0", "Frame_usage_1", "Frame_usage_2", "Total_Cov_Signal"
  )
  m <- fread(metadata_file, showProgress = FALSE)
  keep <- intersect(fields, names(m))
  m <- m[, ..keep]
  for (field in setdiff(fields, names(m))) m[, (field) := NA_character_]
  for (field in setdiff(fields, c("Run", "Frame_usage_0", "Frame_usage_1",
                                  "Frame_usage_2", "Total_Cov_Signal"))) {
    m[, (field) := clean_level(get(field))]
  }
  for (field in c("study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE",
                  "CONDITION", "INHIBITOR", "FRACTION", "Cancer_type",
                  "Sex", "TIMEPOINT")) {
    m[is.na(get(field)) | !nzchar(get(field)), (field) := "MISSING"]
  }
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
  m
}

metadata <- read_metadata()
metadata <- unique(metadata, by = "Run")
if (file.exists(state_file)) {
  states <- fread(state_file, showProgress = FALSE)
  state_cols <- intersect(c("Run", "dominant_state", "dominant_state_score"),
                          names(states))
  if (length(state_cols) == 3L) {
    metadata <- merge(metadata, states[, ..state_cols], by = "Run",
                      all.x = TRUE, sort = FALSE)
  }
}
if (!"dominant_state" %in% names(metadata)) metadata[, dominant_state := NA_character_]
if (!"dominant_state_score" %in% names(metadata)) metadata[, dominant_state_score := NA_real_]
metadata <- unique(metadata, by = "Run")
stage_message("Loaded and deduplicated metadata")

feature_cols <- c(
  "gene_symbol", "tx_id", "feature_id", "feature_type", "category",
  "translon_source", "translon_sources_merged", "feature_bases",
  "cds_overlap_bases", "measured_bases", "measurement_scope", "Run",
  "raw_counts", "fpkm_like", "cds_reference_raw_counts",
  "cds_reference_fpkm_like"
)
feature_cols <- intersect(feature_cols, names(fread(feature_file, nrows = 0)))
fx <- fread(feature_file, select = feature_cols, showProgress = FALSE)
stage_message("Loaded feature-expression table")
for (column in c("raw_counts", "fpkm_like", "feature_bases",
                 "cds_overlap_bases", "measured_bases",
                 "cds_reference_raw_counts", "cds_reference_fpkm_like")) {
  if (column %in% names(fx)) fx[, (column) := safe_num(get(column))]
}
fx[, branch_feature_class := branch_feature_class(feature_type, category)]
fx <- fx[branch_feature_class %chin% c("clean_CDS", "leader_uORF",
                                       "overlapping_uORF", "iORF",
                                       "NTE_NTT")]
fx[, feature_length_bin := as.character(feature_length_bin(feature_bases))]
fx[is.na(feature_length_bin), feature_length_bin := "unknown"]
fx[, is_allocation_feature :=
     branch_feature_class %chin% c("clean_CDS", "leader_uORF",
                                   "overlapping_uORF")]
fx[, is_small_uorf :=
     branch_feature_class %chin% c("leader_uORF", "overlapping_uORF") &
     is.finite(feature_bases) & feature_bases <= 60]

missing_metadata_runs <- setdiff(unique(fx$Run), unique(metadata$Run))
metadata_missing_summary <- data.table(
  metric = c("feature_runs", "metadata_runs", "feature_runs_missing_metadata"),
  value = c(uniqueN(fx$Run), uniqueN(metadata$Run),
            length(missing_metadata_runs))
)
fwrite(metadata_missing_summary,
       file.path(output_dir, "branch_coverage_metadata_match_summary.csv"))

setkey(metadata, Run)
setkey(fx, Run)
fx <- metadata[fx]
setkey(fx, NULL)
stage_message("Joined metadata to feature-expression table")
for (field in c("study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE",
                "CONDITION", "INHIBITOR", "FRACTION", "Cancer_type",
                "Sex", "TIMEPOINT", "design_family", "perturbation_class",
                "dominant_state")) {
  if (!field %in% names(fx)) fx[, (field) := "MISSING"]
  fx[is.na(get(field)) | !nzchar(as.character(get(field))),
     (field) := "MISSING"]
}

replicate_fields <- c(
  "study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE", "CONDITION",
  "INHIBITOR", "FRACTION", "Cancer_type", "Sex", "TIMEPOINT"
)
fx[, merged_replicate_group_id := do.call(
  paste,
  c(.SD, sep = " | ")
), .SDcols = replicate_fields]
stage_message("Built merged replicate group ids")

feature_id_cols <- c(
  "gene_symbol", "tx_id", "feature_id", "feature_type", "category",
  "branch_feature_class", "translon_sources_merged", "feature_bases",
  "cds_overlap_bases", "measured_bases", "measurement_scope",
  "feature_length_bin", "is_allocation_feature", "is_small_uorf"
)
context_cols <- c(
  "study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE", "CONDITION",
  "INHIBITOR", "FRACTION", "Cancer_type", "Sex", "TIMEPOINT",
  "design_family", "perturbation_class", "is_control_condition",
  "dominant_state"
)

fx_model <- fx[is_allocation_feature == TRUE]

replicate_group_runs <- unique(
  fx_model[, .(merged_replicate_group_id, Run)]
)[, .(
  n_runs_merged = .N,
  example_runs = collapse_examples(Run, 6L)
), by = merged_replicate_group_id]
stage_message("Summarized replicate group run membership")

single_units <- fx_model[, .(
  unit_type = "single_run",
  unit_id = Run,
  n_runs_merged = 1L,
  example_runs = Run,
  raw_counts = raw_counts,
  fpkm_like = fpkm_like,
  cds_reference_raw_counts = cds_reference_raw_counts,
  cds_reference_fpkm_like = cds_reference_fpkm_like,
  gene_symbol = gene_symbol,
  tx_id = tx_id,
  feature_id = feature_id,
  feature_type = feature_type,
  category = category,
  branch_feature_class = branch_feature_class,
  translon_sources_merged = translon_sources_merged,
  feature_bases = feature_bases,
  cds_overlap_bases = cds_overlap_bases,
  measured_bases = measured_bases,
  measurement_scope = measurement_scope,
  feature_length_bin = feature_length_bin,
  is_allocation_feature = is_allocation_feature,
  is_small_uorf = is_small_uorf,
  study = study,
  AUTHOR = AUTHOR,
  CELL_LINE = CELL_LINE,
  TISSUE = TISSUE,
  GENE = GENE,
  CONDITION = CONDITION,
  INHIBITOR = INHIBITOR,
  FRACTION = FRACTION,
  Cancer_type = Cancer_type,
  Sex = Sex,
  TIMEPOINT = TIMEPOINT,
  design_family = design_family,
  perturbation_class = perturbation_class,
  is_control_condition = is_control_condition,
  dominant_state = dominant_state
)]
stage_message("Built single-run allocation unit table")

merged_units <- fx_model[, .(
  unit_type = "merged_condition_group",
  unit_id = merged_replicate_group_id[1],
  raw_counts = sum(raw_counts, na.rm = TRUE),
  fpkm_like = sum(fpkm_like, na.rm = TRUE),
  cds_reference_raw_counts = sum(cds_reference_raw_counts, na.rm = TRUE),
  cds_reference_fpkm_like = sum(cds_reference_fpkm_like, na.rm = TRUE),
  feature_type = feature_type[1],
  category = category[1],
  branch_feature_class = branch_feature_class[1],
  translon_sources_merged = translon_sources_merged[1],
  feature_bases = feature_bases[1],
  cds_overlap_bases = cds_overlap_bases[1],
  measured_bases = measured_bases[1],
  measurement_scope = measurement_scope[1],
  feature_length_bin = feature_length_bin[1],
  is_allocation_feature = is_allocation_feature[1],
  is_small_uorf = is_small_uorf[1],
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
  design_family = design_family[1],
  perturbation_class = perturbation_class[1],
  is_control_condition = is_control_condition[1],
  dominant_state = dominant_state[1]
), by = c("merged_replicate_group_id", "gene_symbol", "tx_id", "feature_id")]
merged_units <- merge(
  merged_units,
  replicate_group_runs,
  by = "merged_replicate_group_id",
  all.x = TRUE,
  sort = FALSE
)
merged_units[, merged_replicate_group_id := NULL]
stage_message("Built merged-condition allocation unit table")

units <- rbindlist(list(single_units, merged_units), fill = TRUE)
units[, raw_counts := safe_num(raw_counts)]
units[, fpkm_like := safe_num(fpkm_like)]
stage_message("Combined single and merged allocation unit tables")

allocation_totals <- units[
  is_allocation_feature == TRUE,
  .(
    allocation_total_raw_counts = sum(raw_counts, na.rm = TRUE),
    allocation_total_fpkm_like = sum(fpkm_like, na.rm = TRUE),
    n_allocation_features = .N,
    n_uorf_allocation_features =
      sum(branch_feature_class %chin% c("leader_uORF", "overlapping_uORF"),
          na.rm = TRUE),
    clean_cds_raw_counts = sum(raw_counts[branch_feature_class == "clean_CDS"],
                               na.rm = TRUE),
    min_uorf_raw_counts = {
      x <- raw_counts[branch_feature_class %chin% c("leader_uORF",
                                                    "overlapping_uORF")]
      x <- x[is.finite(x)]
      if (length(x)) min(x) else NA_real_
    },
    max_uorf_raw_counts = {
      x <- raw_counts[branch_feature_class %chin% c("leader_uORF",
                                                    "overlapping_uORF")]
      x <- x[is.finite(x)]
      if (length(x)) max(x) else NA_real_
    },
    n_uorf_features_ge_10 =
      sum(branch_feature_class %chin% c("leader_uORF", "overlapping_uORF") &
            raw_counts >= 10, na.rm = TRUE),
    n_uorf_features_ge_30 =
      sum(branch_feature_class %chin% c("leader_uORF", "overlapping_uORF") &
            raw_counts >= 30, na.rm = TRUE)
  ),
  by = .(unit_type, unit_id, gene_symbol, tx_id)
]
units <- merge(units, allocation_totals,
               by = c("unit_type", "unit_id", "gene_symbol", "tx_id"),
               all.x = TRUE, sort = FALSE)
stage_message("Attached gene-level allocation totals")
units[, relative_use_allocation_counts := fifelse(
  is_allocation_feature == TRUE & allocation_total_raw_counts > 0,
  raw_counts / allocation_total_raw_counts,
  NA_real_
)]
units[, relative_use_allocation_fpkm := fifelse(
  is_allocation_feature == TRUE & allocation_total_fpkm_like > 0,
  fpkm_like / allocation_total_fpkm_like,
  NA_real_
)]

units[, `:=`(
  count_ge_1 = raw_counts >= 1,
  count_ge_3 = raw_counts >= 3,
  count_ge_5 = raw_counts >= 5,
  count_ge_10 = raw_counts >= 10,
  count_ge_30 = raw_counts >= 30,
  count_ge_100 = raw_counts >= 100,
  allocation_total_ge_30 = allocation_total_raw_counts >= 30,
  allocation_total_ge_100 = allocation_total_raw_counts >= 100,
  allocation_total_ge_300 = allocation_total_raw_counts >= 300,
  clean_cds_ge_30 = clean_cds_raw_counts >= 30,
  clean_cds_ge_100 = clean_cds_raw_counts >= 100,
  branch_model_minimal_gate =
    is_allocation_feature == TRUE &
    allocation_total_raw_counts >= 100 &
    clean_cds_raw_counts >= 30 &
    raw_counts >= 10,
  branch_model_strict_gate =
    is_allocation_feature == TRUE &
    allocation_total_raw_counts >= 300 &
    clean_cds_raw_counts >= 100 &
    raw_counts >= 30
)]
units[unit_type == "merged_condition_group",
      unit_type_detail := fifelse(n_runs_merged >= 2,
                                  "merged_condition_group_n2plus",
                                  "merged_condition_group_singleton")]
units[unit_type == "single_run", unit_type_detail := "single_run"]

gate_summary <- units[
  is_allocation_feature == TRUE,
  .(
    rows = .N,
    genes = uniqueN(gene_symbol),
    features = uniqueN(paste(gene_symbol, tx_id, feature_id, sep = "|")),
    units = uniqueN(unit_id),
    median_raw_counts = median_safe(raw_counts),
    q90_raw_counts = quantile_safe(raw_counts, 0.90),
    q95_raw_counts = quantile_safe(raw_counts, 0.95),
    fraction_ge_1 = mean(count_ge_1, na.rm = TRUE),
    fraction_ge_10 = mean(count_ge_10, na.rm = TRUE),
    fraction_ge_30 = mean(count_ge_30, na.rm = TRUE),
    fraction_ge_100 = mean(count_ge_100, na.rm = TRUE),
    fraction_minimal_gate = mean(branch_model_minimal_gate, na.rm = TRUE),
    fraction_strict_gate = mean(branch_model_strict_gate, na.rm = TRUE),
    median_allocation_total_counts =
      median_safe(allocation_total_raw_counts),
    fraction_allocation_total_ge_100 =
      mean(allocation_total_ge_100, na.rm = TRUE),
    fraction_clean_cds_ge_30 = mean(clean_cds_ge_30, na.rm = TRUE),
    median_n_runs_merged = median_safe(n_runs_merged)
  ),
  by = .(unit_type_detail, branch_feature_class, feature_length_bin)
]
setorder(gate_summary, unit_type_detail, branch_feature_class,
         feature_length_bin)
fwrite(gate_summary, file.path(output_dir, "branch_coverage_gate_summary.csv"))

feature_summary <- units[
  is_allocation_feature == TRUE,
  .(
    observations = .N,
    units = uniqueN(unit_id),
    units_ge_10 = sum(count_ge_10, na.rm = TRUE),
    units_ge_30 = sum(count_ge_30, na.rm = TRUE),
    units_ge_100 = sum(count_ge_100, na.rm = TRUE),
    minimal_gate_units = sum(branch_model_minimal_gate, na.rm = TRUE),
    strict_gate_units = sum(branch_model_strict_gate, na.rm = TRUE),
    fraction_ge_10 = mean(count_ge_10, na.rm = TRUE),
    fraction_ge_30 = mean(count_ge_30, na.rm = TRUE),
    fraction_minimal_gate = mean(branch_model_minimal_gate, na.rm = TRUE),
    fraction_strict_gate = mean(branch_model_strict_gate, na.rm = TRUE),
    median_raw_counts = median_safe(raw_counts),
    q90_raw_counts = quantile_safe(raw_counts, 0.90),
    max_raw_counts = max(raw_counts, na.rm = TRUE),
    median_relative_use = median_safe(relative_use_allocation_counts),
    max_relative_use = max(relative_use_allocation_counts, na.rm = TRUE),
    example_units = collapse_examples(unit_id, 4L)
  ),
  by = .(unit_type_detail, gene_symbol, tx_id, feature_id, feature_type,
         category, branch_feature_class, feature_bases, cds_overlap_bases,
         feature_length_bin, is_small_uorf)
]
setorder(feature_summary, unit_type_detail, branch_feature_class,
         feature_bases, gene_symbol, feature_id)
fwrite(feature_summary,
       file.path(output_dir, "branch_coverage_feature_summary.csv"))

gene_summary <- allocation_totals[, .(
  units = .N,
  median_allocation_total_counts = median_safe(allocation_total_raw_counts),
  q90_allocation_total_counts = quantile_safe(allocation_total_raw_counts, 0.90),
  median_clean_cds_counts = median_safe(clean_cds_raw_counts),
  q90_clean_cds_counts = quantile_safe(clean_cds_raw_counts, 0.90),
  median_max_uorf_counts = median_safe(max_uorf_raw_counts),
  q90_max_uorf_counts = quantile_safe(max_uorf_raw_counts, 0.90),
  fraction_allocation_total_ge_100 =
    mean(allocation_total_raw_counts >= 100, na.rm = TRUE),
  fraction_allocation_total_ge_300 =
    mean(allocation_total_raw_counts >= 300, na.rm = TRUE),
  fraction_clean_cds_ge_30 = mean(clean_cds_raw_counts >= 30, na.rm = TRUE),
  fraction_any_uorf_ge_10 = mean(n_uorf_features_ge_10 >= 1, na.rm = TRUE),
  fraction_any_uorf_ge_30 = mean(n_uorf_features_ge_30 >= 1, na.rm = TRUE),
  fraction_two_uorfs_ge_10 = mean(n_uorf_features_ge_10 >= 2, na.rm = TRUE),
  fraction_two_uorfs_ge_30 = mean(n_uorf_features_ge_30 >= 2, na.rm = TRUE)
), by = .(unit_type, gene_symbol, tx_id)]
fwrite(gene_summary,
       file.path(output_dir, "branch_coverage_gene_summary.csv"))

unit_summary <- allocation_totals[, .(
  allocation_gene_rows = .N,
  genes = uniqueN(gene_symbol),
  median_allocation_total_counts = median_safe(allocation_total_raw_counts),
  fraction_allocation_total_ge_100 =
    mean(allocation_total_raw_counts >= 100, na.rm = TRUE),
  fraction_allocation_total_ge_300 =
    mean(allocation_total_raw_counts >= 300, na.rm = TRUE),
  fraction_clean_cds_ge_30 = mean(clean_cds_raw_counts >= 30, na.rm = TRUE),
  fraction_any_uorf_ge_10 = mean(n_uorf_features_ge_10 >= 1, na.rm = TRUE),
  fraction_any_uorf_ge_30 = mean(n_uorf_features_ge_30 >= 1, na.rm = TRUE),
  fraction_two_uorfs_ge_10 = mean(n_uorf_features_ge_10 >= 2, na.rm = TRUE)
), by = .(unit_type)]
fwrite(unit_summary,
       file.path(output_dir, "branch_coverage_unit_summary.csv"))

replicate_group_summary <- merged_units[
  is_allocation_feature == TRUE,
  .(
    allocation_features = .N,
    genes = uniqueN(gene_symbol),
    n_runs_merged = max(n_runs_merged, na.rm = TRUE),
    study = study[1],
    AUTHOR = AUTHOR[1],
    CELL_LINE = CELL_LINE[1],
    TISSUE = TISSUE[1],
    GENE = GENE[1],
    CONDITION = CONDITION[1],
    design_family = design_family[1],
    dominant_state = dominant_state[1]
  ),
  by = unit_id
]
setorder(replicate_group_summary, -n_runs_merged, study, CELL_LINE, CONDITION)
fwrite(replicate_group_summary,
       file.path(output_dir, "branch_coverage_replicate_group_summary.csv"))

small_uorf_feature_wide <- dcast(
  feature_summary[is_small_uorf == TRUE],
  gene_symbol + tx_id + feature_id + feature_type + branch_feature_class +
    feature_bases + cds_overlap_bases + feature_length_bin ~ unit_type_detail,
  value.var = c("units", "units_ge_10", "units_ge_30",
                "minimal_gate_units", "strict_gate_units",
                "median_raw_counts", "q90_raw_counts", "max_raw_counts"),
  fill = 0
)
fwrite(small_uorf_feature_wide,
       file.path(output_dir, "branch_coverage_small_uorf_single_vs_merged.csv"))

recommendation <- data.table(
  criterion = c(
    "single_run_small_uorf_fraction_ge_10",
    "single_run_small_uorf_fraction_ge_30",
    "merged_n2_small_uorf_fraction_ge_10",
    "merged_n2_small_uorf_fraction_ge_30",
    "single_run_allocation_feature_minimal_gate",
    "merged_n2_allocation_feature_minimal_gate",
    "replicate_groups_with_n2plus",
    "recommendation"
  ),
  value = NA_character_
)
single_small <- gate_summary[
  unit_type_detail == "single_run" &
    branch_feature_class %chin% c("leader_uORF", "overlapping_uORF") &
    feature_length_bin %chin% c("<=15 nt", "16-30 nt", "31-60 nt")
]
merged_small <- gate_summary[
  unit_type_detail == "merged_condition_group_n2plus" &
    branch_feature_class %chin% c("leader_uORF", "overlapping_uORF") &
    feature_length_bin %chin% c("<=15 nt", "16-30 nt", "31-60 nt")
]
single_alloc <- gate_summary[unit_type_detail == "single_run"]
merged_alloc <- gate_summary[unit_type_detail == "merged_condition_group_n2plus"]
single_small_ge10 <- weighted.mean(single_small$fraction_ge_10,
                                   single_small$rows, na.rm = TRUE)
single_small_ge30 <- weighted.mean(single_small$fraction_ge_30,
                                   single_small$rows, na.rm = TRUE)
merged_small_ge10 <- weighted.mean(merged_small$fraction_ge_10,
                                   merged_small$rows, na.rm = TRUE)
merged_small_ge30 <- weighted.mean(merged_small$fraction_ge_30,
                                   merged_small$rows, na.rm = TRUE)
single_minimal <- weighted.mean(single_alloc$fraction_minimal_gate,
                                single_alloc$rows, na.rm = TRUE)
merged_minimal <- weighted.mean(merged_alloc$fraction_minimal_gate,
                                merged_alloc$rows, na.rm = TRUE)
n_merged_n2 <- uniqueN(merged_units[n_runs_merged >= 2, unit_id])
recommended_unit <- if (is.finite(merged_small_ge10) &&
                        merged_small_ge10 >= 0.35 &&
                        merged_minimal >= single_minimal * 1.25) {
  "Use merged same-condition replicate groups for branch modeling; do not use single runs as primary units."
} else if (is.finite(merged_small_ge10) && merged_small_ge10 > single_small_ge10) {
  "Merged replicate groups improve small-uORF coverage, but coverage remains marginal; restrict to covered features."
} else {
  "Coverage is insufficient for a general branch-usage model; restrict to manually reviewed high-count genes."
}
recommendation[criterion == "single_run_small_uorf_fraction_ge_10",
               value := sprintf("%.4f", single_small_ge10)]
recommendation[criterion == "single_run_small_uorf_fraction_ge_30",
               value := sprintf("%.4f", single_small_ge30)]
recommendation[criterion == "merged_n2_small_uorf_fraction_ge_10",
               value := sprintf("%.4f", merged_small_ge10)]
recommendation[criterion == "merged_n2_small_uorf_fraction_ge_30",
               value := sprintf("%.4f", merged_small_ge30)]
recommendation[criterion == "single_run_allocation_feature_minimal_gate",
               value := sprintf("%.4f", single_minimal)]
recommendation[criterion == "merged_n2_allocation_feature_minimal_gate",
               value := sprintf("%.4f", merged_minimal)]
recommendation[criterion == "replicate_groups_with_n2plus",
               value := as.character(n_merged_n2)]
recommendation[criterion == "recommendation", value := recommended_unit]
fwrite(recommendation,
       file.path(output_dir, "branch_coverage_recommendation.csv"))

summary_metrics <- data.table(
  metric = c(
    "feature_rows",
    "runs",
    "genes",
    "features",
    "allocation_features",
    "small_uorf_features",
    "merged_condition_groups",
    "merged_condition_groups_n2plus",
    "single_small_uorf_fraction_ge_10",
    "single_small_uorf_fraction_ge_30",
    "merged_n2_small_uorf_fraction_ge_10",
    "merged_n2_small_uorf_fraction_ge_30",
    "single_allocation_feature_minimal_gate_fraction",
    "merged_n2_allocation_feature_minimal_gate_fraction"
  ),
  value = c(
    nrow(fx),
    uniqueN(fx$Run),
    uniqueN(fx$gene_symbol),
    uniqueN(paste(fx$gene_symbol, fx$tx_id, fx$feature_id, sep = "|")),
    uniqueN(paste(fx[is_allocation_feature == TRUE, gene_symbol],
                  fx[is_allocation_feature == TRUE, tx_id],
                  fx[is_allocation_feature == TRUE, feature_id], sep = "|")),
    uniqueN(paste(fx[is_small_uorf == TRUE, gene_symbol],
                  fx[is_small_uorf == TRUE, tx_id],
                  fx[is_small_uorf == TRUE, feature_id], sep = "|")),
    uniqueN(merged_units$unit_id),
    uniqueN(merged_units[n_runs_merged >= 2, unit_id]),
    single_small_ge10,
    single_small_ge30,
    merged_small_ge10,
    merged_small_ge30,
    single_minimal,
    merged_minimal
  )
)
fwrite(summary_metrics,
       file.path(output_dir, "branch_coverage_summary_metrics.csv"))

plot_gate <- gate_summary[
  branch_feature_class %chin% c("leader_uORF", "overlapping_uORF",
                                "clean_CDS") &
    unit_type_detail %chin% c("single_run",
                              "merged_condition_group_n2plus")
]
plot_gate[, feature_length_bin := factor(
  feature_length_bin,
  levels = c("<=15 nt", "16-30 nt", "31-60 nt", "61-120 nt",
             ">120 nt", "unknown")
)]
if (nrow(plot_gate)) {
  p_gate <- ggplot(
    plot_gate,
    aes(x = feature_length_bin, y = fraction_ge_30,
        fill = unit_type_detail,
        text = paste0(
          "Unit: ", unit_type_detail,
          "<br>Feature class: ", branch_feature_class,
          "<br>Length bin: ", feature_length_bin,
          "<br>Rows: ", rows,
          "<br>Fraction >=30 reads: ", round(fraction_ge_30, 3),
          "<br>Fraction >=10 reads: ", round(fraction_ge_10, 3),
          "<br>Median counts: ", round(median_raw_counts, 1)
        ))
  ) +
    geom_col(position = position_dodge(width = 0.72), width = 0.68) +
    facet_wrap(~ branch_feature_class, scales = "free_x") +
    scale_fill_manual(values = c(
      "single_run" = "#6b7280",
      "merged_condition_group_n2plus" = "#047857"
    )) +
    labs(
      title = "RDG Feature Coverage: Single Runs vs Merged Replicate Groups",
      subtitle = "Fraction of feature observations with at least 30 raw reads",
      x = "Feature length bin",
      y = "Fraction >=30 reads",
      fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "bottom"
    )
  ggsave(file.path(figure_dir, "branch_coverage_count_gate_by_length.png"),
         p_gate, width = 10.5, height = 6.8, dpi = 180)
  ggsave(file.path(figure_dir, "branch_coverage_count_gate_by_length.pdf"),
         p_gate, width = 10.5, height = 6.8)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_gate,
      file.path(figure_dir, "branch_coverage_count_gate_by_length.html"),
      title = "RDG branch coverage by length"
    )
  }
}

small_top <- feature_summary[
  is_small_uorf == TRUE &
    unit_type_detail %chin% c("single_run", "merged_condition_group_n2plus")
]
if (nrow(small_top)) {
  small_top[, feature_label := paste0(gene_symbol, " ", feature_id,
                                      " (", feature_bases, " nt)")]
  rank_features <- small_top[
    unit_type_detail == "merged_condition_group_n2plus",
    .(score = max(units_ge_30, na.rm = TRUE)),
    by = feature_label
  ][order(-score)][seq_len(min(.N, 40)), feature_label]
  small_plot <- small_top[feature_label %chin% rank_features]
  small_plot[, feature_label := factor(feature_label, levels = rev(rank_features))]
  p_small <- ggplot(
    small_plot,
    aes(x = units_ge_30, y = feature_label, fill = unit_type_detail,
        text = paste0(
          feature_label,
          "<br>Unit: ", unit_type_detail,
          "<br>Units >=30 reads: ", units_ge_30,
          "<br>Units >=10 reads: ", units_ge_10,
          "<br>Median counts: ", round(median_raw_counts, 1),
          "<br>Q90 counts: ", round(q90_raw_counts, 1)
        ))
  ) +
    geom_col(position = position_dodge(width = 0.72), width = 0.65) +
    scale_fill_manual(values = c(
      "single_run" = "#6b7280",
      "merged_condition_group_n2plus" = "#047857"
    )) +
    labs(
      title = "Covered Small uORFs After Replicate Merging",
      subtitle = "Top small uORF features by merged replicate groups with >=30 reads",
      x = "Number of units with >=30 raw reads",
      y = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "branch_coverage_small_uorf_top_features.png"),
         p_small, width = 10.5, height = 9.0, dpi = 180)
  ggsave(file.path(figure_dir, "branch_coverage_small_uorf_top_features.pdf"),
         p_small, width = 10.5, height = 9.0)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_small,
      file.path(figure_dir, "branch_coverage_small_uorf_top_features.html"),
      title = "Small uORF coverage after replicate merging"
    )
  }
}

message("Saved branch coverage QC outputs in: ", output_dir)
print(recommendation)
