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
output_dir <- file.path(analysis_dir, "dominant_infected_off_frame_cds_scan")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)

message("Infected off-frame CDS peak scan")
message("  1. Reuse cached all-gene CDS frame/bumpiness metrics")
message("  2. Sanitize infected/control labels from sample title/source text")
message("  3. Prioritize off-frame top-peak enrichment over bumpiness alone")
message("  4. Rank gene x infected-context candidates against same-run housekeeping")

run_metrics_file <- file.path(
  analysis_dir, "dominant_cds_frame_bumpiness_qc",
  "cds_frame_bumpiness_run_metrics.csv"
)
housekeeping_file <- file.path(
  analysis_dir, "dominant_cds_frame_bumpiness_qc",
  "housekeeping_baseline_run_qc.csv"
)
if (!file.exists(run_metrics_file)) {
  stop("Missing cached all-gene CDS metrics: ", run_metrics_file,
       call. = FALSE)
}
if (!file.exists(housekeeping_file)) {
  stop("Missing housekeeping run QC: ", housekeeping_file, call. = FALSE)
}

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

collapse_text <- function(x, max_n = 6L) {
  x <- unique(safe_chr(x))
  x <- x[nzchar(x)]
  if (!length(x)) return(NA_character_)
  suffix <- if (length(x) > max_n) ";..." else ""
  paste0(paste(head(x, max_n), collapse = ";"), suffix)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_sum <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(0)
  sum(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

safe_wilcox_p <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2L || length(y) < 2L) return(NA_real_)
  suppressWarnings(stats::wilcox.test(x, y)$p.value)
}

first_existing_cols <- function(file, wanted) {
  available <- names(fread(file, nrows = 0, showProgress = FALSE))
  intersect(wanted, available)
}

add_metadata_frame_qc <- function(dt) {
  usage_cols <- paste0("Frame_usage_", 0:2)
  if ("metadata_possible_misshift" %in% names(dt)) {
    dt[is.na(metadata_possible_misshift), metadata_possible_misshift := FALSE]
    return(dt)
  }
  if (!all(usage_cols %in% names(dt))) {
    dt[, metadata_frame0_fraction := NA_real_]
    dt[, metadata_frame1_fraction := NA_real_]
    dt[, metadata_frame2_fraction := NA_real_]
    dt[, metadata_canonical_frame_fraction := NA_real_]
    dt[, metadata_dominant_frame := NA_character_]
    dt[, metadata_possible_misshift := FALSE]
    return(dt)
  }
  raw0 <- suppressWarnings(as.numeric(dt[["Frame_usage_0"]]))
  raw1 <- suppressWarnings(as.numeric(dt[["Frame_usage_1"]]))
  raw2 <- suppressWarnings(as.numeric(dt[["Frame_usage_2"]]))
  raw_total <- raw0 + raw1 + raw2
  usage_is_percent <- is.finite(raw_total) & raw_total > 1.5
  frame0 <- ifelse(usage_is_percent, raw0 / 100, raw0)
  frame1 <- ifelse(usage_is_percent, raw1 / 100, raw1)
  frame2 <- ifelse(usage_is_percent, raw2 / 100, raw2)
  available <- is.finite(frame0) & is.finite(frame1) & is.finite(frame2) &
    is.finite(raw_total) & raw_total > 0
  dominant_frame <- paste0(
    "frame",
    max.col(cbind(frame0, frame1, frame2), ties.method = "first") - 1L
  )
  dominant_frame[!available] <- NA_character_
  dt[, metadata_frame0_fraction := frame0]
  dt[, metadata_frame1_fraction := frame1]
  dt[, metadata_frame2_fraction := frame2]
  dt[!available, c("metadata_frame0_fraction", "metadata_frame1_fraction",
                   "metadata_frame2_fraction") := .(NA_real_, NA_real_,
                                                    NA_real_)]
  dt[, metadata_canonical_frame_fraction := metadata_frame0_fraction]
  dt[, metadata_dominant_frame := dominant_frame]
  dt[, metadata_possible_misshift :=
       is.finite(metadata_canonical_frame_fraction) &
       metadata_canonical_frame_fraction < 0.45 &
       metadata_dominant_frame != "frame0"]
  dt[is.na(metadata_possible_misshift), metadata_possible_misshift := FALSE]
  dt
}

regex_has <- function(x, pattern) {
  grepl(pattern, x, ignore.case = TRUE, perl = TRUE)
}

classify_pathogen <- function(text) {
  fifelse(
    regex_has(text, "toxoplasma|tachyzoite|\\brh[_ -]?rep|\\brh[_ -]?tox"),
    "Toxoplasma",
    fifelse(
      regex_has(text, "oc43|hcov[-_ ]?oc43"),
      "HCoV-OC43",
      fifelse(
        regex_has(text, "sars[-_ ]?cov[-_ ]?2|covid"),
        "SARS-CoV-2",
        fifelse(
          regex_has(text, "influenza|\\biav\\b|\\bpr8\\b"),
          "Influenza/IAV",
          fifelse(
            regex_has(text, "\\bebv\\b|epstein"),
            "EBV",
            fifelse(
              regex_has(text, "herpes|hsv|cmv"),
              "Herpesvirus",
              fifelse(
                regex_has(text, "hcv|hepatitis c"),
                "HCV",
                fifelse(
                  regex_has(text, "hiv"),
                  "HIV",
                  fifelse(
                    regex_has(text, "virus|viral|infect|vaccinia|zika|dengue|sendai|chikungunya|adenovirus|norovirus"),
                    "infection/viral",
                    "none_or_unspecified"
                  )
                )
              )
            )
          )
        )
      )
    )
  )
}

classify_growth <- function(text) {
  fifelse(
    regex_has(text, "sub[-_ ]?confluent"),
    "subconfluent",
    fifelse(
      regex_has(text, "confluent"),
      "confluent",
      fifelse(
        regex_has(text, "tumou?r|cancer|carcinoma"),
        "tumor/cancer",
        "unspecified"
      )
    )
  )
}

classify_infection_status <- function(text) {
  control_regex <- paste(
    c("uninfected", "mock", "sham", "vehicle", "untreated",
      "\\bctrl\\b", "\\bcontrol\\b", "\\bwt\\b", "wild[-_ ]?type"),
    collapse = "|"
  )
  infection_regex <- paste(
    c("infect", "infection", "virus", "viral", "hcov", "sars", "covid",
      "influenza", "\\biav\\b", "\\bpr8\\b", "toxoplasma", "tachyzoite",
      "vaccinia", "zika", "dengue", "sendai", "chikungunya", "adenovirus",
      "norovirus", "\\bebv\\b", "epstein", "herpes", "hsv", "cmv",
      "hcv", "hiv"),
    collapse = "|"
  )
  is_control <- regex_has(text, control_regex)
  is_infected <- !is_control & regex_has(text, infection_regex)
  fifelse(is_control, "control",
          fifelse(is_infected, "infected", "other"))
}

make_clean_context_value <- function(x, fallback = "NA") {
  x <- safe_chr(x)
  fallback <- rep_len(safe_chr(fallback), length(x))
  empty <- !nzchar(x)
  x[empty] <- fallback[empty]
  x[!nzchar(x)] <- "NA"
  x
}

run_cols <- c(
  "Run", "gene_symbol", "tx_id", "region",
  "cds_positions_used", "cds_total_counts", "cds_mean_coverage",
  "cds_coverage_cv", "cds_zero_fraction", "cds_max_position_share",
  "cds_top1pct_share", "cds_top1pct_signal",
  "frame0_fraction", "frame1_fraction", "frame2_fraction",
  "top1pct_frame0_fraction", "top1pct_frame1_fraction",
  "top1pct_frame2_fraction", "canonical_frame_fraction",
  "off_frame_fraction", "dominant_cds_frame",
  "dominant_cds_frame_fraction", "top1pct_canonical_fraction",
  "top1pct_off_frame_fraction", "possible_library_misshift_for_gene",
  "global_cds_middle_counts", "global_canonical_frame_fraction",
  "global_possible_misshift", "BioProject", "Study_Pubmed_id",
  "AUTHOR", "sample_source", "sample_title", "LIBRARYTYPE",
  "REPLICATE", "CONDITION", "INHIBITOR", "TIMEPOINT", "TISSUE",
  "CELL_LINE", "GENE", "FRACTION", "Cancer_type", "Cell_model",
  "Cell_type", "Organ_system", "Sex", "Life_stage", "study",
  "sample", "sample_id", "Read_Length_Mode", "Pct_Reads_At_Mode_Len",
  "ProteinCoding_Pct", "rRNA_Pct", "Frame_usage_0", "Frame_usage_1",
  "Frame_usage_2", "metadata_frame0_fraction", "metadata_frame1_fraction",
  "metadata_frame2_fraction", "metadata_canonical_frame_fraction",
  "metadata_dominant_frame", "metadata_possible_misshift",
  "top_readlength", "drop_off", "drop_in",
  "tis_peak_strength", "cds_strength_vs_15ntpluss",
  "top1pct_bumpiness_z_vs_depth", "cv_bumpiness_z_vs_depth"
)
run_cols <- first_existing_cols(run_metrics_file, run_cols)
message("Reading selected run metric columns from: ", run_metrics_file)
dt <- fread(run_metrics_file, select = run_cols, showProgress = TRUE)
dt[, gene_symbol := toupper(gene_symbol)]
dt <- add_metadata_frame_qc(dt)

hk_cols <- c(
  "Run", "housekeeping_n_genes", "housekeeping_total_counts",
  "housekeeping_median_canonical_frame_fraction",
  "housekeeping_median_top1pct_share",
  "housekeeping_median_top1pct_off_frame_fraction",
  "housekeeping_median_top1pct_bumpiness_z_vs_depth",
  "housekeeping_median_coverage_cv",
  "housekeeping_n_gene_misshift_candidates",
  "housekeeping_canonical_frame_fraction",
  "housekeeping_possible_misshift",
  "global_canonical_frame_fraction", "global_possible_misshift"
)
hk_cols <- first_existing_cols(housekeeping_file, hk_cols)
message("Reading housekeeping baseline from: ", housekeeping_file)
hk <- fread(housekeeping_file, select = hk_cols, showProgress = FALSE)
hk <- unique(hk, by = "Run")
names(hk) <- sub("^global_", "hk_global_", names(hk))
dt <- merge(dt, hk, by = "Run", all.x = TRUE, sort = FALSE)

for (nm in c("sample_title", "sample_source", "CONDITION", "INHIBITOR",
             "TIMEPOINT", "TISSUE", "CELL_LINE", "GENE", "FRACTION",
             "Cancer_type", "Cell_model", "Cell_type", "Organ_system",
             "Sex", "Life_stage", "study", "BioProject", "AUTHOR",
             "LIBRARYTYPE")) {
  if (!nm %in% names(dt)) dt[, (nm) := NA_character_]
}

dt[, context_text := tolower(paste(
  safe_chr(sample_title), safe_chr(sample_source), safe_chr(CONDITION),
  safe_chr(INHIBITOR), safe_chr(TIMEPOINT), safe_chr(TISSUE),
  safe_chr(CELL_LINE), safe_chr(GENE), safe_chr(FRACTION),
  safe_chr(Cancer_type), safe_chr(Cell_model), safe_chr(Cell_type),
  safe_chr(Organ_system), safe_chr(study), safe_chr(BioProject),
  safe_chr(AUTHOR), safe_chr(LIBRARYTYPE),
  sep = " | "
))]
dt[, infection_status := classify_infection_status(context_text)]
dt[, pathogen := classify_pathogen(context_text)]
dt[, growth_context := classify_growth(context_text)]
dt[, inhibitor_context := make_clean_context_value(INHIBITOR, "none/unspecified")]
dt[, cell_context := make_clean_context_value(CELL_LINE, "unknown_cell")]
dt[, tissue_context := make_clean_context_value(TISSUE, "unknown_tissue")]
dt[, study_context := make_clean_context_value(study, BioProject)]
dt[!nzchar(study_context), study_context := make_clean_context_value(BioProject)]
dt[, bioproject_context := make_clean_context_value(BioProject, "unknown_project")]
dt[, author_context := make_clean_context_value(AUTHOR, "unknown_author")]
dt[, timepoint_context := make_clean_context_value(TIMEPOINT, "unknown_time")]
dt[, match_context_id := paste(
  bioproject_context, author_context, cell_context, tissue_context,
  inhibitor_context, growth_context,
  sep = " | "
)]
dt[, timepoint_context_id := paste(
  match_context_id, timepoint_context,
  sep = " | "
)]

if (!"top1pct_off_frame_fraction" %in% names(dt)) {
  dt[, top1pct_off_frame_fraction :=
       top1pct_frame1_fraction + top1pct_frame2_fraction]
}
if (!"off_frame_fraction" %in% names(dt)) {
  dt[, off_frame_fraction := frame1_fraction + frame2_fraction]
}
dt[, top1pct_off_frame_enrichment :=
     top1pct_off_frame_fraction - off_frame_fraction]
dt[, top1pct_off_frame_ratio :=
     (top1pct_off_frame_fraction + 0.02) / (off_frame_fraction + 0.02)]
dt[, target_minus_housekeeping_bumpiness_z :=
     top1pct_bumpiness_z_vs_depth -
     housekeeping_median_top1pct_bumpiness_z_vs_depth]
dt[, target_minus_housekeeping_canonical_frame :=
     canonical_frame_fraction - housekeeping_canonical_frame_fraction]
dt[, technical_frame_bias_possible :=
     fifelse(is.na(metadata_possible_misshift), FALSE,
             metadata_possible_misshift) |
     fifelse(is.na(global_possible_misshift), FALSE,
             global_possible_misshift) |
     fifelse(is.na(housekeeping_possible_misshift), FALSE,
             housekeeping_possible_misshift)]

min_cds_counts <- 100
min_top1_signal <- 20
dt[, evaluable_for_peak_scan :=
     cds_total_counts >= min_cds_counts &
     cds_top1pct_signal >= min_top1_signal &
     is.finite(top1pct_off_frame_fraction) &
     is.finite(top1pct_off_frame_enrichment) &
     is.finite(target_minus_housekeeping_bumpiness_z)]
dt[, candidate_off_frame_peak_before_technical_filter :=
     evaluable_for_peak_scan &
     top1pct_off_frame_fraction >= 0.55 &
     top1pct_off_frame_enrichment >= 0.20 &
     target_minus_housekeeping_bumpiness_z >= 0.25]
dt[, candidate_off_frame_peak_no_hk_artifact :=
     candidate_off_frame_peak_before_technical_filter &
     !technical_frame_bias_possible]
dt[, candidate_off_frame_peak_frame_risk :=
     candidate_off_frame_peak_before_technical_filter &
     technical_frame_bias_possible]
dt[, candidate_bumpy_no_hk_artifact :=
     evaluable_for_peak_scan &
     target_minus_housekeeping_bumpiness_z >= 1.0 &
     !technical_frame_bias_possible]
dt[, candidate_frame_shift_no_hk_artifact :=
     evaluable_for_peak_scan &
     canonical_frame_fraction < 0.45 &
     housekeeping_canonical_frame_fraction >= 0.50 &
     !technical_frame_bias_possible]
dt[, run_off_frame_peak_score :=
     pmax(0, top1pct_off_frame_enrichment) +
     0.20 * pmax(0, target_minus_housekeeping_bumpiness_z) +
     0.25 * as.numeric(candidate_frame_shift_no_hk_artifact)]

scan_dt <- dt[infection_status %chin% c("infected", "control")]
infected_contexts <- unique(scan_dt[infection_status == "infected",
                                    match_context_id])
scan_dt <- scan_dt[match_context_id %chin% infected_contexts]
message(
  "Rows in infected/control scan table after keeping infected contexts: ",
  nrow(scan_dt), " across ", length(infected_contexts), " contexts"
)

summarize_context <- function(x) {
  case <- x[infection_status == "infected" &
              evaluable_for_peak_scan == TRUE]
  ctrl <- x[infection_status == "control" &
              evaluable_for_peak_scan == TRUE]
  all_evaluable <- x[evaluable_for_peak_scan == TRUE]
  data.table(
    case_n_evaluable = nrow(case),
    control_n_evaluable = nrow(ctrl),
    case_n_off_frame_peak =
      as.integer(case[, sum(candidate_off_frame_peak_no_hk_artifact, na.rm = TRUE)]),
    control_n_off_frame_peak =
      as.integer(ctrl[, sum(candidate_off_frame_peak_no_hk_artifact, na.rm = TRUE)]),
    case_n_bumpy =
      as.integer(case[, sum(candidate_bumpy_no_hk_artifact, na.rm = TRUE)]),
    control_n_bumpy =
      as.integer(ctrl[, sum(candidate_bumpy_no_hk_artifact, na.rm = TRUE)]),
    case_n_frame_shift =
      as.integer(case[, sum(candidate_frame_shift_no_hk_artifact, na.rm = TRUE)]),
    control_n_frame_shift =
      as.integer(ctrl[, sum(candidate_frame_shift_no_hk_artifact, na.rm = TRUE)]),
    case_fraction_off_frame_peak =
      if (nrow(case)) mean(case$candidate_off_frame_peak_no_hk_artifact, na.rm = TRUE) else NA_real_,
    control_fraction_off_frame_peak =
      if (nrow(ctrl)) mean(ctrl$candidate_off_frame_peak_no_hk_artifact, na.rm = TRUE) else NA_real_,
    case_fraction_bumpy =
      if (nrow(case)) mean(case$candidate_bumpy_no_hk_artifact, na.rm = TRUE) else NA_real_,
    control_fraction_bumpy =
      if (nrow(ctrl)) mean(ctrl$candidate_bumpy_no_hk_artifact, na.rm = TRUE) else NA_real_,
    case_median_top1pct_off_frame_fraction =
      as.numeric(case[, safe_median(top1pct_off_frame_fraction)]),
    control_median_top1pct_off_frame_fraction =
      as.numeric(ctrl[, safe_median(top1pct_off_frame_fraction)]),
    case_median_off_frame_enrichment =
      as.numeric(case[, safe_median(top1pct_off_frame_enrichment)]),
    control_median_off_frame_enrichment =
      as.numeric(ctrl[, safe_median(top1pct_off_frame_enrichment)]),
    case_median_target_minus_hk_bumpiness_z =
      as.numeric(case[, safe_median(target_minus_housekeeping_bumpiness_z)]),
    control_median_target_minus_hk_bumpiness_z =
      as.numeric(ctrl[, safe_median(target_minus_housekeeping_bumpiness_z)]),
    case_median_run_off_frame_peak_score =
      as.numeric(case[, safe_median(run_off_frame_peak_score)]),
    control_median_run_off_frame_peak_score =
      as.numeric(ctrl[, safe_median(run_off_frame_peak_score)]),
    case_median_cds_total_counts =
      as.numeric(case[, safe_median(cds_total_counts)]),
    control_median_cds_total_counts =
      as.numeric(ctrl[, safe_median(cds_total_counts)]),
    case_technical_frame_bias_fraction =
      if (nrow(case)) mean(case$technical_frame_bias_possible, na.rm = TRUE) else NA_real_,
    control_technical_frame_bias_fraction =
      if (nrow(ctrl)) mean(ctrl$technical_frame_bias_possible, na.rm = TRUE) else NA_real_,
    off_frame_score_wilcox_p =
      safe_wilcox_p(case$run_off_frame_peak_score,
                    ctrl$run_off_frame_peak_score),
    case_pathogens = collapse_text(case$pathogen, max_n = 4L),
    case_timepoints = collapse_text(case$timepoint_context, max_n = 6L),
    case_sample_titles = collapse_text(case$sample_title, max_n = 4L),
    control_sample_titles = collapse_text(ctrl$sample_title, max_n = 4L),
    top_case_runs = collapse_text(
      case[order(-run_off_frame_peak_score)]$Run,
      max_n = 8L
    ),
    n_evaluable_all = nrow(all_evaluable)
  )
}

context_summary <- scan_dt[
  ,
  summarize_context(.SD),
  by = .(
    gene_symbol, tx_id, match_context_id, bioproject_context,
    author_context, cell_context, tissue_context, inhibitor_context,
    growth_context
  )
]

context_summary[, case_minus_control_fraction_off_frame_peak :=
                  case_fraction_off_frame_peak -
                  fifelse(is.na(control_fraction_off_frame_peak),
                          0, control_fraction_off_frame_peak)]
context_summary[, case_minus_control_median_off_frame_enrichment :=
                  case_median_off_frame_enrichment -
                  control_median_off_frame_enrichment]
context_summary[, case_minus_control_median_off_frame_peak_score :=
                  case_median_run_off_frame_peak_score -
                  control_median_run_off_frame_peak_score]
context_summary[, q_off_frame_score_wilcox :=
                  p.adjust(off_frame_score_wilcox_p, method = "BH")]
context_summary[, context_signal_class := fifelse(
  case_n_evaluable >= 2 & control_n_evaluable >= 2 &
    case_n_off_frame_peak >= 2 &
    case_fraction_off_frame_peak >= control_fraction_off_frame_peak + 0.10 &
    case_median_run_off_frame_peak_score >
      control_median_run_off_frame_peak_score + 0.10 &
    fifelse(is.na(case_technical_frame_bias_fraction), 0,
            case_technical_frame_bias_fraction) < 0.50,
  "strong_matched_infected_off_frame_peak",
  fifelse(
    case_n_evaluable >= 2 & control_n_evaluable >= 2 &
      case_n_off_frame_peak >= 1 &
      case_fraction_off_frame_peak > control_fraction_off_frame_peak &
      case_median_run_off_frame_peak_score >
        control_median_run_off_frame_peak_score,
    "possible_matched_infected_off_frame_peak",
    fifelse(
      case_n_evaluable >= 2 & control_n_evaluable < 2 &
        case_n_off_frame_peak >= 2,
      "unmatched_recurrent_infected_off_frame_peak",
      fifelse(
        control_n_evaluable >= 2 &
          control_n_off_frame_peak >= pmax(2, case_n_off_frame_peak) &
          control_fraction_off_frame_peak >=
            fifelse(is.na(case_fraction_off_frame_peak), 0,
                    case_fraction_off_frame_peak),
        "control_or_baseline_off_frame_peak",
        fifelse(
          fifelse(is.na(case_technical_frame_bias_fraction), 0,
                  case_technical_frame_bias_fraction) >= 0.50 |
            fifelse(is.na(control_technical_frame_bias_fraction), 0,
                    control_technical_frame_bias_fraction) >= 0.50,
          "technical_frame_bias_risk",
          fifelse(
            case_n_evaluable < 2,
            "insufficient_infected_depth",
            "no_clear_infected_signal"
          )
        )
      )
    )
  )
)]

timepoint_summary <- scan_dt[
  infection_status == "infected" & evaluable_for_peak_scan,
  .(
    infected_n_evaluable = .N,
    infected_n_off_frame_peak =
      as.integer(sum(candidate_off_frame_peak_no_hk_artifact, na.rm = TRUE)),
    infected_n_bumpy =
      as.integer(sum(candidate_bumpy_no_hk_artifact, na.rm = TRUE)),
    infected_fraction_off_frame_peak =
      mean(candidate_off_frame_peak_no_hk_artifact, na.rm = TRUE),
    infected_median_top1pct_off_frame_fraction =
      as.numeric(safe_median(top1pct_off_frame_fraction)),
    infected_median_off_frame_enrichment =
      as.numeric(safe_median(top1pct_off_frame_enrichment)),
    infected_median_run_off_frame_peak_score =
      as.numeric(safe_median(run_off_frame_peak_score)),
    infected_median_target_minus_hk_bumpiness_z =
      as.numeric(safe_median(target_minus_housekeeping_bumpiness_z)),
    infected_median_counts = as.numeric(safe_median(cds_total_counts)),
    pathogens = collapse_text(pathogen, max_n = 3L),
    sample_titles = collapse_text(sample_title, max_n = 4L),
    top_runs = collapse_text(Run[order(-run_off_frame_peak_score)], max_n = 6L)
  ),
  by = .(
    gene_symbol, tx_id, timepoint_context_id, bioproject_context,
    author_context, cell_context, tissue_context, inhibitor_context,
    growth_context, timepoint_context
  )
][infected_n_evaluable >= 2]
setorder(timepoint_summary, -infected_n_off_frame_peak,
         -infected_fraction_off_frame_peak,
         -infected_median_run_off_frame_peak_score)

gene_scores <- context_summary[
  ,
  .(
    n_contexts = .N,
    n_evaluable_contexts = sum(case_n_evaluable > 0, na.rm = TRUE),
    n_strong_matched_contexts =
      sum(context_signal_class == "strong_matched_infected_off_frame_peak",
          na.rm = TRUE),
    n_possible_matched_contexts =
      sum(context_signal_class == "possible_matched_infected_off_frame_peak",
          na.rm = TRUE),
    n_unmatched_recurrent_contexts =
      sum(context_signal_class == "unmatched_recurrent_infected_off_frame_peak",
          na.rm = TRUE),
    n_technical_risk_contexts =
      sum(context_signal_class == "technical_frame_bias_risk", na.rm = TRUE),
    total_case_n_evaluable = sum(case_n_evaluable, na.rm = TRUE),
    total_control_n_evaluable = sum(control_n_evaluable, na.rm = TRUE),
    total_case_n_off_frame_peak = sum(case_n_off_frame_peak, na.rm = TRUE),
    total_control_n_off_frame_peak =
      sum(control_n_off_frame_peak, na.rm = TRUE),
    weighted_case_fraction_off_frame_peak =
      safe_sum(case_n_off_frame_peak) / pmax(safe_sum(case_n_evaluable), 1),
    weighted_control_fraction_off_frame_peak =
      safe_sum(control_n_off_frame_peak) /
      pmax(safe_sum(control_n_evaluable), 1),
    max_context_delta_off_frame_peak_fraction =
      safe_max(case_minus_control_fraction_off_frame_peak),
    max_context_delta_off_frame_peak_score =
      safe_max(case_minus_control_median_off_frame_peak_score),
    median_case_top1pct_off_frame_fraction =
      safe_median(case_median_top1pct_off_frame_fraction),
    median_case_off_frame_enrichment =
      safe_median(case_median_off_frame_enrichment),
    median_case_target_minus_hk_bumpiness_z =
      safe_median(case_median_target_minus_hk_bumpiness_z),
    distinct_projects_with_case_signal =
      uniqueN(bioproject_context[case_n_off_frame_peak > 0]),
    distinct_cell_contexts_with_case_signal =
      uniqueN(cell_context[case_n_off_frame_peak > 0]),
    distinct_pathogens_with_case_signal =
      uniqueN(unlist(strsplit(
        paste(case_pathogens[case_n_off_frame_peak > 0], collapse = ";"),
        ";", fixed = TRUE
      ))),
    top_context_label = {
      idx <- order(-case_n_off_frame_peak,
                   -case_minus_control_fraction_off_frame_peak,
                   -case_median_run_off_frame_peak_score)[1L]
      paste(bioproject_context[idx], author_context[idx],
            cell_context[idx], tissue_context[idx],
            inhibitor_context[idx], growth_context[idx], sep = " | ")
    },
    top_context_class = {
      idx <- order(-case_n_off_frame_peak,
                   -case_minus_control_fraction_off_frame_peak,
                   -case_median_run_off_frame_peak_score)[1L]
      context_signal_class[idx]
    },
    top_context_case_runs = {
      idx <- order(-case_n_off_frame_peak,
                   -case_minus_control_fraction_off_frame_peak,
                   -case_median_run_off_frame_peak_score)[1L]
      top_case_runs[idx]
    },
    top_case_pathogens = collapse_text(case_pathogens[case_n_off_frame_peak > 0],
                                      max_n = 5L)
  ),
  by = .(gene_symbol, tx_id)
]

gene_scores[, weighted_delta_fraction_off_frame_peak :=
              weighted_case_fraction_off_frame_peak -
              weighted_control_fraction_off_frame_peak]
gene_scores[, candidate_score :=
              6 * n_strong_matched_contexts +
              3 * n_possible_matched_contexts +
              2 * n_unmatched_recurrent_contexts +
              1.5 * log1p(total_case_n_off_frame_peak) +
              4 * pmax(0, weighted_delta_fraction_off_frame_peak) +
              pmin(3, pmax(0, max_context_delta_off_frame_peak_score)) -
              1.5 * n_technical_risk_contexts]
gene_scores[, review_tier := fifelse(
  !is.finite(candidate_score) | candidate_score <= 0,
  "not_prioritized",
  fifelse(
    n_strong_matched_contexts >= 2 |
      (n_strong_matched_contexts >= 1 &
         total_case_n_off_frame_peak >= 3),
    "high_priority_review",
    fifelse(
      n_strong_matched_contexts >= 1 |
        n_possible_matched_contexts >= 2 |
        (n_possible_matched_contexts >= 1 &
           n_unmatched_recurrent_contexts >= 1),
      "medium_priority_review",
      fifelse(
        n_possible_matched_contexts >= 1 |
          n_unmatched_recurrent_contexts >= 1 |
          total_case_n_off_frame_peak >= 3,
        "low_priority_review",
        "not_prioritized"
      )
    )
  )
)]
setorder(gene_scores, -candidate_score, -total_case_n_off_frame_peak,
         -weighted_delta_fraction_off_frame_peak, gene_symbol)

top_runs <- scan_dt[
  infection_status == "infected" &
    evaluable_for_peak_scan &
    (candidate_off_frame_peak_no_hk_artifact |
       candidate_bumpy_no_hk_artifact |
       candidate_frame_shift_no_hk_artifact),
  .(
    gene_symbol, tx_id, Run, BioProject, study, AUTHOR, CELL_LINE, TISSUE,
    CONDITION, INHIBITOR, TIMEPOINT, pathogen, growth_context, sample_title,
    sample_source, cds_total_counts, cds_top1pct_signal,
    canonical_frame_fraction, off_frame_fraction,
    top1pct_off_frame_fraction, top1pct_off_frame_enrichment,
    top1pct_off_frame_ratio, target_minus_housekeeping_bumpiness_z,
    run_off_frame_peak_score, housekeeping_canonical_frame_fraction,
    housekeeping_median_top1pct_bumpiness_z_vs_depth,
    housekeeping_possible_misshift, metadata_canonical_frame_fraction,
    metadata_dominant_frame, metadata_possible_misshift,
    global_canonical_frame_fraction, global_possible_misshift,
    candidate_off_frame_peak_no_hk_artifact,
    candidate_bumpy_no_hk_artifact, candidate_frame_shift_no_hk_artifact
  )
]
setorder(top_runs, -candidate_off_frame_peak_no_hk_artifact,
         -run_off_frame_peak_score, -cds_total_counts)

frame_risk_runs <- scan_dt[
  infection_status == "infected" &
    candidate_off_frame_peak_frame_risk == TRUE,
  .(
    gene_symbol, tx_id, Run, BioProject, study, AUTHOR, CELL_LINE, TISSUE,
    CONDITION, INHIBITOR, TIMEPOINT, pathogen, growth_context, sample_title,
    sample_source, cds_total_counts, cds_top1pct_signal,
    canonical_frame_fraction, off_frame_fraction,
    top1pct_off_frame_fraction, top1pct_off_frame_enrichment,
    top1pct_off_frame_ratio, target_minus_housekeeping_bumpiness_z,
    run_off_frame_peak_score, housekeeping_canonical_frame_fraction,
    housekeeping_median_top1pct_bumpiness_z_vs_depth,
    housekeeping_possible_misshift, metadata_canonical_frame_fraction,
    metadata_dominant_frame, metadata_possible_misshift,
    global_canonical_frame_fraction, global_possible_misshift,
    candidate_off_frame_peak_before_technical_filter,
    candidate_off_frame_peak_frame_risk
  )
]
setorder(frame_risk_runs, -run_off_frame_peak_score,
         -cds_total_counts, gene_symbol)

frame_risk_gene_scores <- if (nrow(frame_risk_runs) > 0) {
  frame_risk_runs[
    ,
    .(
      n_frame_risk_runs = .N,
      n_projects = uniqueN(BioProject),
      n_cell_contexts = uniqueN(CELL_LINE),
      median_top1pct_off_frame_fraction =
        safe_median(top1pct_off_frame_fraction),
      median_off_frame_enrichment = safe_median(top1pct_off_frame_enrichment),
      median_target_minus_hk_bumpiness_z =
        safe_median(target_minus_housekeeping_bumpiness_z),
      median_metadata_canonical_frame_fraction =
        safe_median(metadata_canonical_frame_fraction),
      metadata_dominant_frames = collapse_text(metadata_dominant_frame,
                                               max_n = 4L),
      run_examples = collapse_text(Run[order(-run_off_frame_peak_score)],
                                   max_n = 8L),
      top_contexts = collapse_text(
        paste(BioProject, CELL_LINE, CONDITION, TIMEPOINT, sep = " | "),
        max_n = 6L
      )
    ),
    by = .(gene_symbol, tx_id)
  ][
    ,
    frame_risk_review_score :=
      log1p(n_frame_risk_runs) +
      pmax(0, median_top1pct_off_frame_fraction - 0.55) +
      pmax(0, median_off_frame_enrichment) +
      pmax(0, median_target_minus_hk_bumpiness_z)
  ][order(-frame_risk_review_score, -n_frame_risk_runs, gene_symbol)]
} else {
  data.table()
}

summary_metrics <- data.table(
  metric = c(
    "n_genes_scanned",
    "n_runs",
    "n_rows_total",
    "n_rows_infected_or_control",
    "n_infected_rows_evaluable",
    "n_control_rows_evaluable",
    "min_cds_counts",
    "min_top1pct_signal",
    "n_high_priority_genes",
    "n_medium_priority_genes",
    "n_low_priority_genes",
    "n_metadata_all_cds_misshift_rows",
    "n_any_technical_frame_bias_rows",
    "n_infected_candidate_off_frame_peak_runs",
    "n_infected_frame_risk_off_frame_peak_runs",
    "n_frame_risk_genes",
    "n_infected_candidate_bumpy_runs",
    "n_infected_candidate_frame_shift_runs"
  ),
  value = as.character(c(
    uniqueN(dt$gene_symbol),
    uniqueN(dt$Run),
    nrow(dt),
    nrow(scan_dt),
    scan_dt[infection_status == "infected" &
              evaluable_for_peak_scan, .N],
    scan_dt[infection_status == "control" &
              evaluable_for_peak_scan, .N],
    min_cds_counts,
    min_top1_signal,
    gene_scores[review_tier == "high_priority_review", .N],
    gene_scores[review_tier == "medium_priority_review", .N],
    gene_scores[review_tier == "low_priority_review", .N],
    dt[metadata_possible_misshift == TRUE, .N],
    dt[technical_frame_bias_possible == TRUE, .N],
    top_runs[candidate_off_frame_peak_no_hk_artifact == TRUE, .N],
    nrow(frame_risk_runs),
    if (nrow(frame_risk_gene_scores)) uniqueN(frame_risk_gene_scores$gene_symbol) else 0L,
    top_runs[candidate_bumpy_no_hk_artifact == TRUE, .N],
    top_runs[candidate_frame_shift_no_hk_artifact == TRUE, .N]
  ))
)

fwrite(gene_scores, file.path(output_dir,
                              "infected_off_frame_cds_candidate_gene_scores.csv"))
fwrite(context_summary, file.path(output_dir,
                                  "infected_off_frame_cds_candidate_contexts.csv"))
fwrite(timepoint_summary, file.path(output_dir,
                                    "infected_off_frame_cds_timepoint_contexts.csv"))
fwrite(top_runs, file.path(output_dir,
                           "infected_off_frame_cds_top_candidate_runs.csv"))
fwrite(frame_risk_runs, file.path(output_dir,
                                  "infected_off_frame_cds_frame_risk_runs.csv"))
fwrite(frame_risk_gene_scores, file.path(output_dir,
                                         "infected_off_frame_cds_frame_risk_gene_scores.csv"))
fwrite(summary_metrics, file.path(output_dir,
                                  "infected_off_frame_cds_summary_metrics.csv"))

top_gene_plot_dt <- head(gene_scores[review_tier != "not_prioritized"], 35L)
if (nrow(top_gene_plot_dt) > 0) {
  top_gene_plot_dt[, gene_label := factor(
    gene_symbol, levels = rev(gene_symbol)
  )]
  p <- ggplot(top_gene_plot_dt,
              aes(x = gene_label, y = candidate_score,
                  fill = review_tier,
                  text = paste0(
                    gene_symbol,
                    "<br>score=", round(candidate_score, 2),
                    "<br>tier=", review_tier,
                    "<br>case off-frame peak runs=",
                    total_case_n_off_frame_peak,
                    "<br>strong contexts=", n_strong_matched_contexts,
                    "<br>possible contexts=", n_possible_matched_contexts,
                    "<br>top context=", top_context_label
                  ))) +
    geom_col(width = 0.72) +
    coord_flip() +
    scale_fill_manual(
      values = c(
        high_priority_review = "#8b0000",
        medium_priority_review = "#d95f02",
        low_priority_review = "#6baed6"
      ),
      drop = FALSE
    ) +
    labs(
      title = "Infected CDS Off-Frame Peak Candidate Ranking",
      subtitle = paste0(
        "Off-frame top peaks after same-run housekeeping/global frame QC; ",
        "min CDS counts=", min_cds_counts,
        ", min top signal=", min_top1_signal
      ),
      x = NULL,
      y = "Candidate score",
      fill = "Review tier"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  ggsave(file.path(figure_dir, "infected_off_frame_cds_top_gene_scores.png"),
         p, width = 9.5, height = 7.5, dpi = 150)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p,
      file.path(figure_dir, "infected_off_frame_cds_top_gene_scores.html"),
      title = "infected_off_frame_cds_top_gene_scores"
    )
  }
}

top_context_plot_dt <- context_summary[
  context_signal_class %chin% c(
    "strong_matched_infected_off_frame_peak",
    "possible_matched_infected_off_frame_peak",
    "unmatched_recurrent_infected_off_frame_peak"
  )
]
if (nrow(top_context_plot_dt) > 0) {
  top_genes <- gene_scores[review_tier != "not_prioritized",
                           head(gene_symbol, 25L)]
  top_context_plot_dt <- top_context_plot_dt[gene_symbol %chin% top_genes]
  if (nrow(top_context_plot_dt) > 0) {
    top_context_plot_dt[, context_label := paste(
      bioproject_context, cell_context, inhibitor_context, growth_context,
      sep = "\n"
    )]
    context_order <- head(top_context_plot_dt[
      ,
      .(score = max(case_minus_control_fraction_off_frame_peak, na.rm = TRUE),
        n = sum(case_n_off_frame_peak, na.rm = TRUE)),
      by = context_label
    ][order(-n, -score)], 18L)$context_label
    top_context_plot_dt <- top_context_plot_dt[context_label %chin% context_order]
    top_context_plot_dt[, gene_label := factor(
      gene_symbol, levels = rev(top_genes)
    )]
    top_context_plot_dt[, context_label := factor(
      context_label, levels = rev(context_order)
    )]
    p2 <- ggplot(
      top_context_plot_dt,
      aes(
        x = context_label,
        y = gene_label,
        size = case_n_off_frame_peak,
        fill = case_minus_control_fraction_off_frame_peak,
        text = paste0(
          gene_symbol,
          "<br>context=", gsub("\n", " | ", as.character(context_label)),
          "<br>class=", context_signal_class,
          "<br>case n=", case_n_evaluable,
          "<br>control n=", control_n_evaluable,
          "<br>case off-frame peak n=", case_n_off_frame_peak,
          "<br>control off-frame peak n=", control_n_off_frame_peak,
          "<br>case-control fraction delta=",
          round(case_minus_control_fraction_off_frame_peak, 3),
          "<br>case pathogens=", case_pathogens,
          "<br>case timepoints=", case_timepoints
        )
      )
    ) +
      geom_point(shape = 21, color = "#333333", alpha = 0.90) +
      scale_fill_gradient2(
        low = "#2166ac", mid = "white", high = "#8b0000",
        midpoint = 0, na.value = "grey90"
      ) +
      scale_size_continuous(range = c(2.5, 10)) +
      labs(
        title = "Infected Contexts with Off-Frame CDS Peak Signals",
        subtitle = "Red means higher off-frame peak fraction in infected samples than matched controls; dot size is infected candidate-run count.",
        x = "Study / cell / inhibitor / growth context",
        y = NULL,
        fill = "Case-control\nfraction delta",
        size = "Case\ncandidate runs"
      ) +
      theme_bw(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                   size = 8),
        plot.title = element_text(face = "bold"),
        legend.position = "right"
      )
    ggsave(file.path(figure_dir, "infected_off_frame_cds_context_matrix.png"),
           p2, width = 13, height = 8.5, dpi = 150)
    if (exists("dominant_save_ggplotly")) {
      dominant_save_ggplotly(
        p2,
        file.path(figure_dir, "infected_off_frame_cds_context_matrix.html"),
        title = "infected_off_frame_cds_context_matrix"
      )
    }
  }
}

top_run_plot_dt <- copy(head(top_runs, 400L))
if (nrow(top_run_plot_dt) > 0) {
  top_plot_genes <- gene_scores[review_tier != "not_prioritized",
                                head(gene_symbol, 14L)]
  top_run_plot_dt[, gene_plot_group := fifelse(
    gene_symbol %chin% top_plot_genes,
    gene_symbol,
    "other"
  )]
  top_run_plot_dt[, gene_plot_group := factor(
    gene_plot_group,
    levels = c(top_plot_genes, "other")
  )]
  p3 <- ggplot(
    top_run_plot_dt,
    aes(
      x = target_minus_housekeeping_bumpiness_z,
      y = top1pct_off_frame_fraction,
      size = pmin(cds_total_counts, 10000),
      color = gene_plot_group,
      text = paste0(
        gene_symbol,
        "<br>Run=", Run,
        "<br>BioProject=", BioProject,
        "<br>cell=", CELL_LINE,
        "<br>condition=", CONDITION,
        "<br>time=", TIMEPOINT,
        "<br>pathogen=", pathogen,
        "<br>counts=", cds_total_counts,
        "<br>top off-frame fraction=",
        round(top1pct_off_frame_fraction, 3),
        "<br>off-frame enrichment=",
        round(top1pct_off_frame_enrichment, 3),
        "<br>target-hk bumpiness z=",
        round(target_minus_housekeeping_bumpiness_z, 3)
      )
    )
  ) +
    geom_hline(yintercept = 0.55, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0.25, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.72) +
    scale_size_continuous(range = c(1.5, 8)) +
    labs(
      title = "Infected Candidate Runs: Off-Frame Peaks vs Bumpiness",
      subtitle = "Upper-right points have both off-frame top peaks and concentrated CDS peaks relative to housekeeping.",
      x = "Target minus housekeeping top-1% bumpiness z",
      y = "Top-1% CDS signal in off-frame positions",
      color = "Gene",
      size = "CDS counts\n(capped)"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")
  ggsave(file.path(figure_dir, "infected_off_frame_cds_top_candidate_runs.png"),
         p3, width = 11, height = 7.5, dpi = 150)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p3,
      file.path(figure_dir, "infected_off_frame_cds_top_candidate_runs.html"),
      title = "infected_off_frame_cds_top_candidate_runs"
    )
  }
}

message("Saved infected off-frame CDS scan outputs to: ", output_dir)
message("Top gene candidates:")
print(head(
  gene_scores[review_tier != "not_prioritized",
              .(gene_symbol, review_tier, candidate_score,
                total_case_n_off_frame_peak,
                n_strong_matched_contexts,
                n_possible_matched_contexts,
                n_unmatched_recurrent_contexts,
                weighted_delta_fraction_off_frame_peak,
                top_context_label)],
  25L
))
