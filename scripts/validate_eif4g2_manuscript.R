#!/usr/bin/env Rscript

# Fail fast when manuscript formatting or the central numerical claims drift
# away from the saved evidence tables. This is intentionally a focused audit,
# not a substitute for rebuilding every source-data panel.

suppressPackageStartupMessages(library(data.table))

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(here, "DESCRIPTION")) &&
        dir.exists(file.path(here, "manuscript"))) {
      return(here)
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find the dominant_cell_states analysis root.", call. = FALSE)
}

analysis_dir <- find_analysis_dir()
article_file <- file.path(analysis_dir, "manuscript", "eif4g2_article_draft.md")
figure_file <- file.path(
  analysis_dir, "manuscript", "eif4g2_figure_plan_and_legends.md"
)
output_dir <- file.path(analysis_dir, "manuscript", "generated")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

article <- readLines(article_file, warn = FALSE, encoding = "UTF-8")
figures <- readLines(figure_file, warn = FALSE, encoding = "UTF-8")
article_text <- gsub("[[:space:]]+", " ", paste(article, collapse = " "))

markdown_words <- function(lines) {
  text <- paste(lines, collapse = " ")
  text <- gsub("https?://[^[:space:]]+", " ", text)
  text <- gsub("[`*_>#|]", " ", text)
  text <- gsub("\\[[^]]*\\]\\([^)]*\\)", " ", text)
  tokens <- strsplit(trimws(text), "[[:space:]]+")[[1L]]
  sum(nzchar(tokens))
}

section_lines <- function(lines, heading, next_heading_pattern = "^## ") {
  start <- match(heading, lines)
  if (is.na(start)) stop("Missing manuscript heading: ", heading, call. = FALSE)
  later <- which(seq_along(lines) > start & grepl(next_heading_pattern, lines))
  end <- if (length(later)) min(later) - 1L else length(lines)
  lines[seq.int(start + 1L, end)]
}

title_line <- article[grepl("^# ", article)][1L]
title <- sub("^# ", "", title_line)
summary_words <- markdown_words(section_lines(article, "## Summary"))
main_start <- match("## A discrete translation branch in the EIF4G2 leader",
                    article)
main_end <- match("## Main references", article) - 1L
if (is.na(main_start) || is.na(main_end) || main_end <= main_start) {
  stop("Could not identify the main-text boundaries.", call. = FALSE)
}
main_words <- markdown_words(article[seq.int(main_start + 1L, main_end)])
reference_count <- sum(grepl("^[0-9]+\\. ", article))

errors <- character()
warnings <- character()
require_true <- function(value, message) {
  if (!isTRUE(value)) errors <<- c(errors, message)
}
warn_unless <- function(value, message) {
  if (!isTRUE(value)) warnings <<- c(warnings, message)
}

require_true(nchar(title, type = "chars") <= 75L,
             "Title exceeds Nature's 75-character working limit.")
require_true(summary_words <= 200L,
             "Summary exceeds the 200-word working limit.")
require_true(main_words <= 4300L,
             "Main text exceeds the 4,300-word expanded Article guide.")
warn_unless(main_words <= 2500L,
            "Main text exceeds the typical six-page 2,500-word guide.")
require_true(reference_count <= 50L,
             "Main reference list exceeds the 50-reference working guide.")

required_phrases <- c(
  "Working-draft status",
  "Collapsing those lanes yielded two tRNA-Glu and two control biological units",
  "descriptive outlier measures rather than z-scores",
  "Direct action at the two GAG codons therefore remains a reporter hypothesis",
  "data do not support a general metastatic-state switch",
  "reporter and immunopeptidomic tests are required to establish causality"
)
for (phrase in required_phrases) {
  require_true(grepl(phrase, article_text, fixed = TRUE),
               paste0("Required claim-boundary text is missing: ", phrase))
}

forbidden_phrases <- c(
  "tRNA-Glu directly decodes the EIF4G2 uORF",
  "the GAG codons sense tRNA-Glu abundance",
  "tRNA-Glu increases DAP5 expression",
  "tRNA-Glu exposes an EIF4G2 neoantigen",
  "tRNA-Glu increases an EIF4G2 neoantigen",
  "a universal collision/recycling sensor",
  "a metastasis-specific switch",
  "replicated six times"
)
for (phrase in forbidden_phrases) {
  require_true(!grepl(phrase, article_text, fixed = TRUE),
               paste0("Unsupported phrase found in article: ", phrase))
}

read_required <- function(relative_path) {
  path <- file.path(analysis_dir, relative_path)
  require_true(file.exists(path), paste0("Missing evidence table: ", relative_path))
  if (!file.exists(path)) return(NULL)
  fread(path)
}

assert_close <- function(observed, expected, label, tolerance = 5e-7) {
  require_true(
    length(observed) == 1L && is.finite(observed) &&
      abs(observed - expected) <= tolerance,
    paste0(label, " drifted: observed ", observed, ", expected ", expected)
  )
}

evidence <- read_required(file.path(
  "results", "dominant_rdg_eif4g2_evidence_packet",
  "eif4g2_contrast_evidence.csv"
))
if (!is.null(evidence)) {
  glu <- evidence[comparison == "tRNA-Glu(UUC) OE"]
  arg <- evidence[comparison == "tRNA-Arg(CCG) OE\nnegative control"]
  abce1 <- evidence[comparison == "ABCE1 KO"]
  assert_close(glu$union_clean_cds_allocation_delta, -0.137059612722506,
               "tRNA-Glu allocation delta")
  assert_close(glu$leave_one_uorf_position_min_abs_delta, 0.0887594823113436,
               "tRNA-Glu single-position bound")
  assert_close(glu$pooled_top3_removed_allocation_delta, -0.0513123055020609,
               "tRNA-Glu top-three-position bound")
  assert_close(arg$union_clean_cds_allocation_delta, 0.0183805139729601,
               "tRNA-Arg allocation delta")
  assert_close(abce1$union_clean_cds_allocation_delta, -0.196742946939377,
               "ABCE1 allocation delta")
}

trna_raw <- read_required(file.path(
  "results", "dominant_rdg_eif4g2_trna_length_resolved",
  "eif4g2_length_contrast_summary.csv"
))
trna_raw_specificity <- read_required(file.path(
  "results", "dominant_rdg_eif4g2_trna_length_resolved",
  "eif4g2_length_specificity_summary.csv"
))
if (!is.null(trna_raw) && !is.null(trna_raw_specificity)) {
  decisive_categories <- c(
    "all_lane_matched_well_phased",
    "all_study_shared_well_phased",
    "canonical_28_32nt_lane_matched_well_phased"
  )
  expected_glu <- c(-0.135631250110818, -0.14146478117896,
                    -0.144531582362685)
  glu_raw <- trna_raw[
    perturbation == "tRNA-Glu" & category %chin% decisive_categories
  ][match(decisive_categories, category)]
  require_true(nrow(glu_raw) == 3L && all(glu_raw$evidence_gate == TRUE),
               "One or more decisive raw tRNA-Glu strata no longer pass.")
  for (i in seq_along(decisive_categories)) {
    assert_close(
      glu_raw$clean_cds_allocation_delta[[i]], expected_glu[[i]],
      paste0("Raw tRNA-Glu ", decisive_categories[[i]], " allocation delta")
    )
  }
  main_glu <- glu_raw[category == "all_lane_matched_well_phased"]
  assert_close(main_glu$clean_cds_vs_uorf1_allocation_delta, -0.134554208571789,
               "Raw tRNA-Glu exact-uORF1 allocation delta")
  require_true(
    main_glu$case_uorf1_counts == 1028L &&
      main_glu$control_uorf1_counts == 1870L &&
      main_glu$target_uorf_phase_gate == TRUE,
    "Raw tRNA-Glu exact-uORF1 counts or phase gate drifted."
  )
  main_specificity <- trna_raw_specificity[
    category == "all_lane_matched_well_phased"
  ]
  assert_close(main_specificity$arg_clean_cds_allocation_delta,
               0.0118184288756015, "Raw tRNA-Arg allocation delta")
  assert_close(main_specificity$glu_minus_arg_clean_cds_allocation_delta,
               -0.14744967898642, "Raw tRNA specificity contrast")
  require_true(nrow(main_specificity) == 1L &&
                 main_specificity$specificity_gate == TRUE,
               "Raw tRNA specificity gate no longer passes.")
  short_glu <- trna_raw[
    perturbation == "tRNA-Glu" &
      category == "short_20_23nt_lane_matched_well_phased"
  ]
  require_true(nrow(short_glu) == 1L &&
                 short_glu$measurement_gate == FALSE &&
                 short_glu$evidence_gate == FALSE,
               "Sparse raw tRNA short-read stratum changed gate status.")
}

rna <- read_required(file.path(
  "results", "dominant_rdg_eif4g2_goodarzi_rna_audit",
  "eif4g2_isoform_rna_audit_summary.csv"
))
if (!is.null(rna)) {
  primary <- rna[tracking_id == "NM_001418"]
  assert_close(primary$total_transcript_log2_glu_vs_control,
               -0.0375986432124639, "NM_001418 total-RNA log2 ratio")
  assert_close(primary$total_transcript_isoform_fraction_delta,
               -0.0254674752217013, "NM_001418 total-RNA isoform delta")
}

specificity <- read_required(file.path(
  "results", "dominant_rdg_trna_codon_branch_audit",
  "eif4g2_specificity_rank_sensitivity.csv"
))
if (!is.null(specificity)) {
  require_true(
    nrow(specificity) == 3L &&
      identical(specificity$minimum_uorf_counts_each_arm, c(5L, 10L, 20L)) &&
      all(specificity$eif4g2_glu_specific_rank == 1L) &&
      identical(specificity$n_eligible_features, c(35L, 27L, 18L)),
    "EIF4G2 feature-wide specificity ranks or eligible counts drifted."
  )
}

uv <- read_required(file.path(
  "results", "dominant_rdg_eif4g2_uv_length_resolved",
  "eif4g2_length_contrast_summary.csv"
))
if (!is.null(uv)) {
  expected_uv <- data.table(
    category = c(
      "all_shared_well_phased",
      "canonical_28_32nt_shared_well_phased",
      "short_20_23nt_shared_well_phased"
    ),
    expected = c(-0.3433283, -0.3410946, -0.2741888)
  )
  checked <- merge(expected_uv, uv, by = "category", all.x = TRUE)
  require_true(nrow(checked) == 3L && all(checked$evidence_gate == TRUE),
               "One or more decisive UV length strata no longer pass.")
  for (i in seq_len(nrow(checked))) {
    assert_close(
      checked$clean_cds_allocation_delta[[i]], checked$expected[[i]],
      paste0("UV ", checked$category[[i]], " allocation delta"),
      tolerance = 5e-7
    )
  }
}

legend_starts <- grep("^### Draft legend", figures)
legend_word_counts <- integer(length(legend_starts))
for (i in seq_along(legend_starts)) {
  start <- legend_starts[[i]] + 1L
  later <- which(seq_along(figures) > legend_starts[[i]] &
                   grepl("^## Figure ", figures))
  end <- if (length(later)) min(later) - 1L else length(figures)
  legend_word_counts[[i]] <- markdown_words(figures[seq.int(start, end)])
}
require_true(all(legend_word_counts <= 250L),
             "At least one main figure legend exceeds 250 words.")

placeholder_count <- sum(grepl("\\[[^]]*(required|author|affiliation|email)",
                               article, ignore.case = TRUE))
warn_unless(placeholder_count == 0L,
            paste0(placeholder_count,
                   " author/submission placeholder line(s) remain, as expected."))

report <- c(
  "EIF4G2 manuscript audit",
  paste0("title_characters: ", nchar(title, type = "chars")),
  paste0("summary_words: ", summary_words),
  paste0("main_text_words: ", main_words),
  paste0("main_references: ", reference_count),
  paste0("figure_legend_words: ", paste(legend_word_counts, collapse = ", ")),
  paste0("errors: ", length(errors)),
  paste0("warnings: ", length(warnings)),
  if (length(warnings)) paste0("WARNING: ", warnings) else "WARNING: none",
  if (length(errors)) paste0("ERROR: ", errors) else "ERROR: none"
)
writeLines(report, file.path(output_dir, "eif4g2_manuscript_audit.txt"))
cat(paste(report, collapse = "\n"), "\n")
if (length(errors)) quit(save = "no", status = 1L)
