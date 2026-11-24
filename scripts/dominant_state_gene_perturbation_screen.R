#!/usr/bin/env Rscript
#
# Systematic gene x dominant-state perturbation screen.
#
# For every gene with a knockdown/knockout/shRNA/siRNA/CRISPR perturbation
# label in the compendium's metadata, compares each of the 9 dominant-state
# scores against a matched control group from the SAME study and cell line,
# then aggregates across studies. This formalizes the ad hoc recipe that
# found the METTL3 (viral mimicry), NAT10 (senescence), and VEGFA (uORF
# gatekeeping) leads in md/dominant_state_pipeline_lessons_2026-09-28.md --
# promoted to a real pipeline step per that note's gap #2, instead of being
# rebuilt from scratch each time.
#
# Control matching (gap #3 from the same note): within a (study, cell_line)
# group, prefer a procedure-matched control (scrambled/non-targeting/mock/
# empty-vector -- i.e. something that went through the same transduction/
# selection process as the perturbation) over a generic untransduced "WT"
# baseline when both exist, since using the wrong control silently changes
# the answer (this exact mistake was caught and fixed for METTL3/HEK293
# earlier this session).

suppressPackageStartupMessages(library(data.table))

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (basename(here) == "dominant_cell_states" &&
        dir.exists(file.path(here, "scripts"))) {
      return(here)
    }
    candidate <- file.path(here, "dominant_cell_states")
    if (dir.exists(file.path(candidate, "scripts"))) {
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
results_dir <- Sys.getenv("DOMINANT_CELL_STATES_RESULTS_DIR",
                          unset = file.path(analysis_dir, "results"))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- function(...) file.path(results_dir, ...)

sample_level_output <- result_file("dominant_state_gene_perturbation_screen_sample_level.csv")
gene_state_output <- result_file("dominant_state_gene_perturbation_screen_gene_state_summary.csv")
top_hits_output <- result_file("dominant_state_gene_perturbation_screen_top_hits.csv")

source(file.path(analysis_dir, "scripts", "dominant_state_module_definitions.R"))
source(file.path(analysis_dir, "scripts", "dominant_column_dictionary.R"))

state_file <- result_file("human_dominant_cell_states_clean_cds_gene_expression.csv")
if (!file.exists(state_file)) {
  stop("Missing ", state_file, ". Run dominant_cell_states_fst_clean_cds.R first.")
}
expr <- fread(state_file)
gene_symbols <- expr$gene_symbol
expr_mat <- as.matrix(expr[, -1, with = FALSE])
rownames(expr_mat) <- gene_symbols
storage.mode(expr_mat) <- "numeric"
expr_norm <- log2(expr_mat + 1)

modules <- load_dominant_state_modules()
baseline_genes <- modules[dominant_state == "Baseline (control)", up_genes[[1]]]
baseline_present <- intersect(baseline_genes, rownames(expr_norm))
baseline <- colMeans(expr_norm[baseline_present, , drop = FALSE], na.rm = TRUE)
gene_z <- sweep(expr_norm, 2, baseline, "-")

score_module <- function(gene_z, up_genes, down_genes = character()) {
  up_present <- intersect(up_genes, rownames(gene_z))
  down_present <- intersect(down_genes, rownames(gene_z))
  up_score <- if (length(up_present) == 0) rep(NA_real_, ncol(gene_z)) else
    colMeans(gene_z[up_present, , drop = FALSE], na.rm = TRUE)
  down_score <- if (length(down_genes) == 0) rep(0, ncol(gene_z)) else if (length(down_present) == 0)
    rep(NA_real_, ncol(gene_z)) else colMeans(gene_z[down_present, , drop = FALSE], na.rm = TRUE)
  up_score - down_score
}

state_scores <- list()
for (i in seq_len(nrow(modules))) {
  state <- modules$dominant_state[[i]]
  state_scores[[state]] <- score_module(gene_z, modules$up_genes[[i]], modules$down_genes[[i]])
}
score_dt <- as.data.table(state_scores)
score_dt[, Run := colnames(gene_z)]
state_cols <- modules$dominant_state

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
if (!file.exists(metadata_file)) stop("Missing metadata file: ", metadata_file)
m <- fread(metadata_file)
keep <- intersect(c("Run", "CELL_LINE", "CONDITION", "GENE", "AUTHOR", "study"), names(m))
m <- unique(m[, ..keep], by = "Run")

clean_level <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "NA", "N/A", "na", "n/a", "NULL", "null", "None", "none")] <- NA_character_
  x
}
for (f in intersect(c("CELL_LINE", "CONDITION", "GENE", "AUTHOR", "study"), names(m))) {
  m[, (f) := clean_level(get(f))]
}

dt <- merge(score_dt, m, by = "Run", all.x = TRUE, sort = FALSE)

exclude_cell_lines <- c("LCL")  # constitutive EBV-driven confound; see
                                # md/lcl_ebv_interferon_translation_note_2026-09-25.md
dt <- dt[!(CELL_LINE %chin% exclude_cell_lines)]

# Two-tier control classification: prefer a procedure-matched control
# (scrambled/non-targeting/mock/empty-vector) within the same study+cell
# line; fall back to a generic WT/control baseline only if no tighter
# control exists in that group.
tight_control_pattern <- "^(scram|scr|scrambled|ntc|non.?targeting|mock|ev|empty.?vector)$"
loose_control_pattern <- "^(wt|control|ctrl|parental)$"
perturb_pattern <- "kd|ko|knockdown|knockout|shrna|sirna|crispr"

dt[, cond_lower := tolower(trimws(CONDITION))]
dt[, is_tight_control := grepl(tight_control_pattern, cond_lower)]
dt[, is_loose_control := grepl(loose_control_pattern, cond_lower)]
dt[, is_perturb := !is.na(GENE) & nzchar(GENE) & grepl(perturb_pattern, cond_lower)]

group_has_tight <- dt[, .(has_tight = any(is_tight_control)), by = .(study, CELL_LINE)]
dt <- merge(dt, group_has_tight, by = c("study", "CELL_LINE"), all.x = TRUE, sort = FALSE)
dt[, is_control := is_tight_control | (is_loose_control & !has_tight)]

pert <- dt[is_perturb == TRUE]
ctrl <- dt[is_control == TRUE]

message("Gene perturbation screen: ", nrow(pert), " candidate perturbation samples, ",
        nrow(ctrl), " candidate control samples, across ",
        uniqueN(dt$study), " studies.")

ctrl_by_group <- ctrl[, c(list(n_ctrl = .N), lapply(.SD, mean, na.rm = TRUE)),
                      by = .(study, CELL_LINE), .SDcols = state_cols]
setnames(ctrl_by_group, state_cols, paste0("ctrl_mean_", make.names(state_cols)))

pert2 <- merge(pert, ctrl_by_group, by = c("study", "CELL_LINE"), all.x = TRUE, sort = FALSE)
pert2 <- pert2[!is.na(n_ctrl) & n_ctrl >= 2]

sample_rows <- rbindlist(lapply(state_cols, function(state_name) {
  ctrl_col <- paste0("ctrl_mean_", make.names(state_name))
  data.table(
    Run = pert2$Run,
    study = pert2$study,
    CELL_LINE = pert2$CELL_LINE,
    GENE = pert2$GENE,
    CONDITION = pert2$CONDITION,
    AUTHOR = pert2$AUTHOR,
    n_ctrl = pert2$n_ctrl,
    state = state_name,
    sample_score = pert2[[state_name]],
    ctrl_mean = pert2[[ctrl_col]]
  )
}), fill = TRUE)
sample_rows[, delta := sample_score - ctrl_mean]
sample_rows <- sample_rows[is.finite(delta)]
fwrite(sample_rows, sample_level_output)
message("Saved: ", sample_level_output)

safe_t_test_p <- function(pert_vals, ctrl_vals) {
  pert_vals <- pert_vals[is.finite(pert_vals)]
  ctrl_vals <- ctrl_vals[is.finite(ctrl_vals)]
  if (length(pert_vals) < 2 || length(ctrl_vals) < 2) return(NA_real_)
  if (length(unique(c(pert_vals, ctrl_vals))) < 2) return(NA_real_)
  tryCatch(t.test(pert_vals, ctrl_vals)$p.value, error = function(e) NA_real_)
}

gene_state_summary <- sample_rows[, {
  by_study <- .SD[, .(
    study_delta = mean(delta, na.rm = TRUE),
    n_pert_in_study = .N
  ), by = study]
  list(
    n_samples = .N,
    n_studies = uniqueN(study),
    n_cell_lines = uniqueN(CELL_LINE),
    mean_delta = mean(delta, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    min_delta = min(delta, na.rm = TRUE),
    max_delta = max(delta, na.rm = TRUE),
    consistent_direction = all(sign(by_study$study_delta) == sign(by_study$study_delta[1])) &&
      by_study$study_delta[1] != 0,
    t_test_p = safe_t_test_p(sample_score, ctrl_mean),
    best_study = study[which.max(abs(delta))],
    best_author = AUTHOR[which.max(abs(delta))]
  )
}, by = .(GENE, state)]
setnames(gene_state_summary, "GENE", "gene_symbol")
gene_state_summary[, p_adj := p.adjust(t_test_p, method = "BH")]
gene_state_summary <- gene_state_summary[order(-abs(mean_delta))]
fwrite(gene_state_summary, gene_state_output)
message("Saved: ", gene_state_output)

top_hits <- gene_state_summary[
  n_studies >= 2 & consistent_direction == TRUE
][order(-abs(mean_delta))]
fwrite(top_hits, top_hits_output)
message("Saved: ", top_hits_output)

message("Top hits (>=2 independent studies, consistent direction):")
if (nrow(top_hits) > 0) {
  print(top_hits[1:min(15, .N), .(gene_symbol, state, n_studies, n_cell_lines, mean_delta, t_test_p)])
} else {
  message("  (none found)")
}

shared_col_descriptions <- c(
  Run = "Sample accession for the perturbed sample.",
  study = "Study/BioProject accession this sample came from.",
  CELL_LINE = "Cell line the perturbation was performed in.",
  GENE = "Gene targeted by the knockdown/knockout/shRNA/siRNA/CRISPR (renamed to `gene_symbol` in the aggregated tables).",
  gene_symbol = "Gene targeted by the perturbation.",
  CONDITION = "Free-text condition label from the source metadata (e.g. 'KD', 'KO', 'shRNA').",
  AUTHOR = "First author of the source study.",
  n_ctrl = "Number of matched control samples found in the same (study, cell line) group.",
  state = "Which of the 9 dominant states (or Baseline) this row's delta refers to.",
  sample_score = "This sample's score for `state` (per-sample gene_z-derived module score).",
  ctrl_mean = "Mean of the matched control group's score for `state`, in the same (study, cell line) group.",
  delta = "sample_score - ctrl_mean. Positive = state score higher in the perturbed sample than its matched control.",
  n_samples = "Number of perturbed samples contributing to this gene x state row.",
  n_studies = "Number of independent studies (BioProjects) contributing.",
  n_cell_lines = "Number of distinct cell lines contributing.",
  mean_delta = "Mean of per-sample deltas across all contributing samples (log2-scale module-score units).",
  median_delta = "Median of per-sample deltas.",
  min_delta = "Smallest per-sample delta (checks whether the direction is consistent, not just the average).",
  max_delta = "Largest per-sample delta.",
  consistent_direction = "TRUE if every contributing study's own mean delta has the same sign (not just the grand mean).",
  t_test_p = "Welch two-sample t-test p-value, sample_score vs. ctrl_mean pooled across all contributing samples. NA if either arm has <2 samples.",
  p_adj = "Benjamini-Hochberg adjusted t_test_p across all gene x state rows in this table.",
  best_study = "The single study contributing this gene x state's largest-magnitude delta.",
  best_author = "First author of best_study."
)

write_column_dictionary(
  sample_level_output, shared_col_descriptions,
  title = "Per-sample perturbation-vs-matched-control deltas, one row per (perturbed sample, dominant state)."
)
write_column_dictionary(
  gene_state_output, shared_col_descriptions,
  title = "Gene x dominant-state summary: matched-control delta aggregated across all contributing studies."
)
write_column_dictionary(
  top_hits_output, shared_col_descriptions,
  title = "Same columns as the gene x state summary, filtered to n_studies >= 2 and consistent_direction == TRUE."
)
