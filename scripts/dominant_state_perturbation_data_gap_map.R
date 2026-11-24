#!/usr/bin/env Rscript
#
# Data-gap map: which biologically interesting gene classes currently have
# thin (0-1 independent study) perturbation coverage in this compendium.
#
# Implements gap #5 from md/dominant_state_pipeline_lessons_2026-09-28.md:
# the systematic gene x state screen (dominant_state_gene_perturbation_screen.R)
# found only 9 genes anywhere in the whole compendium with >=2 independent
# supporting studies -- true cross-lab replication is the exception, not the
# norm. Re-analyzing what's already in hand cannot fix that; this script
# instead reports which curated, biologically-motivated gene classes are
# thinnest, so future data acquisition can be targeted rather than
# opportunistic. This is a "what to go get" report, not a new statistical
# claim about any gene.

suppressPackageStartupMessages(library(data.table))

source(file.path(Sys.getenv("DOMINANT_CELL_STATES_DIR", "."), "scripts",
                 "dominant_script_common.R"), local = TRUE)

analysis_dir <- find_analysis_dir()
results_dir <- Sys.getenv("DOMINANT_CELL_STATES_RESULTS_DIR",
                          unset = file.path(analysis_dir, "results"))
result_file <- function(...) file.path(results_dir, ...)
source(file.path(analysis_dir, "scripts", "dominant_column_dictionary.R"))

output_file <- result_file("dominant_state_perturbation_data_gap_map.csv")
summary_output <- result_file("dominant_state_perturbation_data_gap_map_summary.txt")

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
keep <- intersect(c("Run", "CELL_LINE", "CONDITION", "GENE", "study"), names(m))
m <- unique(m[, ..keep], by = "Run")

clean_level <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "NA", "N/A", "na", "n/a", "NULL", "null", "None", "none")] <- NA_character_
  x
}
for (f in intersect(c("CELL_LINE", "CONDITION", "GENE", "study"), names(m))) {
  m[, (f) := clean_level(get(f))]
}
m[, cond_lower := tolower(trimws(CONDITION))]
perturb_pattern <- "kd|ko|knockdown|knockout|shrna|sirna|crispr"
m[, is_perturb := !is.na(GENE) & nzchar(GENE) & grepl(perturb_pattern, cond_lower)]

pert <- m[is_perturb == TRUE]
per_gene <- pert[, .(n_samples = .N, n_studies = uniqueN(study),
                     n_cell_lines = uniqueN(CELL_LINE)), by = GENE]
setnames(per_gene, "GENE", "gene_symbol")

# Curated gene classes: biologically motivated groups this project has
# repeatedly cared about this session (epitranscriptomic writers/readers/
# erasers, core translation factors, RQC/stress-granule genes, and the
# genome-guardian/viral-mimicry class checked but not found well-represented
# earlier). Extend this list as new classes become interesting -- it is a
# living reference, not a fixed taxonomy.
gene_classes <- list(
  m6a_writer = c("METTL3", "METTL14", "WTAP", "VIRMA", "RBM15", "RBM15B", "ZC3H13"),
  m6a_reader = c("YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", "YTHDC2", "IGF2BP1", "IGF2BP2", "IGF2BP3"),
  m6a_eraser = c("FTO", "ALKBH5"),
  ac4c_pathway = c("NAT10"),
  ribosome_quality_control = c("ZNF598", "GIGYF2", "EIF2A", "PELO", "HBS1L", "ASCC3", "NEMF", "LTN1", "GCN1"),
  stress_granule = c("G3BP1", "G3BP2", "CAPRIN1", "UBAP2L", "DDX3X"),
  isr_core = c("EIF2AK1", "EIF2AK2", "EIF2AK3", "EIF2AK4", "EIF2S1", "PPP1R15A", "PPP1R15B", "ATF4"),
  retroelement_silencing = c("ADAR", "SETDB1", "SETDB2", "TRIM28", "MPHOSPH8", "PPHLN1", "TASOR",
                            "DNMT1", "DNMT3A", "DNMT3B", "TREX1", "SAMHD1"),
  recycling_reinitiation = c("ABCE1", "EIF2D", "MCTS1", "DENR"),
  translation_termination = c("ETF1", "GSPT1", "GSPT2")
)

class_dt <- rbindlist(lapply(names(gene_classes), function(cls) {
  data.table(gene_class = cls, gene_symbol = gene_classes[[cls]])
}))
class_dt <- merge(class_dt, per_gene, by = "gene_symbol", all.x = TRUE)
class_dt[is.na(n_samples), `:=`(n_samples = 0L, n_studies = 0L, n_cell_lines = 0L)]
class_dt[, coverage_tier := fifelse(
  n_studies == 0, "no_coverage",
  fifelse(n_studies == 1, "single_study_only", "replicated")
)]
setorder(class_dt, gene_class, -n_studies)
fwrite(class_dt, output_file)
message("Saved: ", output_file)

class_summary <- class_dt[, .(
  n_genes = .N,
  n_no_coverage = sum(coverage_tier == "no_coverage"),
  n_single_study = sum(coverage_tier == "single_study_only"),
  n_replicated = sum(coverage_tier == "replicated"),
  genes_with_no_coverage = paste(gene_symbol[coverage_tier == "no_coverage"], collapse = ", ")
), by = gene_class]
setorder(class_summary, -n_no_coverage, -n_single_study)

summary_lines <- c(
  sprintf("Perturbation data-gap map (%s)", Sys.Date()),
  sprintf("Gene classes checked: %d, total genes: %d", nrow(class_summary), nrow(class_dt)),
  "",
  "Classes ranked by how thin their coverage is (most genes with zero",
  "perturbation data in this compendium first -- these are where adding a",
  "new dataset would most efficiently open up a new investigation, the",
  "same way METTL3/NAT10 were found for the m6A/ac4c classes):",
  ""
)
for (i in seq_len(nrow(class_summary))) {
  r <- class_summary[i]
  summary_lines <- c(summary_lines, sprintf(
    "%-28s %d genes | %d no coverage | %d single-study | %d replicated%s",
    r$gene_class, r$n_genes, r$n_no_coverage, r$n_single_study, r$n_replicated,
    if (nzchar(r$genes_with_no_coverage)) sprintf(" | missing: %s", r$genes_with_no_coverage) else ""
  ))
}
writeLines(summary_lines, summary_output)
message("Saved: ", summary_output)
message(paste(summary_lines, collapse = "\n"))

write_column_dictionary(
  output_file,
  c(
    gene_class = "Curated, biologically-motivated gene-class label (see gene_classes list in this script) -- extend that list to check new classes.",
    gene_symbol = "Gene in this class.",
    n_samples = "Number of perturbation samples for this gene anywhere in the compendium (0 if never perturbed).",
    n_studies = "Number of independent studies (BioProjects) perturbing this gene.",
    n_cell_lines = "Number of distinct cell lines this gene has been perturbed in.",
    coverage_tier = "no_coverage (0 studies) / single_study_only (1 study) / replicated (>=2 studies)."
  ),
  title = "Perturbation coverage for curated, biologically-motivated gene classes -- a 'what to go acquire next' map, not a statistical claim about any gene."
)
