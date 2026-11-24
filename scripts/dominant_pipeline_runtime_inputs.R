# The megafst index is the update marker for the whole coverage collection.
# Resolving filenames is cheap; the runner signatures check them on each step.
pipeline_runtime_input_files <- function() {
  samples <- tryCatch(ORFik::read.experiment("all_samples-Homo_sapiens", validate = FALSE),
                      error = function(e) NULL)
  if (is.null(samples)) return(list())
  ref <- ORFik::refFolder(samples)
  optimized <- ORFik:::optimized_txdb_path(samples, stop.error = FALSE)
  annotations <- c(samples@txdb, samples@fafile,
    if (!is.null(optimized)) unlist(lapply(c("cds", "leaders", "mrna"),
      function(part) paste0(optimized, "_", part, c(".qs", ".rds")))),
    unlist(lapply(unique(c(ref, dirname(samples@fafile))), function(folder) c(
      file.path(folder, "predicted_translons", paste0("predicted_translons_with_sequence_ranges", c(".qs", ".rds"))),
      file.path(folder, "predicted_translons", "predicted_translons_with_sequence.fst"),
      file.path(folder, "predicted_translons", "TransCode", "predicted_translons_ranges.qs")))))
  index <- file.path(ORFik::resFolder(samples), "collection_tables_indexed", "coverage_index.fst")
  counts <- tryCatch(as.character(ORFik:::countTablePath(samples, "cds")),
                     error = function(e) character())
  lib_sizes <- file.path(ORFik::QCfolder(samples), paste0("totalCounts_mrna", c(".qs", ".rds")))
  experiment_file <- file.path(ORFik::config()[["exp"]], "all_samples-Homo_sapiens.csv")
  merged <- tryCatch(ORFik::read.experiment(
    Sys.getenv("RDG_SDRIVE_EXPERIMENT", "all_merged-Homo_sapiens"), validate = FALSE),
    error = function(e) NULL)
  merged_files <- if (is.null(merged)) character() else c(
    ORFik::filepath(merged, "pshifted"),
    unlist(tryCatch(ORFik::filepath(merged, "bigwig"), error = function(e) character())))
  list(annotations = annotations, index = index, counts = counts,
       sample_info = c(experiment_file, lib_sizes), merged = merged_files)
}

pipeline_step_runtime_inputs <- function(step_name, runtime_inputs) {
  coverage_consumers <- c("rdg_graphs", "cds_frame_bumpiness_qc",
    "infected_off_frame_frame_risk", "infected_off_frame_hotspots",
    "infected_off_frame_frame_diagnostics", "downstream_tis_rescue")
  if (exists("pipeline_optional_coverage_consumers", mode = "function")) {
    coverage_consumers <- c(coverage_consumers, pipeline_optional_coverage_consumers(
      get0("analysis_dir", ifnotfound = NULL)))
  }
  if (step_name %in% coverage_consumers) {
    return(c(runtime_inputs$annotations, runtime_inputs$index,
             runtime_inputs$sample_info))
  }
  switch(step_name,
    manual_translons = c(runtime_inputs$annotations, runtime_inputs$merged),
    fst_pages = c(runtime_inputs$annotations, runtime_inputs$index, runtime_inputs$sample_info),
    clean_cds = c(runtime_inputs$annotations, runtime_inputs$index,
                  runtime_inputs$counts, runtime_inputs$sample_info),
    uorf_regulation = c(runtime_inputs$annotations, runtime_inputs$index,
                       runtime_inputs$sample_info),
    character())
}
