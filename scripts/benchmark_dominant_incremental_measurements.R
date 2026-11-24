# Run after a normal uorf_regulation pass has populated measurement caches.
# Exercises real per-gene computations in an isolated cache; never writes
# production feature tables or alters the curated gene panel.
source("scripts/run_dominant_cell_state_pipeline.R")
analysis_dir <- find_analysis_dir()
setup_runtime_env(analysis_dir)
load_ribocrypt(find_repo_root(analysis_dir))
source(file.path(analysis_dir, "scripts", "dominant_coverage_cache.R"))
source(file.path(analysis_dir, "scripts", "dominant_object_cache.R"))
source(file.path(analysis_dir, "scripts", "dominant_transcript_mapping.R"))
source(file.path(analysis_dir, "scripts", "dominant_internal_orf_signal_helpers.R"))
# Load definitions only, without executing the producer or rewriting its outputs.
producer <- file.path(analysis_dir, "scripts", "dominant_uorf_regulation_inference.R")
for (expr in parse(producer)) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) {
    eval(expr, envir = globalenv())
  }
}
bench_dir <- Sys.getenv("DOMINANT_BENCHMARK_DIR", unset =
  file.path(analysis_dir, "results", paste0("performance_", Sys.Date())))
dir.create(bench_dir, recursive = TRUE, showWarnings = FALSE)
cache_files <- list.files(file.path(analysis_dir, "results", ".measurement_cache"),
                          pattern = "^uorf_.*\\.qs$", full.names = TRUE)
computation <- dominant_computation_signature(measure_uorf_features)
objects <- lapply(cache_files, dominant_read_object)
objects <- Filter(function(x) !is.null(x$value) &&
  identical(x$signature$computation, computation), objects)
if (!length(objects)) stop("Populate current uORF measurement caches first")
genes <- vapply(objects, function(x) x$signature$gene, character(1))
objects <- objects[!duplicated(genes)]
held_out <- which(vapply(objects, function(x) x$signature$gene == "FMR1", logical(1)))
if (length(held_out) != 1L) stop("Expected a current FMR1 cache for the held-out-gene benchmark")
df <- ORFik::read.experiment("all_samples-Homo_sapiens", validate = FALSE)
tx_ids <- vapply(objects, function(x) x$signature$tx_id, character(1))
cds_all <- loadRegion(df, "cds", names.keep = tx_ids)
leader_all <- loadRegion(df, "leaders", names.keep = tx_ids)
fst_index <- objects[[1]]$signature$coverage$fst_index
isolated <- tempfile("incremental-benchmark-", tmpdir = bench_dir)
dir.create(isolated)
compute_calls <- 0L
reference_check <- TRUE
measure <- function(i) {
  object <- objects[[i]]
  signature <- object$signature
  value <- dominant_cached_object(isolated, paste0("uorf_", signature$gene), signature,
    function() {
      compute_calls <<- compute_calls + 1L
      measure_uorf_features(df,
        data.table::data.table(gene_symbol = signature$gene, tx_id = signature$tx_id),
        signature$translons, fst_index, signature$lib_sizes, cds_all, leader_all)
    })
  if (reference_check) stopifnot(isTRUE(all.equal(value, object$value, check.attributes = TRUE)))
  invisible(value)
}
results <- list()
run_case <- function(label, indices) {
  before <- compute_calls
  elapsed <- system.time(lapply(indices, measure))[["elapsed"]]
  results[[length(results) + 1L]] <<- data.frame(
    scenario = label, genes_requested = length(indices),
    genes_computed = compute_calls - before, seconds = elapsed)
}
run_case("initial_panel_without_FMR1", setdiff(seq_along(objects), held_out))
run_case("add_FMR1_to_existing_panel", seq_along(objects))
run_case("unchanged_panel", seq_along(objects))
# Edit a range in memory only; production curated inputs stay untouched.
original <- objects[[held_out]]$signature$translons
edited <- original
first <- edited[[1]]
start(first)[1] <- start(first)[1] + 3L
edited[[1]] <- first
objects[[held_out]]$signature$translons <- edited
reference_check <- FALSE
run_case("simulated_FMR1_range_edit", seq_along(objects))
reference_check <- TRUE
objects[[held_out]]$signature$translons <- original
# An index timestamp change must invalidate every gene, without touching real data.
for (i in seq_along(objects)) {
  objects[[i]]$signature$coverage$fst_index_mtime <-
    objects[[i]]$signature$coverage$fst_index_mtime + 1
}
run_case("simulated_index_rebuild", seq_along(objects))
stopifnot(results[[2]]$genes_computed == 1L,
          results[[3]]$genes_computed == 0L,
          results[[4]]$genes_computed == 1L,
          results[[5]]$genes_computed == length(objects))
data.table::fwrite(data.table::rbindlist(results), file.path(bench_dir, "incremental_gene_benchmark.csv"))
print(data.table::rbindlist(results))
