# Hold independent copies of the original native functions across package
# reloads: an old .Call wrapper alone would resolve the newly installed code.
load_original_mapping_native <- function(directory) {
  root <- file.path(directory, "baseline_native")
  paths <- file.path(root, c("pmapToTranscripts.cpp", "pmapFromTranscripts.cpp"))
  stopifnot(all(file.exists(paths)))
  env <- new.env(parent = globalenv())
  for (i in seq_along(paths)) {
    text <- paste(readLines(paths[i]), collapse = "\n")
    text <- gsub("pmap", "baseline_pmap", text, fixed = TRUE)
    text <- gsub("nExonsCumSums", "baseline_nExonsCumSums", text, fixed = TRUE)
    compiled <- file.path(root, paste0("independent_", basename(paths[i])))
    if (!file.exists(compiled) || !identical(paste(readLines(compiled), collapse = "\n"), text))
      writeLines(text, compiled)
    Rcpp::sourceCpp(compiled, env = env,
      cacheDir = file.path(root, Sys.getenv("ORFIK_MAPPING_BASELINE_CACHE", "rcpp")))
  }
  list(to = env$baseline_pmapToTranscriptsCPP, from = env$baseline_pmapFromTranscriptsCPP)
}
