# One-time conversion of project-owned files; external Bio_data is not traversed.
source("scripts/run_dominant_cell_state_pipeline.R")
analysis_dir <- find_analysis_dir()
setup_runtime_env(analysis_dir)
load_ribocrypt(find_repo_root(analysis_dir))
source("scripts/dominant_serialization.R")
source("scripts/dominant_object_cache.R")
bench_dir <- file.path(analysis_dir, "results", "qs2_migration_2026-10-09")
backup_dir <- file.path(bench_dir, "legacy_rds")
dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(analysis_dir, pattern = "[.]rds$", recursive = TRUE,
                    full.names = TRUE, all.files = TRUE)
files <- files[!grepl("/([.]git|qs2_migration_2026-10-09)/", files)]
old_coverage_md5 <- unname(tools::md5sum(file.path(bench_dir, "source_before", "scripts", "dominant_coverage_cache.R")))
new_coverage_md5 <- unname(tools::md5sum("scripts/dominant_coverage_cache.R"))
# Cache infrastructure is excluded from scientific calculation signatures.
# Rebase only signatures that used the exact pre-migration coverage helper.
infrastructure <- c("dominant_cached_coverage_by_transcript", "dominant_coverage_cache_file",
  "dominant_coverage_cache_dir", "dominant_safe_name", "dominant_hash_object",
  "dominant_display_signature")
records <- list()
# Resume safely after an interruption: archived sources have already passed
# verification. Recheck them and reconstruct any rows since the last checkpoint.
previous_file <- file.path(bench_dir, "conversion_manifest.csv")
previous <- if (file.exists(previous_file)) read.csv(previous_file) else data.frame()
archived <- list.files(backup_dir, "[.]rds$", recursive=TRUE, full.names=TRUE, all.files=TRUE)
for (input in archived) {
  relative <- substring(input,nchar(backup_dir)+2L)
  match_row <- if (nrow(previous)) which(previous$input==relative & previous$status=="converted") else integer()
  if (length(match_row)) {
    records[[length(records)+1L]] <- previous[match_row[1],,drop=FALSE]
    next
  }
  object <- readRDS(input)
  output <- file.path(analysis_dir,sub("[.]rds$",".qs",relative))
  rebased <- FALSE
  if (grepl("/[.]measurement_cache/",paste0("/",relative)) && is.list(object) &&
      is.list(object$signature) && identical(object$signature$cache_helper,old_coverage_md5)) {
    object$signature$cache_helper <- new_coverage_md5
    object$signature$computation <- object$signature$computation[setdiff(names(object$signature$computation),infrastructure)]
    namespace <- sub("__[[:xdigit:]]+[.]rds$","",basename(input))
    output <- file.path(dirname(output),paste0(namespace,"__",digest::digest(object$signature,algo="md5"),".qs"))
    rebased <- TRUE
  }
  restored <- dominant_read_object(output)
  stopifnot(identical(object,restored) || isTRUE(all.equal(object,restored,check.attributes=TRUE)))
  records[[length(records)+1L]] <- data.frame(input=relative,output=substring(output,nchar(analysis_dir)+2L),
    status="converted",signature_rebased=rebased,bytes_before=file.info(input)$size,
    bytes_after=file.info(output)$size,read_seconds=NA_real_,write_seconds=NA_real_,error="")
}
for (i in seq_along(files)) {
  input <- files[[i]]
  relative <- substring(input, nchar(analysis_dir) + 2L)
  output <- sub("[.]rds$", ".qs", input)
  result <- tryCatch({
    read_seconds <- system.time({object <- readRDS(input)}, gcFirst = FALSE)[["elapsed"]]
    rebased <- FALSE
    if (grepl("/[.]measurement_cache/", input) && is.list(object) &&
        is.list(object$signature) && identical(object$signature$cache_helper, old_coverage_md5)) {
      object$signature$cache_helper <- new_coverage_md5
      object$signature$computation <- object$signature$computation[
        setdiff(names(object$signature$computation), infrastructure)]
      # The reusable value is unchanged; only cache infrastructure metadata moves.
      namespace <- sub("__[[:xdigit:]]+[.]rds$", "", basename(input))
      output <- file.path(dirname(input), paste0(namespace, "__",
        digest::digest(object$signature, algo = "md5"), ".qs"))
      rebased <- TRUE
    }
    if (file.exists(output)) {
      stopifnot(isTRUE(all.equal(object, dominant_read_object(output), check.attributes=TRUE)))
      write_seconds <- 0
    } else {
      write_seconds <- system.time(dominant_save_object(object, output), gcFirst = FALSE)[["elapsed"]]
    }
    restored <- dominant_read_object(output)
    if (!identical(object, restored) && !isTRUE(all.equal(object, restored, check.attributes = TRUE))) {
      stop("QS round-trip mismatch")
    }
    bytes_before <- file.info(input)$size
    bytes_after <- file.info(output)$size
    backup <- file.path(backup_dir, relative)
    dir.create(dirname(backup), recursive = TRUE, showWarnings = FALSE)
    if (!file.rename(input, backup)) stop("Could not archive verified RDS: ", input)
    data.frame(input=relative, output=substring(output,nchar(analysis_dir)+2L),
      status="converted", signature_rebased=rebased, bytes_before=bytes_before,
      bytes_after=bytes_after, read_seconds=read_seconds, write_seconds=write_seconds, error="")
  }, error=function(e) data.frame(input=relative,output=substring(output,nchar(analysis_dir)+2L),
      status="error",signature_rebased=FALSE,bytes_before=file.info(input)$size,
      bytes_after=NA_real_,read_seconds=NA_real_,write_seconds=NA_real_,error=conditionMessage(e)))
  records[[length(records)+1L]] <- result
  if (i %% 100L == 0L || result$status == "error") {
    write.csv(do.call(rbind,records),file.path(bench_dir,"conversion_manifest.csv"),row.names=FALSE)
    message("Converted ",i," / ",length(files),"; latest: ",result$status)
  }
}
manifest <- do.call(rbind, records)
write.csv(manifest,file.path(bench_dir,"conversion_manifest.csv"),row.names=FALSE)
print(table(manifest$status))
cat("Bytes:",sum(manifest$bytes_before),"->",sum(manifest$bytes_after,na.rm=TRUE),"\n")
stopifnot(all(manifest$status=="converted"))
