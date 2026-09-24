#!/usr/bin/env Rscript
#
# Reuse massiveNGSpipe's already fetched, adapter/barcode-trimmed, collapsed
# (deduplicated, ORFik "_x<count>" convention) reads as the input for the
# EIF4G2 length-resolved reprocessing (Goodarzi tRNA-Glu/tRNA-Arg RPFs,
# GSE141459 UV RPFs), instead of re-downloading FASTQs from EBI and
# re-trimming with skewer.
#
# Confirmed compatible by direct inspection on this server: massiveNGSpipe
# trimmed these exact runs with the same adapter sequence the skewer-based
# scripts expect, and the resulting collapsed reads already sit mostly
# within the target RPF length window. The one real gap: massiveNGSpipe's
# fastp-based trim does not enforce a hard max-length cutoff, so reads with
# no detected adapter are kept at full raw read length (e.g. 51nt) instead
# of being dropped. skewer's original invocation used "-l 15 -L 40", so that
# same window is enforced here before alignment.
#
# Usage (env-var driven, one call per study accession):
#   REUSE_RUNS="SRR3131003 SRR3131004 ..." \
#   REUSE_STUDY_ACCESSION="PRJNA310107-homo_sapiens" \
#   REUSE_OUTPUT_DIR="/path/to/output" \
#   Rscript scripts/dominant_rdg_reuse_massivengspipe_length_resolved.R

suppressPackageStartupMessages(library(data.table))

min_len <- as.integer(Sys.getenv("REUSE_MIN_LEN", "15"))
max_len <- as.integer(Sys.getenv("REUSE_MAX_LEN", "40"))
threads <- Sys.getenv("REUSE_THREADS", "8")

runs <- strsplit(trimws(Sys.getenv("REUSE_RUNS", "")), "\\s+")[[1]]
runs <- runs[nzchar(runs)]
if (length(runs) == 0) stop("REUSE_RUNS must list one or more SRR run ids.")

study_accession <- Sys.getenv("REUSE_STUDY_ACCESSION", unset = "")
if (!nzchar(study_accession)) {
  stop("REUSE_STUDY_ACCESSION must name the massiveNGSpipe processed_data ",
       "folder for these runs, e.g. 'PRJNA310107-homo_sapiens'.")
}

bam_root <- tryCatch(path.expand(unname(ORFik::config()["bam"])),
                     error = function(e) stop("Could not resolve ORFik ",
                                              "processed_data root: ",
                                              conditionMessage(e)))
collapsed_dir <- file.path(bam_root, study_accession, "trim", "SINGLE")
if (!dir.exists(collapsed_dir)) {
  stop("massiveNGSpipe collapsed-reads directory not found: ", collapsed_dir)
}

output_dir <- Sys.getenv("REUSE_OUTPUT_DIR",
                         unset = file.path("results",
                                          "reuse_length_resolved",
                                          study_accession))
dir.create(file.path(output_dir, "filtered"), recursive = TRUE,
          showWarnings = FALSE)
dir.create(file.path(output_dir, "aligned"), recursive = TRUE,
          showWarnings = FALSE)

star_bin <- Sys.getenv(
  "STAR_BIN",
  unset = path.expand("~/bin/STAR-2.7.4a/bin/Linux_x86_64/STAR")
)
star_index <- Sys.getenv(
  "STAR_INDEX",
  unset = path.expand("~/livemount/Bio_data/references/homo_sapiens/STAR_index/genomeDir")
)
if (!file.exists(star_bin)) stop("STAR executable not found: ", star_bin)
if (!dir.exists(star_index)) stop("STAR index not found: ", star_index)

# The massiveNGSpipe collapsed FASTA is a plain two-line-per-record format:
# ">seq<i>_x<count>" followed by the sequence, already deduplicated. Reading
# it as plain text avoids a Biostrings dependency for this simple format.
filter_collapsed_fasta <- function(input_path, output_path, min_len, max_len) {
  lines <- readLines(gzfile(input_path))
  headers <- lines[c(TRUE, FALSE)]
  seqs <- lines[c(FALSE, TRUE)]
  stopifnot(length(headers) == length(seqs))
  widths <- nchar(seqs)
  keep <- widths >= min_len & widths <= max_len
  con <- gzfile(output_path, "w")
  on.exit(close(con))
  writeLines(as.vector(rbind(headers[keep], seqs[keep])), con)
  list(total = length(seqs), kept = sum(keep))
}

for (run in runs) {
  input_fasta <- file.path(collapsed_dir, paste0("collapsed_trimmed_", run,
                                                 ".fasta.gz"))
  if (!file.exists(input_fasta)) {
    stop("Missing massiveNGSpipe collapsed reads for ", run, ": ",
         input_fasta)
  }

  filtered_fasta <- file.path(output_dir, "filtered",
                              paste0(run, ".filtered_", min_len, "_",
                                    max_len, "nt.fasta.gz"))
  if (!file.exists(filtered_fasta)) {
    message("Filtering ", run, " to ", min_len, "-", max_len, "nt")
    stats <- filter_collapsed_fasta(input_fasta, filtered_fasta, min_len,
                                    max_len)
    message("  kept ", stats$kept, "/", stats$total, " unique sequences")
  } else {
    message("Using existing filtered reads: ", filtered_fasta)
  }

  alignment_dir <- file.path(output_dir, "aligned", run)
  dir.create(alignment_dir, recursive = TRUE, showWarnings = FALSE)
  aligned_bam <- file.path(alignment_dir, "Aligned.sortedByCoord.out.bam")

  alignment_ready <- file.exists(aligned_bam) &&
    system2("samtools", c("quickcheck", aligned_bam)) == 0
  if (!alignment_ready) {
    message("Aligning ", run)
    status <- system2(star_bin, c(
      "--genomeDir", star_index,
      "--readFilesIn", filtered_fasta,
      "--readFilesCommand", "zcat",
      "--runThreadN", threads,
      "--outFileNamePrefix", paste0(alignment_dir, "/"),
      "--outSAMtype", "BAM", "SortedByCoordinate",
      "--outSAMattributes", "NH", "HI", "AS", "nM", "MD",
      "--outFilterIntronMotifs", "RemoveNoncanonicalUnannotated",
      "--outFilterMultimapNmax", "1",
      "--outFilterMismatchNoverLmax", "0.1",
      "--limitBAMsortRAM", "12000000000"
    ))
    if (status != 0) stop("STAR alignment failed for ", run)
    system2("samtools", c("index", "-@", threads, aligned_bam))
  } else {
    message("Using existing alignment: ", aligned_bam)
  }
}

message("Reused-input length-resolved alignment complete: ", output_dir)
