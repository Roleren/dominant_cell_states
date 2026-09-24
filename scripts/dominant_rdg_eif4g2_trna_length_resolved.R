#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicAlignments)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(ORFik)
  library(Rsamtools)
})

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(here, "DESCRIPTION")) &&
        file.exists(file.path(here, "R", "length-resolved-validation.R"))) {
      return(here)
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find dominant_cell_states analysis root.", call. = FALSE)
}

analysis_dir <- find_analysis_dir()
source(file.path(analysis_dir, "R", "length-resolved-validation.R"))

ribocrypt_repo <- Sys.getenv("RIBOCRYPT_REPO", unset = "")
if (!nzchar(ribocrypt_repo)) {
  stop("Set RIBOCRYPT_REPO to the local RiboCrypt checkout.", call. = FALSE)
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required.", call. = FALSE)
}
devtools::load_all(ribocrypt_repo, quiet = TRUE)

input_dir <- Sys.getenv(
  "GOODARZI_TRNA_OUTPUT_DIR",
  "/media/roler/S/data/Bio_data/runtime_inputs/dominant_cell_states/Goodarzi_trna_length_resolved"
)
output_dir <- Sys.getenv(
  "EIF4G2_TRNA_LENGTH_OUTPUT_DIR",
  file.path(analysis_dir, "results", "dominant_rdg_eif4g2_trna_length_resolved")
)
if (!grepl("^/", output_dir)) output_dir <- file.path(analysis_dir, output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
calibrate_only <- identical(
  tolower(Sys.getenv("GOODARZI_TRNA_CALIBRATE_ONLY", unset = "false")),
  "true"
)

run_table <- data.table(
  run = c(
    "SRR3131003", "SRR3131004", "SRR3131019",
    "SRR3131007", "SRR3131008", "SRR3131021",
    "SRR3131011", "SRR3131012", "SRR3131023",
    "SRR3131015", "SRR3131016", "SRR3131025",
    "SRR3130940", "SRR3130941", "SRR3130944", "SRR3130945"
  ),
  perturbation = c(rep("tRNA-Glu", 12L), rep("tRNA-Arg", 4L)),
  condition = c(
    rep("control", 6L), rep("case", 6L),
    "control", "control", "case", "case"
  ),
  replicate = c(
    1L, 1L, 1L, 2L, 2L, 2L,
    1L, 1L, 1L, 2L, 2L, 2L,
    1L, 2L, 1L, 2L
  ),
  lane = c(
    1L, 2L, 3L, 1L, 2L, 3L,
    1L, 2L, 3L, 1L, 2L, 3L,
    1L, 1L, 1L, 1L
  )
)
run_table[, biological_unit := paste(
  perturbation, condition, paste0("replicate_", replicate), sep = "::"
)]
run_table[, bam := file.path(
  input_dir, "aligned", run, "Aligned.sortedByCoord.out.bam"
)]
run_table[, bam_ready :=
  file.exists(bam) & file.info(bam)$size > 0 &
    file.exists(paste0(bam, ".bai"))]
missing_bams <- run_table[bam_ready == FALSE, bam]
if (length(missing_bams)) {
  if (!calibrate_only) {
    stop(
      "Missing aligned BAMs. Run process_goodarzi_trna_length_resolved.sh first: ",
      paste(missing_bams, collapse = ", "), call. = FALSE
    )
  }
  message(
    "Calibration-only mode: skipping ", length(missing_bams),
    " run(s) whose BAMs are not complete."
  )
  run_table <- run_table[bam_ready == TRUE]
  if (!nrow(run_table)) {
    stop("No complete BAMs are available for calibration.", call. = FALSE)
  }
}
run_table[, bam_ready := NULL]

experiment <- ORFik::read.experiment(Sys.getenv("RDG_SDRIVE_EXPERIMENT", unset = "all_merged-Homo_sapiens"))
txdb <- ORFik::loadTxdb(experiment)
target_tx_id <- "ENST00000339995"
target_tx <- ORFik::loadRegion(
  txdb, part = "mrna", names.keep = target_tx_id
)
if (length(target_tx) != 1L) {
  stop("Could not recover the exact EIF4G2 transcript.", call. = FALSE)
}

build_calibration_windows <- function(txdb) {
  candidate_names <- ORFik::filterTranscripts(
    txdb, minFiveUTR = 40L, minCDS = 300L
  )
  candidate_cds <- ORFik::loadRegion(
    txdb, part = "cds", names.keep = candidate_names
  )
  candidate_tx <- ORFik::loadRegion(
    txdb, part = "mrna", names.keep = candidate_names
  )
  common <- intersect(names(candidate_cds), names(candidate_tx))
  rows <- lapply(common, function(tx_name) {
    cds <- candidate_cds[[tx_name]]
    exons <- candidate_tx[[tx_name]]
    if (!length(cds) || !length(exons)) return(NULL)
    first_cds <- cds[1L]
    rank <- as.integer(mcols(first_cds)$exon_rank)
    host <- exons[as.integer(mcols(exons)$exon_rank) == rank]
    if (length(host) != 1L || width(first_cds) < 180L) return(NULL)
    strand_value <- as.character(strand(first_cds))
    upstream_bases <- if (strand_value == "+") {
      start(first_cds) - start(host)
    } else if (strand_value == "-") {
      end(host) - end(first_cds)
    } else {
      return(NULL)
    }
    if (!is.finite(upstream_bases) || upstream_bases < 40L) return(NULL)
    tis <- if (strand_value == "+") start(first_cds) else end(first_cds)
    window_start <- if (strand_value == "+") tis - 40L else tis - 179L
    window_end <- if (strand_value == "+") tis + 179L else tis + 40L
    GRanges(
      seqnames(first_cds), IRanges(window_start, window_end),
      strand = strand_value, tx_id = tx_name, tis = tis
    )
  })
  windows <- do.call(c, rows[!vapply(rows, is.null, logical(1))])
  key <- paste(seqnames(windows), mcols(windows)$tis, strand(windows), sep = ":")
  windows <- windows[!duplicated(key)]
  windows[countOverlaps(windows, windows, ignore.strand = TRUE) == 1L]
}

calibration_cache_file <- file.path(
  output_dir, "calibration_windows_checkpoint.rds"
)
calibration_cache_version <- 1L
calibration_checkpoint <- if (file.exists(calibration_cache_file)) {
  readRDS(calibration_cache_file)
} else {
  NULL
}
if (is.list(calibration_checkpoint) &&
    identical(calibration_checkpoint$version, calibration_cache_version) &&
    inherits(calibration_checkpoint$windows, "GRanges")) {
  calibration_windows <- calibration_checkpoint$windows
  message("Using checkpointed CDS-start calibration windows.")
} else {
  calibration_windows <- build_calibration_windows(txdb)
  if (length(calibration_windows) < 100L) {
    stop("Too few non-overlapping CDS-start calibration windows.",
         call. = FALSE)
  }
  saveRDS(
    list(version = calibration_cache_version, windows = calibration_windows),
    calibration_cache_file
  )
}
calibration_tx_names <- as.character(mcols(calibration_windows)$tx_id)
calibration_tx <- ORFik::loadRegion(
  txdb, part = "tx", names.keep = calibration_tx_names
)

alignment_flag <- scanBamFlag(
  isUnmappedQuery = FALSE,
  isSecondaryAlignment = FALSE,
  isSupplementaryAlignment = FALSE
)
calibration_param <- ScanBamParam(
  flag = alignment_flag,
  which = reduce(granges(calibration_windows), ignore.strand = TRUE)
)
target_span <- range(
  granges(unlist(target_tx, use.names = FALSE)), ignore.strand = TRUE
)
target_param <- ScanBamParam(flag = alignment_flag, which = target_span)

five_prime_points <- function(reads) {
  read_ranges <- granges(reads)
  position <- ifelse(
    as.character(strand(read_ranges)) == "+",
    start(read_ranges), end(read_ranges)
  )
  GRanges(
    seqnames(read_ranges), IRanges(position, width = 1L),
    strand = strand(read_ranges)
  )
}

relative_five_prime_table <- function(reads, windows) {
  points <- five_prime_points(reads)
  hits <- findOverlaps(points, windows, ignore.strand = FALSE)
  query <- queryHits(hits)
  subject <- subjectHits(hits)
  position <- start(points)[query]
  tis <- as.integer(mcols(windows)$tis[subject])
  is_plus <- as.character(strand(windows)[subject]) == "+"
  data.table(
    relative_five_prime = ifelse(is_plus, position - tis, tis - position),
    read_length = ORFik::readWidths(reads, along.reference = FALSE)[query]
  )
}

map_target_psites <- function(reads, shifts, target_tx) {
  shifted <- ORFik::shiftFootprints(reads, shifts, sort = FALSE)
  target_exons <- unlist(target_tx, use.names = FALSE)
  inside_target_exon <- countOverlaps(
    shifted, target_exons, type = "within", ignore.strand = FALSE
  ) > 0L
  shifted <- shifted[inside_target_exon]
  if (!length(shifted)) {
    return(data.table(read_length = integer(), tx_position = integer()))
  }
  mapped <- GenomicFeatures::mapToTranscripts(shifted, target_tx)
  if (!length(mapped)) {
    return(data.table(read_length = integer(), tx_position = integer()))
  }
  data.table(
    read_length = as.integer(mcols(shifted)$size[mcols(mapped)$xHits]),
    tx_position = start(mapped)
  )
}

offset_rows <- list()
calibration_rows <- list()
profile_rows <- list()

for (i in seq_len(nrow(run_table))) {
  run_row <- run_table[i]
  checkpoint_prefix <- file.path(output_dir, paste0("checkpoint_", run_row$run))
  offset_checkpoint <- paste0(checkpoint_prefix, "_offsets.csv")
  calibration_checkpoint <- paste0(checkpoint_prefix, "_calibration.csv")
  profile_checkpoint <- paste0(checkpoint_prefix, "_profile.csv")
  if (all(file.exists(c(
    offset_checkpoint, calibration_checkpoint, profile_checkpoint
  )))) {
    message("Using checkpointed P-site validation for ", run_row$run)
    offset_rows[[run_row$run]] <- fread(offset_checkpoint)
    calibration_rows[[run_row$run]] <- fread(calibration_checkpoint)
    profile_rows[[run_row$run]] <- fread(profile_checkpoint)
    next
  }

  message("Calibrating P-sites for ", run_row$run)
  calibration_reads <- readGAlignments(
    run_row$bam, use.names = FALSE, param = calibration_param
  )
  calibration_reads <- calibration_reads[
    ORFik::readWidths(calibration_reads, along.reference = FALSE) %in% 15:40
  ]
  shifts <- ORFik::detectRibosomeShifts(
    calibration_reads,
    txdb = txdb,
    start = TRUE,
    stop = FALSE,
    top_tx = 10L,
    minFiveUTR = 40L,
    minCDS = 300L,
    txNames = calibration_tx_names,
    firstN = 150L,
    tx = calibration_tx,
    min_reads = 100L,
    min_reads_TIS = 25L,
    accepted.lengths = 15:40,
    must.be.periodic = FALSE,
    strict.fft = FALSE,
    verbose = TRUE
  )
  shifts <- as.data.table(shifts)
  if (!nrow(shifts)) stop("No P-site offsets recovered for ", run_row$run)
  shifts[, physically_plausible_offset :=
           is.finite(fraction) & is.finite(offsets_start) &
           fraction > 0L & offsets_start <= 0L &
           -offsets_start < fraction]
  rejected <- shifts[physically_plausible_offset == FALSE]
  if (nrow(rejected)) {
    warning(
      "Discarding physically impossible P-site offsets for ", run_row$run,
      ": ", paste0(rejected$fraction, "nt=>", rejected$offsets_start,
                    collapse = "; ")
    )
  }
  shifts <- shifts[physically_plausible_offset == TRUE]
  shifts[, physically_plausible_offset := NULL]
  if (!nrow(shifts)) {
    stop("No physically plausible P-site offsets recovered for ", run_row$run)
  }
  shifts[, `:=`(
    run = run_row$run,
    perturbation = run_row$perturbation,
    condition = run_row$condition,
    replicate = run_row$replicate,
    lane = run_row$lane
  )]

  relative <- relative_five_prime_table(
    calibration_reads, calibration_windows
  )
  qc <- as.data.table(dcs_psite_calibration_qc(
    relative$relative_five_prime,
    relative$read_length,
    shifts[, .(fraction, offsets_start)]
  ))
  qc[, `:=`(
    run = run_row$run,
    perturbation = run_row$perturbation,
    condition = run_row$condition,
    replicate = run_row$replicate,
    lane = run_row$lane
  )]

  message("Mapping EIF4G2 P-sites for ", run_row$run)
  target_reads <- readGAlignments(
    run_row$bam, use.names = FALSE, param = target_param
  )
  target_profile <- map_target_psites(
    target_reads, shifts[, .(fraction, offsets_start)], target_tx
  )
  target_profile[, `:=`(
    run = run_row$run,
    perturbation = run_row$perturbation,
    condition = run_row$condition,
    replicate = run_row$replicate,
    lane = run_row$lane
  )]

  offset_rows[[run_row$run]] <- shifts
  calibration_rows[[run_row$run]] <- qc
  profile_rows[[run_row$run]] <- target_profile
  fwrite(shifts, offset_checkpoint)
  fwrite(qc, calibration_checkpoint)
  fwrite(target_profile, profile_checkpoint)
  rm(calibration_reads, target_reads)
  invisible(gc())
}

if (calibrate_only) {
  message(
    "Calibration-only checkpoint complete for ", nrow(run_table),
    " available run(s): ", paste(run_table$run, collapse = ";")
  )
  quit(save = "no", status = 0L)
}

run_table[, uniquely_mapped_reads := vapply(
  bam,
  function(path) sum(Rsamtools::idxstatsBam(path)$mapped),
  numeric(1)
)]

offsets <- rbindlist(offset_rows, fill = TRUE)
calibration <- rbindlist(calibration_rows, fill = TRUE)
profile <- rbindlist(profile_rows, fill = TRUE)

study_run_counts <- run_table[, .(required_runs = uniqueN(run)),
                              by = perturbation]
shared_lengths <- calibration[periodicity_gate == TRUE,
  .(passing_runs = uniqueN(run)), by = .(perturbation, read_length)
]
shared_lengths <- merge(
  shared_lengths, study_run_counts, by = "perturbation", all.x = TRUE
)[passing_runs == required_runs]

lane_run_counts <- run_table[, .(required_runs = uniqueN(run)),
                             by = .(perturbation, lane)]
lane_shared_lengths <- calibration[periodicity_gate == TRUE,
  .(passing_runs = uniqueN(run)), by = .(perturbation, lane, read_length)
]
lane_shared_lengths <- merge(
  lane_shared_lengths, lane_run_counts,
  by = c("perturbation", "lane"), all.x = TRUE
)[passing_runs == required_runs]

perturbation_lane_counts <- run_table[, .(required_lanes = uniqueN(lane)),
                                      by = perturbation]
lane_coverage <- lane_shared_lengths[, .(lanes_with_shared_length = uniqueN(lane)),
                                     by = perturbation]
lane_coverage <- merge(
  perturbation_lane_counts, lane_coverage,
  by = "perturbation", all.x = TRUE, sort = FALSE
)
lane_coverage[is.na(lanes_with_shared_length), lanes_with_shared_length := 0L]
lane_coverage[, all_lanes_have_shared_length :=
                lanes_with_shared_length == required_lanes]

lane_category_membership <- rbindlist(list(
  lane_shared_lengths[, .(
    perturbation, lane,
    category = "all_lane_matched_well_phased"
  )],
  lane_shared_lengths[read_length >= 20L & read_length <= 23L, .(
    perturbation, lane,
    category = "short_20_23nt_lane_matched_well_phased"
  )],
  lane_shared_lengths[read_length >= 28L & read_length <= 32L, .(
    perturbation, lane,
    category = "canonical_28_32nt_lane_matched_well_phased"
  )]
), use.names = TRUE)
lane_category_coverage <- unique(lane_category_membership)[, .(
  category_lanes_with_shared_length = uniqueN(lane)
), by = .(perturbation, category)]
lane_category_coverage <- merge(
  CJ(
    perturbation = unique(run_table$perturbation),
    category = c(
      "all_lane_matched_well_phased",
      "short_20_23nt_lane_matched_well_phased",
      "canonical_28_32nt_lane_matched_well_phased"
    ),
    unique = TRUE
  ),
  lane_category_coverage,
  by = c("perturbation", "category"), all.x = TRUE, sort = FALSE
)
lane_category_coverage <- merge(
  lane_category_coverage, perturbation_lane_counts,
  by = "perturbation", all.x = TRUE, sort = FALSE
)
lane_category_coverage[
  is.na(category_lanes_with_shared_length),
  category_lanes_with_shared_length := 0L
]
lane_category_coverage[, all_category_lanes_have_shared_length :=
                         category_lanes_with_shared_length == required_lanes]

profile[, feature := fcase(
  tx_position >= 23L & tx_position <= 73L, "uORF1",
  tx_position >= 128L & tx_position <= 220L, "uORF2",
  tx_position >= 309L & tx_position <= 3032L, "clean_CDS",
  default = "other_transcript"
)]
profile[, uorf_phase := fcase(
  feature == "uORF1", (tx_position - 23L) %% 3L,
  feature == "uORF2", (tx_position - 128L) %% 3L,
  default = NA_integer_
)]
profile[, length_stratum := dcs_read_length_stratum(read_length)]
profile[, study_shared_well_phased := paste(perturbation, read_length) %chin%
          paste(shared_lengths$perturbation, shared_lengths$read_length)]
profile[, lane_shared_well_phased := paste(perturbation, lane, read_length) %chin%
          paste(
            lane_shared_lengths$perturbation,
            lane_shared_lengths$lane,
            lane_shared_lengths$read_length
          )]

expand_categories <- function(dt) {
  exact <- dt[, .(
    run, perturbation, condition, replicate, lane, read_length,
    tx_position, feature, uorf_phase,
    category = paste0(read_length, "nt")
  )]
  combined <- rbindlist(list(
    exact,
    dt[, .(
      run, perturbation, condition, replicate, lane, read_length,
      tx_position, feature, uorf_phase, category = "all_calibrated"
    )],
    dt[study_shared_well_phased == TRUE, .(
      run, perturbation, condition, replicate, lane, read_length,
      tx_position, feature, uorf_phase,
      category = "all_study_shared_well_phased"
    )],
    dt[study_shared_well_phased == TRUE &
         read_length >= 28L & read_length <= 32L, .(
      run, perturbation, condition, replicate, lane, read_length,
      tx_position, feature, uorf_phase,
      category = "canonical_28_32nt_study_shared_well_phased"
    )],
    dt[lane_shared_well_phased == TRUE, .(
      run, perturbation, condition, replicate, lane, read_length,
      tx_position, feature, uorf_phase,
      category = "all_lane_matched_well_phased"
    )],
    dt[lane_shared_well_phased == TRUE &
         read_length >= 20L & read_length <= 23L, .(
      run, perturbation, condition, replicate, lane, read_length,
      tx_position, feature, uorf_phase,
      category = "short_20_23nt_lane_matched_well_phased"
    )],
    dt[lane_shared_well_phased == TRUE &
         read_length >= 28L & read_length <= 32L, .(
      run, perturbation, condition, replicate, lane, read_length,
      tx_position, feature, uorf_phase,
      category = "canonical_28_32nt_lane_matched_well_phased"
    )]
  ), fill = TRUE)

  lane_specific <- dt[lane_shared_well_phased == TRUE, .(
    run, perturbation, condition, replicate, lane, read_length,
    tx_position, feature, uorf_phase,
    category = paste0("lane", lane, "_shared_well_phased")
  )]
  rbindlist(list(combined, lane_specific), fill = TRUE)
}

expanded <- expand_categories(profile)
feature_levels <- c("uORF1", "uORF2", "clean_CDS")
categories <- sort(unique(expanded$category))
grid <- CJ(run = run_table$run, category = categories,
           feature = feature_levels, unique = TRUE)
feature_counts <- expanded[feature %chin% feature_levels,
                           .(counts = .N),
                           by = .(run, category, feature)]
feature_counts <- merge(
  grid, feature_counts, by = c("run", "category", "feature"),
  all.x = TRUE, sort = FALSE
)
feature_counts[is.na(counts), counts := 0L]
feature_counts <- merge(
  feature_counts,
  run_table[, .(
    run, perturbation, condition, replicate, lane, biological_unit,
    uniquely_mapped_reads
  )],
  by = "run", all.x = TRUE, sort = FALSE
)
run_metrics <- dcast(
  feature_counts,
  run + perturbation + condition + replicate + lane + biological_unit +
    uniquely_mapped_reads + category ~ feature,
  value.var = "counts", fill = 0L
)
run_metrics[, union_uorf := uORF1 + uORF2]
run_metrics[, clean_cds_allocation := clean_CDS /
              pmax(clean_CDS + union_uorf, 1L)]
run_metrics[, `:=`(
  uorf1_cpm = uORF1 / pmax(uniquely_mapped_reads, 1) * 1e6,
  uorf_cpm = union_uorf / pmax(uniquely_mapped_reads, 1) * 1e6,
  clean_cds_cpm = clean_CDS / pmax(uniquely_mapped_reads, 1) * 1e6
)]

uorf_phase_counts <- expanded[
  feature == "uORF1" & is.finite(uorf_phase),
  .(
    uorf_phase0_counts = sum(uorf_phase == 0L),
    uorf_phase1_counts = sum(uorf_phase == 1L),
    uorf_phase2_counts = sum(uorf_phase == 2L)
  ), by = .(run, category)
]
run_metrics <- merge(
  run_metrics, uorf_phase_counts, by = c("run", "category"),
  all.x = TRUE, sort = FALSE
)
for (column in c(
  "uorf_phase0_counts", "uorf_phase1_counts", "uorf_phase2_counts"
)) {
  set(run_metrics, which(is.na(run_metrics[[column]])), column, 0L)
}
run_metrics[, uorf_phase0_fraction := uorf_phase0_counts /
              pmax(uorf_phase0_counts + uorf_phase1_counts +
                     uorf_phase2_counts, 1L)]

unit_metrics <- run_metrics[, .(
  lane_count = uniqueN(lane),
  runs = paste(sort(unique(run)), collapse = ";"),
  uniquely_mapped_reads = sum(uniquely_mapped_reads),
  uORF1 = sum(uORF1),
  uORF2 = sum(uORF2),
  union_uorf = sum(union_uorf),
  clean_CDS = sum(clean_CDS),
  uorf_phase0_counts = sum(uorf_phase0_counts),
  uorf_phase1_counts = sum(uorf_phase1_counts),
  uorf_phase2_counts = sum(uorf_phase2_counts)
), by = .(
  perturbation, condition, replicate, biological_unit, category
)]
unit_metrics[, clean_cds_allocation := clean_CDS /
               pmax(clean_CDS + union_uorf, 1L)]
unit_metrics[, clean_cds_vs_uorf1_allocation := clean_CDS /
               pmax(clean_CDS + uORF1, 1L)]
unit_metrics[, `:=`(
  uorf1_cpm = uORF1 / pmax(uniquely_mapped_reads, 1) * 1e6,
  uorf_cpm = union_uorf / pmax(uniquely_mapped_reads, 1) * 1e6,
  clean_cds_cpm = clean_CDS / pmax(uniquely_mapped_reads, 1) * 1e6
)]
unit_metrics[, uorf_phase0_fraction := uorf_phase0_counts /
               pmax(uorf_phase0_counts + uorf_phase1_counts +
                      uorf_phase2_counts, 1L)]

aggregate_delta <- function(dt, shift = 0L) {
  shifted <- copy(dt)
  shifted[, shifted_position := tx_position + as.integer(shift)]
  shifted[, shifted_feature := fcase(
    shifted_position >= 23L & shifted_position <= 73L, "uORF1",
    shifted_position >= 128L & shifted_position <= 220L, "uORF2",
    shifted_position >= 309L & shifted_position <= 3032L, "clean_CDS",
    default = "other_transcript"
  )]
  x <- shifted[shifted_feature %chin% feature_levels,
               .N, by = .(
                 perturbation, category, condition, shifted_feature
               )]
  x <- dcast(
    x, perturbation + category + condition ~ shifted_feature,
    value.var = "N", fill = 0L
  )
  for (feature in feature_levels) {
    if (!feature %in% names(x)) x[, (feature) := 0L]
  }
  x[, union_uorf := uORF1 + uORF2]
  x[, allocation := clean_CDS / pmax(clean_CDS + union_uorf, 1L)]
  wide <- dcast(
    x, perturbation + category ~ condition, value.var = "allocation"
  )
  if (!all(c("case", "control") %in% names(wide))) {
    return(data.table())
  }
  wide[, .(
    perturbation,
    category,
    coordinate_shift = as.integer(shift),
    clean_cds_allocation_delta = case - control
  )]
}

shift_sensitivity <- rbindlist(lapply(-2:2, function(shift) {
  aggregate_delta(expanded, shift)
}))

observed_position_counts <- expanded[
  feature %chin% feature_levels,
  .(counts = .N),
  by = .(
    perturbation, category, condition,
    region = fifelse(feature == "clean_CDS", "clean_CDS", "uorf"),
    position = tx_position
  )
]
position_grid <- unique(observed_position_counts[, .(
  perturbation, category, region, position
)])[, .(condition = c("case", "control")),
    by = .(perturbation, category, region, position)]
position_counts <- merge(
  position_grid, observed_position_counts,
  by = c("perturbation", "category", "condition", "region", "position"),
  all.x = TRUE, sort = FALSE
)
position_counts[is.na(counts), counts := 0]
position_counts[, audit_category := paste(perturbation, category, sep = "||")]
position_leaveout <- dcs_position_leaveout_audit(
  position_counts[, .(
    category = audit_category, condition, region, position, counts
  )],
  case_label = "case",
  control_label = "control",
  max_greedy_deletions = 3L
)
position_leaveout_summary <- as.data.table(position_leaveout$summary)
position_leaveout_rows <- as.data.table(position_leaveout$leaveout)
setnames(position_leaveout_summary, "category", "audit_category")
position_leaveout_summary[, c("perturbation", "category") := tstrsplit(
  audit_category, "||", fixed = TRUE
)]
position_leaveout_summary[, audit_category := NULL]
setnames(position_leaveout_rows, "category", "audit_category")
position_leaveout_rows[, c("perturbation", "category") := tstrsplit(
  audit_category, "||", fixed = TRUE
)]
position_leaveout_rows[, audit_category := NULL]

contrast <- unit_metrics[, {
  case <- .SD[condition == "case"]
  control <- .SD[condition == "control"]
  case_uorf1 <- sum(case$uORF1)
  control_uorf1 <- sum(control$uORF1)
  case_uorf <- sum(case$union_uorf)
  control_uorf <- sum(control$union_uorf)
  case_cds <- sum(case$clean_CDS)
  control_cds <- sum(control$clean_CDS)
  case_phase <- c(
    sum(case$uorf_phase0_counts), sum(case$uorf_phase1_counts),
    sum(case$uorf_phase2_counts)
  )
  control_phase <- c(
    sum(control$uorf_phase0_counts), sum(control$uorf_phase1_counts),
    sum(control$uorf_phase2_counts)
  )
  case_allocation <- case_cds / pmax(case_cds + case_uorf, 1L)
  control_allocation <- control_cds /
    pmax(control_cds + control_uorf, 1L)
  delta <- case_allocation - control_allocation
  case_uorf1_allocation <- case_cds / pmax(case_cds + case_uorf1, 1L)
  control_uorf1_allocation <- control_cds /
    pmax(control_cds + control_uorf1, 1L)
  uorf1_delta <- case_uorf1_allocation - control_uorf1_allocation
  case_uorf1_cpm <- mean(case$uorf1_cpm)
  control_uorf1_cpm <- mean(control$uorf1_cpm)
  case_uorf_cpm <- mean(case$uorf_cpm)
  control_uorf_cpm <- mean(control$uorf_cpm)
  case_clean_cds_cpm <- mean(case$clean_cds_cpm)
  control_clean_cds_cpm <- mean(control$clean_cds_cpm)
  range_separation <- if (delta < 0) {
    min(control$clean_cds_allocation) - max(case$clean_cds_allocation)
  } else {
    min(case$clean_cds_allocation) - max(control$clean_cds_allocation)
  }
  uorf1_range_separation <- if (uorf1_delta < 0) {
    min(control$clean_cds_vs_uorf1_allocation) -
      max(case$clean_cds_vs_uorf1_allocation)
  } else {
    min(case$clean_cds_vs_uorf1_allocation) -
      max(control$clean_cds_vs_uorf1_allocation)
  }
  list(
    n_case_units = nrow(case),
    n_control_units = nrow(control),
    case_uorf1_counts = case_uorf1,
    control_uorf1_counts = control_uorf1,
    case_union_uorf_counts = case_uorf,
    control_union_uorf_counts = control_uorf,
    case_clean_cds_counts = case_cds,
    control_clean_cds_counts = control_cds,
    case_clean_cds_allocation = case_allocation,
    control_clean_cds_allocation = control_allocation,
    clean_cds_allocation_delta = delta,
    case_clean_cds_vs_uorf1_allocation = case_uorf1_allocation,
    control_clean_cds_vs_uorf1_allocation = control_uorf1_allocation,
    clean_cds_vs_uorf1_allocation_delta = uorf1_delta,
    case_mean_uorf1_cpm = case_uorf1_cpm,
    control_mean_uorf1_cpm = control_uorf1_cpm,
    uorf1_cpm_log2_ratio = log2(
      (case_uorf1_cpm + 0.5) / (control_uorf1_cpm + 0.5)
    ),
    case_mean_uorf_cpm = case_uorf_cpm,
    control_mean_uorf_cpm = control_uorf_cpm,
    uorf_cpm_log2_ratio = log2(
      (case_uorf_cpm + 0.5) / (control_uorf_cpm + 0.5)
    ),
    case_mean_clean_cds_cpm = case_clean_cds_cpm,
    control_mean_clean_cds_cpm = control_clean_cds_cpm,
    clean_cds_cpm_log2_ratio = log2(
      (case_clean_cds_cpm + 0.5) / (control_clean_cds_cpm + 0.5)
    ),
    biological_unit_range_separation = range_separation,
    uorf1_biological_unit_range_separation = uorf1_range_separation,
    case_allocation_min = min(case$clean_cds_allocation),
    case_allocation_max = max(case$clean_cds_allocation),
    control_allocation_min = min(control$clean_cds_allocation),
    control_allocation_max = max(control$clean_cds_allocation),
    case_uorf_phase0_fraction = case_phase[[1L]] / pmax(sum(case_phase), 1L),
    case_uorf_phase1_fraction = case_phase[[2L]] / pmax(sum(case_phase), 1L),
    case_uorf_phase2_fraction = case_phase[[3L]] / pmax(sum(case_phase), 1L),
    control_uorf_phase0_fraction = control_phase[[1L]] /
      pmax(sum(control_phase), 1L),
    control_uorf_phase1_fraction = control_phase[[2L]] /
      pmax(sum(control_phase), 1L),
    control_uorf_phase2_fraction = control_phase[[3L]] /
      pmax(sum(control_phase), 1L)
  )
}, by = .(perturbation, category)]

shift_summary <- shift_sensitivity[, .(
  coordinate_shift_direction_stable =
    length(unique(sign(clean_cds_allocation_delta))) == 1L,
  coordinate_shift_min_abs_delta = min(abs(clean_cds_allocation_delta))
), by = .(perturbation, category)]
contrast <- merge(
  contrast, shift_summary, by = c("perturbation", "category"),
  all.x = TRUE, sort = FALSE
)
contrast <- merge(
  contrast, position_leaveout_summary,
  by = c("perturbation", "category"), all.x = TRUE, sort = FALSE
)

contrast[, read_length := suppressWarnings(as.integer(
  sub("nt$", "", category)
))]
contrast[, all_runs_phase_gate :=
  (!is.na(read_length) & paste(perturbation, read_length) %chin%
     paste(shared_lengths$perturbation, shared_lengths$read_length)) |
  category %chin% c(
    "all_study_shared_well_phased",
    "canonical_28_32nt_study_shared_well_phased"
  )]
contrast <- merge(
  contrast,
  lane_category_coverage[, .(
    perturbation, category, all_category_lanes_have_shared_length
  )],
  by = c("perturbation", "category"), all.x = TRUE, sort = FALSE
)
contrast[is.na(all_category_lanes_have_shared_length),
         all_category_lanes_have_shared_length := FALSE]
contrast[, calibration_design_gate :=
  all_runs_phase_gate == TRUE |
    (all_category_lanes_have_shared_length == TRUE & category %chin% c(
      "all_lane_matched_well_phased",
      "short_20_23nt_lane_matched_well_phased",
      "canonical_28_32nt_lane_matched_well_phased"
    ))]
contrast[, target_uorf_phase_gate :=
  case_uorf_phase0_fraction >= 0.40 &
    case_uorf_phase0_fraction > pmax(
      case_uorf_phase1_fraction, case_uorf_phase2_fraction
    ) &
    control_uorf_phase0_fraction >= 0.40 &
    control_uorf_phase0_fraction > pmax(
      control_uorf_phase1_fraction, control_uorf_phase2_fraction
    )]
contrast[, measurement_gate :=
  calibration_design_gate == TRUE &
    n_case_units == 2L & n_control_units == 2L &
    case_uorf1_counts >= 30L & control_uorf1_counts >= 30L &
    case_union_uorf_counts >= 30L & control_union_uorf_counts >= 30L &
    case_clean_cds_counts >= 100L & control_clean_cds_counts >= 100L &
    target_uorf_phase_gate == TRUE]
contrast[, evidence_gate :=
  measurement_gate == TRUE &
    abs(clean_cds_allocation_delta) >= 0.05 &
    biological_unit_range_separation > 0 &
    sign(clean_cds_vs_uorf1_allocation_delta) ==
      sign(clean_cds_allocation_delta) &
    abs(clean_cds_vs_uorf1_allocation_delta) >= 0.03 &
    uorf1_biological_unit_range_separation > 0 &
    coordinate_shift_direction_stable == TRUE &
    coordinate_shift_min_abs_delta >= 0.03 &
    single_position_direction_stable == TRUE &
    single_position_min_abs_delta >= 0.03 &
    greedy_direction_stable == TRUE &
    abs(greedy_final_clean_cds_allocation_delta) >= 0.03]

specificity <- merge(
  contrast[perturbation == "tRNA-Glu", .(
    category,
    glu_clean_cds_allocation_delta = clean_cds_allocation_delta,
    glu_measurement_gate = measurement_gate,
    glu_evidence_gate = evidence_gate
  )],
  contrast[perturbation == "tRNA-Arg", .(
    category,
    arg_clean_cds_allocation_delta = clean_cds_allocation_delta,
    arg_measurement_gate = measurement_gate,
    arg_evidence_gate = evidence_gate
  )],
  by = "category", all = TRUE, sort = FALSE
)
specificity[, glu_minus_arg_clean_cds_allocation_delta :=
              glu_clean_cds_allocation_delta -
                arg_clean_cds_allocation_delta]
specificity[, specificity_gate :=
  glu_measurement_gate == TRUE & arg_measurement_gate == TRUE &
    glu_evidence_gate == TRUE & abs(arg_clean_cds_allocation_delta) < 0.05 &
    glu_minus_arg_clean_cds_allocation_delta <= -0.05]

lane_length_metrics <- run_metrics[grepl("^[0-9]+nt$", category), {
  case <- .SD[condition == "case"]
  control <- .SD[condition == "control"]
  case_uorf <- sum(case$union_uorf)
  control_uorf <- sum(control$union_uorf)
  case_cds <- sum(case$clean_CDS)
  control_cds <- sum(control$clean_CDS)
  case_allocation <- case_cds / pmax(case_cds + case_uorf, 1L)
  control_allocation <- control_cds /
    pmax(control_cds + control_uorf, 1L)
  list(
    n_case_runs = nrow(case),
    n_control_runs = nrow(control),
    case_union_uorf_counts = case_uorf,
    control_union_uorf_counts = control_uorf,
    case_clean_cds_counts = case_cds,
    control_clean_cds_counts = control_cds,
    clean_cds_allocation_delta = case_allocation - control_allocation
  )
}, by = .(perturbation, lane, category)]
lane_length_metrics[, read_length := as.integer(sub("nt$", "", category))]
lane_length_metrics <- merge(
  lane_length_metrics,
  lane_shared_lengths[, .(
    perturbation, lane, read_length, lane_all_runs_phase_gate = TRUE
  )],
  by = c("perturbation", "lane", "read_length"),
  all.x = TRUE, sort = FALSE
)
lane_length_metrics[is.na(lane_all_runs_phase_gate),
                    lane_all_runs_phase_gate := FALSE]

position_profile <- profile[, .N, by = .(
  run, perturbation, condition, replicate, lane, read_length,
  length_stratum, study_shared_well_phased, lane_shared_well_phased,
  tx_position, feature
)]

fwrite(run_table, file.path(output_dir, "run_manifest.csv"))
fwrite(data.table(
  n_calibration_windows = length(calibration_windows),
  perturbation = c("tRNA-Glu", "tRNA-Arg"),
  shared_well_phased_lengths = vapply(
    c("tRNA-Glu", "tRNA-Arg"),
    function(value) paste(
      shared_lengths[perturbation == value, read_length], collapse = ";"
    ), character(1)
  )
), file.path(output_dir, "calibration_design.csv"))
fwrite(shared_lengths,
       file.path(output_dir, "study_shared_well_phased_lengths.csv"))
fwrite(lane_shared_lengths,
       file.path(output_dir, "lane_shared_well_phased_lengths.csv"))
fwrite(lane_coverage,
       file.path(output_dir, "lane_shared_length_coverage.csv"))
fwrite(lane_category_coverage,
       file.path(output_dir, "lane_category_shared_length_coverage.csv"))
fwrite(offsets, file.path(output_dir, "psite_offsets.csv"))
fwrite(calibration, file.path(output_dir, "psite_offset_qc.csv"))
fwrite(position_profile, file.path(output_dir, "eif4g2_position_profile.csv"))
fwrite(run_metrics, file.path(output_dir, "eif4g2_run_length_metrics.csv"))
fwrite(unit_metrics, file.path(output_dir, "eif4g2_unit_length_metrics.csv"))
fwrite(shift_sensitivity,
       file.path(output_dir, "eif4g2_coordinate_shift_sensitivity.csv"))
fwrite(position_leaveout_rows,
       file.path(output_dir, "eif4g2_position_leaveout_sensitivity.csv"))
fwrite(lane_length_metrics,
       file.path(output_dir, "eif4g2_lane_length_contrast_summary.csv"))
fwrite(contrast,
       file.path(output_dir, "eif4g2_length_contrast_summary.csv"))
fwrite(specificity,
       file.path(output_dir, "eif4g2_length_specificity_summary.csv"))

message("Study-shared well-phased lengths:")
print(shared_lengths)
message("Evidence-gate categories:")
print(contrast[evidence_gate == TRUE])
message("tRNA-Glu-specific evidence-gate categories:")
print(specificity[specificity_gate == TRUE])
message("Saved tRNA length-resolved EIF4G2 validation: ", output_dir)
