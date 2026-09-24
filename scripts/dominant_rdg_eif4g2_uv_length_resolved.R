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
  "GSE141459_OUTPUT_DIR",
  "/media/roler/S/data/Bio_data/runtime_inputs/dominant_cell_states/GSE141459_length_resolved"
)
output_dir <- Sys.getenv(
  "EIF4G2_UV_LENGTH_OUTPUT_DIR",
  file.path(analysis_dir, "results", "dominant_rdg_eif4g2_uv_length_resolved")
)
if (!grepl("^/", output_dir)) output_dir <- file.path(analysis_dir, output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

run_table <- data.table(
  run = paste0("SRR1156908", 0:3),
  condition = c("control", "control", "UV", "UV"),
  replicate = c(1L, 2L, 1L, 2L)
)
run_table[, bam := file.path(
  input_dir, "aligned", run, "Aligned.sortedByCoord.out.bam"
)]
missing_bams <- run_table[!file.exists(bam) | file.info(bam)$size <= 0, bam]
if (length(missing_bams)) {
  stop(
    "Missing aligned BAMs. Run process_gse141459_length_resolved.sh first: ",
    paste(missing_bams, collapse = ", "), call. = FALSE
  )
}

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
    # `candidate_cds[[tx_name]]` is a GRanges; single ranges must be
    # selected with `[`, because GRanges deliberately does not implement `[[`.
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
  windows <- windows[
    countOverlaps(windows, windows, ignore.strand = TRUE) == 1L
  ]
  if (length(windows) < 100L) {
    stop("Too few non-overlapping CDS-start calibration windows.",
         call. = FALSE)
  }
  windows
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
target_span <- range(granges(unlist(target_tx, use.names = FALSE)),
                     ignore.strand = TRUE)
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
    ORFik::readWidths(calibration_reads, along.reference = FALSE) %in% 15:34
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
    accepted.lengths = 15:34,
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
  rejected_shifts <- shifts[physically_plausible_offset == FALSE]
  if (nrow(rejected_shifts)) {
    warning(
      "Discarding physically impossible P-site offsets for ",
      run_row$run, ": ",
      paste0(
        rejected_shifts$fraction, "nt=>", rejected_shifts$offsets_start,
        collapse = "; "
      )
    )
  }
  shifts <- shifts[physically_plausible_offset == TRUE]
  shifts[, physically_plausible_offset := NULL]
  if (!nrow(shifts)) {
    stop("No physically plausible P-site offsets recovered for ",
         run_row$run)
  }
  shifts[, `:=`(
    run = run_row$run,
    condition = run_row$condition,
    replicate = run_row$replicate
  )]
  offset_rows[[run_row$run]] <- shifts

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
    condition = run_row$condition,
    replicate = run_row$replicate
  )]
  calibration_rows[[run_row$run]] <- qc

  message("Mapping EIF4G2 P-sites for ", run_row$run)
  target_reads <- readGAlignments(
    run_row$bam, use.names = FALSE, param = target_param
  )
  target_profile <- map_target_psites(
    target_reads, shifts[, .(fraction, offsets_start)], target_tx
  )
  target_profile[, `:=`(
    run = run_row$run,
    condition = run_row$condition,
    replicate = run_row$replicate
  )]
  profile_rows[[run_row$run]] <- target_profile
  fwrite(shifts, offset_checkpoint)
  fwrite(qc, calibration_checkpoint)
  fwrite(target_profile, profile_checkpoint)
  rm(calibration_reads, target_reads)
  invisible(gc())
}

offsets <- rbindlist(offset_rows, fill = TRUE)
calibration <- rbindlist(calibration_rows, fill = TRUE)
profile <- rbindlist(profile_rows, fill = TRUE)
setcolorder(offsets, c(
  "run", "condition", "replicate", "fraction", "offsets_start"
))
setcolorder(calibration, c(
  "run", "condition", "replicate", "read_length", "offsets_start"
))

shared_well_phased_lengths <- calibration[
  periodicity_gate == TRUE,
  .(n_runs = uniqueN(run)), by = read_length
][n_runs == nrow(run_table), read_length]
if (!length(shared_well_phased_lengths)) {
  warning("No read length passed the phase gate in all four libraries.")
}

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
profile[, shared_well_phased := read_length %in% shared_well_phased_lengths]

expand_categories <- function(dt) {
  rbindlist(list(
    dt[, .(run, condition, replicate, read_length, tx_position, feature,
           uorf_phase,
           category = paste0(read_length, "nt"))],
    dt[, .(run, condition, replicate, read_length, tx_position, feature,
           uorf_phase,
           category = "all_calibrated")],
    dt[, .(run, condition, replicate, read_length, tx_position, feature,
           uorf_phase,
           category = length_stratum)],
    dt[shared_well_phased == TRUE,
       .(run, condition, replicate, read_length, tx_position, feature,
         uorf_phase,
         category = "all_shared_well_phased")],
    dt[shared_well_phased == TRUE & length_stratum == "short_20_23nt",
       .(run, condition, replicate, read_length, tx_position, feature,
         uorf_phase,
         category = "short_20_23nt_shared_well_phased")],
    dt[shared_well_phased == TRUE & length_stratum == "canonical_28_32nt",
       .(run, condition, replicate, read_length, tx_position, feature,
         uorf_phase,
         category = "canonical_28_32nt_shared_well_phased")]
  ), fill = TRUE)
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
  feature_counts, run_table[, .(run, condition, replicate)],
  by = "run", all.x = TRUE, sort = FALSE
)
run_metrics <- dcast(
  feature_counts, run + condition + replicate + category ~ feature,
  value.var = "counts", fill = 0L
)
run_metrics[, union_uorf := uORF1 + uORF2]
run_metrics[, clean_cds_allocation := clean_CDS /
              pmax(clean_CDS + union_uorf, 1L)]
run_metrics[, log2_uorf_to_cds_density := log2(
  ((union_uorf + 0.5) / 144) / ((clean_CDS + 0.5) / 2724)
)]
uorf_phase_counts <- expanded[
  feature %chin% c("uORF1", "uORF2") & is.finite(uorf_phase),
  .(
    uorf_phase0_counts = sum(uorf_phase == 0L),
    uorf_phase1_counts = sum(uorf_phase == 1L),
    uorf_phase2_counts = sum(uorf_phase == 2L)
  ),
  by = .(run, category)
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
               .N, by = .(category, condition, shifted_feature)]
  x <- dcast(x, category + condition ~ shifted_feature,
             value.var = "N", fill = 0L)
  for (feature in feature_levels) {
    if (!feature %in% names(x)) x[, (feature) := 0L]
  }
  x[, union_uorf := uORF1 + uORF2]
  x[, allocation := clean_CDS / pmax(clean_CDS + union_uorf, 1L)]
  wide <- dcast(x, category ~ condition, value.var = "allocation")
  if (!all(c("UV", "control") %in% names(wide))) {
    return(data.table())
  }
  wide[, .(
    category,
    coordinate_shift = as.integer(shift),
    clean_cds_allocation_delta = UV - control
  )]
}

shift_sensitivity <- rbindlist(lapply(-2:2, function(shift) {
  aggregate_delta(expanded, shift)
}))

observed_position_counts <- expanded[
  feature %chin% feature_levels,
  .(counts = .N),
  by = .(
    category, condition,
    region = fifelse(feature == "clean_CDS", "clean_CDS", "uorf"),
    position = tx_position
  )
]
position_grid <- unique(observed_position_counts[, .(
  category, region, position
)])[, .(condition = c("UV", "control")),
    by = .(category, region, position)]
position_counts <- merge(
  position_grid, observed_position_counts,
  by = c("category", "condition", "region", "position"),
  all.x = TRUE, sort = FALSE
)
position_counts[is.na(counts), counts := 0]
position_leaveout <- dcs_position_leaveout_audit(
  position_counts,
  case_label = "UV",
  control_label = "control",
  max_greedy_deletions = 3L
)
position_leaveout_summary <- as.data.table(position_leaveout$summary)
position_leaveout_rows <- as.data.table(position_leaveout$leaveout)

contrast <- run_metrics[, {
  case <- .SD[condition == "UV"]
  control <- .SD[condition == "control"]
  case_uorf <- sum(case$union_uorf)
  control_uorf <- sum(control$union_uorf)
  case_cds <- sum(case$clean_CDS)
  control_cds <- sum(control$clean_CDS)
  case_phase_counts <- colSums(case[, .(
    uorf_phase0_counts, uorf_phase1_counts, uorf_phase2_counts
  )])
  control_phase_counts <- colSums(control[, .(
    uorf_phase0_counts, uorf_phase1_counts, uorf_phase2_counts
  )])
  case_phase_fraction <- case_phase_counts /
    pmax(sum(case_phase_counts), 1L)
  control_phase_fraction <- control_phase_counts /
    pmax(sum(control_phase_counts), 1L)
  aggregate_case <- case_cds / pmax(case_cds + case_uorf, 1L)
  aggregate_control <- control_cds /
    pmax(control_cds + control_uorf, 1L)
  delta <- aggregate_case - aggregate_control
  range_separation <- if (delta < 0) {
    min(control$clean_cds_allocation) - max(case$clean_cds_allocation)
  } else {
    min(case$clean_cds_allocation) - max(control$clean_cds_allocation)
  }
  list(
    n_case_units = nrow(case),
    n_control_units = nrow(control),
    case_union_uorf_counts = case_uorf,
    control_union_uorf_counts = control_uorf,
    case_clean_cds_counts = case_cds,
    control_clean_cds_counts = control_cds,
    case_clean_cds_allocation = aggregate_case,
    control_clean_cds_allocation = aggregate_control,
    clean_cds_allocation_delta = delta,
    biological_unit_range_separation = range_separation,
    case_uorf_phase0_fraction = case_phase_fraction[[1L]],
    case_uorf_phase1_fraction = case_phase_fraction[[2L]],
    case_uorf_phase2_fraction = case_phase_fraction[[3L]],
    control_uorf_phase0_fraction = control_phase_fraction[[1L]],
    control_uorf_phase1_fraction = control_phase_fraction[[2L]],
    control_uorf_phase2_fraction = control_phase_fraction[[3L]],
    case_allocation_min = min(case$clean_cds_allocation),
    case_allocation_max = max(case$clean_cds_allocation),
    control_allocation_min = min(control$clean_cds_allocation),
    control_allocation_max = max(control$clean_cds_allocation)
  )
}, by = category]

shift_summary <- shift_sensitivity[, .(
  coordinate_shift_direction_stable =
    length(unique(sign(clean_cds_allocation_delta))) == 1L,
  coordinate_shift_min_abs_delta = min(abs(clean_cds_allocation_delta))
), by = category]
contrast <- merge(contrast, shift_summary, by = "category", all.x = TRUE,
                  sort = FALSE)
contrast <- merge(
  contrast, position_leaveout_summary,
  by = "category", all.x = TRUE, sort = FALSE
)
contrast[, shared_phase_category := category %chin% c(
  "all_shared_well_phased",
  "short_20_23nt_shared_well_phased",
  "canonical_28_32nt_shared_well_phased"
)]
contrast[, target_uorf_phase_gate :=
  case_uorf_phase0_fraction >= 0.40 &
    case_uorf_phase0_fraction > pmax(
      case_uorf_phase1_fraction, case_uorf_phase2_fraction
    ) &
    control_uorf_phase0_fraction >= 0.40 &
    control_uorf_phase0_fraction > pmax(
      control_uorf_phase1_fraction, control_uorf_phase2_fraction
    )]
contrast[, evidence_gate :=
  shared_phase_category & n_case_units == 2L & n_control_units == 2L &
  case_union_uorf_counts >= 30L & control_union_uorf_counts >= 30L &
  case_clean_cds_counts >= 100L & control_clean_cds_counts >= 100L &
  abs(clean_cds_allocation_delta) >= 0.05 &
  biological_unit_range_separation > 0 &
  target_uorf_phase_gate == TRUE &
  coordinate_shift_direction_stable == TRUE &
  coordinate_shift_min_abs_delta >= 0.03 &
  single_position_direction_stable == TRUE &
  single_position_min_abs_delta >= 0.03 &
  greedy_direction_stable == TRUE &
  abs(greedy_final_clean_cds_allocation_delta) >= 0.03]

position_profile <- profile[, .N, by = .(
  run, condition, replicate, read_length, length_stratum,
  shared_well_phased, tx_position, feature
)]

fwrite(run_table, file.path(output_dir, "run_manifest.csv"))
fwrite(data.table(
  n_calibration_windows = length(calibration_windows),
  shared_well_phased_lengths = paste(shared_well_phased_lengths,
                                     collapse = ";")
), file.path(output_dir, "calibration_design.csv"))
fwrite(offsets, file.path(output_dir, "psite_offsets.csv"))
fwrite(calibration, file.path(output_dir, "psite_offset_qc.csv"))
fwrite(position_profile, file.path(output_dir, "eif4g2_position_profile.csv"))
fwrite(run_metrics, file.path(output_dir, "eif4g2_run_length_metrics.csv"))
fwrite(shift_sensitivity,
       file.path(output_dir, "eif4g2_coordinate_shift_sensitivity.csv"))
fwrite(
  position_leaveout_rows,
  file.path(output_dir, "eif4g2_position_leaveout_sensitivity.csv")
)
fwrite(contrast, file.path(output_dir, "eif4g2_length_contrast_summary.csv"))

message("Shared well-phased lengths: ",
        paste(shared_well_phased_lengths, collapse = ", "))
message("Evidence-gate categories:")
print(contrast[evidence_gate == TRUE])
message("Saved length-resolved EIF4G2 UV validation: ", output_dir)
