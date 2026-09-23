dcs_bind_one_row_metrics <- function(left, right) {
  if (!is.data.frame(left) || !is.data.frame(right) ||
      nrow(left) != 1L || nrow(right) != 1L) {
    stop("Region metric inputs must each be a one-row data frame.",
         call. = FALSE)
  }
  duplicate_names <- intersect(names(left), names(right))
  if (length(duplicate_names)) {
    stop("Region metric names must be unique across inputs: ",
         paste(duplicate_names, collapse = ", "), call. = FALSE)
  }
  cbind(left, right)
}

dcs_manual_identity_fields <- function(candidate_fields,
                                       manual_context_field,
                                       ignored_context_fields = character()) {
  candidate_fields <- unique(as.character(candidate_fields))
  manual_context_field <- as.character(manual_context_field)
  ignored_context_fields <- unique(as.character(ignored_context_fields))
  if (length(manual_context_field) != 1L ||
      is.na(manual_context_field) || !nzchar(manual_context_field)) {
    stop("manual_context_field must be one non-empty string.",
         call. = FALSE)
  }
  candidate_fields[
    !candidate_fields %in% c(manual_context_field, ignored_context_fields)
  ]
}

dcs_biological_unit_ids <- function(run_metadata, group_label) {
  required <- c(
    "Run", "Experiment", "BioSample", "SampleName", "REPLICATE"
  )
  if (!is.data.frame(run_metadata) ||
      !all(required %in% names(run_metadata))) {
    stop(
      "run metadata must contain Run, Experiment, BioSample, SampleName, ",
      "and REPLICATE.", call. = FALSE
    )
  }
  group_label <- as.character(group_label)
  if (length(group_label) != 1L || is.na(group_label) || !nzchar(group_label)) {
    stop("group_label must be one non-empty string.", call. = FALSE)
  }

  clean_identifier <- function(x) {
    x <- trimws(as.character(x))
    missing <- is.na(x) | !nzchar(x) |
      toupper(x) %in% c("MISSING", "NA", "N/A", "UNKNOWN", "NONE")
    x[missing] <- ""
    x
  }

  replicate_id <- clean_identifier(run_metadata$REPLICATE)
  fallback <- clean_identifier(run_metadata$Experiment)
  missing <- !nzchar(fallback)
  fallback[missing] <- clean_identifier(run_metadata$BioSample)[missing]
  missing <- !nzchar(fallback)
  fallback[missing] <- clean_identifier(run_metadata$SampleName)[missing]
  missing <- !nzchar(fallback)
  fallback[missing] <- clean_identifier(run_metadata$Run)[missing]

  # An explicit replicate label is the most conservative biological-unit key.
  # Some SRA submissions assign each sequencing lane a separate Experiment and
  # BioSample even though all lanes share the same biological replicate.
  unit <- ifelse(
    nzchar(replicate_id), paste0("replicate_", replicate_id), fallback
  )
  paste(group_label, unit, sep = "::")
}

dcs_unambiguous_uorf_phase_positions <- function(annotation, n_positions) {
  required <- c("tx_start", "tx_end")
  if (!is.data.frame(annotation) || !all(required %in% names(annotation))) {
    stop("uORF annotation must contain tx_start and tx_end.", call. = FALSE)
  }
  n_positions <- as.integer(n_positions)
  if (length(n_positions) != 1L || is.na(n_positions) || n_positions < 1L) {
    stop("n_positions must be one positive integer.", call. = FALSE)
  }
  if (!nrow(annotation)) {
    return(list(phase0 = integer(), phase1 = integer(), phase2 = integer(),
                informative = integer()))
  }

  phase_rows <- do.call(rbind, lapply(seq_len(nrow(annotation)), function(i) {
    start <- as.integer(annotation$tx_start[[i]])
    end <- as.integer(annotation$tx_end[[i]])
    if (is.na(start) || is.na(end) || end < start) return(NULL)
    position <- seq.int(start, end)
    position <- position[position >= 1L & position <= n_positions]
    if (!length(position)) return(NULL)
    data.frame(position = position, phase = (position - start) %% 3L)
  }))
  if (is.null(phase_rows) || !nrow(phase_rows)) {
    return(list(phase0 = integer(), phase1 = integer(), phase2 = integer(),
                informative = integer()))
  }

  phases <- split(phase_rows$phase, phase_rows$position)
  agreed <- vapply(phases, function(x) {
    x <- unique(x)
    if (length(x) == 1L) x else NA_integer_
  }, integer(1))
  positions <- as.integer(names(agreed))
  informative <- !is.na(agreed)
  positions <- positions[informative]
  agreed <- agreed[informative]
  list(
    phase0 = sort(positions[agreed == 0L]),
    phase1 = sort(positions[agreed == 1L]),
    phase2 = sort(positions[agreed == 2L]),
    informative = sort(positions)
  )
}
