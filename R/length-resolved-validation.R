dcs_read_length_stratum <- function(read_length) {
  read_length <- as.integer(read_length)
  out <- rep("other_15_34nt", length(read_length))
  out[read_length >= 20L & read_length <= 23L] <- "short_20_23nt"
  out[read_length >= 28L & read_length <= 32L] <- "canonical_28_32nt"
  out[is.na(read_length) | read_length < 15L | read_length > 34L] <- NA_character_
  out
}

dcs_psite_calibration_qc <- function(relative_five_prime, read_length,
                                     shift_table, downstream = 149L) {
  required <- c("fraction", "offsets_start")
  if (!is.data.frame(shift_table) ||
      !all(required %in% names(shift_table))) {
    stop("shift_table must contain fraction and offsets_start.",
         call. = FALSE)
  }
  if (length(relative_five_prime) != length(read_length)) {
    stop("relative_five_prime and read_length must have equal length.",
         call. = FALSE)
  }
  downstream <- as.integer(downstream)
  if (length(downstream) != 1L || is.na(downstream) || downstream < 2L) {
    stop("downstream must be one integer of at least 2.", call. = FALSE)
  }

  shift_table <- unique(data.frame(
    read_length = as.integer(shift_table$fraction),
    offsets_start = as.integer(shift_table$offsets_start)
  ))
  if (anyDuplicated(shift_table$read_length)) {
    stop("shift_table has duplicated read lengths.", call. = FALSE)
  }
  if (any(!is.finite(shift_table$read_length)) ||
      any(!is.finite(shift_table$offsets_start)) ||
      any(shift_table$read_length <= 0L) ||
      any(shift_table$offsets_start > 0L) ||
      any(-shift_table$offsets_start >= shift_table$read_length)) {
    stop(
      paste0(
        "Offsets must be finite, non-positive ORFik start offsets whose ",
        "P-site displacement is smaller than the footprint length."
      ),
         call. = FALSE)
  }

  read_length <- as.integer(read_length)
  relative_five_prime <- as.integer(relative_five_prime)
  offset <- -shift_table$offsets_start[
    match(read_length, shift_table$read_length)
  ]
  relative_psite <- relative_five_prime + offset
  keep <- is.finite(relative_psite) & is.finite(read_length)
  if (!any(keep)) {
    return(data.frame(
      read_length = integer(), offsets_start = integer(),
      n_start_window = integer(), n_cds_start = integer(),
      phase0_counts = integer(), phase1_counts = integer(),
      phase2_counts = integer(), phase0_fraction = numeric(),
      phase_margin = numeric(), periodicity_gate = logical()
    ))
  }

  input <- data.frame(
    read_length = read_length[keep],
    relative_psite = relative_psite[keep]
  )
  rows <- lapply(sort(unique(input$read_length)), function(length_value) {
    position <- input$relative_psite[input$read_length == length_value]
    position <- position[position >= -30L & position <= downstream]
    cds_position <- position[position >= 0L]
    phase_counts <- vapply(0:2, function(phase) {
      sum(cds_position %% 3L == phase)
    }, integer(1))
    phase_total <- sum(phase_counts)
    phase_fraction <- if (phase_total) phase_counts / phase_total else {
      rep(NA_real_, 3L)
    }
    phase_margin <- if (phase_total) {
      phase_fraction[[1]] - max(phase_fraction[-1L])
    } else {
      NA_real_
    }
    data.frame(
      read_length = length_value,
      offsets_start = shift_table$offsets_start[
        match(length_value, shift_table$read_length)
      ],
      n_start_window = length(position),
      n_cds_start = phase_total,
      phase0_counts = phase_counts[[1]],
      phase1_counts = phase_counts[[2]],
      phase2_counts = phase_counts[[3]],
      phase0_fraction = phase_fraction[[1]],
      phase_margin = phase_margin,
      periodicity_gate = phase_total >= 100L &&
        is.finite(phase_fraction[[1]]) && phase_fraction[[1]] >= 0.40 &&
        phase_margin >= 0.05
    )
  })
  do.call(rbind, rows)
}

dcs_position_leaveout_audit <- function(
    position_counts, case_label = "case", control_label = "control",
    max_greedy_deletions = 3L) {
  required <- c("category", "condition", "region", "position", "counts")
  if (!is.data.frame(position_counts) ||
      !all(required %in% names(position_counts))) {
    stop(
      "position_counts must contain category, condition, region, position, and counts.",
      call. = FALSE
    )
  }
  max_greedy_deletions <- as.integer(max_greedy_deletions)
  if (length(max_greedy_deletions) != 1L ||
      is.na(max_greedy_deletions) || max_greedy_deletions < 1L) {
    stop("max_greedy_deletions must be a positive integer.",
         call. = FALSE)
  }

  x <- data.frame(
    category = as.character(position_counts$category),
    condition = as.character(position_counts$condition),
    region = as.character(position_counts$region),
    position = as.integer(position_counts$position),
    counts = as.numeric(position_counts$counts),
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(x$position)) || any(!is.finite(x$counts)) ||
      any(x$counts < 0) ||
      !all(x$region %in% c("uorf", "clean_CDS"))) {
    stop("Position counts contain invalid positions, counts, or regions.",
         call. = FALSE)
  }
  x <- stats::aggregate(
    counts ~ category + condition + region + position,
    data = x, FUN = sum
  )
  categories <- unique(x$category)
  summary_rows <- vector("list", length(categories))
  leaveout_rows <- vector("list", length(categories))

  for (i in seq_along(categories)) {
    category_value <- categories[[i]]
    z <- x[x$category == category_value, , drop = FALSE]
    observed_conditions <- unique(z$condition)
    if (!all(c(case_label, control_label) %in% observed_conditions)) {
      stop("Every category must contain both case and control conditions.",
           call. = FALSE)
    }
    positions <- sort(unique(z$position))
    allocation_delta <- function(removed = integer()) {
      kept <- z[!z$position %in% removed, , drop = FALSE]
      region_total <- function(condition_value, region_value) {
        sum(kept$counts[
          kept$condition == condition_value & kept$region == region_value
        ])
      }
      case_uorf <- region_total(case_label, "uorf")
      case_cds <- region_total(case_label, "clean_CDS")
      control_uorf <- region_total(control_label, "uorf")
      control_cds <- region_total(control_label, "clean_CDS")
      case_cds / max(case_cds + case_uorf, 1) -
        control_cds / max(control_cds + control_uorf, 1)
    }

    observed_delta <- allocation_delta()
    single_delta <- vapply(
      positions, function(position) allocation_delta(position), numeric(1)
    )
    worst_single_index <- which.min(abs(single_delta))
    single_rows <- data.frame(
      category = category_value,
      deletion_scope = "single_position",
      deletion_step = 1L,
      removed_position = positions,
      clean_cds_allocation_delta = single_delta,
      stringsAsFactors = FALSE
    )

    removed <- integer()
    greedy_rows <- vector("list", min(max_greedy_deletions,
                                      length(positions)))
    for (step in seq_along(greedy_rows)) {
      candidates <- setdiff(positions, removed)
      candidate_delta <- vapply(
        candidates,
        function(position) allocation_delta(c(removed, position)),
        numeric(1)
      )
      chosen_index <- which.min(abs(candidate_delta))
      removed <- c(removed, candidates[[chosen_index]])
      greedy_rows[[step]] <- data.frame(
        category = category_value,
        deletion_scope = "greedy_worst_case",
        deletion_step = step,
        removed_position = candidates[[chosen_index]],
        clean_cds_allocation_delta = candidate_delta[[chosen_index]],
        stringsAsFactors = FALSE
      )
    }
    greedy_table <- do.call(rbind, greedy_rows)
    final_greedy_delta <- utils::tail(
      greedy_table$clean_cds_allocation_delta, 1L
    )
    summary_rows[[i]] <- data.frame(
      category = category_value,
      observed_clean_cds_allocation_delta = observed_delta,
      single_position_direction_stable =
        all(sign(single_delta) == sign(observed_delta)),
      single_position_min_abs_delta = min(abs(single_delta)),
      worst_single_position = positions[[worst_single_index]],
      greedy_deletions = nrow(greedy_table),
      greedy_removed_positions = paste(removed, collapse = ";"),
      greedy_final_clean_cds_allocation_delta = final_greedy_delta,
      greedy_direction_stable =
        sign(final_greedy_delta) == sign(observed_delta),
      stringsAsFactors = FALSE
    )
    leaveout_rows[[i]] <- rbind(single_rows, greedy_table)
  }

  list(
    summary = do.call(rbind, summary_rows),
    leaveout = do.call(rbind, leaveout_rows)
  )
}
