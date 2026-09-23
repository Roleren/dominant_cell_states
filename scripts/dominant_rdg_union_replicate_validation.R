#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(ORFik)
})

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(
      here, "results", "dominant_rdg_merged_replicate_translon_shift",
      "merged_replicate_translon_shift_review_queue.csv"
    ))) return(here)
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find dominant_cell_states analysis root.", call. = FALSE)
}

clean_text <- function(x, missing = "MISSING") {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x)] <- missing
  x
}

same_text <- function(x, value) {
  clean_text(x) == clean_text(value)
}

split_values <- function(x) {
  x <- unlist(strsplit(as.character(x), ";", fixed = TRUE), use.names = FALSE)
  x <- clean_text(x, missing = "")
  unique(x[nzchar(x)])
}

split_runs <- function(x) {
  x <- unlist(strsplit(as.character(x), ";", fixed = TRUE), use.names = FALSE)
  x <- trimws(x)
  unique(x[nzchar(x) & !is.na(x)])
}

positions_from_rows <- function(dt, n_positions) {
  if (!nrow(dt)) return(integer())
  out <- unlist(Map(seq.int, as.integer(dt$tx_start), as.integer(dt$tx_end)),
                use.names = FALSE)
  sort(unique(out[is.finite(out) & out >= 1L & out <= n_positions]))
}

region_metrics <- function(values, positions) {
  values <- as.numeric(values[positions])
  values[!is.finite(values)] <- 0
  total <- sum(values)
  ordered <- sort(values, decreasing = TRUE)
  data.table(
    counts = total,
    bases = length(positions),
    density = if (length(positions)) total / length(positions) else NA_real_,
    max_position_fraction = if (total > 0) max(values) / total else NA_real_,
    top3_position_fraction = if (total > 0) {
      sum(utils::head(ordered, 3L)) / total
    } else NA_real_,
    nonzero_position_fraction = if (length(positions)) {
      mean(values > 0)
    } else NA_real_
  )
}

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) min(x) else NA_real_
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

aggregate_biological_units <- function(mat, run_metadata, group_label) {
  run_metadata <- copy(run_metadata)
  run_metadata[, biological_unit := dcs_biological_unit_ids(
    run_metadata, group_label
  )]
  units <- unique(run_metadata$biological_unit)
  out <- vapply(units, function(unit) {
    runs <- run_metadata[biological_unit == unit, Run]
    rowSums(mat[, runs, drop = FALSE])
  }, numeric(nrow(mat)))
  if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
  colnames(out) <- units
  out
}

unit_metadata <- function(run_metadata, group_label) {
  run_metadata <- copy(run_metadata)
  run_metadata[, biological_unit := dcs_biological_unit_ids(
    run_metadata, group_label
  )]
  run_metadata[, biological_unit_raw := sub(
    paste0("^", group_label, "::"), "", biological_unit
  )]
  run_metadata[, .(
    group = group_label,
    run_count = uniqueN(Run),
    runs = paste(unique(Run), collapse = ";"),
    sample_names = paste(unique(clean_text(SampleName, missing = "")),
                         collapse = ";"),
    replicate_labels = paste(unique(clean_text(REPLICATE, missing = "")),
                             collapse = ";"),
    total_cov_signal = sum(suppressWarnings(as.numeric(Total_Cov_Signal)),
                           na.rm = TRUE),
    frame0 = safe_mean(Frame_usage_0),
    frame1 = safe_mean(Frame_usage_1),
    frame2 = safe_mean(Frame_usage_2),
    tis_peak_strength = safe_mean(tis_peak_strength),
    top_readlength = safe_median(top_readlength)
  ), by = .(biological_unit, biological_unit_raw)]
}

validate_exact_row <- function(row, annotation, group_table, metadata, df,
                               ignored_context_fields = character()) {
  selector_field <- if ("case_selector_field" %in% names(row)) {
    as.character(row$case_selector_field)
  } else {
    "CONDITION"
  }
  selector_case <- if ("case_selector_value" %in% names(row)) {
    as.character(row$case_selector_value)
  } else {
    as.character(row$case_condition)
  }
  selector_control <- if ("control_selector_value" %in% names(row)) {
    split_values(row$control_selector_value)
  } else {
    split_values(row$control_conditions)
  }
  fields <- c("study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE",
              "INHIBITOR", "FRACTION", "Cancer_type", "Sex", "TIMEPOINT")
  fields <- setdiff(fields, c(ignored_context_fields, selector_field))
  groups <- group_table[gene_symbol == row$gene_symbol & tx_id == row$tx_id]
  for (field in fields) {
    groups <- groups[same_text(get(field), row[[field]])]
  }
  case_group <- groups[same_text(get(selector_field), selector_case)]
  control_group <- groups[
    clean_text(get(selector_field)) %chin% clean_text(selector_control)
  ]
  if (!nrow(case_group) || !nrow(control_group)) {
    stop("Could not recover exact case/control groups.")
  }
  case_runs <- unique(unlist(lapply(case_group$example_runs, split_runs),
                             use.names = FALSE))
  control_runs <- unique(unlist(lapply(control_group$example_runs, split_runs),
                                use.names = FALSE))
  all_runs <- unique(c(case_runs, control_runs))
  md <- metadata[Run %chin% all_runs]
  missing_metadata <- setdiff(all_runs, md$Run)
  if (length(missing_metadata)) {
    stop("Missing metadata for runs: ", paste(missing_metadata, collapse = ", "))
  }

  tx_id <- row$tx_id
  mrna <- ORFik::loadRegion(df, "mrna", names.keep = tx_id)
  if (!length(mrna)) stop("Transcript not found: ", tx_id)
  collection <- RiboCrypt:::collection_path_from_exp(df, tx_id,
                                                      grl_all = mrna)
  reads <- RiboCrypt:::load_collection(collection, grl = mrna,
                                       columns = all_runs)
  mat <- as.matrix(reads)
  storage.mode(mat) <- "numeric"
  if (!all(all_runs %in% colnames(mat))) {
    stop("Coverage matrix is missing requested runs.")
  }

  ann <- annotation[gene_symbol == row$gene_symbol & tx_id == row$tx_id]
  uorf_ann <- ann[feature_type %chin% c("leader_uorf", "overlapping_uorf")]
  cds_ann <- ann[feature_type == "clean_CDS"]
  if (!nrow(uorf_ann) || !nrow(cds_ann)) {
    stop("Missing uORF or clean-CDS annotation.")
  }
  uorf_positions <- positions_from_rows(uorf_ann, nrow(mat))
  uorf_phase_positions <- dcs_unambiguous_uorf_phase_positions(
    uorf_ann, nrow(mat)
  )
  cds_all_positions <- positions_from_rows(cds_ann, nrow(mat))
  clean_positions <- setdiff(cds_all_positions, uorf_positions)
  if (!length(uorf_positions) || !length(clean_positions)) {
    stop("Empty union-uORF or unique clean-CDS region.")
  }

  case_md <- md[Run %chin% case_runs]
  control_md <- md[Run %chin% control_runs]
  case_mat <- aggregate_biological_units(mat[, case_runs, drop = FALSE],
                                         case_md, "case")
  control_mat <- aggregate_biological_units(mat[, control_runs, drop = FALSE],
                                            control_md, "control")
  unit_mat <- cbind(case_mat, control_mat)
  unit_md <- rbind(unit_metadata(case_md, "case"),
                   unit_metadata(control_md, "control"), fill = TRUE)

  feature_unit_rows <- rbindlist(lapply(seq_len(nrow(uorf_ann)), function(j) {
    feature <- uorf_ann[j]
    feature_positions <- positions_from_rows(feature, nrow(mat))
    rbindlist(lapply(seq_len(ncol(unit_mat)), function(i) {
      values <- as.numeric(unit_mat[feature_positions, i])
      total <- sum(values)
      phase <- (feature_positions - as.integer(feature$tx_start)) %% 3L
      phase_counts <- vapply(0:2, function(k) sum(values[phase == k]),
                             numeric(1))
      data.table(
        biological_unit = colnames(unit_mat)[i],
        feature_id = feature$feature_id,
        feature_type = feature$feature_type,
        tx_start = feature$tx_start,
        tx_end = feature$tx_end,
        feature_bases = length(feature_positions),
        start_codon = feature$start_codon,
        start_codon_class = feature$start_codon_class,
        kozak_strength = feature$kozak_strength,
        counts = total,
        density = total / pmax(length(feature_positions), 1L),
        phase0_counts = phase_counts[[1]],
        phase1_counts = phase_counts[[2]],
        phase2_counts = phase_counts[[3]],
        phase0_fraction = if (total > 0) phase_counts[[1]] / total else NA_real_,
        start_9nt_fraction = if (total > 0) {
          sum(utils::head(values, 9L)) / total
        } else {
          NA_real_
        },
        stop_6nt_fraction = if (total > 0) {
          sum(utils::tail(values, 6L)) / total
        } else {
          NA_real_
        },
        max_position_fraction = if (total > 0) max(values) / total else NA_real_,
        nonzero_position_fraction = mean(values > 0)
      )
    }), fill = TRUE)
  }), fill = TRUE)
  feature_unit_rows <- merge(
    feature_unit_rows, unit_md, by = "biological_unit", all.x = TRUE,
    sort = FALSE
  )

  unit_rows <- rbindlist(lapply(seq_len(ncol(unit_mat)), function(i) {
    u <- region_metrics(unit_mat[, i], uorf_positions)
    c <- region_metrics(unit_mat[, i], clean_positions)
    phase_counts <- vapply(
      uorf_phase_positions[c("phase0", "phase1", "phase2")],
      function(position) sum(unit_mat[position, i]), numeric(1)
    )
    phase_total <- sum(phase_counts)
    out <- dcs_bind_one_row_metrics(u[, .(
      uorf_counts = counts,
      uorf_bases = bases,
      uorf_density = density,
      uorf_max_position_fraction = max_position_fraction,
      uorf_top3_position_fraction = top3_position_fraction,
      uorf_nonzero_position_fraction = nonzero_position_fraction
    )], c[, .(
      clean_cds_counts = counts,
      clean_cds_bases = bases,
      clean_cds_density = density,
      clean_cds_max_position_fraction = max_position_fraction,
      clean_cds_top3_position_fraction = top3_position_fraction,
      clean_cds_nonzero_position_fraction = nonzero_position_fraction
    )])
    out[, `:=`(
      uorf_phase0_counts = phase_counts[["phase0"]],
      uorf_phase1_counts = phase_counts[["phase1"]],
      uorf_phase2_counts = phase_counts[["phase2"]],
      uorf_phase0_fraction = if (phase_total > 0) {
        phase_counts[["phase0"]] / phase_total
      } else {
        NA_real_
      },
      uorf_phase_informative_fraction = if (u$counts > 0) {
        phase_total / u$counts
      } else {
        NA_real_
      }
    )]
    out[, biological_unit := colnames(unit_mat)[i]]
    out
  }), fill = TRUE)
  unit_rows <- merge(unit_rows, unit_md, by = "biological_unit", all.x = TRUE,
                     sort = FALSE)
  unit_rows[, clean_cds_allocation := clean_cds_counts /
              pmax(clean_cds_counts + uorf_counts, 1)]
  unit_rows[, log2_uorf_to_clean_density := log2(
    (uorf_density + 0.5 / pmax(uorf_bases, 1)) /
      (clean_cds_density + 0.5 / pmax(clean_cds_bases, 1))
  )]
  unit_rows[, clean_cds_cpm := clean_cds_counts /
              pmax(total_cov_signal, 1) * 1e6]
  unit_rows[, uorf_cpm := uorf_counts /
              pmax(total_cov_signal, 1) * 1e6]

  case_run_lookup <- copy(case_md)
  case_run_lookup[, `:=`(
    group = "case",
    biological_unit = dcs_biological_unit_ids(case_md, "case")
  )]
  control_run_lookup <- copy(control_md)
  control_run_lookup[, `:=`(
    group = "control",
    biological_unit = dcs_biological_unit_ids(control_md, "control")
  )]
  run_lookup <- rbind(
    case_run_lookup, control_run_lookup, fill = TRUE
  )[, .(
    Run, group, biological_unit,
    replicate_label = clean_text(REPLICATE, missing = ""),
    sample_name = clean_text(SampleName, missing = ""),
    fraction = clean_text(FRACTION, missing = ""),
    total_cov_signal = suppressWarnings(as.numeric(Total_Cov_Signal)),
    frame0 = suppressWarnings(as.numeric(Frame_usage_0)),
    frame1 = suppressWarnings(as.numeric(Frame_usage_1)),
    frame2 = suppressWarnings(as.numeric(Frame_usage_2)),
    top_readlength = suppressWarnings(as.numeric(top_readlength))
  )]
  run_rows <- rbindlist(lapply(all_runs, function(run) {
    values <- mat[, run]
    u <- region_metrics(values, uorf_positions)
    c <- region_metrics(values, clean_positions)
    phase_counts <- vapply(
      uorf_phase_positions[c("phase0", "phase1", "phase2")],
      function(position) sum(values[position]), numeric(1)
    )
    phase_total <- sum(phase_counts)
    out <- dcs_bind_one_row_metrics(u[, .(
      uorf_counts = counts,
      uorf_bases = bases,
      uorf_density = density,
      uorf_max_position_fraction = max_position_fraction,
      uorf_top3_position_fraction = top3_position_fraction,
      uorf_nonzero_position_fraction = nonzero_position_fraction
    )], c[, .(
      clean_cds_counts = counts,
      clean_cds_bases = bases,
      clean_cds_density = density,
      clean_cds_max_position_fraction = max_position_fraction,
      clean_cds_top3_position_fraction = top3_position_fraction,
      clean_cds_nonzero_position_fraction = nonzero_position_fraction
    )])
    out[, `:=`(
      Run = run,
      uorf_phase0_counts = phase_counts[["phase0"]],
      uorf_phase1_counts = phase_counts[["phase1"]],
      uorf_phase2_counts = phase_counts[["phase2"]],
      uorf_phase0_fraction = if (phase_total > 0) {
        phase_counts[["phase0"]] / phase_total
      } else {
        NA_real_
      }
    )]
    out
  }), fill = TRUE)
  run_rows <- merge(
    run_rows, run_lookup, by = "Run", all.x = TRUE, sort = FALSE
  )
  run_rows[, clean_cds_allocation := clean_cds_counts /
             pmax(clean_cds_counts + uorf_counts, 1)]
  run_rows[, clean_cds_cpm := clean_cds_counts /
             pmax(total_cov_signal, 1) * 1e6]
  run_rows[, uorf_cpm := uorf_counts /
             pmax(total_cov_signal, 1) * 1e6]

  aggregate_delta_for_shift <- function(shift) {
    up <- uorf_positions + shift
    cp <- clean_positions + shift
    up <- up[up >= 1L & up <= nrow(mat)]
    cp <- cp[cp >= 1L & cp <= nrow(mat)]
    case_u <- sum(case_mat[up, , drop = FALSE])
    case_c <- sum(case_mat[cp, , drop = FALSE])
    control_u <- sum(control_mat[up, , drop = FALSE])
    control_c <- sum(control_mat[cp, , drop = FALSE])
    data.table(
      coordinate_shift = shift,
      clean_cds_allocation_delta =
        case_c / pmax(case_c + case_u, 1) -
        control_c / pmax(control_c + control_u, 1),
      log2_uorf_to_clean_density_delta = log2(
        ((case_u / length(up)) + 0.5 / length(up)) /
          ((case_c / length(cp)) + 0.5 / length(cp))
      ) - log2(
        ((control_u / length(up)) + 0.5 / length(up)) /
          ((control_c / length(cp)) + 0.5 / length(cp))
      )
    )
  }
  shift_rows <- rbindlist(lapply(-2:2, aggregate_delta_for_shift))
  observed <- shift_rows[coordinate_shift == 0]
  direction <- sign(observed$clean_cds_allocation_delta)

  case_units <- unit_rows[group == "case"]
  control_units <- unit_rows[group == "control"]
  case_run_rows <- run_rows[group == "case"]
  control_run_rows <- run_rows[group == "control"]
  median_delta <- safe_median(case_units$clean_cds_allocation) -
    safe_median(control_units$clean_cds_allocation)
  range_separation <- if (is.finite(median_delta) && median_delta >= 0) {
    safe_min(case_units$clean_cds_allocation) -
      safe_max(control_units$clean_cds_allocation)
  } else {
    safe_min(control_units$clean_cds_allocation) -
    safe_max(case_units$clean_cds_allocation)
  }
  run_level_range_separation <- if (
    is.finite(median_delta) && median_delta >= 0
  ) {
    safe_min(case_run_rows$clean_cds_allocation) -
      safe_max(control_run_rows$clean_cds_allocation)
  } else {
    safe_min(control_run_rows$clean_cds_allocation) -
      safe_max(case_run_rows$clean_cds_allocation)
  }
  case_aggregate_uorf <- sum(case_mat[uorf_positions, , drop = FALSE])
  case_aggregate_clean <- sum(case_mat[clean_positions, , drop = FALSE])
  control_aggregate_uorf <- sum(control_mat[uorf_positions, , drop = FALSE])
  control_aggregate_clean <- sum(control_mat[clean_positions, , drop = FALSE])

  allocation_delta_without <- function(removed_positions) {
    retained <- setdiff(uorf_positions, removed_positions)
    case_uorf <- sum(case_mat[retained, , drop = FALSE])
    control_uorf <- sum(control_mat[retained, , drop = FALSE])
    case_aggregate_clean / pmax(case_aggregate_clean + case_uorf, 1) -
      control_aggregate_clean /
        pmax(control_aggregate_clean + control_uorf, 1)
  }
  leaveout_rows <- rbindlist(lapply(uorf_positions, function(position) {
    data.table(
      removal_type = "single_position",
      removed_positions = as.character(position),
      n_removed_positions = 1L,
      removed_pooled_counts = sum(unit_mat[position, , drop = FALSE]),
      clean_cds_allocation_delta = allocation_delta_without(position)
    )
  }))
  pooled_uorf_counts <- rowSums(unit_mat[uorf_positions, , drop = FALSE])
  top3_positions <- uorf_positions[
    head(order(pooled_uorf_counts, decreasing = TRUE), 3L)
  ]
  leaveout_rows <- rbind(
    leaveout_rows,
    data.table(
      removal_type = "pooled_top3_positions",
      removed_positions = paste(top3_positions, collapse = ";"),
      n_removed_positions = length(top3_positions),
      removed_pooled_counts = sum(
        unit_mat[top3_positions, , drop = FALSE]
      ),
      clean_cds_allocation_delta = allocation_delta_without(top3_positions)
    ),
    fill = TRUE
  )
  single_leaveout <- leaveout_rows[removal_type == "single_position"]
  top3_leaveout <- leaveout_rows[removal_type == "pooled_top3_positions"]

  summary <- data.table(
    gene_symbol = row$gene_symbol,
    tx_id = tx_id,
    design_family = row$design_family,
    study = row$study,
    author = row$AUTHOR,
    cell_line = row$CELL_LINE,
    tissue = row$TISSUE,
    case_condition = row$case_condition,
    control_conditions = row$control_conditions,
    timepoint = row$TIMEPOINT,
    inhibitor = row$INHIBITOR,
    fraction = if ("FRACTION" %in% ignored_context_fields) {
      paste(sort(unique(clean_text(groups$FRACTION))), collapse = ";")
    } else {
      row$FRACTION
    },
    n_case_runs = length(case_runs),
    n_control_runs = length(control_runs),
    n_case_biological_units = ncol(case_mat),
    n_control_biological_units = ncol(control_mat),
    n_uorf_features = nrow(uorf_ann),
    summed_uorf_feature_bases = sum(uorf_ann$feature_bases, na.rm = TRUE),
    union_uorf_bases = length(uorf_positions),
    union_uorf_phase_informative_bases =
      length(uorf_phase_positions$informative),
    overlap_inflation_factor = sum(uorf_ann$feature_bases, na.rm = TRUE) /
      length(uorf_positions),
    unique_clean_cds_bases = length(clean_positions),
    case_union_uorf_counts = case_aggregate_uorf,
    control_union_uorf_counts = control_aggregate_uorf,
    case_unique_clean_cds_counts = case_aggregate_clean,
    control_unique_clean_cds_counts = control_aggregate_clean,
    union_clean_cds_allocation_delta =
      observed$clean_cds_allocation_delta,
    median_biological_unit_allocation_delta = median_delta,
    biological_unit_range_separation = range_separation,
    run_level_range_separation = run_level_range_separation,
    run_level_direction_concordant = run_level_range_separation > 0,
    case_allocation_min = safe_min(case_units$clean_cds_allocation),
    case_allocation_max = safe_max(case_units$clean_cds_allocation),
    control_allocation_min = safe_min(control_units$clean_cds_allocation),
    control_allocation_max = safe_max(control_units$clean_cds_allocation),
    log2_uorf_to_clean_density_delta =
      observed$log2_uorf_to_clean_density_delta,
    clean_cds_cpm_log2_ratio = log2(
      (safe_mean(case_units$clean_cds_cpm) + 0.5) /
        (safe_mean(control_units$clean_cds_cpm) + 0.5)
    ),
    uorf_cpm_log2_ratio = log2(
      (safe_mean(case_units$uorf_cpm) + 0.5) /
        (safe_mean(control_units$uorf_cpm) + 0.5)
    ),
    max_uorf_position_fraction = safe_max(unit_rows$uorf_max_position_fraction),
    max_uorf_top3_fraction = safe_max(unit_rows$uorf_top3_position_fraction),
    min_uorf_nonzero_fraction = safe_min(unit_rows$uorf_nonzero_position_fraction),
    case_uorf_phase0_fraction =
      sum(case_units$uorf_phase0_counts) /
      pmax(sum(case_units$uorf_phase0_counts +
                 case_units$uorf_phase1_counts +
                 case_units$uorf_phase2_counts), 1),
    control_uorf_phase0_fraction =
      sum(control_units$uorf_phase0_counts) /
      pmax(sum(control_units$uorf_phase0_counts +
                 control_units$uorf_phase1_counts +
                 control_units$uorf_phase2_counts), 1),
    min_uorf_phase_informative_fraction =
      safe_min(unit_rows$uorf_phase_informative_fraction),
    case_frame0_mean = safe_mean(case_units$frame0),
    control_frame0_mean = safe_mean(control_units$frame0),
    frame0_mean_delta = safe_mean(case_units$frame0) -
      safe_mean(control_units$frame0),
    case_tis_peak_mean = safe_mean(case_units$tis_peak_strength),
    control_tis_peak_mean = safe_mean(control_units$tis_peak_strength),
    coordinate_shift_direction_stable = all(
      sign(shift_rows$clean_cds_allocation_delta) == direction
    ),
    coordinate_shift_min_abs_delta = min(
      abs(shift_rows$clean_cds_allocation_delta), na.rm = TRUE
    ),
    leave_one_uorf_position_direction_stable = all(
      sign(single_leaveout$clean_cds_allocation_delta) == direction
    ),
    leave_one_uorf_position_min_abs_delta = min(
      abs(single_leaveout$clean_cds_allocation_delta), na.rm = TRUE
    ),
    pooled_top3_uorf_positions = top3_leaveout$removed_positions,
    pooled_top3_removed_allocation_delta =
      top3_leaveout$clean_cds_allocation_delta,
    source_review_queue_rank = row$review_queue_rank,
    source_joint_clean_cds_delta = row$joint_delta_clean_CDS,
    source_discovery_priority = row$discovery_priority_score
  )
  summary[, biological_unit_gate :=
            n_case_biological_units >= 2 & n_control_biological_units >= 2]
  summary[, count_gate := case_union_uorf_counts >= 30 &
            control_union_uorf_counts >= 30 &
            case_unique_clean_cds_counts >= 100 &
            control_unique_clean_cds_counts >= 100]
  summary[, shape_gate := max_uorf_position_fraction <= 0.35 &
            max_uorf_top3_fraction <= 0.65 &
            min_uorf_nonzero_fraction >= 0.10]
  summary[, uorf_phase_gate := case_uorf_phase0_fraction >= 0.40 &
            control_uorf_phase0_fraction >= 0.40 &
            min_uorf_phase_informative_fraction >= 0.50]
  summary[, frame_gate := abs(frame0_mean_delta) <= 15 &
            pmin(case_frame0_mean, control_frame0_mean) >= 45]
  summary[, direction_gate := coordinate_shift_direction_stable &
            coordinate_shift_min_abs_delta >= 0.03]
  summary[, peak_robustness_gate :=
            leave_one_uorf_position_direction_stable &
            leave_one_uorf_position_min_abs_delta >= 0.03 &
            sign(pooled_top3_removed_allocation_delta) ==
              sign(union_clean_cds_allocation_delta) &
            abs(pooled_top3_removed_allocation_delta) >= 0.03]
  summary[, exact_union_validation_class := fcase(
    biological_unit_gate & count_gate & shape_gate & uorf_phase_gate & frame_gate &
      direction_gate & biological_unit_range_separation > 0,
    "A_union_validated_range_separated",
    biological_unit_gate & count_gate & shape_gate & uorf_phase_gate & frame_gate &
      direction_gate,
    "B_union_validated_replicate_overlap",
    biological_unit_gate & count_gate & !shape_gate & uorf_phase_gate &
      frame_gate & direction_gate & peak_robustness_gate &
      biological_unit_range_separation > 0,
    "C_shape_flag_leaveout_robust",
    !biological_unit_gate, "D_pseudoreplicate_or_low_biological_n",
    !count_gate, "D_low_union_counts",
    !shape_gate, "D_spike_or_narrow_union_signal",
    !uorf_phase_gate, "D_uorf_phase_inconsistent",
    !frame_gate, "D_frame_imbalance_risk",
    !direction_gate, "D_coordinate_shift_fragile",
    default = "C_review"
  )]

  unit_rows[, `:=`(
    gene_symbol = row$gene_symbol,
    tx_id = tx_id,
    study = row$study,
    case_condition = row$case_condition,
    control_conditions = row$control_conditions,
    timepoint = row$TIMEPOINT,
    source_review_queue_rank = row$review_queue_rank
  )]
  run_rows[, `:=`(
    gene_symbol = row$gene_symbol,
    tx_id = tx_id,
    study = row$study,
    case_condition = row$case_condition,
    control_conditions = row$control_conditions,
    timepoint = row$TIMEPOINT,
    source_review_queue_rank = row$review_queue_rank
  )]
  shift_rows[, `:=`(
    gene_symbol = row$gene_symbol,
    tx_id = tx_id,
    study = row$study,
    case_condition = row$case_condition,
    control_conditions = row$control_conditions,
    timepoint = row$TIMEPOINT,
    source_review_queue_rank = row$review_queue_rank
  )]
  feature_unit_rows[, `:=`(
    gene_symbol = row$gene_symbol,
    tx_id = tx_id,
    study = row$study,
    case_condition = row$case_condition,
    control_conditions = row$control_conditions,
    timepoint = row$TIMEPOINT,
    source_review_queue_rank = row$review_queue_rank
  )]
  leaveout_rows[, `:=`(
    gene_symbol = row$gene_symbol,
    tx_id = tx_id,
    study = row$study,
    case_condition = row$case_condition,
    control_conditions = row$control_conditions,
    timepoint = row$TIMEPOINT,
    source_review_queue_rank = row$review_queue_rank
  )]
  list(
    summary = summary, units = unit_rows, runs = run_rows,
    features = feature_unit_rows,
    shifts = shift_rows, leaveouts = leaveout_rows
  )
}

analysis_dir <- find_analysis_dir()
results_dir <- file.path(analysis_dir, "results")
output_dir <- Sys.getenv(
  "DOMINANT_UNION_VALIDATION_OUTPUT_DIR",
  file.path(results_dir, "dominant_rdg_union_replicate_validation")
)
if (!grepl("^/", output_dir)) output_dir <- file.path(analysis_dir, output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(analysis_dir, "R", "union-validation.R"), local = .GlobalEnv)

ribocrypt_repo <- Sys.getenv("RIBOCRYPT_REPO", unset = "")
if (!nzchar(ribocrypt_repo)) {
  stop("Set RIBOCRYPT_REPO to the local RiboCrypt checkout.", call. = FALSE)
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required.", call. = FALSE)
}
devtools::load_all(ribocrypt_repo, quiet = TRUE)

review_file <- file.path(
  results_dir, "dominant_rdg_merged_replicate_translon_shift",
  "merged_replicate_translon_shift_review_queue.csv"
)
contrast_file <- file.path(
  results_dir, "dominant_rdg_grouped_branch_usage",
  "dominant_rdg_grouped_branch_usage_contrast_table.csv"
)
adapter_file <- file.path(
  results_dir, "dominant_rdg_grouped_branch_usage",
  "dominant_rdg_grouped_branch_usage_adapter.csv"
)
annotation_file <- file.path(
  results_dir, "dominant_rdg_outputs", "rdg_feature_annotation.csv"
)
metadata_file <- Sys.getenv(
  "DOMINANT_METADATA_FILE",
  "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv"
)

validation_source <- Sys.getenv("DOMINANT_UNION_VALIDATION_SOURCE", "review")
if (identical(validation_source, "review")) {
  review <- fread(review_file)
  review <- review[match_scope == "exact_context"]
} else if (identical(validation_source, "contrast")) {
  review <- fread(contrast_file)
  review <- review[
    match_scope == "exact_context" & feature_type == "clean_CDS"
  ]
  review[, `:=`(
    source_abs_delta = abs(grouped_usage_delta),
    review_queue_rank = frank(-abs(grouped_usage_delta), ties.method = "first"),
    joint_delta_clean_CDS = grouped_usage_delta,
    discovery_priority_score = abs(grouped_usage_delta)
  )]
  setorder(review, -source_abs_delta)
} else {
  stop("DOMINANT_UNION_VALIDATION_SOURCE must be 'review' or 'contrast'.",
       call. = FALSE)
}
requested_genes <- split_values(Sys.getenv("DOMINANT_UNION_VALIDATION_GENES", ""))
if (length(requested_genes)) review <- review[gene_symbol %chin% requested_genes]
requested_studies <- split_values(Sys.getenv(
  "DOMINANT_UNION_VALIDATION_STUDIES", ""
))
if (length(requested_studies)) review <- review[study %chin% requested_studies]
manual_context_field <- Sys.getenv(
  "DOMINANT_UNION_VALIDATION_MANUAL_CONTEXT_FIELD", ""
)
manual_case_value <- Sys.getenv(
  "DOMINANT_UNION_VALIDATION_MANUAL_CASE_VALUE", ""
)
manual_control_value <- Sys.getenv(
  "DOMINANT_UNION_VALIDATION_MANUAL_CONTROL_VALUE", ""
)
ignored_context_fields <- split_values(Sys.getenv(
  "DOMINANT_UNION_VALIDATION_IGNORE_CONTEXT_FIELDS", ""
))
allowed_context_fields <- c(
  "study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE", "CONDITION", "INHIBITOR",
  "FRACTION", "Cancer_type", "Sex", "TIMEPOINT"
)
if (nzchar(manual_context_field)) {
  if (!manual_context_field %chin% allowed_context_fields ||
      !nzchar(manual_case_value) || !nzchar(manual_control_value)) {
    stop(
      "Manual context validation requires a valid context field plus ",
      "non-empty case and control values.", call. = FALSE
    )
  }
  if (!length(requested_genes) || !length(requested_studies)) {
    stop(
      "Manual context validation requires requested gene(s) and study/studies.",
      call. = FALSE
    )
  }
  ignored_context_fields <- unique(c(
    ignored_context_fields, manual_context_field
  ))
}
unknown_ignored_fields <- setdiff(ignored_context_fields,
                                  allowed_context_fields)
if (length(unknown_ignored_fields)) {
  stop(
    "Unknown ignored context fields: ",
    paste(unknown_ignored_fields, collapse = ", "), call. = FALSE
  )
}
if (length(ignored_context_fields)) {
  dedup_fields <- c(
    "gene_symbol", "tx_id", "case_condition", "control_conditions",
    "design_family", setdiff(allowed_context_fields, ignored_context_fields)
  )
  dedup_fields <- intersect(dedup_fields, names(review))
  review <- unique(review, by = dedup_fields)
}
max_rows <- suppressWarnings(as.integer(Sys.getenv(
  "DOMINANT_UNION_VALIDATION_MAX_ROWS", "0"
)))
if (is.finite(max_rows) && max_rows > 0) review <- head(review, max_rows)

adapter_columns <- c(
  "group_id", "gene_symbol", "tx_id", "study", "AUTHOR", "CELL_LINE",
  "TISSUE", "GENE", "CONDITION", "INHIBITOR", "FRACTION", "Cancer_type",
  "Sex", "TIMEPOINT", "n_runs_merged", "example_runs"
)
adapter <- fread(adapter_file, select = adapter_columns)
group_table <- unique(adapter)
if (nzchar(manual_context_field)) {
  manual_groups <- group_table[
    gene_symbol %chin% requested_genes & study %chin% requested_studies
  ]
  manual_groups <- manual_groups[
    clean_text(get(manual_context_field)) %chin%
      clean_text(c(manual_case_value, manual_control_value))
  ]
  identity_fields <- dcs_manual_identity_fields(
    c(
      "gene_symbol", "tx_id", "study", "AUTHOR", "CELL_LINE", "TISSUE",
      "GENE", "INHIBITOR", "FRACTION", "Cancer_type", "Sex", "TIMEPOINT"
    ),
    manual_context_field,
    ignored_context_fields
  )
  identity_fields <- intersect(identity_fields, names(manual_groups))
  identities <- unique(manual_groups[, ..identity_fields])
  if (nrow(identities) != 1L) {
    stop(
      "Manual context selector did not resolve to one exact transcript/context ",
      "after excluding the selected axis (found ", nrow(identities), ").",
      call. = FALSE
    )
  }
  observed_values <- unique(clean_text(manual_groups[[manual_context_field]]))
  required_values <- clean_text(c(manual_case_value, manual_control_value))
  if (!all(required_values %chin% observed_values)) {
    stop("Manual case/control context values were not both recovered.",
         call. = FALSE)
  }
  seed <- manual_groups[1]
  case_label <- Sys.getenv(
    "DOMINANT_UNION_VALIDATION_MANUAL_CASE_LABEL",
    paste0(manual_context_field, "=", manual_case_value)
  )
  control_label <- Sys.getenv(
    "DOMINANT_UNION_VALIDATION_MANUAL_CONTROL_LABEL",
    paste0(manual_context_field, "=", manual_control_value)
  )
  review <- data.table(
    gene_symbol = seed$gene_symbol,
    tx_id = seed$tx_id,
    design_family = paste0("manual_", tolower(manual_context_field), "_axis"),
    study = seed$study,
    AUTHOR = seed$AUTHOR,
    CELL_LINE = seed$CELL_LINE,
    TISSUE = seed$TISSUE,
    GENE = seed$GENE,
    INHIBITOR = seed$INHIBITOR,
    FRACTION = seed$FRACTION,
    Cancer_type = seed$Cancer_type,
    Sex = seed$Sex,
    TIMEPOINT = seed$TIMEPOINT,
    case_condition = case_label,
    control_conditions = control_label,
    case_selector_field = manual_context_field,
    case_selector_value = manual_case_value,
    control_selector_value = manual_control_value,
    review_queue_rank = 1L,
    joint_delta_clean_CDS = NA_real_,
    discovery_priority_score = NA_real_
  )
}
if (!nrow(review)) stop("No exact-context rows selected.", call. = FALSE)
annotation <- fread(annotation_file)
metadata <- fread(metadata_file)
needed_metadata <- c(
  "Run", "Experiment", "BioSample", "SampleName", "REPLICATE",
  "Total_Cov_Signal", "Frame_usage_0", "Frame_usage_1", "Frame_usage_2",
  "tis_peak_strength", "top_readlength"
)
missing_metadata <- setdiff(needed_metadata, names(metadata))
if (length(missing_metadata)) {
  stop("Metadata missing columns: ", paste(missing_metadata, collapse = ", "),
       call. = FALSE)
}

df <- ORFik::read.experiment("all_samples-Homo_sapiens", validate = FALSE)
summaries <- list()
units <- list()
runs <- list()
features <- list()
shifts <- list()
leaveouts <- list()
errors <- list()
for (i in seq_len(nrow(review))) {
  row <- review[i]
  message(sprintf(
    "[%d/%d] %s | %s | %s | %s vs %s",
    i, nrow(review), row$gene_symbol, row$study, row$TIMEPOINT,
    row$case_condition, row$control_conditions
  ))
  result <- tryCatch(
    validate_exact_row(
      row, annotation, group_table, metadata, df,
      ignored_context_fields = ignored_context_fields
    ),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    errors[[length(errors) + 1L]] <- data.table(
      gene_symbol = row$gene_symbol,
      tx_id = row$tx_id,
      study = row$study,
      case_condition = row$case_condition,
      control_conditions = row$control_conditions,
      timepoint = row$TIMEPOINT,
      source_review_queue_rank = row$review_queue_rank,
      error = conditionMessage(result)
    )
  } else {
    summaries[[length(summaries) + 1L]] <- result$summary
    units[[length(units) + 1L]] <- result$units
    runs[[length(runs) + 1L]] <- result$runs
    features[[length(features) + 1L]] <- result$features
    shifts[[length(shifts) + 1L]] <- result$shifts
    leaveouts[[length(leaveouts) + 1L]] <- result$leaveouts
  }
}

summary_dt <- rbindlist(summaries, fill = TRUE)
unit_dt <- rbindlist(units, fill = TRUE)
run_dt <- rbindlist(runs, fill = TRUE)
feature_dt <- rbindlist(features, fill = TRUE)
shift_dt <- rbindlist(shifts, fill = TRUE)
leaveout_dt <- rbindlist(leaveouts, fill = TRUE)
error_dt <- rbindlist(errors, fill = TRUE)
if (!ncol(error_dt)) {
  error_dt <- data.table(
    gene_symbol = character(), tx_id = character(), study = character(),
    case_condition = character(), control_conditions = character(),
    timepoint = character(), source_review_queue_rank = integer(),
    error = character()
  )
}
if (nrow(summary_dt)) {
  summary_dt[, union_abs_allocation_delta :=
               abs(union_clean_cds_allocation_delta)]
  setorder(summary_dt, exact_union_validation_class,
           -union_abs_allocation_delta,
           source_review_queue_rank)
  summary_dt[, union_validation_rank := seq_len(.N)]
  setcolorder(summary_dt, c(
    "union_validation_rank", "exact_union_validation_class",
    setdiff(names(summary_dt),
            c("union_validation_rank", "exact_union_validation_class"))
  ))
}

fwrite(summary_dt, file.path(output_dir, "union_exact_contrast_summary.csv"))
fwrite(unit_dt, file.path(output_dir, "union_biological_unit_metrics.csv"))
fwrite(run_dt, file.path(output_dir, "union_run_metrics.csv"))
fwrite(feature_dt,
       file.path(output_dir, "uorf_feature_biological_unit_metrics.csv"))
fwrite(shift_dt, file.path(output_dir, "union_coordinate_shift_sensitivity.csv"))
fwrite(
  leaveout_dt,
  file.path(output_dir, "union_position_leaveout_sensitivity.csv")
)
fwrite(error_dt, file.path(output_dir, "union_validation_errors.csv"))

message("Validated exact rows: ", nrow(summary_dt), "; errors: ", nrow(error_dt))
message("Saved: ", output_dir)
