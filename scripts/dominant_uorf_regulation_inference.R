dominant_loader <- c(
  file.path("scripts", "dominant_ribocrypt_loader.R"),
  file.path("dominant_cell_states", "scripts", "dominant_ribocrypt_loader.R")
)
dominant_loader <- dominant_loader[file.exists(dominant_loader)][1]
if (!is.na(dominant_loader)) {
  source(dominant_loader)
  dominant_load_ribocrypt_if_available()
}

library(data.table)
library(GenomicRanges)
library(IRanges)

analysis_dir <- if (file.exists("dominant_cell_states/human_dominant_cell_states_clean_cds_gene_diagnostics.csv")) {
  "dominant_cell_states"
} else {
  "."
}

metadata_file <- "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv"
diagnostics_file <- file.path(analysis_dir, "human_dominant_cell_states_clean_cds_gene_diagnostics.csv")
feature_output <- file.path(analysis_dir, "dominant_uorf_feature_expression.csv")
sample_model_output <- file.path(analysis_dir, "dominant_uorf_regulation_sample_model.csv")
term_model_output <- file.path(analysis_dir, "dominant_uorf_regulation_term_model.csv")
edge_output <- file.path(analysis_dir, "dominant_uorf_regulation_edges.csv")
allocation_model_output <- file.path(analysis_dir, "dominant_uorf_regulation_allocation_model.csv")
metadata_fields <- c("CELL_LINE", "TISSUE", "GENE", "CONDITION")
min_term_n <- 3L

coverage_cache_file <- file.path(analysis_dir, "scripts", "dominant_coverage_cache.R")
if (!file.exists(coverage_cache_file)) {
  stop("Missing dominant coverage cache helper: ", coverage_cache_file)
}
source(coverage_cache_file)
internal_orf_signal_helper <- file.path(
  analysis_dir,
  "scripts",
  "dominant_internal_orf_signal_helpers.R"
)
if (!file.exists(internal_orf_signal_helper)) {
  stop("Missing internal ORF signal helper: ", internal_orf_signal_helper)
}
source(internal_orf_signal_helper)

clean_level <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "NA", "N/A", "na", "n/a", "NULL", "null", "None", "none")] <- NA_character_
  x
}

range_positions <- function(grl) {
  gr <- unlist(grl, use.names = FALSE)
  if (length(gr) == 0) return(integer())
  unique(unlist(Map(seq.int, start(gr), end(gr)), use.names = FALSE))
}

safe_pmap_positions <- function(query, subject) {
  if (length(query) == 0 || length(subject) == 0) return(integer())
  mapped <- tryCatch(
    suppressWarnings(pmapToTranscriptF(query, subject)),
    error = function(e) GRangesList()
  )
  range_positions(mapped)
}

make_display_region <- function(leader_tx, cds_tx, tx_id) {
  parts <- if (length(leader_tx) > 0 && widthPerGroup(leader_tx) > 0) {
    c(leader_tx, cds_tx)
  } else {
    cds_tx
  }
  display <- GRangesList(reduce(unlistGrl(parts)))
  display <- sortPerGroup(display, quick.rev = TRUE)
  names(display) <- tx_id
  display
}

library_sizes_by_run <- function(df) {
  lib_sizes <- readRDS(get_lib_sizes_file(df))
  if (length(lib_sizes) == nrow(df)) names(lib_sizes) <- runIDs(df)
  lib_sizes
}

load_translon_source <- function(paths, source) {
  path <- unique(paths)
  path <- path[file.exists(path)][1]
  if (is.na(path)) return(GRangesList())
  translons <- read_RDSQS(path)
  mcols(translons)$translon_source <- source
  mcols(translons)$translon_source_path <- path
  translons
}

load_predicted_translons <- function(df) {
  t_paths <- c(
    file.path(refFolder(df), "predicted_translons",
              "predicted_translons_with_sequence_ranges.rds"),
    file.path(dirname(df@fafile), "predicted_translons",
              "predicted_translons_with_sequence_ranges.rds")
  )
  tc_paths <- c(
    file.path(refFolder(df), "predicted_translons", "TransCode",
              "predicted_translons_ranges.qs"),
    file.path(dirname(df@fafile), "predicted_translons", "TransCode",
              "predicted_translons_ranges.qs")
  )
  m_paths <- c(
    file.path(analysis_dir, "manual_translons", "manual_translons_ranges.rds"),
    file.path("dominant_cell_states", "manual_translons", "manual_translons_ranges.rds"),
    file.path("manual_translons", "manual_translons_ranges.rds")
  )
  t <- load_translon_source(t_paths, "T")
  tc <- load_translon_source(tc_paths, "TC")
  m <- load_translon_source(m_paths, "M")
  if (length(t) == 0 && length(tc) == 0 && length(m) == 0) {
    stop("Predicted translon range files not found. Checked: ",
         paste(unique(c(
           t_paths,
           tc_paths,
           m_paths
         )), collapse = "; "))
  }
  message("Loaded translon sources: T=", length(t), ", TC=", length(tc),
          ", M=", length(m))
  c(t, tc, m)
}

collapse_to_tx_iranges <- function(mapped) {
  starts <- vapply(mapped, function(x) {
    if (length(x) == 0) return(NA_integer_)
    min(start(x))
  }, integer(1))
  ends <- vapply(mapped, function(x) {
    if (length(x) == 0) return(NA_integer_)
    max(end(x))
  }, integer(1))
  IRanges(start = starts, end = ends)
}

annotate_translon_categories <- function(translons, cds, mrna) {
  if (length(translons) == 0) return(translons)

  keep <- names(translons) %chin% names(mrna)
  if (!all(keep)) {
    warning("Dropping ", sum(!keep), " predicted translons without mRNA ranges.")
    translons <- translons[keep]
  }
  if (length(translons) == 0) return(translons)

  mapped <- suppressWarnings(pmapToTranscriptF(translons, mrna[names(translons)]))
  tx_iranges <- collapse_to_tx_iranges(mapped)
  valid <- !is.na(start(tx_iranges)) & !is.na(end(tx_iranges))
  categories <- rep(NA_character_, length(translons))

  if (any(valid)) {
    groupings <- match(names(translons)[valid], names(mrna))
    categories[valid] <- ORFik:::categorize_ORFs(
      tx_iranges[valid],
      groupings = as.integer(groupings),
      cds = cds,
      mrna = mrna,
      verbose = FALSE
    )
  }

  mcols(translons)$category <- categories
  mcols(translons)$tx_start <- start(tx_iranges)
  mcols(translons)$tx_end <- end(tx_iranges)
  translons
}

deduplicate_translons_by_tx_coordinates <- function(translons) {
  if (length(translons) == 0) return(translons)
  if (!all(c("tx_start", "tx_end") %chin% names(mcols(translons)))) {
    stop("Translons must be transcript-coordinate annotated before deduplication.")
  }
  key <- paste(names(translons), mcols(translons)$tx_start,
               mcols(translons)$tx_end, sep = "|")
  source <- as.character(mcols(translons)$translon_source)
  source[is.na(source)] <- "unknown"
  merged_source <- vapply(split(source, key), function(x) {
    paste(unique(x), collapse = "+")
  }, character(1))
  mcols(translons)$translon_sources_merged <- unname(merged_source[key])
  keep <- !duplicated(key)
  dropped <- sum(!keep)
  if (dropped > 0) {
    message("Dropped ", dropped,
            " duplicate translons with identical transcript start/stop.")
  }
  translons[keep]
}

is_model_translon_category <- function(category) {
  as.character(category) %chin% c("uORF", "uoORF", "NTE", "NTT", "internal")
}

feature_prefix_from_category <- function(category) {
  category <- as.character(category)
  if (category == "uORF") return("uORF")
  if (category == "uoORF") return("uORF")
  if (category == "NTE") return("NTE")
  if (category == "NTT") return("NTT")
  if (category == "internal") return("iORF")
  "translon"
}

valid_matrix_positions <- function(positions, max_position) {
  positions <- sort(unique(as.integer(positions)))
  positions[
    is.finite(positions) & positions >= 1L & positions <= as.integer(max_position)
  ]
}

normalize_region <- function(mat, positions, lib_sizes) {
  runs <- colnames(mat)
  positions <- valid_matrix_positions(positions, nrow(mat))
  if (length(positions) == 0) {
    return(data.table(
      Run = runs,
      raw_counts = NA_real_,
      fpkm_like = NA_real_
    ))
  }
  raw_sum <- colSums(mat[positions, , drop = FALSE], na.rm = TRUE)
  kb <- length(positions) / 1000
  lib_m <- as.numeric(lib_sizes[runs]) / 1e6
  fpkm_like <- rep(NA_real_, length(runs))
  valid_lib <- is.finite(lib_m) & lib_m > 0
  fpkm_like[valid_lib] <- (raw_sum[valid_lib] / kb) / lib_m[valid_lib]
  data.table(
    Run = runs,
    raw_counts = as.numeric(raw_sum),
    fpkm_like = as.numeric(fpkm_like)
  )
}

add_cds_reference_metrics <- function(expr, reference_expr) {
  expr <- merge(expr, reference_expr, by = "Run", all.x = TRUE, sort = FALSE)
  expr[, feature_vs_cds_mean_ratio := fifelse(
    is.finite(fpkm_like) & is.finite(cds_reference_fpkm_like) &
      cds_reference_fpkm_like > 0,
    fpkm_like / cds_reference_fpkm_like,
    NA_real_
  )]
  expr[]
}

feature_measurement_positions <- function(positions, cds_positions, category) {
  positions <- sort(unique(as.integer(positions)))
  cds_positions <- sort(unique(as.integer(cds_positions)))
  if (as.character(category) == "NTE") {
    return(list(
      positions = setdiff(positions, cds_positions),
      measurement_scope = "unique_non_cds_extension"
    ))
  }
  list(
    positions = positions,
    measurement_scope = fifelse(
      as.character(category) %chin% c("NTT", "internal"),
      "cds_shared_body_diagnostic",
      "full_translon"
    )
  )
}

classify_translon <- function(translon_gr, leader_gr, cds_gr) {
  overlaps_leader <- length(findOverlaps(translon_gr, leader_gr,
                                         ignore.strand = FALSE)) > 0
  overlaps_cds <- length(findOverlaps(translon_gr, cds_gr,
                                      ignore.strand = FALSE)) > 0
  if (overlaps_leader && overlaps_cds) return("overlapping_uorf")
  if (overlaps_leader) return("leader_uorf")
  if (overlaps_cds) return("cds_translon")
  "other_translon"
}

classify_feature <- function(translon_gr, leader_gr, cds_gr, positions,
                             cds_positions, category) {
  feature_type <- classify_translon(translon_gr, leader_gr, cds_gr)
  if (category == "uORF") feature_type <- "leader_uorf"
  if (category == "uoORF") feature_type <- "overlapping_uorf"
  if (category == "NTE") feature_type <- "nte_extension"
  if (category == "NTT") feature_type <- "ntt_truncation"
  if (category == "internal") feature_type <- "internal_orf"
  cds_overlap_bases <- length(intersect(cds_positions, positions))
  list(feature_type = feature_type, cds_overlap_bases = cds_overlap_bases)
}

measure_uorf_features <- function(df, row, translons, fst_index, lib_sizes,
                                  cds_all = NULL, leader_all = NULL) {
  gene <- row$gene_symbol
  tx_id <- row$tx_id
  cds_tx <- if (!is.null(cds_all) && tx_id %in% names(cds_all)) {
    cds_all[tx_id]
  } else {
    loadRegion(df, "cds", names.keep = tx_id)
  }
  leader_tx <- if (!is.null(leader_all) && tx_id %in% names(leader_all)) {
    leader_all[tx_id]
  } else {
    loadRegion(df, "leaders", names.keep = tx_id)
  }
  display <- make_display_region(leader_tx, cds_tx, tx_id)
  coverage <- dominant_cached_coverage_by_transcript(
    display_region = display,
    fst_index = fst_index,
    gene_symbol = gene,
    tx_id = tx_id,
    analysis_dir = analysis_dir
  )
  mat <- as.matrix(coverage)
  storage.mode(mat) <- "numeric"

  cds_positions <- valid_matrix_positions(safe_pmap_positions(cds_tx, display),
                                          nrow(mat))
  leader_gr <- unlist(leader_tx, use.names = FALSE)
  cds_gr <- unlist(cds_tx, use.names = FALSE)

  tx_translons <- translons[names(translons) == tx_id]
  if (length(tx_translons) == 0) return(NULL)
  if (!"category" %in% names(mcols(tx_translons))) {
    stop("Predicted translons must be categorized before uORF modeling.")
  }

  feature_info <- rbindlist(lapply(seq_along(tx_translons), function(i) {
    category <- as.character(mcols(tx_translons)$category[i])
    if (is.na(category)) {
      stop("Could not categorize predicted translon ", i, " for ", tx_id)
    }
    if (!is_model_translon_category(category)) return(NULL)

    translon_gr <- tx_translons[[i]]
    positions <- valid_matrix_positions(
      safe_pmap_positions(GRangesList(translon_gr), display),
      nrow(mat)
    )
    feature_class <- classify_feature(
      translon_gr, leader_gr, cds_gr, positions, cds_positions, category
    )
    measurement <- feature_measurement_positions(
      positions, cds_positions, category
    )
    feature_prefix <- feature_prefix_from_category(category)
    data.table(
      feature_id = paste0(feature_prefix, i),
      feature_prefix = feature_prefix,
      feature_type = feature_class$feature_type,
      category = category,
      translon_source = as.character(mcols(tx_translons)$translon_source[i]),
      translon_sources_merged = as.character(mcols(tx_translons)$translon_sources_merged[i]),
      tx_start = if (length(positions) > 0) min(positions) else NA_integer_,
      tx_end = if (length(positions) > 0) max(positions) else NA_integer_,
      feature_bases = length(positions),
      cds_overlap_bases = feature_class$cds_overlap_bases,
      measured_bases = length(measurement$positions),
      measurement_scope = measurement$measurement_scope,
      positions = list(positions),
      measurement_positions = list(measurement$positions)
    )
  }), fill = TRUE)
  if (is.null(feature_info) || nrow(feature_info) == 0) return(NULL)
  setorder(feature_info, tx_start, tx_end)
  feature_info[, feature_seq := seq_len(.N), by = feature_prefix]
  feature_info[, feature_id := paste0(feature_prefix, feature_seq)]
  feature_info[, c("feature_prefix", "feature_seq") := NULL]

  overlapping_positions <- unique(unlist(
    feature_info[feature_type == "overlapping_uorf", positions],
    use.names = FALSE
  ))
  clean_cds_positions <- setdiff(cds_positions, intersect(cds_positions, overlapping_positions))
  cds_reference_expr <- normalize_region(mat, cds_positions, lib_sizes)
  setnames(cds_reference_expr, c("raw_counts", "fpkm_like"),
           c("cds_reference_raw_counts", "cds_reference_fpkm_like"))

  feature_expr <- rbindlist(lapply(seq_len(nrow(feature_info)), function(i) {
    info <- feature_info[i]
    expr <- normalize_region(mat, info$measurement_positions[[1]], lib_sizes)
    expr <- add_cds_reference_metrics(expr, cds_reference_expr)
    signal_metrics <- dominant_internal_orf_signal_metrics(
      mat = mat,
      feature_positions = info$positions[[1]],
      cds_positions = cds_positions,
      feature_type = info$feature_type,
      category = info$category
    )
    expr <- merge(expr, signal_metrics, by = "Run", all.x = TRUE, sort = FALSE)
    expr[, `:=`(
      gene_symbol = gene,
      tx_id = tx_id,
      feature_id = info$feature_id,
      feature_type = info$feature_type,
      category = info$category,
      translon_source = info$translon_source,
      translon_sources_merged = info$translon_sources_merged,
      feature_bases = info$feature_bases,
      cds_overlap_bases = info$cds_overlap_bases,
      measured_bases = info$measured_bases,
      measurement_scope = info$measurement_scope
    )]
    setcolorder(expr, c("gene_symbol", "tx_id", "feature_id", "feature_type",
                        "category", "translon_source", "translon_sources_merged",
                        "feature_bases", "cds_overlap_bases", "measured_bases",
                        "measurement_scope",
                        "Run", "raw_counts", "fpkm_like",
                        "cds_reference_raw_counts", "cds_reference_fpkm_like",
                        "feature_vs_cds_mean_ratio"))
    expr
  }))

  cds_expr <- normalize_region(mat, clean_cds_positions, lib_sizes)
  cds_expr <- add_cds_reference_metrics(cds_expr, cds_reference_expr)
  cds_expr <- merge(
    cds_expr,
    dominant_internal_orf_empty_signal_metrics(cds_expr$Run),
    by = "Run",
    all.x = TRUE,
    sort = FALSE
  )
  cds_expr[, `:=`(
    gene_symbol = gene,
    tx_id = tx_id,
    feature_id = "clean_CDS",
    feature_type = "clean_CDS",
    category = "clean_CDS",
    translon_source = NA_character_,
    translon_sources_merged = NA_character_,
    feature_bases = length(clean_cds_positions),
    cds_overlap_bases = 0L,
    measured_bases = length(clean_cds_positions),
    measurement_scope = "clean_cds_after_uorf_mask"
  )]
  setcolorder(cds_expr, c("gene_symbol", "tx_id", "feature_id", "feature_type",
                          "category", "translon_source", "translon_sources_merged",
                          "feature_bases", "cds_overlap_bases", "measured_bases",
                          "measurement_scope",
                          "Run", "raw_counts", "fpkm_like",
                          "cds_reference_raw_counts", "cds_reference_fpkm_like",
                          "feature_vs_cds_mean_ratio"))
  feature_expr <- rbind(
    feature_expr,
    cds_expr
  )

  feature_expr[]
}

gene_uorf_columns <- function(wide_dt) {
  feature_cols <- setdiff(names(wide_dt), c("gene_symbol", "Run", "field",
                                            "level", "n", "clean_CDS"))
  feature_cols <- feature_cols[vapply(feature_cols, function(feature) {
    any(is.finite(wide_dt[[feature]]), na.rm = TRUE)
  }, logical(1))]
  feature_order <- function(x) {
    as.integer(gsub("^[A-Za-z_]+", "", x))
  }
    groups <- fifelse(
    grepl("^uORF[0-9]+$", feature_cols), "uORF",
    fifelse(
      grepl("^NTE[0-9]+$", feature_cols), "NTE",
      fifelse(
        grepl("^NTT[0-9]+$", feature_cols), "NTT",
        fifelse(grepl("^iORF[0-9]+$", feature_cols), "iORF", "other")
      )
    )
  )
  feature_cols[order(groups, feature_order(feature_cols), feature_cols)]
}

gene_allocation_columns <- function(wide_dt) {
  gene_uorf_columns(wide_dt)[
    grepl("^uORF[0-9]+$", gene_uorf_columns(wide_dt))
  ]
}

fit_feature_model <- function(wide_dt, gene, level_scope = "sample", weight_col = NULL) {
  feature_cols <- gene_uorf_columns(wide_dt)
  if (length(feature_cols) == 0 || !"clean_CDS" %in% names(wide_dt)) {
    return(data.table())
  }

  model_dt <- copy(wide_dt)
  model_dt[, clean_CDS := log2(clean_CDS + 1)]
  for (feature in feature_cols) model_dt[, (feature) := log2(get(feature) + 1)]
  valid <- is.finite(model_dt$clean_CDS)
  for (feature in feature_cols) valid <- valid & is.finite(model_dt[[feature]])
  model_dt <- model_dt[valid]
  if (nrow(model_dt) < length(feature_cols) + 3) return(data.table())

  for (feature in c("clean_CDS", feature_cols)) {
    s <- sd(model_dt[[feature]], na.rm = TRUE)
    model_dt[, (paste0(feature, "_z")) := if (is.finite(s) && s > 0) {
      (get(feature) - mean(get(feature), na.rm = TRUE)) / s
    } else 0]
  }

  formula <- as.formula(paste(
    "clean_CDS_z ~",
    paste(paste0(feature_cols, "_z"), collapse = " + ")
  ))
  fit <- if (!is.null(weight_col) && weight_col %in% names(model_dt)) {
    lm(formula, data = model_dt, weights = model_dt[[weight_col]])
  } else {
    lm(formula, data = model_dt)
  }
  coefs <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
  coefs <- coefs[term != "(Intercept)"]
  coefs[, feature_id := sub("_z$", "", term)]
  coefs[, `:=`(
    gene_symbol = gene,
    scope = level_scope,
    n_observations = nrow(model_dt),
    cds_response = "clean_CDS",
    coefficient = Estimate,
    p_value = `Pr(>|t|)`
  )]
  coefs[, .(gene_symbol, scope, n_observations, cds_response, feature_id,
            coefficient, p_value)]
}

fit_edges <- function(wide_dt, gene, scope, weight_col = NULL) {
  feature_cols <- gene_uorf_columns(wide_dt)
  ordered_targets <- c(feature_cols[-1], "clean_CDS")
  ordered_sources <- c(feature_cols[seq_len(max(length(feature_cols) - 1, 0))],
                       tail(feature_cols, 1))
  if (length(ordered_sources) == 0) return(data.table())

  rbindlist(lapply(seq_along(ordered_sources), function(i) {
    source <- ordered_sources[[i]]
    target <- ordered_targets[[i]]
    model_dt <- copy(wide_dt)
    model_dt[, source_expr := log2(get(source) + 1)]
    model_dt[, target_expr := log2(get(target) + 1)]
    model_dt <- model_dt[is.finite(source_expr) & is.finite(target_expr)]
    if (nrow(model_dt) < 3 || length(unique(model_dt$source_expr)) < 2) {
      return(NULL)
    }
    source_sd <- sd(model_dt$source_expr)
    target_sd <- sd(model_dt$target_expr)
    if (!is.finite(source_sd) || source_sd == 0 ||
        !is.finite(target_sd) || target_sd == 0) return(NULL)
    model_dt[, source_z := (source_expr - mean(source_expr)) / source_sd]
    model_dt[, target_z := (target_expr - mean(target_expr)) / target_sd]
    fit <- if (!is.null(weight_col) && weight_col %in% names(model_dt)) {
      lm(target_z ~ source_z, data = model_dt, weights = model_dt[[weight_col]])
    } else {
      lm(target_z ~ source_z, data = model_dt)
    }
    co <- summary(fit)$coefficients["source_z", ]
    data.table(
      gene_symbol = gene,
      scope = scope,
      source = source,
      target = target,
      n_observations = nrow(model_dt),
      coefficient = unname(co[["Estimate"]]),
      p_value = unname(co[["Pr(>|t|)"]]),
      relief_effect_if_source_decreases = if (target == "clean_CDS") {
        -unname(co[["Estimate"]])
      } else NA_real_
    )
  }), fill = TRUE)
}

fit_allocation_model <- function(wide_dt, gene, scope, weight_col = NULL) {
  feature_cols <- gene_allocation_columns(wide_dt)
  if (length(feature_cols) == 0 || !"clean_CDS" %in% names(wide_dt)) {
    return(data.table())
  }

  model_dt <- copy(wide_dt)
  model_dt[, total_local_translation := rowSums(.SD, na.rm = TRUE),
           .SDcols = c(feature_cols, "clean_CDS")]
  model_dt <- model_dt[is.finite(total_local_translation) &
                         total_local_translation > 0]
  if (nrow(model_dt) < length(feature_cols) + 3) return(data.table())

  model_dt[, clean_CDS_share := clean_CDS / total_local_translation]
  model_dt[, clean_CDS_share := pmin(pmax(clean_CDS_share, 1e-6), 1 - 1e-6)]
  model_dt[, clean_CDS_share_logit := qlogis(clean_CDS_share)]
  for (feature in feature_cols) {
    model_dt[, (paste0(feature, "_share")) := get(feature) / total_local_translation]
  }
  share_cols <- paste0(feature_cols, "_share")
  valid <- is.finite(model_dt$clean_CDS_share_logit)
  for (feature in share_cols) valid <- valid & is.finite(model_dt[[feature]])
  model_dt <- model_dt[valid]
  if (nrow(model_dt) < length(feature_cols) + 3) return(data.table())

  formula <- as.formula(paste(
    "clean_CDS_share_logit ~",
    paste(share_cols, collapse = " + ")
  ))
  fit <- if (!is.null(weight_col) && weight_col %in% names(model_dt)) {
    lm(formula, data = model_dt, weights = model_dt[[weight_col]])
  } else {
    lm(formula, data = model_dt)
  }
  coefs <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
  coefs <- coefs[term != "(Intercept)"]
  coefs[, feature_id := sub("_share$", "", term)]
  coefs[, `:=`(
    gene_symbol = gene,
    scope = scope,
    n_observations = nrow(model_dt),
    response = "clean_CDS_share_logit",
    coefficient = Estimate,
    p_value = `Pr(>|t|)`,
    interpretation = fifelse(
      Estimate < 0,
      "higher feature share is associated with lower clean CDS share",
      "higher feature share is associated with higher clean CDS share"
    )
  )]
  coefs[, .(gene_symbol, scope, n_observations, response, feature_id,
            coefficient, p_value, interpretation)]
}

df <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
diagnostics <- fread(diagnostics_file)
diagnostics <- diagnostics[expression_source == "fst_clean_cds"]
if (nrow(diagnostics) == 0) {
  stop("No FST-clean marker genes found in ", diagnostics_file)
}

translons <- load_predicted_translons(df)
target_tx <- unique(diagnostics$tx_id)
cds_all <- loadRegion(df, "cds", names.keep = target_tx)
leader_all <- loadRegion(df, "leaders", names.keep = target_tx)
mrna_all <- loadRegion(df, "mrna", names.keep = target_tx)
translons <- translons[names(translons) %chin% target_tx]
translons <- annotate_translon_categories(translons, cds_all, mrna_all)
translons <- deduplicate_translons_by_tx_coordinates(translons)
message("Categorized predicted translons for selected marker transcripts: ",
        paste(names(table(mcols(translons)$category)),
              as.integer(table(mcols(translons)$category)),
              sep = "=", collapse = ", "))
message("Translon source composition after deduplication: ",
        paste(names(table(mcols(translons)$translon_sources_merged)),
              as.integer(table(mcols(translons)$translon_sources_merged)),
              sep = "=", collapse = ", "))
fst_index <- file.path(collection_dir_from_exp(df), "coverage_index.fst")
lib_sizes <- library_sizes_by_run(df)

feature_dt <- rbindlist(lapply(seq_len(nrow(diagnostics)), function(i) {
  message("Measuring uORF features for ", diagnostics$gene_symbol[[i]],
          " (", diagnostics$tx_id[[i]], ")")
  measure_uorf_features(
    df,
    diagnostics[i],
    translons,
    fst_index,
    lib_sizes,
    cds_all = cds_all,
    leader_all = leader_all
  )
}), fill = TRUE)
fwrite(feature_dt, feature_output)
feature_annotation <- unique(feature_dt[
  ,
  .(gene_symbol, feature_id, feature_type, category, feature_bases,
    cds_overlap_bases, translon_source, translon_sources_merged)
])

wide_sample <- dcast(
  feature_dt,
  gene_symbol + Run ~ feature_id,
  value.var = "fpkm_like"
)

sample_models <- rbindlist(lapply(unique(wide_sample$gene_symbol), function(gene) {
  fit_feature_model(wide_sample[gene_symbol == gene], gene, "sample")
}), fill = TRUE)
sample_models <- merge(sample_models, feature_annotation,
                       by = c("gene_symbol", "feature_id"),
                       all.x = TRUE, sort = FALSE)
sample_models[, p_adj := p.adjust(p_value, "BH"), by = gene_symbol]
fwrite(sample_models, sample_model_output)

m <- fread(metadata_file)
m <- m[, .(Run, CELL_LINE, TISSUE, CONDITION, GENE)]
feature_meta <- merge(feature_dt, m, by = "Run", all.x = TRUE, sort = FALSE)

term_feature <- rbindlist(lapply(metadata_fields, function(field) {
  tmp <- copy(feature_meta)
  tmp[, level := clean_level(get(field))]
  tmp <- tmp[!is.na(level)]
	  out <- tmp[
	    ,
	    .(
	      raw_counts = sum(raw_counts, na.rm = TRUE),
	      fpkm_like = mean(fpkm_like, na.rm = TRUE),
	      n = uniqueN(Run)
	    ),
    by = .(level, gene_symbol, feature_id)
  ][n >= min_term_n]
  out[, field := field]
  setcolorder(out, c("field", setdiff(names(out), "field")))
  out
}), fill = TRUE)

wide_term <- dcast(
  term_feature,
  field + level + n + gene_symbol ~ feature_id,
  value.var = "fpkm_like"
)

term_models <- rbindlist(lapply(unique(wide_term$gene_symbol), function(gene) {
  fit_feature_model(wide_term[gene_symbol == gene], gene, "metadata_term",
                    weight_col = "n")
}), fill = TRUE)
term_models <- merge(term_models, feature_annotation,
                     by = c("gene_symbol", "feature_id"),
                     all.x = TRUE, sort = FALSE)
term_models[, p_adj := p.adjust(p_value, "BH"), by = gene_symbol]
fwrite(term_models, term_model_output)

edges <- rbind(
  rbindlist(lapply(unique(wide_sample$gene_symbol), function(gene) {
    fit_edges(wide_sample[gene_symbol == gene], gene, "sample")
  }), fill = TRUE),
  rbindlist(lapply(unique(wide_term$gene_symbol), function(gene) {
    fit_edges(wide_term[gene_symbol == gene], gene, "metadata_term",
              weight_col = "n")
  }), fill = TRUE),
  fill = TRUE
)
edges <- merge(edges, feature_annotation[, .(
  gene_symbol,
  source = feature_id,
  source_type = feature_type,
  source_bases = feature_bases
)], by = c("gene_symbol", "source"), all.x = TRUE, sort = FALSE)
edges <- merge(edges, feature_annotation[, .(
  gene_symbol,
  target = feature_id,
  target_type = feature_type,
  target_bases = feature_bases
)], by = c("gene_symbol", "target"), all.x = TRUE, sort = FALSE)
edges[, p_adj := p.adjust(p_value, "BH"), by = .(gene_symbol, scope)]
fwrite(edges, edge_output)

allocation_models <- rbind(
  rbindlist(lapply(unique(wide_sample$gene_symbol), function(gene) {
    fit_allocation_model(wide_sample[gene_symbol == gene], gene, "sample")
  }), fill = TRUE),
  rbindlist(lapply(unique(wide_term$gene_symbol), function(gene) {
    fit_allocation_model(wide_term[gene_symbol == gene], gene,
                         "metadata_term", weight_col = "n")
  }), fill = TRUE),
  fill = TRUE
)
allocation_models <- merge(allocation_models, feature_annotation,
                           by = c("gene_symbol", "feature_id"),
                           all.x = TRUE, sort = FALSE)
allocation_models[, p_adj := p.adjust(p_value, "BH"),
                  by = .(gene_symbol, scope)]
fwrite(allocation_models, allocation_model_output)

message("Saved: ", feature_output)
message("Saved: ", sample_model_output)
message("Saved: ", term_model_output)
message("Saved: ", edge_output)
message("Saved: ", allocation_model_output)

print(sample_models)
print(term_models)
print(edges)
print(allocation_models)

invisible(list(
  feature_expression = feature_dt,
  sample_models = sample_models,
  term_models = term_models,
  edges = edges,
  allocation_models = allocation_models
))
