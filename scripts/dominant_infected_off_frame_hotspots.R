#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(GenomicRanges)
  library(IRanges)
})

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (basename(here) == "dominant_cell_states" &&
        dir.exists(file.path(here, "scripts"))) {
      return(here)
    }
    candidate <- file.path(here, "dominant_cell_states")
    if (dir.exists(candidate) &&
        dir.exists(file.path(candidate, "scripts"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find dominant_cell_states analysis directory from: ",
       start, call. = FALSE)
}

analysis_dir <- find_analysis_dir()
output_dir <- file.path(analysis_dir, "dominant_infected_off_frame_hotspots")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)
coverage_cache_file <- file.path(analysis_dir, "scripts", "dominant_coverage_cache.R")
if (!file.exists(coverage_cache_file)) {
  stop("Missing dominant coverage cache helper: ", coverage_cache_file,
       call. = FALSE)
}
source(coverage_cache_file)

message("Infected off-frame CDS hotspot localization")
message("  1. Use off-frame candidate gene/run table from previous pass")
message("  2. Load cached transcript coverage for prioritized genes")
message("  3. Localize top-1% CDS positions and keep off-frame hotspots")
message("  4. Compare recurrent infected hotspots to matched controls")

candidate_dir <- file.path(analysis_dir, "dominant_infected_off_frame_cds_scan")
gene_scores_file <- file.path(
  candidate_dir, "infected_off_frame_cds_candidate_gene_scores.csv"
)
top_runs_file <- file.path(
  candidate_dir, "infected_off_frame_cds_top_candidate_runs.csv"
)
run_metrics_file <- file.path(
  analysis_dir, "dominant_cds_frame_bumpiness_qc",
  "cds_frame_bumpiness_run_metrics.csv"
)
housekeeping_file <- file.path(
  analysis_dir, "dominant_cds_frame_bumpiness_qc",
  "housekeeping_baseline_run_qc.csv"
)
diagnostics_file <- file.path(
  analysis_dir, "human_dominant_cell_states_clean_cds_gene_diagnostics.csv"
)

for (file in c(gene_scores_file, top_runs_file, run_metrics_file,
               housekeeping_file, diagnostics_file)) {
  if (!file.exists(file)) stop("Missing required file: ", file, call. = FALSE)
}

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

collapse_text <- function(x, max_n = 6L) {
  x <- unique(safe_chr(x))
  x <- x[nzchar(x)]
  if (!length(x)) return(NA_character_)
  suffix <- if (length(x) > max_n) ";..." else ""
  paste0(paste(head(x, max_n), collapse = ";"), suffix)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

safe_fisher_p <- function(case_hits, case_total, control_hits, control_total) {
  if (control_total < 2L || case_total < 2L) return(NA_real_)
  case_miss <- max(case_total - case_hits, 0L)
  control_miss <- max(control_total - control_hits, 0L)
  suppressWarnings(stats::fisher.test(matrix(
    c(case_hits, case_miss, control_hits, control_miss),
    nrow = 2,
    byrow = TRUE
  ), alternative = "greater")$p.value)
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

middle_cds_positions <- function(cds_positions, min_edge = 45L,
                                 fraction_edge = 0.10) {
  cds_positions <- sort(unique(as.integer(cds_positions)))
  n <- length(cds_positions)
  if (n < 120L) return(cds_positions)
  edge <- max(as.integer(min_edge), floor(n * fraction_edge))
  if (n <= 2L * edge + 30L) {
    edge <- max(15L, floor(n * 0.05))
  }
  cds_positions[(edge + 1L):(n - edge)]
}

regex_has <- function(x, pattern) {
  grepl(pattern, x, ignore.case = TRUE, perl = TRUE)
}

classify_pathogen <- function(text) {
  fifelse(
    regex_has(text, "toxoplasma|tachyzoite|\\brh[_ -]?rep|\\brh[_ -]?tox"),
    "Toxoplasma",
    fifelse(
      regex_has(text, "oc43|hcov[-_ ]?oc43"),
      "HCoV-OC43",
      fifelse(
        regex_has(text, "sars[-_ ]?cov[-_ ]?2|covid"),
        "SARS-CoV-2",
        fifelse(
          regex_has(text, "influenza|\\biav\\b|\\bpr8\\b"),
          "Influenza/IAV",
          fifelse(
            regex_has(text, "\\bebv\\b|epstein"),
            "EBV",
            fifelse(
              regex_has(text, "herpes|hsv|cmv"),
              "Herpesvirus",
              fifelse(
                regex_has(text, "hcv|hepatitis c"),
                "HCV",
                fifelse(
                  regex_has(text, "hiv"),
                  "HIV",
                  fifelse(
                    regex_has(text, "virus|viral|infect|vaccinia|zika|dengue|sendai|chikungunya|adenovirus|norovirus"),
                    "infection/viral",
                    "none_or_unspecified"
                  )
                )
              )
            )
          )
        )
      )
    )
  )
}

classify_growth <- function(text) {
  fifelse(
    regex_has(text, "sub[-_ ]?confluent"),
    "subconfluent",
    fifelse(
      regex_has(text, "confluent"),
      "confluent",
      fifelse(
        regex_has(text, "tumou?r|cancer|carcinoma"),
        "tumor/cancer",
        "unspecified"
      )
    )
  )
}

classify_infection_status <- function(text) {
  control_regex <- paste(
    c("uninfected", "mock", "sham", "vehicle", "untreated",
      "\\bctrl\\b", "\\bcontrol\\b", "\\bwt\\b", "wild[-_ ]?type"),
    collapse = "|"
  )
  infection_regex <- paste(
    c("infect", "infection", "virus", "viral", "hcov", "sars", "covid",
      "influenza", "\\biav\\b", "\\bpr8\\b", "toxoplasma", "tachyzoite",
      "vaccinia", "zika", "dengue", "sendai", "chikungunya", "adenovirus",
      "norovirus", "\\bebv\\b", "epstein", "herpes", "hsv", "cmv",
      "hcv", "hiv"),
    collapse = "|"
  )
  is_control <- regex_has(text, control_regex)
  is_infected <- !is_control & regex_has(text, infection_regex)
  fifelse(is_control, "control",
          fifelse(is_infected, "infected", "other"))
}

make_clean_context_value <- function(x, fallback = "NA") {
  x <- safe_chr(x)
  fallback <- rep_len(safe_chr(fallback), length(x))
  empty <- !nzchar(x)
  x[empty] <- fallback[empty]
  x[!nzchar(x)] <- "NA"
  x
}

first_existing_cols <- function(file, wanted) {
  available <- names(fread(file, nrows = 0, showProgress = FALSE))
  intersect(wanted, available)
}

add_metadata_frame_qc <- function(dt) {
  usage_cols <- paste0("Frame_usage_", 0:2)
  if ("metadata_possible_misshift" %in% names(dt)) {
    dt[is.na(metadata_possible_misshift), metadata_possible_misshift := FALSE]
    return(dt)
  }
  if (!all(usage_cols %in% names(dt))) {
    dt[, metadata_frame0_fraction := NA_real_]
    dt[, metadata_frame1_fraction := NA_real_]
    dt[, metadata_frame2_fraction := NA_real_]
    dt[, metadata_canonical_frame_fraction := NA_real_]
    dt[, metadata_dominant_frame := NA_character_]
    dt[, metadata_possible_misshift := FALSE]
    return(dt)
  }
  raw0 <- suppressWarnings(as.numeric(dt[["Frame_usage_0"]]))
  raw1 <- suppressWarnings(as.numeric(dt[["Frame_usage_1"]]))
  raw2 <- suppressWarnings(as.numeric(dt[["Frame_usage_2"]]))
  raw_total <- raw0 + raw1 + raw2
  usage_is_percent <- is.finite(raw_total) & raw_total > 1.5
  frame0 <- ifelse(usage_is_percent, raw0 / 100, raw0)
  frame1 <- ifelse(usage_is_percent, raw1 / 100, raw1)
  frame2 <- ifelse(usage_is_percent, raw2 / 100, raw2)
  available <- is.finite(frame0) & is.finite(frame1) & is.finite(frame2) &
    is.finite(raw_total) & raw_total > 0
  dominant_frame <- paste0(
    "frame",
    max.col(cbind(frame0, frame1, frame2), ties.method = "first") - 1L
  )
  dominant_frame[!available] <- NA_character_
  dt[, metadata_frame0_fraction := frame0]
  dt[, metadata_frame1_fraction := frame1]
  dt[, metadata_frame2_fraction := frame2]
  dt[!available, c("metadata_frame0_fraction", "metadata_frame1_fraction",
                   "metadata_frame2_fraction") := .(NA_real_, NA_real_,
                                                    NA_real_)]
  dt[, metadata_canonical_frame_fraction := metadata_frame0_fraction]
  dt[, metadata_dominant_frame := dominant_frame]
  dt[, metadata_possible_misshift :=
       is.finite(metadata_canonical_frame_fraction) &
       metadata_canonical_frame_fraction < 0.45 &
       metadata_dominant_frame != "frame0"]
  dt[is.na(metadata_possible_misshift), metadata_possible_misshift := FALSE]
  dt
}

add_context_fields <- function(dt) {
  for (nm in c("sample_title", "sample_source", "CONDITION", "INHIBITOR",
               "TIMEPOINT", "TISSUE", "CELL_LINE", "GENE", "FRACTION",
               "Cancer_type", "Cell_model", "Cell_type", "Organ_system",
               "Sex", "Life_stage", "study", "BioProject", "AUTHOR",
               "LIBRARYTYPE")) {
    if (!nm %in% names(dt)) dt[, (nm) := NA_character_]
  }
  dt[, context_text := tolower(paste(
    safe_chr(sample_title), safe_chr(sample_source), safe_chr(CONDITION),
    safe_chr(INHIBITOR), safe_chr(TIMEPOINT), safe_chr(TISSUE),
    safe_chr(CELL_LINE), safe_chr(GENE), safe_chr(FRACTION),
    safe_chr(Cancer_type), safe_chr(Cell_model), safe_chr(Cell_type),
    safe_chr(Organ_system), safe_chr(study), safe_chr(BioProject),
    safe_chr(AUTHOR), safe_chr(LIBRARYTYPE),
    sep = " | "
  ))]
  dt[, infection_status := classify_infection_status(context_text)]
  dt[, pathogen := classify_pathogen(context_text)]
  dt[, growth_context := classify_growth(context_text)]
  dt[, inhibitor_context :=
       make_clean_context_value(INHIBITOR, "none/unspecified")]
  dt[, cell_context := make_clean_context_value(CELL_LINE, "unknown_cell")]
  dt[, tissue_context := make_clean_context_value(TISSUE, "unknown_tissue")]
  dt[, bioproject_context :=
       make_clean_context_value(BioProject, "unknown_project")]
  dt[, author_context := make_clean_context_value(AUTHOR, "unknown_author")]
  dt[, timepoint_context := make_clean_context_value(TIMEPOINT, "unknown_time")]
  dt[, match_context_id := paste(
    bioproject_context, author_context, cell_context, tissue_context,
    inhibitor_context, growth_context,
    sep = " | "
  )]
  dt[]
}

per_run_top_positions <- function(mat, runs, cds_positions, use_positions,
                                  gene, tx_id, group_lookup, context_lookup,
                                  top_fraction = 0.01,
                                  max_positions_per_run = 30L) {
  runs <- intersect(runs, colnames(mat))
  if (!length(runs)) return(data.table())
  cds_positions <- sort(unique(as.integer(cds_positions)))
  use_positions <- sort(unique(as.integer(use_positions)))
  use_positions <- use_positions[
    use_positions >= 1L & use_positions <= nrow(mat)
  ]
  if (!length(use_positions)) return(data.table())
  full_rel <- match(use_positions, cds_positions)
  frame <- as.integer((use_positions - min(cds_positions)) %% 3L)
  rbindlist(lapply(runs, function(run) {
    counts <- as.numeric(mat[use_positions, run, drop = TRUE])
    total <- sum(counts, na.rm = TRUE)
    if (!is.finite(total) || total <= 0) return(data.table())
    q99 <- as.numeric(stats::quantile(counts, probs = 1 - top_fraction,
                                      na.rm = TRUE, names = FALSE))
    idx <- which(is.finite(counts) & counts > 0 & counts >= q99)
    if (!length(idx)) return(data.table())
    idx <- idx[order(counts[idx], decreasing = TRUE)]
    idx <- head(idx, max_positions_per_run)
    data.table(
      gene_symbol = gene,
      tx_id = tx_id,
      Run = run,
      hotspot_run_group = group_lookup[run],
      match_context_id = context_lookup[run],
      cds_relative_position = full_rel[idx],
      cds_percent = full_rel[idx] / length(cds_positions),
      codon_index = floor((full_rel[idx] - 1L) / 3L) + 1L,
      transcript_position = use_positions[idx],
      cds_frame = paste0("frame", frame[idx]),
      is_off_frame = frame[idx] != 0L,
      position_count = counts[idx],
      position_count_share = counts[idx] / total,
      hotspot_rank_in_run = seq_along(idx),
      run_cds_total_counts = total,
      top_threshold_count = q99
    )
  }), fill = TRUE)
}

aggregate_profile <- function(mat, runs, cds_positions, use_positions,
                              gene, tx_id, group_name) {
  runs <- intersect(runs, colnames(mat))
  if (!length(runs)) return(data.table())
  cds_positions <- sort(unique(as.integer(cds_positions)))
  use_positions <- sort(unique(as.integer(use_positions)))
  use_positions <- use_positions[
    use_positions >= 1L & use_positions <= nrow(mat)
  ]
  if (!length(use_positions)) return(data.table())
  counts <- as.numeric(rowSums(mat[use_positions, runs, drop = FALSE],
                               na.rm = TRUE))
  total <- sum(counts, na.rm = TRUE)
  rel <- match(use_positions, cds_positions)
  data.table(
    gene_symbol = gene,
    tx_id = tx_id,
    profile_group = group_name,
    n_runs = length(runs),
    cds_relative_position = rel,
    cds_percent = rel / length(cds_positions),
    codon_index = floor((rel - 1L) / 3L) + 1L,
    transcript_position = use_positions,
    cds_frame = paste0("frame", (use_positions - min(cds_positions)) %% 3L),
    counts = counts,
    normalized_density = if (is.finite(total) && total > 0) {
      counts / total * 1000
    } else {
      NA_real_
    }
  )
}

gene_scores <- fread(gene_scores_file, showProgress = FALSE)
gene_scores[, gene_symbol := toupper(gene_symbol)]
candidate_genes <- gene_scores[
  review_tier != "not_prioritized" &
    is.finite(candidate_score) &
    candidate_score > 0,
  gene_symbol
]
if (!length(candidate_genes)) stop("No prioritized genes found.", call. = FALSE)

top_runs <- fread(top_runs_file, showProgress = FALSE)
top_runs[, gene_symbol := toupper(gene_symbol)]
off_frame_candidate_runs <- top_runs[
  gene_symbol %chin% candidate_genes &
    candidate_off_frame_peak_no_hk_artifact == TRUE
]
if (!nrow(off_frame_candidate_runs)) {
  stop("No infected off-frame candidate runs found.", call. = FALSE)
}

run_cols <- c(
  "Run", "gene_symbol", "tx_id", "cds_total_counts",
  "cds_top1pct_signal", "frame0_fraction", "frame1_fraction",
  "frame2_fraction", "top1pct_frame0_fraction",
  "top1pct_frame1_fraction", "top1pct_frame2_fraction",
  "canonical_frame_fraction", "off_frame_fraction",
  "top1pct_off_frame_fraction", "global_possible_misshift",
  "Frame_usage_0", "Frame_usage_1", "Frame_usage_2",
  "metadata_frame0_fraction", "metadata_frame1_fraction",
  "metadata_frame2_fraction", "metadata_canonical_frame_fraction",
  "metadata_dominant_frame", "metadata_possible_misshift",
  "BioProject", "Study_Pubmed_id", "AUTHOR", "sample_source",
  "sample_title", "LIBRARYTYPE", "REPLICATE", "CONDITION",
  "INHIBITOR", "TIMEPOINT", "TISSUE", "CELL_LINE", "GENE",
  "FRACTION", "Cancer_type", "Cell_model", "Cell_type",
  "Organ_system", "Sex", "Life_stage", "study",
  "top1pct_bumpiness_z_vs_depth"
)
run_cols <- first_existing_cols(run_metrics_file, run_cols)
message("Reading run metrics needed for candidate contexts.")
run_dt <- fread(run_metrics_file, select = run_cols, showProgress = FALSE)
run_dt[, gene_symbol := toupper(gene_symbol)]
run_dt <- run_dt[gene_symbol %chin% candidate_genes]
run_dt <- add_metadata_frame_qc(run_dt)

hk_cols <- c(
  "Run", "housekeeping_median_top1pct_bumpiness_z_vs_depth",
  "housekeeping_canonical_frame_fraction",
  "housekeeping_possible_misshift",
  "global_possible_misshift"
)
hk_cols <- first_existing_cols(housekeeping_file, hk_cols)
hk <- fread(housekeeping_file, select = hk_cols, showProgress = FALSE)
hk <- unique(hk, by = "Run")
names(hk) <- sub("^global_", "hk_global_", names(hk))
run_dt <- merge(run_dt, hk, by = "Run", all.x = TRUE, sort = FALSE)
run_dt <- add_context_fields(run_dt)
if (!"top1pct_off_frame_fraction" %in% names(run_dt)) {
  run_dt[, top1pct_off_frame_fraction :=
           top1pct_frame1_fraction + top1pct_frame2_fraction]
}
if (!"off_frame_fraction" %in% names(run_dt)) {
  run_dt[, off_frame_fraction := frame1_fraction + frame2_fraction]
}
run_dt[, top1pct_off_frame_enrichment :=
         top1pct_off_frame_fraction - off_frame_fraction]
run_dt[, target_minus_housekeeping_bumpiness_z :=
         top1pct_bumpiness_z_vs_depth -
         housekeeping_median_top1pct_bumpiness_z_vs_depth]
run_dt[, technical_frame_bias_possible :=
         fifelse(is.na(metadata_possible_misshift), FALSE,
                 metadata_possible_misshift) |
         fifelse(is.na(global_possible_misshift), FALSE,
                 global_possible_misshift) |
         fifelse(is.na(housekeeping_possible_misshift), FALSE,
                 housekeeping_possible_misshift)]
run_dt[, evaluable_for_peak_scan :=
         cds_total_counts >= 100 &
         cds_top1pct_signal >= 20 &
         is.finite(top1pct_off_frame_fraction) &
         is.finite(top1pct_off_frame_enrichment) &
         is.finite(target_minus_housekeeping_bumpiness_z)]

off_frame_candidate_runs <- add_context_fields(copy(off_frame_candidate_runs))
candidate_run_keys <- unique(off_frame_candidate_runs[, .(gene_symbol, Run)])
run_dt <- merge(
  run_dt,
  candidate_run_keys[, candidate_infected_off_frame_run := TRUE,
                     by = .(gene_symbol, Run)],
  by = c("gene_symbol", "Run"),
  all.x = TRUE,
  sort = FALSE
)
run_dt[is.na(candidate_infected_off_frame_run),
       candidate_infected_off_frame_run := FALSE]

candidate_contexts <- unique(off_frame_candidate_runs[
  ,
  .(gene_symbol, match_context_id)
])
run_dt <- merge(
  run_dt,
  candidate_contexts[, in_candidate_context := TRUE,
                     by = .(gene_symbol, match_context_id)],
  by = c("gene_symbol", "match_context_id"),
  all.x = TRUE,
  sort = FALSE
)
run_dt[is.na(in_candidate_context), in_candidate_context := FALSE]
run_dt[, hotspot_run_group := fifelse(
  candidate_infected_off_frame_run,
  "candidate_infected",
  fifelse(
    in_candidate_context & infection_status == "control" &
      evaluable_for_peak_scan,
    "matched_control",
    fifelse(
      in_candidate_context & infection_status == "infected" &
        evaluable_for_peak_scan,
      "same_context_infected_noncandidate",
      "not_used"
    )
  )
)]
run_group_dt <- run_dt[
  hotspot_run_group != "not_used",
  .(
    gene_symbol, tx_id, Run, hotspot_run_group, match_context_id,
    BioProject, study, AUTHOR, CELL_LINE, TISSUE, CONDITION, INHIBITOR,
    TIMEPOINT, pathogen, growth_context, sample_title, sample_source,
    cds_total_counts, top1pct_off_frame_fraction,
    top1pct_off_frame_enrichment, target_minus_housekeeping_bumpiness_z,
    metadata_canonical_frame_fraction, metadata_dominant_frame,
    metadata_possible_misshift, global_possible_misshift,
    housekeeping_possible_misshift, technical_frame_bias_possible
  )
]
run_group_dt <- unique(run_group_dt, by = c("gene_symbol", "Run",
                                           "hotspot_run_group"))

df <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
run_order <- runIDs(df)
fst_index <- file.path(collection_dir_from_exp(df), "coverage_index.fst")
if (!file.exists(fst_index)) stop("Missing FST coverage index: ", fst_index)

diagnostics <- fread(diagnostics_file, showProgress = FALSE)
diagnostics[, gene_symbol := toupper(gene_symbol)]
diagnostics <- diagnostics[
  expression_source == "fst_clean_cds" &
    gene_symbol %chin% candidate_genes &
    !is.na(tx_id) & nzchar(tx_id)
]
selected_tx <- unique(diagnostics[, .(gene_symbol, tx_id)])
selected_tx <- selected_tx[match(candidate_genes, gene_symbol), ]
selected_tx <- selected_tx[!is.na(gene_symbol)]

cds_all <- loadRegion(df, part = "cds", names.keep = selected_tx$tx_id)
leader_all <- loadRegion(df, part = "leaders", names.keep = selected_tx$tx_id)

all_hotspots <- list()
all_profiles <- list()
gene_run_summary <- list()

for (i in seq_len(nrow(selected_tx))) {
  gene <- selected_tx$gene_symbol[[i]]
  tx_id <- selected_tx$tx_id[[i]]
  gene_runs <- run_group_dt[gene_symbol == gene]
  if (!nrow(gene_runs)) next
  if (!tx_id %in% names(cds_all)) next
  message("Localizing hotspots for ", gene, " (", tx_id, ")")
  cds_tx <- cds_all[tx_id]
  leader_tx <- if (tx_id %in% names(leader_all)) {
    leader_all[tx_id]
  } else {
    GRangesList()
  }
  display_region <- make_display_region(leader_tx, cds_tx, tx_id)
  cds_positions <- safe_pmap_positions(cds_tx, display_region)
  use_positions <- middle_cds_positions(cds_positions)
  if (!length(use_positions)) next

  coverage <- dominant_cached_coverage_by_transcript(
    display_region = display_region,
    fst_index = fst_index,
    gene_symbol = gene,
    tx_id = tx_id,
    analysis_dir = analysis_dir
  )
  mat <- as.matrix(coverage)
  storage.mode(mat) <- "numeric"
  if (!identical(colnames(mat), run_order)) {
    mat <- mat[, run_order, drop = FALSE]
  }

  group_lookup <- setNames(gene_runs$hotspot_run_group, gene_runs$Run)
  context_lookup <- setNames(gene_runs$match_context_id, gene_runs$Run)
  run_set <- unique(gene_runs$Run)
  all_hotspots[[gene]] <- per_run_top_positions(
    mat = mat,
    runs = run_set,
    cds_positions = cds_positions,
    use_positions = use_positions,
    gene = gene,
    tx_id = tx_id,
    group_lookup = group_lookup,
    context_lookup = context_lookup
  )

  all_profiles[[gene]] <- rbindlist(lapply(
    c("candidate_infected", "matched_control",
      "same_context_infected_noncandidate"),
    function(group_name) {
      aggregate_profile(
        mat = mat,
        runs = gene_runs[hotspot_run_group == group_name, unique(Run)],
        cds_positions = cds_positions,
        use_positions = use_positions,
        gene = gene,
        tx_id = tx_id,
        group_name = group_name
      )
    }
  ), fill = TRUE)

  gene_run_summary[[gene]] <- gene_runs[
    ,
    .(
      n_runs = uniqueN(Run),
      median_cds_total_counts = as.numeric(safe_median(cds_total_counts)),
      median_top1pct_off_frame_fraction =
        as.numeric(safe_median(top1pct_off_frame_fraction)),
      median_target_minus_hk_bumpiness_z =
        as.numeric(safe_median(target_minus_housekeeping_bumpiness_z)),
      contexts = collapse_text(match_context_id, max_n = 4L),
      pathogens = collapse_text(pathogen, max_n = 4L),
      top_runs = collapse_text(
        Run[order(-top1pct_off_frame_fraction,
                  -target_minus_housekeeping_bumpiness_z)],
        max_n = 8L
      )
    ),
    by = .(gene_symbol, tx_id, hotspot_run_group)
  ]
}

hotspots <- rbindlist(all_hotspots, fill = TRUE)
profiles <- rbindlist(all_profiles, fill = TRUE)
run_summary <- rbindlist(gene_run_summary, fill = TRUE)

fwrite(hotspots, file.path(output_dir,
                           "infected_off_frame_hotspot_run_positions.csv"))
fwrite(profiles, file.path(output_dir,
                           "infected_off_frame_hotspot_aggregate_profiles.csv"))
fwrite(run_summary, file.path(output_dir,
                              "infected_off_frame_hotspot_group_run_summary.csv"))

off_hotspots <- hotspots[is_off_frame == TRUE]
group_sizes <- run_group_dt[
  hotspot_run_group != "not_used",
  .(group_n_runs = uniqueN(Run)),
  by = .(gene_symbol, hotspot_run_group)
]
position_group <- off_hotspots[
  ,
  .(
    n_runs_with_hotspot = uniqueN(Run),
    total_position_counts = sum(position_count, na.rm = TRUE),
    median_position_count_share = safe_median(position_count_share),
    max_position_count_share = safe_max(position_count_share),
    best_rank_in_run = min(hotspot_rank_in_run, na.rm = TRUE),
    run_examples = collapse_text(Run[order(hotspot_rank_in_run)], max_n = 8L),
    contexts = collapse_text(match_context_id, max_n = 4L)
  ),
  by = .(
    gene_symbol, tx_id, cds_relative_position, cds_percent, codon_index,
    transcript_position, cds_frame, hotspot_run_group
  )
]
position_group <- merge(
  position_group,
  group_sizes,
  by = c("gene_symbol", "hotspot_run_group"),
  all.x = TRUE,
  sort = FALSE
)
position_group[, fraction_runs_with_hotspot :=
                 n_runs_with_hotspot / pmax(group_n_runs, 1)]

hotspot_summary <- dcast(
  position_group,
  gene_symbol + tx_id + cds_relative_position + cds_percent + codon_index +
    transcript_position + cds_frame ~ hotspot_run_group,
  value.var = c(
    "n_runs_with_hotspot", "group_n_runs", "fraction_runs_with_hotspot",
    "total_position_counts", "median_position_count_share",
    "max_position_count_share", "best_rank_in_run", "run_examples",
    "contexts"
  ),
  fill = 0
)

ensure_col <- function(dt, nm, value) {
  if (!nm %in% names(dt)) dt[, (nm) := value]
  invisible(dt)
}
for (nm in c("n_runs_with_hotspot_candidate_infected",
             "n_runs_with_hotspot_matched_control",
             "n_runs_with_hotspot_same_context_infected_noncandidate",
             "group_n_runs_candidate_infected",
             "group_n_runs_matched_control",
             "group_n_runs_same_context_infected_noncandidate",
             "fraction_runs_with_hotspot_candidate_infected",
             "fraction_runs_with_hotspot_matched_control",
             "fraction_runs_with_hotspot_same_context_infected_noncandidate",
             "total_position_counts_candidate_infected",
             "median_position_count_share_candidate_infected",
             "max_position_count_share_candidate_infected",
             "best_rank_in_run_candidate_infected")) {
  ensure_col(hotspot_summary, nm, 0)
}

hotspot_summary <- hotspot_summary[
  n_runs_with_hotspot_candidate_infected > 0
]
hotspot_summary[, case_minus_control_fraction :=
                  fraction_runs_with_hotspot_candidate_infected -
                  fraction_runs_with_hotspot_matched_control]
hotspot_summary[, case_minus_same_context_infected_fraction :=
                  fraction_runs_with_hotspot_candidate_infected -
                  fraction_runs_with_hotspot_same_context_infected_noncandidate]
hotspot_summary[, fisher_p_matched_control := mapply(
  safe_fisher_p,
  n_runs_with_hotspot_candidate_infected,
  group_n_runs_candidate_infected,
  n_runs_with_hotspot_matched_control,
  group_n_runs_matched_control
)]
hotspot_summary[, q_matched_control :=
                  p.adjust(fisher_p_matched_control, method = "BH")]
hotspot_summary[, hotspot_score :=
                  2.0 * fraction_runs_with_hotspot_candidate_infected +
                  log1p(n_runs_with_hotspot_candidate_infected) +
                  pmax(0, case_minus_control_fraction) +
                  0.5 * (group_n_runs_matched_control >= 2 &
                           n_runs_with_hotspot_matched_control == 0) +
                  pmax(0, median_position_count_share_candidate_infected) -
                  1.5 * pmax(0, fraction_runs_with_hotspot_matched_control -
                               fraction_runs_with_hotspot_candidate_infected)]
hotspot_summary[, hotspot_review_class := fifelse(
  group_n_runs_matched_control >= 2 &
    n_runs_with_hotspot_candidate_infected >= 2 &
    fraction_runs_with_hotspot_candidate_infected >= 0.25 &
    case_minus_control_fraction >= 0.50 &
    fraction_runs_with_hotspot_matched_control <= 0.25,
  "matched_specific_recurrent_infected_hotspot",
  fifelse(
    group_n_runs_matched_control >= 2 &
      n_runs_with_hotspot_candidate_infected >= 2 &
      fraction_runs_with_hotspot_candidate_infected >
        fraction_runs_with_hotspot_matched_control,
    "matched_enriched_shared_hotspot",
    fifelse(
    group_n_runs_matched_control < 2 &
      n_runs_with_hotspot_candidate_infected >= 2 &
      fraction_runs_with_hotspot_candidate_infected >= 0.15,
    "unmatched_recurrent_infected_hotspot",
    fifelse(
      n_runs_with_hotspot_candidate_infected == 1,
      "single_run_infected_hotspot",
      fifelse(
        fraction_runs_with_hotspot_matched_control >=
          fraction_runs_with_hotspot_candidate_infected,
        "shared_or_control_hotspot",
        "weak_recurrent_infected_hotspot"
      )
    )
    )
  )
)]
setorder(hotspot_summary, -hotspot_score,
         -n_runs_with_hotspot_candidate_infected,
         gene_symbol, cds_relative_position)

gene_hotspot_summary <- hotspot_summary[
  ,
  .(
    n_candidate_hotspots = .N,
    n_matched_recurrent_hotspots =
      sum(hotspot_review_class ==
            "matched_specific_recurrent_infected_hotspot"),
    n_matched_enriched_shared_hotspots =
      sum(hotspot_review_class == "matched_enriched_shared_hotspot"),
    n_unmatched_recurrent_hotspots =
      sum(hotspot_review_class == "unmatched_recurrent_infected_hotspot"),
    n_weak_recurrent_hotspots =
      sum(hotspot_review_class == "weak_recurrent_infected_hotspot"),
    n_single_run_hotspots =
      sum(hotspot_review_class == "single_run_infected_hotspot"),
    best_hotspot_score = max(hotspot_score, na.rm = TRUE),
    best_hotspot_class = hotspot_review_class[which.max(hotspot_score)],
    best_codon_index = codon_index[which.max(hotspot_score)],
    best_cds_frame = cds_frame[which.max(hotspot_score)],
    best_fraction_runs_with_hotspot =
      fraction_runs_with_hotspot_candidate_infected[which.max(hotspot_score)],
    best_case_minus_control_fraction =
      case_minus_control_fraction[which.max(hotspot_score)],
    best_run_examples =
      run_examples_candidate_infected[which.max(hotspot_score)]
  ),
  by = .(gene_symbol, tx_id)
]
gene_hotspot_summary <- merge(
  gene_hotspot_summary,
  gene_scores[
    ,
    .(gene_symbol, candidate_score, review_tier,
      total_case_n_off_frame_peak, top_context_label)
  ],
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)
setorder(gene_hotspot_summary, -n_matched_recurrent_hotspots,
         -n_matched_enriched_shared_hotspots,
         -n_unmatched_recurrent_hotspots, -best_hotspot_score)

fwrite(hotspot_summary, file.path(output_dir,
                                  "infected_off_frame_hotspot_recurrence.csv"))
fwrite(gene_hotspot_summary, file.path(output_dir,
                                       "infected_off_frame_hotspot_gene_summary.csv"))

summary_metrics <- data.table(
  metric = c(
    "n_prioritized_genes",
    "n_genes_with_hotspots",
    "n_candidate_infected_runs",
    "n_matched_control_runs",
    "n_same_context_noncandidate_infected_runs",
    "n_off_frame_hotspot_positions",
    "n_matched_recurrent_hotspots",
    "n_matched_enriched_shared_hotspots",
    "n_unmatched_recurrent_hotspots",
    "n_weak_recurrent_hotspots"
  ),
  value = as.character(c(
    length(candidate_genes),
    uniqueN(hotspot_summary$gene_symbol),
    run_group_dt[hotspot_run_group == "candidate_infected", uniqueN(Run)],
    run_group_dt[hotspot_run_group == "matched_control", uniqueN(Run)],
    run_group_dt[
      hotspot_run_group == "same_context_infected_noncandidate",
      uniqueN(Run)
    ],
    nrow(hotspot_summary),
    hotspot_summary[
      hotspot_review_class ==
        "matched_specific_recurrent_infected_hotspot", .N
    ],
    hotspot_summary[
      hotspot_review_class == "matched_enriched_shared_hotspot", .N
    ],
    hotspot_summary[
      hotspot_review_class == "unmatched_recurrent_infected_hotspot", .N
    ],
    hotspot_summary[
      hotspot_review_class == "weak_recurrent_infected_hotspot", .N
    ]
  ))
)
fwrite(summary_metrics, file.path(output_dir,
                                  "infected_off_frame_hotspot_summary_metrics.csv"))

top_hotspots <- head(hotspot_summary, 80L)
if (nrow(top_hotspots) > 0) {
  gene_order <- gene_hotspot_summary$gene_symbol
  top_hotspots[, gene_label := factor(gene_symbol, levels = rev(gene_order))]
  p <- ggplot(
    top_hotspots,
    aes(
      x = cds_percent,
      y = gene_label,
      size = n_runs_with_hotspot_candidate_infected,
      fill = case_minus_control_fraction,
      text = paste0(
        gene_symbol,
        "<br>codon=", codon_index,
        "<br>frame=", cds_frame,
        "<br>class=", hotspot_review_class,
        "<br>case hits=", n_runs_with_hotspot_candidate_infected,
        " / ", group_n_runs_candidate_infected,
        "<br>control hits=", n_runs_with_hotspot_matched_control,
        " / ", group_n_runs_matched_control,
        "<br>case-control fraction delta=",
        round(case_minus_control_fraction, 3),
        "<br>score=", round(hotspot_score, 3),
        "<br>runs=", run_examples_candidate_infected
      )
    )
  ) +
    geom_point(shape = 21, color = "#303030", alpha = 0.92) +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#8b0000",
      midpoint = 0, na.value = "grey90"
    ) +
    scale_size_continuous(range = c(2.5, 9)) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "Recurrent Infected Off-Frame CDS Hotspots",
      subtitle = "Dots are off-frame top-1% CDS positions in infected candidate runs; x is relative CDS position.",
      x = "Relative CDS position",
      y = NULL,
      fill = "Case-control\nfraction delta",
      size = "Candidate\nrun hits"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  ggsave(file.path(figure_dir,
                   "infected_off_frame_hotspot_recurrence.png"),
         p, width = 11, height = 7.5, dpi = 150)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p,
      file.path(figure_dir,
                "infected_off_frame_hotspot_recurrence.html"),
      title = "infected_off_frame_hotspot_recurrence"
    )
  }
}

top_profile_genes <- head(gene_hotspot_summary$gene_symbol, 8L)
profile_plot_dt <- profiles[
  gene_symbol %chin% top_profile_genes &
    profile_group %chin% c("candidate_infected", "matched_control")
]
if (nrow(profile_plot_dt) > 0) {
  profile_plot_dt[, gene_symbol := factor(
    gene_symbol, levels = top_profile_genes
  )]
  p2 <- ggplot(
    profile_plot_dt,
    aes(
      x = cds_percent,
      y = normalized_density,
      color = profile_group,
      group = profile_group,
      text = paste0(
        as.character(gene_symbol),
        "<br>group=", profile_group,
        "<br>codon=", codon_index,
        "<br>frame=", cds_frame,
        "<br>density=", round(normalized_density, 3),
        "<br>n runs=", n_runs
      )
    )
  ) +
    geom_line(linewidth = 0.35, alpha = 0.85, na.rm = TRUE) +
    facet_wrap(~ gene_symbol, scales = "free_y", ncol = 2) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_color_manual(
      values = c(candidate_infected = "#8b0000",
                 matched_control = "#2166ac"),
      drop = FALSE
    ) +
    labs(
      title = "Aggregate Middle-CDS Profiles for Off-Frame Hotspot Candidates",
      subtitle = "Candidate infected runs are compared to matched controls when available.",
      x = "Relative CDS position",
      y = "Normalized density per 1000 CDS reads",
      color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  ggsave(file.path(figure_dir,
                   "infected_off_frame_hotspot_aggregate_profiles.png"),
         p2, width = 11, height = 8, dpi = 150)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p2,
      file.path(figure_dir,
                "infected_off_frame_hotspot_aggregate_profiles.html"),
      title = "infected_off_frame_hotspot_aggregate_profiles"
    )
  }
}

message("Saved infected off-frame hotspot outputs to: ", output_dir)
message("Top hotspot gene summary:")
print(head(gene_hotspot_summary[
  ,
    .(gene_symbol, review_tier, candidate_score,
    n_matched_recurrent_hotspots, n_matched_enriched_shared_hotspots,
    n_unmatched_recurrent_hotspots,
    best_hotspot_class, best_hotspot_score, best_codon_index,
    best_cds_frame, best_fraction_runs_with_hotspot)
], 25L))
