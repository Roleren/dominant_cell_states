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

analysis_dir <- if (dir.exists("dominant_cell_states")) {
  "dominant_cell_states"
} else {
  "."
}

manual_dir <- file.path(analysis_dir, "manual_translons")
dir.create(manual_dir, recursive = TRUE, showWarnings = FALSE)

manual_ranges_file <- file.path(manual_dir, "manual_translons_ranges.rds")
candidate_manifest_file <- file.path(manual_dir, "manual_translon_candidates.csv")
selected_manifest_file <- file.path(manual_dir, "manual_translons_manifest.csv")
manual_additions_file <- file.path(manual_dir, "translons_to_be_added_to_manual.txt")

discovered_manual_targets <- data.table(
  gene_symbol = c("SNAI2", "MMP2"),
  tx_id = c("ENST00000020945", "ENST00000219070")
)

fixed_manual_targets <- data.table(
  gene_symbol = "LMNB1",
  tx_id = "ENST00000261366",
  tx_start = 313L,
  tx_end = 795L,
  start_codon = "CTG",
  stop_codon = "TGA",
  selection_reason = paste(
    "user-specified most-upstream active LMNB1 CTG overlapping uORF;",
    "alternative starts exist, but tx 313-795 is used as manual M"
  ),
  manual_type = "uoORF",
  strict_category_match = FALSE
)

normalize_manual_type <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "", x)
  fifelse(
    x %chin% c("ouorf", "uoorf", "uo", "overlappinguorf", "overlapuorf"),
    "uoORF",
    fifelse(
      x %chin% c("uorf", "leaderuorf"),
      "uORF",
      fifelse(
        x %chin% c("iorf", "internal", "internalorf", "cdsorf"),
        "internal",
        x
      )
    )
  )
}

expected_category_from_manual_type <- function(x) {
  x <- normalize_manual_type(x)
  fifelse(
    x == "uoORF",
    "uoORF",
    fifelse(
      x == "uORF",
      "uORF",
      fifelse(x == "internal", "internal", NA_character_)
    )
  )
}

manual_category_override_allowed <- function(manual_type, category) {
  manual_type <- normalize_manual_type(manual_type)
  category <- as.character(category)
  manual_type == "internal" & category %chin% c("internal", "uORF", "uoORF")
}

read_manual_additions <- function(path, analysis_dir) {
  if (!file.exists(path)) {
    return(data.table())
  }
  additions <- fread(path, strip.white = TRUE)
  setnames(additions, tolower(gsub("[^A-Za-z0-9]+", "_", names(additions))))
  required <- c("gene", "start", "stop", "type")
  if (!all(required %chin% names(additions))) {
    stop("Manual additions file must contain columns: gene, start, stop, type")
  }
  additions <- additions[!is.na(gene) & !is.na(start) & !is.na(stop)]
  if (nrow(additions) == 0) return(data.table())

  additions[, gene_symbol := toupper(trimws(gene))]
  additions[, tx_start := as.integer(start)]
  additions[, tx_end := as.integer(stop)]
  additions[, manual_type := normalize_manual_type(type)]
  additions[, `:=`(gene = NULL, start = NULL, stop = NULL, type = NULL)]
  additions <- additions[!is.na(tx_start) & !is.na(tx_end)]
  reverse_rows <- which(additions$tx_start > additions$tx_end)
  if (length(reverse_rows) > 0) {
    start_original <- additions$tx_start[reverse_rows]
    additions$tx_start[reverse_rows] <- additions$tx_end[reverse_rows]
    additions$tx_end[reverse_rows] <- start_original
  }

  tx_lookup <- data.table()
  diagnostics_file <- file.path(
    analysis_dir,
    "human_dominant_cell_states_clean_cds_gene_diagnostics.csv"
  )
  rdg_annotation_file <- file.path(
    analysis_dir,
    "dominant_rdg_outputs",
    "rdg_feature_annotation.csv"
  )
  if (file.exists(diagnostics_file)) {
    tx_lookup <- fread(diagnostics_file)[
      ,
      .(gene_symbol = toupper(gene_symbol), tx_id, lookup_source = "clean_cds_diagnostics")
    ]
  }
  if (nrow(tx_lookup) == 0 && file.exists(rdg_annotation_file)) {
    tx_lookup <- unique(fread(rdg_annotation_file)[
      ,
      .(gene_symbol = toupper(gene_symbol), tx_id, lookup_source = "rdg_feature_annotation")
    ])
  }
  if (nrow(tx_lookup) == 0) {
    stop(
      "Could not resolve manual addition gene symbols to transcript IDs. ",
      "Rerun clean-CDS diagnostics first or add tx_id support to the manual file."
    )
  }

  additions <- merge(additions, tx_lookup, by = "gene_symbol", all.x = TRUE, sort = FALSE)
  unresolved <- additions[is.na(tx_id), gene_symbol]
  if (length(unresolved) > 0) {
    stop("Manual additions could not be resolved to selected transcripts: ",
         paste(unique(unresolved), collapse = ", "))
  }
  unique(additions, by = c("gene_symbol", "tx_id", "tx_start", "tx_end"))
}

manual_file_targets <- read_manual_additions(manual_additions_file, analysis_dir)
if (nrow(manual_file_targets) > 0) {
  manual_file_targets[, `:=`(
    start_codon = NA_character_,
    stop_codon = NA_character_,
    strict_category_match = TRUE,
    selection_reason = paste0(
      "manual translon from ",
      basename(manual_additions_file),
      " (type=", manual_type, ", tx coordinates ",
      tx_start, "-", tx_end, ")"
    )
  )]
  fixed_manual_targets <- rbind(
    fixed_manual_targets,
    manual_file_targets[, .(
      gene_symbol, tx_id, tx_start, tx_end, start_codon, stop_codon,
      selection_reason, manual_type, strict_category_match
    )],
    fill = TRUE
  )
  fixed_manual_targets <- unique(
    fixed_manual_targets,
    by = c("gene_symbol", "tx_id", "tx_start", "tx_end")
  )
}

manual_targets <- unique(rbind(
  discovered_manual_targets,
  fixed_manual_targets[, .(gene_symbol, tx_id)]
))

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

load_translon_source <- function(paths, source) {
  path <- unique(paths)
  path <- path[file.exists(path)][1]
  if (is.na(path)) return(GRangesList())
  translons <- read_RDSQS(path)
  mcols(translons)$translon_source <- source
  translons
}

load_t_tc_translons <- function(df) {
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
  c(load_translon_source(t_paths, "T"), load_translon_source(tc_paths, "TC"))
}

annotate_tx_coordinates <- function(translons, cds, mrna) {
  if (length(translons) == 0) return(translons)

  keep <- names(translons) %chin% names(mrna)
  translons <- translons[keep]
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

genomic_string <- function(gr) {
  paste0(
    as.character(seqnames(gr)), ":", start(gr), "-", end(gr), ":",
    as.character(strand(gr)),
    collapse = ";"
  )
}

score_with_merged_rfp <- function(df_activity, candidates, mrna) {
  libs <- outputLibs(df_activity, type = "pshifted", output.mode = "envirlist")
  if (length(libs) < 1) {
    stop("No pshifted libraries were returned by outputLibs() for human_all_merged_l50.")
  }
  rfp <- libs[[1]]
  start_gr <- startSites(candidates, TRUE, TRUE, TRUE)
  stop_gr <- stopSites(candidates, TRUE, TRUE, TRUE)
  upstream_gr <- windowPerGroup(start_gr, mrna[names(start_gr)], 20, -2)
  downstream_gr <- windowPerGroup(stop_gr, mrna[names(stop_gr)], -2, 20)

  orf_stats <- ORFik:::coveragePerORFStatistics(candidates, rfp)
  upstream_stats <- ORFik:::coveragePerORFStatistics(upstream_gr, rfp)
  downstream_stats <- ORFik:::coveragePerORFStatistics(downstream_gr, rfp)
  start_reads <- countOverlapsW(start_gr, rfp, weight = "score")
  stop_reads <- countOverlapsW(stop_gr, rfp, weight = "score")

  predicted <- (orf_stats$mean > upstream_stats$mean * 1.3) &
    orf_stats$ORFScores > 2.5 &
    ((start_reads + 3) > orf_stats$median)

  data.table(
    start_reads = start_reads,
    stop_reads = stop_reads,
    predicted_by_detect_rule = predicted,
    orf_mean = orf_stats$mean,
    orf_median = orf_stats$median,
    orfscore = orf_stats$ORFScores,
    upstream_mean = upstream_stats$mean,
    downstream_mean = downstream_stats$mean
  )
}

df <- ORFik::read.experiment("all_samples-Homo_sapiens", validate = FALSE)
target_tx <- manual_targets$tx_id
cds <- loadRegion(df, "cds", names.keep = target_tx)
mrna <- loadRegion(df, "mrna", names.keep = target_tx)

message("Finding manual uoORF candidates with ORFik findORFs + categorize_and_filter_ORFs")
discovered_tx <- discovered_manual_targets$tx_id
discovered_cds <- cds[discovered_tx]
discovered_mrna <- mrna[discovered_tx]
orf_candidates <- findORFs(
  seqs = txSeqsFromFa(discovered_mrna, df, TRUE),
  longestORF = FALSE,
  startCodon = startDefinition(1),
  stopCodon = stopDefinition(1),
  minimumLength = 0
)
uo_candidates <- ORFik:::categorize_and_filter_ORFs(
  orf_candidates,
  ORF_categories_to_keep = "uoORF",
  cds = discovered_cds,
  mrna = discovered_mrna,
  map_to_gr = TRUE
)
names(uo_candidates) <- sub("_[0-9]+$", "", names(uo_candidates))
uo_candidates <- uo_candidates[names(uo_candidates) %chin% discovered_tx]
uo_candidates <- annotate_tx_coordinates(uo_candidates, discovered_cds, discovered_mrna)

if (length(uo_candidates) == 0) {
  stop("No candidate uoORFs were found for manual targets.")
}

message("Checking which candidate uoORFs are absent from T and TC")
t_tc <- load_t_tc_translons(df)
t_tc <- t_tc[names(t_tc) %chin% target_tx]
t_tc <- annotate_tx_coordinates(t_tc, cds, mrna)
t_tc_key <- if (length(t_tc) > 0) {
  paste(names(t_tc), mcols(t_tc)$tx_start, mcols(t_tc)$tx_end, sep = "|")
} else {
  character()
}

candidate_key <- paste(
  names(uo_candidates),
  mcols(uo_candidates)$tx_start,
  mcols(uo_candidates)$tx_end,
  sep = "|"
)

candidate_dt <- data.table(
  candidate_index = seq_along(uo_candidates),
  manual_input_type = "orfik_findORFs_candidate",
  gene_symbol = discovered_manual_targets$gene_symbol[
    match(names(uo_candidates), discovered_manual_targets$tx_id)
  ],
  tx_id = names(uo_candidates),
  category = as.character(mcols(uo_candidates)$category),
  tx_start = as.integer(mcols(uo_candidates)$tx_start),
  tx_end = as.integer(mcols(uo_candidates)$tx_end),
  width = as.integer(widthPerGroup(uo_candidates, FALSE)),
  genomic = vapply(seq_along(uo_candidates), function(i) {
    genomic_string(uo_candidates[[i]])
  }, character(1)),
  missing_from_T_TC = !candidate_key %chin% t_tc_key
)

message("Scoring manual uoORF candidates on human_all_merged_l50 pshifted coverage")
df_activity <- ORFik::read.experiment("human_all_merged_l50", validate = FALSE)
activity_scores <- score_with_merged_rfp(df_activity, uo_candidates, mrna)
candidate_dt <- cbind(candidate_dt, activity_scores)

candidate_dt[
  missing_from_T_TC == TRUE,
  manual_rank_by_start_reads := frank(-start_reads, ties.method = "first"),
  by = gene_symbol
]
candidate_dt[, manual_selected := missing_from_T_TC == TRUE & manual_rank_by_start_reads == 1]
candidate_dt[, selection_reason := fifelse(
  manual_selected,
  "highest merged-l50 start-site support among missing uoORF candidates",
  NA_character_
)]

fixed_iranges <- IRangesList(lapply(seq_len(nrow(fixed_manual_targets)), function(i) {
  IRanges(fixed_manual_targets$tx_start[[i]], fixed_manual_targets$tx_end[[i]])
}))
names(fixed_iranges) <- as.character(seq_len(nrow(fixed_manual_targets)))
fixed_translons <- pmapFromTranscriptF(
  fixed_iranges,
  mrna[fixed_manual_targets$tx_id],
  removeEmpty = TRUE
)
names(fixed_translons) <- fixed_manual_targets$tx_id
fixed_translons <- annotate_tx_coordinates(fixed_translons, cds, mrna)

fixed_key <- paste(
  names(fixed_translons),
  mcols(fixed_translons)$tx_start,
  mcols(fixed_translons)$tx_end,
  sep = "|"
)

fixed_dt <- data.table(
  candidate_index = NA_integer_,
  manual_input_type = "fixed_user_tx_coordinate",
  gene_symbol = fixed_manual_targets$gene_symbol,
  tx_id = names(fixed_translons),
  category = as.character(mcols(fixed_translons)$category),
  tx_start = as.integer(mcols(fixed_translons)$tx_start),
  tx_end = as.integer(mcols(fixed_translons)$tx_end),
  width = as.integer(widthPerGroup(fixed_translons, FALSE)),
  genomic = vapply(seq_along(fixed_translons), function(i) {
    genomic_string(fixed_translons[[i]])
  }, character(1)),
  missing_from_T_TC = !fixed_key %chin% t_tc_key,
  start_reads = NA_integer_,
  stop_reads = NA_integer_,
  predicted_by_detect_rule = NA,
  orf_mean = NA_real_,
  orf_median = NA_real_,
  orfscore = NA_real_,
  upstream_mean = NA_real_,
  downstream_mean = NA_real_,
  manual_rank_by_start_reads = NA_integer_,
  manual_type = fixed_manual_targets$manual_type,
  strict_category_match = fixed_manual_targets$strict_category_match,
  expected_category = fifelse(
    is.na(fixed_manual_targets$manual_type),
    NA_character_,
    expected_category_from_manual_type(fixed_manual_targets$manual_type)
  ),
  selection_reason = fixed_manual_targets$selection_reason,
  start_codon = fixed_manual_targets$start_codon,
  stop_codon = fixed_manual_targets$stop_codon
)
fixed_dt[, category_match := is.na(expected_category) | category == expected_category]
fixed_dt[, category_override_allowed :=
           manual_category_override_allowed(manual_type, category)]
fixed_dt[, manual_selected :=
           strict_category_match != TRUE |
             category_match == TRUE |
             category_override_allowed == TRUE]
fixed_dt[
  manual_selected == TRUE & category_match == FALSE & category_override_allowed == TRUE,
  selection_reason := paste0(
    selection_reason,
    "; selected with reviewed manual iORF override although ORFik categorized ",
    "the coordinates as ", category
  )
]
fixed_dt[manual_selected == FALSE, selection_reason := paste0(
  "not selected: manual type ", manual_type, " expects category ",
  expected_category, " but ORFik categorized tx coordinates ",
  tx_start, "-", tx_end, " as ", category, "; check coordinate system"
)]

candidate_dt[, start_codon := NA_character_]
candidate_dt[, stop_codon := NA_character_]
candidate_dt <- rbind(candidate_dt, fixed_dt, fill = TRUE)
setorder(candidate_dt, gene_symbol, tx_start, tx_end)

selected_indices <- candidate_dt[
  manual_selected == TRUE & manual_input_type == "orfik_findORFs_candidate",
  candidate_index
]
selected_discovered_translons <- uo_candidates[selected_indices]
names(selected_discovered_translons) <- names(uo_candidates)[selected_indices]
selected_fixed_indices <- which(fixed_dt$manual_selected == TRUE)
selected_fixed_translons <- fixed_translons[selected_fixed_indices]
manual_translons <- c(selected_discovered_translons, selected_fixed_translons)
selected_dt <- rbind(
  candidate_dt[manual_selected == TRUE & manual_input_type == "orfik_findORFs_candidate"],
  fixed_dt[manual_selected == TRUE],
  fill = TRUE
)

mcols(manual_translons)$manual_gene_symbol <- selected_dt$gene_symbol
mcols(manual_translons)$manual_tx_start <- selected_dt$tx_start
mcols(manual_translons)$manual_tx_end <- selected_dt$tx_end
mcols(manual_translons)$manual_selection_reason <- selected_dt$selection_reason
mcols(manual_translons)$manual_start_codon <- selected_dt$start_codon
mcols(manual_translons)$manual_stop_codon <- selected_dt$stop_codon
mcols(manual_translons)$manual_type <- selected_dt$manual_type
mcols(manual_translons)$manual_category_match <- selected_dt$category_match
mcols(manual_translons)$manual_category_override_allowed <- selected_dt$category_override_allowed
mcols(manual_translons)$manual_source_short <- "M"
mcols(manual_translons)$manual_source_name <- "manual"

saveRDS(manual_translons, manual_ranges_file)
fwrite(candidate_dt, candidate_manifest_file)
setorder(selected_dt, gene_symbol, tx_start, tx_end)
fwrite(selected_dt, selected_manifest_file)

message("Saved manual translon ranges: ", manual_ranges_file)
message("Saved manual candidate manifest: ", candidate_manifest_file)
message("Saved selected manual manifest: ", selected_manifest_file)
message("Selected manual translons: ",
        paste(selected_dt$gene_symbol, selected_dt$tx_start, selected_dt$tx_end,
              sep = ":", collapse = ", "))
