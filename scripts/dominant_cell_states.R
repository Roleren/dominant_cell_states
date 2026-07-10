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

analysis_dir <- if (dir.exists("dominant_cell_states")) "dominant_cell_states" else "."

message("Dominant-state QC annotation pseudo code:")
message("  Step 0: expr_norm = log2(count + 1)")
message("  Step 1: baseline = mean(expr_norm[control_genes])")
message("  Step 1: gene_z = expr_norm - baseline")
message("  Step 2: score = mean(gene_z[up_genes])")
message("  Step 2: signed_score = mean(gene_z[up_genes]) - mean(gene_z[down_genes])")

module_definition_file <- file.path(
  analysis_dir,
  "scripts",
  "dominant_state_module_definitions.R"
)
if (!file.exists(module_definition_file)) {
  stop("Missing dominant-state module definition file: ", module_definition_file)
}
source(module_definition_file)
dominant_state_modules <- load_dominant_state_modules(analysis_dir)

df <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
counts <- countTable(df, "cds", "fpkm")
symbols <- symbols(df)

strip_version <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

gene_list_text <- function(x) {
  vapply(x, paste, collapse = ", ", FUN.VALUE = character(1))
}

tx_column <- function(symbols_dt) {
  preferred <- c("ensembl_tx_name", "tx_name", "transcript", "transcript_id")
  preferred <- preferred[preferred %in% names(symbols_dt)]
  if (length(preferred) > 0) return(preferred[1])

  candidates <- grep("tx|transcript|value", names(symbols_dt),
                     ignore.case = TRUE, value = TRUE)
  if (length(candidates) != 1) {
    stop("Could not uniquely identify transcript column in symbols table. ",
         "Columns: ", paste(names(symbols_dt), collapse = ", "))
  }
  candidates
}

symbol_column <- function(symbols_dt) {
  preferred <- c("external_gene_name", "gene_symbol", "symbol", "hgnc_symbol",
                 "gene_name", "label")
  preferred <- preferred[preferred %in% names(symbols_dt)]
  if (length(preferred) > 0) return(preferred[1])

  candidates <- grep("symbol|gene.*name|external|label", names(symbols_dt),
                     ignore.case = TRUE, value = TRUE)
  if (length(candidates) < 1) {
    stop("Could not identify gene-symbol column in symbols table. Columns: ",
         paste(names(symbols_dt), collapse = ", "))
  }
  candidates[1]
}

counts_to_gene_matrix <- function(count_table, symbols_dt, target_genes) {
  target_genes <- unique(toupper(target_genes))
  count_dt <- as.data.table(count_table)
  feature_id <- rownames(count_table)

  if (is.null(feature_id) || length(feature_id) != nrow(count_dt) ||
      all(grepl("^[0-9]+$", feature_id))) {
    candidate_id_cols <- intersect(
      c("id", "transcript", "transcript_id", "tx_name", "ensembl_tx_name"),
      names(count_dt)
    )
    if (length(candidate_id_cols) == 0) {
      stop("Counts table has no usable row names or transcript ID column.")
    }
    feature_id <- count_dt[[candidate_id_cols[1]]]
    count_dt[, (candidate_id_cols[1]) := NULL]
  }

  sample_names <- names(count_dt)
  count_matrix <- as.matrix(count_dt)
  storage.mode(count_matrix) <- "numeric"
  if (ncol(count_matrix) == 0) stop("Counts table has no sample columns.")
  if (is.null(sample_names)) {
    sample_names <- paste0("sample_", seq_len(ncol(count_matrix)))
  }
  colnames(count_matrix) <- sample_names

  feature_map <- data.table(
    row_index = seq_along(feature_id),
    feature_id = as.character(feature_id)
  )
  feature_map[, feature_id_no_version := strip_version(feature_id)]
  feature_map[, gene_symbol := NA_character_]

  if (nrow(symbols_dt) > 0) {
    tx_col <- tx_column(symbols_dt)
    sym_col <- symbol_column(symbols_dt)
    symbol_map <- as.data.table(symbols_dt)[
      !is.na(get(tx_col)) & !is.na(get(sym_col)) &
        nzchar(get(sym_col)) & toupper(get(sym_col)) %in% target_genes,
      .(
        feature_id = as.character(get(tx_col)),
        gene_symbol = toupper(as.character(get(sym_col)))
      )
    ]
    symbol_map <- unique(symbol_map)
    symbol_map <- symbol_map[!duplicated(feature_id)]

    feature_map[, gene_symbol := symbol_map$gene_symbol[
      match(feature_id, symbol_map$feature_id)
    ]]

    missing_symbol <- is.na(feature_map$gene_symbol) | !nzchar(feature_map$gene_symbol)
    if (any(missing_symbol)) {
      symbol_map[, feature_id_no_version := strip_version(feature_id)]
      symbol_map <- unique(symbol_map[, .(feature_id_no_version, gene_symbol)])
      symbol_map <- symbol_map[!duplicated(feature_id_no_version)]
      feature_map[missing_symbol, gene_symbol := symbol_map$gene_symbol[
        match(feature_id_no_version, symbol_map$feature_id_no_version)
      ]]
    }
  }

  direct_gene_rows <- toupper(feature_map$feature_id) %in% target_genes
  feature_map[is.na(gene_symbol) & direct_gene_rows,
              gene_symbol := toupper(feature_id)]

  feature_map <- feature_map[!is.na(gene_symbol) & nzchar(gene_symbol)]
  if (nrow(feature_map) == 0) {
    stop("No count rows could be mapped to gene symbols.")
  }

  marker_matrix <- count_matrix[feature_map$row_index, , drop = FALSE]
  gene_groups <- split(seq_len(nrow(marker_matrix)), feature_map$gene_symbol)
  gene_matrix <- do.call(rbind, lapply(gene_groups, function(i) {
    colMeans(marker_matrix[i, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(gene_matrix) <- names(gene_groups)
  colnames(gene_matrix) <- sample_names
  gene_matrix
}

score_module <- function(gene_z, up_genes, down_genes = character()) {
  up_present <- intersect(up_genes, rownames(gene_z))
  down_present <- intersect(down_genes, rownames(gene_z))

  if (length(up_present) == 0) {
    up_score <- rep(NA_real_, ncol(gene_z))
  } else {
    up_score <- colMeans(gene_z[up_present, , drop = FALSE], na.rm = TRUE)
  }

  if (length(down_genes) == 0) {
    down_score <- rep(0, ncol(gene_z))
  } else if (length(down_present) == 0) {
    down_score <- rep(NA_real_, ncol(gene_z))
  } else {
    down_score <- colMeans(gene_z[down_present, , drop = FALSE], na.rm = TRUE)
  }

  up_score - down_score
}

target_genes <- dominant_state_target_genes(
  dominant_state_modules,
  analysis_dir,
  include_postviral_auxiliary = TRUE
)
gene_expr <- counts_to_gene_matrix(counts, symbols, target_genes)
sample_names <- colnames(gene_expr)
run_ids <- runIDs(df)
if (length(run_ids) != ncol(gene_expr)) {
  stop("runIDs(df) length (", length(run_ids),
       ") does not match count-table sample columns (", ncol(gene_expr), ").")
}
sample_runs <- data.table(
  sample_index = seq_len(ncol(gene_expr)),
  Run = as.character(run_ids),
  sample = sample_names
)

expr_norm <- log2(gene_expr + 1)
baseline_genes <- dominant_state_modules[
  dominant_state == "Baseline (control)", up_genes[[1]]
]
baseline_present <- intersect(baseline_genes, rownames(expr_norm))

if (length(baseline_present) < 5) {
  stop("Too few baseline genes were found in the count table: ",
       length(baseline_present), " / ", length(baseline_genes))
}

baseline <- colMeans(expr_norm[baseline_present, , drop = FALSE], na.rm = TRUE)
gene_z <- sweep(expr_norm, 2, baseline, "-")

dominant_state_scores <- rbindlist(lapply(seq_len(nrow(dominant_state_modules)), function(i) {
  state <- dominant_state_modules[i]
  score <- score_module(gene_z, state$up_genes[[1]], state$down_genes[[1]])
  data.table(
    sample_index = seq_along(score),
    Run = as.character(run_ids),
    sample = sample_names,
    dominant_state = state$dominant_state,
    score = as.numeric(score),
    up_genes_present = length(intersect(state$up_genes[[1]], rownames(gene_z))),
    up_genes_total = length(state$up_genes[[1]]),
    down_genes_present = length(intersect(state$down_genes[[1]], rownames(gene_z))),
    down_genes_total = length(state$down_genes[[1]])
  )
}))

dominant_state_qc_annotations <- dominant_state_scores[
  dominant_state != "Baseline (control)",
  .SD[which.max(score)],
  by = sample_index
][
  ,
  .(sample_index, Run, sample, dominant_state, dominant_state_score = score)
]

dominant_state_scores_wide <- dcast(
  dominant_state_scores[, .(sample_index, dominant_state, score)],
  sample_index ~ dominant_state,
  value.var = "score"
)

dominant_state_qc_annotations <- merge(
  dominant_state_qc_annotations,
  dominant_state_scores_wide,
  by = "sample_index",
  sort = FALSE
)
setorder(dominant_state_qc_annotations, sample_index)
dominant_state_qc_annotations[, sample_index := NULL]
setcolorder(
  dominant_state_qc_annotations,
  c("Run", setdiff(names(dominant_state_qc_annotations), "Run"))
)

marker_coverage <- dominant_state_modules[
  ,
  .(
    dominant_state,
    up_genes = gene_list_text(up_genes),
    down_genes = gene_list_text(down_genes),
    up_genes_present = vapply(up_genes, function(x) {
      length(intersect(x, rownames(gene_z)))
    }, integer(1)),
    up_genes_total = lengths(up_genes),
    down_genes_present = vapply(down_genes, function(x) {
      length(intersect(x, rownames(gene_z)))
    }, integer(1)),
    down_genes_total = lengths(down_genes)
  )
]

stopifnot(nrow(gene_expr) > 0)
stopifnot(ncol(gene_expr) > 0)
stopifnot(identical(dominant_state_qc_annotations$Run,
                    as.character(run_ids)))
stopifnot(nrow(dominant_state_scores) ==
            nrow(dominant_state_modules) * ncol(gene_expr))
stopifnot(all(is.finite(dominant_state_scores[
  dominant_state != "Baseline (control)", score
])))
stopifnot(all(marker_coverage[
  dominant_state != "Baseline (control)", up_genes_present
] > 0))
stopifnot(all(marker_coverage[
  down_genes_total > 0, down_genes_present
] > 0))

message("\nMarker definition / coverage table:")
print(marker_coverage)

message("\nDominant-state QC annotation result table:")
print(dominant_state_qc_annotations)

invisible(list(
  modules = dominant_state_modules,
  marker_coverage = marker_coverage,
  gene_counts = gene_expr,
  expr_norm = expr_norm,
  gene_z = gene_z,
  scores = dominant_state_scores,
  sample_runs = sample_runs,
  qc_annotations = dominant_state_qc_annotations
))
