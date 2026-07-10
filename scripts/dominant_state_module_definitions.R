suppressPackageStartupMessages({
  library(data.table)
})

dominant_state_analysis_dir <- function() {
  if (dir.exists("dominant_cell_states")) "dominant_cell_states" else "."
}

dominant_state_results_dir <- function(analysis_dir = dominant_state_analysis_dir()) {
  env_results <- Sys.getenv("DOMINANT_CELL_STATES_RESULTS_DIR", unset = "")
  if (nzchar(env_results)) {
    if (grepl("^(/|~)", env_results)) {
      return(normalizePath(env_results, mustWork = FALSE))
    }
    return(normalizePath(file.path(analysis_dir, env_results),
                         mustWork = FALSE))
  }
  file.path(analysis_dir, "results")
}

dominant_state_curated_input_file <- function(analysis_dir, ...) {
  file.path(dominant_state_results_dir(analysis_dir), "curated_inputs", ...)
}

base_dominant_state_modules <- function() {
  data.table(
    dominant_state = c(
      "Proliferation / Cell cycle",
      "mTOR / Translation capacity (TOP program)",
      "ISR / ER stress (ATF4 axis)",
      "Interferon / antiviral",
      "Hypoxia / HIF1A program",
      "EMT / Mesenchymal shift",
      "Senescence / SASP",
      "OXPHOS (mitochondrial)",
      "Baseline (control)"
    ),
    up_genes = list(
      c("MKI67", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6",
        "MCM7", "CDK1", "PLK1", "RRM1", "RRM2", "FEN1", "PRIM2",
        "POLA1"),
      c("RPS6", "RPLP0", "RPL32", "RPL13A", "RPS3", "EIF3E", "EIF4B",
        "EEF2", "PABPC1", "NPM1", "NCL", "FBL", "MYC"),
      c("ATF4", "DDIT3", "ATF3", "ASNS", "PPP1R15A", "DDIT4", "CEBPB",
        "HERPUD1", "HSPA5", "XBP1"),
      c("STAT1", "ISG15", "MX1", "OAS1", "IFIT1", "IFIT3", "IFI44",
        "CXCL10"),
      c("CA9", "VEGFA", "EGLN3", "SLC2A1", "BNIP3", "LDHA", "PDK1",
        "NDRG1"),
      c("VIM", "FN1", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "CDH2", "MMP2",
        "THY1", "POSTN"),
      c("CDKN1A", "CDKN2A", "SERPINE1", "IGFBP7", "IGFBP5", "TIMP1",
        "TGFBI", "CXCL8", "IL6", "CCND1"),
      c("NDUFA1", "NDUFS1", "SDHB", "UQCRC1", "COX5A", "ATP5F1A",
        "ATP5MC1", "PHB2"),
      c("ACTB", "GAPDH", "HPRT1", "TBP", "PPIA", "UBC", "YWHAZ",
        "EEF1A1", "PSMB4", "RPS18", "RPL27", "RPL30", "RPS27A", "TPT1",
        "RPL13A", "B2M")
    ),
    down_genes = list(
      character(),
      character(),
      character(),
      character(),
      character(),
      c("CDH1", "EPCAM", "KRT8", "KRT18", "KRT19"),
      "LMNB1",
      c("SLC2A1", "LDHA", "PDK1"),
      character()
    )
  )
}

load_expanded_oxphos_gene_set <- function(analysis_dir = dominant_state_analysis_dir()) {
  gene_file <- dominant_state_curated_input_file(
    analysis_dir,
    "expanded_oxphos_marker_genes.csv"
  )
  if (!file.exists(gene_file)) {
    return(NULL)
  }

  genes <- fread(gene_file, showProgress = FALSE)
  required <- c("gene_symbol", "score_role")
  if (!all(required %chin% names(genes))) {
    stop("Expanded OXPHOS gene set is missing columns: ",
         paste(setdiff(required, names(genes)), collapse = ", "),
         call. = FALSE)
  }

  genes[, gene_symbol := toupper(trimws(gene_symbol))]
  genes[, score_role := tolower(trimws(score_role))]
  genes <- genes[
    nzchar(gene_symbol) &
      score_role %chin% c("oxphos_capacity_up", "glycolysis_hypoxia_down")
  ]

  up_genes <- sort(unique(genes[score_role == "oxphos_capacity_up",
                                gene_symbol]))
  down_genes <- sort(unique(genes[score_role == "glycolysis_hypoxia_down",
                                  gene_symbol]))
  if (length(up_genes) == 0 && length(down_genes) == 0) {
    return(NULL)
  }

  list(
    up_genes = up_genes,
    down_genes = down_genes,
    source_file = gene_file
  )
}

apply_expanded_oxphos_module <- function(modules,
                                         analysis_dir = dominant_state_analysis_dir()) {
  expanded <- load_expanded_oxphos_gene_set(analysis_dir)
  if (is.null(expanded)) {
    return(modules)
  }

  idx <- which(modules$dominant_state == "OXPHOS (mitochondrial)")
  if (length(idx) != 1) {
    stop("Expected exactly one OXPHOS module, found ", length(idx),
         call. = FALSE)
  }

  modules <- copy(modules)
  modules[idx, `:=`(
    up_genes = list(expanded$up_genes),
    down_genes = list(expanded$down_genes)
  )]
  attr(modules, "expanded_oxphos_marker_file") <- expanded$source_file
  attr(modules, "expanded_oxphos_up_genes") <- length(expanded$up_genes)
  attr(modules, "expanded_oxphos_down_genes") <- length(expanded$down_genes)
  modules
}

load_postviral_fatigue_gene_set <- function(analysis_dir = dominant_state_analysis_dir()) {
  gene_file <- dominant_state_curated_input_file(
    analysis_dir,
    "postviral_fatigue_dominant_state_genes.csv"
  )
  if (!file.exists(gene_file)) {
    return(data.table(
      State = character(),
      Gene = character(),
      Direction = character(),
      Module = character(),
      Tier = character(),
      Rationale = character()
    ))
  }

  genes <- fread(gene_file)
  required <- c("State", "Gene", "Direction", "Module", "Tier")
  if (!all(required %chin% names(genes))) {
    stop("Post-viral fatigue gene set is missing columns: ",
         paste(setdiff(required, names(genes)), collapse = ", "),
         call. = FALSE)
  }

  genes[, Gene := toupper(trimws(Gene))]
  genes[, Direction := tolower(trimws(Direction))]
  genes[, Tier := tolower(trimws(Tier))]
  genes <- genes[nzchar(Gene) & Direction %chin% c("up", "down")]
  unique(genes, by = c("Gene", "Direction", "Tier"))
}

dominant_state_symbol_aliases <- function() {
  c(CHOP = "DDIT3", DDX58 = "RIGI")
}

resolve_dominant_state_gene_symbols <- function(genes) {
  genes <- toupper(trimws(as.character(genes)))
  aliases <- dominant_state_symbol_aliases()
  out <- genes
  alias_hits <- out %chin% names(aliases)
  out[alias_hits] <- unname(aliases[out[alias_hits]])
  out[nzchar(out)]
}

load_additional_prior_gene_set <- function(analysis_dir = dominant_state_analysis_dir()) {
  priority_file <- dominant_state_curated_input_file(
    analysis_dir,
    "curated_additional_gene_priorities.csv"
  )
  if (!file.exists(priority_file)) {
    return(character())
  }

  genes <- fread(priority_file, showProgress = FALSE)
  if (!"gene_symbol" %chin% names(genes)) {
    stop("Additional prior gene table is missing `gene_symbol`: ",
         priority_file, call. = FALSE)
  }

  sort(unique(resolve_dominant_state_gene_symbols(genes$gene_symbol)))
}

load_dominant_state_modules <- function(analysis_dir = dominant_state_analysis_dir(),
                                        include_postviral = TRUE,
                                        postviral_tiers = "conservative") {
  modules <- base_dominant_state_modules()
  modules <- apply_expanded_oxphos_module(modules, analysis_dir)
  if (!include_postviral) return(modules)

  fatigue_genes <- load_postviral_fatigue_gene_set(analysis_dir)
  fatigue_score_genes <- fatigue_genes[Tier %chin% tolower(postviral_tiers)]
  if (nrow(fatigue_score_genes) == 0) {
    warning("No post-viral fatigue genes found for scoring tiers: ",
            paste(postviral_tiers, collapse = ", "))
    return(modules)
  }

  fatigue_module <- data.table(
    dominant_state = "Post-viral fatigue / ribosome stress",
    up_genes = list(sort(unique(fatigue_score_genes[Direction == "up", Gene]))),
    down_genes = list(sort(unique(fatigue_score_genes[Direction == "down", Gene])))
  )

  rbind(
    modules[dominant_state != "Baseline (control)"],
    fatigue_module,
    modules[dominant_state == "Baseline (control)"],
    fill = TRUE
  )
}

dominant_state_target_genes <- function(modules,
                                        analysis_dir = dominant_state_analysis_dir(),
                                        include_postviral_auxiliary = TRUE,
                                        include_additional_priors = TRUE) {
  target_genes <- unique(unlist(c(modules$up_genes, modules$down_genes),
                                use.names = FALSE))
  if (include_postviral_auxiliary) {
    fatigue_genes <- load_postviral_fatigue_gene_set(analysis_dir)
    target_genes <- unique(c(target_genes, fatigue_genes$Gene))
  }
  if (include_additional_priors) {
    target_genes <- unique(c(target_genes, load_additional_prior_gene_set(analysis_dir)))
  }
  sort(unique(resolve_dominant_state_gene_symbols(target_genes)))
}
