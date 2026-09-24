dominant_loader <- c(
  file.path("scripts", "dominant_ribocrypt_loader.R"),
  file.path("dominant_cell_states", "scripts", "dominant_ribocrypt_loader.R")
)
dominant_loader <- dominant_loader[file.exists(dominant_loader)][1]
if (!is.na(dominant_loader)) {
  source(dominant_loader)
  dominant_load_ribocrypt_if_available()
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

analysis_dir <- if (dir.exists("dominant_cell_states")) "dominant_cell_states" else "."
output_dir <- file.path(analysis_dir, "dominant_state_auxiliary_modules")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

module_definition_file <- file.path(
  analysis_dir, "scripts", "dominant_state_module_definitions.R"
)
if (!file.exists(module_definition_file)) {
  stop("Missing module definition file: ", module_definition_file, call. = FALSE)
}
source(module_definition_file)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)

state_file <- file.path(analysis_dir, "human_dominant_cell_states_clean_cds.csv")
expression_file <- file.path(
  analysis_dir, "human_dominant_cell_states_clean_cds_gene_expression.csv"
)
metadata_file <- local({
  candidates <- c(
    Sys.getenv("DOMINANT_METADATA_FILE", unset = NA_character_),
    "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv",
    path.expand("~/livemount/Bio_data/NGS_pipeline/metadata_done_samples_extended_qc.csv")
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) existing[[1]] else candidates[[1]]
})
atlas_file <- file.path(analysis_dir, "dominant_rdg_atlas",
                        "dominant_rdg_atlas_cards.csv")

if (!file.exists(state_file) || !file.exists(expression_file)) {
  stop("Auxiliary module screen requires clean-CDS state/expression files.",
       call. = FALSE)
}

candidate_marker_modules <- function() {
  rbindlist(list(
    data.table(
      candidate_state = "Proteostasis / Heat shock",
      gene_symbol = c("HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1", "HSPH1",
                      "DNAJA1", "DNAJB1", "HSPB1", "BAG3", "HSF1",
                      "HSPA5", "HERPUD1"),
      direction = "up",
      marker_group = c(rep("cytosolic_chaperone", 10), "er_chaperone",
                       "er_chaperone"),
      rationale = "Heat-shock/proteostasis marker; expected to reshape translation during proteotoxic stress."
    ),
    data.table(
      candidate_state = "Autophagy / Lysosome (TFEB)",
      gene_symbol = c("SQSTM1", "MAP1LC3B", "GABARAPL1", "WIPI1", "ATG5",
                      "ATG7", "BECN1", "LAMP1", "LAMP2", "CTSB", "CTSD",
                      "TFEB", "TFE3", "DDIT4", "SESN2", "RPS6", "MYC"),
      direction = c(rep("up", 15), "down", "down"),
      marker_group = c(rep("autophagy_lysosome", 13),
                       rep("stress_mtor_inhibition", 2),
                       rep("growth_translation_capacity", 2)),
      rationale = "Autophagy/lysosome and low-mTOR state; useful for nutrient, proteostasis, and post-viral metabolism contexts."
    ),
    data.table(
      candidate_state = "NRF2 / Oxidative-ferroptosis defense",
      gene_symbol = c("HMOX1", "NQO1", "GCLC", "GCLM", "SLC7A11",
                      "TXNRD1", "SRXN1", "FTH1", "FTL", "GPX4",
                      "SQSTM1", "ATF4", "DDIT4", "SESN2", "GDF15"),
      direction = "up",
      marker_group = c(rep("nrf2_redox", 11), rep("isr_redox_cross_talk", 4)),
      rationale = "Oxidative stress, antioxidant capacity, and ferroptosis-defense marker set."
    ),
    data.table(
      candidate_state = "DNA damage / p53 checkpoint",
      gene_symbol = c("CDKN1A", "MDM2", "GADD45A", "DDB2", "BBC3",
                      "PMAIP1", "FAS", "RRM2B", "BTG2", "SESN2", "DDIT4",
                      "CCND1", "MYC", "CDK1", "PLK1", "MKI67"),
      direction = c(rep("up", 11), rep("down", 5)),
      marker_group = c(rep("p53_checkpoint", 11), rep("proliferation_down", 5)),
      rationale = "p53/DNA-damage arrest; should separate checkpoint arrest from simple low proliferation."
    ),
    data.table(
      candidate_state = "Antigen presentation / immune composition",
      gene_symbol = c("HLA-DRA", "HLA-DRB1", "CD74", "B2M", "TAP1",
                      "TAP2", "PSMB8", "PSMB9", "CIITA", "HLA-A",
                      "HLA-B", "HLA-C", "CD79A", "MS4A1"),
      direction = "up",
      marker_group = c(rep("antigen_presentation", 12), rep("b_cell", 2)),
      rationale = "MHC/antigen-presentation and immune-composition state; important confounder for tissue, tumor, and infection contrasts."
    ),
    data.table(
      candidate_state = "Secretory / ER proteostasis",
      gene_symbol = c("HSPA5", "XBP1", "HERPUD1", "ATF4", "DDIT3",
                      "DDIT4", "PDIA4", "CALR", "CANX", "DNAJB9",
                      "SEL1L", "EDEM1", "ERN1", "ATF6"),
      direction = "up",
      marker_group = c(rep("upr_er_chaperone", 3),
                       rep("isr_er_cross_talk", 3),
                       rep("erad_secretory", 8)),
      rationale = "ER protein-folding and secretory-load state, separated from the narrower ATF4-axis ISR score."
    ),
    data.table(
      candidate_state = "Mitochondrial UPR / mitonuclear stress",
      gene_symbol = c("ATF4", "DDIT3", "TRIB3", "GDF15", "DDIT4", "SESN2",
                      "HSPD1", "HSPE1", "HSPA9", "LONP1", "CLPP",
                      "NDUFS1", "NDUFA9", "SDHB", "UQCRC1", "COX5A",
                      "ATP5F1A", "TFAM", "PPARGC1A"),
      direction = c(rep("up", 11), rep("down", 8)),
      marker_group = c(rep("mitochondrial_stress_response", 11),
                       rep("mitochondrial_capacity_down", 8)),
      rationale = "Mitonuclear stress program; a more mechanistic split of the broad OXPHOS/post-viral fatigue signal."
    ),
    data.table(
      candidate_state = "ISR kinase / initiation control",
      gene_symbol = c("EIF2AK1", "EIF2AK2", "EIF2AK3", "EIF2AK4", "GCN1",
                      "PPP1R15B", "EIF2S1", "EIF2B1", "EIF2B2", "EIF2B3",
                      "EIF2B4", "EIF2B5"),
      direction = c(rep("up", 6), rep("down", 6)),
      marker_group = c(rep("stress_kinase_feedback", 6),
                       rep("ternary_complex_capacity", 6)),
      rationale = "Attempts to separate upstream ISR kinase/initiation machinery from the downstream ATF4 marker output."
    ),
    data.table(
      candidate_state = "RNA sensing / RNase L / stress granule",
      gene_symbol = c("RIGI", "MAVS", "TBK1", "IRF3", "OAS2", "OAS3",
                      "RNASEL", "IFIT2", "ZC3HAV1", "DDX3X", "G3BP1",
                      "G3BP2", "TIA1", "TIAL1", "IFIT1", "ISG15"),
      direction = "up",
      marker_group = c(rep("rna_sensing", 4), rep("rnasel_oas", 3),
                       rep("antiviral_decay", 2), rep("stress_granule", 5),
                       rep("isg_anchor", 2)),
      rationale = "Viral RNA sensing, RNase-L/OAS, antiviral RNA decay, and stress-granule context missing from the original ISG-only antiviral state."
    ),
    data.table(
      candidate_state = "Host mRNA export / initiation shutoff",
      gene_symbol = c("IFIT1", "ISG15", "ATF4", "DDIT3", "NUP98", "RAE1",
                      "NXF1", "ALYREF", "EIF4G1", "EIF4G2", "EIF4E",
                      "EIF4A1", "EIF4A2", "PABPN1", "RPS6", "RPLP0",
                      "MYC"),
      direction = c(rep("up", 4), rep("down", 13)),
      marker_group = c(rep("antiviral_isr_anchor", 4),
                       rep("mrna_export_capacity", 4),
                       rep("initiation_capacity", 6),
                       rep("translation_growth_capacity", 3)),
      rationale = "Host-shutoff-like state combining antiviral/ISR pressure with low mRNA export, initiation, and translation-growth capacity."
    ),
    data.table(
      candidate_state = "Ribosome quality control / collision",
      gene_symbol = c("ZNF598", "GIGYF2", "EIF2A", "PELO", "HBS1L",
                      "ASCC3", "NEMF", "LTN1", "GCN1", "EIF2AK4"),
      direction = "up",
      marker_group = c(rep("collision_sensor_decay", 2),
                       rep("alternate_initiation_rescue", 1),
                       rep("ribosome_rescue", 2),
                       rep("quality_control_complex", 3),
                       rep("gcn2_collision_axis", 2)),
      rationale = "Collision/no-go/ribosome-rescue machinery that could explain branch-specific accumulation or stalled coverage."
    ),
    data.table(
      candidate_state = "DecodeME GWAS-prior expression layer",
      gene_symbol = c("BTN2A2", "OLFM4", "RABGAP1L", "ZNFX1", "FBXL4",
                      "CA10", "ARFGEF2", "CSE1L"),
      direction = "up",
      marker_group = "gwas_prior_gene",
      rationale = "ME/CFS GWAS-prior genes tracked as an exploratory expression/context layer, not as a mechanistic dominant-state claim."
    ),
    data.table(
      candidate_state = "Host shutoff / antiviral translation repression",
      gene_symbol = c("IFIT1", "IFIT3", "IFITM1", "ISG15", "MX1", "OAS1",
                      "STAT1", "IRF7", "ATF4", "DDIT3", "RPS6", "RPLP0",
                      "RPL32", "RPS3", "EIF4B", "EEF2", "PABPC1", "MYC"),
      direction = c(rep("up", 10), rep("down", 8)),
      marker_group = c(rep("antiviral_isr", 10),
                       rep("translation_capacity_down", 8)),
      rationale = "Viral-host-shutoff-like score combining antiviral induction with reduced translation-capacity markers."
    )
  ), fill = TRUE)[
    ,
    `:=`(
      gene_symbol = toupper(gene_symbol),
      direction = tolower(direction)
    )
  ]
}

read_gene_expression_matrix <- function(path) {
  dt <- fread(path, showProgress = FALSE)
  if (!"gene_symbol" %chin% names(dt)) {
    stop("Expression table lacks gene_symbol: ", path, call. = FALSE)
  }
  gene_symbols <- toupper(as.character(dt$gene_symbol))
  expr <- as.matrix(dt[, setdiff(names(dt), "gene_symbol"), with = FALSE])
  storage.mode(expr) <- "numeric"
  rownames(expr) <- gene_symbols
  expr
}

mad_or_sd <- function(x) {
  spread <- mad(x, center = median(x, na.rm = TRUE), na.rm = TRUE)
  if (!is.finite(spread) || spread <= 0) spread <- sd(x, na.rm = TRUE)
  if (!is.finite(spread) || spread <= 0) 1 else spread
}

percentile_z <- function(x) {
  out <- rep(NA_real_, length(x))
  valid <- is.finite(x)
  if (!any(valid)) return(out)
  p <- (rank(x[valid], ties.method = "average") - 0.5) / sum(valid)
  out[valid] <- qnorm(pmin(pmax(p, 1e-6), 1 - 1e-6))
  out
}

score_candidate_module <- function(module_genes, gene_z, run_ids, sample_names) {
  up <- module_genes[direction == "up", unique(gene_symbol)]
  down <- module_genes[direction == "down", unique(gene_symbol)]
  up_present <- intersect(up, rownames(gene_z))
  down_present <- intersect(down, rownames(gene_z))
  up_score <- if (length(up_present) > 0) {
    colMeans(gene_z[up_present, , drop = FALSE], na.rm = TRUE)
  } else {
    rep(NA_real_, ncol(gene_z))
  }
  down_score <- if (length(down) > 0 && length(down_present) > 0) {
    colMeans(gene_z[down_present, , drop = FALSE], na.rm = TRUE)
  } else if (length(down) == 0) {
    rep(0, ncol(gene_z))
  } else {
    rep(NA_real_, ncol(gene_z))
  }
  signed_score <- as.numeric(up_score - down_score)
  center <- median(signed_score, na.rm = TRUE)
  spread <- mad_or_sd(signed_score)
  data.table(
    Run = run_ids,
    sample = sample_names,
    candidate_state = module_genes$candidate_state[1],
    signed_score = signed_score,
    robust_z = (signed_score - center) / spread,
    percentile_z = percentile_z(signed_score),
    percentile = pnorm(percentile_z(signed_score)),
    up_arm_score = as.numeric(up_score),
    down_arm_score = if (length(down) > 0) as.numeric(down_score) else NA_real_,
    up_genes_present = length(up_present),
    up_genes_total = length(up),
    down_genes_present = length(down_present),
    down_genes_total = length(down)
  )
}

cohen_d <- function(x, g) {
  x1 <- x[g]
  x0 <- x[!g]
  n1 <- length(x1)
  n0 <- length(x0)
  if (n1 < 2 || n0 < 2) return(NA_real_)
  s1 <- var(x1, na.rm = TRUE)
  s0 <- var(x0, na.rm = TRUE)
  pooled <- sqrt(((n1 - 1) * s1 + (n0 - 1) * s0) / (n1 + n0 - 2))
  if (!is.finite(pooled) || pooled <= 0) return(NA_real_)
  (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / pooled
}

safe_fisher_p <- function(a, b, c, d) {
  out <- tryCatch(
    fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value,
    error = function(e) NA_real_
  )
  as.numeric(out)
}

save_plot <- function(plot, name, width, height) {
  ggsave(file.path(figure_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 220, bg = "white")
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), plot,
         width = width, height = height, bg = "white")
  if (exists("dominant_save_ggplotly", mode = "function") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      plot, file.path(figure_dir, paste0(name, ".html")), title = name
    )
  }
}

message("Scoring auxiliary dominant-state modules")
expr <- read_gene_expression_matrix(expression_file)
state_dt <- fread(state_file, showProgress = FALSE)
if (!all(c("Run", "sample") %chin% names(state_dt))) {
  stop("Clean-CDS state table lacks Run/sample columns.", call. = FALSE)
}
if (!all(state_dt$Run %chin% colnames(expr))) {
  stop("Expression matrix is missing Run columns from the state table.",
       call. = FALSE)
}
expr <- expr[, state_dt$Run, drop = FALSE]

modules <- load_dominant_state_modules(analysis_dir)
baseline_genes <- modules[dominant_state == "Baseline (control)", up_genes[[1]]]
baseline_present <- intersect(toupper(baseline_genes), rownames(expr))
if (length(baseline_present) < 5) {
  stop("Too few baseline genes available for auxiliary scoring.", call. = FALSE)
}

expr_norm <- log2(expr + 1)
baseline <- colMeans(expr_norm[baseline_present, , drop = FALSE], na.rm = TRUE)
gene_z <- sweep(expr_norm, 2, baseline, "-")

candidate_genes <- candidate_marker_modules()
available_genes <- rownames(expr)
candidate_genes[, available_in_clean_cds := gene_symbol %chin% available_genes]
candidate_genes[
  ,
  `:=`(
    state_genes_present = sum(available_in_clean_cds),
    state_genes_total = .N
  ),
  by = candidate_state
]
candidate_genes[, state_gene_fraction_present :=
                  state_genes_present / pmax(state_genes_total, 1)]

scores <- rbindlist(
  lapply(split(candidate_genes, candidate_genes$candidate_state),
         score_candidate_module, gene_z = gene_z,
         run_ids = state_dt$Run, sample_names = state_dt$sample),
  fill = TRUE
)

existing_state_cols <- intersect(modules$dominant_state, names(state_dt))
existing_state_cols <- setdiff(existing_state_cols, "Baseline (control)")
existing_long <- melt(
  state_dt[, c("Run", existing_state_cols), with = FALSE],
  id.vars = "Run",
  variable.name = "existing_state",
  value.name = "existing_score",
  variable.factor = FALSE
)
correlations <- merge(
  scores[, .(Run, candidate_state, signed_score)],
  existing_long,
  by = "Run",
  allow.cartesian = TRUE
)[
  ,
  .(
    spearman_rho = suppressWarnings(cor(signed_score, existing_score,
                                        method = "spearman",
                                        use = "pairwise.complete.obs")),
    pearson_r = suppressWarnings(cor(signed_score, existing_score,
                                     method = "pearson",
                                     use = "pairwise.complete.obs"))
  ),
  by = .(candidate_state, existing_state)
]
correlations[!is.finite(spearman_rho), spearman_rho := NA_real_]
correlations[!is.finite(pearson_r), pearson_r := NA_real_]

top_existing <- correlations[
  ,
  .SD[which.max(abs(spearman_rho))],
  by = candidate_state
][
  ,
  .(
    candidate_state,
    nearest_existing_state = existing_state,
    max_abs_existing_spearman = abs(spearman_rho),
    nearest_existing_spearman = spearman_rho
  )
]

score_summary <- scores[
  ,
  .(
    up_genes_present = first(up_genes_present),
    up_genes_total = first(up_genes_total),
    down_genes_present = first(down_genes_present),
    down_genes_total = first(down_genes_total),
    genes_present = first(up_genes_present + down_genes_present),
    genes_total = first(up_genes_total + down_genes_total),
    gene_fraction_present = first((up_genes_present + down_genes_present) /
                                    pmax(up_genes_total + down_genes_total, 1)),
    score_median = median(signed_score, na.rm = TRUE),
    score_mad = mad_or_sd(signed_score),
    score_q95 = as.numeric(quantile(signed_score, 0.95, na.rm = TRUE)),
    score_max = max(signed_score, na.rm = TRUE),
    high_percentile_samples = sum(percentile >= 0.95, na.rm = TRUE),
    top_decile_samples = sum(percentile >= 0.90, na.rm = TRUE)
  ),
  by = candidate_state
]
score_summary <- merge(score_summary, top_existing, by = "candidate_state",
                       all.x = TRUE, sort = FALSE)
score_summary[, novelty_score :=
                pmax(0, 1 - fifelse(is.na(max_abs_existing_spearman), 0,
                                     max_abs_existing_spearman))]

metadata_enrichment <- data.table()
numeric_metadata <- data.table()
if (file.exists(metadata_file)) {
  metadata <- fread(metadata_file, showProgress = FALSE)
  if ("Run" %chin% names(metadata)) {
    merged <- merge(
      scores,
      metadata,
      by = "Run",
      all.x = TRUE,
      sort = FALSE
    )
    categorical_axes <- intersect(
      c("CELL_LINE", "TISSUE", "CONDITION", "GENE", "Cancer_type",
        "Cell_model", "Cell_type", "Organ_system", "Sex", "Life_stage",
        "AUTHOR", "INHIBITOR", "FRACTION", "LIBRARYTYPE"),
      names(merged)
    )
    metadata_enrichment <- rbindlist(lapply(categorical_axes, function(axis) {
      dt <- merged[
        !is.na(get(axis)) & nzchar(trimws(as.character(get(axis)))) &
          is.finite(signed_score)
      ]
      if (nrow(dt) == 0) return(data.table())
      dt[, axis_value := trimws(as.character(get(axis)))]
      dt <- dt[axis_value != "MISSING"]
      if (nrow(dt) == 0) return(data.table())
      state_total <- dt[, .(
        total_n = .N,
        total_top = sum(percentile >= 0.90, na.rm = TRUE)
      ), by = candidate_state]
      out <- dt[
        ,
        .(
          n_level = .N,
          n_top_decile = sum(percentile >= 0.90, na.rm = TRUE),
          mean_score = mean(signed_score, na.rm = TRUE),
          median_score = median(signed_score, na.rm = TRUE)
        ),
        by = .(candidate_state, axis_value)
      ][n_level >= 5]
      if (nrow(out) == 0) return(data.table())
      out <- merge(out, state_total, by = "candidate_state",
                   all.x = TRUE, sort = FALSE)
      out[
        ,
        `:=`(
          n_other = total_n - n_level,
          n_other_top = total_top - n_top_decile
        )
      ]
      out <- out[n_other >= 5]
      if (nrow(out) == 0) return(data.table())
      out[
        ,
        `:=`(
          axis = axis,
          mean_other = NA_real_,
          mean_delta = NA_real_,
          cohen_d = NA_real_,
          top_fraction = n_top_decile / pmax(n_level, 1),
          top_fraction_other = n_other_top / pmax(n_other, 1),
          fisher_p = NA_real_
        )
      ]
      for (i in seq_len(nrow(out))) {
        sub <- dt[candidate_state == out$candidate_state[i]]
        in_level <- sub$axis_value == out$axis_value[i]
        out$mean_other[i] <- mean(sub$signed_score[!in_level], na.rm = TRUE)
        out$mean_delta[i] <- out$mean_score[i] - out$mean_other[i]
        out$cohen_d[i] <- cohen_d(sub$signed_score, in_level)
        out$fisher_p[i] <- safe_fisher_p(
          out$n_top_decile[i],
          out$n_level[i] - out$n_top_decile[i],
          out$n_other_top[i],
          out$n_other[i] - out$n_other_top[i]
        )
      }
      out[
        ,
        .(
          candidate_state,
          axis,
          axis_value,
          n_level,
          n_top_decile,
          top_fraction,
          top_fraction_other,
          top_or = ((n_top_decile + 0.5) /
                      (n_level - n_top_decile + 0.5)) /
            ((n_other_top + 0.5) /
               (n_other - n_other_top + 0.5)),
          mean_score,
          mean_other,
          mean_delta,
          cohen_d,
          fisher_p
        )
      ]
    }), fill = TRUE)
    if (nrow(metadata_enrichment) > 0) {
      metadata_enrichment[, q := p.adjust(fisher_p, method = "BH"),
                          by = candidate_state]
      metadata_enrichment[, abs_cohen_d := abs(cohen_d)]
      setorder(metadata_enrichment, candidate_state, q, -abs_cohen_d)
    }

    numeric_axes <- names(merged)[vapply(merged, is.numeric, logical(1))]
    numeric_axes <- setdiff(
      intersect(
        numeric_axes,
        c("YEAR", "avgLength", "spots", "bases", "size_MB",
          "total mapped reads %", "Uniquely mapped reads %",
          "ProteinCoding_Pct", "rRNA_Pct", "Active_PC_Genes_min10",
          "Active_PC_Genes_min100", "Median_Expressed_PC_Level",
          "Sense_Mapping_Pct", "Frame_usage_0", "Frame_usage_1",
          "Frame_usage_2", "top_readlength", "drop_off", "drop_in",
          "drop_in_off_ratio", "tis_peak_strength", "tts_peak_strength",
          "pp_peak_strength")
      ),
      c("signed_score", "robust_z", "percentile_z", "percentile")
    )
    numeric_metadata <- rbindlist(lapply(numeric_axes, function(axis) {
      merged[
        is.finite(signed_score) & is.finite(get(axis)),
        .(
          n = .N,
          spearman_rho = suppressWarnings(cor(signed_score, get(axis),
                                              method = "spearman",
                                              use = "complete.obs")),
          pearson_r = suppressWarnings(cor(signed_score, get(axis),
                                           method = "pearson",
                                           use = "complete.obs"))
        ),
        by = candidate_state
      ][, metadata_variable := axis]
    }), fill = TRUE)
    if (nrow(numeric_metadata) > 0) {
      numeric_metadata[, abs_spearman := abs(spearman_rho)]
      setorder(numeric_metadata, candidate_state, -abs_spearman)
    }
  }
}

top_metadata <- if (nrow(metadata_enrichment) > 0) {
  metadata_enrichment[
    is.finite(cohen_d),
    .SD[which.max(abs(cohen_d))],
    by = candidate_state
  ][
    ,
    .(
      candidate_state,
      top_metadata_axis = axis,
      top_metadata_value = axis_value,
      top_metadata_n = n_level,
      top_metadata_cohen_d = cohen_d,
      top_metadata_q = q
    )
  ]
} else {
  data.table(candidate_state = character())
}
score_summary <- merge(score_summary, top_metadata, by = "candidate_state",
                       all.x = TRUE, sort = FALSE)
score_summary[is.na(top_metadata_cohen_d), top_metadata_cohen_d := 0]
score_summary[
  ,
  recommendation := fcase(
    genes_present < 4,
    "needs_fst_marker_expansion",
    gene_fraction_present < 0.40 & genes_present < 6,
    "needs_fst_marker_expansion",
    novelty_score < 0.12,
    "mostly_existing_state_alias",
    abs(top_metadata_cohen_d) >= 0.45 & novelty_score >= 0.20,
    "candidate_dominant_state_for_pipeline",
    default = "exploratory_screen_only"
  )
]
recommendation_labels <- c(
  candidate_dominant_state_for_pipeline = "Pipeline candidate",
  exploratory_screen_only = "Exploratory only",
  mostly_existing_state_alias = "Existing-state alias",
  needs_fst_marker_expansion = "Needs FST markers"
)
score_summary[
  ,
  recommendation_label := recommendation_labels[recommendation]
]
score_summary[is.na(recommendation_label), recommendation_label := recommendation]
score_summary[
  ,
  missing_priority_genes := vapply(candidate_state, function(state) {
    paste(candidate_genes[
      candidate_state == state & !available_in_clean_cds,
      head(unique(gene_symbol), 12)
    ], collapse = ";")
  }, character(1))
]
score_summary[, abs_top_metadata_cohen_d := abs(top_metadata_cohen_d)]
setorder(score_summary, -novelty_score, -abs_top_metadata_cohen_d,
         -genes_present)

rdg_gene_links <- data.table()
if (file.exists(atlas_file)) {
  atlas <- fread(atlas_file, showProgress = FALSE)
  if ("gene_symbol" %chin% names(atlas)) {
    atlas[, gene_symbol := toupper(gene_symbol)]
    atlas_summary <- atlas[
      ,
      .(
        atlas_rows = .N,
        top_design_family = design_family[which.max(
          fifelse(is.finite(atlas_browser_review_priority_score),
                  atlas_browser_review_priority_score, -Inf)
        )],
        max_atlas_review_priority = max(atlas_browser_review_priority_score,
                                        na.rm = TRUE),
        top_atlas_evidence_class = atlas_evidence_class[which.max(
          fifelse(is.finite(atlas_browser_review_priority_score),
                  atlas_browser_review_priority_score, -Inf)
        )],
        top_model_agreement_class = model_agreement_class[which.max(
          fifelse(is.finite(atlas_browser_review_priority_score),
                  atlas_browser_review_priority_score, -Inf)
        )],
        top_warnings = warnings[which.max(
          fifelse(is.finite(atlas_browser_review_priority_score),
                  atlas_browser_review_priority_score, -Inf)
        )]
      ),
      by = gene_symbol
    ]
    rdg_gene_links <- merge(
      unique(candidate_genes[, .(candidate_state, gene_symbol, direction,
                                 marker_group, available_in_clean_cds)]),
      atlas_summary,
      by = "gene_symbol",
      all.x = TRUE,
      sort = FALSE
    )
    rdg_gene_links[is.na(atlas_rows), atlas_rows := 0L]
  }
}

summary_metrics <- rbindlist(list(
  data.table(metric = "candidate_modules", value = uniqueN(candidate_genes$candidate_state)),
  data.table(metric = "candidate_module_genes", value = nrow(candidate_genes)),
  data.table(metric = "genes_available_in_current_clean_cds",
             value = uniqueN(candidate_genes[
               available_in_clean_cds == TRUE, gene_symbol
             ])),
  data.table(metric = "pipeline_ready_candidate_states",
             value = score_summary[
               recommendation == "candidate_dominant_state_for_pipeline", .N
             ]),
  data.table(metric = "states_needing_fst_marker_expansion",
             value = score_summary[
               recommendation == "needs_fst_marker_expansion", .N
             ]),
  data.table(metric = "mostly_existing_state_aliases",
             value = score_summary[
               recommendation == "mostly_existing_state_alias", .N
             ])
), fill = TRUE)

fwrite(candidate_genes,
       file.path(output_dir, "dominant_state_auxiliary_module_genes.csv"))
fwrite(scores,
       file.path(output_dir, "dominant_state_auxiliary_module_scores.csv"))
fwrite(dcast(scores, Run + sample ~ candidate_state, value.var = "signed_score"),
       file.path(output_dir, "dominant_state_auxiliary_module_scores_wide.csv"))
fwrite(correlations,
       file.path(output_dir,
                 "dominant_state_auxiliary_existing_state_correlations.csv"))
fwrite(score_summary,
       file.path(output_dir, "dominant_state_auxiliary_module_summary.csv"))
fwrite(metadata_enrichment,
       file.path(output_dir, "dominant_state_auxiliary_metadata_enrichment.csv"))
fwrite(numeric_metadata,
       file.path(output_dir,
                 "dominant_state_auxiliary_numeric_metadata_correlations.csv"))
fwrite(rdg_gene_links,
       file.path(output_dir, "dominant_state_auxiliary_rdg_gene_links.csv"))
fwrite(summary_metrics,
       file.path(output_dir, "dominant_state_auxiliary_summary_metrics.csv"))

if (nrow(score_summary) > 0) {
  p_summary <- ggplot(
    score_summary,
    aes(x = reorder(candidate_state, novelty_score), y = novelty_score,
        fill = recommendation_label)
  ) +
    geom_col(width = 0.72) +
    geom_point(aes(y = gene_fraction_present), shape = 21, size = 2.5,
               fill = "white", colour = "black", stroke = 0.25) +
    coord_flip() +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
    labs(
      x = NULL,
      y = "Novelty vs existing states (bars); marker availability (points)",
      fill = "Recommendation",
      title = "Auxiliary Dominant-State Screen"
    ) +
    scale_fill_manual(
      values = c(
        "Pipeline candidate" = "#4C78A8",
        "Exploratory only" = "#72B7B2",
        "Existing-state alias" = "#B279A2",
        "Needs FST markers" = "#F58518"
      )
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 9)
    )
  save_plot(p_summary, "auxiliary_state_novelty_and_availability", 10.5, 5.2)
}

if (nrow(correlations) > 0) {
  p_cor <- ggplot(
    correlations,
    aes(x = existing_state, y = candidate_state, fill = spearman_rho)
  ) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0, limits = c(-1, 1), name = "Spearman"
    ) +
    labs(x = "Existing dominant state", y = "Auxiliary candidate state",
         title = "Auxiliary State Redundancy Map") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_cor, "auxiliary_state_existing_state_correlations", 9.5, 5.2)
}

if (nrow(metadata_enrichment) > 0) {
  top_plot <- metadata_enrichment[
    is.finite(cohen_d),
    .SD[order(-abs_cohen_d)[seq_len(min(.N, 4))]],
    by = candidate_state
  ]
  top_plot[, label := paste(axis, axis_value, sep = "=")]
  p_meta <- ggplot(
    top_plot,
    aes(x = reorder(label, cohen_d), y = cohen_d, fill = candidate_state)
  ) +
    geom_col(width = 0.72, show.legend = FALSE) +
    facet_wrap(~ candidate_state, scales = "free_y", ncol = 2) +
    coord_flip() +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.25) +
    labs(x = NULL, y = "Cohen d for auxiliary score",
         title = "Top Metadata Associations For Auxiliary States") +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.y = element_blank())
  save_plot(p_meta, "auxiliary_state_top_metadata_effects", 10.5, 8.5)
}

message("Saved auxiliary dominant-state screen to: ", output_dir)
