#!/usr/bin/env Rscript
#
# Additive, standalone follow-up analysis implementing the two sub-score
# splits proposed in md/dominant_state_gene_set_revision_plan_2026-09-25.md:
#
#   ISR / ER stress (ATF4 axis) ->
#     isr_atf4_direct   (ATF4, ATF3, ASNS -- direct ATF4 transcriptional
#                        targets, plausibly dominated by ATF4's own
#                        uORF-translation-efficiency switch)
#     isr_terminal_upr  (DDIT3, PPP1R15A, DDIT4, HERPUD1, HSPA5, XBP1,
#                        CEBPB -- the terminal/growth-arrest UPR genes)
#
#   Senescence / SASP ->
#     senescence_acute_sasp (SERPINE1, IGFBP7, TGFBI, TIMP1, CCND1, IL6,
#                            CDKN1A, CXCL8 -- the cohesive acute SASP core)
#     senescence_commitment (CDKN2A up, LMNB1 down -- the slow,
#                            chromatin-level commitment/cell-cycle-arrest
#                            axis, which showed near-zero/negative cohesion
#                            with the acute SASP core)
#
# This does not modify dominant_state_module_definitions.R or the main
# clean_cds scoring pipeline -- it recomputes the same score_module() logic
# (gene_z <- log2(expr+1) recentred on the Baseline mean) directly from the
# already-produced per-sample gene expression matrix, using the corrected
# (B2M/UBC removed) Baseline gene list, so it reflects the current module
# definitions without requiring this to be wired into the pipeline's
# dependency graph yet.

suppressPackageStartupMessages(library(data.table))

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (basename(here) == "dominant_cell_states" &&
        dir.exists(file.path(here, "scripts"))) {
      return(here)
    }
    candidate <- file.path(here, "dominant_cell_states")
    if (dir.exists(file.path(candidate, "scripts"))) {
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
results_dir <- Sys.getenv("DOMINANT_CELL_STATES_RESULTS_DIR",
                          unset = file.path(analysis_dir, "results"))
result_file <- function(...) file.path(results_dir, ...)

source(file.path(analysis_dir, "scripts", "dominant_state_module_definitions.R"))
modules <- load_dominant_state_modules()
baseline_genes <- modules[dominant_state == "Baseline (control)", up_genes[[1]]]

expr_file <- result_file("human_dominant_cell_states_clean_cds_gene_expression.csv")
if (!file.exists(expr_file)) {
  stop("Missing ", expr_file, " -- run dominant_cell_states_fst_clean_cds.R first.")
}
expr <- fread(expr_file)
gene_symbols <- expr$gene_symbol
expr_mat <- as.matrix(expr[, -1, with = FALSE])
rownames(expr_mat) <- gene_symbols
storage.mode(expr_mat) <- "numeric"

expr_norm <- log2(expr_mat + 1)
baseline_present <- intersect(baseline_genes, rownames(expr_norm))
if (length(baseline_present) < 5) {
  stop("Too few baseline genes present: ", length(baseline_present))
}
baseline <- colMeans(expr_norm[baseline_present, , drop = FALSE], na.rm = TRUE)
gene_z <- sweep(expr_norm, 2, baseline, "-")

score_module <- function(gene_z, up_genes, down_genes = character()) {
  up_present <- intersect(up_genes, rownames(gene_z))
  down_present <- intersect(down_genes, rownames(gene_z))
  up_score <- if (length(up_present) == 0) rep(NA_real_, ncol(gene_z)) else
    colMeans(gene_z[up_present, , drop = FALSE], na.rm = TRUE)
  down_score <- if (length(down_genes) == 0) rep(0, ncol(gene_z)) else if (length(down_present) == 0)
    rep(NA_real_, ncol(gene_z)) else colMeans(gene_z[down_present, , drop = FALSE], na.rm = TRUE)
  up_score - down_score
}

isr_atf4_direct_genes <- c("ATF4", "ATF3", "ASNS")
isr_terminal_upr_genes <- c("DDIT3", "PPP1R15A", "DDIT4", "HERPUD1", "HSPA5", "XBP1", "CEBPB")
senescence_acute_sasp_genes <- c("SERPINE1", "IGFBP7", "TGFBI", "TIMP1", "CCND1", "IL6", "CDKN1A", "CXCL8")
senescence_commitment_up_genes <- c("CDKN2A")
senescence_commitment_down_genes <- c("LMNB1")

out <- data.table(
  Run = colnames(gene_z),
  isr_atf4_direct = score_module(gene_z, isr_atf4_direct_genes),
  isr_terminal_upr = score_module(gene_z, isr_terminal_upr_genes),
  senescence_acute_sasp = score_module(gene_z, senescence_acute_sasp_genes),
  senescence_commitment = score_module(gene_z, senescence_commitment_up_genes,
                                       senescence_commitment_down_genes)
)

out_file <- result_file("dominant_state_isr_senescence_subscores.csv")
fwrite(out, out_file)
message("Saved: ", out_file)

cor_isr <- cor(out$isr_atf4_direct, out$isr_terminal_upr, use = "pairwise.complete.obs")
cor_senesc <- cor(out$senescence_acute_sasp, out$senescence_commitment, use = "pairwise.complete.obs")
message(sprintf("ISR sub-score correlation (atf4_direct vs terminal_upr): r = %.3f", cor_isr))
message(sprintf("Senescence sub-score correlation (acute_sasp vs commitment): r = %.3f", cor_senesc))

summary_dt <- data.table(
  sub_score = c("isr_atf4_direct", "isr_terminal_upr",
               "senescence_acute_sasp", "senescence_commitment"),
  mean = c(mean(out$isr_atf4_direct, na.rm = TRUE),
          mean(out$isr_terminal_upr, na.rm = TRUE),
          mean(out$senescence_acute_sasp, na.rm = TRUE),
          mean(out$senescence_commitment, na.rm = TRUE)),
  sd = c(sd(out$isr_atf4_direct, na.rm = TRUE),
        sd(out$isr_terminal_upr, na.rm = TRUE),
        sd(out$senescence_acute_sasp, na.rm = TRUE),
        sd(out$senescence_commitment, na.rm = TRUE)),
  cor_with_partner = c(cor_isr, cor_isr, cor_senesc, cor_senesc)
)
fwrite(summary_dt, result_file("dominant_state_isr_senescence_subscores_summary.csv"))
print(summary_dt)
