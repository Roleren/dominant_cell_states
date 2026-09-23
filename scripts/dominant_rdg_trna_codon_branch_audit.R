#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Biostrings)
  library(data.table)
  library(GenomicRanges)
  library(Rsamtools)
})

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(here, "DESCRIPTION")) &&
        dir.exists(file.path(here, "results"))) return(here)
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find dominant_cell_states analysis root.", call. = FALSE)
}

analysis_dir <- find_analysis_dir()
results_dir <- file.path(analysis_dir, "results")
output_dir <- file.path(results_dir, "dominant_rdg_trna_codon_branch_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

adapter <- fread(file.path(
  results_dir, "dominant_rdg_grouped_branch_usage",
  "dominant_rdg_grouped_branch_usage_adapter.csv"
))
annotation <- fread(file.path(
  results_dir, "dominant_rdg_outputs", "rdg_feature_annotation.csv"
))
transcript_fasta <- Sys.getenv(
  "DOMINANT_TRANSCRIPT_FASTA",
  "/media/roler/S/data/Bio_data/references/Homo_sapiens_GRCh38_101/Ribo-code-transcripts/transcripts_sequence.fa"
)
if (!file.exists(transcript_fasta)) {
  stop("Transcript FASTA not found: ", transcript_fasta, call. = FALSE)
}

feature_contrast <- function(study_id, case_value) {
  x <- adapter[
    study == study_id & CONDITION %chin% c("WT", case_value) &
      feature_type %chin% c("leader_uorf", "clean_CDS")
  ]
  # FRACTION is a sequencing-lane label in the tRNA-Glu study. Summing these
  # rows reconstructs each two-replicate arm without treating lanes as n.
  x <- x[, .(
    counts = sum(raw_counts),
    n_runs = sum(n_runs_merged),
    lane_labels = paste(sort(unique(FRACTION)), collapse = ";")
  ), by = .(
    gene_symbol, tx_id, feature_id, feature_type, CONDITION,
    feature_bases, start_codon = fifelse(
      feature_type == "clean_CDS", NA_character_, "uORF"
    )
  )]
  cds <- x[feature_type == "clean_CDS", .(
    gene_symbol, tx_id, CONDITION, clean_cds_counts = counts
  )]
  uorf <- merge(
    x[feature_type == "leader_uorf"], cds,
    by = c("gene_symbol", "tx_id", "CONDITION"), all = FALSE
  )
  uorf[, uorf_allocation := counts / pmax(counts + clean_cds_counts, 1L)]
  wide <- dcast(
    uorf,
    gene_symbol + tx_id + feature_id + feature_bases ~ CONDITION,
    value.var = c("counts", "clean_cds_counts", "uorf_allocation"),
    fill = 0
  )
  required <- paste0(
    c("counts_", "clean_cds_counts_", "uorf_allocation_"),
    rep(c(case_value, "WT"), each = 3L)
  )
  if (!all(required %in% names(wide))) {
    stop("Failed to recover both arms for ", study_id, call. = FALSE)
  }
  setnames(
    wide,
    c(
      paste0("counts_", case_value), "counts_WT",
      paste0("clean_cds_counts_", case_value), "clean_cds_counts_WT",
      paste0("uorf_allocation_", case_value), "uorf_allocation_WT"
    ),
    c(
      "case_uorf_counts", "control_uorf_counts",
      "case_clean_cds_counts", "control_clean_cds_counts",
      "case_uorf_allocation", "control_uorf_allocation"
    )
  )
  wide[, uorf_allocation_delta :=
         case_uorf_allocation - control_uorf_allocation]
  wide[, `:=`(study = study_id, case_condition = case_value)]
  wide
}

glu <- feature_contrast("PRJNA310107-homo_sapiens", "glu")
arg <- feature_contrast("PRJNA310096-homo_sapiens", "Arg")

feature_annotation <- unique(annotation[
  feature_type == "leader_uorf" & tx_start >= 1L & tx_end >= tx_start,
  .(
    gene_symbol, tx_id, feature_id, tx_start, tx_end, feature_bases,
    annotated_start_codon = start_codon,
    start_codon_class, kozak_strength
  )
])
cds_annotation <- unique(annotation[
  feature_type == "clean_CDS" & tx_start >= 1L & tx_end >= tx_start,
  .(gene_symbol, tx_id, cds_start = tx_start, cds_end = tx_end)
])

mark_overlap_free <- function(x) {
  x[, uorf_overlap_free := {
    ranges <- IRanges(tx_start, tx_end)
    countOverlaps(ranges, ranges) == 1L
  }, by = tx_id]
  x
}
feature_annotation <- mark_overlap_free(feature_annotation)

sequence_metrics <- function(dt, start_col, end_col, prefix) {
  ranges <- GRanges(
    dt$tx_id,
    IRanges(as.integer(dt[[start_col]]), as.integer(dt[[end_col]]))
  )
  sequence <- as.character(scanFa(FaFile(transcript_fasta), param = ranges))
  metric <- rbindlist(lapply(sequence, function(value) {
    n <- nchar(value)
    codon_start <- seq.int(1L, n - 2L, by = 3L)
    codons <- substring(value, codon_start, codon_start + 2L)
    data.table(
      sequence = value,
      codon_count = length(codons),
      glu_codon_count = sum(codons %chin% c("GAA", "GAG")),
      gag_count = sum(codons == "GAG"),
      gaa_count = sum(codons == "GAA"),
      arg_cgg_count = sum(codons == "CGG"),
      gc_fraction = sum(strsplit(value, "", fixed = TRUE)[[1L]] %chin%
                          c("G", "C")) / n
    )
  }))
  metric[, glu_codon_fraction := glu_codon_count / pmax(codon_count, 1L)]
  metric[, arg_cgg_fraction := arg_cgg_count / pmax(codon_count, 1L)]
  setnames(metric, names(metric), paste0(prefix, names(metric)))
  cbind(copy(dt), metric)
}

feature_annotation <- sequence_metrics(
  feature_annotation, "tx_start", "tx_end", "uorf_"
)
cds_annotation <- sequence_metrics(
  cds_annotation, "cds_start", "cds_end", "cds_"
)
sequence_table <- merge(
  feature_annotation,
  cds_annotation[, .(
    gene_symbol, tx_id, cds_codon_count, cds_glu_codon_count,
    cds_glu_codon_fraction, cds_arg_cgg_count, cds_arg_cgg_fraction,
    cds_gc_fraction
  )],
  by = c("gene_symbol", "tx_id"), all = FALSE
)
sequence_table[, relative_glu_codon_fraction :=
                 uorf_glu_codon_fraction - cds_glu_codon_fraction]
sequence_table[, relative_arg_cgg_fraction :=
                 uorf_arg_cgg_fraction - cds_arg_cgg_fraction]

joint <- merge(
  glu[, .(
    gene_symbol, tx_id, feature_id,
    glu_case_uorf_counts = case_uorf_counts,
    glu_control_uorf_counts = control_uorf_counts,
    glu_case_clean_cds_counts = case_clean_cds_counts,
    glu_control_clean_cds_counts = control_clean_cds_counts,
    glu_uorf_allocation_delta = uorf_allocation_delta,
    glu_control_uorf_allocation = control_uorf_allocation
  )],
  arg[, .(
    gene_symbol, tx_id, feature_id,
    arg_case_uorf_counts = case_uorf_counts,
    arg_control_uorf_counts = control_uorf_counts,
    arg_case_clean_cds_counts = case_clean_cds_counts,
    arg_control_clean_cds_counts = control_clean_cds_counts,
    arg_uorf_allocation_delta = uorf_allocation_delta,
    arg_control_uorf_allocation = control_uorf_allocation
  )],
  by = c("gene_symbol", "tx_id", "feature_id"), all = FALSE
)
joint <- merge(joint, sequence_table,
               by = c("gene_symbol", "tx_id", "feature_id"), all = FALSE)
joint[, glu_specific_allocation_delta :=
        glu_uorf_allocation_delta - arg_uorf_allocation_delta]
joint[, count_gate :=
        pmin(glu_case_uorf_counts, glu_control_uorf_counts,
             arg_case_uorf_counts, arg_control_uorf_counts) >= 20L &
        pmin(glu_case_clean_cds_counts, glu_control_clean_cds_counts,
             arg_case_clean_cds_counts, arg_control_clean_cds_counts) >= 100L]

analysis <- joint[
  count_gate == TRUE & uorf_overlap_free == TRUE &
    uorf_codon_count >= 5L & is.finite(relative_glu_codon_fraction)
]

rank_test <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 10L) return(list(rho = NA_real_, p = NA_real_, n = sum(keep)))
  test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman",
                                    exact = FALSE))
  list(rho = unname(test$estimate), p = test$p.value, n = sum(keep))
}

tests <- rbindlist(list(
  cbind(data.table(
    test = "Glu OE delta vs relative Glu-codon fraction"
  ), as.data.table(rank_test(
    analysis$relative_glu_codon_fraction,
    analysis$glu_uorf_allocation_delta
  ))),
  cbind(data.table(
    test = "Arg OE delta vs relative Glu-codon fraction (negative control)"
  ), as.data.table(rank_test(
    analysis$relative_glu_codon_fraction,
    analysis$arg_uorf_allocation_delta
  ))),
  cbind(data.table(
    test = "Glu-specific delta vs relative Glu-codon fraction"
  ), as.data.table(rank_test(
    analysis$relative_glu_codon_fraction,
    analysis$glu_specific_allocation_delta
  ))),
  cbind(data.table(
    test = "Arg OE delta vs relative CGG-codon fraction"
  ), as.data.table(rank_test(
    analysis$relative_arg_cgg_fraction,
    analysis$arg_uorf_allocation_delta
  )))
), fill = TRUE)

rank_sensitivity <- rbindlist(lapply(c(5L, 10L, 20L), function(minimum) {
  x <- joint[
    pmin(glu_case_uorf_counts, glu_control_uorf_counts,
         arg_case_uorf_counts, arg_control_uorf_counts) >= minimum &
      pmin(glu_case_clean_cds_counts, glu_control_clean_cds_counts,
           arg_case_clean_cds_counts, arg_control_clean_cds_counts) >= 100L &
      uorf_overlap_free == TRUE & uorf_codon_count >= 5L
  ]
  x[, glu_delta_rank := frank(-glu_uorf_allocation_delta,
                              ties.method = "min")]
  x[, glu_specific_rank := frank(-glu_specific_allocation_delta,
                                 ties.method = "min")]
  eif <- x[gene_symbol == "EIF4G2" & feature_id == "uORF1"]
  specificity_median <- median(x$glu_specific_allocation_delta, na.rm = TRUE)
  specificity_raw_mad <- mad(
    x$glu_specific_allocation_delta, center = specificity_median,
    constant = 1, na.rm = TRUE
  )
  data.table(
    minimum_uorf_counts_each_arm = minimum,
    n_eligible_features = nrow(x),
    eif4g2_glu_delta_rank = if (nrow(eif)) eif$glu_delta_rank else NA_integer_,
    eif4g2_glu_specific_rank = if (nrow(eif)) {
      eif$glu_specific_rank
    } else {
      NA_integer_
    },
    eif4g2_glu_delta = if (nrow(eif)) {
      eif$glu_uorf_allocation_delta
    } else {
      NA_real_
    },
    eif4g2_arg_delta = if (nrow(eif)) {
      eif$arg_uorf_allocation_delta
    } else {
      NA_real_
    },
    eif4g2_glu_specific_delta = if (nrow(eif)) {
      eif$glu_specific_allocation_delta
    } else {
      NA_real_
    },
    background_glu_specific_median = specificity_median,
    background_glu_specific_raw_mad = specificity_raw_mad,
    eif4g2_glu_specific_robust_z = if (
      nrow(eif) && is.finite(specificity_raw_mad) &&
        specificity_raw_mad > 0
    ) {
      (eif$glu_specific_allocation_delta - specificity_median) /
        specificity_raw_mad
    } else {
      NA_real_
    },
    empirical_top_rank_fraction = if (nrow(eif) && nrow(x)) {
      eif$glu_specific_rank / nrow(x)
    } else {
      NA_real_
    }
  )
}))

analysis[, glu_specific_rank := frank(-glu_specific_allocation_delta,
                                      ties.method = "min")]
setorder(analysis, glu_specific_rank)

fwrite(joint, file.path(output_dir, "trna_codon_feature_audit.csv"))
fwrite(analysis, file.path(output_dir, "trna_codon_analysis_set.csv"))
fwrite(tests, file.path(output_dir, "trna_codon_rank_tests.csv"))
fwrite(rank_sensitivity,
       file.path(output_dir, "eif4g2_specificity_rank_sensitivity.csv"))

png(file.path(output_dir, "trna_glu_codon_branch_audit.png"),
    width = 1800, height = 900, res = 180)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
plot(
  analysis$relative_glu_codon_fraction,
  analysis$glu_uorf_allocation_delta,
  pch = 16, col = grDevices::adjustcolor("#0072B2", 0.45),
  xlab = "uORF minus CDS Glu-codon fraction",
  ylab = "uORF allocation delta (tRNA-Glu OE - WT)",
  main = "tRNA-Glu perturbation"
)
abline(h = 0, v = 0, lty = 2, col = "grey60")
eif <- analysis[gene_symbol == "EIF4G2"]
if (nrow(eif)) {
  points(eif$relative_glu_codon_fraction, eif$glu_uorf_allocation_delta,
         pch = 21, bg = "#D55E00", col = "white", cex = 1.5)
  text(eif$relative_glu_codon_fraction, eif$glu_uorf_allocation_delta,
       labels = eif$feature_id, pos = 4, cex = 0.8)
}
plot(
  analysis$relative_glu_codon_fraction,
  analysis$arg_uorf_allocation_delta,
  pch = 16, col = grDevices::adjustcolor("#666666", 0.45),
  xlab = "uORF minus CDS Glu-codon fraction",
  ylab = "uORF allocation delta (tRNA-Arg OE - WT)",
  main = "Matched tRNA-Arg control"
)
abline(h = 0, v = 0, lty = 2, col = "grey60")
if (nrow(eif)) {
  points(eif$relative_glu_codon_fraction, eif$arg_uorf_allocation_delta,
         pch = 21, bg = "#D55E00", col = "white", cex = 1.5)
  text(eif$relative_glu_codon_fraction, eif$arg_uorf_allocation_delta,
       labels = eif$feature_id, pos = 4, cex = 0.8)
}
dev.off()

png(file.path(output_dir, "trna_glu_specificity_outlier.png"),
    width = 1800, height = 900, res = 180)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
limits <- range(c(
  analysis$arg_uorf_allocation_delta,
  analysis$glu_uorf_allocation_delta
), finite = TRUE)
plot(
  analysis$arg_uorf_allocation_delta,
  analysis$glu_uorf_allocation_delta,
  pch = 16, col = grDevices::adjustcolor("#666666", 0.55),
  xlim = limits, ylim = limits,
  xlab = "uORF allocation delta: tRNA-Arg OE",
  ylab = "uORF allocation delta: tRNA-Glu OE",
  main = "Matched tRNA perturbations"
)
abline(a = 0, b = 1, lty = 2, col = "grey55")
eif <- analysis[gene_symbol == "EIF4G2" & feature_id == "uORF1"]
if (nrow(eif)) {
  points(eif$arg_uorf_allocation_delta, eif$glu_uorf_allocation_delta,
         pch = 21, bg = "#D55E00", col = "white", cex = 1.7)
  text(eif$arg_uorf_allocation_delta, eif$glu_uorf_allocation_delta,
       labels = "EIF4G2 uORF1", pos = 4, cex = 0.85)
}
plot(
  analysis$glu_specific_rank,
  analysis$glu_specific_allocation_delta,
  pch = 16, col = grDevices::adjustcolor("#0072B2", 0.65),
  xlab = "Rank among strict count-supported leader uORFs",
  ylab = "tRNA-Glu delta minus tRNA-Arg delta",
  main = "Perturbation-specific branch response"
)
abline(h = 0, lty = 2, col = "grey55")
if (nrow(eif)) {
  points(eif$glu_specific_rank, eif$glu_specific_allocation_delta,
         pch = 21, bg = "#D55E00", col = "white", cex = 1.7)
  text(eif$glu_specific_rank, eif$glu_specific_allocation_delta,
       labels = "EIF4G2 uORF1", pos = 4, cex = 0.85)
}
dev.off()

print(tests)
print(rank_sensitivity)
message("Analysis features: ", nrow(analysis))
message("Saved tRNA codon/branch audit: ", output_dir)
