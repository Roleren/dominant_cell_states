#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

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
output_dir <- file.path(
  analysis_dir, "results", "dominant_rdg_eif4g2_goodarzi_rna_audit"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

samples <- CJ(
  condition = c("ctrl", "glu"),
  replicate = 1:2,
  assay = c("RPF", "TT"),
  unique = TRUE
)
samples[, sample := sprintf("MDA-%s-%d-%s", condition, replicate, assay)]
samples[, file_name := sprintf(
  "GSE77347_%s_isoforms.fpkm_tracking.gz", sample
)]
samples[, url := sprintf(
  paste0(
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE77347&file=%s",
    "&format=file"
  ),
  file_name
)]

read_sample <- function(row) {
  destination <- file.path(tempdir(), row$file_name)
  if (!file.exists(destination) || file.info(destination)$size <= 0) {
    download.file(row$url, destination, mode = "wb", quiet = TRUE,
                  method = "libcurl")
  }
  x <- fread(destination)
  x <- x[gene_short_name == "EIF4G2"]
  if (!nrow(x)) stop("EIF4G2 missing from ", row$file_name, call. = FALSE)
  x[, `:=`(
    condition = row$condition,
    replicate = row$replicate,
    assay = row$assay,
    source_file = row$file_name
  )]
  x
}

isoforms <- rbindlist(lapply(seq_len(nrow(samples)), function(i) {
  read_sample(samples[i])
}), fill = TRUE)

paired <- dcast(
  isoforms,
  tracking_id + gene_short_name + condition + replicate ~ assay,
  value.var = "FPKM"
)
if (!all(c("RPF", "TT") %in% names(paired))) {
  stop("Failed to recover paired RPF and total-transcript FPKM.",
       call. = FALSE)
}
paired[, translation_efficiency_proxy := RPF / pmax(TT, 1e-8)]
paired[, primary_eif4g2_isoform := tracking_id == "NM_001418"]

summary_long <- paired[, .(
  mean_rpf_fpkm = mean(RPF),
  mean_total_transcript_fpkm = mean(TT),
  mean_translation_efficiency_proxy = mean(translation_efficiency_proxy)
), by = .(tracking_id, condition)]
summary_long[, total_transcript_isoform_fraction :=
  mean_total_transcript_fpkm / sum(mean_total_transcript_fpkm),
  by = condition
]
summary_long[, rpf_isoform_fraction :=
  mean_rpf_fpkm / sum(mean_rpf_fpkm),
  by = condition
]
summary <- dcast(
  summary_long,
  tracking_id ~ condition,
  value.var = c(
    "mean_rpf_fpkm", "mean_total_transcript_fpkm",
    "mean_translation_efficiency_proxy",
    "total_transcript_isoform_fraction", "rpf_isoform_fraction"
  )
)
summary[, rpf_log2_glu_vs_control := log2(
  mean_rpf_fpkm_glu / mean_rpf_fpkm_ctrl
)]
summary[, total_transcript_log2_glu_vs_control := log2(
  mean_total_transcript_fpkm_glu / mean_total_transcript_fpkm_ctrl
)]
summary[, translation_efficiency_log2_glu_vs_control := log2(
  mean_translation_efficiency_proxy_glu /
    mean_translation_efficiency_proxy_ctrl
)]
summary[, total_transcript_isoform_fraction_delta :=
  total_transcript_isoform_fraction_glu -
    total_transcript_isoform_fraction_ctrl]
summary[, rpf_isoform_fraction_delta :=
  rpf_isoform_fraction_glu - rpf_isoform_fraction_ctrl]
summary[, total_transcript_stability_gate :=
          abs(total_transcript_log2_glu_vs_control) <= 0.25]
summary[, isoform_composition_stability_gate :=
          abs(total_transcript_isoform_fraction_delta) <= 0.05]
summary[, primary_eif4g2_isoform := tracking_id == "NM_001418"]

fwrite(samples, file.path(output_dir, "geo_processed_file_manifest.csv"))
fwrite(isoforms, file.path(output_dir, "eif4g2_isoform_fpkm.csv"))
fwrite(paired, file.path(output_dir, "eif4g2_paired_rpf_total_transcript.csv"))
fwrite(summary, file.path(output_dir, "eif4g2_isoform_rna_audit_summary.csv"))

print(summary[primary_eif4g2_isoform == TRUE])
message("Saved Goodarzi EIF4G2 RNA audit: ", output_dir)
