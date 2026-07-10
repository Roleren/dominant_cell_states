dcs_legacy_result_link_paths <- function(root = dcs_project_root()) {
  generated_names <- c(
    "dominant_pipeline_last_run_summary.csv",
    "htmlwidget_libs",
    "postviral_fatigue_outputs",
    "readlength_inputs",
    "required_fst_page_basenames.txt",
    "required_fst_pages.csv",
    "required_fst_pages.txt",
    "required_fst_pages_by_gene.csv",
    "required_fst_pages_summary.txt",
    "required_fst_selected_transcripts.csv"
  )

  results_dir <- dcs_results_dir(root, create = FALSE)
  if (dir.exists(results_dir)) {
    result_names <- list.files(
      results_dir,
      full.names = FALSE,
      recursive = FALSE,
      all.files = FALSE,
      no.. = TRUE
    )
    generated_names <- unique(c(generated_names, result_names))
  }

  generated_names <- generated_names[
    grepl("^dominant_", generated_names) |
      grepl("^human_dominant_", generated_names) |
      generated_names %in% c(
        "htmlwidget_libs",
        "postviral_fatigue_outputs",
        "readlength_inputs",
        "required_fst_page_basenames.txt",
        "required_fst_pages.csv",
        "required_fst_pages.txt",
        "required_fst_pages_by_gene.csv",
        "required_fst_pages_summary.txt",
        "required_fst_selected_transcripts.csv"
      )
  ]

  manual_generated <- file.path(
    "manual_translons",
    c(
      "manual_translon_candidates.csv",
      "manual_translon_requested_gene_uoorf_candidates.csv",
      "manual_translons_manifest.csv",
      "manual_translons_ranges.rds"
    )
  )

  root_entries <- list.files(
    root,
    full.names = TRUE,
    recursive = FALSE,
    all.files = FALSE,
    no.. = TRUE
  )
  root_readlinks <- Sys.readlink(root_entries)
  root_links <- root_entries[!is.na(root_readlinks) & nzchar(root_readlinks)]

  manual_paths <- file.path(root, manual_generated)
  manual_readlinks <- Sys.readlink(manual_paths)
  manual_links <- manual_paths[
    !is.na(manual_readlinks) & nzchar(manual_readlinks)
  ]

  unique(c(file.path(root, generated_names), root_links, manual_paths, manual_links))
}

dcs_remove_legacy_result_links <- function(root = dcs_project_root()) {
  paths <- dcs_legacy_result_link_paths(root)
  readlinks <- Sys.readlink(paths)
  links <- paths[!is.na(readlinks) & nzchar(readlinks)]
  if (length(links)) {
    unlink(links)
  }
  invisible(links)
}
