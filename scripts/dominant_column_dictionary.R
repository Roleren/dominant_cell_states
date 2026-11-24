# Shared helper: write a small, paired column dictionary next to a derived
# output table.
#
# Implements gap #4 from md/dominant_state_pipeline_lessons_2026-09-28.md:
# derived-table columns were not self-documenting, and understanding a
# non-obvious one (e.g. `branch_class` in
# dominant_rdg_context_interaction_branch_context_effects.csv) required
# grepping the source script directly -- this cost real time and
# contributed to an initial wrong conclusion (see
# papers/to_be_explored_more_later/VEGFA_uorf_viral_infection/README.md).
#
# Usage from any script that writes a results CSV:
#
#   source(file.path(analysis_dir, "scripts", "dominant_column_dictionary.R"))
#   write_column_dictionary(
#     output_csv_path,
#     descriptions = c(
#       gene_symbol = "HGNC gene symbol.",
#       my_score    = "Positive = more of X; negative = less. Units: log2."
#     ),
#     title = "My output table"
#   )
#
# This writes `<output_csv_path minus .csv>_columns.md` alongside the table.
# `descriptions` only needs entries for columns that are NOT self-explanatory
# from their name alone (obvious ones like `Run`, `gene_symbol`, `n_samples`
# can be omitted) -- the dictionary always lists every actual column in the
# table, flagging any without a supplied description as "(undocumented)"
# rather than silently omitting them, so the gap itself stays visible.

write_column_dictionary <- function(csv_path, descriptions = character(),
                                    title = NULL, notes = character()) {
  if (!file.exists(csv_path)) {
    warning("write_column_dictionary: csv does not exist yet, skipping: ", csv_path)
    return(invisible(NULL))
  }
  header <- tryCatch(
    names(data.table::fread(csv_path, nrows = 0)),
    error = function(e) NULL
  )
  if (is.null(header)) {
    warning("write_column_dictionary: could not read header of ", csv_path)
    return(invisible(NULL))
  }

  dict_path <- sub("\\.csv$", "_columns.md", csv_path)
  title <- if (is.null(title)) basename(csv_path) else title
  title <- sub("[.]+$", "", trimws(title))

  lines <- c(
    paste0("# Column dictionary: ", basename(csv_path)),
    "",
    paste0(title, ". Generated alongside the table it describes -- if the ",
          "table's own columns change, re-running the producing step ",
          "regenerates this file too, so drift between the two is visible ",
          "(an undocumented new column shows up as \"(undocumented)\" below ",
          "rather than silently missing)."),
    ""
  )
  if (length(notes) > 0) {
    lines <- c(lines, notes, "")
  }
  lines <- c(lines, "| column | meaning |", "|---|---|")
  for (col in header) {
    desc <- if (col %in% names(descriptions)) descriptions[[col]] else "(undocumented)"
    desc <- gsub("\\|", "\\\\|", desc)
    lines <- c(lines, sprintf("| `%s` | %s |", col, desc))
  }

  n_undocumented <- sum(!header %in% names(descriptions))
  if (n_undocumented > 0) {
    lines <- c(lines, "", sprintf(
      "%d of %d columns are undocumented above -- worth filling in if this table gets used again.",
      n_undocumented, length(header)
    ))
  }

  writeLines(lines, dict_path)
  message("Saved: ", dict_path,
          if (n_undocumented > 0) sprintf(" (%d undocumented columns)", n_undocumented) else "")
  invisible(dict_path)
}
