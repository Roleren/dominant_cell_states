# Shared helpers for pipeline step scripts.
#
# Step scripts source this file with
#   source(file.path(Sys.getenv("DOMINANT_CELL_STATES_DIR", "."), "scripts",
#                    "dominant_script_common.R"), local = TRUE)
# The runner sets DOMINANT_CELL_STATES_DIR for every step; standalone runs
# start from the project root.

# Project root: DOMINANT_CELL_STATES_DIR when set, otherwise the nearest
# directory (or dominant_cell_states child) that holds scripts/ and the
# scientific note. Unlike the per-script copies this replaced, the checkout
# directory may have any name (e.g. a git worktree).
find_analysis_dir <- function(start = getwd()) {
  is_root <- function(path) {
    dir.exists(file.path(path, "scripts")) &&
      file.exists(file.path(path, "md", "dominant_cell_state_scientific_note.md"))
  }
  env_root <- Sys.getenv("DOMINANT_CELL_STATES_DIR", unset = "")
  if (nzchar(env_root)) {
    env_root <- normalizePath(env_root, mustWork = FALSE)
    if (is_root(env_root)) return(env_root)
    stop("DOMINANT_CELL_STATES_DIR is not a dominant_cell_states root: ",
         env_root, call. = FALSE)
  }
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (is_root(here)) return(here)
    candidate <- file.path(here, "dominant_cell_states")
    if (is_root(candidate)) return(normalizePath(candidate, mustWork = TRUE))
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find dominant_cell_states analysis directory from: ",
       start, call. = FALSE)
}
