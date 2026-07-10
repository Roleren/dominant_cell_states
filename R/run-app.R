#' Run the RDG atlas browser Shiny app
#'
#' Launches the dominant-cell-state ribosome decision graph atlas browser using
#' the local analysis checkout as the data root. Generated app inputs are read
#' from `results/`.
#'
#' @param root Dominant-cell-states project root. Defaults to
#'   [dcs_project_root()].
#' @param host Host interface passed to [shiny::runApp()].
#' @param port Port passed to [shiny::runApp()]. `NULL` lets Shiny choose.
#' @param launch.browser Whether Shiny should open a browser.
#' @param ... Additional arguments passed to [shiny::runApp()].
#'
#' @return The value returned by [shiny::runApp()].
#' @export
run_dcs_app <- function(root = dcs_project_root(),
                        host = Sys.getenv("RDG_ATLAS_BROWSER_HOST",
                                          "127.0.0.1"),
                        port = {
                          value <- Sys.getenv("RDG_ATLAS_BROWSER_PORT", "")
                          if (nzchar(value)) as.integer(value) else NULL
                        },
                        launch.browser = interactive(),
                        ...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to run the RDG atlas browser.",
         call. = FALSE)
  }

  root <- normalizePath(root, mustWork = TRUE)
  app_dir <- system.file(
    "shiny",
    "rdg_atlas_browser",
    package = "dominantCellStates"
  )
  if (!nzchar(app_dir)) {
    app_dir <- file.path(root, "inst", "shiny", "rdg_atlas_browser")
  }
  if (!file.exists(file.path(app_dir, "app.R"))) {
    stop("Missing RDG atlas browser app.R under: ", app_dir, call. = FALSE)
  }

  old_root <- Sys.getenv("DOMINANT_CELL_STATES_DIR", unset = NA_character_)
  old_results <- Sys.getenv(
    "DOMINANT_CELL_STATES_RESULTS_DIR",
    unset = NA_character_
  )
  Sys.setenv(
    DOMINANT_CELL_STATES_DIR = root,
    DOMINANT_CELL_STATES_RESULTS_DIR = dcs_results_dir(root, create = FALSE)
  )
  on.exit({
    if (is.na(old_root)) {
      Sys.unsetenv("DOMINANT_CELL_STATES_DIR")
    } else {
      Sys.setenv(DOMINANT_CELL_STATES_DIR = old_root)
    }
    if (is.na(old_results)) {
      Sys.unsetenv("DOMINANT_CELL_STATES_RESULTS_DIR")
    } else {
      Sys.setenv(DOMINANT_CELL_STATES_RESULTS_DIR = old_results)
    }
  }, add = TRUE)

  args <- list(
    appDir = app_dir,
    host = host,
    launch.browser = launch.browser,
    ...
  )
  if (!is.null(port)) {
    args$port <- port
  }
  do.call(shiny::runApp, args)
}
