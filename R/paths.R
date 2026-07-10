dcs_is_project_root <- function(path) {
  file.exists(file.path(path, "md", "dominant_cell_state_scientific_note.md")) &&
    dir.exists(file.path(path, "scripts")) &&
    dir.exists(file.path(path, "inst", "shiny", "rdg_atlas_browser"))
}

dcs_project_root <- function(start = getwd()) {
  env_root <- Sys.getenv("DOMINANT_CELL_STATES_DIR", unset = "")
  if (nzchar(env_root)) {
    env_root <- normalizePath(env_root, mustWork = FALSE)
    if (dcs_is_project_root(env_root)) {
      return(env_root)
    }
    stop("DOMINANT_CELL_STATES_DIR is not a dominant_cell_states root: ",
         env_root, call. = FALSE)
  }

  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (dcs_is_project_root(here)) {
      return(here)
    }

    candidate <- file.path(here, "dominant_cell_states")
    if (dcs_is_project_root(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }

    parent <- dirname(here)
    if (identical(parent, here)) {
      break
    }
    here <- parent
  }

  stop("Could not find dominant_cell_states project root from: ", start,
       call. = FALSE)
}

dcs_source_path <- function(..., root = dcs_project_root()) {
  file.path(root, ...)
}

dcs_results_dir <- function(root = dcs_project_root(), create = FALSE) {
  env_results <- Sys.getenv("DOMINANT_CELL_STATES_RESULTS_DIR", unset = "")
  results_dir <- if (nzchar(env_results)) {
    if (grepl("^(/|~)", env_results)) env_results else file.path(root, env_results)
  } else {
    file.path(root, "results")
  }
  results_dir <- normalizePath(results_dir, mustWork = FALSE)
  if (isTRUE(create)) {
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  }
  results_dir
}

dcs_result_path <- function(..., root = dcs_project_root(),
                            create_parent = FALSE) {
  path <- file.path(dcs_results_dir(root, create = create_parent), ...)
  if (isTRUE(create_parent)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  }
  path
}

dcs_description_package <- function(path) {
  desc_file <- file.path(path, "DESCRIPTION")
  if (!file.exists(desc_file)) {
    return(NA_character_)
  }
  desc <- tryCatch(utils::read.dcf(desc_file), error = function(e) NULL)
  if (is.null(desc) || !"Package" %in% colnames(desc)) {
    return(NA_character_)
  }
  unname(desc[1, "Package"])
}

dcs_ribocrypt_repo <- function(start = getwd(), required = TRUE) {
  env_repo <- Sys.getenv("RIBOCRYPT_REPO", unset = "")
  if (nzchar(env_repo)) {
    env_repo <- normalizePath(env_repo, mustWork = FALSE)
    if (identical(dcs_description_package(env_repo), "RiboCrypt")) {
      return(env_repo)
    }
    if (isTRUE(required)) {
      stop("RIBOCRYPT_REPO does not point to the RiboCrypt package root: ",
           env_repo, call. = FALSE)
    }
    return("")
  }

  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (identical(dcs_description_package(here), "RiboCrypt")) {
      return(here)
    }
    parent <- dirname(here)
    if (identical(parent, here)) {
      break
    }
    here <- parent
  }

  if (isTRUE(required)) {
    stop("Could not find local RiboCrypt. Set RIBOCRYPT_REPO=/path/to/RiboCrypt.",
         call. = FALSE)
  }
  ""
}

dcs_load_ribocrypt <- function(repo = dcs_ribocrypt_repo(required = TRUE),
                               quiet = TRUE) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Package 'devtools' is required to load local RiboCrypt.",
         call. = FALSE)
  }
  devtools::load_all(repo, quiet = quiet)
  invisible(repo)
}
