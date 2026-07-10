dominant_description_package <- function(path) {
  desc_file <- file.path(path, "DESCRIPTION")
  if (!file.exists(desc_file)) return(NA_character_)
  desc <- tryCatch(read.dcf(desc_file), error = function(e) NULL)
  if (is.null(desc) || !"Package" %in% colnames(desc)) return(NA_character_)
  unname(desc[1, "Package"])
}

dominant_find_ribocrypt_repo <- function(start = getwd()) {
  env_repo <- Sys.getenv("RIBOCRYPT_REPO", unset = "")
  if (nzchar(env_repo)) {
    env_repo <- normalizePath(env_repo, mustWork = FALSE)
    if (identical(dominant_description_package(env_repo), "RiboCrypt")) {
      return(env_repo)
    }
    return("")
  }

  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (identical(dominant_description_package(here), "RiboCrypt")) {
      return(here)
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  ""
}

dominant_load_ribocrypt_if_available <- function(start = getwd()) {
  repo_root <- dominant_find_ribocrypt_repo(start)
  if (!nzchar(repo_root)) return(invisible(FALSE))
  if (!requireNamespace("devtools", quietly = TRUE)) return(invisible(FALSE))
  devtools::load_all(repo_root, quiet = TRUE)
  invisible(TRUE)
}
