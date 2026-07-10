\name{dcs_pipeline_runner}
\alias{dcs_pipeline_runner}
\alias{dcs_load_pipeline}
\alias{dcs_run_pipeline}
\title{Load and run the dominant-cell-states pipeline}
\usage{
dcs_pipeline_runner(root = dcs_project_root())

dcs_load_pipeline(root = dcs_project_root(), envir = parent.frame())

dcs_run_pipeline(..., root = dcs_project_root())
}
\arguments{
\item{root}{Dominant-cell-states project root.}

\item{envir}{Environment where the pipeline runner should be sourced.}

\item{...}{Arguments passed to \code{run_dominant_cell_state_pipeline()}.}
}
\value{
\code{dcs_pipeline_runner()} returns the runner path,
\code{dcs_load_pipeline()} invisibly returns \code{envir}, and
\code{dcs_run_pipeline()} returns the pipeline result.
}
\description{
Convenience wrappers around \code{scripts/run_dominant_cell_state_pipeline.R}.
The pipeline loads local RiboCrypt through \code{RIBOCRYPT_REPO}.
}
